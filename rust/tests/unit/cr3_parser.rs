use piex::cr3::{
    extract_cr3_orientation, extract_cr3_preview_data, extract_cr3_preview_data_into,
    JpegBoxOffsets, PRVW_JPEG_OFFSETS, THMB_JPEG_OFFSETS, UUID_MOOV, UUID_PRVW,
};
use piex::stream::{ByteSliceStream, Endian};
use piex::tiff::directory::TiffType;
use piex::tiff::tags::tiff::*;
use piex::types::{ImageFormat, PiexError, PreviewImageData, Rational};

// ============================================================================
// Helper Utilities
// ============================================================================

struct BoxBuilder {
    buf: Vec<u8>,
}

impl BoxBuilder {
    fn new() -> Self {
        Self { buf: Vec::new() }
    }

    fn add_box(&mut self, tag: &[u8; 4], payload: &[u8]) {
        let total_size = u32::try_from(8 + payload.len()).expect("box payload fits in u32");
        self.buf.extend_from_slice(&total_size.to_be_bytes());
        self.buf.extend_from_slice(tag);
        self.buf.extend_from_slice(payload);
    }
}

struct TiffWriter {
    buf: Vec<u8>,
    endian: Endian,
}

impl TiffWriter {
    fn new(endian: Endian) -> Self {
        Self { buf: Vec::new(), endian }
    }

    fn write_header(&mut self, first_ifd_offset: u32) {
        match self.endian {
            Endian::Little => {
                self.buf.extend_from_slice(b"II");
                self.buf.extend_from_slice(&42u16.to_le_bytes());
                self.buf.extend_from_slice(&first_ifd_offset.to_le_bytes());
            }
            Endian::Big => {
                self.buf.extend_from_slice(b"MM");
                self.buf.extend_from_slice(&42u16.to_be_bytes());
                self.buf.extend_from_slice(&first_ifd_offset.to_be_bytes());
            }
        }
    }

    fn pack_u16(&self, val: u16) -> [u8; 2] {
        match self.endian {
            Endian::Little => val.to_le_bytes(),
            Endian::Big => val.to_be_bytes(),
        }
    }

    fn pack_u32(&self, val: u32) -> [u8; 4] {
        match self.endian {
            Endian::Little => val.to_le_bytes(),
            Endian::Big => val.to_be_bytes(),
        }
    }

    fn write_entry_inline(&mut self, tag: u16, type_: TiffType, count: u32, inline_bytes: &[u8]) {
        self.buf.extend_from_slice(&self.pack_u16(tag));
        self.buf.extend_from_slice(&self.pack_u16(type_ as u16));
        self.buf.extend_from_slice(&self.pack_u32(count));
        let mut value_field = [0u8; 4];
        let copy_len = inline_bytes.len().min(4);
        value_field[..copy_len].copy_from_slice(&inline_bytes[..copy_len]);
        self.buf.extend_from_slice(&value_field);
    }

    fn write_entry_offset(&mut self, tag: u16, type_: TiffType, count: u32, data_offset: u32) {
        self.buf.extend_from_slice(&self.pack_u16(tag));
        self.buf.extend_from_slice(&self.pack_u16(type_ as u16));
        self.buf.extend_from_slice(&self.pack_u32(count));
        self.buf.extend_from_slice(&self.pack_u32(data_offset));
    }
}

fn create_prvw_payload(width: u16, height: u16, jpeg_data: &[u8]) -> Vec<u8> {
    // PRVW payload relative offsets from box header (8 bytes header):
    // payload offset 6 (width), 8 (height), 12 (jpeg_size), 16 (jpeg_data)
    let mut payload = vec![0u8; 16 + jpeg_data.len()];
    payload[6..8].copy_from_slice(&width.to_be_bytes());
    payload[8..10].copy_from_slice(&height.to_be_bytes());
    let jpeg_len = u32::try_from(jpeg_data.len()).expect("JPEG length fits u32");
    payload[12..16].copy_from_slice(&jpeg_len.to_be_bytes());
    payload[16..].copy_from_slice(jpeg_data);
    payload
}

fn create_thmb_payload(width: u16, height: u16, jpeg_data: &[u8]) -> Vec<u8> {
    // THMB payload relative offsets from box header (8 bytes header):
    // payload offset 4 (width), 6 (height), 8 (jpeg_size), 16 (jpeg_data)
    let mut payload = vec![0u8; 16 + jpeg_data.len()];
    payload[4..6].copy_from_slice(&width.to_be_bytes());
    payload[6..8].copy_from_slice(&height.to_be_bytes());
    let jpeg_len = u32::try_from(jpeg_data.len()).expect("JPEG length fits u32");
    payload[8..12].copy_from_slice(&jpeg_len.to_be_bytes());
    payload[16..].copy_from_slice(jpeg_data);
    payload
}

// ============================================================================
// 1. PRVW Box & JpegBoxOffsets Unit Tests
// ============================================================================

#[test]
fn test_prvw_jpeg_offsets_constants() {
    assert_eq!(
        PRVW_JPEG_OFFSETS,
        JpegBoxOffsets { width: 14, height: 16, jpeg_size: 20, jpeg_data: 24 }
    );
}

#[test]
fn test_prvw_box_extraction() {
    let jpeg_payload = [0xFF, 0xD8, 0xFF, 0xE0, 0x00, 0x10, 0x4A, 0x46]; // JPEG SOI header
    let prvw_payload = create_prvw_payload(1920, 1080, &jpeg_payload);

    let mut builder = BoxBuilder::new();
    builder.add_box(b"PRVW", &prvw_payload);
    builder.add_box(b"mdat", &[]);

    let mut stream = ByteSliceStream::new(&builder.buf);
    let preview_data =
        extract_cr3_preview_data(&mut stream).expect("extract_cr3_preview_data failed");

    assert_eq!(preview_data.preview.width, 1920);
    assert_eq!(preview_data.preview.height, 1080);
    assert_eq!(preview_data.preview.length, jpeg_payload.len() as u32);
    assert_eq!(preview_data.preview.offset, 24); // 8 bytes header + 16 bytes offset
    assert_eq!(preview_data.preview.format, ImageFormat::JpegCompressed);
}

// ============================================================================
// 2. THMB Box Unit Tests
// ============================================================================

#[test]
fn test_thmb_jpeg_offsets_constants() {
    assert_eq!(
        THMB_JPEG_OFFSETS,
        JpegBoxOffsets { width: 12, height: 14, jpeg_size: 16, jpeg_data: 24 }
    );
}

#[test]
fn test_thmb_box_thumbnail_extraction() {
    let prvw_jpeg = [0xFF, 0xD8, 0xFF, 0xE0, 0x00, 0x10];
    let thmb_jpeg = [0xFF, 0xD8, 0xFF, 0xDB, 0x00, 0x43];

    let prvw_payload = create_prvw_payload(1600, 1200, &prvw_jpeg);
    let thmb_payload = create_thmb_payload(160, 120, &thmb_jpeg);

    let mut builder = BoxBuilder::new();
    builder.add_box(b"PRVW", &prvw_payload);
    let prvw_box_len = builder.buf.len();
    builder.add_box(b"THMB", &thmb_payload);
    builder.add_box(b"mdat", &[]);

    let mut stream = ByteSliceStream::new(&builder.buf);
    let preview_data =
        extract_cr3_preview_data(&mut stream).expect("extract_cr3_preview_data failed");

    // Main preview
    assert_eq!(preview_data.preview.width, 1600);
    assert_eq!(preview_data.preview.height, 1200);

    // Thumbnail
    assert_eq!(preview_data.thumbnail.width, 160);
    assert_eq!(preview_data.thumbnail.height, 120);
    assert_eq!(preview_data.thumbnail.length, thmb_jpeg.len() as u32);
    assert_eq!(preview_data.thumbnail.offset, (prvw_box_len + 24) as u32);
    assert_eq!(preview_data.thumbnail.format, ImageFormat::JpegCompressed);
}

// ============================================================================
// 3. UUID Box Container Navigation Tests
// ============================================================================

#[test]
fn test_uuid_moov_and_uuid_prvw_container_navigation() {
    let jpeg_payload = [0xFF, 0xD8, 0xFF, 0xE0, 0x00, 0x10];
    let prvw_payload = create_prvw_payload(1280, 720, &jpeg_payload);

    // Build inner PRVW box
    let mut prvw_builder = BoxBuilder::new();
    prvw_builder.add_box(b"PRVW", &prvw_payload);

    // Build UUID_PRVW box wrapping the inner PRVW box
    // UUID_PRVW header is 16 bytes UUID + 8 bytes padding = 24 bytes payload prefix
    let mut uuid_prvw_payload = Vec::new();
    uuid_prvw_payload.extend_from_slice(&UUID_PRVW);
    uuid_prvw_payload.extend_from_slice(&[0u8; 8]); // 8 bytes padding
    uuid_prvw_payload.extend_from_slice(&prvw_builder.buf);

    let mut uuid_prvw_builder = BoxBuilder::new();
    uuid_prvw_builder.add_box(b"uuid", &uuid_prvw_payload);

    // Build UUID_MOOV box wrapping the UUID_PRVW box
    // UUID_MOOV header is 16 bytes UUID payload prefix
    let mut uuid_moov_payload = Vec::new();
    uuid_moov_payload.extend_from_slice(&UUID_MOOV);
    uuid_moov_payload.extend_from_slice(&uuid_prvw_builder.buf);

    let mut outer_builder = BoxBuilder::new();

    // Add an unknown UUID box before UUID_MOOV to verify skipping
    let mut unknown_uuid_payload = vec![0xABu8; 16]; // Custom unknown UUID
    unknown_uuid_payload.extend_from_slice(&[0x12; 32]); // Dummy payload
    outer_builder.add_box(b"uuid", &unknown_uuid_payload);

    // Add UUID_MOOV box
    outer_builder.add_box(b"uuid", &uuid_moov_payload);
    outer_builder.add_box(b"mdat", &[]);

    let mut stream = ByteSliceStream::new(&outer_builder.buf);
    let preview_data =
        extract_cr3_preview_data(&mut stream).expect("extract_cr3_preview_data failed");

    assert_eq!(preview_data.preview.width, 1280);
    assert_eq!(preview_data.preview.height, 720);
    assert_eq!(preview_data.preview.format, ImageFormat::JpegCompressed);
}

// ============================================================================
// 4. STBL Sub-tree & CRAW/STSZ/CO64 Parsing Tests
// ============================================================================

fn build_stbl_payload(content_type: u16, sample_size_direct: u32, chunk_offset: u64) -> Vec<u8> {
    // 1. CRAW box inside STSD
    // CRAW box size = 100 bytes
    let mut craw_payload = vec![0u8; 92]; // 100 total - 8 header
    craw_payload[24..26].copy_from_slice(&3840u16.to_be_bytes()); // width at +32 (craw offset + 32,
                                                                  // i.e., payload + 24)
    craw_payload[26..28].copy_from_slice(&2160u16.to_be_bytes()); // height at +34 (craw offset +
                                                                  // 34, i.e., payload + 26)
    craw_payload[78..80].copy_from_slice(&content_type.to_be_bytes()); // content_type at +86 (craw
                                                                       // offset + 86, i.e., payload
                                                                       // + 78)

    let mut craw_box = BoxBuilder::new();
    craw_box.add_box(b"CRAW", &craw_payload);

    // STSD box payload: 8 bytes STSD header (version/flags/entry count) + CRAW box
    let mut stsd_payload = vec![0u8; 8];
    stsd_payload.extend_from_slice(&craw_box.buf);

    let mut stsd_box = BoxBuilder::new();
    stsd_box.add_box(b"stsd", &stsd_payload);

    // 2. STSZ box
    // Payload: 4 bytes version/flags, 4 bytes sample_size, 4 bytes count, [4 bytes if count > 0]
    let mut stsz_payload = vec![0u8; 12]; // 4 version/flags + 4 sample_size + 4 count
    stsz_payload[4..8].copy_from_slice(&sample_size_direct.to_be_bytes());
    if sample_size_direct == 0 {
        stsz_payload[8..12].copy_from_slice(&1u32.to_be_bytes()); // count = 1
        stsz_payload.extend_from_slice(&10000u32.to_be_bytes()); // sample size = 10000
    }

    let mut stsz_box = BoxBuilder::new();
    stsz_box.add_box(b"stsz", &stsz_payload);

    // 3. CO64 box
    // Header offset + 4 (count = 1), Header offset + 8 (chunk_offset)
    let mut co64_payload = vec![0u8; 12]; // 4 bytes version/flags, 4 bytes count, 4 bytes
                                          // chunk_offset
    co64_payload[4..8].copy_from_slice(&1u32.to_be_bytes()); // count = 1 at header_offset + 4
    co64_payload[8..12].copy_from_slice(&(chunk_offset as u32).to_be_bytes()); // chunk_offset at
                                                                               // header_offset + 8

    let mut co64_box = BoxBuilder::new();
    co64_box.add_box(b"co64", &co64_payload);

    // Combine into STBL payload
    let mut stbl_payload = Vec::new();
    stbl_payload.extend_from_slice(&stsd_box.buf);
    stbl_payload.extend_from_slice(&stsz_box.buf);
    stbl_payload.extend_from_slice(&co64_box.buf);

    stbl_payload
}

#[test]
fn test_stbl_subtree_mdat_image_calculation() {
    let stbl_payload = build_stbl_payload(3 /* JPEG */, 20480, 16);

    let mut stbl_builder = BoxBuilder::new();
    stbl_builder.add_box(b"stbl", &stbl_payload);

    // Put STBL inside container boxes: MOOV -> TRAK -> MDIA -> MINF -> STBL
    let mut minf = BoxBuilder::new();
    minf.add_box(b"minf", &stbl_builder.buf);

    let mut mdia = BoxBuilder::new();
    mdia.add_box(b"mdia", &minf.buf);

    let mut trak = BoxBuilder::new();
    trak.add_box(b"trak", &mdia.buf);

    let mut moov = BoxBuilder::new();
    moov.add_box(b"moov", &trak.buf);

    let mut outer_builder = BoxBuilder::new();
    outer_builder.buf.extend_from_slice(&moov.buf);

    // After MOOV, place MDAT box
    let mdat_offset = outer_builder.buf.len();
    // Payload of MDAT contains JPEG payload starting at chunk_offset = 16
    let mut mdat_payload = vec![0u8; 16]; // 16 bytes offset prefix inside MDAT
    mdat_payload.extend_from_slice(&[0xFF, 0xD8, 0xFF, 0xE0]); // JPEG SOI
    mdat_payload.resize(16 + 20480, 0u8);

    outer_builder.add_box(b"mdat", &mdat_payload);

    let mut stream = ByteSliceStream::new(&outer_builder.buf);
    let preview_data =
        extract_cr3_preview_data(&mut stream).expect("extract_cr3_preview_data failed");

    assert_eq!(preview_data.preview.width, 3840);
    assert_eq!(preview_data.preview.height, 2160);
    assert_eq!(preview_data.preview.length, 20480);

    // mdat_image offset = mdat header offset (mdat_offset + 8) + chunk_offset (16)
    let expected_offset = u32::try_from(mdat_offset + 8 + 16).unwrap();
    assert_eq!(preview_data.preview.offset, expected_offset);
    assert_eq!(preview_data.preview.format, ImageFormat::JpegCompressed);
}

// ============================================================================
// 5. Preview Image Selection Logic (JPEG vs HEVC Fallback)
// ============================================================================

#[test]
fn test_preview_selection_invalid_jpeg_magic_bytes_falls_back_to_prvw() {
    // 1. STBL specifying JPEG content_type = 3
    let stbl_payload = build_stbl_payload(3 /* JPEG */, 15000, 0);

    let mut stbl_builder = BoxBuilder::new();
    stbl_builder.add_box(b"stbl", &stbl_payload);

    let mut moov = BoxBuilder::new();
    moov.add_box(b"moov", &stbl_builder.buf);

    let mut outer_builder = BoxBuilder::new();
    outer_builder.buf.extend_from_slice(&moov.buf);

    // 2. PRVW box containing valid JPEG
    let jpeg_bytes = [0xFF, 0xD8, 0xFF, 0xE0, 0x00, 0x10];
    let prvw_payload = create_prvw_payload(1920, 1080, &jpeg_bytes);
    outer_builder.add_box(b"PRVW", &prvw_payload);

    // 3. MDAT box containing corrupted non-JPEG magic bytes
    let mdat_payload = vec![0x00, 0x00, 0x00, 0x00, 0x00, 0x00]; // Invalid magic bytes
    outer_builder.add_box(b"mdat", &mdat_payload);

    let mut stream = ByteSliceStream::new(&outer_builder.buf);
    let preview_data =
        extract_cr3_preview_data(&mut stream).expect("extract_cr3_preview_data failed");

    // Must fallback to PRVW image
    assert_eq!(preview_data.preview.width, 1920);
    assert_eq!(preview_data.preview.height, 1080);
    assert_eq!(preview_data.preview.format, ImageFormat::JpegCompressed);
}

/// Verifies that HEVC mdat image IS selected as preview image with format
/// ImageFormat::HevcCompressed.
#[test]
fn test_preview_selection_hevc_mdat_selected_as_preview() {
    // 1. STBL specifying HEVC (content_type = 4)
    let stbl_payload = build_stbl_payload(4 /* HEVC */, 15000, 0);

    let mut stbl_builder = BoxBuilder::new();
    stbl_builder.add_box(b"stbl", &stbl_payload);

    let mut moov = BoxBuilder::new();
    moov.add_box(b"moov", &stbl_builder.buf);

    let mut outer_builder = BoxBuilder::new();
    outer_builder.buf.extend_from_slice(&moov.buf);

    // 2. PRVW box containing valid JPEG
    let jpeg_bytes = [0xFF, 0xD8, 0xFF, 0xE0, 0x00, 0x10];
    let prvw_payload = create_prvw_payload(1920, 1080, &jpeg_bytes);
    outer_builder.add_box(b"PRVW", &prvw_payload);

    // 3. MDAT box containing HEVC non-JPEG payload
    let mdat_payload = vec![0x00, 0x00, 0x00, 0x01, 0x40, 0x01]; // HEVC NAL unit header
    outer_builder.add_box(b"mdat", &mdat_payload);

    let mut stream = ByteSliceStream::new(&outer_builder.buf);
    let preview_data =
        extract_cr3_preview_data(&mut stream).expect("extract_cr3_preview_data failed");

    // HEVC mdat image IS selected as preview (with width 3840, height 2160, and HevcCompressed
    // format).
    assert_eq!(preview_data.preview.width, 3840);
    assert_eq!(preview_data.preview.height, 2160);
    assert_eq!(preview_data.preview.format, ImageFormat::HevcCompressed);
}

// ============================================================================
// 6. EXIF Metadata Parsing (CMT1, CMT2 & extract_cr3_orientation)
// ============================================================================

fn create_cmt1_payload() -> Vec<u8> {
    let mut writer = TiffWriter::new(Endian::Little);
    writer.write_header(8);

    // IFD0 with 5 entries
    writer.buf.extend_from_slice(&writer.pack_u16(5));

    // Make ("Canon\0") at offset 100
    writer.write_entry_offset(TIFF_TAG_MAKE as u16, TiffType::Ascii, 6, 100);
    // Model ("EOS R5\0") at offset 120
    writer.write_entry_offset(TIFF_TAG_MODEL as u16, TiffType::Ascii, 7, 120);
    // Orientation = 6 (Rotate 90 CW)
    writer.write_entry_inline(TIFF_TAG_ORIENTATION as u16, TiffType::Short, 1, &6u16.to_le_bytes());
    // Full width = 8192
    writer.write_entry_inline(
        TIFF_TAG_IMAGE_WIDTH as u16,
        TiffType::Long,
        1,
        &8192u32.to_le_bytes(),
    );
    // Full height = 5464
    writer.write_entry_inline(
        TIFF_TAG_IMAGE_LENGTH as u16,
        TiffType::Long,
        1,
        &5464u32.to_le_bytes(),
    );

    // Next IFD = 0
    writer.buf.extend_from_slice(&writer.pack_u32(0));

    if writer.buf.len() < 106 {
        writer.buf.resize(106, 0);
    }
    writer.buf[100..106].copy_from_slice(b"Canon\0");

    if writer.buf.len() < 127 {
        writer.buf.resize(127, 0);
    }
    writer.buf[120..127].copy_from_slice(b"EOS R5\0");

    writer.buf
}

fn create_cmt2_payload() -> Vec<u8> {
    let mut writer = TiffWriter::new(Endian::Little);
    writer.write_header(8);

    // EXIF IFD with 5 entries
    writer.buf.extend_from_slice(&writer.pack_u16(5));

    // DateTimeOriginal ("2024:01:15 10:20:30\0") at offset 100
    writer.write_entry_offset(EXIF_TAG_DATE_TIME_ORIGINAL as u16, TiffType::Ascii, 20, 100);
    // ExposureTime (1/500) at offset 130
    writer.write_entry_offset(EXIF_TAG_EXPOSURE_TIME as u16, TiffType::Rational, 1, 130);
    // FNumber (40/10 = 4.0) at offset 140
    writer.write_entry_offset(EXIF_TAG_FNUMBER as u16, TiffType::Rational, 1, 140);
    // FocalLength (85/1) at offset 150
    writer.write_entry_offset(EXIF_TAG_FOCAL_LENGTH as u16, TiffType::Rational, 1, 150);
    // ISO (400)
    writer.write_entry_inline(EXIF_TAG_ISO_SPEED as u16, TiffType::Long, 1, &400u32.to_le_bytes());

    // Next IFD = 0
    writer.buf.extend_from_slice(&writer.pack_u32(0));

    if writer.buf.len() < 120 {
        writer.buf.resize(120, 0);
    }
    writer.buf[100..120].copy_from_slice(b"2024:01:15 10:20:30\0");

    if writer.buf.len() < 138 {
        writer.buf.resize(138, 0);
    }
    writer.buf[130..134].copy_from_slice(&1u32.to_le_bytes());
    writer.buf[134..138].copy_from_slice(&500u32.to_le_bytes());

    if writer.buf.len() < 148 {
        writer.buf.resize(148, 0);
    }
    writer.buf[140..144].copy_from_slice(&40u32.to_le_bytes());
    writer.buf[144..148].copy_from_slice(&10u32.to_le_bytes());

    if writer.buf.len() < 158 {
        writer.buf.resize(158, 0);
    }
    writer.buf[150..154].copy_from_slice(&85u32.to_le_bytes());
    writer.buf[154..158].copy_from_slice(&1u32.to_le_bytes());

    writer.buf
}

#[test]
fn test_cmt1_and_cmt2_exif_metadata_extraction() {
    let cmt1_bytes = create_cmt1_payload();
    let cmt2_bytes = create_cmt2_payload();

    let jpeg_bytes = [0xFF, 0xD8, 0xFF, 0xE0, 0x00, 0x10];
    let prvw_payload = create_prvw_payload(1920, 1080, &jpeg_bytes);

    let mut builder = BoxBuilder::new();
    builder.add_box(b"CMT1", &cmt1_bytes);
    builder.add_box(b"CMT2", &cmt2_bytes);
    builder.add_box(b"PRVW", &prvw_payload);
    builder.add_box(b"mdat", &[]);

    let mut stream = ByteSliceStream::new(&builder.buf);
    let preview_data =
        extract_cr3_preview_data(&mut stream).expect("extract_cr3_preview_data failed");

    // CMT1 metadata
    assert_eq!(preview_data.maker, b"Canon\0");
    assert_eq!(preview_data.model, b"EOS R5\0");
    assert_eq!(preview_data.exif_orientation, 6);
    assert_eq!(preview_data.full_width, 8192);
    assert_eq!(preview_data.full_height, 5464);

    // CMT2 metadata
    assert_eq!(preview_data.date_time, b"2024:01:15 10:20:30\0");
    assert_eq!(preview_data.exposure_time, Rational::new(1, 500));
    assert_eq!(preview_data.fnumber, Rational::new(40, 10));
    assert_eq!(preview_data.focal_length, Rational::new(85, 1));
    assert_eq!(preview_data.iso, 400);
}

#[test]
fn test_extract_cr3_orientation_standalone() {
    let cmt1_bytes = create_cmt1_payload();

    let mut builder = BoxBuilder::new();
    builder.add_box(b"CMT1", &cmt1_bytes);

    let mut stream = ByteSliceStream::new(&builder.buf);
    let orientation = extract_cr3_orientation(&mut stream).expect("extract_cr3_orientation failed");

    assert_eq!(orientation, 6);
}

// ============================================================================
// Parity Remediation Unit Tests (DIV-04, DIV-10, DIV-11)
// ============================================================================

/// DIV-04: Verifies that missing mdat box causes extraction failure without speculative PRVW stream
/// traversal retry.
#[test]
fn test_cr3_missing_mdat_box_fails_parity() {
    let jpeg_bytes = [0xFF, 0xD8, 0xFF, 0xE0, 0x00, 0x10];
    let prvw_payload = create_prvw_payload(1920, 1080, &jpeg_bytes);

    let mut builder = BoxBuilder::new();
    builder.add_box(b"moov", &[0u8; 8]);
    builder.add_box(b"PRVW", &prvw_payload);
    // Note: No mdat box is added

    let mut stream = ByteSliceStream::new(&builder.buf);
    let res = extract_cr3_preview_data(&mut stream);
    assert!(res.is_err(), "Extraction must fail when mdat box is missing");
}

/// DIV-10: Verifies single box header read per box during sequential traversal.
#[test]
fn test_div10_single_box_header_read_per_box() {
    let jpeg_bytes = [0xFF, 0xD8, 0xFF, 0xE0, 0x00, 0x10];
    let prvw_payload = create_prvw_payload(1920, 1080, &jpeg_bytes);

    let mut builder = BoxBuilder::new();
    builder.add_box(b"ftyp", b"crx \0\0\0\x01crx \0\0\0\0");
    builder.add_box(b"PRVW", &prvw_payload);
    builder.add_box(b"mdat", &[]);

    let mut stream = ByteSliceStream::new(&builder.buf);
    let preview_data =
        extract_cr3_preview_data(&mut stream).expect("extract_cr3_preview_data failed");
    assert_eq!(preview_data.preview.width, 1920);
    assert_eq!(preview_data.preview.height, 1080);
}

/// DIV-10: Verifies that a corrupted box header at non-zero offset propagates the error.
#[test]
fn test_div10_corrupted_box_header_at_nonzero_offset_propagates_error() {
    let mut builder = BoxBuilder::new();
    builder.add_box(b"ftyp", b"crx \0\0\0\x01crx \0\0\0\0"); // 24 bytes
                                                             // Append 4 truncated bytes (incomplete
                                                             // box header at offset 24)
    builder.buf.extend_from_slice(&32u32.to_be_bytes());

    let mut stream = ByteSliceStream::new(&builder.buf);
    let res = extract_cr3_preview_data(&mut stream);
    assert!(
        matches!(res, Err(PiexError::InvalidBoxHeader { offset: 24 })),
        "Expected InvalidBoxHeader at offset 24, got {:?}",
        res
    );
}

/// DIV-11: Verifies that corrupted CMT1 with corrupted/empty TIFF directory returns error.
#[test]
fn test_cmt1_corrupted_empty_tiff_directory_returns_err() {
    // Corrupted CMT1 payload with entry count 0xFFFF
    let mut cmt1_payload = Vec::new();
    cmt1_payload.extend_from_slice(b"MM\x00\x2A\x00\x00\x00\x08");
    cmt1_payload.extend_from_slice(&[0xFF, 0xFF]);

    let jpeg_bytes = [0xFF, 0xD8, 0xFF, 0xE0, 0x00, 0x10];
    let prvw_payload = create_prvw_payload(1920, 1080, &jpeg_bytes);

    let mut builder = BoxBuilder::new();
    builder.add_box(b"CMT1", &cmt1_payload);
    builder.add_box(b"PRVW", &prvw_payload);
    builder.add_box(b"mdat", &[]);

    let mut stream = ByteSliceStream::new(&builder.buf);
    let res = extract_cr3_preview_data(&mut stream);
    assert!(res.is_err(), "Expected error on corrupted CMT1 TIFF directory");

    let orient_res = extract_cr3_orientation(&mut stream);
    assert!(orient_res.is_err(), "Expected error on orientation with corrupted CMT1");
}

/// DIV-11: Verifies that corrupted CMT2 with corrupted/empty TIFF directory returns error.
#[test]
fn test_cmt2_corrupted_empty_tiff_directory_returns_err() {
    let cmt1_bytes = create_cmt1_payload();

    // Corrupted CMT2 payload with entry count 0xFFFF
    let mut cmt2_payload = Vec::new();
    cmt2_payload.extend_from_slice(b"MM\x00\x2A\x00\x00\x00\x08");
    cmt2_payload.extend_from_slice(&[0xFF, 0xFF]);

    let jpeg_bytes = [0xFF, 0xD8, 0xFF, 0xE0, 0x00, 0x10];
    let prvw_payload = create_prvw_payload(1920, 1080, &jpeg_bytes);

    let mut builder = BoxBuilder::new();
    builder.add_box(b"CMT1", &cmt1_bytes);
    builder.add_box(b"CMT2", &cmt2_payload);
    builder.add_box(b"PRVW", &prvw_payload);
    builder.add_box(b"mdat", &[]);

    let mut stream = ByteSliceStream::new(&builder.buf);
    let res = extract_cr3_preview_data(&mut stream);
    assert!(res.is_err(), "Expected error on corrupted CMT2 TIFF directory");
}

/// DIV-11: Verifies that missing stsd in stbl returns CorruptData error.
#[test]
fn test_stbl_missing_stsd_box_returns_err() {
    let mut stbl_builder = BoxBuilder::new();
    stbl_builder.add_box(b"stbl", &[]); // Empty stbl, no stsd

    let mut minf = BoxBuilder::new();
    minf.add_box(b"minf", &stbl_builder.buf);
    let mut mdia = BoxBuilder::new();
    mdia.add_box(b"mdia", &minf.buf);
    let mut trak = BoxBuilder::new();
    trak.add_box(b"trak", &mdia.buf);
    let mut moov = BoxBuilder::new();
    moov.add_box(b"moov", &trak.buf);

    let jpeg_bytes = [0xFF, 0xD8, 0xFF, 0xE0, 0x00, 0x10];
    let prvw_payload = create_prvw_payload(1920, 1080, &jpeg_bytes);

    let mut outer_builder = BoxBuilder::new();
    outer_builder.buf.extend_from_slice(&moov.buf);
    outer_builder.add_box(b"PRVW", &prvw_payload);
    outer_builder.add_box(b"mdat", &[]);

    let mut stream = ByteSliceStream::new(&outer_builder.buf);
    let res = extract_cr3_preview_data(&mut stream);
    assert!(res.is_err(), "Expected error when stsd is missing from stbl");
}

/// DIV-11: Verifies that corrupted/missing box inside stsd returns CorruptData error.
#[test]
fn test_stbl_missing_craw_box_inside_stsd_returns_err() {
    let mut stsd_box = BoxBuilder::new();
    // 8 bytes stsd header + invalid child box header with size = 2 (< 8 bytes minimum)
    stsd_box.add_box(b"stsd", &[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, b'C', b'R', b'A', b'W']);

    let mut stbl_builder = BoxBuilder::new();
    stbl_builder.add_box(b"stbl", &stsd_box.buf);

    let mut minf = BoxBuilder::new();
    minf.add_box(b"minf", &stbl_builder.buf);
    let mut mdia = BoxBuilder::new();
    mdia.add_box(b"mdia", &minf.buf);
    let mut trak = BoxBuilder::new();
    trak.add_box(b"trak", &mdia.buf);
    let mut moov = BoxBuilder::new();
    moov.add_box(b"moov", &trak.buf);

    let jpeg_bytes = [0xFF, 0xD8, 0xFF, 0xE0, 0x00, 0x10];
    let prvw_payload = create_prvw_payload(1920, 1080, &jpeg_bytes);

    let mut outer_builder = BoxBuilder::new();
    outer_builder.buf.extend_from_slice(&moov.buf);
    outer_builder.add_box(b"PRVW", &prvw_payload);
    outer_builder.add_box(b"mdat", &[]);

    let mut stream = ByteSliceStream::new(&outer_builder.buf);
    let res = extract_cr3_preview_data(&mut stream);
    assert!(res.is_err(), "Expected error when craw is missing/corrupted inside stsd");
}

/// DIV-11: Verifies that missing stsz in stbl returns CorruptData error.
#[test]
fn test_stbl_missing_stsz_box_returns_err() {
    let mut craw_payload = vec![0u8; 92];
    craw_payload[24..26].copy_from_slice(&3840u16.to_be_bytes());
    craw_payload[26..28].copy_from_slice(&2160u16.to_be_bytes());
    craw_payload[78..80].copy_from_slice(&3u16.to_be_bytes()); // JPEG

    let mut craw_box = BoxBuilder::new();
    craw_box.add_box(b"CRAW", &craw_payload);

    let mut stsd_payload = vec![0u8; 8];
    stsd_payload.extend_from_slice(&craw_box.buf);

    let mut stsd_box = BoxBuilder::new();
    stsd_box.add_box(b"stsd", &stsd_payload);

    // stbl contains only stsd (no stsz)
    let mut stbl_builder = BoxBuilder::new();
    stbl_builder.add_box(b"stbl", &stsd_box.buf);

    let mut minf = BoxBuilder::new();
    minf.add_box(b"minf", &stbl_builder.buf);
    let mut mdia = BoxBuilder::new();
    mdia.add_box(b"mdia", &minf.buf);
    let mut trak = BoxBuilder::new();
    trak.add_box(b"trak", &mdia.buf);
    let mut moov = BoxBuilder::new();
    moov.add_box(b"moov", &trak.buf);

    let jpeg_bytes = [0xFF, 0xD8, 0xFF, 0xE0, 0x00, 0x10];
    let prvw_payload = create_prvw_payload(1920, 1080, &jpeg_bytes);

    let mut outer_builder = BoxBuilder::new();
    outer_builder.buf.extend_from_slice(&moov.buf);
    outer_builder.add_box(b"PRVW", &prvw_payload);
    outer_builder.add_box(b"mdat", &[]);

    let mut stream = ByteSliceStream::new(&outer_builder.buf);
    let res = extract_cr3_preview_data(&mut stream);
    assert!(res.is_err(), "Expected error when stsz is missing in stbl");
}

/// DIV-11: Verifies that missing co64 in stbl returns CorruptData error.
#[test]
fn test_stbl_missing_co64_box_returns_err() {
    let mut craw_payload = vec![0u8; 92];
    craw_payload[24..26].copy_from_slice(&3840u16.to_be_bytes());
    craw_payload[26..28].copy_from_slice(&2160u16.to_be_bytes());
    craw_payload[78..80].copy_from_slice(&3u16.to_be_bytes()); // JPEG

    let mut craw_box = BoxBuilder::new();
    craw_box.add_box(b"CRAW", &craw_payload);

    let mut stsd_payload = vec![0u8; 8];
    stsd_payload.extend_from_slice(&craw_box.buf);

    let mut stsd_box = BoxBuilder::new();
    stsd_box.add_box(b"stsd", &stsd_payload);

    let mut stsz_payload = vec![0u8; 12];
    stsz_payload[4..8].copy_from_slice(&10000u32.to_be_bytes());
    let mut stsz_box = BoxBuilder::new();
    stsz_box.add_box(b"stsz", &stsz_payload);

    // stbl contains stsd and stsz, but no co64
    let mut stbl_payload = Vec::new();
    stbl_payload.extend_from_slice(&stsd_box.buf);
    stbl_payload.extend_from_slice(&stsz_box.buf);

    let mut stbl_builder = BoxBuilder::new();
    stbl_builder.add_box(b"stbl", &stbl_payload);

    let mut minf = BoxBuilder::new();
    minf.add_box(b"minf", &stbl_builder.buf);
    let mut mdia = BoxBuilder::new();
    mdia.add_box(b"mdia", &minf.buf);
    let mut trak = BoxBuilder::new();
    trak.add_box(b"trak", &mdia.buf);
    let mut moov = BoxBuilder::new();
    moov.add_box(b"moov", &trak.buf);

    let jpeg_bytes = [0xFF, 0xD8, 0xFF, 0xE0, 0x00, 0x10];
    let prvw_payload = create_prvw_payload(1920, 1080, &jpeg_bytes);

    let mut outer_builder = BoxBuilder::new();
    outer_builder.buf.extend_from_slice(&moov.buf);
    outer_builder.add_box(b"PRVW", &prvw_payload);
    outer_builder.add_box(b"mdat", &[]);

    let mut stream = ByteSliceStream::new(&outer_builder.buf);
    let res = extract_cr3_preview_data(&mut stream);
    assert!(res.is_err(), "Expected error when co64 is missing in stbl");
}

// ============================================================================
// 7. Error Handling & Malformed Box Structure Tests
// ============================================================================

#[test]
fn test_error_handling_empty_stream() {
    let buf = [];
    let mut stream = ByteSliceStream::new(&buf);
    let res = extract_cr3_preview_data(&mut stream);
    assert!(matches!(res, Err(PiexError::GenericFail(_))));
}

#[test]
fn test_error_handling_truncated_box_header() {
    let buf = [0x00, 0x00, 0x00, 0x20, 0x43]; // 5 bytes < 8 bytes box header
    let mut stream = ByteSliceStream::new(&buf);
    let res = extract_cr3_preview_data(&mut stream);
    assert!(matches!(res, Err(PiexError::InvalidBoxHeader { offset: 0 })));
}

#[test]
fn test_error_handling_stsz_invalid_sample_count() {
    // STSZ with sample_size = 0, but count != 1 (count = 2)
    let mut stsz_payload = vec![0u8; 12];
    stsz_payload[4..8].copy_from_slice(&0u32.to_be_bytes()); // sample_size = 0
    stsz_payload[8..12].copy_from_slice(&2u32.to_be_bytes()); // count = 2 (Invalid!)

    let mut stsz_box = BoxBuilder::new();
    stsz_box.add_box(b"stsz", &stsz_payload);

    let mut craw_payload = vec![0u8; 92];
    craw_payload[24..26].copy_from_slice(&100u16.to_be_bytes());
    craw_payload[26..28].copy_from_slice(&100u16.to_be_bytes());
    craw_payload[78..80].copy_from_slice(&3u16.to_be_bytes()); // JPEG

    let mut craw_box = BoxBuilder::new();
    craw_box.add_box(b"CRAW", &craw_payload);

    let mut stsd_payload = vec![0u8; 8];
    stsd_payload.extend_from_slice(&craw_box.buf);

    let mut stsd_box = BoxBuilder::new();
    stsd_box.add_box(b"stsd", &stsd_payload);

    let mut co64_payload = vec![0u8; 16];
    co64_payload[4..8].copy_from_slice(&1u32.to_be_bytes());
    co64_payload[8..16].copy_from_slice(&0u64.to_be_bytes());

    let mut co64_box = BoxBuilder::new();
    co64_box.add_box(b"co64", &co64_payload);

    let mut stbl_payload = Vec::new();
    stbl_payload.extend_from_slice(&stsd_box.buf);
    stbl_payload.extend_from_slice(&stsz_box.buf);
    stbl_payload.extend_from_slice(&co64_box.buf);

    let mut stbl_builder = BoxBuilder::new();
    stbl_builder.add_box(b"stbl", &stbl_payload);

    let mut moov = BoxBuilder::new();
    moov.add_box(b"moov", &stbl_builder.buf);

    let mut stream = ByteSliceStream::new(&moov.buf);
    let res = extract_cr3_preview_data(&mut stream);

    assert!(matches!(res, Err(PiexError::CorruptData(_))));
}

#[test]
fn test_error_handling_co64_invalid_entry_count() {
    // CO64 with count != 1 (count = 0)
    let stsz_payload = vec![0u8; 12]; // sample_size = 0, count = 0 -> direct sample_size 0

    let mut stsz_box = BoxBuilder::new();
    stsz_box.add_box(b"stsz", &stsz_payload);

    let mut craw_payload = vec![0u8; 92];
    craw_payload[24..26].copy_from_slice(&100u16.to_be_bytes());
    craw_payload[26..28].copy_from_slice(&100u16.to_be_bytes());
    craw_payload[78..80].copy_from_slice(&3u16.to_be_bytes());

    let mut craw_box = BoxBuilder::new();
    craw_box.add_box(b"CRAW", &craw_payload);

    let mut stsd_payload = vec![0u8; 8];
    stsd_payload.extend_from_slice(&craw_box.buf);

    let mut stsd_box = BoxBuilder::new();
    stsd_box.add_box(b"stsd", &stsd_payload);

    // CO64 with count = 0 at payload offset +4
    let mut co64_payload = vec![0u8; 16];
    co64_payload[4..8].copy_from_slice(&0u32.to_be_bytes()); // count = 0 (Invalid!)

    let mut co64_box = BoxBuilder::new();
    co64_box.add_box(b"co64", &co64_payload);

    let mut stbl_payload = Vec::new();
    stbl_payload.extend_from_slice(&stsd_box.buf);
    stbl_payload.extend_from_slice(&stsz_box.buf);
    stbl_payload.extend_from_slice(&co64_box.buf);

    let mut stbl_builder = BoxBuilder::new();
    stbl_builder.add_box(b"stbl", &stbl_payload);

    let mut moov = BoxBuilder::new();
    moov.add_box(b"moov", &stbl_builder.buf);

    let mut stream = ByteSliceStream::new(&moov.buf);
    let res = extract_cr3_preview_data(&mut stream);

    assert!(matches!(res, Err(PiexError::CorruptData(_))));
}

fn create_valid_jpeg_payload() -> Vec<u8> {
    vec![0xFF, 0xD8, 0xFF, 0xE0, 0x00, 0x10, 0x4A, 0x46, 0x00, 0x01, 0x01, 0x01]
}

fn build_stbl_tree(
    content_type: u16,
    sample_size: u32,
    chunk_offset: u32,
    craw_tag: &[u8; 4],
) -> Vec<u8> {
    let mut craw_payload = vec![0u8; 92];
    craw_payload[24..26].copy_from_slice(&3840u16.to_be_bytes());
    craw_payload[26..28].copy_from_slice(&2160u16.to_be_bytes());
    craw_payload[78..80].copy_from_slice(&content_type.to_be_bytes());

    let mut craw_builder = BoxBuilder::new();
    craw_builder.add_box(craw_tag, &craw_payload);

    let mut stsd_payload = vec![0u8; 8];
    stsd_payload.extend_from_slice(&craw_builder.buf);
    let mut stsd_builder = BoxBuilder::new();
    stsd_builder.add_box(b"stsd", &stsd_payload);

    let mut stsz_payload = vec![0u8; 12];
    stsz_payload[4..8].copy_from_slice(&sample_size.to_be_bytes());
    let mut stsz_builder = BoxBuilder::new();
    stsz_builder.add_box(b"stsz", &stsz_payload);

    let mut co64_payload = vec![0u8; 12];
    co64_payload[4..8].copy_from_slice(&1u32.to_be_bytes()); // count = 1
    co64_payload[8..12].copy_from_slice(&chunk_offset.to_be_bytes());
    let mut co64_builder = BoxBuilder::new();
    co64_builder.add_box(b"co64", &co64_payload);

    let mut stbl_payload = Vec::new();
    stbl_payload.extend_from_slice(&stsd_builder.buf);
    stbl_payload.extend_from_slice(&stsz_builder.buf);
    stbl_payload.extend_from_slice(&co64_builder.buf);

    let mut stbl_builder = BoxBuilder::new();
    stbl_builder.add_box(b"stbl", &stbl_payload);

    let mut minf_builder = BoxBuilder::new();
    minf_builder.add_box(b"minf", &stbl_builder.buf);

    let mut mdia_builder = BoxBuilder::new();
    mdia_builder.add_box(b"mdia", &minf_builder.buf);

    let mut trak_builder = BoxBuilder::new();
    trak_builder.add_box(b"trak", &mdia_builder.buf);

    let mut moov_builder = BoxBuilder::new();
    moov_builder.add_box(b"moov", &trak_builder.buf);

    moov_builder.buf
}

// ============================================================================
// 1. Truncated Box Headers & Box Length Boundary Cases
// ============================================================================

#[test]
fn test_truncated_box_headers_all_sub_8_lengths() {
    for len in 0..8 {
        let buf = vec![0x41u8; len];
        let mut stream = ByteSliceStream::new(&buf);
        let res = extract_cr3_preview_data(&mut stream);
        assert!(res.is_err(), "Stream of length {len} must fail extraction");
    }
}

#[test]
fn test_truncated_64bit_extended_box_headers() {
    for ext_len in 8..16 {
        let mut buf = vec![0x00, 0x00, 0x00, 0x01, b'm', b'd', b'a', b't'];
        buf.resize(ext_len, 0x00);
        let mut stream = ByteSliceStream::new(&buf);
        let res = extract_cr3_preview_data(&mut stream);
        assert!(res.is_err(), "Extended box of length {ext_len} must fail");
    }
}

#[test]
fn test_invalid_box_sizes_zero_and_underflow() {
    // 1. Length 0 box
    let buf_zero = vec![0, 0, 0, 0, b'f', b't', b'y', b'p'];
    let mut stream = ByteSliceStream::new(&buf_zero);
    let res = extract_cr3_preview_data(&mut stream);
    assert!(res.is_err());

    // 2. Length < 8 (e.g. 4)
    let buf_four = vec![0, 0, 0, 4, b'f', b't', b'y', b'p'];
    let mut stream = ByteSliceStream::new(&buf_four);
    let res = extract_cr3_preview_data(&mut stream);
    assert!(res.is_err());

    // 3. Extended box with length < 16 (e.g. 10)
    let mut buf_ext_short = vec![0, 0, 0, 1, b'm', b'd', b'a', b't'];
    buf_ext_short.extend_from_slice(&10u64.to_be_bytes());
    let mut stream = ByteSliceStream::new(&buf_ext_short);
    let res = extract_cr3_preview_data(&mut stream);
    assert!(res.is_err());
}

// ============================================================================
// 2. Partial / Truncated Box Payloads
// ============================================================================

#[test]
fn test_partial_prvw_box_payloads() {
    let jpeg = create_valid_jpeg_payload();
    let full_prvw = create_prvw_payload(1920, 1080, &jpeg);

    // Truncate at end of stream where no subsequent box can be aliased
    for cut in [0, 4, 8, 14, 16, 20, 23] {
        let mut builder = BoxBuilder::new();
        builder.add_box(b"ftyp", b"crx \0\0\0\x01crx \0\0\0\0");
        builder.add_box(b"PRVW", &full_prvw[..cut.min(full_prvw.len())]);

        let mut stream = ByteSliceStream::new(&builder.buf);
        let res = extract_cr3_preview_data(&mut stream);
        assert!(res.is_err(), "Truncated PRVW at {cut} must fail");
    }
}

#[test]
fn test_partial_thmb_box_payloads() {
    let jpeg = create_valid_jpeg_payload();
    let thmb_payload = create_thmb_payload(160, 120, &jpeg);

    // Truncate at end of stream
    for cut in [0, 4, 8, 12, 14, 16, 20, 23] {
        let mut builder = BoxBuilder::new();
        builder.add_box(b"ftyp", b"crx \0\0\0\x01crx \0\0\0\0");
        builder.add_box(b"THMB", &thmb_payload[..cut.min(thmb_payload.len())]);

        let mut stream = ByteSliceStream::new(&builder.buf);
        let res = extract_cr3_preview_data(&mut stream);
        assert!(res.is_err(), "Truncated THMB at {cut} must fail");
    }
}

#[test]
fn test_partial_cmt1_and_cmt2_box_payloads() {
    let cmt1 = create_cmt1_payload();

    // Corrupted CMT1 with invalid entry count at end of stream
    for cut in [0, 4, 8, 12, 20, 40, 80] {
        let mut builder = BoxBuilder::new();
        builder.add_box(b"ftyp", b"crx \0\0\0\x01crx \0\0\0\0");
        builder.add_box(b"CMT1", &cmt1[..cut.min(cmt1.len())]);

        let mut stream = ByteSliceStream::new(&builder.buf);
        let res = extract_cr3_orientation(&mut stream);
        assert!(res.is_err(), "Truncated CMT1 at {cut} must fail orientation");
    }

    // Corrupted CMT2 with entry count 0xFFFF
    let mut bad_cmt2 = Vec::new();
    bad_cmt2.extend_from_slice(b"MM\x00\x2A\x00\x00\x00\x08");
    bad_cmt2.extend_from_slice(&[0xFF, 0xFF]);

    let mut builder = BoxBuilder::new();
    builder.add_box(b"ftyp", b"crx \0\0\0\x01crx \0\0\0\0");
    builder.add_box(b"CMT1", &cmt1);
    builder.add_box(b"CMT2", &bad_cmt2);
    builder.add_box(b"PRVW", &create_prvw_payload(1920, 1080, &create_valid_jpeg_payload()));
    builder.add_box(b"mdat", &create_valid_jpeg_payload());

    let mut stream = ByteSliceStream::new(&builder.buf);
    let res = extract_cr3_preview_data(&mut stream);
    assert!(res.is_err(), "Corrupted CMT2 must fail extraction");
}

// ============================================================================
// 3. Corrupted UUID Boxes
// ============================================================================

#[test]
fn test_corrupted_uuid_moov_and_uuid_prvw() {
    let jpeg = create_valid_jpeg_payload();

    // UUID_MOOV with truncated UUID payload (< 16 bytes)
    for len in [0, 4, 8, 15] {
        let mut builder = BoxBuilder::new();
        builder.add_box(b"ftyp", b"crx \0\0\0\x01crx \0\0\0\0");
        builder.add_box(b"uuid", &UUID_MOOV[..len]);
        builder.add_box(b"mdat", &jpeg);

        let mut stream = ByteSliceStream::new(&builder.buf);
        let res = extract_cr3_preview_data(&mut stream);
        assert!(res.is_err(), "UUID_MOOV truncated to {len} must fail");
    }

    // UUID_PRVW with truncated UUID payload (< 24 bytes)
    for len in [0, 8, 15, 16, 20, 23] {
        let mut payload = Vec::new();
        payload.extend_from_slice(&UUID_PRVW);
        payload.extend_from_slice(&[0u8; 8]);

        let mut builder = BoxBuilder::new();
        builder.add_box(b"ftyp", b"crx \0\0\0\x01crx \0\0\0\0");
        builder.add_box(b"uuid", &payload[..len.min(payload.len())]);
        builder.add_box(b"mdat", &jpeg);

        let mut stream = ByteSliceStream::new(&builder.buf);
        let res = extract_cr3_preview_data(&mut stream);
        assert!(res.is_err(), "UUID_PRVW truncated to {len} must fail");
    }

    // Unknown UUID box is safely skipped
    let mut unknown_payload = vec![0xEFu8; 16];
    unknown_payload.extend_from_slice(b"arbitrary_vendor_payload");
    let prvw_payload = create_prvw_payload(1920, 1080, &jpeg);

    let mut builder = BoxBuilder::new();
    builder.add_box(b"ftyp", b"crx \0\0\0\x01crx \0\0\0\0");
    builder.add_box(b"uuid", &unknown_payload);
    builder.add_box(b"PRVW", &prvw_payload);
    builder.add_box(b"mdat", &jpeg);

    let mut stream = ByteSliceStream::new(&builder.buf);
    let data = extract_cr3_preview_data(&mut stream).expect("Unknown UUID must be skipped");
    assert_eq!(data.preview.width, 1920);
    assert_eq!(data.preview.height, 1080);
}

// ============================================================================
// 4. Malformed Sample Tables (STBL)
// ============================================================================

#[test]
fn test_stbl_non_craw_tag_skipped_gracefully() {
    let jpeg = create_valid_jpeg_payload();
    let prvw_payload = create_prvw_payload(1920, 1080, &jpeg);

    // STBL with CTMD (Canon Timed Metadata) instead of CRAW
    let moov = build_stbl_tree(3, 10000, 0, b"CTMD");

    let mut builder = BoxBuilder::new();
    builder.add_box(b"ftyp", b"crx \0\0\0\x01crx \0\0\0\0");
    builder.buf.extend_from_slice(&moov);
    builder.add_box(b"PRVW", &prvw_payload);
    builder.add_box(b"mdat", &jpeg);

    let mut stream = ByteSliceStream::new(&builder.buf);
    let data = extract_cr3_preview_data(&mut stream).expect("CTMD in stsd must fall back to PRVW");
    assert_eq!(data.preview.width, 1920);
    assert_eq!(data.preview.height, 1080);
}

#[test]
fn test_stbl_non_jpeg_non_hevc_content_types_fallback_to_prvw() {
    let jpeg = create_valid_jpeg_payload();
    let prvw_payload = create_prvw_payload(1920, 1080, &jpeg);

    for invalid_ct in [0, 1, 2, 5, 6, 100, 999] {
        let moov = build_stbl_tree(invalid_ct, 10000, 0, b"CRAW");

        let mut builder = BoxBuilder::new();
        builder.add_box(b"ftyp", b"crx \0\0\0\x01crx \0\0\0\0");
        builder.buf.extend_from_slice(&moov);
        builder.add_box(b"PRVW", &prvw_payload);
        builder.add_box(b"mdat", &jpeg);

        let mut stream = ByteSliceStream::new(&builder.buf);
        let data = extract_cr3_preview_data(&mut stream)
            .expect("Invalid content type must fallback to PRVW");
        assert_eq!(data.preview.width, 1920);
        assert_eq!(data.preview.height, 1080);
    }
}

#[test]
fn test_stsz_and_co64_entry_count_validation() {
    let jpeg = create_valid_jpeg_payload();
    let prvw_payload = create_prvw_payload(1920, 1080, &jpeg);

    // STSZ sample_size = 0, count = 0
    let mut stsz_zero = vec![0u8; 12];
    stsz_zero[4..8].copy_from_slice(&0u32.to_be_bytes());
    stsz_zero[8..12].copy_from_slice(&0u32.to_be_bytes());

    let mut craw = vec![0u8; 92];
    craw[24..26].copy_from_slice(&3840u16.to_be_bytes());
    craw[26..28].copy_from_slice(&2160u16.to_be_bytes());
    craw[78..80].copy_from_slice(&3u16.to_be_bytes());

    let mut craw_box = BoxBuilder::new();
    craw_box.add_box(b"CRAW", &craw);

    let mut stsd_box = BoxBuilder::new();
    let mut stsd_p = vec![0u8; 8];
    stsd_p.extend_from_slice(&craw_box.buf);
    stsd_box.add_box(b"stsd", &stsd_p);

    let mut stsz_box = BoxBuilder::new();
    stsz_box.add_box(b"stsz", &stsz_zero);

    let mut co64_box = BoxBuilder::new();
    let mut co64_p = vec![0u8; 12];
    co64_p[4..8].copy_from_slice(&1u32.to_be_bytes());
    co64_box.add_box(b"co64", &co64_p);

    let mut stbl_p = Vec::new();
    stbl_p.extend_from_slice(&stsd_box.buf);
    stbl_p.extend_from_slice(&stsz_box.buf);
    stbl_p.extend_from_slice(&co64_box.buf);

    let mut stbl_box = BoxBuilder::new();
    stbl_box.add_box(b"stbl", &stbl_p);
    let mut minf = BoxBuilder::new();
    minf.add_box(b"minf", &stbl_box.buf);
    let mut mdia = BoxBuilder::new();
    mdia.add_box(b"mdia", &minf.buf);
    let mut trak = BoxBuilder::new();
    trak.add_box(b"trak", &mdia.buf);
    let mut moov = BoxBuilder::new();
    moov.add_box(b"moov", &trak.buf);

    let mut builder = BoxBuilder::new();
    builder.add_box(b"ftyp", b"crx \0\0\0\x01crx \0\0\0\0");
    builder.buf.extend_from_slice(&moov.buf);
    builder.add_box(b"PRVW", &prvw_payload);
    builder.add_box(b"mdat", &jpeg);

    let mut stream = ByteSliceStream::new(&builder.buf);
    let res = extract_cr3_preview_data(&mut stream);
    assert!(res.is_err(), "stsz count == 0 when sample_size == 0 must fail");
}

// ============================================================================
// 5. Missing MDAT Streams & PRVW Selection (DIV-04)
// ============================================================================

#[test]
fn test_div04_missing_mdat_fails_even_with_prvw() {
    let jpeg = create_valid_jpeg_payload();
    let prvw_payload = create_prvw_payload(1920, 1080, &jpeg);

    let mut builder = BoxBuilder::new();
    builder.add_box(b"ftyp", b"crx \0\0\0\x01crx \0\0\0\0");
    builder.add_box(b"moov", &[0u8; 8]);
    builder.add_box(b"PRVW", &prvw_payload);
    // Note: No mdat box

    let mut stream = ByteSliceStream::new(&builder.buf);
    let res = extract_cr3_preview_data(&mut stream);
    assert!(res.is_err(), "Missing mdat must fail without speculative retry");
}

#[test]
fn test_mdat_corrupt_jpeg_falls_back_to_prvw() {
    let jpeg = create_valid_jpeg_payload();
    let prvw_payload = create_prvw_payload(1920, 1080, &jpeg);
    let moov = build_stbl_tree(3, 10000, 0, b"CRAW");

    let mut builder = BoxBuilder::new();
    builder.add_box(b"ftyp", b"crx \0\0\0\x01crx \0\0\0\0");
    builder.buf.extend_from_slice(&moov);
    builder.add_box(b"PRVW", &prvw_payload);
    // Corrupt mdat (no JPEG SOI)
    builder.add_box(b"mdat", &[0x00, 0x00, 0x00, 0x00, 0x00, 0x00]);

    let mut stream = ByteSliceStream::new(&builder.buf);
    let data =
        extract_cr3_preview_data(&mut stream).expect("Corrupt mdat JPEG must fallback to PRVW");
    assert_eq!(data.preview.width, 1920);
    assert_eq!(data.preview.height, 1080);
    assert_eq!(data.preview.format, ImageFormat::JpegCompressed);
}

#[test]
fn test_mdat_hevc_selected_without_soi_check() {
    let jpeg = create_valid_jpeg_payload();
    let prvw_payload = create_prvw_payload(1920, 1080, &jpeg);
    let moov = build_stbl_tree(4 /* HEVC */, 10000, 0, b"CRAW");

    let mut builder = BoxBuilder::new();
    builder.add_box(b"ftyp", b"crx \0\0\0\x01crx \0\0\0\0");
    builder.buf.extend_from_slice(&moov);
    builder.add_box(b"PRVW", &prvw_payload);
    // HEVC NAL data in mdat
    builder.add_box(b"mdat", &[0x00, 0x00, 0x00, 0x01, 0x40, 0x01]);

    let mut stream = ByteSliceStream::new(&builder.buf);
    let data = extract_cr3_preview_data(&mut stream).expect("HEVC mdat must be selected");
    assert_eq!(data.preview.width, 3840);
    assert_eq!(data.preview.height, 2160);
    assert_eq!(data.preview.format, ImageFormat::HevcCompressed);
}

// ============================================================================
// 6. Systematic Fuzz Sweeps (Truncation & Bit-flips)
// ============================================================================

#[test]
fn test_systematic_truncation_sweep() {
    let jpeg = create_valid_jpeg_payload();
    let prvw_payload = create_prvw_payload(1920, 1080, &jpeg);
    let thmb_payload = create_thmb_payload(160, 120, &jpeg);
    let cmt1 = create_cmt1_payload();
    let cmt2 = create_cmt2_payload();
    let moov = build_stbl_tree(3, jpeg.len() as u32, 0, b"CRAW");

    let mut builder = BoxBuilder::new();
    builder.add_box(b"ftyp", b"crx \0\0\0\x01crx \0\0\0\0");
    builder.add_box(b"CMT1", &cmt1);
    builder.add_box(b"CMT2", &cmt2);
    builder.buf.extend_from_slice(&moov);
    builder.add_box(b"PRVW", &prvw_payload);
    builder.add_box(b"THMB", &thmb_payload);
    builder.add_box(b"mdat", &jpeg);

    let full_stream = builder.buf;

    // Truncate at every 5th byte from 0 to full_stream.len() - 1
    for cut in (0..full_stream.len() - 1).step_by(5) {
        let mut stream = ByteSliceStream::new(&full_stream[..cut]);
        let _ = extract_cr3_preview_data(&mut stream); // Must not panic
        let _ = extract_cr3_orientation(&mut stream); // Must not panic
    }
}

#[test]
fn test_regression_crash_03_cr3_truncated_trailer_box() {
    // Truncated ISOBMFF stream after ftyp box
    let truncated_cr3 = b"\x00\x00\x00\x18ftypcrx \x00\x00\x00\x00crx kDNG\x00\x00\x00\x08moov";
    let stream = ByteSliceStream::new(truncated_cr3);

    let res = extract_cr3_preview_data(&stream);
    // Verified graceful failure / error handling without panicking on truncated trailing boxes
    assert!(res.is_err());
}

// ============================================================================
// 7. Declared Box-Size Enforcement on Canon UUID Containers
// ============================================================================

fn append_uuid_box_with_declared_payload(
    buf: &mut Vec<u8>,
    extended_64: bool,
    declared_payload_len: usize,
    physical_payload: &[u8],
) {
    if extended_64 {
        buf.extend_from_slice(&1u32.to_be_bytes());
        buf.extend_from_slice(b"uuid");
        let total = 16u64 + declared_payload_len as u64;
        buf.extend_from_slice(&total.to_be_bytes());
    } else {
        let total = u32::try_from(8 + declared_payload_len).expect("compact box fits u32");
        buf.extend_from_slice(&total.to_be_bytes());
        buf.extend_from_slice(b"uuid");
    }
    buf.extend_from_slice(physical_payload);
}

#[test]
fn test_uuid_prvw_rejects_short_declared_payload_even_with_valid_trailing_boxes() {
    let jpeg = create_valid_jpeg_payload();
    let prvw_payload = create_prvw_payload(1280, 720, &jpeg);

    let mut prvw_box = BoxBuilder::new();
    prvw_box.add_box(b"PRVW", &prvw_payload);

    let mut mdat_box = BoxBuilder::new();
    mdat_box.add_box(b"mdat", &jpeg);

    let mut full_tail = Vec::new();
    full_tail.extend_from_slice(&UUID_PRVW);
    full_tail.extend_from_slice(&[0u8; 8]);
    full_tail.extend_from_slice(&prvw_box.buf);
    full_tail.extend_from_slice(&mdat_box.buf);

    for extended_64 in [false, true] {
        for declared_payload_len in 0..24 {
            let mut builder = BoxBuilder::new();
            builder.add_box(b"ftyp", b"crx \0\0\0\x01crx \0\0\0\0");
            let uuid_offset = builder.buf.len() as u64;
            append_uuid_box_with_declared_payload(
                &mut builder.buf,
                extended_64,
                declared_payload_len,
                &full_tail,
            );

            let stream = ByteSliceStream::new(&builder.buf);
            let res = extract_cr3_preview_data(&stream);
            assert_eq!(
                res,
                Err(PiexError::InvalidBoxHeader { offset: uuid_offset }),
                "UUID_PRVW with extended_64={extended_64} declared_payload_len={declared_payload_len} must fail",
            );
        }

        for declared_payload_len in [24, 24 + prvw_box.buf.len()] {
            let mut builder = BoxBuilder::new();
            builder.add_box(b"ftyp", b"crx \0\0\0\x01crx \0\0\0\0");
            append_uuid_box_with_declared_payload(
                &mut builder.buf,
                extended_64,
                declared_payload_len,
                &full_tail,
            );

            let stream = ByteSliceStream::new(&builder.buf);
            let data = extract_cr3_preview_data(&stream).unwrap_or_else(|err| {
                panic!(
                    "UUID_PRVW with extended_64={extended_64} declared_payload_len={declared_payload_len} failed: {err:?}"
                )
            });
            assert_eq!(data.preview.width, 1280);
            assert_eq!(data.preview.height, 720);
            assert_eq!(data.preview.format, ImageFormat::JpegCompressed);
        }
    }
}

#[test]
fn test_uuid_moov_rejects_short_declared_payload_for_preview_and_orientation() {
    let jpeg = create_valid_jpeg_payload();
    let cmt1_payload = create_cmt1_payload();
    let prvw_payload = create_prvw_payload(1280, 720, &jpeg);

    let mut children = BoxBuilder::new();
    children.add_box(b"CMT1", &cmt1_payload);
    children.add_box(b"PRVW", &prvw_payload);

    let mut mdat_box = BoxBuilder::new();
    mdat_box.add_box(b"mdat", &jpeg);

    let mut full_tail = Vec::new();
    full_tail.extend_from_slice(&UUID_MOOV);
    full_tail.extend_from_slice(&children.buf);
    full_tail.extend_from_slice(&mdat_box.buf);

    for extended_64 in [false, true] {
        for declared_payload_len in 0..16 {
            let mut builder = BoxBuilder::new();
            builder.add_box(b"ftyp", b"crx \0\0\0\x01crx \0\0\0\0");
            let uuid_offset = builder.buf.len() as u64;
            append_uuid_box_with_declared_payload(
                &mut builder.buf,
                extended_64,
                declared_payload_len,
                &full_tail,
            );

            let stream = ByteSliceStream::new(&builder.buf);
            assert_eq!(
                extract_cr3_preview_data(&stream),
                Err(PiexError::InvalidBoxHeader { offset: uuid_offset }),
                "UUID_MOOV preview with extended_64={extended_64} declared_payload_len={declared_payload_len}",
            );
            assert_eq!(
                extract_cr3_orientation(&stream),
                Err(PiexError::InvalidBoxHeader { offset: uuid_offset }),
                "UUID_MOOV orientation with extended_64={extended_64} declared_payload_len={declared_payload_len}",
            );
        }

        for declared_payload_len in [16, 16 + children.buf.len()] {
            let mut builder = BoxBuilder::new();
            builder.add_box(b"ftyp", b"crx \0\0\0\x01crx \0\0\0\0");
            append_uuid_box_with_declared_payload(
                &mut builder.buf,
                extended_64,
                declared_payload_len,
                &full_tail,
            );

            let stream = ByteSliceStream::new(&builder.buf);
            let data = extract_cr3_preview_data(&stream).unwrap_or_else(|err| {
                panic!(
                    "UUID_MOOV preview with extended_64={extended_64} declared_payload_len={declared_payload_len} failed: {err:?}"
                )
            });
            assert_eq!(data.preview.width, 1280);
            assert_eq!(data.preview.height, 720);
            assert_eq!(data.exif_orientation, 6);
            assert_eq!(extract_cr3_orientation(&stream), Ok(6));
        }
    }
}

#[test]
fn test_extract_cr3_preview_data_into_preserves_prior_boxes_on_short_uuid_error() {
    let jpeg = create_valid_jpeg_payload();
    let cmt1_payload = create_cmt1_payload();
    let thmb_payload = create_thmb_payload(160, 120, &jpeg);
    let prvw_payload = create_prvw_payload(1280, 720, &jpeg);

    let mut prvw_and_mdat = BoxBuilder::new();
    prvw_and_mdat.add_box(b"PRVW", &prvw_payload);
    prvw_and_mdat.add_box(b"mdat", &jpeg);

    let mut prvw_tail = Vec::new();
    prvw_tail.extend_from_slice(&UUID_PRVW);
    prvw_tail.extend_from_slice(&[0u8; 8]);
    prvw_tail.extend_from_slice(&prvw_and_mdat.buf);

    let mut builder = BoxBuilder::new();
    builder.add_box(b"ftyp", b"crx \0\0\0\x01crx \0\0\0\0");
    builder.add_box(b"CMT1", &cmt1_payload);
    builder.add_box(b"THMB", &thmb_payload);
    let short_uuid_offset = builder.buf.len() as u64;
    append_uuid_box_with_declared_payload(&mut builder.buf, false, 23, &prvw_tail);

    let stream = ByteSliceStream::new(&builder.buf);
    let mut out =
        PreviewImageData { iso: 3200, date_time: b"keep\0".to_vec(), ..Default::default() };

    let res = extract_cr3_preview_data_into(&stream, &mut out);
    assert_eq!(res, Err(PiexError::InvalidBoxHeader { offset: short_uuid_offset }));

    assert_eq!(out.maker, b"Canon\0");
    assert_eq!(out.model, b"EOS R5\0");
    assert_eq!(out.exif_orientation, 6);
    assert_eq!(out.full_width, 8192);
    assert_eq!(out.full_height, 5464);
    assert_eq!(out.thumbnail.width, 160);
    assert_eq!(out.thumbnail.height, 120);
    assert_eq!(out.iso, 3200);
    assert_eq!(out.date_time, b"keep\0");
    assert_eq!(out.preview.length, 0);
    assert_eq!(out.preview.width, 0);
}

#[test]
fn test_short_uuid_moov_preserves_prior_box_writes_across_preview_and_orientation() {
    let jpeg = create_valid_jpeg_payload();
    let cmt1_payload = create_cmt1_payload();
    let cmt2_payload = create_cmt2_payload();
    let thmb_payload = create_thmb_payload(160, 120, &jpeg);
    let prvw_payload = create_prvw_payload(1280, 720, &jpeg);

    // Subcase A: CMT2 + THMB precede short UUID_MOOV (declared payload = 15), with CMT1 + PRVW +
    // mdat trailing. Both preview and orientation fail at UUID_MOOV; prior CMT2 + THMB writes and
    // untouched CMT1 caller fields remain intact.
    let mut post_moov_a = BoxBuilder::new();
    post_moov_a.add_box(b"CMT1", &cmt1_payload);
    post_moov_a.add_box(b"PRVW", &prvw_payload);
    post_moov_a.add_box(b"mdat", &jpeg);
    let mut tail_a = Vec::new();
    tail_a.extend_from_slice(&UUID_MOOV);
    tail_a.extend_from_slice(&post_moov_a.buf);

    for extended_64 in [false, true] {
        let mut builder = BoxBuilder::new();
        builder.add_box(b"ftyp", b"crx \0\0\0\x01crx \0\0\0\0");
        builder.add_box(b"CMT2", &cmt2_payload);
        builder.add_box(b"THMB", &thmb_payload);
        let uuid_offset = builder.buf.len() as u64;
        append_uuid_box_with_declared_payload(&mut builder.buf, extended_64, 15, &tail_a);

        let stream = ByteSliceStream::new(&builder.buf);
        let mut out = PreviewImageData {
            maker: b"keep_make\0".to_vec(),
            model: b"keep_model\0".to_vec(),
            exif_orientation: 8,
            full_width: 777,
            full_height: 888,
            ..Default::default()
        };
        assert_eq!(
            extract_cr3_preview_data_into(&stream, &mut out),
            Err(PiexError::InvalidBoxHeader { offset: uuid_offset }),
        );
        assert_eq!(out.iso, 400);
        assert_eq!(out.date_time, b"2024:01:15 10:20:30\0");
        assert_eq!(out.exposure_time, Rational::new(1, 500));
        assert_eq!(out.fnumber, Rational::new(40, 10));
        assert_eq!(out.focal_length, Rational::new(85, 1));
        assert_eq!(out.thumbnail.width, 160);
        assert_eq!(out.thumbnail.height, 120);
        assert_eq!(out.maker, b"keep_make\0");
        assert_eq!(out.model, b"keep_model\0");
        assert_eq!(out.exif_orientation, 8);
        assert_eq!(out.full_width, 777);
        assert_eq!(out.full_height, 888);
        assert_eq!(out.preview.length, 0);
        assert_eq!(
            extract_cr3_orientation(&stream),
            Err(PiexError::InvalidBoxHeader { offset: uuid_offset }),
        );
    }

    // Subcase B: CMT1 + THMB precede short UUID_MOOV (declared payload = 15), with CMT2 + PRVW +
    // mdat trailing. Preview fails preserving CMT1 + THMB and caller CMT2 sentinels, whereas
    // orientation stops at CMT1 and returns Ok(6) before reaching UUID_MOOV.
    let mut post_moov_b = BoxBuilder::new();
    post_moov_b.add_box(b"CMT2", &cmt2_payload);
    post_moov_b.add_box(b"PRVW", &prvw_payload);
    post_moov_b.add_box(b"mdat", &jpeg);
    let mut tail_b = Vec::new();
    tail_b.extend_from_slice(&UUID_MOOV);
    tail_b.extend_from_slice(&post_moov_b.buf);

    for extended_64 in [false, true] {
        let mut builder = BoxBuilder::new();
        builder.add_box(b"ftyp", b"crx \0\0\0\x01crx \0\0\0\0");
        builder.add_box(b"CMT1", &cmt1_payload);
        builder.add_box(b"THMB", &thmb_payload);
        let uuid_offset = builder.buf.len() as u64;
        append_uuid_box_with_declared_payload(&mut builder.buf, extended_64, 15, &tail_b);

        let stream = ByteSliceStream::new(&builder.buf);
        let mut out =
            PreviewImageData { iso: 3200, date_time: b"keep\0".to_vec(), ..Default::default() };
        assert_eq!(
            extract_cr3_preview_data_into(&stream, &mut out),
            Err(PiexError::InvalidBoxHeader { offset: uuid_offset }),
        );
        assert_eq!(out.maker, b"Canon\0");
        assert_eq!(out.model, b"EOS R5\0");
        assert_eq!(out.exif_orientation, 6);
        assert_eq!(out.full_width, 8192);
        assert_eq!(out.full_height, 5464);
        assert_eq!(out.thumbnail.width, 160);
        assert_eq!(out.iso, 3200);
        assert_eq!(out.date_time, b"keep\0");
        assert_eq!(out.preview.length, 0);
        assert_eq!(extract_cr3_orientation(&stream), Ok(6));
    }
}

#[test]
fn test_unknown_uuid_short_declared_payload_skips_to_next_box_offset() {
    let jpeg = create_valid_jpeg_payload();
    let prvw_payload = create_prvw_payload(1280, 720, &jpeg);

    let mut builder = BoxBuilder::new();
    builder.add_box(b"ftyp", b"crx \0\0\0\x01crx \0\0\0\0");
    append_uuid_box_with_declared_payload(&mut builder.buf, false, 8, &[0xAAu8; 8]);
    let expected_prvw_jpeg_offset = u32::try_from(builder.buf.len() + 24).expect("fits u32");
    builder.add_box(b"PRVW", &prvw_payload);
    builder.add_box(b"mdat", &jpeg);

    let stream = ByteSliceStream::new(&builder.buf);
    let data =
        extract_cr3_preview_data(&stream).expect("short unknown UUID with trailing box skips");
    assert_eq!(data.preview.width, 1280);
    assert_eq!(data.preview.height, 720);
    assert_eq!(data.preview.offset, expected_prvw_jpeg_offset);
}
