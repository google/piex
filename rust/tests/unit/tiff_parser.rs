use piex::stream::{ByteSliceStream, Endian};
use piex::tiff::directory::{TiffContent, TiffDirectory, TiffError, TiffType};
use piex::tiff::parser::TiffParser;
use piex::tiff::tags::gps::*;
use piex::tiff::tags::panasonic::*;
use piex::tiff::tags::tiff::*;
use piex::types::{ColorSpace, Image, ImageFormat, PreviewImageData, Rational};

// ============================================================================
// Binary Helper Utility
// ============================================================================

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

// ============================================================================
// 1. TiffParser::get_endianness & Constructors
// ============================================================================

#[test]
fn test_get_endianness_valid_little_endian() {
    let mut header = Vec::new();
    header.extend_from_slice(b"II");
    header.extend_from_slice(&42u16.to_le_bytes());
    header.extend_from_slice(&8u32.to_le_bytes());

    let stream = ByteSliceStream::new(&header);
    let endian = TiffParser::get_endianness(&stream, 0).unwrap();
    assert_eq!(endian, Endian::Little);

    let parser = TiffParser::new(&stream).unwrap();
    assert_eq!(parser.endian(), Endian::Little);
    assert_eq!(parser.tiff_offset(), 0);
}

#[test]
fn test_get_endianness_valid_big_endian() {
    let mut header = Vec::new();
    header.extend_from_slice(b"MM");
    header.extend_from_slice(&42u16.to_be_bytes());
    header.extend_from_slice(&8u32.to_be_bytes());

    let stream = ByteSliceStream::new(&header);
    let endian = TiffParser::get_endianness(&stream, 0).unwrap();
    assert_eq!(endian, Endian::Big);

    let parser = TiffParser::new(&stream).unwrap();
    assert_eq!(parser.endian(), Endian::Big);
    assert_eq!(parser.tiff_offset(), 0);
}

#[test]
fn test_get_endianness_invalid_byte_order() {
    let header = b"XX\x2A\x00\x08\x00\x00\x00";
    let stream = ByteSliceStream::new(header);
    let res = TiffParser::get_endianness(&stream, 0);
    assert_eq!(res, Err(TiffError::InvalidByteOrder));
}

#[test]
fn test_get_endianness_accepts_arbitrary_magic() {
    // LE prefix with magic 43 (0x002B)
    let header_le = b"II\x2B\x00\x08\x00\x00\x00";
    let stream_le = ByteSliceStream::new(header_le);
    assert_eq!(TiffParser::get_endianness(&stream_le, 0), Ok(Endian::Little));

    // BE prefix with magic 0x1234
    let header_be = b"MM\x12\x34\x00\x00\x00\x08";
    let stream_be = ByteSliceStream::new(header_be);
    assert_eq!(TiffParser::get_endianness(&stream_be, 0), Ok(Endian::Big));
}

#[test]
fn test_get_endianness_truncated_stream() {
    let header = b"II\x2A"; // Only 3 bytes
    let stream = ByteSliceStream::new(header);
    let res = TiffParser::get_endianness(&stream, 0);
    assert!(matches!(res, Err(TiffError::StreamError(_))));
}

#[test]
fn test_constructor_with_offset() {
    let mut buf = vec![0u8; 16]; // 16 bytes padding
    let tiff_offset = 16u32;
    buf.extend_from_slice(b"II");
    buf.extend_from_slice(&42u16.to_le_bytes());
    buf.extend_from_slice(&8u32.to_le_bytes());

    let stream = ByteSliceStream::new(&buf);
    let parser = TiffParser::with_offset(&stream, tiff_offset).unwrap();
    assert_eq!(parser.endian(), Endian::Little);
    assert_eq!(parser.tiff_offset(), 16);
}

// ============================================================================
// 2. parse_directory Tests
// ============================================================================

#[test]
fn test_parse_directory_inline_and_offset_entries() {
    for endian in [Endian::Little, Endian::Big] {
        let mut writer = TiffWriter::new(endian);
        writer.write_header(8);

        // IFD starts at byte 8
        // Number of entries = 2
        writer.buf.extend_from_slice(&writer.pack_u16(2));

        // Entry 1: Inline Short tag (width = 1920) -> tag 0x0100, type 3, count 1
        let width_bytes = writer.pack_u16(1920);
        writer.write_entry_inline(TIFF_TAG_IMAGE_WIDTH as u16, TiffType::Short, 1, &width_bytes);

        // Entry 2: Out-of-line ASCII tag (Make = "Canon\0") -> tag 0x010F, type 2, count 6
        // Data payload offset at byte 100
        let payload_offset = 100u32;
        writer.write_entry_offset(TIFF_TAG_MAKE as u16, TiffType::Ascii, 6, payload_offset);

        // Next IFD offset = 0
        writer.buf.extend_from_slice(&writer.pack_u32(0));

        // Pad up to payload_offset
        if writer.buf.len() < payload_offset as usize {
            writer.buf.resize(payload_offset as usize, 0);
        }
        writer.buf.extend_from_slice(b"Canon\0");

        let stream = ByteSliceStream::new(&writer.buf);
        let parser = TiffParser::new(&stream).unwrap();

        let desired_tags = [TIFF_TAG_IMAGE_WIDTH, TIFF_TAG_MAKE];

        let mut next_ifd_offset = 999u32;
        let dir = parser.parse_directory(8, &desired_tags, &mut next_ifd_offset).unwrap();

        assert_eq!(next_ifd_offset, 0);
        assert_eq!(dir.get_u32(TIFF_TAG_IMAGE_WIDTH).unwrap(), 1920);
        assert_eq!(dir.get_string(TIFF_TAG_MAKE).unwrap(), "Canon");
    }
}

#[test]
fn test_parse_directory_desired_tags_filtering() {
    let mut writer = TiffWriter::new(Endian::Little);
    writer.write_header(8);
    writer.buf.extend_from_slice(&writer.pack_u16(2));

    // Entry 1: tag 0x0100 (IMAGE_WIDTH)
    writer.write_entry_inline(
        TIFF_TAG_IMAGE_WIDTH as u16,
        TiffType::Short,
        1,
        &100u16.to_le_bytes(),
    );
    // Entry 2: tag 0x0101 (IMAGE_LENGTH)
    writer.write_entry_inline(
        TIFF_TAG_IMAGE_LENGTH as u16,
        TiffType::Short,
        1,
        &200u16.to_le_bytes(),
    );
    writer.buf.extend_from_slice(&0u32.to_le_bytes());

    let stream = ByteSliceStream::new(&writer.buf);
    let parser = TiffParser::new(&stream).unwrap();

    // Filter only tag 0x0100
    let desired = [TIFF_TAG_IMAGE_WIDTH];

    let mut next_ifd = 0;
    let dir = parser.parse_directory(8, &desired, &mut next_ifd).unwrap();

    assert!(dir.has(TIFF_TAG_IMAGE_WIDTH));
    assert!(!dir.has(TIFF_TAG_IMAGE_LENGTH));
}

#[test]
fn test_parse_directory_corrupted_entry_handling() {
    // 1. Invalid type ID
    let mut writer = TiffWriter::new(Endian::Little);
    writer.write_header(8);
    writer.buf.extend_from_slice(&1u16.to_le_bytes()); // 1 entry
    writer.buf.extend_from_slice(&0x0100u16.to_le_bytes()); // tag
    writer.buf.extend_from_slice(&99u16.to_le_bytes()); // invalid type ID 99
    writer.buf.extend_from_slice(&1u32.to_le_bytes()); // count
    writer.buf.extend_from_slice(&0u32.to_le_bytes()); // value
    writer.buf.extend_from_slice(&0u32.to_le_bytes()); // next offset

    let stream = ByteSliceStream::new(&writer.buf);
    let parser = TiffParser::new(&stream).unwrap();
    let desired = [0x0100];

    let mut next_ifd = 0;
    let dir = parser.parse_directory(8, &desired, &mut next_ifd).unwrap();
    assert!(!dir.has(0x0100));

    // 2. Element size count multiplication overflow
    let mut writer_overflow = TiffWriter::new(Endian::Little);
    writer_overflow.write_header(8);
    writer_overflow.buf.extend_from_slice(&1u16.to_le_bytes());
    writer_overflow.buf.extend_from_slice(&0x0100u16.to_le_bytes());
    writer_overflow.buf.extend_from_slice(&(TiffType::Long as u16).to_le_bytes()); // 4 bytes/elem
    writer_overflow.buf.extend_from_slice(&(u32::MAX / 2).to_le_bytes()); // Overflow count
    writer_overflow.buf.extend_from_slice(&0u32.to_le_bytes());
    writer_overflow.buf.extend_from_slice(&0u32.to_le_bytes());

    let stream_overflow = ByteSliceStream::new(&writer_overflow.buf);
    let parser_overflow = TiffParser::new(&stream_overflow).unwrap();
    let res_overflow = parser_overflow.parse_directory(8, &desired, &mut next_ifd);
    assert!(matches!(res_overflow, Err(TiffError::CorruptedEntry { .. })));
}

// ============================================================================
// 3. parse & parse_ifd_chain Multi-IFD / Cycle Detection / Limit Tests
// ============================================================================

#[test]
fn test_parse_multi_ifd_chain_traversal() {
    let mut writer = TiffWriter::new(Endian::Little);
    writer.write_header(8);

    // IFD 0 at offset 8 -> points to IFD 1 at offset 32
    writer.buf.extend_from_slice(&1u16.to_le_bytes());
    writer.write_entry_inline(
        TIFF_TAG_IMAGE_WIDTH as u16,
        TiffType::Short,
        1,
        &100u16.to_le_bytes(),
    );
    writer.buf.extend_from_slice(&32u32.to_le_bytes()); // Next IFD offset = 32

    // Pad to offset 32
    writer.buf.resize(32, 0);

    // IFD 1 at offset 32 -> next IFD = 0
    writer.buf.extend_from_slice(&1u16.to_le_bytes());
    writer.write_entry_inline(
        TIFF_TAG_IMAGE_LENGTH as u16,
        TiffType::Short,
        1,
        &200u16.to_le_bytes(),
    );
    writer.buf.extend_from_slice(&0u32.to_le_bytes()); // End of chain

    let stream = ByteSliceStream::new(&writer.buf);
    let parser = TiffParser::new(&stream).unwrap();

    let desired = [TIFF_TAG_IMAGE_WIDTH, TIFF_TAG_IMAGE_LENGTH];

    let content = parser.parse(&desired, 10).unwrap();
    assert_eq!(content.tiff_directories.len(), 2);
    assert_eq!(content.tiff_directories[0].get_u32(TIFF_TAG_IMAGE_WIDTH).unwrap(), 100);
    assert_eq!(content.tiff_directories[1].get_u32(TIFF_TAG_IMAGE_LENGTH).unwrap(), 200);
}

#[test]
fn test_parse_ifd_chain_cycle_detection_guard() {
    let mut writer = TiffWriter::new(Endian::Little);
    writer.write_header(8);

    // IFD 0 at offset 8 -> points to IFD 1 at offset 32
    writer.buf.extend_from_slice(&1u16.to_le_bytes());
    writer.write_entry_inline(
        TIFF_TAG_IMAGE_WIDTH as u16,
        TiffType::Short,
        1,
        &100u16.to_le_bytes(),
    );
    writer.buf.extend_from_slice(&32u32.to_le_bytes());

    // IFD 1 at offset 32 -> points BACK to IFD 0 at offset 8 (Cycle!)
    writer.buf.resize(32, 0);
    writer.buf.extend_from_slice(&1u16.to_le_bytes());
    writer.write_entry_inline(
        TIFF_TAG_IMAGE_LENGTH as u16,
        TiffType::Short,
        1,
        &200u16.to_le_bytes(),
    );
    writer.buf.extend_from_slice(&8u32.to_le_bytes()); // Cycle back to offset 8!

    let stream = ByteSliceStream::new(&writer.buf);
    let parser = TiffParser::new(&stream).unwrap();

    let desired = [TIFF_TAG_IMAGE_WIDTH, TIFF_TAG_IMAGE_LENGTH];

    // Should terminate gracefully without infinite loop and collect unique visited IFDs
    let content = parser.parse(&desired, 10).unwrap();
    assert_eq!(content.tiff_directories.len(), 2);
}

#[test]
fn test_parse_max_number_ifds_limit() {
    let mut writer = TiffWriter::new(Endian::Little);
    writer.write_header(8);

    // IFD 0 at offset 8 -> points to IFD 1 at offset 32
    writer.buf.extend_from_slice(&1u16.to_le_bytes());
    writer.write_entry_inline(
        TIFF_TAG_IMAGE_WIDTH as u16,
        TiffType::Short,
        1,
        &100u16.to_le_bytes(),
    );
    writer.buf.extend_from_slice(&32u32.to_le_bytes());

    // IFD 1 at offset 32
    writer.buf.resize(32, 0);
    writer.buf.extend_from_slice(&1u16.to_le_bytes());
    writer.write_entry_inline(
        TIFF_TAG_IMAGE_LENGTH as u16,
        TiffType::Short,
        1,
        &200u16.to_le_bytes(),
    );
    writer.buf.extend_from_slice(&0u32.to_le_bytes());

    let stream = ByteSliceStream::new(&writer.buf);
    let parser = TiffParser::new(&stream).unwrap();

    let desired = [TIFF_TAG_IMAGE_WIDTH, TIFF_TAG_IMAGE_LENGTH];

    // Set max_number_ifds = 1 -> should stop after 1st IFD
    let res = parser.parse(&desired, 1).unwrap();
    assert_eq!(res.tiff_directories.len(), 1);
}

// ============================================================================
// 4. get_jpeg_dimensions Tests
// ============================================================================

fn create_jpeg_stream(sof_marker: u16, width: u16, height: u16) -> Vec<u8> {
    let mut buf = Vec::new();
    // Prepend 8-byte valid TIFF header so TiffParser constructor succeeds
    buf.extend_from_slice(b"II\x2A\x00\x08\x00\x00\x00");

    // JPEG SOI at offset 8
    buf.extend_from_slice(&[0xFF, 0xD8]);

    // Extra APP0 marker
    buf.extend_from_slice(&[0xFF, 0xE0]);
    buf.extend_from_slice(&16u16.to_be_bytes()); // length = 16
    buf.extend_from_slice(&[0u8; 14]); // payload

    // SOF marker
    buf.extend_from_slice(&sof_marker.to_be_bytes());
    buf.extend_from_slice(&17u16.to_be_bytes()); // SOF segment length = 17
    buf.push(8); // Precision = 8
    buf.extend_from_slice(&height.to_be_bytes());
    buf.extend_from_slice(&width.to_be_bytes());
    buf.push(3); // 3 components
    buf.extend_from_slice(&[0u8; 9]);

    // SOS marker
    buf.extend_from_slice(&[0xFF, 0xDA]);

    buf
}

#[test]
fn test_get_jpeg_dimensions_accepts_sof0_only() {
    let jpeg_buf = create_jpeg_stream(0xFFC0, 1920, 1080);
    let stream = ByteSliceStream::new(&jpeg_buf);
    let parser = TiffParser::new(&stream).unwrap();
    let (w, h) = parser.get_jpeg_dimensions(8).unwrap();
    assert_eq!(w, 1920);
    assert_eq!(h, 1080);
}

#[test]
fn test_get_jpeg_dimensions_rejects_sof1_and_sof2() {
    let non_sof0_markers = [0xFFC1u16, 0xFFC2u16, 0xFFC3u16];
    for &sof in &non_sof0_markers {
        let jpeg_buf = create_jpeg_stream(sof, 1920, 1080);
        let stream = ByteSliceStream::new(&jpeg_buf);
        let parser = TiffParser::new(&stream).unwrap();
        let res = parser.get_jpeg_dimensions(8);
        assert_eq!(res, Err(TiffError::TagNotFound { tag: 0xFFC0 }));
    }
}

#[test]
fn test_get_jpeg_dimensions_corrupted_cases() {
    // 1. SOS before SOF
    let mut jpeg_sos_first = Vec::new();
    jpeg_sos_first.extend_from_slice(b"II\x2A\x00\x08\x00\x00\x00"); // TIFF header
    jpeg_sos_first.extend_from_slice(&[0xFF, 0xD8]); // SOI at offset 8
    jpeg_sos_first.extend_from_slice(&[0xFF, 0xDA]); // SOS
    let stream_sos = ByteSliceStream::new(&jpeg_sos_first);
    let parser_sos = TiffParser::new(&stream_sos).unwrap();
    assert_eq!(parser_sos.get_jpeg_dimensions(8), Err(TiffError::TagNotFound { tag: 0xFFC0 }));

    // 2. Marker with invalid length < 2
    let mut jpeg_bad_len = Vec::new();
    jpeg_bad_len.extend_from_slice(b"II\x2A\x00\x08\x00\x00\x00"); // TIFF header
    jpeg_bad_len.extend_from_slice(&[0xFF, 0xD8]); // SOI at offset 8
    jpeg_bad_len.extend_from_slice(&[0xFF, 0xE0]); // APP0
    jpeg_bad_len.extend_from_slice(&1u16.to_be_bytes()); // invalid length 1 < 2
    let stream_bad_len = ByteSliceStream::new(&jpeg_bad_len);
    let parser_bad_len = TiffParser::new(&stream_bad_len).unwrap();
    assert!(matches!(parser_bad_len.get_jpeg_dimensions(8), Err(TiffError::CorruptedEntry { .. })));
}

// ============================================================================
// 5. get_preview_image_data & Metadata Extraction Tests
// ============================================================================

#[test]
fn test_get_preview_image_data_full_metadata_extraction() {
    let mut main_dir = TiffDirectory::new(Endian::Little);
    main_dir.add_entry(TIFF_TAG_MAKE, TiffType::Ascii, 6, 0, b"Canon\0".to_vec());
    main_dir.add_entry(TIFF_TAG_MODEL, TiffType::Ascii, 9, 0, b"EOS 80D\0".to_vec());
    main_dir.add_entry(TIFF_TAG_ORIENTATION, TiffType::Short, 1, 0, 6u16.to_le_bytes().to_vec());
    main_dir.add_entry(EXIF_TAG_COLOR_SPACE, TiffType::Short, 1, 0, 1u16.to_le_bytes().to_vec()); // sRGB
    main_dir.add_entry(EXIF_TAG_WIDTH, TiffType::Long, 1, 0, 6000u32.to_le_bytes().to_vec());
    main_dir.add_entry(EXIF_TAG_HEIGHT, TiffType::Long, 1, 0, 4000u32.to_le_bytes().to_vec());
    main_dir.add_entry(EXIF_TAG_EXIF_IFD, TiffType::Long, 1, 0, 1000u32.to_le_bytes().to_vec());
    main_dir.add_entry(EXIF_TAG_GPS, TiffType::Long, 1, 0, 2000u32.to_le_bytes().to_vec());

    let mut exif_dir = TiffDirectory::new(Endian::Little);
    exif_dir.add_entry(
        EXIF_TAG_DATE_TIME_ORIGINAL,
        TiffType::Ascii,
        20,
        0,
        b"2023:10:25 14:30:00\0".to_vec(),
    );
    exif_dir.add_entry(EXIF_TAG_ISO_SPEED, TiffType::Long, 1, 0, 800u32.to_le_bytes().to_vec());

    // Exposure time: 1/250 -> rational
    let mut exp_bytes = Vec::new();
    exp_bytes.extend_from_slice(&1u32.to_le_bytes());
    exp_bytes.extend_from_slice(&250u32.to_le_bytes());
    exif_dir.add_entry(EXIF_TAG_EXPOSURE_TIME, TiffType::Rational, 1, 0, exp_bytes);

    // FNumber: 28/10 = 2.8
    let mut fn_bytes = Vec::new();
    fn_bytes.extend_from_slice(&28u32.to_le_bytes());
    fn_bytes.extend_from_slice(&10u32.to_le_bytes());
    exif_dir.add_entry(EXIF_TAG_FNUMBER, TiffType::Rational, 1, 0, fn_bytes);

    // Focal length: 50/1
    let mut fl_bytes = Vec::new();
    fl_bytes.extend_from_slice(&50u32.to_le_bytes());
    fl_bytes.extend_from_slice(&1u32.to_le_bytes());
    exif_dir.add_entry(EXIF_TAG_FOCAL_LENGTH, TiffType::Rational, 1, 0, fl_bytes);

    // GPS directory setup
    let mut gps_dir = TiffDirectory::new(Endian::Little);
    gps_dir.add_entry(GPS_TAG_LATITUDE_REF, TiffType::Ascii, 2, 0, b"N\0".to_vec());
    gps_dir.add_entry(GPS_TAG_LONGITUDE_REF, TiffType::Ascii, 2, 0, b"W\0".to_vec());

    // Lat: 37/1, 46/1, 29/1
    let mut lat_bytes = Vec::new();
    for &(num, den) in &[(37u32, 1u32), (46, 1), (29, 1)] {
        lat_bytes.extend_from_slice(&num.to_le_bytes());
        lat_bytes.extend_from_slice(&den.to_le_bytes());
    }
    gps_dir.add_entry(GPS_TAG_LATITUDE, TiffType::Rational, 3, 0, lat_bytes);

    // Lon: 122/1, 25/1, 10/1
    let mut lon_bytes = Vec::new();
    for &(num, den) in &[(122u32, 1u32), (25, 1), (10, 1)] {
        lon_bytes.extend_from_slice(&num.to_le_bytes());
        lon_bytes.extend_from_slice(&den.to_le_bytes());
    }
    gps_dir.add_entry(GPS_TAG_LONGITUDE, TiffType::Rational, 3, 0, lon_bytes);

    // Time stamp: 10/1, 15/1, 30/1
    let mut time_bytes = Vec::new();
    for &(num, den) in &[(10u32, 1u32), (15, 1), (30, 1)] {
        time_bytes.extend_from_slice(&num.to_le_bytes());
        time_bytes.extend_from_slice(&den.to_le_bytes());
    }
    gps_dir.add_entry(GPS_TAG_TIME_STAMP, TiffType::Rational, 3, 0, time_bytes);
    gps_dir.add_entry(GPS_TAG_DATE_STAMP, TiffType::Ascii, 11, 0, b"2023:10:25\0".to_vec());

    let mut content = TiffContent::new();
    content.tiff_directories.push(main_dir);
    content.exif_directory = Some(exif_dir);
    content.gps_directory = Some(gps_dir);

    // Create a dummy stream for TiffParser
    let mut header = Vec::new();
    header.extend_from_slice(b"II\x2A\x00\x08\x00\x00\x00");
    let stream = ByteSliceStream::new(&header);
    let parser = TiffParser::new(&stream).unwrap();

    let mut preview_data = PreviewImageData::default();
    parser.get_preview_image_data(&content, &mut preview_data).unwrap();

    assert_eq!(preview_data.maker, b"Canon\0");
    assert_eq!(preview_data.model, b"EOS 80D\0");
    assert_eq!(preview_data.exif_orientation, 6);
    assert_eq!(preview_data.color_space, ColorSpace::Srgb);
    assert_eq!(preview_data.full_width, 6000);
    assert_eq!(preview_data.full_height, 4000);
    assert_eq!(preview_data.date_time, b"2023:10:25 14:30:00\0");
    assert_eq!(preview_data.iso, 800);
    assert_eq!(preview_data.exposure_time, Rational::new(1, 250));
    assert_eq!(preview_data.fnumber, Rational::new(28, 10));
    assert_eq!(preview_data.focal_length, Rational::new(50, 1));

    // Check GPS
    assert!(preview_data.gps.is_valid);
    assert_eq!(preview_data.gps.latitude_ref, b'N');
    assert_eq!(preview_data.gps.longitude_ref, b'W');
    assert_eq!(preview_data.gps.latitude[0], Rational::new(37, 1));
    assert_eq!(preview_data.gps.longitude[0], Rational::new(122, 1));
    assert_eq!(preview_data.gps.date_stamp, b"2023:10:25");
}

#[test]
fn test_thumbnail_vs_preview_image_selection() {
    let mut header = Vec::new();
    header.extend_from_slice(b"II\x2A\x00\x08\x00\x00\x00");
    let stream = ByteSliceStream::new(&header);
    let parser = TiffParser::new(&stream).unwrap();

    // 1. Thumbnail image (width 320, height 240 <= 512)
    let thumb_image = Image::new(320, 240, 5000, 1000, ImageFormat::JpegCompressed);
    assert!(parser.is_thumbnail(&thumb_image));

    // 2. Preview image (width 1920, height 1080 > 512)
    let preview_image = Image::new(1920, 1080, 50000, 10000, ImageFormat::JpegCompressed);
    assert!(!parser.is_thumbnail(&preview_image));

    // Test thumbnail boundary exactly 512px
    let boundary_image = Image::new(512, 512, 10000, 2000, ImageFormat::JpegCompressed);
    assert!(parser.is_thumbnail(&boundary_image));

    // 513px is preview
    let over_boundary_image = Image::new(513, 512, 10000, 2000, ImageFormat::JpegCompressed);
    assert!(!parser.is_thumbnail(&over_boundary_image));
}

#[test]
fn test_panasonic_iso_tag_fallback() {
    let mut dir = TiffDirectory::new(Endian::Little);
    dir.add_entry(PANASONIC_TAG_ISO, TiffType::Short, 1, 0, 400u16.to_le_bytes().to_vec());

    let mut content = TiffContent::new();
    content.tiff_directories.push(dir);

    let mut header = Vec::new();
    header.extend_from_slice(b"II\x2A\x00\x08\x00\x00\x00");
    let stream = ByteSliceStream::new(&header);
    let parser = TiffParser::new(&stream).unwrap();

    let mut preview_data = PreviewImageData::default();
    assert_eq!(preview_data.iso, 0);

    parser.get_preview_image_data(&content, &mut preview_data).unwrap();
    assert_eq!(preview_data.iso, 400);
}

// ============================================================================
// 6. get_image_data & Strip Offsets
// ============================================================================

#[test]
fn test_get_image_data_jpeg_offset_and_byte_count() {
    let mut dir = TiffDirectory::new(Endian::Little);
    dir.add_entry(TIFF_TAG_JPEG_OFFSET, TiffType::Long, 1, 0, 2048u32.to_le_bytes().to_vec());
    dir.add_entry(TIFF_TAG_JPEG_BYTE_COUNT, TiffType::Long, 1, 0, 10240u32.to_le_bytes().to_vec());

    let mut header = Vec::new();
    header.extend_from_slice(b"II\x2A\x00\x08\x00\x00\x00");
    let stream = ByteSliceStream::new(&header);
    let parser = TiffParser::new(&stream).unwrap();

    let img = parser.get_image_data(&dir).unwrap();
    assert_eq!(img.offset, 2048);
    assert_eq!(img.length, 10240);
    assert_eq!(img.format, ImageFormat::JpegCompressed);
}

#[test]
fn test_get_image_data_consecutive_vs_non_consecutive_strips() {
    let mut header = Vec::new();
    header.extend_from_slice(b"II\x2A\x00\x08\x00\x00\x00");
    let stream = ByteSliceStream::new(&header);
    let parser = TiffParser::new(&stream).unwrap();

    // Consecutive strips: strip 0 at 100 (len 50), strip 1 at 150 (len 50)
    let mut dir_consecutive = TiffDirectory::new(Endian::Little);
    let mut offsets_bytes = Vec::new();
    offsets_bytes.extend_from_slice(&100u32.to_le_bytes());
    offsets_bytes.extend_from_slice(&150u32.to_le_bytes());
    dir_consecutive.add_entry(TIFF_TAG_STRIP_OFFSETS, TiffType::Long, 2, 0, offsets_bytes);

    let mut counts_bytes = Vec::new();
    counts_bytes.extend_from_slice(&50u32.to_le_bytes());
    counts_bytes.extend_from_slice(&50u32.to_le_bytes());
    dir_consecutive.add_entry(TIFF_TAG_STRIP_BYTE_COUNTS, TiffType::Long, 2, 0, counts_bytes);
    dir_consecutive.add_entry(
        TIFF_TAG_COMPRESSION,
        TiffType::Short,
        1,
        0,
        1u16.to_le_bytes().to_vec(),
    ); // Uncompressed

    let img = parser.get_image_data(&dir_consecutive).unwrap();
    assert_eq!(img.offset, 100);
    assert_eq!(img.length, 100); // 50 + 50
    assert_eq!(img.format, ImageFormat::UncompressedRgb);

    // Non-consecutive strips: strip 0 at 100 (len 50), strip 1 at 200 (gap!)
    let mut dir_gap = TiffDirectory::new(Endian::Little);
    let mut gap_offsets = Vec::new();
    gap_offsets.extend_from_slice(&100u32.to_le_bytes());
    gap_offsets.extend_from_slice(&200u32.to_le_bytes());
    dir_gap.add_entry(TIFF_TAG_STRIP_OFFSETS, TiffType::Long, 2, 0, gap_offsets);

    let mut gap_counts = Vec::new();
    gap_counts.extend_from_slice(&50u32.to_le_bytes());
    gap_counts.extend_from_slice(&50u32.to_le_bytes());
    dir_gap.add_entry(TIFF_TAG_STRIP_BYTE_COUNTS, TiffType::Long, 2, 0, gap_counts);
    dir_gap.add_entry(TIFF_TAG_COMPRESSION, TiffType::Short, 1, 0, 1u16.to_le_bytes().to_vec());

    assert_eq!(
        parser.get_image_data(&dir_gap),
        Err(TiffError::TagNotFound { tag: TIFF_TAG_STRIP_OFFSETS })
    );
}

// ============================================================================
// 6. Exif, SubIFD, and GPS IFD Error Propagation Tests (DIV-07 & DIV-08)
// ============================================================================

#[test]
fn test_corrupted_exif_ifd_fails_parse() {
    let mut buf = Vec::new();
    buf.extend_from_slice(b"II\x2A\x00\x08\x00\x00\x00"); // Header: LE, IFD0 @ 8
    buf.extend_from_slice(&1u16.to_le_bytes()); // 1 entry
    buf.extend_from_slice(&(EXIF_TAG_EXIF_IFD as u16).to_le_bytes()); // Tag 0x8769
    buf.extend_from_slice(&(TiffType::Long as u16).to_le_bytes()); // Type 4
    buf.extend_from_slice(&1u32.to_le_bytes()); // Count 1
    buf.extend_from_slice(&99999u32.to_le_bytes()); // Offset 99999 (out of bounds)
    buf.extend_from_slice(&0u32.to_le_bytes()); // Next IFD = 0

    let stream = ByteSliceStream::new(&buf);
    let parser = TiffParser::new(&stream).unwrap();
    let desired = [EXIF_TAG_EXIF_IFD];
    assert!(parser.parse(&desired, 1).is_err());
}

#[test]
fn test_corrupted_sub_ifd_fails_parse() {
    let mut buf = Vec::new();
    buf.extend_from_slice(b"II\x2A\x00\x08\x00\x00\x00"); // Header: LE, IFD0 @ 8
    buf.extend_from_slice(&1u16.to_le_bytes()); // 1 entry
    buf.extend_from_slice(&(EXIF_TAG_SUB_IFDS as u16).to_le_bytes()); // Tag 0x014A
    buf.extend_from_slice(&(TiffType::Long as u16).to_le_bytes()); // Type 4
    buf.extend_from_slice(&1u32.to_le_bytes()); // Count 1
    buf.extend_from_slice(&99999u32.to_le_bytes()); // Offset 99999 (out of bounds)
    buf.extend_from_slice(&0u32.to_le_bytes()); // Next IFD = 0

    let stream = ByteSliceStream::new(&buf);
    let parser = TiffParser::new(&stream).unwrap();
    let desired = [EXIF_TAG_SUB_IFDS];
    assert!(parser.parse(&desired, 1).is_err());
}

#[test]
fn test_corrupted_gps_ifd_fails_parse() {
    let mut buf = Vec::new();
    buf.extend_from_slice(b"II\x2A\x00\x08\x00\x00\x00"); // Header: LE, IFD0 @ 8
    buf.extend_from_slice(&1u16.to_le_bytes()); // 1 entry
    buf.extend_from_slice(&(EXIF_TAG_GPS as u16).to_le_bytes()); // Tag 0x8825
    buf.extend_from_slice(&(TiffType::Long as u16).to_le_bytes()); // Type 4
    buf.extend_from_slice(&1u32.to_le_bytes()); // Count 1
    buf.extend_from_slice(&99999u32.to_le_bytes()); // Offset 99999 (out of bounds)
    buf.extend_from_slice(&0u32.to_le_bytes()); // Next IFD = 0

    let stream = ByteSliceStream::new(&buf);
    let parser = TiffParser::new(&stream).unwrap();
    let desired = [EXIF_TAG_GPS];
    assert!(parser.parse(&desired, 1).is_err());
}

#[test]
fn test_valid_exif_with_corrupted_gps_ifd_fails_parse() {
    let mut writer = TiffWriter::new(Endian::Little);
    writer.write_header(8);

    // IFD0 with 2 entries: EXIF_TAG_EXIF_IFD -> 38, EXIF_TAG_GPS -> 99999
    writer.buf.extend_from_slice(&writer.pack_u16(2));
    writer.write_entry_offset(EXIF_TAG_EXIF_IFD as u16, TiffType::Long, 1, 38);
    writer.write_entry_offset(EXIF_TAG_GPS as u16, TiffType::Long, 1, 99999);
    writer.buf.extend_from_slice(&writer.pack_u32(0)); // Next IFD = 0

    // Pad up to offset 38
    if writer.buf.len() < 38 {
        writer.buf.resize(38, 0);
    }
    // Exif IFD at offset 38: 0 entries, next IFD 0
    writer.buf.extend_from_slice(&writer.pack_u16(0));
    writer.buf.extend_from_slice(&writer.pack_u32(0));

    let stream = ByteSliceStream::new(&writer.buf);
    let parser = TiffParser::new(&stream).unwrap();
    let desired = [EXIF_TAG_EXIF_IFD, EXIF_TAG_GPS];
    assert!(parser.parse(&desired, 1).is_err());
}

// ============================================================================
// 7. get_full_dimension32 DefaultCropSize Tests (DIV-09)
// ============================================================================

#[test]
fn test_get_full_dimension32_default_crop_size_invalid_does_not_fall_through() {
    let mut header = Vec::new();
    header.extend_from_slice(b"II\x2A\x00\x08\x00\x00\x00");
    let stream = ByteSliceStream::new(&header);
    let parser = TiffParser::new(&stream).unwrap();

    let mut dir = TiffDirectory::new(Endian::Little);
    // Add fallback tags that would succeed if fallthrough occurred
    dir.add_entry(EXIF_TAG_WIDTH, TiffType::Long, 1, 0, 6000u32.to_le_bytes().to_vec());
    dir.add_entry(EXIF_TAG_HEIGHT, TiffType::Long, 1, 0, 4000u32.to_le_bytes().to_vec());
    dir.add_entry(TIFF_TAG_IMAGE_WIDTH, TiffType::Long, 1, 0, 6000u32.to_le_bytes().to_vec());
    dir.add_entry(TIFF_TAG_IMAGE_LENGTH, TiffType::Long, 1, 0, 4000u32.to_le_bytes().to_vec());

    // 1. Invalid DefaultCropSize: count = 1 (expected 2)
    dir.add_entry(EXIF_TAG_DEFAULT_CROP_SIZE, TiffType::Long, 1, 0, 1800u32.to_le_bytes().to_vec());
    let res = parser.get_full_dimension32(&dir);
    assert!(res.is_err(), "Expected error when DefaultCropSize is invalid, but got {:?}", res);

    // 2. Invalid DefaultCropSize: count = 2 but width is 0
    let mut dir_zero = TiffDirectory::new(Endian::Little);
    dir_zero.add_entry(EXIF_TAG_WIDTH, TiffType::Long, 1, 0, 6000u32.to_le_bytes().to_vec());
    dir_zero.add_entry(EXIF_TAG_HEIGHT, TiffType::Long, 1, 0, 4000u32.to_le_bytes().to_vec());
    let mut zero_crop = Vec::new();
    zero_crop.extend_from_slice(&0u32.to_le_bytes());
    zero_crop.extend_from_slice(&1200u32.to_le_bytes());
    dir_zero.add_entry(EXIF_TAG_DEFAULT_CROP_SIZE, TiffType::Long, 2, 0, zero_crop);
    let res_zero = parser.get_full_dimension32(&dir_zero);
    assert!(
        res_zero.is_err(),
        "Expected error when DefaultCropSize has 0 dimension, but got {:?}",
        res_zero
    );
}

#[test]
fn test_get_full_dimension32_default_crop_size_valid_overrides_exif_dimensions() {
    let mut header = Vec::new();
    header.extend_from_slice(b"II\x2A\x00\x08\x00\x00\x00");
    let stream = ByteSliceStream::new(&header);
    let parser = TiffParser::new(&stream).unwrap();

    let mut dir = TiffDirectory::new(Endian::Little);
    dir.add_entry(EXIF_TAG_WIDTH, TiffType::Long, 1, 0, 6000u32.to_le_bytes().to_vec());
    dir.add_entry(EXIF_TAG_HEIGHT, TiffType::Long, 1, 0, 4000u32.to_le_bytes().to_vec());

    let mut crop_bytes = Vec::new();
    crop_bytes.extend_from_slice(&1800u32.to_le_bytes());
    crop_bytes.extend_from_slice(&1200u32.to_le_bytes());
    dir.add_entry(EXIF_TAG_DEFAULT_CROP_SIZE, TiffType::Long, 2, 0, crop_bytes);

    let (w, h) = parser.get_full_dimension32(&dir).unwrap();
    assert_eq!(w, 1800);
    assert_eq!(h, 1200);
}

#[test]
fn test_get_full_dimension32_default_crop_size_rational_valid() {
    let mut header = Vec::new();
    header.extend_from_slice(b"II\x2A\x00\x08\x00\x00\x00");
    let stream = ByteSliceStream::new(&header);
    let parser = TiffParser::new(&stream).unwrap();

    let mut dir = TiffDirectory::new(Endian::Little);
    dir.add_entry(EXIF_TAG_WIDTH, TiffType::Long, 1, 0, 6000u32.to_le_bytes().to_vec());
    dir.add_entry(EXIF_TAG_HEIGHT, TiffType::Long, 1, 0, 4000u32.to_le_bytes().to_vec());

    let mut rat_bytes = Vec::new();
    rat_bytes.extend_from_slice(&3600u32.to_le_bytes()); // num1
    rat_bytes.extend_from_slice(&2u32.to_le_bytes()); // den1 -> 1800
    rat_bytes.extend_from_slice(&2400u32.to_le_bytes()); // num2
    rat_bytes.extend_from_slice(&2u32.to_le_bytes()); // den2 -> 1200
    dir.add_entry(EXIF_TAG_DEFAULT_CROP_SIZE, TiffType::Rational, 2, 0, rat_bytes);

    let (w, h) = parser.get_full_dimension32(&dir).unwrap();
    assert_eq!(w, 1800);
    assert_eq!(h, 1200);
}

#[test]
fn test_regression_crash_04_tiff_parser_truncated_ifd_bounds() {
    // Truncated TIFF header with IFD count stating 100 entries but stream ending prematurely
    let truncated_tiff = b"II*\x00\x08\x00\x00\x00\x64\x00\x00\x00";
    let stream = ByteSliceStream::new(truncated_tiff);

    let parser = TiffParser::new(&stream);
    assert!(parser.is_ok());

    let parser = parser.unwrap();
    let desired_tags: [u32; 0] = [];
    let res = parser.parse(&desired_tags, 1);

    // Should return error for truncated IFD entries rather than panicking or out of bounds access
    assert!(res.is_err());
}

#[test]
fn test_parse_sub_ifd_with_tiff_ifd_type_single_and_multiple() {
    for endian in [Endian::Little, Endian::Big] {
        // 1. Single inline SubIFD encoded with TiffType::Ifd (13).
        let mut writer = TiffWriter::new(endian);
        writer.write_header(8);
        let sub_ifd_offset = 26u32;
        writer.buf.extend_from_slice(&writer.pack_u16(1));
        writer.write_entry_offset(EXIF_TAG_SUB_IFDS as u16, TiffType::Ifd, 1, sub_ifd_offset);
        writer.buf.extend_from_slice(&writer.pack_u32(0));
        assert_eq!(writer.buf.len(), sub_ifd_offset as usize);

        writer.buf.extend_from_slice(&writer.pack_u16(4));
        writer.write_entry_offset(TIFF_TAG_IMAGE_WIDTH as u16, TiffType::Long, 1, 4000);
        writer.write_entry_offset(TIFF_TAG_IMAGE_LENGTH as u16, TiffType::Long, 1, 3000);
        writer.write_entry_offset(TIFF_TAG_JPEG_OFFSET as u16, TiffType::Long, 1, 7);
        writer.write_entry_offset(TIFF_TAG_JPEG_BYTE_COUNT as u16, TiffType::Long, 1, 6);
        writer.buf.extend_from_slice(&writer.pack_u32(0));

        let stream = ByteSliceStream::new(&writer.buf);
        let parser = TiffParser::new(&stream).unwrap();
        let desired = [
            EXIF_TAG_SUB_IFDS,
            TIFF_TAG_IMAGE_WIDTH,
            TIFF_TAG_IMAGE_LENGTH,
            TIFF_TAG_JPEG_OFFSET,
            TIFF_TAG_JPEG_BYTE_COUNT,
        ];
        let content = parser.parse(&desired, 2).unwrap();
        assert_eq!(content.tiff_directories.len(), 1);
        assert_eq!(content.tiff_directories[0].sub_directories().len(), 1);

        let mut preview_data = PreviewImageData::default();
        assert!(parser.get_preview_image_data(&content, &mut preview_data).is_ok());
        assert_eq!(preview_data.full_width, 4000);
        assert_eq!(preview_data.full_height, 3000);
        assert_eq!(preview_data.thumbnail.offset, 7);
        assert_eq!(preview_data.thumbnail.length, 6);

        // 2. Multiple out-of-line SubIFDs encoded with TiffType::Ifd (13).
        let mut multi = TiffWriter::new(endian);
        multi.write_header(8);
        let offsets_array_pos = 26u32;
        let sub_ifd_0_offset = 34u32;
        let sub_ifd_1_offset = sub_ifd_0_offset + 18;
        multi.buf.extend_from_slice(&multi.pack_u16(1));
        multi.write_entry_offset(EXIF_TAG_SUB_IFDS as u16, TiffType::Ifd, 2, offsets_array_pos);
        multi.buf.extend_from_slice(&multi.pack_u32(0));
        assert_eq!(multi.buf.len(), offsets_array_pos as usize);

        let off0 = multi.pack_u32(sub_ifd_0_offset);
        let off1 = multi.pack_u32(sub_ifd_1_offset);
        multi.buf.extend_from_slice(&off0);
        multi.buf.extend_from_slice(&off1);
        assert_eq!(multi.buf.len(), sub_ifd_0_offset as usize);

        multi.buf.extend_from_slice(&multi.pack_u16(1));
        multi.write_entry_offset(TIFF_TAG_SUB_FILE_TYPE as u16, TiffType::Long, 1, 1);
        multi.buf.extend_from_slice(&multi.pack_u32(0));
        assert_eq!(multi.buf.len(), sub_ifd_1_offset as usize);

        multi.buf.extend_from_slice(&multi.pack_u16(3));
        multi.write_entry_offset(TIFF_TAG_SUB_FILE_TYPE as u16, TiffType::Long, 1, 0);
        multi.write_entry_offset(TIFF_TAG_IMAGE_WIDTH as u16, TiffType::Long, 1, 6000);
        multi.write_entry_offset(TIFF_TAG_IMAGE_LENGTH as u16, TiffType::Long, 1, 4000);
        multi.buf.extend_from_slice(&multi.pack_u32(0));

        let multi_stream = ByteSliceStream::new(&multi.buf);
        let multi_parser = TiffParser::new(&multi_stream).unwrap();
        let multi_desired = [
            EXIF_TAG_SUB_IFDS,
            TIFF_TAG_SUB_FILE_TYPE,
            TIFF_TAG_IMAGE_WIDTH,
            TIFF_TAG_IMAGE_LENGTH,
        ];
        let multi_content = multi_parser.parse(&multi_desired, 2).unwrap();
        assert_eq!(multi_content.tiff_directories.len(), 1);
        assert_eq!(multi_content.tiff_directories[0].sub_directories().len(), 2);

        let mut multi_preview = PreviewImageData::default();
        assert!(multi_parser.get_preview_image_data(&multi_content, &mut multi_preview).is_ok());
        assert_eq!(multi_preview.full_width, 6000);
        assert_eq!(multi_preview.full_height, 4000);

        // 3. TiffType::Short (3) for EXIF_TAG_SUB_IFDS is ignored (matching C++ ParseSubIfds).
        let mut short_writer = TiffWriter::new(endian);
        short_writer.write_header(8);
        short_writer.buf.extend_from_slice(&short_writer.pack_u16(1));
        let packed_short = short_writer.pack_u16(26);
        short_writer.write_entry_inline(
            EXIF_TAG_SUB_IFDS as u16,
            TiffType::Short,
            1,
            &packed_short,
        );
        short_writer.buf.extend_from_slice(&short_writer.pack_u32(0));
        let short_stream = ByteSliceStream::new(&short_writer.buf);
        let short_parser = TiffParser::new(&short_stream).unwrap();
        let short_content = short_parser.parse(&[EXIF_TAG_SUB_IFDS], 2).unwrap();
        assert!(short_content.tiff_directories[0].sub_directories().is_empty());

        // 4. Non-zero tiff_offset resolves SubIFDs relative to tiff_offset and rejects u32
        //    overflow.
        let tiff_offset = 16u32;
        let mut embedded = vec![0xAAu8; tiff_offset as usize];
        embedded.extend_from_slice(&writer.buf);
        let embedded_stream = ByteSliceStream::new(&embedded);
        let embedded_parser = TiffParser::with_offset(&embedded_stream, tiff_offset).unwrap();
        let embedded_content = embedded_parser
            .parse(&[EXIF_TAG_SUB_IFDS, TIFF_TAG_IMAGE_WIDTH, TIFF_TAG_IMAGE_LENGTH], 2)
            .unwrap();
        assert_eq!(embedded_content.tiff_directories[0].sub_directories().len(), 1);

        // 5. Cyclic next_ifd_offset pointing back to IFD0 stops after 1 IFD.
        let mut cyclic = TiffWriter::new(endian);
        cyclic.write_header(8);
        cyclic.buf.extend_from_slice(&cyclic.pack_u16(1));
        cyclic.write_entry_offset(TIFF_TAG_IMAGE_WIDTH as u16, TiffType::Long, 1, 3200);
        cyclic.buf.extend_from_slice(&cyclic.pack_u32(8));
        let cyclic_stream = ByteSliceStream::new(&cyclic.buf);
        let cyclic_parser = TiffParser::new(&cyclic_stream).unwrap();
        let cyclic_content = cyclic_parser.parse(&[TIFF_TAG_IMAGE_WIDTH], 2).unwrap();
        assert_eq!(cyclic_content.tiff_directories.len(), 1);
    }
}

#[test]
fn test_cfa_pattern_dim_set_only_when_two_element_tag_present() {
    let mut header = Vec::new();
    header.extend_from_slice(b"II\x2A\x00\x08\x00\x00\x00");
    let stream = ByteSliceStream::new(&header);

    // 1. Tag absent -> cfa_pattern_dim_set remains false and caller value untouched.
    {
        let dir = TiffDirectory::new(Endian::Little);
        let mut content = TiffContent::new();
        content.tiff_directories.push(dir);
        let parser = TiffParser::new(&stream).unwrap();
        let mut preview_data = PreviewImageData::default();
        preview_data.cfa_pattern_dim = [9, 8];
        parser.get_preview_image_data(&content, &mut preview_data).unwrap();
        assert!(!preview_data.cfa_pattern_dim_set);
        assert_eq!(preview_data.cfa_pattern_dim, [9, 8]);
    }

    // 2. Tag present with 1 element -> cfa_pattern_dim_set remains false.
    {
        let mut dir = TiffDirectory::new(Endian::Little);
        dir.add_entry(TIFF_TAG_CFA_PATTERN_DIM, TiffType::Short, 1, 0, 2u16.to_le_bytes().to_vec());
        let mut content = TiffContent::new();
        content.tiff_directories.push(dir);
        let parser = TiffParser::new(&stream).unwrap();
        let mut preview_data = PreviewImageData::default();
        preview_data.cfa_pattern_dim = [9, 8];
        parser.get_preview_image_data(&content, &mut preview_data).unwrap();
        assert!(!preview_data.cfa_pattern_dim_set);
        assert_eq!(preview_data.cfa_pattern_dim, [9, 8]);
    }

    // 3. Tag present with 2 elements -> cfa_pattern_dim_set is true and values updated.
    {
        let mut dir = TiffDirectory::new(Endian::Little);
        let mut cfa_bytes = Vec::new();
        cfa_bytes.extend_from_slice(&2u16.to_le_bytes());
        cfa_bytes.extend_from_slice(&4u16.to_le_bytes());
        dir.add_entry(TIFF_TAG_CFA_PATTERN_DIM, TiffType::Short, 2, 0, cfa_bytes);
        let mut content = TiffContent::new();
        content.tiff_directories.push(dir);
        let parser = TiffParser::new(&stream).unwrap();
        let mut preview_data = PreviewImageData::default();
        parser.get_preview_image_data(&content, &mut preview_data).unwrap();
        assert!(preview_data.cfa_pattern_dim_set);
        assert_eq!(preview_data.cfa_pattern_dim, [2, 4]);
    }
}
