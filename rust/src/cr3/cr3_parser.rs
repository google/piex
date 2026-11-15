#![forbid(unsafe_code)]

//! Canon CR3 RAW image container parser.
//!
//! Parses ISOBMFF structures in Canon CR3 files to extract JPEG/HEVC preview images,
//! thumbnails, and EXIF metadata (camera model, make, exposure parameters, orientation).

use crate::cr3::box_reader::{
    find_next_box_with_tag, read_next_box, BoxHeader, BoxType, UUID_MOOV, UUID_PRVW,
};
use crate::stream::StreamInterface;
use crate::tiff::tags::tiff::*;
use crate::tiff::{TiffError, TiffParser};
use crate::types::{Image, ImageFormat, PiexError, PreviewImageData};

/// Offsets for width, height, JPEG size, and JPEG data within PRVW and THMB boxes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct JpegBoxOffsets {
    /// Byte offset of the 16-bit width value from box start.
    pub width: u64,
    /// Byte offset of the 16-bit height value from box start.
    pub height: u64,
    /// Byte offset of the 32-bit JPEG byte count from box start.
    pub jpeg_size: u64,
    /// Byte offset of the start of JPEG payload data from box start.
    pub jpeg_data: u64,
}

/// Offsets for `PRVW` preview image box.
pub const PRVW_JPEG_OFFSETS: JpegBoxOffsets =
    JpegBoxOffsets { width: 14, height: 16, jpeg_size: 20, jpeg_data: 24 };

/// Offsets for `THMB` thumbnail image box.
pub const THMB_JPEG_OFFSETS: JpegBoxOffsets =
    JpegBoxOffsets { width: 12, height: 14, jpeg_size: 16, jpeg_data: 24 };

/// Internal parsing context tracking candidate embedded images found in CR3 streams.
#[derive(Default, Debug, Clone)]
struct Cr3ProcessContext {
    mdat_image: Option<Image>,
    prvw_image: Option<Image>,
}

fn read_u16_be<R: StreamInterface + ?Sized>(stream: &R, offset: u64) -> Result<u16, PiexError> {
    let offset_usize = usize::try_from(offset).map_err(|_| PiexError::OutOfBounds {
        requested_offset: usize::MAX,
        requested_size: 2,
        total_len: stream.len(),
    })?;
    let mut buf = [0u8; 2];
    stream.get_data(offset_usize, &mut buf)?;
    Ok(u16::from_be_bytes(buf))
}

fn read_u32_be<R: StreamInterface + ?Sized>(stream: &R, offset: u64) -> Result<u32, PiexError> {
    let offset_usize = usize::try_from(offset).map_err(|_| PiexError::OutOfBounds {
        requested_offset: usize::MAX,
        requested_size: 4,
        total_len: stream.len(),
    })?;
    let mut buf = [0u8; 4];
    stream.get_data(offset_usize, &mut buf)?;
    Ok(u32::from_be_bytes(buf))
}

#[allow(dead_code)]
fn read_u64_be<R: StreamInterface + ?Sized>(stream: &R, offset: u64) -> Result<u64, PiexError> {
    let offset_usize = usize::try_from(offset).map_err(|_| PiexError::OutOfBounds {
        requested_offset: usize::MAX,
        requested_size: 8,
        total_len: stream.len(),
    })?;
    let mut buf = [0u8; 8];
    stream.get_data(offset_usize, &mut buf)?;
    Ok(u64::from_be_bytes(buf))
}

/// Returns `true` if child boxes within the given box tag should be processed recursively.
fn should_process_children(tag: BoxType) -> bool {
    matches!(tag, BoxType::TRAK | BoxType::MOOV | BoxType::MDIA | BoxType::MINF)
}

/// Extracts embedded JPEG properties (dimensions, offset, length) from a `PRVW` or `THMB` box.
fn process_jpeg_box<R: StreamInterface + ?Sized>(
    stream: &R,
    box_header: &BoxHeader,
    offsets: &JpegBoxOffsets,
) -> Result<Image, PiexError> {
    let width_offset = box_header
        .offset
        .checked_add(offsets.width)
        .ok_or_else(|| PiexError::CorruptData("JPEG box width offset overflow".to_string()))?;
    let width = read_u16_be(stream, width_offset)?;

    let height_offset = box_header
        .offset
        .checked_add(offsets.height)
        .ok_or_else(|| PiexError::CorruptData("JPEG box height offset overflow".to_string()))?;
    let height = read_u16_be(stream, height_offset)?;

    let jpeg_size_offset = box_header
        .offset
        .checked_add(offsets.jpeg_size)
        .ok_or_else(|| PiexError::CorruptData("JPEG box size offset overflow".to_string()))?;
    let jpeg_size = read_u32_be(stream, jpeg_size_offset)?;

    let jpeg_data_offset = box_header
        .offset
        .checked_add(offsets.jpeg_data)
        .ok_or_else(|| PiexError::CorruptData("JPEG box data offset overflow".to_string()))?;

    let image_offset = u32::try_from(jpeg_data_offset).map_err(|_| {
        PiexError::CorruptData("JPEG image offset overflowed 32-bit integer".to_string())
    })?;

    Ok(Image {
        format: ImageFormat::JpegCompressed,
        width,
        height,
        offset: image_offset,
        length: jpeg_size,
    })
}

/// Byte offset of the 16-bit image width field from the start of a `CRAW` box.
const CRAW_WIDTH_OFFSET: u64 = 32;
/// Byte offset of the 16-bit image height field from the start of a `CRAW` box.
const CRAW_HEIGHT_OFFSET: u64 = 34;
/// Byte offset of the 16-bit content type field from the start of a `CRAW` box.
const CRAW_CONTENT_TYPE_OFFSET: u64 = 86;
/// Minimum total byte length of a `CRAW` box covering width, height, and content type.
const CRAW_MIN_BOX_LEN: u64 = CRAW_CONTENT_TYPE_OFFSET + std::mem::size_of::<u16>() as u64;

/// Reads image width, height, and content type from a Canon RAW (`CRAW`) box.
fn process_craw_box<R: StreamInterface + ?Sized>(
    stream: &R,
    craw_box: &BoxHeader,
) -> Result<(u16, u16, u16), PiexError> {
    let box_len = craw_box
        .next_box_offset
        .checked_sub(craw_box.offset)
        .ok_or(PiexError::InvalidBoxHeader { offset: craw_box.offset })?;
    if box_len < CRAW_MIN_BOX_LEN {
        return Err(PiexError::CorruptData(
            "CRAW box too small to contain image metadata".to_string(),
        ));
    }

    let width_offset = craw_box
        .offset
        .checked_add(CRAW_WIDTH_OFFSET)
        .ok_or_else(|| PiexError::CorruptData("CRAW box width offset overflow".to_string()))?;
    let width = read_u16_be(stream, width_offset)?;

    let height_offset = craw_box
        .offset
        .checked_add(CRAW_HEIGHT_OFFSET)
        .ok_or_else(|| PiexError::CorruptData("CRAW box height offset overflow".to_string()))?;
    let height = read_u16_be(stream, height_offset)?;

    let content_type_offset =
        craw_box.offset.checked_add(CRAW_CONTENT_TYPE_OFFSET).ok_or_else(|| {
            PiexError::CorruptData("CRAW box content_type offset overflow".to_string())
        })?;
    let content_type = read_u16_be(stream, content_type_offset)?;

    Ok((width, height, content_type))
}

/// Reads sample size from a Sample Size (`stsz`) box.
fn process_stsz_box<R: StreamInterface + ?Sized>(
    stream: &R,
    stsz_box: &BoxHeader,
) -> Result<u32, PiexError> {
    let sample_size_offset = stsz_box.offset.checked_add(12).ok_or_else(|| {
        PiexError::CorruptData("stsz box sample_size offset overflow".to_string())
    })?;
    let sample_size = read_u32_be(stream, sample_size_offset)?;
    if sample_size > 0 {
        return Ok(sample_size);
    }

    let count_offset = stsz_box
        .offset
        .checked_add(16)
        .ok_or_else(|| PiexError::CorruptData("stsz box count offset overflow".to_string()))?;
    let count = read_u32_be(stream, count_offset)?;
    if count != 1 {
        return Err(PiexError::CorruptData(
            "Expected stsz sample count == 1 when sample_size is 0".to_string(),
        ));
    }

    let entry0_offset = stsz_box
        .offset
        .checked_add(20)
        .ok_or_else(|| PiexError::CorruptData("stsz box entry0 offset overflow".to_string()))?;
    read_u32_be(stream, entry0_offset)
}

/// Reads 64-bit chunk offset from a 64-bit Chunk Offset (`co64`) box.
fn process_co64_box<R: StreamInterface + ?Sized>(
    stream: &R,
    co64_box: &BoxHeader,
) -> Result<u32, PiexError> {
    let count_offset = co64_box
        .header_offset
        .checked_add(4)
        .ok_or_else(|| PiexError::CorruptData("co64 box count offset overflow".to_string()))?;
    let count = read_u32_be(stream, count_offset)?;
    if count != 1 {
        return Err(PiexError::CorruptData("Expected co64 entry count == 1".to_string()));
    }

    let offset_offset = co64_box
        .header_offset
        .checked_add(8)
        .ok_or_else(|| PiexError::CorruptData("co64 box offset offset overflow".to_string()))?;
    read_u32_be(stream, offset_offset)
}

/// Navigates `stsd` -> `CRAW` -> `stsz` -> `co64` inside a Sample Table (`stbl`) box
/// to locate image parameters.
fn process_stbl_box<R: StreamInterface + ?Sized>(
    stream: &R,
    stbl_box: &BoxHeader,
    ctx: &mut Cr3ProcessContext,
) -> Result<(), PiexError> {
    let stsd_box =
        match find_next_box_with_tag(stream, stbl_box.header_offset, stbl_box, BoxType::STSD)? {
            Some(b) => b,
            None => return Err(PiexError::CorruptData("stsd box not found in stbl".to_string())),
        };

    let stsd_payload_offset = stsd_box
        .header_offset
        .checked_add(8)
        .ok_or_else(|| PiexError::CorruptData("stsd box payload offset overflow".to_string()))?;
    let craw_box = match read_next_box(stream, stsd_payload_offset)? {
        Some(b) => b,
        None => return Err(PiexError::CorruptData("Invalid or missing box in stsd".to_string())),
    };

    if craw_box.box_type != BoxType::CRAW {
        return Ok(());
    }
    if craw_box.next_box_offset > stsd_box.next_box_offset {
        return Err(PiexError::CorruptData(
            "CRAW sample entry exceeds enclosing stsd box".to_string(),
        ));
    }

    let (width, height, content_type) = process_craw_box(stream, &craw_box)?;

    const JPEG_CONTENT_TYPE: u16 = 3;
    const HEVC_CONTENT_TYPE: u16 = 4;

    if content_type != JPEG_CONTENT_TYPE && content_type != HEVC_CONTENT_TYPE {
        return Ok(());
    }

    let format = if content_type == JPEG_CONTENT_TYPE {
        ImageFormat::JpegCompressed
    } else {
        ImageFormat::HevcCompressed
    };

    let stsz_box =
        match find_next_box_with_tag(stream, stsd_box.next_box_offset, stbl_box, BoxType::STSZ)? {
            Some(b) => b,
            None => return Err(PiexError::CorruptData("stsz box not found in stbl".to_string())),
        };
    let sample_size = process_stsz_box(stream, &stsz_box)?;

    let co64_box =
        match find_next_box_with_tag(stream, stsz_box.next_box_offset, stbl_box, BoxType::CO64)? {
            Some(b) => b,
            None => return Err(PiexError::CorruptData("co64 box not found in stbl".to_string())),
        };
    let chunk_offset = process_co64_box(stream, &co64_box)?;

    ctx.mdat_image =
        Some(Image { format, width, height, offset: chunk_offset, length: sample_size });

    Ok(())
}

/// Parses IFD0 metadata (camera make, model, orientation, dimensions) from `CMT1` box payload.
fn parse_exif_ifd0<R: StreamInterface + ?Sized>(
    stream: &R,
    box_header: &BoxHeader,
    data: &mut PreviewImageData,
) -> Result<(), PiexError> {
    let tiff_offset = u32::try_from(box_header.header_offset)
        .map_err(|_| PiexError::CorruptData("TIFF offset overflow in CMT1 box".to_string()))?;

    let parser = TiffParser::with_offset(stream, tiff_offset).map_err(|e| match e {
        TiffError::StreamError(se) => se,
        other => PiexError::CorruptData(other.to_string()),
    })?;

    const CMT1_TAGS: [u32; 5] = [
        TIFF_TAG_MODEL,
        TIFF_TAG_MAKE,
        TIFF_TAG_ORIENTATION,
        TIFF_TAG_IMAGE_WIDTH,
        TIFF_TAG_IMAGE_LENGTH,
    ];

    let content = parser.parse(&CMT1_TAGS, 1).map_err(|e| match e {
        TiffError::StreamError(se) => se,
        other => PiexError::CorruptData(other.to_string()),
    })?;

    if content.tiff_directories.len() != 1 {
        return Err(PiexError::CorruptData(
            "Expected exactly 1 TIFF directory in CMT1".to_string(),
        ));
    }

    let dir = &content.tiff_directories[0];
    if let Ok(model) = dir.get_ascii_bytes(TIFF_TAG_MODEL) {
        data.model = model.to_vec();
    }
    if let Ok(maker) = dir.get_ascii_bytes(TIFF_TAG_MAKE) {
        data.maker = maker.to_vec();
    }
    if let Ok(orientation) = dir.get_u32(TIFF_TAG_ORIENTATION) {
        data.exif_orientation = orientation;
    }
    if let Ok(width) = dir.get_u32(TIFF_TAG_IMAGE_WIDTH) {
        data.full_width = width;
    }
    if let Ok(height) = dir.get_u32(TIFF_TAG_IMAGE_LENGTH) {
        data.full_height = height;
    }

    Ok(())
}

/// Parses EXIF IFD metadata (date/time, exposure parameters, ISO rating) from `CMT2` box payload.
fn parse_exif_exif_ifd<R: StreamInterface + ?Sized>(
    stream: &R,
    box_header: &BoxHeader,
    data: &mut PreviewImageData,
) -> Result<(), PiexError> {
    let tiff_offset = u32::try_from(box_header.header_offset)
        .map_err(|_| PiexError::CorruptData("TIFF offset overflow in CMT2 box".to_string()))?;

    let parser = TiffParser::with_offset(stream, tiff_offset).map_err(|e| match e {
        TiffError::StreamError(se) => se,
        other => PiexError::CorruptData(other.to_string()),
    })?;

    const CMT2_TAGS: [u32; 5] = [
        EXIF_TAG_DATE_TIME_ORIGINAL,
        EXIF_TAG_EXPOSURE_TIME,
        EXIF_TAG_FNUMBER,
        EXIF_TAG_FOCAL_LENGTH,
        EXIF_TAG_ISO_SPEED,
    ];

    let content = parser.parse(&CMT2_TAGS, 1).map_err(|e| match e {
        TiffError::StreamError(se) => se,
        other => PiexError::CorruptData(other.to_string()),
    })?;

    if content.tiff_directories.len() != 1 {
        return Err(PiexError::CorruptData(
            "Expected exactly 1 TIFF directory in CMT2".to_string(),
        ));
    }

    let dir = &content.tiff_directories[0];
    if let Ok(date_time) = dir.get_ascii_bytes(EXIF_TAG_DATE_TIME_ORIGINAL) {
        data.date_time = date_time.to_vec();
    }
    if let Ok(exposure_time) = dir.get_rational(EXIF_TAG_EXPOSURE_TIME) {
        data.exposure_time = exposure_time;
    }
    if let Ok(fnumber) = dir.get_rational(EXIF_TAG_FNUMBER) {
        data.fnumber = fnumber;
    }
    if let Ok(focal_length) = dir.get_rational(EXIF_TAG_FOCAL_LENGTH) {
        data.focal_length = focal_length;
    }
    if let Ok(iso) = dir.get_u32(EXIF_TAG_ISO_SPEED) {
        data.iso = iso;
    }

    Ok(())
}

/// Payload prefix byte count in a Canon `uuid` MOOV container box (16-byte UUID).
const UUID_MOOV_PREFIX_LEN: u64 = UUID_MOOV.len() as u64;

/// Payload prefix byte count in a Canon `uuid` PRVW container box
/// (16-byte UUID plus 8 bytes before the nested `PRVW` child box).
const UUID_PRVW_PREFIX_LEN: u64 = UUID_PRVW.len() as u64 + 8;

/// Returns the first child-box offset inside `box_header` after a `prefix_len`-byte
/// UUID container prefix, rejecting headers that do not enclose the full prefix.
fn uuid_child_offset(box_header: &BoxHeader, prefix_len: u64) -> Result<u64, PiexError> {
    let payload_len = box_header
        .next_box_offset
        .checked_sub(box_header.header_offset)
        .ok_or(PiexError::InvalidBoxHeader { offset: box_header.offset })?;
    if payload_len < prefix_len {
        return Err(PiexError::InvalidBoxHeader { offset: box_header.offset });
    }
    box_header
        .header_offset
        .checked_add(prefix_len)
        .ok_or(PiexError::InvalidBoxHeader { offset: box_header.offset })
}

/// Reads and dispatches processing for an individual ISOBMFF box at `box_header`.
fn process_box<R: StreamInterface + ?Sized>(
    stream: &R,
    box_header: &BoxHeader,
    data: &mut PreviewImageData,
    ctx: &mut Cr3ProcessContext,
) -> Result<u64, PiexError> {
    if box_header.box_type == BoxType::UUID {
        let mut uuid = [0u8; 16];
        let header_offset_usize =
            usize::try_from(box_header.header_offset).map_err(|_| PiexError::OutOfBounds {
                requested_offset: usize::MAX,
                requested_size: 16,
                total_len: stream.len(),
            })?;
        stream.get_data(header_offset_usize, &mut uuid)?;

        return match uuid {
            UUID_MOOV => uuid_child_offset(box_header, UUID_MOOV_PREFIX_LEN),
            UUID_PRVW => uuid_child_offset(box_header, UUID_PRVW_PREFIX_LEN),
            _ => Ok(box_header.next_box_offset),
        };
    }

    if should_process_children(box_header.box_type) {
        return Ok(box_header.header_offset);
    }

    if box_header.box_type == BoxType::MDAT {
        if let Some(ref mut mdat_img) = ctx.mdat_image {
            let hdr_off_32 = u32::try_from(box_header.header_offset).map_err(|_| {
                PiexError::CorruptData("MDAT box header offset overflowed 32-bit".to_string())
            })?;
            mdat_img.offset = mdat_img
                .offset
                .checked_add(hdr_off_32)
                .ok_or_else(|| PiexError::CorruptData("MDAT offset overflow".to_string()))?;
        }
    } else if box_header.box_type == BoxType::STBL {
        process_stbl_box(stream, box_header, ctx)?;
    } else if box_header.box_type == BoxType::PRVW {
        let img = process_jpeg_box(stream, box_header, &PRVW_JPEG_OFFSETS)?;
        ctx.prvw_image = Some(img);
    } else if box_header.box_type == BoxType::THMB {
        let img = process_jpeg_box(stream, box_header, &THMB_JPEG_OFFSETS)?;
        data.thumbnail = img;
    } else if box_header.box_type == BoxType::CMT1 {
        parse_exif_ifd0(stream, box_header, data)?;
    } else if box_header.box_type == BoxType::CMT2 {
        parse_exif_exif_ifd(stream, box_header, data)?;
    }

    Ok(box_header.next_box_offset)
}

/// Sequentially traverses ISOBMFF boxes from `offset = 0` until `last_chunk` is reached.
fn process_stream<R: StreamInterface + ?Sized>(
    stream: &R,
    last_chunk: BoxType,
    data: &mut PreviewImageData,
    ctx: &mut Cr3ProcessContext,
) -> Result<(), PiexError> {
    let total_len = stream.len() as u64;
    let mut offset = 0u64;

    while offset < total_len {
        let box_header = match read_next_box(stream, offset)? {
            Some(h) => h,
            None => break,
        };

        if box_header.next_box_offset <= offset {
            return Err(PiexError::InvalidBoxHeader { offset });
        }

        let current_tag = box_header.box_type;
        let new_offset = process_box(stream, &box_header, data, ctx)?;

        if new_offset <= offset {
            return Err(PiexError::InvalidBoxHeader { offset });
        }

        if current_tag == last_chunk {
            return Ok(());
        }

        offset = new_offset;
    }

    if last_chunk != BoxType::new(b"    ") {
        return Err(PiexError::CorruptData("Target box not found in CR3 stream".to_string()));
    }

    Ok(())
}

/// Validates if `image` is non-empty and has valid magic header (if JPEG) or is non-JPEG format.
fn is_valid_image<R: StreamInterface + ?Sized>(stream: &R, image: &Image) -> bool {
    if image.length == 0 {
        return false;
    }
    if image.format != ImageFormat::JpegCompressed {
        // Non-JPEG images (e.g. HEVC) are valid preview images. Pass responsibility to caller.
        return true;
    }
    let offset_usize = match usize::try_from(image.offset) {
        Ok(v) => v,
        Err(_) => return false,
    };
    let mut header = [0u8; 3];
    if stream.get_data(offset_usize, &mut header).is_err() {
        return false;
    }
    header == [0xFF, 0xD8, 0xFF]
}

/// Parses CR3 ISOBMFF stream and extracts preview image, thumbnail, and metadata
/// into `preview_image_data`.
pub fn extract_cr3_preview_data_into<R: StreamInterface + ?Sized>(
    stream: &R,
    preview_image_data: &mut PreviewImageData,
) -> Result<(), PiexError> {
    if stream.len() == 0 {
        return Err(PiexError::GenericFail("Empty stream provided to CR3 parser".to_string()));
    }
    let mut ctx = Cr3ProcessContext::default();

    process_stream(stream, BoxType::MDAT, preview_image_data, &mut ctx)?;

    if let Some(ref mdat) = ctx.mdat_image {
        if mdat.length > 0 && is_valid_image(stream, mdat) {
            preview_image_data.preview = mdat.clone();
            return Ok(());
        }
    }

    if let Some(ref prvw) = ctx.prvw_image {
        if prvw.length > 0 && is_valid_image(stream, prvw) {
            preview_image_data.preview = prvw.clone();
            return Ok(());
        }
    }

    Err(PiexError::GenericFail("No valid preview image found in CR3 stream".to_string()))
}

/// Parses CR3 ISOBMFF stream and returns extracted `PreviewImageData`.
pub fn extract_cr3_preview_data<R: StreamInterface + ?Sized>(
    stream: &R,
) -> Result<PreviewImageData, PiexError> {
    let mut preview_image_data = PreviewImageData::default();
    extract_cr3_preview_data_into(stream, &mut preview_image_data)?;
    Ok(preview_image_data)
}

/// Extracts EXIF orientation value from CR3 stream.
///
/// Scans stream up to `CMT1` box and returns `exif_orientation`.
///
/// # Errors
/// Returns `PiexError` if the stream is corrupted or lacks orientation metadata.
pub fn extract_cr3_orientation<R: StreamInterface + ?Sized>(stream: &R) -> Result<u32, PiexError> {
    let mut preview_image_data = PreviewImageData::default();
    let mut ctx = Cr3ProcessContext::default();

    process_stream(stream, BoxType::CMT1, &mut preview_image_data, &mut ctx)?;

    Ok(preview_image_data.exif_orientation)
}
