#![forbid(unsafe_code)]

//! Core RAW preview extraction pipeline and helper utilities.
//!
//! Provides functions to identify RAW file streams, check format support requirements,
//! retrieve supported file extension lists, and extract embedded preview images
//! and EXIF/TIFF metadata for standard TIFF-based formats (ARW, CR2, NEF, PEF, SRW)
//! and complex formats (DNG, CR3, ORF, RAF, RW2).

use std::io::{Read, Seek, SeekFrom};

use crate::format_recognizer::{
    get_number_of_bytes_for_is_raw_lite, recognize_raw_image_type_lite,
};
use crate::stream::{Endian, StreamInterface};
use crate::tiff::tags::olympus::{
    OLYMPUS_TAG_ASPECT_FRAME, OLYMPUS_TAG_CAMERA_SETTINGS, OLYMPUS_TAG_RAW_PROCESSING,
};
use crate::tiff::tags::panasonic::{
    PANASONIC_TAG_BOTTOM_BORDER, PANASONIC_TAG_ISO, PANASONIC_TAG_JPEG_IMAGE,
    PANASONIC_TAG_LEFT_BORDER, PANASONIC_TAG_RIGHT_BORDER, PANASONIC_TAG_TOP_BORDER,
};
use crate::tiff::tags::pentax::PENTAX_TAG_COLOR_SPACE;
use crate::tiff::tags::tiff::{
    EXIF_TAG_COLOR_SPACE, EXIF_TAG_DATE_TIME_ORIGINAL, EXIF_TAG_DEFAULT_CROP_SIZE,
    EXIF_TAG_EXIF_IFD, EXIF_TAG_EXPOSURE_TIME, EXIF_TAG_FNUMBER, EXIF_TAG_FOCAL_LENGTH,
    EXIF_TAG_GPS, EXIF_TAG_HEIGHT, EXIF_TAG_ISO_SPEED, EXIF_TAG_MAKERNOTES, EXIF_TAG_SUB_IFDS,
    EXIF_TAG_WIDTH, TIFF_TAG_CFA_PATTERN_DIM, TIFF_TAG_COMPRESSION, TIFF_TAG_DATE_TIME,
    TIFF_TAG_IMAGE_LENGTH, TIFF_TAG_IMAGE_WIDTH, TIFF_TAG_JPEG_BYTE_COUNT, TIFF_TAG_JPEG_OFFSET,
    TIFF_TAG_MAKE, TIFF_TAG_MODEL, TIFF_TAG_ORIENTATION, TIFF_TAG_PHOTOMETRIC,
    TIFF_TAG_STRIP_BYTE_COUNTS, TIFF_TAG_STRIP_OFFSETS, TIFF_TAG_SUB_FILE_TYPE,
};
use crate::tiff::{TiffContent, TiffError, TiffParser, TiffType};
pub use crate::types::{
    ColorSpace, DngInformation, Image, ImageFormat, PiexError, PreviewImageData, RawImageType,
};

/// Returns the minimum byte header length required by `is_raw`.
pub fn bytes_required_for_is_raw() -> usize {
    get_number_of_bytes_for_is_raw_lite()
}

/// Returns `true` if `stream` contains a recognized RAW file format header.
///
/// The `stream` position will be restored to its initial offset after inspecting the header.
pub fn is_raw<R: Read + Seek>(stream: &mut R) -> bool {
    let req_bytes = bytes_required_for_is_raw();
    let mut header = [0u8; 10_000];
    let to_read = req_bytes.min(header.len());
    let current_pos = match stream.stream_position() {
        Ok(pos) => pos,
        Err(_) => return false,
    };

    let read_bytes = match stream.read(&mut header[..to_read]) {
        Ok(n) => n,
        Err(_) => return false,
    };
    let _ = stream.seek(SeekFrom::Start(current_pos));

    if read_bytes == 0 {
        return false;
    }

    recognize_raw_image_type_lite(&header[..read_bytes]) != RawImageType::NonRaw
}

/// Returns a static slice of supported uppercase RAW file extension strings.
pub fn supported_extensions() -> &'static [&'static str] {
    &["ARW", "CR2", "CR3", "DNG", "NEF", "NRW", "ORF", "PEF", "RAF", "RW2", "SRW"]
}

/// Extracts the preview image and metadata into `preview_data`.
///
/// Returns `Ok(RawImageType)` on success, or `PiexError` on failure.
pub fn get_preview_image_data_into<S: StreamInterface + ?Sized>(
    stream: &S,
    preview_data: &mut PreviewImageData,
) -> Result<RawImageType, PiexError> {
    let req_bytes = bytes_required_for_is_raw();
    let header_len = req_bytes.min(stream.len());
    if header_len == 0 {
        return Err(PiexError::StreamReadFailed {
            offset: 0,
            expected: req_bytes,
            available: stream.len(),
        });
    }
    let mut header = [0u8; 10_000];
    let to_read = header_len.min(header.len());
    stream.get_data(0, &mut header[..to_read])?;

    let raw_type = recognize_raw_image_type_lite(&header[..to_read]);
    match raw_type {
        RawImageType::Arw => extract_arw(stream, preview_data)?,
        RawImageType::Cr2 => extract_cr2(stream, preview_data)?,
        RawImageType::Cr3 => extract_cr3(stream, preview_data)?,
        RawImageType::Dng => extract_dng(stream, preview_data)?,
        RawImageType::Nef | RawImageType::Nrw => extract_nef(stream, preview_data)?,
        RawImageType::Orf => extract_orf(stream, preview_data)?,
        RawImageType::Pef => extract_pef(stream, preview_data)?,
        RawImageType::Raf => extract_raf(stream, preview_data)?,
        RawImageType::Rw2 => extract_rw2(stream, preview_data)?,
        RawImageType::Srw => extract_srw(stream, preview_data)?,
        _ => return Err(PiexError::UnsupportedFormat),
    }

    Ok(raw_type)
}

/// Extracts the preview image and metadata from `stream`.
///
/// Returns `Ok((PreviewImageData, RawImageType))` on success, or `PiexError` on failure.
pub fn get_preview_image_data<S: StreamInterface + ?Sized>(
    stream: &S,
) -> Result<(PreviewImageData, RawImageType), PiexError> {
    let mut preview_data = PreviewImageData::default();
    let raw_type = get_preview_image_data_into(stream, &mut preview_data)?;
    Ok((preview_data, raw_type))
}

/// Retrieves the EXIF orientation for the image in `stream`.
///
/// Returns `Ok(orientation)` on success, or `PiexError` on failure.
pub fn get_orientation<S: StreamInterface + ?Sized>(stream: &S) -> Result<u32, PiexError> {
    // max(GetNumberOfBytesForIsOfType(kRafImage), GetNumberOfBytesForIsOfType(kCr3Image))
    let min_bytes = 12;
    let header_len = min_bytes.min(stream.len());
    if header_len == 0 {
        return Err(PiexError::StreamReadFailed {
            offset: 0,
            expected: min_bytes,
            available: stream.len(),
        });
    }
    let mut header = [0u8; 12];
    stream.get_data(0, &mut header[..header_len])?;

    if crate::format_recognizer::is_of_type(&header[..header_len], RawImageType::Raf) {
        raf_get_orientation(stream)
    } else if crate::format_recognizer::is_of_type(&header[..header_len], RawImageType::Cr3) {
        crate::cr3::extract_cr3_orientation(stream)
    } else {
        get_exif_orientation_internal(stream)
    }
}

fn get_exif_orientation_internal<S: StreamInterface + ?Sized>(
    stream: &S,
) -> Result<u32, PiexError> {
    let parser = match TiffParser::new(stream) {
        Ok(p) => p,
        Err(TiffError::StreamError(se)) => return Err(se),
        Err(e) => return Err(PiexError::CorruptData(e.to_string())),
    };
    const ORIENT_TAGS: [u32; 1] = [TIFF_TAG_ORIENTATION];
    let content = match parser.parse(&ORIENT_TAGS, 1) {
        Ok(c) => c,
        Err(TiffError::StreamError(se)) => return Err(se),
        Err(e) => return Err(PiexError::CorruptData(e.to_string())),
    };
    for dir in &content.tiff_directories {
        if dir.has(TIFF_TAG_ORIENTATION) {
            if let Ok(orient) = dir.get_u32(TIFF_TAG_ORIENTATION) {
                return Ok(orient);
            }
        }
    }
    Err(PiexError::CorruptData("TIFF_TAG_ORIENTATION not found in TIFF directories".into()))
}

const BASE_DESIRED_TAGS: [u32; 15] = [
    EXIF_TAG_COLOR_SPACE,
    EXIF_TAG_DATE_TIME_ORIGINAL,
    EXIF_TAG_EXPOSURE_TIME,
    EXIF_TAG_FNUMBER,
    EXIF_TAG_FOCAL_LENGTH,
    EXIF_TAG_GPS,
    EXIF_TAG_ISO_SPEED,
    TIFF_TAG_COMPRESSION,
    TIFF_TAG_DATE_TIME,
    EXIF_TAG_EXIF_IFD,
    TIFF_TAG_CFA_PATTERN_DIM,
    TIFF_TAG_MAKE,
    TIFF_TAG_MODEL,
    TIFF_TAG_ORIENTATION,
    TIFF_TAG_PHOTOMETRIC,
];

/// Internal helper function to parse TIFF metadata and fill preview image data.
fn get_preview_data<S: StreamInterface + ?Sized>(
    extended_tags: &[u32],
    tiff_offset: u32,
    number_of_ifds: u16,
    stream: &S,
    preview_data: &mut PreviewImageData,
) -> Result<TiffContent, PiexError> {
    let mut desired_tags = [0u32; 32];
    desired_tags[..BASE_DESIRED_TAGS.len()].copy_from_slice(&BASE_DESIRED_TAGS);
    let mut len = BASE_DESIRED_TAGS.len();
    for &tag in extended_tags {
        if len < desired_tags.len() && !desired_tags[..len].contains(&tag) {
            desired_tags[len] = tag;
            len += 1;
        }
    }

    let parser = TiffParser::with_offset(stream, tiff_offset)
        .map_err(|e| PiexError::CorruptData(e.to_string()))?;
    let content = parser
        .parse(&desired_tags[..len], number_of_ifds)
        .map_err(|e| PiexError::CorruptData(e.to_string()))?;

    if content.tiff_directories.is_empty() {
        return Err(PiexError::CorruptData("No TIFF directories found".into()));
    }

    parser
        .get_preview_image_data(&content, preview_data)
        .map_err(|e| PiexError::CorruptData(e.to_string()))?;

    Ok(content)
}

/// Queries TIFF tags (`TIFF_TAG_JPEG_BYTE_COUNT`, `TIFF_TAG_JPEG_OFFSET` plus `extended_tags`)
/// using `TiffParser` and sets `preview_data.thumbnail` offset and length if found.
fn get_thumbnail_offset_and_length<S: StreamInterface + ?Sized>(
    extended_tags: &[u32],
    stream: &S,
    preview_data: &mut PreviewImageData,
) -> Result<(), PiexError> {
    let mut desired_tags = [0u32; 16];
    desired_tags[0] = TIFF_TAG_JPEG_BYTE_COUNT;
    desired_tags[1] = TIFF_TAG_JPEG_OFFSET;
    let mut len = 2;
    for &tag in extended_tags {
        if len < desired_tags.len() && !desired_tags[..len].contains(&tag) {
            desired_tags[len] = tag;
            len += 1;
        }
    }

    let mut thumbnail_data = PreviewImageData::default();
    if get_preview_data(&desired_tags[..len], 0, 2, stream, &mut thumbnail_data).is_ok() {
        preview_data.thumbnail = thumbnail_data.thumbnail;
    }
    Ok(())
}

/// Parses Pentax PEF Makernote tag `0x0037`
/// (`0 = ColorSpace::Srgb`, non-zero = `ColorSpace::AdobeRgb`).
///
/// Returns `Ok(Some(ColorSpace))` if the tag is present and parsed,
/// `Ok(None)` if Makernote is parsed successfully but tag is absent,
/// or `Err(PiexError)` if TIFF, Exif IFD, Makernote offset, or Makernote IFD parsing fails.
fn pef_get_color_space<S: StreamInterface + ?Sized>(
    stream: &S,
) -> Result<Option<ColorSpace>, PiexError> {
    let parser = TiffParser::new(stream).map_err(|e| PiexError::CorruptData(e.to_string()))?;
    const PEF_EXIF_TAGS: [u32; 2] = [EXIF_TAG_EXIF_IFD, EXIF_TAG_MAKERNOTES];
    let content =
        parser.parse(&PEF_EXIF_TAGS, 1).map_err(|e| PiexError::CorruptData(e.to_string()))?;
    let exif_dir = content
        .exif_directory
        .as_ref()
        .ok_or_else(|| PiexError::CorruptData("Missing Exif IFD in PEF".into()))?;
    let (makernote_offset, _) = exif_dir
        .get_offset_and_length(EXIF_TAG_MAKERNOTES, TiffType::Undefined)
        .map_err(|e| PiexError::CorruptData(e.to_string()))?;

    let skip_offset = 6u32;
    const MN_TAGS: [u32; 1] = [PENTAX_TAG_COLOR_SPACE];
    let mut next_offset = 0u32;
    let mn_parser = TiffParser::with_offset_and_endian(stream, makernote_offset, parser.endian());
    let mn_ifd_offset = makernote_offset
        .checked_add(skip_offset)
        .ok_or_else(|| PiexError::CorruptData("Makernote offset overflow".into()))?;
    let mn_dir = mn_parser
        .parse_directory(mn_ifd_offset, &MN_TAGS, &mut next_offset)
        .map_err(|e| PiexError::CorruptData(e.to_string()))?;

    if mn_dir.has(PENTAX_TAG_COLOR_SPACE) {
        let color_space_val = mn_dir
            .get_u32(PENTAX_TAG_COLOR_SPACE)
            .map_err(|e| PiexError::CorruptData(e.to_string()))?;
        if color_space_val == 0 {
            Ok(Some(ColorSpace::Srgb))
        } else {
            Ok(Some(ColorSpace::AdobeRgb))
        }
    } else {
        Ok(None)
    }
}

/// Sony ARW extractor.
///
/// Queries thumbnail offset/length, sets up tag set (`EXIF_TAG_HEIGHT`, `EXIF_TAG_WIDTH`,
/// `TIFF_TAG_JPEG_BYTE_COUNT`, `TIFF_TAG_JPEG_OFFSET`, `EXIF_TAG_SUB_IFDS`),
/// and calls `TiffParser`.
pub fn extract_arw<S: StreamInterface + ?Sized>(
    stream: &S,
    preview_data: &mut PreviewImageData,
) -> Result<(), PiexError> {
    const ARW_TAGS: [u32; 5] = [
        EXIF_TAG_HEIGHT,
        EXIF_TAG_WIDTH,
        TIFF_TAG_JPEG_BYTE_COUNT,
        TIFF_TAG_JPEG_OFFSET,
        EXIF_TAG_SUB_IFDS,
    ];

    let _ = get_thumbnail_offset_and_length(&[], stream, preview_data);
    get_preview_data(&ARW_TAGS, 0, 1, stream, preview_data)?;
    Ok(())
}

/// Canon CR2 extractor.
///
/// Queries thumbnail data via IFD1, sets up strip tags (`EXIF_TAG_HEIGHT`, `EXIF_TAG_WIDTH`,
/// `TIFF_TAG_STRIP_BYTE_COUNTS`, `TIFF_TAG_STRIP_OFFSETS`), and calls `TiffParser`.
pub fn extract_cr2<S: StreamInterface + ?Sized>(
    stream: &S,
    preview_data: &mut PreviewImageData,
) -> Result<(), PiexError> {
    const CR2_TAGS: [u32; 4] =
        [EXIF_TAG_HEIGHT, EXIF_TAG_WIDTH, TIFF_TAG_STRIP_BYTE_COUNTS, TIFF_TAG_STRIP_OFFSETS];

    let _ = get_thumbnail_offset_and_length(&[], stream, preview_data);
    get_preview_data(&CR2_TAGS, 0, 1, stream, preview_data)?;
    Ok(())
}

/// Nikon NEF/NRW extractor.
///
/// Uses `TiffParser` to extract primary preview metadata.
/// Validates JPEG dimensions (`get_jpeg_dimensions`) against full sensor dimensions
/// (`EPSILON = 0.9` ratio test).
pub fn extract_nef<S: StreamInterface + ?Sized>(
    stream: &S,
    preview_data: &mut PreviewImageData,
) -> Result<(), PiexError> {
    const NEF_TAGS: [u32; 7] = [
        TIFF_TAG_IMAGE_WIDTH,
        TIFF_TAG_IMAGE_LENGTH,
        TIFF_TAG_JPEG_BYTE_COUNT,
        TIFF_TAG_JPEG_OFFSET,
        TIFF_TAG_STRIP_BYTE_COUNTS,
        TIFF_TAG_STRIP_OFFSETS,
        EXIF_TAG_SUB_IFDS,
    ];

    get_preview_data(&NEF_TAGS, 0, 2, stream, preview_data)?;

    if preview_data.thumbnail.length == 0 {
        let mut thumbnail_data = PreviewImageData::default();
        let _ = get_thumbnail_offset_and_length(&[], stream, &mut thumbnail_data);
        preview_data.thumbnail = thumbnail_data.thumbnail;
    }

    if preview_data.preview.length > 0 {
        const EPSILON: f32 = 0.9;
        let parser = TiffParser::new(stream).map_err(|e| PiexError::CorruptData(e.to_string()))?;
        let (width, height) = parser
            .get_jpeg_dimensions(preview_data.preview.offset)
            .map_err(|_| PiexError::UnsupportedFormat)?;

        if preview_data.full_width == 0 || preview_data.full_height == 0 {
            return Err(PiexError::UnsupportedFormat);
        }

        if (width as f32 / preview_data.full_width as f32) > EPSILON
            || (height as f32 / preview_data.full_height as f32) > EPSILON
        {
            preview_data.full_width = width as u32;
            preview_data.full_height = height as u32;
        }
    }

    Ok(())
}

/// Pentax PEF extractor.
///
/// Extracts preview with extended tags (`TIFF_TAG_IMAGE_WIDTH`, `TIFF_TAG_IMAGE_LENGTH`,
/// `TIFF_TAG_JPEG_BYTE_COUNT`, `TIFF_TAG_JPEG_OFFSET`, `EXIF_TAG_SUB_IFDS`),
/// sets `preview_data.color_space` via `pef_get_color_space`, and fetches thumbnail.
pub fn extract_pef<S: StreamInterface + ?Sized>(
    stream: &S,
    preview_data: &mut PreviewImageData,
) -> Result<(), PiexError> {
    const PEF_TAGS: [u32; 5] = [
        TIFF_TAG_IMAGE_WIDTH,
        TIFF_TAG_IMAGE_LENGTH,
        TIFF_TAG_JPEG_BYTE_COUNT,
        TIFF_TAG_JPEG_OFFSET,
        EXIF_TAG_SUB_IFDS,
    ];

    get_preview_data(&PEF_TAGS, 0, 3, stream, preview_data)?;

    if let Some(cs) = pef_get_color_space(stream)? {
        preview_data.color_space = cs;
    }

    let mut thumbnail_data = PreviewImageData::default();
    let _ = get_thumbnail_offset_and_length(&[], stream, &mut thumbnail_data);
    preview_data.thumbnail = thumbnail_data.thumbnail;

    Ok(())
}

/// Samsung SRW extractor.
///
/// Processes SubIFD containers for embedded preview JPEGs.
pub fn extract_srw<S: StreamInterface + ?Sized>(
    stream: &S,
    preview_data: &mut PreviewImageData,
) -> Result<(), PiexError> {
    const SUBIFD_TAGS: [u32; 1] = [EXIF_TAG_SUB_IFDS];
    let _ = get_thumbnail_offset_and_length(&SUBIFD_TAGS, stream, preview_data);

    const SRW_TAGS: [u32; 5] = [
        EXIF_TAG_WIDTH,
        EXIF_TAG_HEIGHT,
        TIFF_TAG_JPEG_BYTE_COUNT,
        TIFF_TAG_JPEG_OFFSET,
        EXIF_TAG_SUB_IFDS,
    ];

    get_preview_data(&SRW_TAGS, 0, 1, stream, preview_data)?;
    Ok(())
}

fn get_dng_info_internal<S: StreamInterface + ?Sized>(
    extended_tags: &[u32],
    stream: &S,
    width: &mut u32,
    height: &mut u32,
    cfa: &mut [u32; 2],
) -> bool {
    const DNG_INFO_BASE_TAGS: [u32; 4] = [
        EXIF_TAG_DEFAULT_CROP_SIZE,
        TIFF_TAG_CFA_PATTERN_DIM,
        EXIF_TAG_EXIF_IFD,
        TIFF_TAG_SUB_FILE_TYPE,
    ];
    let mut desired_tags = [0u32; 8];
    desired_tags[..DNG_INFO_BASE_TAGS.len()].copy_from_slice(&DNG_INFO_BASE_TAGS);
    let mut len = DNG_INFO_BASE_TAGS.len();
    for &tag in extended_tags {
        if len < desired_tags.len() && !desired_tags[..len].contains(&tag) {
            desired_tags[len] = tag;
            len += 1;
        }
    }

    let Ok(parser) = TiffParser::new(stream) else {
        return false;
    };
    let Ok(content) = parser.parse(&desired_tags[..len], 1) else {
        return false;
    };
    let Some(first_dir) = content.tiff_directories.first() else {
        return false;
    };

    let target_dir = if first_dir.sub_directories().is_empty() {
        first_dir
    } else {
        &first_dir.sub_directories()[0]
    };

    if let Ok((w, h)) = parser.get_full_dimension32(target_dir) {
        if w > 0 && h > 0 {
            *width = w;
            *height = h;
        } else {
            return false;
        }
    } else {
        return false;
    }

    if let Ok(cfa_vec) = target_dir.get_u32_vec(TIFF_TAG_CFA_PATTERN_DIM) {
        if cfa_vec.len() == 2 {
            *cfa = [cfa_vec[0], cfa_vec[1]];
            return true;
        }
    }

    false
}

/// Extracts full DNG width, height, and 2-element `cfa_pattern_dim` into out parameters.
pub fn get_dng_information_out<S: StreamInterface + ?Sized>(
    stream: &S,
    width: &mut u32,
    height: &mut u32,
    cfa: &mut [u32; 2],
) -> bool {
    if !get_dng_info_internal(&[], stream, width, height, cfa) {
        const SUBIFD_TAGS: [u32; 1] = [EXIF_TAG_SUB_IFDS];
        get_dng_info_internal(&SUBIFD_TAGS, stream, width, height, cfa)
    } else {
        true
    }
}

/// Extracts full DNG width, height, and 2-element `cfa_pattern_dim`.
/// First attempts with empty tags; if not found, checks SubIFD.
pub fn get_dng_information<S: StreamInterface + ?Sized>(stream: &S) -> Option<DngInformation> {
    let mut width = 0;
    let mut height = 0;
    let mut cfa = [0, 0];
    if get_dng_information_out(stream, &mut width, &mut height, &mut cfa) {
        Some(DngInformation { width, height, cfa_pattern_dim: cfa })
    } else {
        None
    }
}

/// Adobe DNG extractor logic.
///
/// Iterates primary IFD and SubIFDs, categorizes images <= 512px as thumbnails
/// (`THUMBNAIL_MAX_DIMENSION = 512`), and tracks the largest preview JPEG.
pub fn extract_dng<S: StreamInterface + ?Sized>(
    stream: &S,
    preview_data: &mut PreviewImageData,
) -> Result<(), PiexError> {
    const DNG_THUMBNAIL_MAX_DIMENSION: u32 = 512;

    const DNG_TAGS: [u32; 7] = [
        EXIF_TAG_DEFAULT_CROP_SIZE,
        TIFF_TAG_SUB_FILE_TYPE,
        TIFF_TAG_IMAGE_WIDTH,
        TIFF_TAG_IMAGE_LENGTH,
        TIFF_TAG_STRIP_BYTE_COUNTS,
        TIFF_TAG_STRIP_OFFSETS,
        EXIF_TAG_SUB_IFDS,
    ];

    let content = get_preview_data(&DNG_TAGS, 0, 3, stream, preview_data)?;
    let tiff_directory = content
        .tiff_directories
        .first()
        .ok_or_else(|| PiexError::CorruptData("No TIFF directory found".into()))?;

    if tiff_directory.has(EXIF_TAG_DEFAULT_CROP_SIZE) {
        if let Ok(crop) = tiff_directory.get_u32_vec(EXIF_TAG_DEFAULT_CROP_SIZE) {
            if crop.len() == 2 && crop[0] > 0 && crop[1] > 0 {
                preview_data.full_width = crop[0];
                preview_data.full_height = crop[1];
            } else {
                return Err(PiexError::CorruptData("Invalid default crop size".into()));
            }
        } else if let Ok(crop_rat) = tiff_directory.get_rational_vec(EXIF_TAG_DEFAULT_CROP_SIZE) {
            if crop_rat.len() == 2
                && crop_rat[0].numerator > 0
                && crop_rat[0].denominator > 0
                && crop_rat[1].numerator > 0
                && crop_rat[1].denominator > 0
            {
                preview_data.full_width = crop_rat[0].numerator / crop_rat[0].denominator;
                preview_data.full_height = crop_rat[1].numerator / crop_rat[1].denominator;
            } else {
                return Err(PiexError::CorruptData("Invalid default crop size rational".into()));
            }
        } else {
            return Err(PiexError::CorruptData("Invalid default crop size tag".into()));
        }
    }

    let parser = TiffParser::new(stream).map_err(|e| PiexError::CorruptData(e.to_string()))?;
    let mut preview = Image::default();
    let mut thumbnail = Image::default();

    let is_thumbnail_img = |img: &Image| -> bool {
        u32::from(img.width) <= DNG_THUMBNAIL_MAX_DIMENSION
            && u32::from(img.height) <= DNG_THUMBNAIL_MAX_DIMENSION
    };

    if let Ok(temp_image) = parser.get_image_data(tiff_directory) {
        if is_thumbnail_img(&temp_image) {
            thumbnail = temp_image;
        } else if temp_image.format == ImageFormat::JpegCompressed {
            preview = temp_image;
        }
    }

    for sub_ifd in tiff_directory.sub_directories() {
        if let Ok(temp_image) = parser.get_image_data(sub_ifd) {
            if is_thumbnail_img(&temp_image) {
                if temp_image.is_larger_than(&thumbnail) {
                    thumbnail = temp_image;
                }
            } else if temp_image.format == ImageFormat::JpegCompressed
                && temp_image.is_larger_than(&preview)
            {
                preview = temp_image;
            }
        }
    }

    preview_data.preview = preview;
    preview_data.thumbnail = thumbnail;

    Ok(())
}

/// Canon CR3 extractor.
///
/// Delegates preview extraction to `crate::cr3::extract_cr3_preview_data`.
pub fn extract_cr3<S: StreamInterface + ?Sized>(
    stream: &S,
    preview_data: &mut PreviewImageData,
) -> Result<(), PiexError> {
    crate::cr3::extract_cr3_preview_data_into(stream, preview_data)
}

/// Locates EXIF IFD -> MakerNote IFD -> CameraSettings IFD (`0x0100` thumbnail,
/// `0x0101` preview offset, `0x0102` preview length) -> RawProcessing IFD
/// (`OLYMPUS_TAG_ASPECT_FRAME = 0x1113`) to set `full_width` and `full_height`.
pub fn get_olympus_preview_image<S: StreamInterface + ?Sized>(
    stream: &S,
    preview_data: &mut PreviewImageData,
) -> Result<(), PiexError> {
    let parser = TiffParser::new(stream).map_err(|e| PiexError::CorruptData(e.to_string()))?;
    let mut header_bytes = [0u8; 4];
    stream.get_data(4, &mut header_bytes)?;
    let initial_ifd_offset = match parser.endian() {
        Endian::Little => u32::from_le_bytes(header_bytes),
        Endian::Big => u32::from_be_bytes(header_bytes),
    };

    const OLYMPUS_MAIN_TAGS: [u32; 1] = [EXIF_TAG_EXIF_IFD];
    let mut next_offset = 0u32;
    let tiff_ifd = parser
        .parse_directory(initial_ifd_offset, &OLYMPUS_MAIN_TAGS, &mut next_offset)
        .map_err(|e| PiexError::CorruptData(e.to_string()))?;

    let exif_offset =
        tiff_ifd.get_u32(EXIF_TAG_EXIF_IFD).map_err(|e| PiexError::CorruptData(e.to_string()))?;

    const MN_TAGS: [u32; 1] = [EXIF_TAG_MAKERNOTES];
    let exif_dir = parser
        .parse_directory(exif_offset, &MN_TAGS, &mut next_offset)
        .map_err(|e| PiexError::CorruptData(e.to_string()))?;

    let (makernote_offset, _) = exif_dir
        .get_offset_and_length(EXIF_TAG_MAKERNOTES, TiffType::Undefined)
        .map_err(|e| PiexError::CorruptData(e.to_string()))?;

    let mn_parser = TiffParser::with_offset_and_endian(stream, makernote_offset, parser.endian());
    const MN_DESIRED_TAGS: [u32; 3] = [
        0x0100, // thumbnail tag in Olympus makernote
        OLYMPUS_TAG_CAMERA_SETTINGS,
        OLYMPUS_TAG_RAW_PROCESSING,
    ];

    let makernote_ifd = mn_parser
        .parse_directory(
            makernote_offset
                .checked_add(12)
                .ok_or_else(|| PiexError::CorruptData("Makernote offset overflow".into()))?,
            &MN_DESIRED_TAGS,
            &mut next_offset,
        )
        .map_err(|e| PiexError::CorruptData(e.to_string()))?;

    if makernote_ifd.has(0x0100) {
        let (off, len) = makernote_ifd
            .get_offset_and_length(0x0100, TiffType::Undefined)
            .map_err(|e| PiexError::CorruptData(e.to_string()))?;
        preview_data.thumbnail.offset = off;
        preview_data.thumbnail.length = len;
    }

    let mut has_preview = false;
    if makernote_ifd.has(OLYMPUS_TAG_CAMERA_SETTINGS) {
        if let Ok((camera_settings_ptr, _)) =
            makernote_ifd.get_offset_and_length(OLYMPUS_TAG_CAMERA_SETTINGS, TiffType::Ifd)
        {
            let mut cs_ptr_bytes = [0u8; 4];
            if stream.get_data(camera_settings_ptr as usize, &mut cs_ptr_bytes).is_ok() {
                let camera_settings_offset = match parser.endian() {
                    Endian::Little => u32::from_le_bytes(cs_ptr_bytes),
                    Endian::Big => u32::from_be_bytes(cs_ptr_bytes),
                };
                const CS_DESIRED_TAGS: [u32; 2] = [0x0101, 0x0102];
                if let Ok(camera_settings_ifd) = mn_parser.parse_directory(
                    makernote_offset.saturating_add(camera_settings_offset),
                    &CS_DESIRED_TAGS,
                    &mut next_offset,
                ) {
                    if camera_settings_ifd.has(0x0101) && camera_settings_ifd.has(0x0102) {
                        has_preview = true;
                        if let Ok(p_off) = camera_settings_ifd.get_u32(0x0101) {
                            preview_data.preview.offset = p_off;
                        }
                        preview_data.preview.offset =
                            preview_data.preview.offset.saturating_add(makernote_offset);
                        if let Ok(p_len) = camera_settings_ifd.get_u32(0x0102) {
                            preview_data.preview.length = p_len;
                        }
                    }
                }
            }
        }
    }
    if !has_preview {
        return Err(PiexError::CorruptData("No preview in Olympus camera settings".into()));
    }

    if makernote_ifd.has(OLYMPUS_TAG_RAW_PROCESSING) {
        if let Ok((raw_proc_ptr, _)) =
            makernote_ifd.get_offset_and_length(OLYMPUS_TAG_RAW_PROCESSING, TiffType::Ifd)
        {
            let mut rp_ptr_bytes = [0u8; 4];
            if stream.get_data(raw_proc_ptr as usize, &mut rp_ptr_bytes).is_ok() {
                let raw_proc_offset = match parser.endian() {
                    Endian::Little => u32::from_le_bytes(rp_ptr_bytes),
                    Endian::Big => u32::from_be_bytes(rp_ptr_bytes),
                };
                const RP_DESIRED_TAGS: [u32; 1] = [OLYMPUS_TAG_ASPECT_FRAME];
                if let Ok(raw_processing_ifd) = mn_parser.parse_directory(
                    makernote_offset.saturating_add(raw_proc_offset),
                    &RP_DESIRED_TAGS,
                    &mut next_offset,
                ) {
                    if raw_processing_ifd.has(OLYMPUS_TAG_ASPECT_FRAME) {
                        if let Ok(aspect_frame) =
                            raw_processing_ifd.get_u32_vec(OLYMPUS_TAG_ASPECT_FRAME)
                        {
                            if aspect_frame.len() == 4
                                && aspect_frame[2] > aspect_frame[0]
                                && aspect_frame[3] > aspect_frame[1]
                            {
                                preview_data.full_width = aspect_frame[2] - aspect_frame[0] + 1;
                                preview_data.full_height = aspect_frame[3] - aspect_frame[1] + 1;
                                if preview_data.full_width < preview_data.full_height {
                                    std::mem::swap(
                                        &mut preview_data.full_width,
                                        &mut preview_data.full_height,
                                    );
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    Ok(())
}

/// Olympus ORF extractor.
///
/// Extracts primary EXIF preview metadata and calls `get_olympus_preview_image`.
pub fn extract_orf<S: StreamInterface + ?Sized>(
    stream: &S,
    preview_data: &mut PreviewImageData,
) -> Result<(), PiexError> {
    const ORF_EXTENDED_TAGS: [u32; 2] = [TIFF_TAG_JPEG_BYTE_COUNT, TIFF_TAG_JPEG_OFFSET];
    get_preview_data(&ORF_EXTENDED_TAGS, 0, 2, stream, preview_data)?;
    let _ = get_olympus_preview_image(stream, preview_data);
    Ok(())
}

/// Reads 32-bit big-endian preview offset at offset 84, reads EXIF orientation tag
/// from `preview_offset + 12`.
pub fn raf_get_orientation<S: StreamInterface + ?Sized>(stream: &S) -> Result<u32, PiexError> {
    let mut buf = [0u8; 4];
    stream.get_data(84, &mut buf)?;
    let preview_offset = u32::from_be_bytes(buf);
    let exif_offset = preview_offset.checked_add(12).ok_or_else(|| PiexError::OutOfBounds {
        requested_offset: preview_offset as usize,
        requested_size: 12,
        total_len: stream.len(),
    })?;

    let parser = TiffParser::with_offset(stream, exif_offset)
        .map_err(|e| PiexError::CorruptData(e.to_string()))?;
    const RAF_TAGS: [u32; 1] = [TIFF_TAG_ORIENTATION];
    let content = parser.parse(&RAF_TAGS, 1).map_err(|e| PiexError::CorruptData(e.to_string()))?;

    for dir in &content.tiff_directories {
        if dir.has(TIFF_TAG_ORIENTATION) {
            if let Ok(orient) = dir.get_u32(TIFF_TAG_ORIENTATION) {
                return Ok(orient);
            }
        }
    }
    Err(PiexError::CorruptData("TIFF_TAG_ORIENTATION not found in RAF EXIF".into()))
}

/// Reads CFA header table offset at byte position 92, iterates entries for tag `0x0111`
/// to extract width and height.
pub fn raf_get_dimension<S: StreamInterface + ?Sized>(stream: &S) -> Result<(u32, u32), PiexError> {
    let mut buf4 = [0u8; 4];
    stream.get_data(92, &mut buf4)?;
    let cfa_header_offset = u32::from_be_bytes(buf4);

    stream.get_data(cfa_header_offset as usize, &mut buf4)?;
    let cfa_header_entries = u32::from_be_bytes(buf4);

    let mut index = cfa_header_offset
        .checked_add(4)
        .ok_or_else(|| PiexError::CorruptData("CFA header offset overflow".into()))?;

    for _ in 0..cfa_header_entries {
        let mut entry_buf = [0u8; 8];
        stream.get_data(index as usize, &mut entry_buf)?;

        let id = u16::from_be_bytes([entry_buf[0], entry_buf[1]]);
        let length = u16::from_be_bytes([entry_buf[2], entry_buf[3]]);

        if id == 0x0111 {
            let tmp_height = u16::from_be_bytes([entry_buf[4], entry_buf[5]]);
            let tmp_width = u16::from_be_bytes([entry_buf[6], entry_buf[7]]);
            return Ok((u32::from(tmp_width), u32::from(tmp_height)));
        }

        index = index
            .checked_add(4 + u32::from(length))
            .ok_or_else(|| PiexError::CorruptData("CFA entry index overflow".into()))?;
    }

    Err(PiexError::CorruptData("CFA crop dimensions tag 0x0111 not found".into()))
}

/// Fuji RAF extractor.
///
/// Reads RAF header preview offset/length, computes sensor crop dimensions,
/// and extracts EXIF metadata.
pub fn extract_raf<S: StreamInterface + ?Sized>(
    stream: &S,
    preview_data: &mut PreviewImageData,
) -> Result<(), PiexError> {
    let mut buf = [0u8; 8];
    stream.get_data(84, &mut buf)?;
    let preview_offset = u32::from_be_bytes([buf[0], buf[1], buf[2], buf[3]]);
    let preview_length = u32::from_be_bytes([buf[4], buf[5], buf[6], buf[7]]);

    let (fw, fh) = raf_get_dimension(stream)?;
    preview_data.full_width = fw;
    preview_data.full_height = fh;

    if preview_length > 0 {
        let exif_offset = preview_offset
            .checked_add(12)
            .ok_or_else(|| PiexError::CorruptData("RAF EXIF offset overflow".into()))?;
        const RAF_EXTENDED_TAGS: [u32; 2] = [TIFF_TAG_JPEG_BYTE_COUNT, TIFF_TAG_JPEG_OFFSET];
        get_preview_data(&RAF_EXTENDED_TAGS, exif_offset, 2, stream, preview_data)?;
    }

    preview_data.thumbnail.offset = preview_data.thumbnail.offset.saturating_add(160);
    preview_data.preview.offset = preview_offset;
    preview_data.preview.length = preview_length;

    Ok(())
}

/// Panasonic RW2 extractor.
///
/// Parses Panasonic border tags (`PANASONIC_TAG_TOP_BORDER`, `LEFT_BORDER`,
/// `BOTTOM_BORDER`, `RIGHT_BORDER`), `PANASONIC_TAG_ISO`, `PANASONIC_TAG_JPEG_IMAGE`,
/// and relative EXIF thumbnail offsets (+12 byte shift).
pub fn extract_rw2<S: StreamInterface + ?Sized>(
    stream: &S,
    preview_data: &mut PreviewImageData,
) -> Result<(), PiexError> {
    const RW2_EXTENDED_TAGS: [u32; 8] = [
        PANASONIC_TAG_TOP_BORDER,
        PANASONIC_TAG_LEFT_BORDER,
        PANASONIC_TAG_BOTTOM_BORDER,
        PANASONIC_TAG_RIGHT_BORDER,
        PANASONIC_TAG_ISO,
        PANASONIC_TAG_JPEG_IMAGE,
        TIFF_TAG_JPEG_BYTE_COUNT,
        TIFF_TAG_JPEG_OFFSET,
    ];

    let mut raw_preview_data = PreviewImageData::default();
    get_preview_data(&RW2_EXTENDED_TAGS, 0, 1, stream, &mut raw_preview_data)?;

    if raw_preview_data.preview.length > 0 {
        let exif_offset = raw_preview_data
            .preview
            .offset
            .checked_add(12)
            .ok_or_else(|| PiexError::CorruptData("RW2 EXIF offset overflow".into()))?;
        const EXIF_TAGS: [u32; 2] = [TIFF_TAG_JPEG_BYTE_COUNT, TIFF_TAG_JPEG_OFFSET];
        get_preview_data(&EXIF_TAGS, exif_offset, 2, stream, preview_data)?;
        preview_data.thumbnail.offset = preview_data.thumbnail.offset.saturating_add(exif_offset);
    }

    preview_data.preview = raw_preview_data.preview;
    preview_data.iso = raw_preview_data.iso;
    preview_data.full_width = raw_preview_data.full_width;
    preview_data.full_height = raw_preview_data.full_height;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_bytes_required_for_is_raw() {
        assert_eq!(bytes_required_for_is_raw(), 10_000);
    }

    #[test]
    fn test_is_raw_with_empty_stream() {
        let mut cursor = Cursor::new(vec![]);
        assert!(!is_raw(&mut cursor));
        assert_eq!(cursor.position(), 0);
    }

    #[test]
    fn test_is_raw_with_small_stream() {
        let mut cursor_non_raw = Cursor::new(b"NON_RAW!".to_vec());
        assert!(!is_raw(&mut cursor_non_raw));
        assert_eq!(cursor_non_raw.position(), 0);

        let mut cursor_raf = Cursor::new(b"FUJIFILM".to_vec());
        assert!(is_raw(&mut cursor_raf));
        assert_eq!(cursor_raf.position(), 0);
    }

    #[test]
    fn test_is_raw_with_valid_raf() {
        let mut header = vec![0u8; 10_000];
        header[..8].copy_from_slice(b"FUJIFILM");
        let mut cursor = Cursor::new(header);
        cursor.set_position(0);

        assert!(is_raw(&mut cursor));
        assert_eq!(cursor.position(), 0);
    }

    #[test]
    fn test_is_raw_restores_stream_position() {
        let mut header = vec![0u8; 12_000];
        header[1000..1008].copy_from_slice(b"FUJIFILM");
        let mut cursor = Cursor::new(header);
        cursor.set_position(1000);

        assert!(is_raw(&mut cursor));
        assert_eq!(cursor.position(), 1000);
    }

    #[test]
    fn test_is_raw_with_non_raw() {
        let mut cursor = Cursor::new(vec![0u8; 10_000]);
        assert!(!is_raw(&mut cursor));
        assert_eq!(cursor.position(), 0);
    }

    #[test]
    fn test_supported_extensions() {
        let exts = supported_extensions();
        assert_eq!(exts.len(), 11);
        assert!(exts.contains(&"ARW"));
        assert!(exts.contains(&"CR2"));
        assert!(exts.contains(&"CR3"));
        assert!(exts.contains(&"DNG"));
        assert!(exts.contains(&"NEF"));
        assert!(exts.contains(&"NRW"));
        assert!(exts.contains(&"ORF"));
        assert!(exts.contains(&"PEF"));
        assert!(exts.contains(&"RAF"));
        assert!(exts.contains(&"RW2"));
        assert!(exts.contains(&"SRW"));

        // All uppercase requirement check
        assert!(exts.iter().all(|ext| ext.to_uppercase() == *ext));
    }

    #[test]
    fn test_reexported_types() {
        let dng = DngInformation::default();
        assert_eq!(dng.width, 0);

        let preview_data = PreviewImageData::default();
        assert_eq!(preview_data.full_width, 0);

        let err = PiexError::GenericFail("test".into());
        assert_eq!(err, PiexError::GenericFail("test".into()));
    }
}
