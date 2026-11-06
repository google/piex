//! C-compatible ABI type definitions and conversion traits for PIEX.

use crate::ffi_ext::{c_buf_to_bytes, copy_bytes_to_c_buf};
use crate::types::{
    ColorSpace, GpsData, Image, ImageFormat, PiexError, PreviewImageData, Rational, RawImageType,
};
use core::ffi::c_char;

/// C-compatible status/error code return values for PIEX FFI API functions.
#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
#[allow(non_camel_case_types)]
pub enum piex_error_t {
    /// Operation succeeded without error.
    #[default]
    PIEX_OK = 0,
    /// General or stream operation failure.
    PIEX_FAIL = 1,
    /// Format or feature is unsupported.
    PIEX_UNSUPPORTED = 2,
    /// Invalid argument passed to FFI function.
    PIEX_INVALID_ARGUMENT = 3,
    /// Panic caught at FFI boundary.
    PIEX_PANIC = 4,
}

impl From<PiexError> for piex_error_t {
    fn from(err: PiexError) -> Self {
        match err {
            PiexError::UnsupportedFormat => piex_error_t::PIEX_UNSUPPORTED,
            PiexError::OutOfBounds { .. } => piex_error_t::PIEX_FAIL,
            PiexError::StreamReadFailed { .. } => piex_error_t::PIEX_FAIL,
            PiexError::InvalidHeader(_) => piex_error_t::PIEX_FAIL,
            PiexError::InvalidBoxHeader { .. } => piex_error_t::PIEX_FAIL,
            PiexError::CorruptData(_) => piex_error_t::PIEX_FAIL,
            PiexError::IoError(_) => piex_error_t::PIEX_FAIL,
            PiexError::GenericFail(_) => piex_error_t::PIEX_FAIL,
        }
    }
}

impl From<&PiexError> for piex_error_t {
    fn from(err: &PiexError) -> Self {
        match err {
            PiexError::UnsupportedFormat => piex_error_t::PIEX_UNSUPPORTED,
            _ => piex_error_t::PIEX_FAIL,
        }
    }
}

/// C-compatible enum representing RAW image camera formats.
#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
#[allow(non_camel_case_types)]
pub enum piex_raw_image_type_t {
    /// Non-RAW or unknown image format.
    #[default]
    PIEX_RAW_TYPE_NON_RAW = 0,
    /// Sony ARW format.
    PIEX_RAW_TYPE_ARW = 1,
    /// Canon CR2 format.
    PIEX_RAW_TYPE_CR2 = 2,
    /// Canon CR3 format.
    PIEX_RAW_TYPE_CR3 = 3,
    /// Canon CRW format.
    PIEX_RAW_TYPE_CRW = 4,
    /// Kodak DCR format.
    PIEX_RAW_TYPE_DCR = 5,
    /// Digital Negative (DNG) format.
    PIEX_RAW_TYPE_DNG = 6,
    /// Epson ERF format.
    PIEX_RAW_TYPE_ERF = 7,
    /// Kodak KDC format.
    PIEX_RAW_TYPE_KDC = 8,
    /// Leaf MOS format.
    PIEX_RAW_TYPE_MOS = 9,
    /// Minolta MRW format.
    PIEX_RAW_TYPE_MRW = 10,
    /// Nikon NEF format.
    PIEX_RAW_TYPE_NEF = 11,
    /// Nikon NRW format.
    PIEX_RAW_TYPE_NRW = 12,
    /// Olympus ORF format.
    PIEX_RAW_TYPE_ORF = 13,
    /// Pentax PEF format.
    PIEX_RAW_TYPE_PEF = 14,
    /// Phase One RAW format.
    PIEX_RAW_TYPE_PHASE_ONE = 15,
    /// Logitech Fotoman PXN format.
    PIEX_RAW_TYPE_PXN = 16,
    /// Apple QTK format.
    PIEX_RAW_TYPE_QTK = 17,
    /// Fuji RAF format.
    PIEX_RAW_TYPE_RAF = 18,
    /// Hasselblad 3FR format.
    PIEX_RAW_TYPE_RAW_3FR = 19,
    /// Contax N RAW format.
    PIEX_RAW_TYPE_RAW_CONTAX_N = 20,
    /// Rollei RDC RAW format.
    PIEX_RAW_TYPE_RDC = 21,
    /// Panasonic RW2 format.
    PIEX_RAW_TYPE_RW2 = 22,
    /// Samsung SRW format.
    PIEX_RAW_TYPE_SRW = 23,
    /// Sigma / Polaroid X3F format.
    PIEX_RAW_TYPE_X3F = 24,
}

impl From<RawImageType> for piex_raw_image_type_t {
    fn from(ty: RawImageType) -> Self {
        match ty {
            RawImageType::NonRaw => piex_raw_image_type_t::PIEX_RAW_TYPE_NON_RAW,
            RawImageType::Arw => piex_raw_image_type_t::PIEX_RAW_TYPE_ARW,
            RawImageType::Cr2 => piex_raw_image_type_t::PIEX_RAW_TYPE_CR2,
            RawImageType::Cr3 => piex_raw_image_type_t::PIEX_RAW_TYPE_CR3,
            RawImageType::Crw => piex_raw_image_type_t::PIEX_RAW_TYPE_CRW,
            RawImageType::Dcr => piex_raw_image_type_t::PIEX_RAW_TYPE_DCR,
            RawImageType::Dng => piex_raw_image_type_t::PIEX_RAW_TYPE_DNG,
            RawImageType::Erf => piex_raw_image_type_t::PIEX_RAW_TYPE_ERF,
            RawImageType::Kdc => piex_raw_image_type_t::PIEX_RAW_TYPE_KDC,
            RawImageType::Mos => piex_raw_image_type_t::PIEX_RAW_TYPE_MOS,
            RawImageType::Mrw => piex_raw_image_type_t::PIEX_RAW_TYPE_MRW,
            RawImageType::Nef => piex_raw_image_type_t::PIEX_RAW_TYPE_NEF,
            RawImageType::Nrw => piex_raw_image_type_t::PIEX_RAW_TYPE_NRW,
            RawImageType::Orf => piex_raw_image_type_t::PIEX_RAW_TYPE_ORF,
            RawImageType::Pef => piex_raw_image_type_t::PIEX_RAW_TYPE_PEF,
            RawImageType::PhaseOne => piex_raw_image_type_t::PIEX_RAW_TYPE_PHASE_ONE,
            RawImageType::Pxn => piex_raw_image_type_t::PIEX_RAW_TYPE_PXN,
            RawImageType::Qtk => piex_raw_image_type_t::PIEX_RAW_TYPE_QTK,
            RawImageType::Raf => piex_raw_image_type_t::PIEX_RAW_TYPE_RAF,
            RawImageType::Raw3fr => piex_raw_image_type_t::PIEX_RAW_TYPE_RAW_3FR,
            RawImageType::RawContaxN => piex_raw_image_type_t::PIEX_RAW_TYPE_RAW_CONTAX_N,
            RawImageType::Rdc => piex_raw_image_type_t::PIEX_RAW_TYPE_RDC,
            RawImageType::Rw2 => piex_raw_image_type_t::PIEX_RAW_TYPE_RW2,
            RawImageType::Srw => piex_raw_image_type_t::PIEX_RAW_TYPE_SRW,
            RawImageType::X3f => piex_raw_image_type_t::PIEX_RAW_TYPE_X3F,
        }
    }
}

impl From<piex_raw_image_type_t> for RawImageType {
    fn from(ty: piex_raw_image_type_t) -> Self {
        RawImageType::from_u32(ty as u32).unwrap_or(RawImageType::NonRaw)
    }
}

/// C-compatible enum for preview image data formatting options.
#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
#[allow(non_camel_case_types)]
pub enum piex_image_format_t {
    /// JPEG compressed image format.
    #[default]
    PIEX_FORMAT_JPEG_COMPRESSED = 0,
    /// Uncompressed RGB image format.
    PIEX_FORMAT_UNCOMPRESSED_RGB = 1,
    /// HEVC compressed image format.
    PIEX_FORMAT_HEVC_COMPRESSED = 2,
}

impl From<ImageFormat> for piex_image_format_t {
    fn from(format: ImageFormat) -> Self {
        match format {
            ImageFormat::JpegCompressed => piex_image_format_t::PIEX_FORMAT_JPEG_COMPRESSED,
            ImageFormat::UncompressedRgb => piex_image_format_t::PIEX_FORMAT_UNCOMPRESSED_RGB,
            ImageFormat::HevcCompressed => piex_image_format_t::PIEX_FORMAT_HEVC_COMPRESSED,
        }
    }
}

impl From<piex_image_format_t> for ImageFormat {
    fn from(format: piex_image_format_t) -> Self {
        match format {
            piex_image_format_t::PIEX_FORMAT_JPEG_COMPRESSED => ImageFormat::JpegCompressed,
            piex_image_format_t::PIEX_FORMAT_UNCOMPRESSED_RGB => ImageFormat::UncompressedRgb,
            piex_image_format_t::PIEX_FORMAT_HEVC_COMPRESSED => ImageFormat::HevcCompressed,
        }
    }
}

/// C-compatible enum for preview image color space.
#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
#[allow(non_camel_case_types)]
pub enum piex_color_space_t {
    /// sRGB color space.
    #[default]
    PIEX_COLOR_SPACE_SRGB = 0,
    /// Adobe RGB color space.
    PIEX_COLOR_SPACE_ADOBE_RGB = 1,
}

impl From<ColorSpace> for piex_color_space_t {
    fn from(cs: ColorSpace) -> Self {
        match cs {
            ColorSpace::Srgb => piex_color_space_t::PIEX_COLOR_SPACE_SRGB,
            ColorSpace::AdobeRgb => piex_color_space_t::PIEX_COLOR_SPACE_ADOBE_RGB,
        }
    }
}

impl From<piex_color_space_t> for ColorSpace {
    fn from(cs: piex_color_space_t) -> Self {
        match cs {
            piex_color_space_t::PIEX_COLOR_SPACE_SRGB => ColorSpace::Srgb,
            piex_color_space_t::PIEX_COLOR_SPACE_ADOBE_RGB => ColorSpace::AdobeRgb,
        }
    }
}

/// C-compatible preview image descriptor structure.
#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[allow(non_camel_case_types)]
pub struct piex_image_t {
    /// Width of the image in pixels.
    pub width: u16,
    /// Height of the image in pixels.
    pub height: u16,
    /// Length of image data in bytes.
    pub length: u32,
    /// Offset of image data from the beginning of the stream.
    pub offset: u32,
    /// Image encoding/compression format.
    pub format: piex_image_format_t,
}

impl From<Image> for piex_image_t {
    fn from(img: Image) -> Self {
        Self {
            width: img.width,
            height: img.height,
            length: img.length,
            offset: img.offset,
            format: img.format.into(),
        }
    }
}

impl From<piex_image_t> for Image {
    fn from(img: piex_image_t) -> Self {
        Self {
            width: img.width,
            height: img.height,
            length: img.length,
            offset: img.offset,
            format: img.format.into(),
        }
    }
}

/// C-compatible unsigned rational number structure.
#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[allow(non_camel_case_types)]
pub struct piex_rational_t {
    /// Numerator value.
    pub numerator: u32,
    /// Denominator value.
    pub denominator: u32,
}

impl From<Rational> for piex_rational_t {
    fn from(r: Rational) -> Self {
        Self { numerator: r.numerator, denominator: r.denominator }
    }
}

impl From<piex_rational_t> for Rational {
    fn from(r: piex_rational_t) -> Self {
        Self { numerator: r.numerator, denominator: r.denominator }
    }
}

/// C-compatible GPS metadata descriptor structure.
#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[allow(non_camel_case_types)]
pub struct piex_gps_t {
    /// Indicates whether GPS data is valid.
    pub is_valid: bool,
    /// Latitude reference ('N' or 'S').
    pub latitude_ref: c_char,
    /// Latitude degrees, minutes, and seconds.
    pub latitude: [piex_rational_t; 3],
    /// Longitude reference ('E' or 'W').
    pub longitude_ref: c_char,
    /// Longitude degrees, minutes, and seconds.
    pub longitude: [piex_rational_t; 3],
    /// Altitude reference (true for above sea level, false for below).
    pub altitude_ref: bool,
    /// Altitude rational value.
    pub altitude: piex_rational_t,
    /// GPS timestamp (hours, minutes, seconds).
    pub time_stamp: [piex_rational_t; 3],
    /// GPS date stamp string buffer (e.g. "YYYY:MM:DD\0").
    pub date_stamp: [c_char; 12],
}

impl Default for piex_gps_t {
    fn default() -> Self {
        Self {
            is_valid: false,
            latitude_ref: b'N' as c_char,
            latitude: [piex_rational_t::default(); 3],
            longitude_ref: b'E' as c_char,
            longitude: [piex_rational_t::default(); 3],
            altitude_ref: false,
            altitude: piex_rational_t::default(),
            time_stamp: [piex_rational_t::default(); 3],
            date_stamp: [0; 12],
        }
    }
}

impl From<&GpsData> for piex_gps_t {
    fn from(gps: &GpsData) -> Self {
        let mut date_stamp = [0; 12];
        copy_bytes_to_c_buf(&gps.date_stamp, &mut date_stamp);

        Self {
            is_valid: gps.is_valid,
            latitude_ref: gps.latitude_ref as c_char,
            latitude: gps.latitude.map(piex_rational_t::from),
            longitude_ref: gps.longitude_ref as c_char,
            longitude: gps.longitude.map(piex_rational_t::from),
            altitude_ref: gps.altitude_ref,
            altitude: piex_rational_t::from(gps.altitude),
            time_stamp: gps.time_stamp.map(piex_rational_t::from),
            date_stamp,
        }
    }
}

impl From<GpsData> for piex_gps_t {
    fn from(gps: GpsData) -> Self {
        piex_gps_t::from(&gps)
    }
}

impl From<&piex_gps_t> for GpsData {
    fn from(gps: &piex_gps_t) -> Self {
        Self {
            is_valid: gps.is_valid,
            latitude_ref: gps.latitude_ref as u8,
            latitude: gps.latitude.map(Rational::from),
            longitude_ref: gps.longitude_ref as u8,
            longitude: gps.longitude.map(Rational::from),
            altitude_ref: gps.altitude_ref,
            altitude: Rational::from(gps.altitude),
            time_stamp: gps.time_stamp.map(Rational::from),
            date_stamp: c_buf_to_bytes(&gps.date_stamp),
        }
    }
}

impl From<piex_gps_t> for GpsData {
    fn from(gps: piex_gps_t) -> Self {
        GpsData::from(&gps)
    }
}

/// C-compatible preview image data and EXIF metadata container structure.
#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[allow(non_camel_case_types)]
pub struct piex_preview_image_data_t {
    /// Preview image descriptor.
    pub preview: piex_image_t,
    /// Thumbnail image descriptor.
    pub thumbnail: piex_image_t,
    /// EXIF orientation value.
    pub exif_orientation: u32,
    /// Color space of the preview.
    pub color_space: piex_color_space_t,
    /// Full RAW image width in pixels.
    pub full_width: u32,
    /// Full RAW image height in pixels.
    pub full_height: u32,
    /// Camera manufacturer string buffer (NUL-terminated).
    pub maker: [c_char; 64],
    /// Camera model string buffer (NUL-terminated).
    pub model: [c_char; 64],
    /// Date and time string buffer (NUL-terminated).
    pub date_time: [c_char; 32],
    /// ISO speed rating.
    pub iso: u32,
    /// Exposure time rational value.
    pub exposure_time: piex_rational_t,
    /// Aperture F-number rational value.
    pub fnumber: piex_rational_t,
    /// Focal length rational value in millimeters.
    pub focal_length: piex_rational_t,
    /// GPS location metadata structure.
    pub gps: piex_gps_t,
    /// Color Filter Array (CFA) pattern dimensions [width, height].
    pub cfa_pattern_dim: [u32; 2],
}

impl Default for piex_preview_image_data_t {
    fn default() -> Self {
        Self {
            preview: piex_image_t::default(),
            thumbnail: piex_image_t::default(),
            exif_orientation: 1,
            color_space: piex_color_space_t::default(),
            full_width: 0,
            full_height: 0,
            maker: [0; 64],
            model: [0; 64],
            date_time: [0; 32],
            iso: 0,
            exposure_time: piex_rational_t::default(),
            fnumber: piex_rational_t::default(),
            focal_length: piex_rational_t::default(),
            gps: piex_gps_t::default(),
            cfa_pattern_dim: [0, 0],
        }
    }
}

impl From<&PreviewImageData> for piex_preview_image_data_t {
    fn from(data: &PreviewImageData) -> Self {
        let mut maker = [0; 64];
        let mut model = [0; 64];
        let mut date_time = [0; 32];

        copy_bytes_to_c_buf(&data.maker, &mut maker);
        copy_bytes_to_c_buf(&data.model, &mut model);
        copy_bytes_to_c_buf(&data.date_time, &mut date_time);

        Self {
            preview: piex_image_t::from(data.preview),
            thumbnail: piex_image_t::from(data.thumbnail),
            exif_orientation: data.exif_orientation,
            color_space: piex_color_space_t::from(data.color_space),
            full_width: data.full_width,
            full_height: data.full_height,
            maker,
            model,
            date_time,
            iso: data.iso,
            exposure_time: piex_rational_t::from(data.exposure_time),
            fnumber: piex_rational_t::from(data.fnumber),
            focal_length: piex_rational_t::from(data.focal_length),
            gps: piex_gps_t::from(&data.gps),
            cfa_pattern_dim: data.cfa_pattern_dim,
        }
    }
}

impl From<PreviewImageData> for piex_preview_image_data_t {
    fn from(data: PreviewImageData) -> Self {
        piex_preview_image_data_t::from(&data)
    }
}

impl From<&piex_preview_image_data_t> for PreviewImageData {
    fn from(data: &piex_preview_image_data_t) -> Self {
        Self {
            preview: Image::from(data.preview),
            thumbnail: Image::from(data.thumbnail),
            exif_orientation: data.exif_orientation,
            color_space: ColorSpace::from(data.color_space),
            full_width: data.full_width,
            full_height: data.full_height,
            maker: c_buf_to_bytes(&data.maker),
            model: c_buf_to_bytes(&data.model),
            date_time: c_buf_to_bytes(&data.date_time),
            iso: data.iso,
            exposure_time: Rational::from(data.exposure_time),
            fnumber: Rational::from(data.fnumber),
            focal_length: Rational::from(data.focal_length),
            gps: GpsData::from(&data.gps),
            cfa_pattern_dim: data.cfa_pattern_dim,
            cfa_pattern_dim_set: false,
        }
    }
}

impl From<piex_preview_image_data_t> for PreviewImageData {
    fn from(data: piex_preview_image_data_t) -> Self {
        PreviewImageData::from(&data)
    }
}

/// Callback function pointer type for reading raw stream bytes from C callers.
///
/// Returns 0 (`PIEX_OK`) on success, or non-zero on failure.
#[allow(non_camel_case_types)]
pub type piex_read_fn_t = Option<
    unsafe extern "C" fn(
        user_data: *mut core::ffi::c_void,
        offset: usize,
        length: usize,
        out_data: *mut u8,
    ) -> i32,
>;

/// C-compatible stream handle containing read callback pointer and opaque user context pointer.
#[repr(C)]
#[derive(Debug, Clone, Copy)]
#[allow(non_camel_case_types)]
pub struct piex_stream_t {
    /// C function pointer invoked to read stream data bytes.
    pub read_fn: piex_read_fn_t,
    /// Opaque caller-provided user context pointer.
    pub user_data: *mut core::ffi::c_void,
}

impl Default for piex_stream_t {
    fn default() -> Self {
        Self { read_fn: None, user_data: core::ptr::null_mut() }
    }
}
