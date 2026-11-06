#![forbid(unsafe_code)]

//! Types and error definitions for the PIEX RAW preview extraction library.

use thiserror::Error;

/// Error types returned by PIEX operations.
#[derive(Debug, Error, PartialEq, Eq, Clone)]
pub enum PiexError {
    /// Out of bounds data access.
    #[error(
        "Out of bounds read: requested offset {requested_offset} + size {requested_size}, but total length is {total_len}"
    )]
    OutOfBounds {
        /// Requested offset in bytes.
        requested_offset: usize,
        /// Requested size in bytes.
        requested_size: usize,
        /// Total length of the container or stream.
        total_len: usize,
    },

    /// Stream read operation failed.
    #[error(
        "Stream read failure at offset {offset}: expected {expected} bytes, available {available} bytes"
    )]
    StreamReadFailed {
        /// Offset where read was attempted.
        offset: usize,
        /// Expected byte count.
        expected: usize,
        /// Available byte count.
        available: usize,
    },

    /// Invalid header encountered.
    #[error("Invalid header magic or signature: {0}")]
    InvalidHeader(String),

    /// Invalid ISOBMFF box header encountered.
    #[error("Invalid box header at offset {offset}")]
    InvalidBoxHeader {
        /// Offset where the invalid box header was found.
        offset: u64,
    },

    /// Corrupt data encountered.
    #[error("Corrupt image data structure: {0}")]
    CorruptData(String),

    /// Format unsupported.
    #[error("Unsupported RAW format or feature")]
    UnsupportedFormat,

    /// Underlying I/O error message.
    #[error("I/O failure: {0}")]
    IoError(String),

    /// Generic failure message.
    #[error("Extraction failed: {0}")]
    GenericFail(String),
}

impl From<std::io::Error> for PiexError {
    fn from(err: std::io::Error) -> Self {
        PiexError::IoError(err.to_string())
    }
}

/// Supported image compression and encoding formats.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
#[repr(C)]
pub enum ImageFormat {
    /// JPEG compressed image data.
    #[default]
    JpegCompressed = 0,
    /// Uncompressed RGB image data.
    UncompressedRgb = 1,
    /// HEVC compressed image data.
    HevcCompressed = 2,
}

/// Supported RAW image types and non-RAW identifier.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
#[repr(C)]
pub enum RawImageType {
    /// Non-RAW or unknown image format.
    #[default]
    NonRaw = 0,
    /// Sony ARW format.
    Arw,
    /// Canon CR2 format.
    Cr2,
    /// Canon CR3 format.
    Cr3,
    /// Canon CRW format.
    Crw,
    /// Kodak DCR format.
    Dcr,
    /// Digital Negative (DNG) format.
    Dng,
    /// Epson ERF format.
    Erf,
    /// Kodak KDC format.
    Kdc,
    /// Leaf MOS format.
    Mos,
    /// Minolta MRW format.
    Mrw,
    /// Nikon NEF format.
    Nef,
    /// Nikon NRW format.
    Nrw,
    /// Olympus ORF format.
    Orf,
    /// Pentax PEF format.
    Pef,
    /// Phase One RAW format.
    PhaseOne,
    /// Logitech Fotoman PXN format.
    Pxn,
    /// Apple QTK format.
    Qtk,
    /// Fuji RAF format.
    Raf,
    /// Hasselblad 3FR format.
    Raw3fr,
    /// Contax N RAW format.
    RawContaxN,
    /// Rollei RDC RAW format.
    Rdc,
    /// Panasonic RW2 format.
    Rw2,
    /// Samsung SRW format.
    Srw,
    /// Sigma / Polaroid X3F format.
    X3f,
}

impl RawImageType {
    /// Converts a raw `u32` discriminant to `RawImageType`, returning `None` if out of range.
    pub fn from_u32(val: u32) -> Option<Self> {
        match val {
            0 => Some(Self::NonRaw),
            1 => Some(Self::Arw),
            2 => Some(Self::Cr2),
            3 => Some(Self::Cr3),
            4 => Some(Self::Crw),
            5 => Some(Self::Dcr),
            6 => Some(Self::Dng),
            7 => Some(Self::Erf),
            8 => Some(Self::Kdc),
            9 => Some(Self::Mos),
            10 => Some(Self::Mrw),
            11 => Some(Self::Nef),
            12 => Some(Self::Nrw),
            13 => Some(Self::Orf),
            14 => Some(Self::Pef),
            15 => Some(Self::PhaseOne),
            16 => Some(Self::Pxn),
            17 => Some(Self::Qtk),
            18 => Some(Self::Raf),
            19 => Some(Self::Raw3fr),
            20 => Some(Self::RawContaxN),
            21 => Some(Self::Rdc),
            22 => Some(Self::Rw2),
            23 => Some(Self::Srw),
            24 => Some(Self::X3f),
            _ => None,
        }
    }
}

/// Defines the properties of an embedded preview or thumbnail image.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[repr(C)]
pub struct Image {
    /// Width of the image in pixels.
    pub width: u16,
    /// Height of the image in pixels.
    pub height: u16,
    /// Length of the image data in bytes.
    pub length: u32,
    /// Offset of the image data in bytes from the start of the stream.
    pub offset: u32,
    /// Format of the image data.
    pub format: ImageFormat,
}

impl Image {
    /// Creates a new `Image` instance with the specified parameters.
    pub fn new(width: u16, height: u16, length: u32, offset: u32, format: ImageFormat) -> Self {
        Self { width, height, length, offset, format }
    }

    /// Checks if the image is valid (i.e. length is greater than 0).
    pub fn is_valid(&self) -> bool {
        self.length > 0
    }

    /// Returns `true` if both width and height are strictly greater than `rhs`'s dimensions.
    pub fn is_larger_than(&self, rhs: &Image) -> bool {
        self.width > rhs.width && self.height > rhs.height
    }
}

/// Color space of the preview image.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
#[repr(C)]
pub enum ColorSpace {
    /// sRGB color space.
    #[default]
    Srgb = 0,
    /// Adobe RGB color space.
    AdobeRgb = 1,
}

/// Represents a fractional number with a 32-bit unsigned numerator and denominator.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(C)]
pub struct Rational {
    /// Numerator value.
    pub numerator: u32,
    /// Denominator value.
    pub denominator: u32,
}

impl Rational {
    /// Creates a new `Rational` fraction.
    pub fn new(numerator: u32, denominator: u32) -> Self {
        Self { numerator, denominator }
    }

    /// Converts the rational fraction to a 64-bit floating point number.
    /// Returns 0.0 if denominator is 0.
    pub fn as_f64(&self) -> f64 {
        if self.denominator == 0 {
            0.0
        } else {
            self.numerator as f64 / self.denominator as f64
        }
    }
}

impl Default for Rational {
    fn default() -> Self {
        Self { numerator: 0, denominator: 1 }
    }
}

/// Represents a fractional number with a 32-bit signed numerator and denominator.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct SRational {
    /// Numerator value.
    pub numerator: i32,
    /// Denominator value.
    pub denominator: i32,
}

impl SRational {
    /// Creates a new `SRational` fraction.
    pub fn new(numerator: i32, denominator: i32) -> Self {
        Self { numerator, denominator }
    }

    /// Converts the rational fraction to a 64-bit floating point number.
    /// Returns 0.0 if denominator is 0.
    pub fn as_f64(&self) -> f64 {
        if self.denominator == 0 {
            0.0
        } else {
            self.numerator as f64 / self.denominator as f64
        }
    }
}

impl Default for SRational {
    fn default() -> Self {
        Self { numerator: 0, denominator: 1 }
    }
}

pub use crate::cxx_bridge::ffi::{GpsData, PreviewImageData};

impl Default for GpsData {
    fn default() -> Self {
        Self {
            is_valid: false,
            latitude_ref: b'N',
            latitude: [Rational::default(); 3],
            longitude_ref: b'E',
            longitude: [Rational::default(); 3],
            altitude_ref: false,
            altitude: Rational::default(),
            time_stamp: [Rational::default(); 3],
            date_stamp: Vec::new(),
        }
    }
}

/// Structural dimensions and color filter array (CFA) metadata for Adobe DNG images.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Default)]
pub struct DngInformation {
    /// Full DNG image width.
    pub width: u32,
    /// Full DNG image height.
    pub height: u32,
    /// Color Filter Array (CFA) pattern dimensions (width, height).
    pub cfa_pattern_dim: [u32; 2],
}

impl Default for PreviewImageData {
    fn default() -> Self {
        Self {
            preview: Image::default(),
            thumbnail: Image::default(),
            exif_orientation: 1,
            color_space: ColorSpace::Srgb,
            full_width: 0,
            full_height: 0,
            maker: Vec::new(),
            model: Vec::new(),
            date_time: Vec::new(),
            iso: 0,
            exposure_time: Rational::default(),
            fnumber: Rational::default(),
            focal_length: Rational::default(),
            gps: GpsData::default(),
            cfa_pattern_dim: [0, 0],
            cfa_pattern_dim_set: false,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_raw_image_type_default() {
        assert_eq!(RawImageType::default(), RawImageType::NonRaw);
    }

    #[test]
    fn test_raw_image_type_equality_and_copy() {
        let t1 = RawImageType::Arw;
        let t2 = t1;
        assert_eq!(t1, t2);
        assert_ne!(RawImageType::Arw, RawImageType::Dng);
    }

    #[test]
    fn test_raw_image_type_debug_formatting() {
        let formatted = format!("{:?}", RawImageType::Cr2);
        assert_eq!(formatted, "Cr2");
    }

    #[test]
    fn test_raw_image_type_c_repr_values() {
        assert_eq!(RawImageType::NonRaw as i32, 0);
        assert_eq!(RawImageType::Arw as i32, 1);
        assert_eq!(RawImageType::Erf as i32, 7);
        assert_eq!(RawImageType::PhaseOne as i32, 15);
        assert_eq!(RawImageType::Pxn as i32, 16);
        assert_eq!(RawImageType::Rdc as i32, 21);
        assert_eq!(RawImageType::X3f as i32, 24);
    }

    #[test]
    fn test_dng_information() {
        let dng = DngInformation { width: 4000, height: 3000, cfa_pattern_dim: [2, 2] };
        assert_eq!(dng.width, 4000);
        assert_eq!(dng.height, 3000);
        assert_eq!(dng.cfa_pattern_dim, [2, 2]);

        let default_dng = DngInformation::default();
        assert_eq!(default_dng.width, 0);
        assert_eq!(default_dng.height, 0);
        assert_eq!(default_dng.cfa_pattern_dim, [0, 0]);

        assert_eq!(dng.clone(), dng);
        let debug = format!("{:?}", dng);
        assert!(debug.contains("4000"));
    }
}
