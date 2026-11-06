use std::panic::{catch_unwind, AssertUnwindSafe};
use std::pin::Pin;

use crate::cr3;
use crate::extractor;
use crate::format_recognizer;
use crate::stream::StreamInterface;
use crate::types::{ColorSpace, Image, ImageFormat, PiexError, Rational, RawImageType};

unsafe impl cxx::ExternType for ImageFormat {
    type Id = cxx::type_id!("piex_bridge::ImageFormat");
    type Kind = cxx::kind::Trivial;
}

unsafe impl cxx::ExternType for Image {
    type Id = cxx::type_id!("piex_bridge::Image");
    type Kind = cxx::kind::Trivial;
}

unsafe impl cxx::ExternType for ColorSpace {
    type Id = cxx::type_id!("piex_bridge::ColorSpace");
    type Kind = cxx::kind::Trivial;
}

unsafe impl cxx::ExternType for Rational {
    type Id = cxx::type_id!("piex_bridge::Rational");
    type Kind = cxx::kind::Trivial;
}

#[cxx::bridge(namespace = "piex_bridge")]
pub mod ffi {
    /// GPS metadata associated with an image.
    #[derive(Debug, Clone, PartialEq, Eq)]
    pub struct GpsData {
        /// Indicates if the GPS data is valid to use.
        pub is_valid: bool,
        /// Latitude reference ('N' or 'S').
        pub latitude_ref: u8,
        /// Latitude degrees, minutes, seconds.
        pub latitude: [Rational; 3],
        /// Longitude reference ('E' or 'W').
        pub longitude_ref: u8,
        /// Longitude degrees, minutes, seconds.
        pub longitude: [Rational; 3],
        /// Altitude reference (true for above sea level, false for below).
        pub altitude_ref: bool,
        /// Altitude value.
        pub altitude: Rational,
        /// Time stamp (hours, minutes, seconds).
        pub time_stamp: [Rational; 3],
        /// Raw date stamp bytes, normally "YYYY:MM:DD". Not necessarily UTF-8.
        pub date_stamp: Vec<u8>,
    }

    /// Contains extracted preview images and metadata.
    #[derive(Debug, Clone, PartialEq, Eq)]
    pub struct PreviewImageData {
        /// Main preview image information.
        pub preview: Image,
        /// Thumbnail image information.
        pub thumbnail: Image,
        /// EXIF orientation value (1 is default horizontal).
        pub exif_orientation: u32,
        /// Color space of the preview.
        pub color_space: ColorSpace,
        /// Full RAW image width in pixels.
        pub full_width: u32,
        /// Full RAW image height in pixels.
        pub full_height: u32,
        /// Raw camera manufacturer bytes (TIFF Make), including any trailing NUL.
        /// Not necessarily UTF-8.
        pub maker: Vec<u8>,
        /// Raw camera model bytes (TIFF Model), including any trailing NUL.
        /// Not necessarily UTF-8.
        pub model: Vec<u8>,
        /// Raw date/time bytes (EXIF DateTimeOriginal), including any trailing NUL.
        /// Not necessarily UTF-8.
        pub date_time: Vec<u8>,
        /// ISO speed rating.
        pub iso: u32,
        /// Exposure time in seconds.
        pub exposure_time: Rational,
        /// F-number aperture value.
        pub fnumber: Rational,
        /// Focal length in millimeters.
        pub focal_length: Rational,
        /// GPS coordinates and timestamp.
        pub gps: GpsData,
        /// Color Filter Array (CFA) pattern dimensions (width, height).
        pub cfa_pattern_dim: [u32; 2],
    }

    unsafe extern "C++" {
        include!("rust/bridge/piex_cxx_glue.h");
        include!("src/piex_types.h");

        type ImageFormat = crate::types::ImageFormat;
        type Image = crate::types::Image;
        type ColorSpace = crate::types::ColorSpace;
        type Rational = crate::types::Rational;
        type StreamInterface;

        unsafe fn StreamGetData(
            stream: Pin<&mut StreamInterface>,
            offset: usize,
            length: usize,
            data: *mut u8,
        ) -> i32;
    }

    extern "Rust" {
        fn rust_bytes_to_vec(bytes: &[u8]) -> Vec<u8>;
        fn rust_bytes_required_for_is_raw() -> usize;
        fn rust_is_raw(stream: Pin<&mut StreamInterface>) -> bool;
        fn rust_get_preview_image_data(
            stream: Pin<&mut StreamInterface>,
            out_data: &mut PreviewImageData,
            out_type: &mut u32,
        ) -> i32;
        fn rust_get_dng_information(
            stream: Pin<&mut StreamInterface>,
            out_width: &mut u32,
            out_height: &mut u32,
            out_cfa: &mut [u32; 2],
        ) -> bool;
        fn rust_get_orientation(
            stream: Pin<&mut StreamInterface>,
            out_orientation: &mut u32,
        ) -> bool;
        fn rust_cr3_get_orientation(
            stream: Pin<&mut StreamInterface>,
            out_orientation: &mut u32,
        ) -> bool;
        fn rust_cr3_get_preview_data(
            stream: Pin<&mut StreamInterface>,
            out_data: &mut PreviewImageData,
        ) -> i32;
        fn rust_is_raw_lite(data: &[u8]) -> bool;
        fn rust_recognize_raw_image_type_lite(data: &[u8]) -> u32;
        fn rust_is_raw_type(raw_type: u32) -> bool;
        fn rust_is_of_type(data: &[u8], raw_type: u32) -> bool;
        fn rust_get_number_of_bytes_for_is_raw_lite() -> usize;
        fn rust_get_number_of_bytes_for_is_of_type(raw_type: u32) -> usize;
    }
}

pub struct CxxStreamAdapter<'a> {
    stream: std::sync::Mutex<Pin<&'a mut ffi::StreamInterface>>,
    length: usize,
}

impl<'a> CxxStreamAdapter<'a> {
    pub fn new(stream: Pin<&'a mut ffi::StreamInterface>) -> Self {
        Self { stream: std::sync::Mutex::new(stream), length: usize::MAX }
    }
}

impl StreamInterface for CxxStreamAdapter<'_> {
    fn get_data(&self, offset: usize, data: &mut [u8]) -> Result<(), PiexError> {
        let req_len = data.len();
        if offset > self.length {
            return Err(PiexError::OutOfBounds {
                requested_offset: offset,
                requested_size: req_len,
                total_len: self.length,
            });
        }
        if req_len == 0 {
            return Ok(());
        }

        let _end = offset.checked_add(req_len).ok_or_else(|| PiexError::OutOfBounds {
            requested_offset: offset,
            requested_size: req_len,
            total_len: self.length,
        })?;

        let mut stream_guard = self.stream.lock().map_err(|_| {
            PiexError::GenericFail("Mutex poisoned during CXX stream read".to_string())
        })?;

        // SAFETY:
        // - `stream_guard` contains a valid pinned C++ StreamInterface.
        // - `data` is a mutable slice of length `req_len`.
        // - `data.as_mut_ptr()` points to at least `req_len` bytes.
        let status = unsafe {
            ffi::StreamGetData(stream_guard.as_mut(), offset, req_len, data.as_mut_ptr())
        };

        if status == 0 {
            Ok(())
        } else {
            Err(PiexError::StreamReadFailed {
                offset,
                expected: req_len,
                available: self.length.saturating_sub(offset),
            })
        }
    }

    fn len(&self) -> usize {
        self.length
    }
}

/// Copies `bytes` into a new `Vec<u8>` with a single memcpy. The C++ bridge
/// uses this because `rust::Vec` has no public way to set its length.
fn rust_bytes_to_vec(bytes: &[u8]) -> Vec<u8> {
    bytes.to_vec()
}

fn rust_bytes_required_for_is_raw() -> usize {
    catch_unwind(AssertUnwindSafe(|| format_recognizer::get_number_of_bytes_for_is_raw_lite()))
        .unwrap_or(0)
}

fn rust_is_raw(stream: Pin<&mut ffi::StreamInterface>) -> bool {
    catch_unwind(AssertUnwindSafe(|| {
        let adapter = CxxStreamAdapter::new(stream);
        let req_bytes = extractor::bytes_required_for_is_raw();
        let mut header = vec![0u8; 10_000];
        let to_read = req_bytes.min(header.len());
        if adapter.get_data(0, &mut header[..to_read]).is_err() {
            return false;
        }
        format_recognizer::recognize_raw_image_type_lite(&header[..to_read]) != RawImageType::NonRaw
    }))
    .unwrap_or(false)
}

fn rust_get_preview_image_data(
    stream: Pin<&mut ffi::StreamInterface>,
    out_data: &mut ffi::PreviewImageData,
    out_type: &mut u32,
) -> i32 {
    catch_unwind(AssertUnwindSafe(|| {
        let adapter = CxxStreamAdapter::new(stream);
        let req_bytes = extractor::bytes_required_for_is_raw();
        let mut header = vec![0u8; 10_000];
        let to_read = req_bytes.min(header.len());
        if adapter.get_data(0, &mut header[..to_read]).is_err() {
            return 1;
        }
        let recognized = format_recognizer::recognize_raw_image_type_lite(&header[..to_read]);
        *out_type = recognized as u32;
        match extractor::get_preview_image_data_into(&adapter, out_data) {
            Ok(raw_type) => {
                *out_type = raw_type as u32;
                0
            }
            Err(PiexError::UnsupportedFormat) => 2,
            Err(_) => 1,
        }
    }))
    .unwrap_or(1)
}

fn rust_get_dng_information(
    stream: Pin<&mut ffi::StreamInterface>,
    out_width: &mut u32,
    out_height: &mut u32,
    out_cfa: &mut [u32; 2],
) -> bool {
    catch_unwind(AssertUnwindSafe(|| {
        let adapter = CxxStreamAdapter::new(stream);
        extractor::get_dng_information_out(&adapter, out_width, out_height, out_cfa)
    }))
    .unwrap_or(false)
}

fn rust_get_orientation(stream: Pin<&mut ffi::StreamInterface>, out_orientation: &mut u32) -> bool {
    catch_unwind(AssertUnwindSafe(|| {
        let adapter = CxxStreamAdapter::new(stream);
        match extractor::get_orientation(&adapter) {
            Ok(orient) => {
                *out_orientation = orient;
                true
            }
            Err(_) => false,
        }
    }))
    .unwrap_or(false)
}

fn rust_cr3_get_orientation(
    stream: Pin<&mut ffi::StreamInterface>,
    out_orientation: &mut u32,
) -> bool {
    catch_unwind(AssertUnwindSafe(|| {
        let adapter = CxxStreamAdapter::new(stream);
        match cr3::extract_cr3_orientation(&adapter) {
            Ok(orient) => {
                *out_orientation = orient;
                true
            }
            Err(_) => false,
        }
    }))
    .unwrap_or(false)
}

fn rust_cr3_get_preview_data(
    stream: Pin<&mut ffi::StreamInterface>,
    out_data: &mut ffi::PreviewImageData,
) -> i32 {
    catch_unwind(AssertUnwindSafe(|| {
        let adapter = CxxStreamAdapter::new(stream);
        match cr3::extract_cr3_preview_data_into(&adapter, out_data) {
            Ok(()) => 0,
            Err(PiexError::UnsupportedFormat) => 2,
            Err(_) => 1,
        }
    }))
    .unwrap_or(1)
}

fn rust_is_raw_lite(data: &[u8]) -> bool {
    catch_unwind(AssertUnwindSafe(|| format_recognizer::is_raw_lite(data))).unwrap_or(false)
}

fn rust_recognize_raw_image_type_lite(data: &[u8]) -> u32 {
    catch_unwind(AssertUnwindSafe(|| format_recognizer::recognize_raw_image_type_lite(data) as u32))
        .unwrap_or(0)
}

fn rust_is_raw_type(raw_type: u32) -> bool {
    catch_unwind(AssertUnwindSafe(
        || matches!(RawImageType::from_u32(raw_type), Some(t) if t != RawImageType::NonRaw),
    ))
    .unwrap_or(false)
}

fn rust_is_of_type(data: &[u8], raw_type: u32) -> bool {
    catch_unwind(AssertUnwindSafe(|| {
        RawImageType::from_u32(raw_type).is_some_and(|t| format_recognizer::is_of_type(data, t))
    }))
    .unwrap_or(false)
}

fn rust_get_number_of_bytes_for_is_raw_lite() -> usize {
    catch_unwind(AssertUnwindSafe(|| format_recognizer::get_number_of_bytes_for_is_raw_lite()))
        .unwrap_or(0)
}

fn rust_get_number_of_bytes_for_is_of_type(raw_type: u32) -> usize {
    catch_unwind(AssertUnwindSafe(|| {
        RawImageType::from_u32(raw_type)
            .map_or(0, format_recognizer::get_number_of_bytes_for_is_of_type)
    }))
    .unwrap_or(0)
}
