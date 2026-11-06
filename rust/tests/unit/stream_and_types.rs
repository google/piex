use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};
use std::io::{Cursor, ErrorKind, Seek, SeekFrom};

use piex::stream::{ByteSliceStream, Endian, EndianReader, StreamAdapter, StreamInterface};
use piex::types::{
    ColorSpace, GpsData, Image, ImageFormat, PiexError, PreviewImageData, Rational, RawImageType,
};

// ============================================================================
// 1. Endian Reader Tests (Primitives: u8, u16, i16, u32, i32, rational)
// ============================================================================

#[test]
fn test_read_u8() {
    let data = [0x00, 0x42, 0xFF, 0x12];
    assert_eq!(data.read_u8(0).unwrap(), 0x00);
    assert_eq!(data.read_u8(1).unwrap(), 0x42);
    assert_eq!(data.read_u8(2).unwrap(), 0xFF);
    assert_eq!(data.read_u8(3).unwrap(), 0x12);

    // Out of bounds
    assert!(matches!(data.read_u8(4), Err(PiexError::OutOfBounds { .. })));
    assert!(matches!(data.read_u8(usize::MAX), Err(PiexError::OutOfBounds { .. })));
}

#[test]
fn test_read_u16() {
    let data_le = [0x34, 0x12, 0x00, 0x00];
    assert_eq!(data_le.read_u16(0, Endian::Little).unwrap(), 0x1234);

    let data_be = [0x12, 0x34, 0x00, 0x00];
    assert_eq!(data_be.read_u16(0, Endian::Big).unwrap(), 0x1234);

    // Boundary check: offset 3 leaves only 1 byte available
    assert!(matches!(
        data_le.read_u16(3, Endian::Little),
        Err(PiexError::OutOfBounds { requested_offset: 3, requested_size: 2, total_len: 4 })
    ));
    assert!(matches!(
        data_be.read_u16(4, Endian::Big),
        Err(PiexError::OutOfBounds { requested_offset: 4, requested_size: 2, total_len: 4 })
    ));
}

#[test]
fn test_read_i16() {
    // -2 in 16-bit 2's complement is 0xFFFE
    let data_le = [0xFE, 0xFF, 0x00, 0x80]; // 0x8000 = -32768 in LE is [0x00, 0x80]
    assert_eq!(data_le.read_i16(0, Endian::Little).unwrap(), -2);
    assert_eq!(data_le.read_i16(2, Endian::Little).unwrap(), -32768);

    let data_be = [0xFF, 0xFE, 0x80, 0x00];
    assert_eq!(data_be.read_i16(0, Endian::Big).unwrap(), -2);
    assert_eq!(data_be.read_i16(2, Endian::Big).unwrap(), -32768);

    // Positive i16
    let pos_data = [0xD2, 0x04]; // 1234 in LE [0xD2, 0x04]
    assert_eq!(pos_data.read_i16(0, Endian::Little).unwrap(), 1234);

    // Out of bounds
    assert!(matches!(pos_data.read_i16(1, Endian::Little), Err(PiexError::OutOfBounds { .. })));
}

#[test]
fn test_read_u32() {
    let data_le = [0x78, 0x56, 0x34, 0x12, 0x00];
    assert_eq!(data_le.read_u32(0, Endian::Little).unwrap(), 0x12345678);

    let data_be = [0x12, 0x34, 0x56, 0x78, 0x00];
    assert_eq!(data_be.read_u32(0, Endian::Big).unwrap(), 0x12345678);

    // Out of bounds when requesting 4 bytes with only 3 remaining
    assert!(matches!(
        data_le.read_u32(2, Endian::Little),
        Err(PiexError::OutOfBounds { requested_offset: 2, requested_size: 4, total_len: 5 })
    ));
}

#[test]
fn test_read_i32() {
    // -1 = 0xFFFFFFFF
    let neg_one = [0xFF, 0xFF, 0xFF, 0xFF];
    assert_eq!(neg_one.read_i32(0, Endian::Little).unwrap(), -1);
    assert_eq!(neg_one.read_i32(0, Endian::Big).unwrap(), -1);

    // i32::MIN = 0x80000000 -> LE: [0x00, 0x00, 0x00, 0x80], BE: [0x80, 0x00, 0x00, 0x00]
    let min_le = [0x00, 0x00, 0x00, 0x80];
    assert_eq!(min_le.read_i32(0, Endian::Little).unwrap(), i32::MIN);

    let min_be = [0x80, 0x00, 0x00, 0x00];
    assert_eq!(min_be.read_i32(0, Endian::Big).unwrap(), i32::MIN);

    // Out of bounds
    assert!(matches!(min_le.read_i32(1, Endian::Little), Err(PiexError::OutOfBounds { .. })));
}

#[test]
fn test_read_rational() {
    // Rational: 2 uint32s (num, den). 10 / 2
    let rat_le = [
        0x0A, 0x00, 0x00, 0x00, // Num = 10
        0x02, 0x00, 0x00, 0x00, // Den = 2
    ];
    let r_le = rat_le.read_rational(0, Endian::Little).unwrap();
    assert_eq!(r_le, Rational::new(10, 2));

    let rat_be = [
        0x00, 0x00, 0x00, 0x0A, // Num = 10
        0x00, 0x00, 0x00, 0x02, // Den = 2
    ];
    let r_be = rat_be.read_rational(0, Endian::Big).unwrap();
    assert_eq!(r_be, Rational::new(10, 2));

    // 7 bytes available -> OutOfBounds
    let short_data = [0u8; 7];
    assert!(matches!(
        short_data.read_rational(0, Endian::Little),
        Err(PiexError::OutOfBounds { .. })
    ));
}

// ============================================================================
// 2. Boundary Validation, Slice Reading & Offset Overflows
// ============================================================================

#[test]
fn test_read_bytes_boundaries() {
    let data = b"Hello, World!";

    // Zero-length reads within bounds
    assert_eq!(data.read_bytes(0, 0).unwrap(), b"");
    assert_eq!(data.read_bytes(5, 0).unwrap(), b"");
    assert_eq!(data.read_bytes(data.len(), 0).unwrap(), b"");
    // Zero-length read out of bounds (offset > len)
    assert!(matches!(data.read_bytes(data.len() + 1, 0), Err(PiexError::OutOfBounds { .. })));

    // Full slice read
    assert_eq!(data.read_bytes(0, data.len()).unwrap(), data);

    // Partial slice read
    assert_eq!(data.read_bytes(7, 5).unwrap(), b"World");

    // Out of bounds reads
    assert!(matches!(
        data.read_bytes(0, data.len() + 1),
        Err(PiexError::OutOfBounds { requested_offset: 0, requested_size: 14, total_len: 13 })
    ));
    assert!(matches!(
        data.read_bytes(10, 5),
        Err(PiexError::OutOfBounds { requested_offset: 10, requested_size: 5, total_len: 13 })
    ));

    // Arithmetic overflow in offset + length calculation
    assert!(matches!(data.read_bytes(usize::MAX, 1), Err(PiexError::OutOfBounds { .. })));
    assert!(matches!(data.read_bytes(1, usize::MAX), Err(PiexError::OutOfBounds { .. })));
    assert!(matches!(data.read_bytes(usize::MAX - 1, 2), Err(PiexError::OutOfBounds { .. })));
}

// ============================================================================
// 3. String Extraction & Null Termination Handling
// ============================================================================

#[test]
fn test_read_string() {
    // Null-terminated string within buffer
    let data1 = b"Canon\0EOS 5D";
    assert_eq!(data1.read_string(0, 12).unwrap(), "Canon");
    assert_eq!(data1.read_string(6, 6).unwrap(), "EOS 5D");

    // String without null byte
    let data2 = b"Nikon";
    assert_eq!(data2.read_string(0, 5).unwrap(), "Nikon");

    // Empty string (zero length requested)
    assert_eq!(data2.read_string(0, 0).unwrap(), "");

    // Leading null byte
    let data3 = b"\0Camera";
    assert_eq!(data3.read_string(0, 7).unwrap(), "");

    // String with lossy UTF-8 characters
    let data4 = b"PIEX\xFF\xFE\0REST";
    let s = data4.read_string(0, 8).unwrap();
    assert!(s.starts_with("PIEX"));

    // String read exceeding slice bounds
    assert!(matches!(data2.read_string(0, 10), Err(PiexError::OutOfBounds { .. })));
}

// ============================================================================
// 4. Rational Type & Division-by-Zero Handling
// ============================================================================

#[test]
fn test_rational_operations() {
    let r1 = Rational::new(1, 2);
    assert_eq!(r1.numerator, 1);
    assert_eq!(r1.denominator, 2);
    assert_eq!(r1.as_f64(), 0.5);

    let r2 = Rational::new(10, 3);
    assert!((r2.as_f64() - 3.3333333333333335).abs() < 1e-10);

    let r_zero_num = Rational::new(0, 5);
    assert_eq!(r_zero_num.as_f64(), 0.0);

    // Division by zero denominator
    let r_zero_den = Rational::new(100, 0);
    assert_eq!(r_zero_den.as_f64(), 0.0);

    // Default rational
    let r_def = Rational::default();
    assert_eq!(r_def.numerator, 0);
    assert_eq!(r_def.denominator, 1);
    assert_eq!(r_def.as_f64(), 0.0);

    // Copy & Clone
    let r_copy = r1;
    assert_eq!(r1, r_copy);

    // Debug representation
    let dbg = format!("{:?}", r1);
    assert!(dbg.contains("Rational"));

    // Hash implementation
    let mut hasher1 = DefaultHasher::new();
    let mut hasher2 = DefaultHasher::new();
    r1.hash(&mut hasher1);
    r_copy.hash(&mut hasher2);
    assert_eq!(hasher1.finish(), hasher2.finish());
}

// ============================================================================
// 5. ByteSliceStream & SubStream Slicing / Nesting
// ============================================================================

#[test]
fn test_byte_slice_stream_operations() {
    let buffer = b"0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    let stream = ByteSliceStream::new(buffer);

    assert_eq!(stream.len(), buffer.len());
    assert!(!stream.is_empty());

    // Reading data
    let mut out = [0u8; 10];
    stream.get_data(0, &mut out).unwrap();
    assert_eq!(&out, b"0123456789");

    stream.get_data(10, &mut out).unwrap();
    assert_eq!(&out, b"ABCDEFGHIJ");

    // Zero-length read
    let mut empty_buf = [0u8; 0];
    stream.get_data(5, &mut empty_buf).unwrap();

    // Read extending beyond stream boundary
    let mut big_out = [0u8; 50];
    assert!(matches!(
        stream.get_data(30, &mut big_out),
        Err(PiexError::OutOfBounds { requested_offset: 30, requested_size: 50, total_len: 36 })
    ));

    // Overflow check
    assert!(matches!(
        stream.get_data(usize::MAX, &mut big_out),
        Err(PiexError::OutOfBounds { .. })
    ));
}

#[test]
fn test_byte_slice_sub_stream_nesting() {
    let buffer = b"ABC------------------XYZ";
    let stream = ByteSliceStream::new(buffer);

    // Sub-stream 1
    let sub1 = stream.sub_stream(3, 18).unwrap();
    assert_eq!(sub1.len(), 18);

    // Sub-stream 2 nested inside sub1
    let sub2 = sub1.sub_stream(2, 10).unwrap();
    assert_eq!(sub2.len(), 10);

    // Sub-stream with 0 length
    let empty_sub = stream.sub_stream(5, 0).unwrap();
    assert_eq!(empty_sub.len(), 0);
    assert!(empty_sub.is_empty());

    // Sub-stream out of bounds
    assert!(matches!(
        stream.sub_stream(20, 10),
        Err(PiexError::OutOfBounds { requested_offset: 20, requested_size: 10, total_len: 24 })
    ));
    assert!(matches!(stream.sub_stream(usize::MAX, 1), Err(PiexError::OutOfBounds { .. })));
}

#[test]
fn test_byte_slice_sub_stream_zero_length_out_of_bounds() {
    let buffer = b"0123456789";
    let stream = ByteSliceStream::new(buffer);

    // Zero-length within bounds
    let empty_sub_0 = stream.sub_stream(0, 0).unwrap();
    assert_eq!(empty_sub_0.len(), 0);
    let empty_sub_5 = stream.sub_stream(5, 0).unwrap();
    assert_eq!(empty_sub_5.len(), 0);

    // Zero-length at exact end boundary
    let empty_sub_end = stream.sub_stream(10, 0).unwrap();
    assert_eq!(empty_sub_end.len(), 0);

    // Zero-length out of bounds (offset > len)
    let res_oob = stream.sub_stream(11, 0);
    assert!(
        matches!(
            res_oob,
            Err(PiexError::OutOfBounds { requested_offset: 11, requested_size: 0, total_len: 10 })
        ),
        "Expected OutOfBounds for zero-length sub-stream past end, got {:?}",
        res_oob
    );

    let res_max = stream.sub_stream(usize::MAX, 0);
    assert!(
        matches!(res_max, Err(PiexError::OutOfBounds { .. })),
        "Expected OutOfBounds for zero-length sub-stream at usize::MAX, got {:?}",
        res_max
    );
}

// ============================================================================
// 6. StreamAdapter (Cursor wrapping & seeking behavior)
// ============================================================================

#[test]
fn test_stream_adapter_cursor() {
    let buffer = vec![0u8, 1, 2, 3, 4, 5, 6, 7, 8, 9];
    let cursor = Cursor::new(buffer.clone());
    let adapter = StreamAdapter::new(cursor).unwrap();

    assert_eq!(adapter.len(), 10);
    assert!(!adapter.is_empty());

    // Get data at offset 0
    let mut out = [0u8; 4];
    adapter.get_data(0, &mut out).unwrap();
    assert_eq!(&out, &[0, 1, 2, 3]);

    // Get data at non-zero offset
    adapter.get_data(5, &mut out).unwrap();
    assert_eq!(&out, &[5, 6, 7, 8]);

    // Verify repeated access at offset 0 works (seeking behavior)
    adapter.get_data(0, &mut out[..2]).unwrap();
    assert_eq!(&out[..2], &[0, 1]);

    // Zero-length read
    adapter.get_data(2, &mut out[..0]).unwrap();

    // Read past length
    let mut large = [0u8; 10];
    assert!(matches!(
        adapter.get_data(8, &mut large),
        Err(PiexError::OutOfBounds { requested_offset: 8, requested_size: 10, total_len: 10 })
    ));

    // Overflow check
    assert!(matches!(adapter.get_data(usize::MAX, &mut large), Err(PiexError::OutOfBounds { .. })));
}

#[test]
fn test_stream_adapter_empty_and_positions() {
    let empty_cursor = Cursor::new(Vec::<u8>::new());
    let adapter = StreamAdapter::new(empty_cursor).unwrap();

    assert_eq!(adapter.len(), 0);
    assert!(adapter.is_empty());

    // Cursor initialized at offset non-zero
    let mut cursor = Cursor::new(vec![10, 20, 30, 40]);
    cursor.seek(SeekFrom::Start(2)).unwrap();
    let adapter2 = StreamAdapter::new(cursor).unwrap();

    // The length should still be calculated as 4
    assert_eq!(adapter2.len(), 4);
    let mut out = [0u8; 2];
    adapter2.get_data(0, &mut out).unwrap();
    assert_eq!(&out, &[10, 20]);
}

struct CustomLenReader {
    len: u64,
    pos: u64,
}

impl std::io::Read for CustomLenReader {
    fn read(&mut self, _buf: &mut [u8]) -> std::io::Result<usize> {
        Ok(0)
    }
}

impl std::io::Seek for CustomLenReader {
    fn seek(&mut self, pos: SeekFrom) -> std::io::Result<u64> {
        match pos {
            SeekFrom::Start(p) => {
                self.pos = p;
                Ok(self.pos)
            }
            SeekFrom::End(p) => {
                let new_pos = (self.len as i64).saturating_add(p);
                self.pos = new_pos as u64;
                Ok(self.pos)
            }
            SeekFrom::Current(p) => {
                let new_pos = (self.pos as i64).saturating_add(p);
                self.pos = new_pos as u64;
                Ok(self.pos)
            }
        }
    }
}

#[test]
fn test_stream_adapter_boundary_and_overflow() {
    // Valid stream size fits in usize
    let valid_reader = CustomLenReader { len: 1024, pos: 0 };
    let adapter = StreamAdapter::new(valid_reader).unwrap();
    assert_eq!(adapter.len(), 1024);

    // Stream size exceeding usize capacity (if target pointer width < 64-bit)
    if (u64::MAX as u128) > (usize::MAX as u128) {
        let overflow_reader = CustomLenReader { len: (usize::MAX as u64) + 1, pos: 0 };
        let res = StreamAdapter::new(overflow_reader);
        assert!(matches!(res, Err(PiexError::GenericFail(ref msg)) if msg.contains("capacity")));
    }
}

// ============================================================================
// 7. Image, ImageFormat, ColorSpace, GpsData & PreviewImageData Integrity
// ============================================================================

#[test]
fn test_image_struct() {
    let img = Image::new(1920, 1080, 50000, 1024, ImageFormat::JpegCompressed);
    assert_eq!(img.width, 1920);
    assert_eq!(img.height, 1080);
    assert_eq!(img.length, 50000);
    assert_eq!(img.offset, 1024);
    assert_eq!(img.format, ImageFormat::JpegCompressed);
    assert!(img.is_valid());

    let invalid_img = Image::new(100, 100, 0, 500, ImageFormat::UncompressedRgb);
    assert!(!invalid_img.is_valid());

    // Dimension comparisons
    let small_img = Image::new(800, 600, 1000, 0, ImageFormat::JpegCompressed);
    assert!(img.is_larger_than(&small_img));
    assert!(!small_img.is_larger_than(&img));

    // Edge case: width larger, height equal or smaller
    let same_height = Image::new(2000, 1080, 1000, 0, ImageFormat::JpegCompressed);
    assert!(!same_height.is_larger_than(&img));

    // Default Image
    let def_img = Image::default();
    assert_eq!(def_img.width, 0);
    assert_eq!(def_img.height, 0);
    assert_eq!(def_img.length, 0);
    assert_eq!(def_img.offset, 0);
    assert_eq!(def_img.format, ImageFormat::JpegCompressed);
}

#[test]
fn test_enums_and_traits() {
    // ImageFormat defaults & equality
    assert_eq!(ImageFormat::default(), ImageFormat::JpegCompressed);
    assert_ne!(ImageFormat::JpegCompressed, ImageFormat::HevcCompressed);

    // ColorSpace defaults & equality
    assert_eq!(ColorSpace::default(), ColorSpace::Srgb);
    assert_ne!(ColorSpace::Srgb, ColorSpace::AdobeRgb);

    // RawImageType defaults & equality
    assert_eq!(RawImageType::default(), RawImageType::NonRaw);
    assert_ne!(RawImageType::Arw, RawImageType::Nef);

    // Hash checks for enums
    let mut h1 = DefaultHasher::new();
    let mut h2 = DefaultHasher::new();
    ImageFormat::JpegCompressed.hash(&mut h1);
    ImageFormat::JpegCompressed.hash(&mut h2);
    assert_eq!(h1.finish(), h2.finish());

    let mut h3 = DefaultHasher::new();
    let mut h4 = DefaultHasher::new();
    RawImageType::Arw.hash(&mut h3);
    RawImageType::Arw.hash(&mut h4);
    assert_eq!(h3.finish(), h4.finish());
}

#[test]
fn test_gps_data_defaults() {
    let gps = GpsData::default();
    assert!(!gps.is_valid);
    assert_eq!(gps.latitude_ref, b'N');
    assert_eq!(gps.longitude_ref, b'E');
    assert!(!gps.altitude_ref);
    assert_eq!(gps.latitude, [Rational::default(); 3]);
    assert_eq!(gps.longitude, [Rational::default(); 3]);
    assert_eq!(gps.altitude, Rational::default());
    assert_eq!(gps.time_stamp, [Rational::default(); 3]);
    assert!(gps.date_stamp.is_empty());

    let gps_clone = gps.clone();
    assert_eq!(gps, gps_clone);
}

#[test]
fn test_preview_image_data_defaults() {
    let pid = PreviewImageData::default();
    assert_eq!(pid.preview, Image::default());
    assert_eq!(pid.thumbnail, Image::default());
    assert_eq!(pid.exif_orientation, 1);
    assert_eq!(pid.color_space, ColorSpace::Srgb);
    assert_eq!(pid.full_width, 0);
    assert_eq!(pid.full_height, 0);
    assert!(pid.maker.is_empty());
    assert!(pid.model.is_empty());
    assert!(pid.date_time.is_empty());
    assert_eq!(pid.iso, 0);
    assert_eq!(pid.exposure_time, Rational::default());
    assert_eq!(pid.fnumber, Rational::default());
    assert_eq!(pid.focal_length, Rational::default());
    assert_eq!(pid.gps, GpsData::default());
    assert_eq!(pid.cfa_pattern_dim, [0, 0]);

    let pid_clone = pid.clone();
    assert_eq!(pid, pid_clone);
}

#[test]
fn test_piex_error_display_and_conversion() {
    let e_out =
        PiexError::OutOfBounds { requested_offset: 100, requested_size: 50, total_len: 120 };
    assert_eq!(
        format!("{}", e_out),
        "Out of bounds read: requested offset 100 + size 50, but total length is 120"
    );

    let e_stream = PiexError::StreamReadFailed { offset: 10, expected: 20, available: 5 };
    assert_eq!(
        format!("{}", e_stream),
        "Stream read failure at offset 10: expected 20 bytes, available 5 bytes"
    );

    let e_header = PiexError::InvalidHeader("Bad magic".to_string());
    assert_eq!(format!("{}", e_header), "Invalid header magic or signature: Bad magic");

    let e_corrupt = PiexError::CorruptData("Corrupt IFD".to_string());
    assert_eq!(format!("{}", e_corrupt), "Corrupt image data structure: Corrupt IFD");

    let e_unsupported = PiexError::UnsupportedFormat;
    assert_eq!(format!("{}", e_unsupported), "Unsupported RAW format or feature");

    let e_fail = PiexError::GenericFail("General failure".to_string());
    assert_eq!(format!("{}", e_fail), "Extraction failed: General failure");

    let io_err = std::io::Error::new(ErrorKind::UnexpectedEof, "EOF reached");
    let piex_io: PiexError = io_err.into();
    assert_eq!(format!("{}", piex_io), "I/O failure: EOF reached");
}

#[test]
fn test_byte_slice_stream_empty_slice_boundaries() {
    let empty_data: &[u8] = &[];
    let stream = ByteSliceStream::new(empty_data);
    assert_eq!(stream.len(), 0);
    assert!(stream.is_empty());

    // offset = 0, len = 0 (valid: offset <= len && len <= len - offset)
    let sub = stream.sub_stream(0, 0).expect("0,0 on empty should succeed");
    assert_eq!(sub.len(), 0);
    assert!(sub.is_empty());

    // offset = 0, len = 1 (OOB)
    assert!(matches!(
        stream.sub_stream(0, 1),
        Err(PiexError::OutOfBounds { requested_offset: 0, requested_size: 1, total_len: 0 })
    ));

    // offset = 0, len = usize::MAX (OOB)
    assert!(matches!(stream.sub_stream(0, usize::MAX), Err(PiexError::OutOfBounds { .. })));

    // offset = 1, len = 0 (offset > len -> OOB)
    assert!(matches!(
        stream.sub_stream(1, 0),
        Err(PiexError::OutOfBounds { requested_offset: 1, requested_size: 0, total_len: 0 })
    ));

    // offset = 1, len = 1 (OOB)
    assert!(matches!(stream.sub_stream(1, 1), Err(PiexError::OutOfBounds { .. })));

    // offset = usize::MAX, len = 0 (OOB)
    assert!(matches!(stream.sub_stream(usize::MAX, 0), Err(PiexError::OutOfBounds { .. })));

    // offset = usize::MAX, len = 1 (OOB)
    assert!(matches!(stream.sub_stream(usize::MAX, 1), Err(PiexError::OutOfBounds { .. })));

    // offset = usize::MAX, len = usize::MAX (OOB / arithmetic overflow)
    assert!(matches!(
        stream.sub_stream(usize::MAX, usize::MAX),
        Err(PiexError::OutOfBounds { .. })
    ));
}

#[test]
fn test_byte_slice_stream_single_byte_slice_boundaries() {
    let data: &[u8] = &[0x5a];
    let stream = ByteSliceStream::new(data);
    assert_eq!(stream.len(), 1);
    assert!(!stream.is_empty());

    // offset = 0, len = 0
    let sub = stream.sub_stream(0, 0).expect("0,0 on 1-byte slice should succeed");
    assert_eq!(sub.len(), 0);

    // offset = 0, len = 1
    let sub1 = stream.sub_stream(0, 1).expect("0,1 on 1-byte slice should succeed");
    assert_eq!(sub1.len(), 1);
    let mut buf = [0u8; 1];
    sub1.get_data(0, &mut buf).unwrap();
    assert_eq!(buf[0], 0x5a);

    // offset = 0, len = 2 (OOB)
    assert!(matches!(
        stream.sub_stream(0, 2),
        Err(PiexError::OutOfBounds { requested_offset: 0, requested_size: 2, total_len: 1 })
    ));

    // offset = 0, len = usize::MAX (OOB)
    assert!(matches!(stream.sub_stream(0, usize::MAX), Err(PiexError::OutOfBounds { .. })));

    // offset = 1 (offset == data.len()), len = 0 (valid empty sub-stream at end boundary)
    let sub_end = stream.sub_stream(1, 0).expect("1,0 at end boundary should succeed");
    assert_eq!(sub_end.len(), 0);

    // offset = 1 (offset == data.len()), len = 1 (OOB)
    assert!(matches!(
        stream.sub_stream(1, 1),
        Err(PiexError::OutOfBounds { requested_offset: 1, requested_size: 1, total_len: 1 })
    ));

    // offset = 2 (offset == data.len() + 1), len = 0 (offset > data.len() -> OOB)
    assert!(matches!(
        stream.sub_stream(2, 0),
        Err(PiexError::OutOfBounds { requested_offset: 2, requested_size: 0, total_len: 1 })
    ));

    // offset = 2, len = 1 (OOB)
    assert!(matches!(stream.sub_stream(2, 1), Err(PiexError::OutOfBounds { .. })));

    // offset = usize::MAX, len = 0
    assert!(matches!(stream.sub_stream(usize::MAX, 0), Err(PiexError::OutOfBounds { .. })));

    // offset = usize::MAX, len = usize::MAX
    assert!(matches!(
        stream.sub_stream(usize::MAX, usize::MAX),
        Err(PiexError::OutOfBounds { .. })
    ));
}

#[test]
fn test_byte_slice_stream_exhaustive_grid() {
    let raw_bytes: Vec<u8> = (0..32).collect();
    let stream = ByteSliceStream::new(&raw_bytes);
    let n = raw_bytes.len();

    // Exhaustive grid for offsets 0..=40 and lengths 0..=40
    for offset in 0..=40 {
        for len in 0..=40 {
            let res = stream.sub_stream(offset, len);
            // C++ SoT oracle logic:
            // Valid iff offset <= n AND len <= n - offset
            let expected_valid = offset <= n && len <= (n - offset);
            if expected_valid {
                assert!(
                    res.is_ok(),
                    "Expected Ok for offset={}, len={}, total={}, got {:?}",
                    offset,
                    len,
                    n,
                    res
                );
                let sub = res.unwrap();
                assert_eq!(sub.len(), len);
            } else {
                assert!(
                    res.is_err(),
                    "Expected Err for offset={}, len={}, total={}, got {:?}",
                    offset,
                    len,
                    n,
                    res
                );
            }
        }
    }
}

#[test]
fn test_byte_slice_stream_overflow_arithmetic() {
    let data = b"hello world";
    let stream = ByteSliceStream::new(data);

    let overflow_cases = [
        (usize::MAX, 0),
        (usize::MAX, 1),
        (usize::MAX, usize::MAX),
        (usize::MAX - 1, 2),
        (usize::MAX - 10, 11),
        (usize::MAX / 2 + 1, usize::MAX / 2 + 1),
        (1, usize::MAX),
        (0, usize::MAX),
        (data.len(), usize::MAX),
        (data.len() + 1, 0),
    ];

    for (offset, len) in overflow_cases {
        let res = stream.sub_stream(offset, len);
        assert!(res.is_err(), "Expected error for offset={}, len={}, got {:?}", offset, len, res);
    }
}

#[test]
fn test_stream_adapter_boundary_values() {
    // 1. Empty reader
    let empty_cursor = Cursor::new(Vec::<u8>::new());
    let adapter = StreamAdapter::new(empty_cursor).expect("Empty adapter creation");
    assert_eq!(adapter.len(), 0);
    assert!(adapter.is_empty());

    let mut buf0 = [];
    assert!(adapter.get_data(0, &mut buf0).is_ok());

    let mut buf1 = [0u8; 1];
    assert!(adapter.get_data(0, &mut buf1).is_err());
    assert!(adapter.get_data(1, &mut buf1).is_err());
    assert!(adapter.get_data(usize::MAX, &mut buf1).is_err());

    // 2. 1-byte reader
    let one_cursor = Cursor::new(vec![0x77]);
    let adapter1 = StreamAdapter::new(one_cursor).expect("1-byte adapter creation");
    assert_eq!(adapter1.len(), 1);

    let mut read_buf = [0u8; 1];
    adapter1.get_data(0, &mut read_buf).expect("Read 1 byte at 0");
    assert_eq!(read_buf[0], 0x77);

    let mut read_two = [0u8; 2];
    assert!(adapter1.get_data(0, &mut read_two).is_err());
    assert!(adapter1.get_data(1, &mut read_buf).is_err());
    assert!(adapter1.get_data(2, &mut read_buf).is_err());
    assert!(adapter1.get_data(usize::MAX, &mut read_buf).is_err());

    // Overflow requested size
    let mut overflow_buf = [0u8; 1];
    assert!(adapter1.get_data(usize::MAX, &mut overflow_buf).is_err());
}

#[test]
fn test_endian_reader_extreme_boundaries() {
    let data = [0x01, 0x02, 0x03, 0x04];

    // Out of bounds read_u8
    assert!(data.read_u8(0).is_ok());
    assert!(data.read_u8(3).is_ok());
    assert!(data.read_u8(4).is_err());
    assert!(data.read_u8(usize::MAX).is_err());

    // Out of bounds read_u16
    assert!(data.read_u16(0, Endian::Little).is_ok());
    assert!(data.read_u16(2, Endian::Little).is_ok());
    assert!(data.read_u16(3, Endian::Little).is_err());
    assert!(data.read_u16(4, Endian::Little).is_err());
    assert!(data.read_u16(usize::MAX, Endian::Little).is_err());

    // Out of bounds read_u32
    assert!(data.read_u32(0, Endian::Little).is_ok());
    assert!(data.read_u32(1, Endian::Little).is_err());
    assert!(data.read_u32(4, Endian::Little).is_err());
    assert!(data.read_u32(usize::MAX, Endian::Little).is_err());

    // Out of bounds read_rational (requires 8 bytes)
    assert!(data.read_rational(0, Endian::Little).is_err());
    assert!(data.read_rational(usize::MAX, Endian::Little).is_err());

    // Out of bounds read_bytes
    assert!(data.read_bytes(0, 0).is_ok());
    assert!(data.read_bytes(4, 0).is_ok());
    assert!(data.read_bytes(5, 0).is_err(), "offset 5 > len 4 on read_bytes(5, 0) must fail");
    assert!(data.read_bytes(usize::MAX, 0).is_err(), "offset usize::MAX on read_bytes must fail");
    assert!(data.read_bytes(0, 4).is_ok());
    assert!(data.read_bytes(0, 5).is_err());
    assert!(data.read_bytes(1, 4).is_err());
    assert!(data.read_bytes(usize::MAX, 1).is_err());
    assert!(data.read_bytes(usize::MAX - 1, 2).is_err());
}

#[test]
fn test_zero_length_reads_at_invalid_offsets() {
    let buf = [1u8, 2, 3, 4, 5];
    let slice_stream = ByteSliceStream::new(&buf);

    let mut empty_buf = [];
    assert!(slice_stream.get_data(0, &mut empty_buf).is_ok());
    assert!(slice_stream.get_data(5, &mut empty_buf).is_ok());
    assert!(slice_stream.get_data(6, &mut empty_buf).is_err());
    assert!(slice_stream.get_data(usize::MAX, &mut empty_buf).is_err());

    let cursor = Cursor::new(buf.to_vec());
    let adapter = StreamAdapter::new(cursor).expect("create adapter");
    assert!(adapter.get_data(0, &mut empty_buf).is_ok());
    assert!(adapter.get_data(5, &mut empty_buf).is_ok());
    assert!(adapter.get_data(6, &mut empty_buf).is_err());
    assert!(adapter.get_data(usize::MAX, &mut empty_buf).is_err());
}

#[test]
fn test_byte_slice_sub_stream_after_later_reads_clamps_bounds() {
    let data: [u8; 16] = core::array::from_fn(|i| 0xA0 + i as u8);
    let root = ByteSliceStream::new(&data);

    for &(prior_read_offset, sub_pos, sub_len) in &[
        (0usize, 6usize, 4usize),
        (0, 2, 6),
        (4, 5, 2),
        (4, 2, 8),
        (8, 2, 7),
        (8, 2, 6),
        (8, 2, 5),
        (12, 2, 4),
        (12, 2, 0),
    ] {
        let mut warm = [0u8; 4];
        root.get_data(prior_read_offset, &mut warm).unwrap();
        let sub = root.sub_stream(sub_pos, sub_len).unwrap();
        assert_eq!(sub.len(), sub_len);

        for i in 0..=(data.len() - sub_pos) {
            let mut probe = [0u8; 1];
            let res = sub.get_data(i, &mut probe);
            if i < sub_len {
                assert_eq!(res, Ok(()));
                assert_eq!(probe[0], data[sub_pos + i]);
            } else {
                assert!(matches!(res, Err(PiexError::OutOfBounds { .. })));
            }
        }
    }
}
