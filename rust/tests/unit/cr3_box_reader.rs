use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};
use std::io::Cursor;

use piex::cr3::{
    extract_cr3_orientation, extract_cr3_preview_data, find_next_box_with_tag, read_next_box,
    BoxHeader, BoxType, Uuid, UUID_MOOV, UUID_PRVW,
};
use piex::stream::{ByteSliceStream, StreamAdapter, StreamInterface};
use piex::types::PiexError;

// ============================================================================
// 1. BoxType Constants, Construction, Derives & Formatting Tests
// ============================================================================

#[test]
fn test_box_type_all_constants() {
    assert_eq!(BoxType::FTYP.as_bytes(), b"ftyp");
    assert_eq!(BoxType::MOOV.as_bytes(), b"moov");
    assert_eq!(BoxType::UUID.as_bytes(), b"uuid");
    assert_eq!(BoxType::PRVW.as_bytes(), b"PRVW");
    assert_eq!(BoxType::THMB.as_bytes(), b"THMB");
    assert_eq!(BoxType::CMT1.as_bytes(), b"CMT1");
    assert_eq!(BoxType::CMT2.as_bytes(), b"CMT2");
    assert_eq!(BoxType::TRAK.as_bytes(), b"trak");
    assert_eq!(BoxType::TKHD.as_bytes(), b"tkhd");
    assert_eq!(BoxType::MDIA.as_bytes(), b"mdia");
    assert_eq!(BoxType::MINF.as_bytes(), b"minf");
    assert_eq!(BoxType::STBL.as_bytes(), b"stbl");
    assert_eq!(BoxType::STSD.as_bytes(), b"stsd");
    assert_eq!(BoxType::CRAW.as_bytes(), b"CRAW");
    assert_eq!(BoxType::STSZ.as_bytes(), b"stsz");
    assert_eq!(BoxType::CO64.as_bytes(), b"co64");
    assert_eq!(BoxType::MDAT.as_bytes(), b"mdat");
}

#[test]
fn test_box_type_construction_and_derives() {
    let tag = BoxType::new(b"test");
    assert_eq!(tag.as_bytes(), b"test");
    assert_eq!(tag, BoxType([b't', b'e', b's', b't']));

    // Copy and Clone
    let tag_copy = tag;
    let tag_clone = tag.clone();
    assert_eq!(tag, tag_copy);
    assert_eq!(tag, tag_clone);

    // Hash trait
    let mut hasher1 = DefaultHasher::new();
    let mut hasher2 = DefaultHasher::new();
    tag.hash(&mut hasher1);
    tag_copy.hash(&mut hasher2);
    assert_eq!(hasher1.finish(), hasher2.finish());
}

#[test]
fn test_box_type_display_and_debug_formatting() {
    // Printable ASCII FourCCs
    let ftyp = BoxType::FTYP;
    assert_eq!(format!("{ftyp}"), "ftyp");
    assert_eq!(format!("{ftyp:?}"), "BoxType(ftyp)");

    let craw = BoxType::CRAW;
    assert_eq!(format!("{craw}"), "CRAW");
    assert_eq!(format!("{craw:?}"), "BoxType(CRAW)");

    let space_padded = BoxType::new(b"a b ");
    assert_eq!(format!("{space_padded}"), "a b ");
    assert_eq!(format!("{space_padded:?}"), "BoxType(a b )");

    // Non-ASCII or non-graphic FourCCs
    let non_ascii_zero = BoxType([0x00, 0x01, 0x02, 0x03]);
    assert_eq!(format!("{non_ascii_zero}"), "[0x00, 0x01, 0x02, 0x03]");
    assert_eq!(format!("{non_ascii_zero:?}"), "BoxType([0, 1, 2, 3])");

    let non_ascii_high = BoxType([0xFF, 0xFE, 0xFD, 0xFC]);
    assert_eq!(format!("{non_ascii_high}"), "[0xff, 0xfe, 0xfd, 0xfc]");
    assert_eq!(format!("{non_ascii_high:?}"), "BoxType([255, 254, 253, 252])");

    let non_graphic_control = BoxType([0x1B, 0x07, 0x08, 0x7F]);
    assert_eq!(format!("{non_graphic_control}"), "[0x1b, 0x07, 0x08, 0x7f]");
}

#[test]
fn test_uuid_constants_verification() {
    let expected_moov: Uuid = [
        0x85, 0xc0, 0xb6, 0x87, 0x82, 0x0f, 0x11, 0xe0, 0x81, 0x11, 0xf4, 0xce, 0x46, 0x2b, 0x6a,
        0x48,
    ];
    let expected_prvw: Uuid = [
        0xea, 0xf4, 0x2b, 0x5e, 0x1c, 0x98, 0x4b, 0x88, 0xb9, 0xfb, 0xb7, 0xdc, 0x40, 0x6e, 0x4d,
        0x16,
    ];

    assert_eq!(UUID_MOOV, expected_moov);
    assert_eq!(UUID_PRVW, expected_prvw);
}

// ============================================================================
// 2. BoxHeader Struct Methods & Boundary Checks
// ============================================================================

#[test]
fn test_box_header_methods_and_validity() {
    // Standard 32-bit box header
    let header = BoxHeader::new(BoxType::FTYP, 0, 8, 24);
    assert_eq!(header.header_length(), 8);
    assert_eq!(header.total_length(), 24);
    assert_eq!(header.payload_length(), 16);
    assert!(header.is_valid());

    // Extended 64-bit box header
    let ext_header = BoxHeader::new(BoxType::MDAT, 100, 116, 200);
    assert_eq!(ext_header.header_length(), 16);
    assert_eq!(ext_header.total_length(), 100);
    assert_eq!(ext_header.payload_length(), 84);
    assert!(ext_header.is_valid());

    // Zero-payload box header
    let empty_header = BoxHeader::new(BoxType::MOOV, 50, 58, 58);
    assert_eq!(empty_header.header_length(), 8);
    assert_eq!(empty_header.total_length(), 8);
    assert_eq!(empty_header.payload_length(), 0);
    assert!(empty_header.is_valid());

    // Invalid box headers
    let invalid_next_less_than_offset = BoxHeader::new(BoxType::FTYP, 10, 18, 5);
    assert!(!invalid_next_less_than_offset.is_valid());

    let invalid_header_less_than_offset = BoxHeader::new(BoxType::FTYP, 10, 8, 20);
    assert!(!invalid_header_less_than_offset.is_valid());

    let invalid_next_less_than_header = BoxHeader::new(BoxType::FTYP, 10, 20, 15);
    assert!(!invalid_next_less_than_header.is_valid());

    let invalid_zero_length = BoxHeader::new(BoxType::FTYP, 10, 10, 10);
    assert!(!invalid_zero_length.is_valid());
}

// ============================================================================
// 3. Compact 32-bit Box Header Reading Tests
// ============================================================================

#[test]
fn test_read_next_box_compact_header_various_sizes() {
    let sizes: &[u32] = &[8, 16, 32, 100, 1024];

    for &size in sizes {
        let mut buf = Vec::new();
        buf.extend_from_slice(&size.to_be_bytes());
        buf.extend_from_slice(b"ftyp");
        buf.resize(size as usize, 0u8);

        let mut stream = ByteSliceStream::new(&buf);
        let header = read_next_box(&mut stream, 0).unwrap().unwrap();

        assert_eq!(header.box_type, BoxType::FTYP);
        assert_eq!(header.offset, 0);
        assert_eq!(header.header_offset, 8);
        assert_eq!(header.next_box_offset, u64::from(size));
        assert_eq!(header.header_length(), 8);
        assert_eq!(header.total_length(), u64::from(size));
        assert_eq!(header.payload_length(), u64::from(size) - 8);
        assert!(header.is_valid());
    }
}

#[test]
fn test_read_next_box_compact_header_at_nonzero_offset() {
    let mut buf = vec![0xAAu8; 100]; // Prefix padding
    let box_offset = 100u64;
    let box_size = 24u32;

    buf.extend_from_slice(&box_size.to_be_bytes());
    buf.extend_from_slice(b"moov");
    buf.resize(buf.len() + 16, 0u8);

    let mut stream = ByteSliceStream::new(&buf);
    let header = read_next_box(&mut stream, box_offset).unwrap().unwrap();

    assert_eq!(header.box_type, BoxType::MOOV);
    assert_eq!(header.offset, 100);
    assert_eq!(header.header_offset, 108);
    assert_eq!(header.next_box_offset, 124);
    assert_eq!(header.header_length(), 8);
    assert_eq!(header.total_length(), 24);
    assert_eq!(header.payload_length(), 16);
}

// ============================================================================
// 4. Extended 64-bit Box Header Reading Tests
// ============================================================================

#[test]
fn test_read_next_box_extended_64bit_header() {
    let mut buf = Vec::new();
    // 32-bit length = 1 triggers 64-bit size reading
    buf.extend_from_slice(&1u32.to_be_bytes());
    buf.extend_from_slice(b"mdat");
    // 64-bit size = 64
    buf.extend_from_slice(&64u64.to_be_bytes());
    buf.resize(64, 0u8);

    let mut stream = ByteSliceStream::new(&buf);
    let header = read_next_box(&mut stream, 0).unwrap().unwrap();

    assert_eq!(header.box_type, BoxType::MDAT);
    assert_eq!(header.offset, 0);
    assert_eq!(header.header_offset, 16);
    assert_eq!(header.next_box_offset, 64);
    assert_eq!(header.header_length(), 16);
    assert_eq!(header.total_length(), 64);
    assert_eq!(header.payload_length(), 48);
}

#[test]
fn test_read_next_box_extended_header_at_nonzero_offset() {
    let mut buf = vec![0xFFu8; 50];
    let offset = 50u64;

    buf.extend_from_slice(&1u32.to_be_bytes());
    buf.extend_from_slice(b"uuid");
    buf.extend_from_slice(&100u64.to_be_bytes()); // total size 100
    buf.resize(150, 0u8);

    let mut stream = ByteSliceStream::new(&buf);
    let header = read_next_box(&mut stream, offset).unwrap().unwrap();

    assert_eq!(header.box_type, BoxType::UUID);
    assert_eq!(header.offset, 50);
    assert_eq!(header.header_offset, 66);
    assert_eq!(header.next_box_offset, 150);
    assert_eq!(header.header_length(), 16);
    assert_eq!(header.total_length(), 100);
    assert_eq!(header.payload_length(), 84);
}

// ============================================================================
// 5. Box Extending to EOF (`length_32 == 0`)
// ============================================================================

#[test]
fn test_read_next_box_length_zero_extends_to_eof() {
    let mut buf = Vec::new();
    // length_32 = 0 indicates box extends to the end of the stream
    buf.extend_from_slice(&0u32.to_be_bytes());
    buf.extend_from_slice(b"mdat");
    buf.extend_from_slice(&[0x12; 32]); // 40 bytes total stream length

    let mut stream = ByteSliceStream::new(&buf);
    let header = read_next_box(&mut stream, 0).unwrap().unwrap();

    assert_eq!(header.box_type, BoxType::MDAT);
    assert_eq!(header.offset, 0);
    assert_eq!(header.header_offset, 8);
    assert_eq!(header.next_box_offset, 40);
    assert_eq!(header.total_length(), 40);
    assert_eq!(header.payload_length(), 32);
}

#[test]
fn test_read_next_box_length_zero_at_nonzero_offset() {
    let mut buf = vec![0x00u8; 20];
    buf.extend_from_slice(&0u32.to_be_bytes());
    buf.extend_from_slice(b"trak");
    buf.extend_from_slice(&[0x34; 80]); // Total stream length 108

    let mut stream = ByteSliceStream::new(&buf);
    let header = read_next_box(&mut stream, 20).unwrap().unwrap();

    assert_eq!(header.box_type, BoxType::TRAK);
    assert_eq!(header.offset, 20);
    assert_eq!(header.header_offset, 28);
    assert_eq!(header.next_box_offset, 108);
    assert_eq!(header.total_length(), 88);
    assert_eq!(header.payload_length(), 80);
}

// ============================================================================
// 6. Clean EOF Handling
// ============================================================================

#[test]
fn test_read_next_box_clean_eof_handling() {
    let buf = vec![0u8; 16];
    let mut stream = ByteSliceStream::new(&buf);

    // Exact end of stream
    assert_eq!(read_next_box(&mut stream, 16).unwrap(), None);

    // Beyond end of stream
    assert_eq!(read_next_box(&mut stream, 17).unwrap(), None);
    assert_eq!(read_next_box(&mut stream, 1000).unwrap(), None);
    assert_eq!(read_next_box(&mut stream, u64::MAX - 10).unwrap(), None);
}

// ============================================================================
// 7. Error Handling on Truncated Streams & Malformed Sizes
// ============================================================================

#[test]
fn test_read_next_box_truncated_header_fewer_than_8_bytes() {
    for len in 0..8 {
        let buf = vec![0u8; len];
        let mut stream = ByteSliceStream::new(&buf);
        let res = read_next_box(&mut stream, 0);
        if len == 0 {
            assert_eq!(res.unwrap(), None);
        } else {
            assert!(
                matches!(res, Err(PiexError::InvalidBoxHeader { offset: 0 })),
                "Expected InvalidBoxHeader for length {len}"
            );
        }
    }
}

#[test]
fn test_read_next_box_truncated_extended_header_fewer_than_16_bytes() {
    for len in 8..16 {
        let mut buf = Vec::new();
        buf.extend_from_slice(&1u32.to_be_bytes()); // length_32 = 1
        buf.extend_from_slice(b"uuid");
        buf.resize(len, 0u8);

        let mut stream = ByteSliceStream::new(&buf);
        let res = read_next_box(&mut stream, 0);
        assert!(
            matches!(res, Err(PiexError::InvalidBoxHeader { offset: 0 })),
            "Expected InvalidBoxHeader for extended header with buffer length {len}"
        );
    }
}

#[test]
fn test_read_next_box_invalid_compact_sizes() {
    // Sizes 2..8 are invalid because header length is 8
    for size in 2u32..8u32 {
        let mut buf = Vec::new();
        buf.extend_from_slice(&size.to_be_bytes());
        buf.extend_from_slice(b"ftyp");
        buf.resize(16, 0u8);

        let mut stream = ByteSliceStream::new(&buf);
        let res = read_next_box(&mut stream, 0);
        assert!(
            matches!(res, Err(PiexError::InvalidBoxHeader { offset: 0 })),
            "Expected InvalidBoxHeader for invalid compact size {size}"
        );
    }
}

#[test]
fn test_read_next_box_invalid_extended_sizes() {
    // Sizes < 16 for extended headers are invalid because extended header is 16 bytes
    let invalid_sizes: &[u64] = &[0, 1, 8, 15];

    for &size in invalid_sizes {
        let mut buf = Vec::new();
        buf.extend_from_slice(&1u32.to_be_bytes());
        buf.extend_from_slice(b"uuid");
        buf.extend_from_slice(&size.to_be_bytes());
        buf.resize(32, 0u8);

        let mut stream = ByteSliceStream::new(&buf);
        let res = read_next_box(&mut stream, 0);
        assert!(
            matches!(res, Err(PiexError::InvalidBoxHeader { offset: 0 })),
            "Expected InvalidBoxHeader for invalid extended size {size}"
        );
    }
}

// ============================================================================
// 8. Bounded Scanning & `find_next_box_with_tag` Tests
// ============================================================================

#[test]
fn test_find_next_box_with_tag_sequential_children() {
    let mut buf = Vec::new();

    // Outer container: moov box (length = 56)
    buf.extend_from_slice(&56u32.to_be_bytes());
    buf.extend_from_slice(b"moov");

    // Child 1: ftyp (length = 16)
    buf.extend_from_slice(&16u32.to_be_bytes());
    buf.extend_from_slice(b"ftyp");
    buf.extend_from_slice(&[0u8; 8]);

    // Child 2: trak (length = 16)
    buf.extend_from_slice(&16u32.to_be_bytes());
    buf.extend_from_slice(b"trak");
    buf.extend_from_slice(&[0u8; 8]);

    // Child 3: uuid (length = 16)
    buf.extend_from_slice(&16u32.to_be_bytes());
    buf.extend_from_slice(b"uuid");
    buf.extend_from_slice(&[0u8; 8]);

    let mut stream = ByteSliceStream::new(&buf);
    let enclosing = read_next_box(&mut stream, 0).unwrap().unwrap();
    assert_eq!(enclosing.box_type, BoxType::MOOV);

    // 1. Search for first child (ftyp)
    let ftyp_box =
        find_next_box_with_tag(&mut stream, enclosing.header_offset, &enclosing, BoxType::FTYP)
            .unwrap()
            .unwrap();
    assert_eq!(ftyp_box.box_type, BoxType::FTYP);
    assert_eq!(ftyp_box.offset, 8);
    assert_eq!(ftyp_box.next_box_offset, 24);

    // 2. Search for middle child (trak) starting from enclosing header offset
    let trak_box =
        find_next_box_with_tag(&mut stream, enclosing.header_offset, &enclosing, BoxType::TRAK)
            .unwrap()
            .unwrap();
    assert_eq!(trak_box.box_type, BoxType::TRAK);
    assert_eq!(trak_box.offset, 24);
    assert_eq!(trak_box.next_box_offset, 40);

    // 3. Search for last child (uuid)
    let uuid_box =
        find_next_box_with_tag(&mut stream, enclosing.header_offset, &enclosing, BoxType::UUID)
            .unwrap()
            .unwrap();
    assert_eq!(uuid_box.box_type, BoxType::UUID);
    assert_eq!(uuid_box.offset, 40);
    assert_eq!(uuid_box.next_box_offset, 56);

    // 4. Search for non-existent child tag
    let not_found =
        find_next_box_with_tag(&mut stream, enclosing.header_offset, &enclosing, BoxType::MDAT)
            .unwrap();
    assert!(not_found.is_none());
}

#[test]
fn test_find_next_box_with_tag_enclosing_boundary_overrun() {
    let mut buf = Vec::new();

    // Outer container: moov box (length = 32)
    buf.extend_from_slice(&32u32.to_be_bytes());
    buf.extend_from_slice(b"moov");

    // Child box says length = 40, which extends past enclosing boundary (32)
    buf.extend_from_slice(&40u32.to_be_bytes());
    buf.extend_from_slice(b"trak");
    buf.resize(48, 0u8);

    let mut stream = ByteSliceStream::new(&buf);
    let enclosing = read_next_box(&mut stream, 0).unwrap().unwrap();

    let res =
        find_next_box_with_tag(&mut stream, enclosing.header_offset, &enclosing, BoxType::TRAK)
            .unwrap();
    // Expect Ok(None) because child box extends past enclosing box boundary
    assert!(res.is_none());
}

#[test]
fn test_find_next_box_with_tag_starting_at_or_after_boundary() {
    let mut buf = Vec::new();

    // Outer container: moov box (length = 24)
    buf.extend_from_slice(&24u32.to_be_bytes());
    buf.extend_from_slice(b"moov");
    buf.resize(24, 0u8);

    let mut stream = ByteSliceStream::new(&buf);
    let enclosing = read_next_box(&mut stream, 0).unwrap().unwrap();

    // Search starting at boundary
    let at_bound =
        find_next_box_with_tag(&mut stream, enclosing.next_box_offset, &enclosing, BoxType::FTYP)
            .unwrap();
    assert!(at_bound.is_none());

    // Search starting past boundary
    let past_bound = find_next_box_with_tag(
        &mut stream,
        enclosing.next_box_offset + 10,
        &enclosing,
        BoxType::FTYP,
    )
    .unwrap();
    assert!(past_bound.is_none());
}

#[test]
fn test_find_next_box_with_tag_detects_non_advancing_progress() {
    let mut buf = Vec::new();

    // Outer box: moov (length 40)
    buf.extend_from_slice(&40u32.to_be_bytes());
    buf.extend_from_slice(b"moov");

    // Corrupted child box with non-advancing box specification
    // Offset 8: child box has size 0, but total stream is 40 -> total_length becomes 32, but offset
    // doesn't advance if size is invalid Let's create invalid short length 4
    buf.extend_from_slice(&4u32.to_be_bytes());
    buf.extend_from_slice(b"bad_");
    buf.resize(40, 0u8);

    let mut stream = ByteSliceStream::new(&buf);
    let enclosing = read_next_box(&mut stream, 0).unwrap().unwrap();

    let res =
        find_next_box_with_tag(&mut stream, enclosing.header_offset, &enclosing, BoxType::TRAK);
    assert!(matches!(res, Err(PiexError::InvalidBoxHeader { .. })));
}

// ============================================================================
// 9. Multi-level Nested Container Hierarchy Test
// ============================================================================

#[test]
fn test_find_next_box_nested_hierarchy() {
    // Structure:
    // moov (outer, len 88)
    //   trak (child, len 80)
    //     mdia (child, len 72)
    //       minf (child, len 64)
    //         stbl (child, len 56)
    //           stsd (leaf, len 48)

    let mut buf = Vec::new();

    // moov: offset 0, len 88
    buf.extend_from_slice(&88u32.to_be_bytes());
    buf.extend_from_slice(b"moov");

    // trak: offset 8, len 80
    buf.extend_from_slice(&80u32.to_be_bytes());
    buf.extend_from_slice(b"trak");

    // mdia: offset 16, len 72
    buf.extend_from_slice(&72u32.to_be_bytes());
    buf.extend_from_slice(b"mdia");

    // minf: offset 24, len 64
    buf.extend_from_slice(&64u32.to_be_bytes());
    buf.extend_from_slice(b"minf");

    // stbl: offset 32, len 56
    buf.extend_from_slice(&56u32.to_be_bytes());
    buf.extend_from_slice(b"stbl");

    // stsd: offset 40, len 48
    buf.extend_from_slice(&48u32.to_be_bytes());
    buf.extend_from_slice(b"stsd");
    buf.resize(88, 0u8);

    let mut stream = ByteSliceStream::new(&buf);

    let moov = read_next_box(&mut stream, 0).unwrap().unwrap();
    assert_eq!(moov.box_type, BoxType::MOOV);

    let trak = find_next_box_with_tag(&mut stream, moov.header_offset, &moov, BoxType::TRAK)
        .unwrap()
        .unwrap();
    assert_eq!(trak.box_type, BoxType::TRAK);

    let mdia = find_next_box_with_tag(&mut stream, trak.header_offset, &trak, BoxType::MDIA)
        .unwrap()
        .unwrap();
    assert_eq!(mdia.box_type, BoxType::MDIA);

    let minf = find_next_box_with_tag(&mut stream, mdia.header_offset, &mdia, BoxType::MINF)
        .unwrap()
        .unwrap();
    assert_eq!(minf.box_type, BoxType::MINF);

    let stbl = find_next_box_with_tag(&mut stream, minf.header_offset, &minf, BoxType::STBL)
        .unwrap()
        .unwrap();
    assert_eq!(stbl.box_type, BoxType::STBL);

    let stsd = find_next_box_with_tag(&mut stream, stbl.header_offset, &stbl, BoxType::STSD)
        .unwrap()
        .unwrap();
    assert_eq!(stsd.box_type, BoxType::STSD);
    assert_eq!(stsd.offset, 40);
    assert_eq!(stsd.next_box_offset, 88);
}

// ============================================================================
// 10. StreamAdapter Compatibility Test
// ============================================================================

#[test]
fn test_box_reader_with_stream_adapter() {
    let mut buf = Vec::new();
    buf.extend_from_slice(&24u32.to_be_bytes());
    buf.extend_from_slice(b"CRAW");
    buf.extend_from_slice(&[0xFFu8; 16]);

    let cursor = Cursor::new(&buf);
    let mut adapter = StreamAdapter::new(cursor).unwrap();

    let header = read_next_box(&mut adapter, 0).unwrap().unwrap();
    assert_eq!(header.box_type, BoxType::CRAW);
    assert_eq!(header.offset, 0);
    assert_eq!(header.header_offset, 8);
    assert_eq!(header.next_box_offset, 24);
}

struct SparseBoxStream {
    chunks: Vec<(usize, Vec<u8>)>,
    total_len: usize,
}

impl StreamInterface for SparseBoxStream {
    fn get_data(&self, offset: usize, data: &mut [u8]) -> Result<(), PiexError> {
        let end = offset.checked_add(data.len()).ok_or(PiexError::OutOfBounds {
            requested_offset: offset,
            requested_size: data.len(),
            total_len: self.total_len,
        })?;
        for (base, bytes) in &self.chunks {
            if offset >= *base && end <= base.saturating_add(bytes.len()) {
                data.copy_from_slice(&bytes[offset - base..end - base]);
                return Ok(());
            }
        }
        Err(PiexError::OutOfBounds {
            requested_offset: offset,
            requested_size: data.len(),
            total_len: self.total_len,
        })
    }

    fn len(&self) -> usize {
        self.total_len
    }
}

fn ext_box_header(tag: &[u8; 4], len64: u64) -> Vec<u8> {
    [&1u32.to_be_bytes()[..], tag, &len64.to_be_bytes()].concat()
}

#[test]
fn test_read_next_box_extended_64bit_upper_32_bits_preserved() {
    let offset = 24u64;
    let valid_lengths = [
        u64::from(u32::MAX),
        1u64 << 32,
        (1u64 << 32) | 16,
        (0x1234_5678u64 << 32) | 0x9ABC_DEF0u64,
        u64::MAX - offset,
    ];
    for len64 in valid_lengths {
        let hdr = ext_box_header(b"mdat", len64);
        let stream =
            SparseBoxStream { chunks: vec![(offset as usize, hdr)], total_len: usize::MAX };
        let parsed = read_next_box(&stream, offset).unwrap().unwrap();
        assert_eq!(parsed.box_type, BoxType::MDAT);
        assert_eq!(parsed.offset, offset);
        assert_eq!(parsed.header_offset, offset + 16);
        assert_eq!(parsed.next_box_offset, offset + len64);
        assert_eq!(parsed.total_length(), len64);
    }

    for overflow_len64 in [u64::MAX - offset + 1, u64::MAX] {
        let hdr = ext_box_header(b"mdat", overflow_len64);
        let stream =
            SparseBoxStream { chunks: vec![(offset as usize, hdr)], total_len: usize::MAX };
        assert!(matches!(
            read_next_box(&stream, offset),
            Err(PiexError::InvalidBoxHeader { offset: 24 })
        ));
    }

    // Verify > 4 GiB sibling box traversal via find_next_box_with_tag.
    let child1_len = (1u64 << 32) + 16;
    let child2_off = 16 + child1_len;
    let moov_len = child2_off + 16;
    let stream = SparseBoxStream {
        chunks: vec![
            (0, ext_box_header(b"moov", moov_len)),
            (16, ext_box_header(b"free", child1_len)),
            (child2_off as usize, ext_box_header(b"trak", 16)),
        ],
        total_len: (moov_len as usize) + 16,
    };
    let moov = read_next_box(&stream, 0).unwrap().unwrap();
    let trak =
        find_next_box_with_tag(&stream, moov.header_offset, &moov, BoxType::TRAK).unwrap().unwrap();
    assert_eq!(trak.offset, child2_off);
    assert_eq!(trak.next_box_offset, moov_len);
}

// Same sparse file as `Cr3SkippedBoxWithMixed64BitLengthEndsWhereClaimed` in corpus_test.cc: a
// skipped box at offset 90 whose length ends on a 4 GiB multiple, then a 12-byte `mdat`.
#[test]
fn test_skipped_box_with_mixed_64bit_length_ends_where_claimed() {
    let far_len = 0x0102_0304_FFFF_FFA6u64;
    let head = [
        &12u32.to_be_bytes()[..],
        b"ftypcrx ",
        &ext_box_header(b"moov", 16),
        &34u32.to_be_bytes(),
        b"CMT1MM\0\x2A\0\0\0\x08\0\x01\x01\x12\0\x03\0\0\0\x01\0\x06\0\0\0\0\0\0",
        &28u32.to_be_bytes(),
        b"PRVW\0\0\0\0\0\0\0\xA0\0\x78\0\0\0\0\0\x04\xFF\xD8\xFF\xD9",
        &ext_box_header(b"free", far_len),
    ]
    .concat();
    assert_eq!(head.len(), 106);
    let free_end = 90 + far_len;
    let mdat = [&12u32.to_be_bytes()[..], b"mdat\0\0\0\0"].concat();
    for (mdat_offset, preview_found) in
        [(free_end - 1, false), (free_end, true), (free_end + 1, false)]
    {
        let stream = SparseBoxStream {
            chunks: vec![(0, head.clone()), (mdat_offset as usize, mdat.clone())],
            total_len: usize::MAX,
        };
        assert_eq!(extract_cr3_orientation(&stream).unwrap(), 6, "mdat at {mdat_offset}");
        let preview =
            extract_cr3_preview_data(&stream).map(|d| (d.preview.offset, d.preview.length));
        assert_eq!(preview.ok(), preview_found.then_some((86, 4)), "mdat at {mdat_offset}");
    }
}
