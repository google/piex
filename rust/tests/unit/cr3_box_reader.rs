use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};
use std::io::Cursor;

use piex::cr3::{
    find_next_box_with_tag, read_next_box, BoxHeader, BoxType, Uuid, UUID_MOOV, UUID_PRVW,
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

// ============================================================================
// 11. 64-bit Extended Box Bit-Lanes, Overflow Lattice & Multi-GiB CR3 Streams
// ============================================================================

struct SparseVirtualStream {
    segments: Vec<(usize, Vec<u8>)>,
    logical_len: usize,
}

impl StreamInterface for SparseVirtualStream {
    fn get_data(&self, offset: usize, data: &mut [u8]) -> Result<(), PiexError> {
        if data.is_empty() {
            return if offset <= self.logical_len {
                Ok(())
            } else {
                Err(PiexError::OutOfBounds {
                    requested_offset: offset,
                    requested_size: 0,
                    total_len: self.logical_len,
                })
            };
        }
        let end = offset.checked_add(data.len()).ok_or(PiexError::OutOfBounds {
            requested_offset: offset,
            requested_size: data.len(),
            total_len: self.logical_len,
        })?;
        if end > self.logical_len {
            return Err(PiexError::OutOfBounds {
                requested_offset: offset,
                requested_size: data.len(),
                total_len: self.logical_len,
            });
        }
        data.fill(0);
        for (seg_off, seg_bytes) in &self.segments {
            let seg_end = seg_off.saturating_add(seg_bytes.len());
            let copy_start = offset.max(*seg_off);
            let copy_end = end.min(seg_end);
            if copy_start < copy_end {
                data[copy_start - offset..copy_end - offset]
                    .copy_from_slice(&seg_bytes[copy_start - seg_off..copy_end - seg_off]);
            }
        }
        Ok(())
    }

    fn len(&self) -> usize {
        self.logical_len
    }
}

unsafe extern "C" fn sparse_ffi_read(
    user_data: *mut core::ffi::c_void,
    offset: usize,
    length: usize,
    out_data: *mut u8,
) -> i32 {
    if user_data.is_null() || (length > 0 && out_data.is_null()) {
        return -1;
    }
    // SAFETY: `user_data` points to a live `SparseVirtualStream` for the duration of each FFI call,
    // and when `length > 0`, `out_data` points to a writable buffer of `length` bytes allocated by
    // `CStreamAdapter` (when `length == 0`, we pass `&mut []` without dereferencing `out_data`).
    let (stream, dst): (&SparseVirtualStream, &mut [u8]) = unsafe {
        (
            &*(user_data as *const SparseVirtualStream),
            if length == 0 { &mut [] } else { core::slice::from_raw_parts_mut(out_data, length) },
        )
    };
    if stream.get_data(offset, dst).is_ok() {
        0
    } else {
        -1
    }
}

fn encode_extended_header(tag: &[u8; 4], len64: u64) -> [u8; 16] {
    let mut hdr = [0u8; 16];
    hdr[0..4].copy_from_slice(&1u32.to_be_bytes());
    hdr[4..8].copy_from_slice(tag);
    hdr[8..16].copy_from_slice(&len64.to_be_bytes());
    hdr
}

fn encode_compact_box(tag: &[u8; 4], payload: &[u8]) -> Vec<u8> {
    let total = u32::try_from(8 + payload.len()).unwrap();
    let mut out = Vec::with_capacity(8 + payload.len());
    out.extend_from_slice(&total.to_be_bytes());
    out.extend_from_slice(tag);
    out.extend_from_slice(payload);
    out
}

#[test]
fn test_read_next_box_extended_64bit_all_bit_and_byte_lanes() {
    let base_offset = 24usize;
    let mut sizes = Vec::new();
    for bit in 0..4 {
        sizes.push(16u64 | (1u64 << bit));
    }
    for bit in 4..64 {
        let candidate = 1u64 << bit;
        if candidate <= u64::MAX - (base_offset as u64) {
            sizes.push(candidate);
        }
    }
    for lo in 0u64..=15u64 {
        sizes.push((1u64 << 32) | lo);
    }
    sizes.extend_from_slice(&[
        0x0000_0000_8000_0010,
        0x8000_0000_0000_0000,
        0x8000_0000_8000_0010,
        0x0102_0304_0506_0708,
        0x8040_2010_8040_2010,
        0xFEDC_BA98_7654_3210,
        u64::from(u32::MAX),
        u64::from(u32::MAX) + 1,
    ]);

    for len64 in sizes {
        let hdr = encode_extended_header(b"mdat", len64);
        let stream = SparseVirtualStream {
            segments: vec![(base_offset, hdr.to_vec())],
            logical_len: base_offset + 64,
        };
        let parsed = read_next_box(&stream, base_offset as u64)
            .unwrap_or_else(|e| panic!("failed on len64={len64:#x}: {e:?}"))
            .expect("expected Some(BoxHeader)");
        assert_eq!(parsed.box_type, BoxType::MDAT);
        assert_eq!(parsed.offset, base_offset as u64);
        assert_eq!(parsed.header_offset, (base_offset as u64) + 16);
        assert_eq!(parsed.total_length(), len64, "len64={len64:#x}");
        assert_eq!(parsed.next_box_offset, (base_offset as u64) + len64);
        assert_eq!(parsed.payload_length(), len64 - 16);
    }

    // Skipped child box inside an enclosing container whose low 32 bits (16) fit inside
    // the container while its upper 32 bits overrun `enclosing.next_box_offset`.
    for hi_bit in [32u32, 47, 62, 63] {
        let overrun_len = (1u64 << hi_bit) | 16u64;
        let mut child_bytes = Vec::new();
        child_bytes.extend_from_slice(&encode_extended_header(b"free", overrun_len));
        child_bytes.extend_from_slice(&encode_compact_box(b"trak", &[0u8; 8]));
        let moov = encode_compact_box(b"moov", &child_bytes);

        let stream = ByteSliceStream::new(&moov);
        let enclosing = read_next_box(&stream, 0).unwrap().unwrap();
        let found =
            find_next_box_with_tag(&stream, enclosing.header_offset, &enclosing, BoxType::TRAK)
                .unwrap();
        assert!(found.is_none(), "expected overrun child at bit {hi_bit} to stop scan");
    }
}

#[test]
fn test_read_next_box_extended_64bit_overflow_and_max_offset_lattice() {
    let mut offsets = vec![0usize, 24, 512];
    if let Ok(near_4g) = usize::try_from(0xFFFF_FFE0u64) {
        offsets.push(near_4g);
    }
    if usize::MAX >= 16 {
        offsets.push(usize::MAX - 16);
    }

    for offset in offsets {
        let offset_u64 = offset as u64;
        let exact_max = u64::MAX - offset_u64;

        for short_len in [0u64, 1, 8, 15] {
            let hdr = encode_extended_header(b"uuid", short_len);
            let stream = SparseVirtualStream {
                segments: vec![(offset, hdr.to_vec())],
                logical_len: usize::MAX,
            };
            assert!(
                matches!(
                    read_next_box(&stream, offset_u64),
                    Err(PiexError::InvalidBoxHeader { offset: o }) if o == offset_u64
                ),
                "offset={offset:#x} short_len={short_len}"
            );
        }

        let mut valid_max_lattice = vec![16u64];
        if exact_max > 16 {
            valid_max_lattice.push(exact_max - 1);
        }
        if exact_max >= 16 {
            valid_max_lattice.push(exact_max);
        }
        for len64 in valid_max_lattice {
            let hdr = encode_extended_header(b"mdat", len64);
            let stream = SparseVirtualStream {
                segments: vec![(offset, hdr.to_vec())],
                logical_len: usize::MAX,
            };
            let parsed = read_next_box(&stream, offset_u64).unwrap().unwrap();
            assert_eq!(parsed.total_length(), len64);
            assert_eq!(parsed.next_box_offset, offset_u64 + len64);
        }

        if offset_u64 > 0 {
            let mut overflow_lattice = vec![u64::MAX];
            if let Some(plus_one) = exact_max.checked_add(1) {
                overflow_lattice.push(plus_one);
            }
            if u64::MAX - 1 > exact_max {
                overflow_lattice.push(u64::MAX - 1);
            }
            for len64 in overflow_lattice {
                let hdr = encode_extended_header(b"mdat", len64);
                let stream = SparseVirtualStream {
                    segments: vec![(offset, hdr.to_vec())],
                    logical_len: usize::MAX,
                };
                assert!(
                    matches!(
                        read_next_box(&stream, offset_u64),
                        Err(PiexError::InvalidBoxHeader { offset: o }) if o == offset_u64
                    ),
                    "offset={offset:#x} overflow len64={len64:#x}"
                );
            }
        }
    }

    let empty_at_max = SparseVirtualStream { segments: Vec::new(), logical_len: usize::MAX };
    assert_eq!(read_next_box(&empty_at_max, usize::MAX as u64).unwrap(), None);
    for tail_rem in 1usize..16usize {
        let off = usize::MAX - tail_rem;
        let hdr = encode_extended_header(b"mdat", 16);
        let stream = SparseVirtualStream {
            segments: vec![(off, hdr[..tail_rem].to_vec())],
            logical_len: usize::MAX,
        };
        assert!(matches!(
            read_next_box(&stream, off as u64),
            Err(PiexError::InvalidBoxHeader { offset: o }) if o == off as u64
        ));
    }
}

#[cfg(target_pointer_width = "64")]
#[test]
fn test_sparse_cr3_multi_gib_preview_and_ffi_extended_boxes() {
    use piex::c_types::{
        piex_error_t, piex_preview_image_data_t, piex_raw_image_type_t, piex_stream_t,
    };
    use piex::cr3::{extract_cr3_orientation, extract_cr3_preview_data};
    use piex::ffi::{
        piex_cr3_get_orientation, piex_cr3_get_preview_data, piex_get_preview_image_data,
    };

    let mut cmt1_tiff = Vec::new();
    cmt1_tiff.extend_from_slice(b"II\x2a\x00");
    cmt1_tiff.extend_from_slice(&8u32.to_le_bytes());
    cmt1_tiff.extend_from_slice(&1u16.to_le_bytes());
    cmt1_tiff.extend_from_slice(&0x0112u16.to_le_bytes());
    cmt1_tiff.extend_from_slice(&3u16.to_le_bytes());
    cmt1_tiff.extend_from_slice(&1u32.to_le_bytes());
    cmt1_tiff.extend_from_slice(&6u16.to_le_bytes());
    cmt1_tiff.extend_from_slice(&0u16.to_le_bytes());
    cmt1_tiff.extend_from_slice(&0u32.to_le_bytes());

    let mut uuid_moov_payload = UUID_MOOV.to_vec();
    uuid_moov_payload.extend_from_slice(&encode_compact_box(b"CMT1", &cmt1_tiff));
    let moov_children = encode_compact_box(b"uuid", &uuid_moov_payload);
    let mut moov = encode_extended_header(b"moov", (16 + moov_children.len()) as u64).to_vec();
    moov.extend_from_slice(&moov_children);

    let jpeg = b"\xFF\xD8\xFF\xE0\x00\x04PX\xFF\xD9";
    let mut prvw_payload = vec![0u8; 6];
    prvw_payload.extend_from_slice(&1620u16.to_be_bytes());
    prvw_payload.extend_from_slice(&1080u16.to_be_bytes());
    prvw_payload.extend_from_slice(&[0u8; 2]);
    prvw_payload.extend_from_slice(&(jpeg.len() as u32).to_be_bytes());
    prvw_payload.extend_from_slice(jpeg);
    let mut uuid_prvw_payload = UUID_PRVW.to_vec();
    uuid_prvw_payload.extend_from_slice(&[0u8; 8]);
    uuid_prvw_payload.extend_from_slice(&encode_compact_box(b"PRVW", &prvw_payload));
    let uuid_prvw = encode_compact_box(b"uuid", &uuid_prvw_payload);

    let mut prefix = encode_compact_box(b"ftyp", b"crx \x00\x00\x00\x01crx isom");
    prefix.extend_from_slice(&moov);
    prefix.extend_from_slice(&uuid_prvw);

    let free_offset = prefix.len();
    let skip_past_4gib: u64 = (1u64 << 32) | 64u64;
    let mut seg0 = prefix.clone();
    seg0.extend_from_slice(&encode_extended_header(b"free", skip_past_4gib));
    let real_mdat = encode_compact_box(b"mdat", &[0xAA; 16]);
    let target_mdat_off = free_offset + (skip_past_4gib as usize);

    // 1. Positive >4 GiB sparse seek lands cleanly on `mdat` across Rust and C-FFI.
    let mut valid_stream = SparseVirtualStream {
        segments: vec![(0, seg0.clone()), (target_mdat_off, real_mdat.clone())],
        logical_len: target_mdat_off + real_mdat.len(),
    };
    assert_eq!(extract_cr3_orientation(&valid_stream).unwrap(), 6);
    let preview = extract_cr3_preview_data(&valid_stream).unwrap();
    assert_eq!(preview.exif_orientation, 6);
    assert_eq!(preview.preview.width, 1620);
    assert_eq!(preview.preview.height, 1080);

    let c_stream = piex_stream_t {
        read_fn: Some(sparse_ffi_read),
        user_data: (&mut valid_stream as *mut SparseVirtualStream).cast(),
    };
    let mut ffi_orient = 0u32;
    let mut ffi_data = piex_preview_image_data_t::default();
    let mut ffi_type = piex_raw_image_type_t::PIEX_RAW_TYPE_NON_RAW;
    // SAFETY: `c_stream` wraps live `valid_stream`, and all output references are valid locals.
    unsafe {
        assert!(piex_cr3_get_orientation(Some(&c_stream), Some(&mut ffi_orient)));
        assert_eq!(ffi_orient, 6);
        assert_eq!(
            piex_cr3_get_preview_data(Some(&c_stream), Some(&mut ffi_data)),
            piex_error_t::PIEX_OK
        );
        assert_eq!(ffi_data.preview.width, 1620);
        assert_eq!(
            piex_get_preview_image_data(Some(&c_stream), Some(&mut ffi_data), Some(&mut ffi_type)),
            piex_error_t::PIEX_OK
        );
        assert_eq!(ffi_type, piex_raw_image_type_t::PIEX_RAW_TYPE_CR3);
    }

    // 2. Decoy `mdat` placed at low-32-bit offset `free_offset + 64` is never visited when the
    //    stream ends before `free_offset + ((1 << 32) | 64)`.
    let decoy_off = free_offset + 64;
    let mut decoy_stream = SparseVirtualStream {
        segments: vec![(0, seg0), (decoy_off, real_mdat)],
        logical_len: 16_384,
    };
    assert!(extract_cr3_preview_data(&decoy_stream).is_err());
    let decoy_c_stream = piex_stream_t {
        read_fn: Some(sparse_ffi_read),
        user_data: (&mut decoy_stream as *mut SparseVirtualStream).cast(),
    };
    // SAFETY: `decoy_c_stream` wraps live `decoy_stream` and `ffi_data` is a valid local.
    unsafe {
        assert_eq!(
            piex_cr3_get_preview_data(Some(&decoy_c_stream), Some(&mut ffi_data)),
            piex_error_t::PIEX_FAIL
        );
    }

    // 3. Overflowing `mdat` extended length (`u64::MAX`) at non-zero offset is rejected.
    let mut overflow_seg = prefix;
    overflow_seg.extend_from_slice(&encode_extended_header(b"mdat", u64::MAX));
    let mut overflow_stream =
        SparseVirtualStream { segments: vec![(0, overflow_seg)], logical_len: 16_384 };
    assert!(extract_cr3_preview_data(&overflow_stream).is_err());
    let overflow_c_stream = piex_stream_t {
        read_fn: Some(sparse_ffi_read),
        user_data: (&mut overflow_stream as *mut SparseVirtualStream).cast(),
    };
    // SAFETY: `overflow_c_stream` wraps live `overflow_stream` and `ffi_data` is a valid local.
    unsafe {
        assert_eq!(
            piex_cr3_get_preview_data(Some(&overflow_c_stream), Some(&mut ffi_data)),
            piex_error_t::PIEX_FAIL
        );
    }
}
