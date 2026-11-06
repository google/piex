#![forbid(unsafe_code)]

//! ISOBMFF box header parsing and tag scanning utilities.

use std::fmt;

use crate::stream::StreamInterface;
use crate::types::PiexError;

/// Represents a 4-byte ISOBMFF box tag (e.g. `ftyp`, `moov`, `uuid`).
#[repr(transparent)]
#[derive(Copy, Clone, PartialEq, Eq, Hash)]
pub struct BoxType(pub [u8; 4]);

impl BoxType {
    /// File Type box tag (`ftyp`).
    pub const FTYP: BoxType = BoxType(*b"ftyp");
    /// Movie box tag (`moov`).
    pub const MOOV: BoxType = BoxType(*b"moov");
    /// Extended User Data box tag (`uuid`).
    pub const UUID: BoxType = BoxType(*b"uuid");
    /// Preview image box tag (`PRVW`).
    pub const PRVW: BoxType = BoxType(*b"PRVW");
    /// Thumbnail image box tag (`THMB`).
    pub const THMB: BoxType = BoxType(*b"THMB");
    /// Comment 1 metadata box tag (`CMT1`).
    pub const CMT1: BoxType = BoxType(*b"CMT1");
    /// Comment 2 metadata box tag (`CMT2`).
    pub const CMT2: BoxType = BoxType(*b"CMT2");
    /// Track box tag (`trak`).
    pub const TRAK: BoxType = BoxType(*b"trak");
    /// Track Header box tag (`tkhd`).
    pub const TKHD: BoxType = BoxType(*b"tkhd");
    /// Media box tag (`mdia`).
    pub const MDIA: BoxType = BoxType(*b"mdia");
    /// Media Information box tag (`minf`).
    pub const MINF: BoxType = BoxType(*b"minf");
    /// Sample Table box tag (`stbl`).
    pub const STBL: BoxType = BoxType(*b"stbl");
    /// Sample Description box tag (`stsd`).
    pub const STSD: BoxType = BoxType(*b"stsd");
    /// Canon RAW box tag (`CRAW`).
    pub const CRAW: BoxType = BoxType(*b"CRAW");
    /// Sample Size box tag (`stsz`).
    pub const STSZ: BoxType = BoxType(*b"stsz");
    /// Chunk Offset 64-bit box tag (`co64`).
    pub const CO64: BoxType = BoxType(*b"co64");
    /// Media Data box tag (`mdat`).
    pub const MDAT: BoxType = BoxType(*b"mdat");

    /// Creates a new `BoxType` from a 4-byte array reference.
    pub const fn new(s: &[u8; 4]) -> Self {
        Self(*s)
    }

    /// Returns a reference to the underlying 4-byte array.
    pub const fn as_bytes(&self) -> &[u8; 4] {
        &self.0
    }
}

impl fmt::Display for BoxType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.0.iter().all(|&b| b.is_ascii_graphic() || b == b' ') {
            if let Ok(s) = std::str::from_utf8(&self.0) {
                return write!(f, "{s}");
            }
        }
        write!(
            f,
            "[{:#04x}, {:#04x}, {:#04x}, {:#04x}]",
            self.0[0], self.0[1], self.0[2], self.0[3]
        )
    }
}

impl fmt::Debug for BoxType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.0.iter().all(|&b| b.is_ascii_graphic() || b == b' ') {
            if let Ok(s) = std::str::from_utf8(&self.0) {
                return write!(f, "BoxType({s})");
            }
        }
        write!(f, "BoxType({:?})", self.0)
    }
}

/// 16-byte UUID for extended ISOBMFF box headers.
pub type Uuid = [u8; 16];

/// UUID identifying the `moov` extended box.
pub const UUID_MOOV: Uuid = [
    0x85, 0xc0, 0xb6, 0x87, 0x82, 0x0f, 0x11, 0xe0, 0x81, 0x11, 0xf4, 0xce, 0x46, 0x2b, 0x6a, 0x48,
];

/// UUID identifying the `PRVW` extended box.
pub const UUID_PRVW: Uuid = [
    0xea, 0xf4, 0x2b, 0x5e, 0x1c, 0x98, 0x4b, 0x88, 0xb9, 0xfb, 0xb7, 0xdc, 0x40, 0x6e, 0x4d, 0x16,
];

/// Represents the location and header parameters of an ISOBMFF box.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BoxHeader {
    /// Four-character box tag.
    pub box_type: BoxType,
    /// Absolute starting byte offset of the box in the stream.
    pub offset: u64,
    /// Absolute byte offset where the box payload starts.
    pub header_offset: u64,
    /// Absolute byte offset where the next box starts.
    pub next_box_offset: u64,
}

impl BoxHeader {
    /// Creates a new `BoxHeader` with specified attributes.
    pub const fn new(
        box_type: BoxType,
        offset: u64,
        header_offset: u64,
        next_box_offset: u64,
    ) -> Self {
        Self { box_type, offset, header_offset, next_box_offset }
    }

    /// Returns the length of the box header in bytes.
    pub const fn header_length(&self) -> u64 {
        self.header_offset.saturating_sub(self.offset)
    }

    /// Returns the total length of the box (header + payload) in bytes.
    pub const fn total_length(&self) -> u64 {
        self.next_box_offset.saturating_sub(self.offset)
    }

    /// Returns the payload length of the box in bytes.
    pub const fn payload_length(&self) -> u64 {
        self.next_box_offset.saturating_sub(self.header_offset)
    }

    /// Returns `true` if the box header boundaries are non-empty and well-ordered.
    pub fn is_valid(&self) -> bool {
        self.next_box_offset > self.offset
            && self.header_offset >= self.offset
            && self.next_box_offset >= self.header_offset
    }
}

/// Reads the next ISOBMFF box header from `stream` at the specified `offset`.
///
/// Returns `Ok(None)` cleanly when `offset` is at or beyond the stream boundary.
///
/// # Errors
/// Returns `PiexError::InvalidBoxHeader` if the box header is malformed, truncated,
/// or specifies invalid lengths or integer overflow.
pub fn read_next_box<R: StreamInterface + ?Sized>(
    stream: &R,
    offset: u64,
) -> Result<Option<BoxHeader>, PiexError> {
    let stream_len = stream.len();
    let offset_usize = match usize::try_from(offset) {
        Ok(val) => val,
        Err(_) => return Err(PiexError::InvalidBoxHeader { offset }),
    };

    if offset_usize >= stream_len {
        return Ok(None);
    }

    let end_header = match offset_usize.checked_add(8) {
        Some(end) => end,
        None => return Err(PiexError::InvalidBoxHeader { offset }),
    };

    if end_header > stream_len {
        return Err(PiexError::InvalidBoxHeader { offset });
    }

    let mut header_buf = [0u8; 8];
    stream
        .get_data(offset_usize, &mut header_buf)
        .map_err(|_| PiexError::InvalidBoxHeader { offset })?;

    let length_32 =
        u32::from_be_bytes([header_buf[0], header_buf[1], header_buf[2], header_buf[3]]);
    let box_type = BoxType([header_buf[4], header_buf[5], header_buf[6], header_buf[7]]);

    let (total_length, header_offset) = if length_32 == 1 {
        let ext_offset = match offset_usize.checked_add(8) {
            Some(val) => val,
            None => return Err(PiexError::InvalidBoxHeader { offset }),
        };
        let end_ext = match ext_offset.checked_add(8) {
            Some(val) => val,
            None => return Err(PiexError::InvalidBoxHeader { offset }),
        };
        if end_ext > stream_len {
            return Err(PiexError::InvalidBoxHeader { offset });
        }

        let mut len_buf = [0u8; 8];
        stream
            .get_data(ext_offset, &mut len_buf)
            .map_err(|_| PiexError::InvalidBoxHeader { offset })?;

        let length_64 = u64::from_be_bytes(len_buf);
        let header_off = match offset.checked_add(16) {
            Some(val) => val,
            None => return Err(PiexError::InvalidBoxHeader { offset }),
        };
        (length_64, header_off)
    } else {
        // C++ GetNextBox treats `length_32 == 0` as a literal 0-byte compact length,
        // which is then rejected by the `total_length < header_length` check below.
        let header_off = match offset.checked_add(8) {
            Some(val) => val,
            None => return Err(PiexError::InvalidBoxHeader { offset }),
        };
        (u64::from(length_32), header_off)
    };

    let header_length = match header_offset.checked_sub(offset) {
        Some(len) => len,
        None => return Err(PiexError::InvalidBoxHeader { offset }),
    };

    if total_length < header_length {
        return Err(PiexError::InvalidBoxHeader { offset });
    }

    let next_box_offset = match offset.checked_add(total_length) {
        Some(val) => val,
        None => return Err(PiexError::InvalidBoxHeader { offset }),
    };

    let box_header = BoxHeader::new(box_type, offset, header_offset, next_box_offset);
    if !box_header.is_valid() {
        return Err(PiexError::InvalidBoxHeader { offset });
    }

    Ok(Some(box_header))
}

/// Scans for the next child box matching `expected_tag` within `enclosing_box`.
///
/// Bounded scanning begins at `offset` and continues while `offset <
/// enclosing_box.next_box_offset`.
///
/// # Errors
/// Returns `PiexError::InvalidBoxHeader` if a child box header is corrupt or lacks monotonic
/// progress.
pub fn find_next_box_with_tag<R: StreamInterface + ?Sized>(
    stream: &R,
    mut offset: u64,
    enclosing_box: &BoxHeader,
    expected_tag: BoxType,
) -> Result<Option<BoxHeader>, PiexError> {
    while offset < enclosing_box.next_box_offset {
        let box_header = match read_next_box(stream, offset)? {
            Some(header) => header,
            None => return Ok(None),
        };

        if box_header.next_box_offset <= offset {
            return Err(PiexError::InvalidBoxHeader { offset });
        }

        if box_header.next_box_offset > enclosing_box.next_box_offset {
            return Ok(None);
        }

        if box_header.box_type == expected_tag {
            return Ok(Some(box_header));
        }

        offset = box_header.next_box_offset;
    }

    Ok(None)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::stream::ByteSliceStream;

    #[test]
    fn test_box_type_constants_and_display() {
        assert_eq!(BoxType::FTYP.as_bytes(), b"ftyp");
        assert_eq!(BoxType::MOOV.as_bytes(), b"moov");
        assert_eq!(BoxType::UUID.as_bytes(), b"uuid");
        assert_eq!(BoxType::PRVW.as_bytes(), b"PRVW");

        let b = BoxType::new(b"ftyp");
        assert_eq!(format!("{b}"), "ftyp");
        assert_eq!(format!("{b:?}"), "BoxType(ftyp)");

        let non_ascii = BoxType([0xff, 0x00, 0x01, 0x02]);
        assert_eq!(format!("{non_ascii}"), "[0xff, 0x00, 0x01, 0x02]");
    }

    #[test]
    fn test_box_header_methods() {
        let header = BoxHeader::new(BoxType::FTYP, 0, 8, 24);
        assert_eq!(header.header_length(), 8);
        assert_eq!(header.total_length(), 24);
        assert_eq!(header.payload_length(), 16);
        assert!(header.is_valid());

        let invalid = BoxHeader::new(BoxType::FTYP, 10, 8, 5);
        assert!(!invalid.is_valid());
    }

    #[test]
    fn test_read_next_box_standard_length() {
        let mut buf = Vec::new();
        // 32-bit length = 16
        buf.extend_from_slice(&16u32.to_be_bytes());
        buf.extend_from_slice(b"ftyp");
        buf.extend_from_slice(&[0u8; 8]);

        let mut stream = ByteSliceStream::new(&buf);
        let header = read_next_box(&mut stream, 0).unwrap().unwrap();

        assert_eq!(header.box_type, BoxType::FTYP);
        assert_eq!(header.offset, 0);
        assert_eq!(header.header_offset, 8);
        assert_eq!(header.next_box_offset, 16);
        assert_eq!(header.total_length(), 16);
        assert_eq!(header.payload_length(), 8);
    }

    #[test]
    fn test_read_next_box_extended_length() {
        let mut buf = Vec::new();
        // 32-bit length = 1 (indicates 64-bit extended size)
        buf.extend_from_slice(&1u32.to_be_bytes());
        buf.extend_from_slice(b"mdat");
        // 64-bit length = 32
        buf.extend_from_slice(&32u64.to_be_bytes());
        buf.extend_from_slice(&[0u8; 16]);

        let mut stream = ByteSliceStream::new(&buf);
        let header = read_next_box(&mut stream, 0).unwrap().unwrap();

        assert_eq!(header.box_type, BoxType::MDAT);
        assert_eq!(header.offset, 0);
        assert_eq!(header.header_offset, 16);
        assert_eq!(header.next_box_offset, 32);
        assert_eq!(header.total_length(), 32);
        assert_eq!(header.payload_length(), 16);
    }

    #[test]
    fn test_read_next_box_eof() {
        let buf = [0u8; 8];
        let mut stream = ByteSliceStream::new(&buf);
        let res = read_next_box(&mut stream, 8).unwrap();
        assert!(res.is_none());
    }

    #[test]
    fn test_read_next_box_invalid_short_length() {
        let mut buf = Vec::new();
        // 32-bit length = 4 (less than header size 8)
        buf.extend_from_slice(&4u32.to_be_bytes());
        buf.extend_from_slice(b"ftyp");

        let mut stream = ByteSliceStream::new(&buf);
        let res = read_next_box(&mut stream, 0);
        assert!(matches!(res, Err(PiexError::InvalidBoxHeader { offset: 0 })));
    }

    #[test]
    fn test_find_next_box_with_tag() {
        let mut buf = Vec::new();
        // Outer box: moov (length 40)
        buf.extend_from_slice(&40u32.to_be_bytes());
        buf.extend_from_slice(b"moov");

        // Child box 1: ftyp (length 16)
        buf.extend_from_slice(&16u32.to_be_bytes());
        buf.extend_from_slice(b"ftyp");
        buf.extend_from_slice(&[0u8; 8]);

        // Child box 2: trak (length 16)
        buf.extend_from_slice(&16u32.to_be_bytes());
        buf.extend_from_slice(b"trak");
        buf.extend_from_slice(&[0u8; 8]);

        let mut stream = ByteSliceStream::new(&buf);
        let enclosing = read_next_box(&mut stream, 0).unwrap().unwrap();
        assert_eq!(enclosing.box_type, BoxType::MOOV);

        let found =
            find_next_box_with_tag(&mut stream, enclosing.header_offset, &enclosing, BoxType::TRAK)
                .unwrap()
                .unwrap();

        assert_eq!(found.box_type, BoxType::TRAK);
        assert_eq!(found.offset, 24);
        assert_eq!(found.next_box_offset, 40);

        // Searching for non-existent tag
        let not_found =
            find_next_box_with_tag(&mut stream, enclosing.header_offset, &enclosing, BoxType::MDAT)
                .unwrap();
        assert!(not_found.is_none());
    }
}
