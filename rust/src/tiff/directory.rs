#![forbid(unsafe_code)]

//! TIFF directory parsing types and Image File Directory (IFD) structures.

use thiserror::Error;

use crate::stream::{Endian, StreamError};
use crate::types::{Rational, SRational};

/// Errors encountered during TIFF parsing and directory operations.
#[derive(Debug, Error, PartialEq, Eq, Clone)]
pub enum TiffError {
    /// Stream or I/O error occurred during TIFF parsing.
    #[error("Stream error: {0}")]
    StreamError(#[from] StreamError),

    /// Header byte order marker is invalid.
    #[error("Invalid byte order prefix")]
    InvalidByteOrder,

    /// TIFF magic number is invalid (expected 0x002A).
    #[error("Invalid TIFF magic number: {0:#06x}, expected 0x002a")]
    InvalidMagicNumber(u16),

    /// Requested tag was not found in the TIFF directory.
    #[error("Tag {tag:#06x} not found in TIFF directory")]
    TagNotFound {
        /// Tag ID.
        tag: u32,
    },

    /// Tag data type does not match expected type.
    #[error("Type mismatch for tag {tag:#06x}: expected {expected:?}, actual {actual:?}")]
    TypeMismatch {
        /// Tag ID.
        tag: u32,
        /// Expected TIFF type.
        expected: TiffType,
        /// Actual TIFF type.
        actual: TiffType,
    },

    /// Number of elements does not match expected count.
    #[error("Invalid count for tag {tag:#06x}: expected {expected}, actual {actual}")]
    InvalidCount {
        /// Tag ID.
        tag: u32,
        /// Expected count.
        expected: usize,
        /// Actual count.
        actual: usize,
    },

    /// Rational denominator is zero.
    #[error("Zero denominator encountered in rational for tag {tag:#06x}")]
    ZeroDenominator {
        /// Tag ID.
        tag: u32,
    },

    /// Maximum number of IFDs limit was reached.
    #[error("Maximum IFD limit reached: {0}")]
    MaxIfdLimitReached(u16),

    /// Directory entry is corrupted at offset.
    #[error("Corrupted directory entry at offset {offset:#010x}")]
    CorruptedEntry {
        /// Byte offset of corrupted entry.
        offset: u32,
    },

    /// Type ID in directory entry is invalid.
    #[error("Invalid TIFF type ID: {0}")]
    InvalidTypeId(u16),

    /// ASCII string tag contains invalid UTF-8 encoding.
    #[error("Invalid string encoding for tag {tag:#06x}")]
    InvalidStringEncoding {
        /// Tag ID.
        tag: u32,
    },
}

/// Standard TIFF data types.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(u16)]
pub enum TiffType {
    /// 8-bit unsigned integer.
    Byte = 1,
    /// 8-bit byte array containing a NUL-terminated ASCII string.
    Ascii = 2,
    /// 16-bit unsigned integer.
    Short = 3,
    /// 32-bit unsigned integer.
    Long = 4,
    /// Two 32-bit unsigned integers: numerator and denominator.
    Rational = 5,
    /// 8-bit signed integer.
    SByte = 6,
    /// 8-bit byte array containing undefined data.
    Undefined = 7,
    /// 16-bit signed integer.
    SShort = 8,
    /// 32-bit signed integer.
    SLong = 9,
    /// Two 32-bit signed integers: numerator and denominator.
    SRational = 10,
    /// 32-bit IEEE floating point.
    Float = 11,
    /// 64-bit IEEE floating point.
    Double = 12,
    /// 32-bit unsigned integer offset to sub-IFD.
    Ifd = 13,
}

impl TiffType {
    /// Returns the element size in bytes for a single element of this type.
    pub fn element_size(self) -> usize {
        match self {
            TiffType::Byte | TiffType::Ascii | TiffType::SByte | TiffType::Undefined => 1,
            TiffType::Short | TiffType::SShort => 2,
            TiffType::Long | TiffType::SLong | TiffType::Float | TiffType::Ifd => 4,
            TiffType::Rational | TiffType::SRational | TiffType::Double => 8,
        }
    }

    /// Converts a 16-bit integer type code into a `TiffType`.
    pub fn from_u16(val: u16) -> Option<Self> {
        match val {
            1 => Some(TiffType::Byte),
            2 => Some(TiffType::Ascii),
            3 => Some(TiffType::Short),
            4 => Some(TiffType::Long),
            5 => Some(TiffType::Rational),
            6 => Some(TiffType::SByte),
            7 => Some(TiffType::Undefined),
            8 => Some(TiffType::SShort),
            9 => Some(TiffType::SLong),
            10 => Some(TiffType::SRational),
            11 => Some(TiffType::Float),
            12 => Some(TiffType::Double),
            13 => Some(TiffType::Ifd),
            _ => None,
        }
    }
}

/// Representation of TIFF tag data: inline for small values (<= 8 bytes), heap for large arrays.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum EntryValue {
    /// Small data stored directly inside the enum (0 to 8 bytes) with zero heap allocation.
    Inline([u8; 8], u8),
    /// Data larger than 8 bytes stored on the heap.
    Heap(Vec<u8>),
}

impl EntryValue {
    /// Returns a slice over the entry bytes.
    #[inline]
    pub fn as_slice(&self) -> &[u8] {
        match self {
            EntryValue::Inline(buf, len) => &buf[..(*len as usize)],
            EntryValue::Heap(vec) => vec.as_slice(),
        }
    }

    /// Returns length in bytes.
    #[inline]
    pub fn len(&self) -> usize {
        match self {
            EntryValue::Inline(_, len) => *len as usize,
            EntryValue::Heap(vec) => vec.len(),
        }
    }

    /// Returns true if empty.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

impl std::ops::Deref for EntryValue {
    type Target = [u8];

    #[inline]
    fn deref(&self) -> &Self::Target {
        self.as_slice()
    }
}

/// Information for a single directory entry within an IFD.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DirectoryEntry {
    /// Data type of the tag values.
    pub type_: TiffType,
    /// Number of values of type_ (not byte count).
    pub count: u32,
    /// Byte offset of entry's data in stream, or 0 if contained in entry or unassigned.
    pub offset: u32,
    /// Byte representation of entry value(s) stored inline if <= 8 bytes.
    pub value: EntryValue,
}

impl DirectoryEntry {
    /// Creates a new `DirectoryEntry` instance from a `Vec<u8>`.
    pub fn new(type_: TiffType, count: u32, offset: u32, value: Vec<u8>) -> Self {
        let val = if value.len() <= 8 {
            let mut buf = [0u8; 8];
            buf[..value.len()].copy_from_slice(&value);
            EntryValue::Inline(buf, value.len() as u8)
        } else {
            EntryValue::Heap(value)
        };
        Self { type_, count, offset, value: val }
    }

    /// Creates a new `DirectoryEntry` instance from a byte slice,
    /// storing inline if `<= 8` bytes to avoid heap allocations.
    pub fn new_slice(type_: TiffType, count: u32, offset: u32, value: &[u8]) -> Self {
        let val = if value.len() <= 8 {
            let mut buf = [0u8; 8];
            buf[..value.len()].copy_from_slice(value);
            EntryValue::Inline(buf, value.len() as u8)
        } else {
            EntryValue::Heap(value.to_vec())
        };
        Self { type_, count, offset, value: val }
    }

    /// Returns a byte slice of the entry value.
    #[inline]
    pub fn value_bytes(&self) -> &[u8] {
        &self.value
    }
}

/// Represents a TIFF Image File Directory (IFD).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TiffDirectory {
    directory_entries: Vec<(u32, DirectoryEntry)>,
    tag_order: Vec<u32>,
    sub_directories: Vec<TiffDirectory>,
    endian: Endian,
}

impl TiffDirectory {
    /// Creates a new empty `TiffDirectory` with specified byte order.
    pub fn new(endian: Endian) -> Self {
        Self {
            directory_entries: Vec::new(),
            tag_order: Vec::new(),
            sub_directories: Vec::new(),
            endian,
        }
    }

    /// Returns the byte endianness for this directory.
    pub fn endian(&self) -> Endian {
        self.endian
    }

    fn find_entry(&self, tag: u32) -> Option<&DirectoryEntry> {
        self.directory_entries.iter().find(|(t, _)| *t == tag).map(|(_, e)| e)
    }

    /// Checks if tag is present in directory.
    pub fn has(&self, tag: u32) -> bool {
        self.find_entry(tag).is_some()
    }

    /// Gets raw byte slice for tag if type is `Byte` or `Undefined` without cloning.
    pub fn get_bytes_ref(&self, tag: u32) -> Result<&[u8], TiffError> {
        let entry = self.find_entry(tag).ok_or(TiffError::TagNotFound { tag })?;
        if entry.type_ != TiffType::Byte && entry.type_ != TiffType::Undefined {
            return Err(TiffError::TypeMismatch {
                tag,
                expected: TiffType::Byte,
                actual: entry.type_,
            });
        }
        Ok(&entry.value)
    }

    /// Gets raw byte value for tag if type is `Byte` or `Undefined`.
    pub fn get_bytes(&self, tag: u32) -> Result<Vec<u8>, TiffError> {
        let entry = self.find_entry(tag).ok_or(TiffError::TagNotFound { tag })?;
        if entry.type_ != TiffType::Byte && entry.type_ != TiffType::Undefined {
            return Err(TiffError::TypeMismatch {
                tag,
                expected: TiffType::Byte,
                actual: entry.type_,
            });
        }
        Ok(entry.value.to_vec())
    }

    /// Gets ASCII string value for tag.
    pub fn get_string(&self, tag: u32) -> Result<String, TiffError> {
        let entry = self.find_entry(tag).ok_or(TiffError::TagNotFound { tag })?;
        if entry.type_ != TiffType::Ascii {
            return Err(TiffError::TypeMismatch {
                tag,
                expected: TiffType::Ascii,
                actual: entry.type_,
            });
        }
        let val = &entry.value[..];
        let end = val.iter().position(|&b| b == 0).unwrap_or(val.len());
        Ok(String::from_utf8_lossy(&val[..end]).into_owned())
    }

    /// Gets the complete ASCII value bytes for tag, including trailing and
    /// embedded NUL bytes, exactly like C++ `TiffDirectory::Get(tag, std::string*)`.
    pub fn get_ascii_bytes(&self, tag: u32) -> Result<&[u8], TiffError> {
        let entry = self.find_entry(tag).ok_or(TiffError::TagNotFound { tag })?;
        if entry.type_ != TiffType::Ascii {
            return Err(TiffError::TypeMismatch {
                tag,
                expected: TiffType::Ascii,
                actual: entry.type_,
            });
        }
        Ok(&entry.value)
    }

    /// Gets single u32 value for tag (supports `Short` and `Long`).
    pub fn get_u32(&self, tag: u32) -> Result<u32, TiffError> {
        let entry = self.find_entry(tag).ok_or(TiffError::TagNotFound { tag })?;
        if entry.count != 1 {
            return Err(TiffError::InvalidCount { tag, expected: 1, actual: entry.count as usize });
        }
        let vec = self.get_u32_vec(tag)?;
        vec.first().copied().ok_or(TiffError::InvalidCount { tag, expected: 1, actual: 0 })
    }

    /// Gets vector of u32 values for tag (supports `Short` and `Long`).
    pub fn get_u32_vec(&self, tag: u32) -> Result<Vec<u32>, TiffError> {
        let entry = self.find_entry(tag).ok_or(TiffError::TagNotFound { tag })?;
        let elem_size = match entry.type_ {
            TiffType::Short => 2,
            TiffType::Long => 4,
            _ => {
                return Err(TiffError::TypeMismatch {
                    tag,
                    expected: TiffType::Long,
                    actual: entry.type_,
                })
            }
        };
        let required_bytes = (entry.count as usize).saturating_mul(elem_size);
        let val_bytes = &entry.value[..];
        if val_bytes.len() < required_bytes {
            return Err(TiffError::InvalidCount {
                tag,
                expected: entry.count as usize,
                actual: val_bytes.len() / elem_size,
            });
        }
        match entry.type_ {
            TiffType::Short => Ok(val_bytes[..required_bytes]
                .chunks_exact(2)
                .map(|b| match self.endian {
                    Endian::Little => u16::from_le_bytes([b[0], b[1]]) as u32,
                    Endian::Big => u16::from_be_bytes([b[0], b[1]]) as u32,
                })
                .collect()),
            TiffType::Long => Ok(val_bytes[..required_bytes]
                .chunks_exact(4)
                .map(|b| match self.endian {
                    Endian::Little => u32::from_le_bytes([b[0], b[1], b[2], b[3]]),
                    Endian::Big => u32::from_be_bytes([b[0], b[1], b[2], b[3]]),
                })
                .collect()),
            _ => unreachable!(),
        }
    }

    /// Gets single `Rational` value for tag (supports `Short`, `Long`, `Rational`).
    pub fn get_rational(&self, tag: u32) -> Result<Rational, TiffError> {
        let entry = self.find_entry(tag).ok_or(TiffError::TagNotFound { tag })?;
        if entry.count != 1 {
            return Err(TiffError::InvalidCount { tag, expected: 1, actual: entry.count as usize });
        }
        let vec = self.get_rational_vec(tag)?;
        vec.first().copied().ok_or(TiffError::InvalidCount { tag, expected: 1, actual: 0 })
    }

    /// Gets vector of `Rational` values for tag (supports `Short`, `Long`, `Rational`).
    pub fn get_rational_vec(&self, tag: u32) -> Result<Vec<Rational>, TiffError> {
        let entry = self.find_entry(tag).ok_or(TiffError::TagNotFound { tag })?;
        match entry.type_ {
            TiffType::Short | TiffType::Long => {
                let u32_vec = self.get_u32_vec(tag)?;
                Ok(u32_vec.into_iter().map(|v| Rational::new(v, 1)).collect())
            }
            TiffType::Rational => {
                let elem_size = 8;
                let required_bytes = (entry.count as usize).saturating_mul(elem_size);
                let val_bytes = &entry.value[..];
                if val_bytes.len() < required_bytes {
                    return Err(TiffError::InvalidCount {
                        tag,
                        expected: entry.count as usize,
                        actual: val_bytes.len() / elem_size,
                    });
                }
                val_bytes[..required_bytes]
                    .chunks_exact(8)
                    .map(|b| {
                        let (num, den) = match self.endian {
                            Endian::Little => (
                                u32::from_le_bytes([b[0], b[1], b[2], b[3]]),
                                u32::from_le_bytes([b[4], b[5], b[6], b[7]]),
                            ),
                            Endian::Big => (
                                u32::from_be_bytes([b[0], b[1], b[2], b[3]]),
                                u32::from_be_bytes([b[4], b[5], b[6], b[7]]),
                            ),
                        };
                        if den == 0 {
                            Err(TiffError::ZeroDenominator { tag })
                        } else {
                            Ok(Rational::new(num, den))
                        }
                    })
                    .collect()
            }
            _ => Err(TiffError::TypeMismatch {
                tag,
                expected: TiffType::Rational,
                actual: entry.type_,
            }),
        }
    }

    /// Gets single `SRational` value for tag (supports `SShort`, `SLong`, `SRational`).
    pub fn get_srational(&self, tag: u32) -> Result<SRational, TiffError> {
        let entry = self.find_entry(tag).ok_or(TiffError::TagNotFound { tag })?;
        if entry.count != 1 {
            return Err(TiffError::InvalidCount { tag, expected: 1, actual: entry.count as usize });
        }
        let vec = self.get_srational_vec(tag)?;
        vec.first().copied().ok_or(TiffError::InvalidCount { tag, expected: 1, actual: 0 })
    }

    /// Gets vector of `SRational` values for tag (supports `SShort`, `SLong`, `SRational`).
    pub fn get_srational_vec(&self, tag: u32) -> Result<Vec<SRational>, TiffError> {
        let entry = self.find_entry(tag).ok_or(TiffError::TagNotFound { tag })?;
        let elem_size = match entry.type_ {
            TiffType::SShort => 2,
            TiffType::SLong => 4,
            TiffType::SRational => 8,
            _ => {
                return Err(TiffError::TypeMismatch {
                    tag,
                    expected: TiffType::SRational,
                    actual: entry.type_,
                })
            }
        };
        let required_bytes = (entry.count as usize).saturating_mul(elem_size);
        let val_bytes = &entry.value[..];
        if val_bytes.len() < required_bytes {
            return Err(TiffError::InvalidCount {
                tag,
                expected: entry.count as usize,
                actual: val_bytes.len() / elem_size,
            });
        }
        match entry.type_ {
            TiffType::SShort => Ok(val_bytes[..required_bytes]
                .chunks_exact(2)
                .map(|b| {
                    let val = match self.endian {
                        Endian::Little => i16::from_le_bytes([b[0], b[1]]) as i32,
                        Endian::Big => i16::from_be_bytes([b[0], b[1]]) as i32,
                    };
                    SRational::new(val, 1)
                })
                .collect()),
            TiffType::SLong => Ok(val_bytes[..required_bytes]
                .chunks_exact(4)
                .map(|b| {
                    let val = match self.endian {
                        Endian::Little => i32::from_le_bytes([b[0], b[1], b[2], b[3]]),
                        Endian::Big => i32::from_be_bytes([b[0], b[1], b[2], b[3]]),
                    };
                    SRational::new(val, 1)
                })
                .collect()),
            TiffType::SRational => val_bytes[..required_bytes]
                .chunks_exact(8)
                .map(|b| {
                    let (num, den) = match self.endian {
                        Endian::Little => (
                            i32::from_le_bytes([b[0], b[1], b[2], b[3]]),
                            i32::from_le_bytes([b[4], b[5], b[6], b[7]]),
                        ),
                        Endian::Big => (
                            i32::from_be_bytes([b[0], b[1], b[2], b[3]]),
                            i32::from_be_bytes([b[4], b[5], b[6], b[7]]),
                        ),
                    };
                    if den == 0 {
                        Err(TiffError::ZeroDenominator { tag })
                    } else {
                        Ok(SRational::new(num, den))
                    }
                })
                .collect(),
            _ => unreachable!(),
        }
    }

    /// Gets offset and byte length of tag data if entry matches expected type.
    pub fn get_offset_and_length(
        &self,
        tag: u32,
        expected_type: TiffType,
    ) -> Result<(u32, u32), TiffError> {
        let entry = self.find_entry(tag).ok_or(TiffError::TagNotFound { tag })?;
        if entry.type_ != expected_type {
            return Err(TiffError::TypeMismatch {
                tag,
                expected: expected_type,
                actual: entry.type_,
            });
        }
        let length = entry.count.saturating_mul(expected_type.element_size() as u32);
        Ok((entry.offset, length))
    }

    /// Adds an entry to directory with given parameters.
    pub fn add_entry(
        &mut self,
        tag: u32,
        type_: TiffType,
        count: u32,
        offset: u32,
        value: Vec<u8>,
    ) {
        self.add_directory_entry(tag, DirectoryEntry::new(type_, count, offset, value));
    }

    /// Adds a `DirectoryEntry` to directory.
    pub fn add_directory_entry(&mut self, tag: u32, entry: DirectoryEntry) {
        if let Some(pos) = self.directory_entries.iter().position(|(t, _)| *t == tag) {
            self.directory_entries[pos] = (tag, entry);
        } else {
            self.tag_order.push(tag);
            self.directory_entries.push((tag, entry));
        }
    }

    /// Adds a sub-directory to this directory.
    pub fn add_sub_directory(&mut self, sub_directory: TiffDirectory) {
        self.sub_directories.push(sub_directory);
    }

    /// Returns slice of sub-directories contained in this directory.
    pub fn sub_directories(&self) -> &[TiffDirectory] {
        &self.sub_directories
    }
}

/// Container for all directories extracted from a TIFF file.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TiffContent {
    /// Main image file directories.
    pub tiff_directories: Vec<TiffDirectory>,
    /// EXIF sub-directory, if present.
    pub exif_directory: Option<TiffDirectory>,
    /// GPS sub-directory, if present.
    pub gps_directory: Option<TiffDirectory>,
}

impl TiffContent {
    /// Creates a new empty `TiffContent`.
    pub fn new() -> Self {
        Self { tiff_directories: Vec::new(), exif_directory: None, gps_directory: None }
    }
}

impl Default for TiffContent {
    fn default() -> Self {
        Self::new()
    }
}
