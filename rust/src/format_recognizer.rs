#![forbid(unsafe_code)]

//! RAW image format identification engine and format checkers.

use crate::stream::Endian;
use crate::types::RawImageType;

/// Maximum byte size inspected for format identification.
pub const MAX_HEADER_INSPECTION_SIZE: usize = 10_000;

/// Retrieves a 16-bit unsigned integer from `source` at `offset` using `endian`.
pub fn get_u16(source: &[u8], offset: usize, endian: Endian) -> Option<u16> {
    let end = offset.checked_add(2)?;
    let bytes = source.get(offset..end)?;
    let arr: [u8; 2] = bytes.try_into().ok()?;
    Some(match endian {
        Endian::Little => u16::from_le_bytes(arr),
        Endian::Big => u16::from_be_bytes(arr),
    })
}

/// Retrieves a 32-bit unsigned integer from `source` at `offset` using `endian`.
pub fn get_u32(source: &[u8], offset: usize, endian: Endian) -> Option<u32> {
    let end = offset.checked_add(4)?;
    let bytes = source.get(offset..end)?;
    let arr: [u8; 4] = bytes.try_into().ok()?;
    Some(match endian {
        Endian::Little => u32::from_le_bytes(arr),
        Endian::Big => u32::from_be_bytes(arr),
    })
}

/// Checks if `needle` matches `source` starting at `offset`.
pub fn is_signature_matched(source: &[u8], offset: usize, needle: &[u8]) -> bool {
    if needle.is_empty() {
        return true;
    }
    let end = match offset.checked_add(needle.len()) {
        Some(e) => e,
        None => return false,
    };
    if source.len() < end {
        return false;
    }
    &source[offset..end] == needle
}

/// Checks if `needle` is found in `source` within the window `[search_offset, search_offset +
/// search_range]`.
///
/// Matches C++ `IsSignatureFound` loop bound:
///   `i < search_offset + search_range - signature.size()`.
pub fn is_signature_found(
    source: &[u8],
    search_offset: usize,
    search_range: usize,
    needle: &[u8],
) -> bool {
    if needle.is_empty() {
        return true;
    }
    let end_bound = match search_offset.checked_add(search_range) {
        Some(end) => end,
        None => return false,
    };
    if source.len() < end_bound || needle.len() >= search_range {
        return false;
    }
    let max_start = search_offset + search_range - needle.len();
    (search_offset..max_start).any(|i| &source[i..i + needle.len()] == needle)
}

/// Trait for format identification checkers.
pub trait FormatChecker: Send + Sync {
    /// Returns the `RawImageType` associated with this checker.
    fn format_type(&self) -> RawImageType;

    /// Returns the requested data size (in bytes) needed for identification.
    fn requested_size(&self) -> usize;

    /// Returns `true` if `source` matches this RAW image format.
    fn is_my_type(&self, source: &[u8]) -> bool;
}

fn limit_source(source: &[u8], requested_size: usize) -> Option<&[u8]> {
    if source.len() < requested_size {
        None
    } else {
        Some(&source[..requested_size])
    }
}

/// Sony ARW format checker.
#[derive(Debug, Clone, Copy, Default)]
pub struct ArwChecker;

impl FormatChecker for ArwChecker {
    fn format_type(&self) -> RawImageType {
        RawImageType::Arw
    }

    fn requested_size(&self) -> usize {
        10_000
    }

    fn is_my_type(&self, source: &[u8]) -> bool {
        let limited = match limit_source(source, self.requested_size()) {
            Some(l) => l,
            None => return false,
        };
        let endian = match Endian::from_tiff_prefix(limited) {
            Some(e) => e,
            None => return false,
        };

        if get_u16(limited, 2, endian) != Some(0x2A) || get_u32(limited, 4, endian) != Some(8) {
            return false;
        }

        if !is_signature_found(limited, 0, limited.len(), b"SONY") {
            return false;
        }

        let section = b"\x00\xb0\x01\x00\x04\x00\x00\x00";
        let versions: [&[u8]; 7] = [
            b"\x02\x00",
            b"\x03\x00",
            b"\x03\x01",
            b"\x03\x02",
            b"\x03\x03",
            b"\x04\x00",
            b"\x05\x00",
        ];

        versions.iter().any(|ver| {
            let mut full_sig = [0u8; 10];
            full_sig[..8].copy_from_slice(section);
            full_sig[8..10].copy_from_slice(ver);
            is_signature_found(limited, 0, limited.len(), &full_sig)
        })
    }
}

/// Canon CR3 format checker.
#[derive(Debug, Clone, Copy, Default)]
pub struct Cr3Checker;

impl FormatChecker for Cr3Checker {
    fn format_type(&self) -> RawImageType {
        RawImageType::Cr3
    }

    fn requested_size(&self) -> usize {
        12
    }

    fn is_my_type(&self, source: &[u8]) -> bool {
        let limited = match limit_source(source, self.requested_size()) {
            Some(l) => l,
            None => return false,
        };
        is_signature_matched(limited, 4, b"ftypcrx ")
    }
}

/// Canon CR2 format checker.
#[derive(Debug, Clone, Copy, Default)]
pub struct Cr2Checker;

impl FormatChecker for Cr2Checker {
    fn format_type(&self) -> RawImageType {
        RawImageType::Cr2
    }

    fn requested_size(&self) -> usize {
        16
    }

    fn is_my_type(&self, source: &[u8]) -> bool {
        let limited = match limit_source(source, self.requested_size()) {
            Some(l) => l,
            None => return false,
        };
        let endian = match Endian::from_tiff_prefix(limited) {
            Some(e) => e,
            None => return false,
        };

        if get_u16(limited, 2, endian) != Some(42) {
            return false;
        }

        is_signature_matched(limited, 8, b"CR\x02\x00")
    }
}

/// Canon CRW format checker.
#[derive(Debug, Clone, Copy, Default)]
pub struct CrwChecker;

impl FormatChecker for CrwChecker {
    fn format_type(&self) -> RawImageType {
        RawImageType::Crw
    }

    fn requested_size(&self) -> usize {
        14
    }

    fn is_my_type(&self, source: &[u8]) -> bool {
        let limited = match limit_source(source, self.requested_size()) {
            Some(l) => l,
            None => return false,
        };
        let endian = match Endian::from_tiff_prefix(limited) {
            Some(e) => e,
            None => return false,
        };

        let sig = match endian {
            Endian::Big => b"\x00\x10\xba\xb0\xac\xbb\x00\x02".as_slice(),
            Endian::Little => b"HEAPCCDR".as_slice(),
        };

        is_signature_matched(limited, 6, sig)
    }
}

/// Kodak DCR format checker.
#[derive(Debug, Clone, Copy, Default)]
pub struct DcrChecker;

impl FormatChecker for DcrChecker {
    fn format_type(&self) -> RawImageType {
        RawImageType::Dcr
    }

    fn requested_size(&self) -> usize {
        5000
    }

    fn is_my_type(&self, source: &[u8]) -> bool {
        let limited = match limit_source(source, self.requested_size()) {
            Some(l) => l,
            None => return false,
        };
        let endian = match Endian::from_tiff_prefix(limited) {
            Some(e) => e,
            None => return false,
        };

        let sig = b"KODAK           ";
        if is_signature_matched(limited, 16, sig) {
            return true;
        }

        let (tag0, tag1) = match endian {
            Endian::Big => (b"\x03\xe9\x00\x02".as_slice(), b"\x0c\xe5\x00\x02".as_slice()),
            Endian::Little => (b"\xe9\x03\x02\x00".as_slice(), b"\xe5\x0c\x02\x00".as_slice()),
        };

        is_signature_found(limited, 0, limited.len(), tag0)
            && is_signature_found(limited, 0, limited.len(), tag1)
    }
}

/// Digital Negative (DNG) format checker.
#[derive(Debug, Clone, Copy, Default)]
pub struct DngChecker;

impl FormatChecker for DngChecker {
    fn format_type(&self) -> RawImageType {
        RawImageType::Dng
    }

    fn requested_size(&self) -> usize {
        1024
    }

    fn is_my_type(&self, source: &[u8]) -> bool {
        let limited = match limit_source(source, self.requested_size()) {
            Some(l) => l,
            None => return false,
        };
        let endian = match Endian::from_tiff_prefix(limited) {
            Some(e) => e,
            None => return false,
        };

        let dng_tags: [&[u8]; 5] = match endian {
            Endian::Big => [
                b"\xc6\x12\x00\x01\x00\x00\x00\x04",
                b"\xc6\x13\x00\x01\x00\x00\x00\x04",
                b"\xc6\x14\x00\x02",
                b"\xc6\x20",
                b"\xc6\x2d\x00\x04\x00\x00\x00\x01",
            ],
            Endian::Little => [
                b"\x12\xc6\x01\x00\x04\x00\x00\x00",
                b"\x13\xc6\x01\x00\x04\x00\x00\x00",
                b"\x14\xc6\x02\x00",
                b"\x20\xc6",
                b"\x2d\xc6\x04\x00\x01\x00\x00\x00",
            ],
        };

        let tags_found = dng_tags
            .iter()
            .filter(|tag| is_signature_found(limited, 0, limited.len(), tag))
            .count();

        tags_found >= 2
    }
}

/// Kodak KDC format checker.
#[derive(Debug, Clone, Copy, Default)]
pub struct KdcChecker;

impl FormatChecker for KdcChecker {
    fn format_type(&self) -> RawImageType {
        RawImageType::Kdc
    }

    fn requested_size(&self) -> usize {
        5000
    }

    fn is_my_type(&self, source: &[u8]) -> bool {
        let limited = match limit_source(source, self.requested_size()) {
            Some(l) => l,
            None => return false,
        };
        let endian = match Endian::from_tiff_prefix(limited) {
            Some(e) => e,
            None => return false,
        };

        let (tag0, tag1) = match endian {
            Endian::Big => (b"\xfa\x0d\x00\x01".as_slice(), b"\xfa\x00\x00\x02".as_slice()),
            Endian::Little => (b"\x0d\xfa\x01\x00".as_slice(), b"\x00\xfa\x02\x00".as_slice()),
        };

        is_signature_found(limited, 0, limited.len(), tag0)
            && is_signature_found(limited, 0, limited.len(), tag1)
    }
}

/// Leaf MOS format checker.
#[derive(Debug, Clone, Copy, Default)]
pub struct MosChecker;

impl FormatChecker for MosChecker {
    fn format_type(&self) -> RawImageType {
        RawImageType::Mos
    }

    fn requested_size(&self) -> usize {
        5000
    }

    fn is_my_type(&self, source: &[u8]) -> bool {
        let limited = match limit_source(source, self.requested_size()) {
            Some(l) => l,
            None => return false,
        };
        if Endian::from_tiff_prefix(limited).is_none() {
            return false;
        }

        let sig = b"PKTS\x00\x00\x00\x01";
        is_signature_found(limited, 0, limited.len(), sig)
    }
}

/// Minolta MRW format checker.
#[derive(Debug, Clone, Copy, Default)]
pub struct MrwChecker;

impl FormatChecker for MrwChecker {
    fn format_type(&self) -> RawImageType {
        RawImageType::Mrw
    }

    fn requested_size(&self) -> usize {
        4
    }

    fn is_my_type(&self, source: &[u8]) -> bool {
        let limited = match limit_source(source, self.requested_size()) {
            Some(l) => l,
            None => return false,
        };
        is_signature_matched(limited, 0, b"\x00MRM")
    }
}

fn contains_nrw_signature(source: &[u8]) -> bool {
    is_signature_found(source, 0, source.len(), b"NRW   ")
}

fn contains_nikon_signatures(source: &[u8], endian: Endian) -> bool {
    let sig_nikon = b"NIKON";
    let ref_black_white = match endian {
        Endian::Big => b"\x02\x14\x00\x05".as_slice(),
        Endian::Little => b"\x14\x02\x05\x00".as_slice(),
    };
    is_signature_found(source, 0, source.len(), sig_nikon)
        && is_signature_found(source, 0, source.len(), ref_black_white)
}

/// Nikon NEF format checker.
#[derive(Debug, Clone, Copy, Default)]
pub struct NefChecker;

impl FormatChecker for NefChecker {
    fn format_type(&self) -> RawImageType {
        RawImageType::Nef
    }

    fn requested_size(&self) -> usize {
        4000
    }

    fn is_my_type(&self, source: &[u8]) -> bool {
        let limited = match limit_source(source, self.requested_size()) {
            Some(l) => l,
            None => return false,
        };
        let endian = match Endian::from_tiff_prefix(limited) {
            Some(e) => e,
            None => return false,
        };

        if get_u16(limited, 2, endian) != Some(0x2A) {
            return false;
        }

        contains_nikon_signatures(limited, endian) && !contains_nrw_signature(limited)
    }
}

/// Nikon NRW format checker.
#[derive(Debug, Clone, Copy, Default)]
pub struct NrwChecker;

impl FormatChecker for NrwChecker {
    fn format_type(&self) -> RawImageType {
        RawImageType::Nrw
    }

    fn requested_size(&self) -> usize {
        4000
    }

    fn is_my_type(&self, source: &[u8]) -> bool {
        let limited = match limit_source(source, self.requested_size()) {
            Some(l) => l,
            None => return false,
        };
        let endian = match Endian::from_tiff_prefix(limited) {
            Some(e) => e,
            None => return false,
        };

        if get_u16(limited, 2, endian) != Some(0x2A) || get_u32(limited, 4, endian) != Some(8) {
            return false;
        }

        contains_nikon_signatures(limited, endian) && contains_nrw_signature(limited)
    }
}

/// Epson ERF format checker.
#[derive(Debug, Clone, Copy, Default)]
pub struct ErfChecker;

impl FormatChecker for ErfChecker {
    fn format_type(&self) -> RawImageType {
        RawImageType::Erf
    }

    fn requested_size(&self) -> usize {
        1024
    }

    fn is_my_type(&self, source: &[u8]) -> bool {
        let limited = match limit_source(source, self.requested_size()) {
            Some(l) => l,
            None => return false,
        };
        let endian = match Endian::from_tiff_prefix(limited) {
            Some(e) => e,
            None => return false,
        };

        if get_u16(limited, 2, endian) != Some(0x2A) {
            return false;
        }

        is_signature_found(limited, 0, limited.len(), b"EPSON")
    }
}

/// Olympus ORF format checker.
#[derive(Debug, Clone, Copy, Default)]
pub struct OrfChecker;

impl FormatChecker for OrfChecker {
    fn format_type(&self) -> RawImageType {
        RawImageType::Orf
    }

    fn requested_size(&self) -> usize {
        3000
    }

    fn is_my_type(&self, source: &[u8]) -> bool {
        let limited = match limit_source(source, self.requested_size()) {
            Some(l) => l,
            None => return false,
        };
        let endian = match Endian::from_tiff_prefix(limited) {
            Some(e) => e,
            None => return false,
        };

        let tag = match get_u16(limited, 2, endian) {
            Some(t) => t,
            None => return false,
        };

        if tag != 0x4F52 && tag != 0x5352 {
            return false;
        }

        is_signature_found(limited, 0, limited.len(), b"OLYMP")
            || is_signature_found(limited, 0, limited.len(), b"OM Digital")
    }
}

/// Pentax PEF format checker.
#[derive(Debug, Clone, Copy, Default)]
pub struct PefChecker;

impl FormatChecker for PefChecker {
    fn format_type(&self) -> RawImageType {
        RawImageType::Pef
    }

    fn requested_size(&self) -> usize {
        1280
    }

    fn is_my_type(&self, source: &[u8]) -> bool {
        let limited = match limit_source(source, self.requested_size()) {
            Some(l) => l,
            None => return false,
        };
        let endian = match Endian::from_tiff_prefix(limited) {
            Some(e) => e,
            None => return false,
        };

        if get_u16(limited, 2, endian) != Some(0x2A) || get_u32(limited, 4, endian) != Some(8) {
            return false;
        }

        let sig_aoc = b"\x41\x4f\x43\x00\x4d\x4d";
        let sig_pentax = b"\x50\x45\x4e\x54\x41\x58\x20\x00";

        is_signature_found(limited, 0, limited.len(), sig_aoc)
            || is_signature_found(limited, 0, limited.len(), sig_pentax)
    }
}

/// Phase One RAW format checker.
#[derive(Debug, Clone, Copy, Default)]
pub struct PhaseOneChecker;

impl FormatChecker for PhaseOneChecker {
    fn format_type(&self) -> RawImageType {
        RawImageType::PhaseOne
    }

    fn requested_size(&self) -> usize {
        16
    }

    fn is_my_type(&self, source: &[u8]) -> bool {
        let limited = match limit_source(source, self.requested_size()) {
            Some(l) => l,
            None => return false,
        };
        is_signature_matched(limited, 0, b"MMMMRawT")
            || is_signature_matched(limited, 8, b"MMMMRawT")
    }
}

/// Logitech Fotoman PXN format checker.
#[derive(Debug, Clone, Copy, Default)]
pub struct PxnChecker;

impl FormatChecker for PxnChecker {
    fn format_type(&self) -> RawImageType {
        RawImageType::Pxn
    }

    fn requested_size(&self) -> usize {
        4
    }

    fn is_my_type(&self, source: &[u8]) -> bool {
        let limited = match limit_source(source, self.requested_size()) {
            Some(l) => l,
            None => return false,
        };
        is_signature_matched(limited, 0, b"PXN\x00")
    }
}

/// Apple QTK format checker.
#[derive(Debug, Clone, Copy, Default)]
pub struct QtkChecker;

impl FormatChecker for QtkChecker {
    fn format_type(&self) -> RawImageType {
        RawImageType::Qtk
    }

    fn requested_size(&self) -> usize {
        8
    }

    fn is_my_type(&self, source: &[u8]) -> bool {
        let limited = match limit_source(source, self.requested_size()) {
            Some(l) => l,
            None => return false,
        };
        is_signature_matched(limited, 0, b"qktk\x00\x00\x00\x08")
            || is_signature_matched(limited, 0, b"qktn\x00\x00\x00\x08")
    }
}

/// Fuji RAF format checker.
#[derive(Debug, Clone, Copy, Default)]
pub struct RafChecker;

impl FormatChecker for RafChecker {
    fn format_type(&self) -> RawImageType {
        RawImageType::Raf
    }

    fn requested_size(&self) -> usize {
        8
    }

    fn is_my_type(&self, source: &[u8]) -> bool {
        let limited = match limit_source(source, self.requested_size()) {
            Some(l) => l,
            None => return false,
        };
        is_signature_matched(limited, 0, b"FUJIFILM")
    }
}

/// Hasselblad 3FR format checker.
#[derive(Debug, Clone, Copy, Default)]
pub struct Raw3frChecker;

impl FormatChecker for Raw3frChecker {
    fn format_type(&self) -> RawImageType {
        RawImageType::Raw3fr
    }

    fn requested_size(&self) -> usize {
        1024
    }

    fn is_my_type(&self, source: &[u8]) -> bool {
        let limited = match limit_source(source, self.requested_size()) {
            Some(l) => l,
            None => return false,
        };
        let endian = match Endian::from_tiff_prefix(limited) {
            Some(e) => e,
            None => return false,
        };

        if get_u16(limited, 2, endian) != Some(0x2A) {
            return false;
        }

        is_signature_found(limited, 0, limited.len(), b"Hasselblad")
    }
}

/// Contax N RAW format checker.
#[derive(Debug, Clone, Copy, Default)]
pub struct RawContaxNChecker;

impl FormatChecker for RawContaxNChecker {
    fn format_type(&self) -> RawImageType {
        RawImageType::RawContaxN
    }

    fn requested_size(&self) -> usize {
        36
    }

    fn is_my_type(&self, source: &[u8]) -> bool {
        let limited = match limit_source(source, self.requested_size()) {
            Some(l) => l,
            None => return false,
        };
        is_signature_matched(limited, 25, b"ARECOYK")
    }
}

/// Rollei RDC RAW format checker.
#[derive(Debug, Clone, Copy, Default)]
pub struct RdcChecker;

impl FormatChecker for RdcChecker {
    fn format_type(&self) -> RawImageType {
        RawImageType::Rdc
    }

    fn requested_size(&self) -> usize {
        16
    }

    fn is_my_type(&self, source: &[u8]) -> bool {
        let limited = match limit_source(source, self.requested_size()) {
            Some(l) => l,
            None => return false,
        };
        is_signature_matched(limited, 0, b"DSC-Image <C> by")
    }
}

/// Panasonic RW2 format checker.
#[derive(Debug, Clone, Copy, Default)]
pub struct Rw2Checker;

impl FormatChecker for Rw2Checker {
    fn format_type(&self) -> RawImageType {
        RawImageType::Rw2
    }

    fn requested_size(&self) -> usize {
        4
    }

    fn is_my_type(&self, source: &[u8]) -> bool {
        let limited = match limit_source(source, self.requested_size()) {
            Some(l) => l,
            None => return false,
        };
        let endian = match Endian::from_tiff_prefix(limited) {
            Some(e) => e,
            None => return false,
        };

        get_u16(limited, 2, endian) == Some(0x55)
    }
}

/// Samsung SRW format checker.
#[derive(Debug, Clone, Copy, Default)]
pub struct SrwChecker;

impl FormatChecker for SrwChecker {
    fn format_type(&self) -> RawImageType {
        RawImageType::Srw
    }

    fn requested_size(&self) -> usize {
        256
    }

    fn is_my_type(&self, source: &[u8]) -> bool {
        let limited = match limit_source(source, self.requested_size()) {
            Some(l) => l,
            None => return false,
        };
        let endian = match Endian::from_tiff_prefix(limited) {
            Some(e) => e,
            None => return false,
        };

        if get_u16(limited, 2, endian) != Some(0x2A) || get_u32(limited, 4, endian) != Some(8) {
            return false;
        }

        is_signature_found(limited, 0, limited.len(), b"SAMSUNG")
    }
}

/// Sigma / Polaroid X3F format checker.
#[derive(Debug, Clone, Copy, Default)]
pub struct X3fChecker;

impl FormatChecker for X3fChecker {
    fn format_type(&self) -> RawImageType {
        RawImageType::X3f
    }

    fn requested_size(&self) -> usize {
        4
    }

    fn is_my_type(&self, source: &[u8]) -> bool {
        let limited = match limit_source(source, self.requested_size()) {
            Some(l) => l,
            None => return false,
        };
        is_signature_matched(limited, 0, b"FOVb")
    }
}

/// Static registry of all RAW format checkers in exact C++ `TypeCheckerList` stable_sort order.
pub static ALL_CHECKERS: &[&'static (dyn FormatChecker + 'static)] = &[
    &MrwChecker,
    &PxnChecker,
    &Rw2Checker,
    &X3fChecker,
    &QtkChecker,
    &RafChecker,
    &Cr3Checker,
    &CrwChecker,
    &Cr2Checker,
    &PhaseOneChecker,
    &RdcChecker,
    &RawContaxNChecker,
    &SrwChecker,
    &ErfChecker,
    &Raw3frChecker,
    &DngChecker,
    &PefChecker,
    &OrfChecker,
    &NefChecker,
    &NrwChecker,
    &DcrChecker,
    &KdcChecker,
    &MosChecker,
    &ArwChecker,
];

/// Checks if `image_type` is a RAW image format (i.e. not `RawImageType::NonRaw`).
pub fn is_raw(image_type: RawImageType) -> bool {
    image_type != RawImageType::NonRaw
}

/// Identifies the RAW image format of `source`.
///
/// Returns `RawImageType::NonRaw` if no format matches.
pub fn recognize_raw_image_type_lite(source: &[u8]) -> RawImageType {
    ALL_CHECKERS
        .iter()
        .find(|checker| checker.is_my_type(source))
        .map(|checker| checker.format_type())
        .unwrap_or(RawImageType::NonRaw)
}

/// Returns `true` if `source` matches `target_type`.
pub fn is_of_type(source: &[u8], target_type: RawImageType) -> bool {
    if target_type == RawImageType::NonRaw {
        return false;
    }
    ALL_CHECKERS
        .iter()
        .find(|checker| checker.format_type() == target_type)
        .map(|checker| checker.is_my_type(source))
        .unwrap_or(false)
}

/// Returns the maximum number of header bytes required to identify any RAW format.
pub fn get_number_of_bytes_for_is_raw_lite() -> usize {
    ALL_CHECKERS.iter().map(|checker| checker.requested_size()).max().unwrap_or(0)
}

/// Returns the smallest checker window size strictly greater than `after`,
/// or `0` if no checker requires more than `after` bytes.
pub fn next_recognition_window_size(after: usize) -> usize {
    ALL_CHECKERS
        .iter()
        .map(|checker| checker.requested_size())
        .filter(|&size| size > after)
        .min()
        .unwrap_or(0)
}

/// Identifies the RAW image format of `source` using only checkers whose
/// `requested_size()` equals `window_size`, preserving `ALL_CHECKERS` order.
pub fn recognize_raw_image_type_for_window(source: &[u8], window_size: usize) -> RawImageType {
    ALL_CHECKERS
        .iter()
        .filter(|checker| checker.requested_size() == window_size)
        .find(|checker| checker.is_my_type(source))
        .map(|checker| checker.format_type())
        .unwrap_or(RawImageType::NonRaw)
}

/// Returns the number of header bytes required to identify `target_type`.
pub fn get_number_of_bytes_for_is_of_type(target_type: RawImageType) -> usize {
    ALL_CHECKERS
        .iter()
        .find(|checker| checker.format_type() == target_type)
        .map(|checker| checker.requested_size())
        .unwrap_or(0)
}

/// Returns `true` if `source` is recognized as any known RAW image type.
pub fn is_raw_lite(source: &[u8]) -> bool {
    is_raw(recognize_raw_image_type_lite(source))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_raw() {
        assert!(!is_raw(RawImageType::NonRaw));
        assert!(is_raw(RawImageType::Arw));
        assert!(is_raw(RawImageType::Cr2));
        assert!(is_raw(RawImageType::Cr3));
        assert!(is_raw(RawImageType::Crw));
        assert!(is_raw(RawImageType::Dcr));
        assert!(is_raw(RawImageType::Dng));
        assert!(is_raw(RawImageType::Kdc));
        assert!(is_raw(RawImageType::Mos));
        assert!(is_raw(RawImageType::Mrw));
        assert!(is_raw(RawImageType::Nef));
        assert!(is_raw(RawImageType::Nrw));
        assert!(is_raw(RawImageType::Orf));
        assert!(is_raw(RawImageType::Pef));
        assert!(is_raw(RawImageType::Qtk));
        assert!(is_raw(RawImageType::Raf));
        assert!(is_raw(RawImageType::Raw3fr));
        assert!(is_raw(RawImageType::RawContaxN));
        assert!(is_raw(RawImageType::Rw2));
        assert!(is_raw(RawImageType::Srw));
        assert!(is_raw(RawImageType::X3f));
    }

    #[test]
    fn test_recognize_mrw() {
        let data = b"\x00MRM_padding";
        assert_eq!(recognize_raw_image_type_lite(data), RawImageType::Mrw);
        assert!(is_of_type(data, RawImageType::Mrw));
        assert!(is_raw_lite(data));
    }

    #[test]
    fn test_recognize_rw2() {
        let mut data = vec![0x49, 0x49, 0x55, 0x00];
        assert_eq!(recognize_raw_image_type_lite(&data), RawImageType::Rw2);
        assert!(is_of_type(&data, RawImageType::Rw2));

        data = vec![0x4D, 0x4D, 0x00, 0x55];
        assert_eq!(recognize_raw_image_type_lite(&data), RawImageType::Rw2);
    }

    #[test]
    fn test_recognize_x3f() {
        let data = b"FOVb_extra_data";
        assert_eq!(recognize_raw_image_type_lite(data), RawImageType::X3f);
        assert!(is_of_type(data, RawImageType::X3f));
    }

    #[test]
    fn test_recognize_cr3() {
        let data = b"\x00\x00\x00\x18ftypcrx \x00\x00\x00\x00";
        assert_eq!(recognize_raw_image_type_lite(data), RawImageType::Cr3);
        assert!(is_of_type(data, RawImageType::Cr3));
    }

    #[test]
    fn test_recognize_qtk() {
        let data1 = b"qktk\x00\x00\x00\x08extra";
        let data2 = b"qktn\x00\x00\x00\x08extra";
        assert_eq!(recognize_raw_image_type_lite(data1), RawImageType::Qtk);
        assert_eq!(recognize_raw_image_type_lite(data2), RawImageType::Qtk);
        assert!(is_of_type(data1, RawImageType::Qtk));
    }

    #[test]
    fn test_recognize_raf() {
        let data = b"FUJIFILM_RAW_HEADER";
        assert_eq!(recognize_raw_image_type_lite(data), RawImageType::Raf);
        assert!(is_of_type(data, RawImageType::Raf));
    }

    #[test]
    fn test_recognize_cr2() {
        // II + 42 + offset 8 + "CR\x02\x00"
        let mut data = vec![0u8; 16];
        data[0..2].copy_from_slice(b"II");
        data[2..4].copy_from_slice(&42u16.to_le_bytes());
        data[8..12].copy_from_slice(b"CR\x02\x00");
        assert_eq!(recognize_raw_image_type_lite(&data), RawImageType::Cr2);
        assert!(is_of_type(&data, RawImageType::Cr2));
    }

    #[test]
    fn test_recognize_crw() {
        let mut data = vec![0u8; 14];
        data[0..2].copy_from_slice(b"II");
        data[6..14].copy_from_slice(b"HEAPCCDR");
        assert_eq!(recognize_raw_image_type_lite(&data), RawImageType::Crw);
        assert!(is_of_type(&data, RawImageType::Crw));
    }

    #[test]
    fn test_recognize_raw_contax_n() {
        let mut data = vec![0u8; 36];
        data[25..32].copy_from_slice(b"ARECOYK");
        assert_eq!(recognize_raw_image_type_lite(&data), RawImageType::RawContaxN);
        assert!(is_of_type(&data, RawImageType::RawContaxN));
    }

    #[test]
    fn test_recognize_srw() {
        let mut buf = vec![0u8; 256];
        buf[0..2].copy_from_slice(b"MM");
        buf[2..4].copy_from_slice(&0x2Au16.to_be_bytes());
        buf[4..8].copy_from_slice(&8u32.to_be_bytes());
        buf[50..57].copy_from_slice(b"SAMSUNG");
        assert_eq!(recognize_raw_image_type_lite(&buf), RawImageType::Srw);
        assert!(is_of_type(&buf, RawImageType::Srw));
    }

    #[test]
    fn test_recognize_raw_3fr() {
        let mut buf = vec![0u8; 1024];
        buf[0..2].copy_from_slice(b"II");
        buf[2..4].copy_from_slice(&0x2Au16.to_le_bytes());
        buf[100..110].copy_from_slice(b"Hasselblad");
        assert_eq!(recognize_raw_image_type_lite(&buf), RawImageType::Raw3fr);
        assert!(is_of_type(&buf, RawImageType::Raw3fr));
    }

    #[test]
    fn test_recognize_pef() {
        let mut buf = vec![0u8; 1280];
        buf[0..2].copy_from_slice(b"MM");
        buf[2..4].copy_from_slice(&0x2Au16.to_be_bytes());
        buf[4..8].copy_from_slice(&8u32.to_be_bytes());
        buf[20..28].copy_from_slice(b"PENTAX \x00");
        assert_eq!(recognize_raw_image_type_lite(&buf), RawImageType::Pef);
        assert!(is_of_type(&buf, RawImageType::Pef));
    }

    #[test]
    fn test_recognize_orf() {
        let mut buf = vec![0u8; 3000];
        buf[0..2].copy_from_slice(b"II");
        buf[2..4].copy_from_slice(&0x4F52u16.to_le_bytes());
        buf[50..55].copy_from_slice(b"OLYMP");
        assert_eq!(recognize_raw_image_type_lite(&buf), RawImageType::Orf);
        assert!(is_of_type(&buf, RawImageType::Orf));
    }

    #[test]
    fn test_recognize_nef_and_nrw() {
        // NEF: II + 0x2A + NIKON + RefBlackWhite + NO NRW
        let mut nef_buf = vec![0u8; 4000];
        nef_buf[0..2].copy_from_slice(b"II");
        nef_buf[2..4].copy_from_slice(&0x2Au16.to_le_bytes());
        nef_buf[10..15].copy_from_slice(b"NIKON");
        nef_buf[20..24].copy_from_slice(b"\x14\x02\x05\x00");
        assert_eq!(recognize_raw_image_type_lite(&nef_buf), RawImageType::Nef);

        // NRW: II + 0x2A + offset 8 + NIKON + RefBlackWhite + NRW
        let mut nrw_buf = nef_buf.clone();
        nrw_buf[4..8].copy_from_slice(&8u32.to_le_bytes());
        nrw_buf[30..36].copy_from_slice(b"NRW   ");
        assert_eq!(recognize_raw_image_type_lite(&nrw_buf), RawImageType::Nrw);
    }

    #[test]
    fn test_recognize_dcr() {
        let mut buf = vec![0u8; 5000];
        buf[0..2].copy_from_slice(b"MM");
        buf[16..32].copy_from_slice(b"KODAK           ");
        assert_eq!(recognize_raw_image_type_lite(&buf), RawImageType::Dcr);

        // Tag search fallback
        let mut buf2 = vec![0u8; 5000];
        buf2[0..2].copy_from_slice(b"MM");
        buf2[10..14].copy_from_slice(b"\x03\xe9\x00\x02");
        buf2[20..24].copy_from_slice(b"\x0c\xe5\x00\x02");
        assert_eq!(recognize_raw_image_type_lite(&buf2), RawImageType::Dcr);
    }

    #[test]
    fn test_recognize_kdc() {
        let mut buf = vec![0u8; 5000];
        buf[0..2].copy_from_slice(b"II");
        buf[10..14].copy_from_slice(b"\x0d\xfa\x01\x00");
        buf[20..24].copy_from_slice(b"\x00\xfa\x02\x00");
        assert_eq!(recognize_raw_image_type_lite(&buf), RawImageType::Kdc);
    }

    #[test]
    fn test_recognize_mos() {
        let mut buf = vec![0u8; 5000];
        buf[0..2].copy_from_slice(b"II");
        buf[10..18].copy_from_slice(b"PKTS\x00\x00\x00\x01");
        assert_eq!(recognize_raw_image_type_lite(&buf), RawImageType::Mos);
    }

    #[test]
    fn test_recognize_arw() {
        let mut buf = vec![0u8; 10_000];
        buf[0..2].copy_from_slice(b"II");
        buf[2..4].copy_from_slice(&0x2Au16.to_le_bytes());
        buf[4..8].copy_from_slice(&8u32.to_le_bytes());
        buf[100..104].copy_from_slice(b"SONY");
        let sig = b"\x00\xb0\x01\x00\x04\x00\x00\x00\x02\x00";
        buf[200..200 + sig.len()].copy_from_slice(sig);
        assert_eq!(recognize_raw_image_type_lite(&buf), RawImageType::Arw);
    }

    #[test]
    fn test_recognize_dng() {
        let mut buf = vec![0u8; 1024];
        buf[0..2].copy_from_slice(b"II");
        buf[10..18].copy_from_slice(b"\x12\xc6\x01\x00\x04\x00\x00\x00");
        buf[30..38].copy_from_slice(b"\x13\xc6\x01\x00\x04\x00\x00\x00");
        assert_eq!(recognize_raw_image_type_lite(&buf), RawImageType::Dng);
    }

    #[test]
    fn test_dng_precedence() {
        // Raw3fr and DNG tags present: Raw3fr should match first
        let mut buf = vec![0u8; 1024];
        buf[0..2].copy_from_slice(b"II");
        buf[2..4].copy_from_slice(&0x2Au16.to_le_bytes());
        buf[10..18].copy_from_slice(b"\x12\xc6\x01\x00\x04\x00\x00\x00");
        buf[30..38].copy_from_slice(b"\x13\xc6\x01\x00\x04\x00\x00\x00");
        buf[100..110].copy_from_slice(b"Hasselblad");

        assert_eq!(recognize_raw_image_type_lite(&buf), RawImageType::Raw3fr);
    }

    #[test]
    fn test_truncated_headers() {
        assert_eq!(recognize_raw_image_type_lite(&[]), RawImageType::NonRaw);
        assert_eq!(recognize_raw_image_type_lite(b"II"), RawImageType::NonRaw);
        assert_eq!(recognize_raw_image_type_lite(b"\x00MR"), RawImageType::NonRaw);
        assert!(!is_raw_lite(&[]));
    }

    #[test]
    fn test_non_raw_fallback() {
        let jpeg_header = b"\xFF\xD8\xFF\xE0\x00\x10JFIF\x00\x01\x01\x01\x00\x48\x00\x48\x00\x00";
        assert_eq!(recognize_raw_image_type_lite(jpeg_header), RawImageType::NonRaw);
        assert!(!is_raw_lite(jpeg_header));
        assert!(!is_of_type(jpeg_header, RawImageType::Arw));
    }

    #[test]
    fn test_byte_counts() {
        assert_eq!(get_number_of_bytes_for_is_raw_lite(), 10_000);
        assert_eq!(get_number_of_bytes_for_is_of_type(RawImageType::Mrw), 4);
        assert_eq!(get_number_of_bytes_for_is_of_type(RawImageType::Cr3), 12);
        assert_eq!(get_number_of_bytes_for_is_of_type(RawImageType::Arw), 10_000);
        assert_eq!(get_number_of_bytes_for_is_of_type(RawImageType::NonRaw), 0);
    }
}
