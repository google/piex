use std::io::Cursor;

use piex::format_recognizer::{
    get_number_of_bytes_for_is_of_type, get_number_of_bytes_for_is_raw_lite, get_u16, get_u32,
    is_of_type, is_raw, is_raw_lite, is_signature_found, is_signature_matched,
    recognize_raw_image_type_lite, MAX_HEADER_INSPECTION_SIZE,
};
use piex::stream::Endian;
use piex::types::RawImageType;

// ============================================================================
// Synthetic Header Generators
// ============================================================================

fn make_arw_header(little_endian: bool) -> Vec<u8> {
    let mut buf = vec![0u8; 10_000];
    if little_endian {
        buf[0..2].copy_from_slice(b"II");
        buf[2..4].copy_from_slice(&0x2Au16.to_le_bytes());
        buf[4..8].copy_from_slice(&8u32.to_le_bytes());
    } else {
        buf[0..2].copy_from_slice(b"MM");
        buf[2..4].copy_from_slice(&0x2Au16.to_be_bytes());
        buf[4..8].copy_from_slice(&8u32.to_be_bytes());
    }
    buf[100..104].copy_from_slice(b"SONY");
    let sig = b"\x00\xb0\x01\x00\x04\x00\x00\x00\x02\x00";
    buf[200..200 + sig.len()].copy_from_slice(sig);
    buf
}

fn make_cr2_header(little_endian: bool) -> Vec<u8> {
    let mut buf = vec![0u8; 16];
    if little_endian {
        buf[0..2].copy_from_slice(b"II");
        buf[2..4].copy_from_slice(&42u16.to_le_bytes());
    } else {
        buf[0..2].copy_from_slice(b"MM");
        buf[2..4].copy_from_slice(&42u16.to_be_bytes());
    }
    buf[8..12].copy_from_slice(b"CR\x02\x00");
    buf
}

fn make_cr3_header() -> Vec<u8> {
    let mut buf = vec![0u8; 16];
    buf[0..4].copy_from_slice(&16u32.to_be_bytes());
    buf[4..12].copy_from_slice(b"ftypcrx ");
    buf
}

fn make_crw_header(little_endian: bool) -> Vec<u8> {
    let mut buf = vec![0u8; 14];
    if little_endian {
        buf[0..2].copy_from_slice(b"II");
        buf[6..14].copy_from_slice(b"HEAPCCDR");
    } else {
        buf[0..2].copy_from_slice(b"MM");
        buf[6..14].copy_from_slice(b"\x00\x10\xba\xb0\xac\xbb\x00\x02");
    }
    buf
}

fn make_dcr_header(by_signature: bool, little_endian: bool) -> Vec<u8> {
    let mut buf = vec![0u8; 5000];
    if little_endian {
        buf[0..2].copy_from_slice(b"II");
    } else {
        buf[0..2].copy_from_slice(b"MM");
    }

    if by_signature {
        buf[16..32].copy_from_slice(b"KODAK           ");
    } else if little_endian {
        buf[10..14].copy_from_slice(b"\xe9\x03\x02\x00");
        buf[20..24].copy_from_slice(b"\xe5\x0c\x02\x00");
    } else {
        buf[10..14].copy_from_slice(b"\x03\xe9\x00\x02");
        buf[20..24].copy_from_slice(b"\x0c\xe5\x00\x02");
    }
    buf
}

fn make_dng_header(little_endian: bool) -> Vec<u8> {
    let mut buf = vec![0u8; 1024];
    if little_endian {
        buf[0..2].copy_from_slice(b"II");
        buf[10..18].copy_from_slice(b"\x12\xc6\x01\x00\x04\x00\x00\x00");
        buf[30..38].copy_from_slice(b"\x13\xc6\x01\x00\x04\x00\x00\x00");
    } else {
        buf[0..2].copy_from_slice(b"MM");
        buf[10..18].copy_from_slice(b"\xc6\x12\x00\x01\x00\x00\x00\x04");
        buf[30..38].copy_from_slice(b"\xc6\x13\x00\x01\x00\x00\x00\x04");
    }
    buf
}

fn make_erf_header(little_endian: bool) -> Vec<u8> {
    let mut buf = vec![0u8; 1024];
    if little_endian {
        buf[0..2].copy_from_slice(b"II");
        buf[2..4].copy_from_slice(&0x2Au16.to_le_bytes());
    } else {
        buf[0..2].copy_from_slice(b"MM");
        buf[2..4].copy_from_slice(&0x2Au16.to_be_bytes());
    }
    buf[100..105].copy_from_slice(b"EPSON");
    buf
}

fn make_kdc_header(little_endian: bool) -> Vec<u8> {
    let mut buf = vec![0u8; 5000];
    if little_endian {
        buf[0..2].copy_from_slice(b"II");
        buf[10..14].copy_from_slice(b"\x0d\xfa\x01\x00");
        buf[20..24].copy_from_slice(b"\x00\xfa\x02\x00");
    } else {
        buf[0..2].copy_from_slice(b"MM");
        buf[10..14].copy_from_slice(b"\xfa\x0d\x00\x01");
        buf[20..24].copy_from_slice(b"\xfa\x00\x00\x02");
    }
    buf
}

fn make_mos_header(little_endian: bool) -> Vec<u8> {
    let mut buf = vec![0u8; 5000];
    if little_endian {
        buf[0..2].copy_from_slice(b"II");
    } else {
        buf[0..2].copy_from_slice(b"MM");
    }
    buf[10..18].copy_from_slice(b"PKTS\x00\x00\x00\x01");
    buf
}

fn make_mrw_header() -> Vec<u8> {
    b"\x00MRM_synthetic_minolta_header".to_vec()
}

fn make_nef_header(little_endian: bool) -> Vec<u8> {
    let mut buf = vec![0u8; 4000];
    if little_endian {
        buf[0..2].copy_from_slice(b"II");
        buf[2..4].copy_from_slice(&0x2Au16.to_le_bytes());
        buf[10..15].copy_from_slice(b"NIKON");
        buf[20..24].copy_from_slice(b"\x14\x02\x05\x00");
    } else {
        buf[0..2].copy_from_slice(b"MM");
        buf[2..4].copy_from_slice(&0x2Au16.to_be_bytes());
        buf[10..15].copy_from_slice(b"NIKON");
        buf[20..24].copy_from_slice(b"\x02\x14\x00\x05");
    }
    buf
}

fn make_nrw_header(little_endian: bool) -> Vec<u8> {
    let mut buf = make_nef_header(little_endian);
    if little_endian {
        buf[4..8].copy_from_slice(&8u32.to_le_bytes());
    } else {
        buf[4..8].copy_from_slice(&8u32.to_be_bytes());
    }
    buf[30..36].copy_from_slice(b"NRW   ");
    buf
}

fn make_orf_header(little_endian: bool, is_or: bool) -> Vec<u8> {
    let mut buf = vec![0u8; 3000];
    let tag = if is_or { 0x4F52u16 } else { 0x5352u16 };
    if little_endian {
        buf[0..2].copy_from_slice(b"II");
        buf[2..4].copy_from_slice(&tag.to_le_bytes());
    } else {
        buf[0..2].copy_from_slice(b"MM");
        buf[2..4].copy_from_slice(&tag.to_be_bytes());
    }
    buf[50..55].copy_from_slice(b"OLYMP");
    buf
}

fn make_orf_om_digital_header(little_endian: bool, is_or: bool) -> Vec<u8> {
    let mut buf = vec![0u8; 3000];
    let tag = if is_or { 0x4F52u16 } else { 0x5352u16 };
    if little_endian {
        buf[0..2].copy_from_slice(b"II");
        buf[2..4].copy_from_slice(&tag.to_le_bytes());
    } else {
        buf[0..2].copy_from_slice(b"MM");
        buf[2..4].copy_from_slice(&tag.to_be_bytes());
    }
    buf[50..60].copy_from_slice(b"OM Digital");
    buf
}

fn make_pef_header(little_endian: bool, use_pentax_sig: bool) -> Vec<u8> {
    let mut buf = vec![0u8; 1280];
    if little_endian {
        buf[0..2].copy_from_slice(b"II");
        buf[2..4].copy_from_slice(&0x2Au16.to_le_bytes());
        buf[4..8].copy_from_slice(&8u32.to_le_bytes());
    } else {
        buf[0..2].copy_from_slice(b"MM");
        buf[2..4].copy_from_slice(&0x2Au16.to_be_bytes());
        buf[4..8].copy_from_slice(&8u32.to_be_bytes());
    }
    if use_pentax_sig {
        buf[20..28].copy_from_slice(b"\x50\x45\x4e\x54\x41\x58\x20\x00");
    } else {
        buf[20..26].copy_from_slice(b"\x41\x4f\x43\x00\x4d\x4d");
    }
    buf
}

fn make_phase_one_header(offset: usize) -> Vec<u8> {
    let mut buf = vec![0u8; 16];
    buf[offset..offset + 8].copy_from_slice(b"MMMMRawT");
    buf
}

fn make_pxn_header() -> Vec<u8> {
    b"PXN\x00".to_vec()
}

fn make_qtk_header(variant_n: bool) -> Vec<u8> {
    let mut buf = vec![0u8; 16];
    if variant_n {
        buf[0..8].copy_from_slice(b"qktn\x00\x00\x00\x08");
    } else {
        buf[0..8].copy_from_slice(b"qktk\x00\x00\x00\x08");
    }
    buf
}

fn make_raf_header() -> Vec<u8> {
    b"FUJIFILM_RAW_HEADER_DATA_STREAM".to_vec()
}

fn make_raw_3fr_header(little_endian: bool) -> Vec<u8> {
    let mut buf = vec![0u8; 1024];
    if little_endian {
        buf[0..2].copy_from_slice(b"II");
        buf[2..4].copy_from_slice(&0x2Au16.to_le_bytes());
    } else {
        buf[0..2].copy_from_slice(b"MM");
        buf[2..4].copy_from_slice(&0x2Au16.to_be_bytes());
    }
    buf[100..110].copy_from_slice(b"Hasselblad");
    buf
}

fn make_raw_contax_n_header() -> Vec<u8> {
    let mut buf = vec![0u8; 40];
    buf[25..32].copy_from_slice(b"ARECOYK");
    buf
}

fn make_rdc_header() -> Vec<u8> {
    b"DSC-Image <C> by".to_vec()
}

fn make_rw2_header(little_endian: bool) -> Vec<u8> {
    let mut buf = vec![0u8; 8];
    if little_endian {
        buf[0..2].copy_from_slice(b"II");
        buf[2..4].copy_from_slice(&0x55u16.to_le_bytes());
    } else {
        buf[0..2].copy_from_slice(b"MM");
        buf[2..4].copy_from_slice(&0x55u16.to_be_bytes());
    }
    buf
}

fn make_srw_header(little_endian: bool) -> Vec<u8> {
    let mut buf = vec![0u8; 256];
    if little_endian {
        buf[0..2].copy_from_slice(b"II");
        buf[2..4].copy_from_slice(&0x2Au16.to_le_bytes());
        buf[4..8].copy_from_slice(&8u32.to_le_bytes());
    } else {
        buf[0..2].copy_from_slice(b"MM");
        buf[2..4].copy_from_slice(&0x2Au16.to_be_bytes());
        buf[4..8].copy_from_slice(&8u32.to_be_bytes());
    }
    buf[50..57].copy_from_slice(b"SAMSUNG");
    buf
}

fn make_x3f_header() -> Vec<u8> {
    b"FOVb_sigma_foveon_header".to_vec()
}

// ============================================================================
// Tests
// ============================================================================

const ALL_24_TYPES: [RawImageType; 24] = [
    RawImageType::Arw,
    RawImageType::Cr2,
    RawImageType::Cr3,
    RawImageType::Crw,
    RawImageType::Dcr,
    RawImageType::Dng,
    RawImageType::Erf,
    RawImageType::Kdc,
    RawImageType::Mos,
    RawImageType::Mrw,
    RawImageType::Nef,
    RawImageType::Nrw,
    RawImageType::Orf,
    RawImageType::Pef,
    RawImageType::PhaseOne,
    RawImageType::Pxn,
    RawImageType::Qtk,
    RawImageType::Raf,
    RawImageType::Raw3fr,
    RawImageType::RawContaxN,
    RawImageType::Rdc,
    RawImageType::Rw2,
    RawImageType::Srw,
    RawImageType::X3f,
];

#[test]
fn test_all_24_raw_formats_positive_recognition() {
    let test_cases: [(RawImageType, Vec<u8>); 24] = [
        (RawImageType::Arw, make_arw_header(true)),
        (RawImageType::Cr2, make_cr2_header(true)),
        (RawImageType::Cr3, make_cr3_header()),
        (RawImageType::Crw, make_crw_header(true)),
        (RawImageType::Dcr, make_dcr_header(true, true)),
        (RawImageType::Dng, make_dng_header(true)),
        (RawImageType::Erf, make_erf_header(true)),
        (RawImageType::Kdc, make_kdc_header(true)),
        (RawImageType::Mos, make_mos_header(true)),
        (RawImageType::Mrw, make_mrw_header()),
        (RawImageType::Nef, make_nef_header(true)),
        (RawImageType::Nrw, make_nrw_header(true)),
        (RawImageType::Orf, make_orf_header(true, true)),
        (RawImageType::Pef, make_pef_header(true, true)),
        (RawImageType::PhaseOne, make_phase_one_header(0)),
        (RawImageType::Pxn, make_pxn_header()),
        (RawImageType::Qtk, make_qtk_header(false)),
        (RawImageType::Raf, make_raf_header()),
        (RawImageType::Raw3fr, make_raw_3fr_header(true)),
        (RawImageType::RawContaxN, make_raw_contax_n_header()),
        (RawImageType::Rdc, make_rdc_header()),
        (RawImageType::Rw2, make_rw2_header(true)),
        (RawImageType::Srw, make_srw_header(true)),
        (RawImageType::X3f, make_x3f_header()),
    ];

    for (expected_type, header) in &test_cases {
        let recognized = recognize_raw_image_type_lite(header);
        assert_eq!(recognized, *expected_type, "Failed recognizing {:?}", expected_type);
        assert!(is_of_type(header, *expected_type), "is_of_type failed for {:?}", expected_type);
        assert!(is_raw_lite(header), "is_raw_lite returned false for {:?}", expected_type);
        assert!(
            is_raw(*expected_type),
            "is_raw returned false for enum variant {:?}",
            expected_type
        );

        // Disambiguation check: ensure it does not match other RAW types
        for other_type in ALL_24_TYPES {
            if other_type != *expected_type {
                assert!(
                    !is_of_type(header, other_type),
                    "Header for {:?} incorrectly matched {:?}",
                    expected_type,
                    other_type
                );
            }
        }
    }
}

#[test]
fn test_endian_and_signature_variants() {
    // ARW BE
    let arw_be = make_arw_header(false);
    assert_eq!(recognize_raw_image_type_lite(&arw_be), RawImageType::Arw);

    // CR2 BE
    let cr2_be = make_cr2_header(false);
    assert_eq!(recognize_raw_image_type_lite(&cr2_be), RawImageType::Cr2);

    // CRW BE
    let crw_be = make_crw_header(false);
    assert_eq!(recognize_raw_image_type_lite(&crw_be), RawImageType::Crw);

    // DCR Tag search variant
    let dcr_tags = make_dcr_header(false, false);
    assert_eq!(recognize_raw_image_type_lite(&dcr_tags), RawImageType::Dcr);

    // DNG BE
    let dng_be = make_dng_header(false);
    assert_eq!(recognize_raw_image_type_lite(&dng_be), RawImageType::Dng);

    // ERF BE
    let erf_be = make_erf_header(false);
    assert_eq!(recognize_raw_image_type_lite(&erf_be), RawImageType::Erf);

    // KDC BE
    let kdc_be = make_kdc_header(false);
    assert_eq!(recognize_raw_image_type_lite(&kdc_be), RawImageType::Kdc);

    // MOS BE
    let mos_be = make_mos_header(false);
    assert_eq!(recognize_raw_image_type_lite(&mos_be), RawImageType::Mos);

    // NEF BE
    let nef_be = make_nef_header(false);
    assert_eq!(recognize_raw_image_type_lite(&nef_be), RawImageType::Nef);

    // NRW BE
    let nrw_be = make_nrw_header(false);
    assert_eq!(recognize_raw_image_type_lite(&nrw_be), RawImageType::Nrw);

    // ORF variants
    let orf_sr = make_orf_header(true, false);
    assert_eq!(recognize_raw_image_type_lite(&orf_sr), RawImageType::Orf);

    let orf_om = make_orf_om_digital_header(true, true);
    assert_eq!(recognize_raw_image_type_lite(&orf_om), RawImageType::Orf);

    let orf_om_be = make_orf_om_digital_header(false, false);
    assert_eq!(recognize_raw_image_type_lite(&orf_om_be), RawImageType::Orf);

    // PEF AOC signature
    let pef_aoc = make_pef_header(true, false);
    assert_eq!(recognize_raw_image_type_lite(&pef_aoc), RawImageType::Pef);

    // Phase One at offset 8
    let phase_one_offset8 = make_phase_one_header(8);
    assert_eq!(recognize_raw_image_type_lite(&phase_one_offset8), RawImageType::PhaseOne);

    // QTK qktn variant
    let qtk_n = make_qtk_header(true);
    assert_eq!(recognize_raw_image_type_lite(&qtk_n), RawImageType::Qtk);

    // RAW 3FR BE
    let raw3fr_be = make_raw_3fr_header(false);
    assert_eq!(recognize_raw_image_type_lite(&raw3fr_be), RawImageType::Raw3fr);

    // RW2 BE
    let rw2_be = make_rw2_header(false);
    assert_eq!(recognize_raw_image_type_lite(&rw2_be), RawImageType::Rw2);

    // SRW BE
    let srw_be = make_srw_header(false);
    assert_eq!(recognize_raw_image_type_lite(&srw_be), RawImageType::Srw);
}

#[test]
fn test_orf_om_digital_signature_recognition() {
    let orf_om_header = make_orf_om_digital_header(true, true);
    assert_eq!(recognize_raw_image_type_lite(&orf_om_header), RawImageType::Orf);
    assert!(is_of_type(&orf_om_header, RawImageType::Orf));
    assert!(is_raw_lite(&orf_om_header));
}

#[test]
fn test_negative_cases_non_raw_media_and_arbitrary_streams() {
    let png_header = b"\x89PNG\r\n\x1a\n\x00\x00\x00\rIHDR\x00\x00\x00\x10\x00\x00\x00\x10\x08\x06";
    let jpeg_header = b"\xFF\xD8\xFF\xE0\x00\x10JFIF\x00\x01\x01\x01\x00\x48\x00\x48\x00\x00";
    let gif87a_header = b"GIF87a\x10\x00\x10\x00\xf0\x00\x00";
    let gif89a_header = b"GIF89a\x10\x00\x10\x00\xf0\x00\x00";
    let bmp_header = b"BM\x36\x00\x00\x00\x00\x00\x00\x00\x36\x00\x00\x00";
    let webp_header = b"RIFF\x38\x00\x00\x00WEBPVP8 \x2c\x00\x00\x00";
    let generic_tiff_header = b"II\x2A\x00\x08\x00\x00\x00\x00\x00\x00\x00\x00\x00";

    let non_raw_samples: &[(&str, &[u8])] = &[
        ("PNG", png_header),
        ("JPEG", jpeg_header),
        ("GIF87a", gif87a_header),
        ("GIF89a", gif89a_header),
        ("BMP", bmp_header),
        ("WEBP", webp_header),
        ("Generic TIFF", generic_tiff_header),
    ];

    for (label, sample) in non_raw_samples {
        assert_eq!(
            recognize_raw_image_type_lite(sample),
            RawImageType::NonRaw,
            "Sample {} was expected to be NonRaw",
            label
        );
        assert!(!is_raw_lite(sample), "is_raw_lite should return false for {}", label);
        assert!(
            !is_of_type(sample, RawImageType::NonRaw),
            "is_of_type(..., NonRaw) should be false for {}",
            label
        );
    }

    // Zeroed buffers
    for size in [4, 16, 256, 1024, 5000, 10000, 20000] {
        let buf = vec![0u8; size];
        assert_eq!(recognize_raw_image_type_lite(&buf), RawImageType::NonRaw);
        assert!(!is_raw_lite(&buf));
    }

    // 0xFF filled buffers
    for size in [4, 16, 256, 1024, 5000, 10000] {
        let buf = vec![0xFFu8; size];
        assert_eq!(recognize_raw_image_type_lite(&buf), RawImageType::NonRaw);
        assert!(!is_raw_lite(&buf));
    }

    // Arbitrary pattern buffer (0x01, 0x02, 0x03, ...)
    let pattern_buf: Vec<u8> = (0..10_000).map(|i| (i % 251) as u8).collect();
    assert_eq!(recognize_raw_image_type_lite(&pattern_buf), RawImageType::NonRaw);
    assert!(!is_raw_lite(&pattern_buf));

    assert!(!is_raw(RawImageType::NonRaw));
}

#[test]
fn test_format_priority_and_dng_precedence() {
    // 1. Hasselblad 3FR with DNG tags present (3FR precedes DNG)
    let mut raw3fr_dng = make_raw_3fr_header(true);
    raw3fr_dng[10..18].copy_from_slice(b"\x12\xc6\x01\x00\x04\x00\x00\x00");
    raw3fr_dng[30..38].copy_from_slice(b"\x13\xc6\x01\x00\x04\x00\x00\x00");
    assert_eq!(recognize_raw_image_type_lite(&raw3fr_dng), RawImageType::Raw3fr);

    // 2. Epson ERF with DNG tags present (ERF precedes DNG)
    let mut erf_dng = make_erf_header(true);
    erf_dng[10..18].copy_from_slice(b"\x12\xc6\x01\x00\x04\x00\x00\x00");
    erf_dng[30..38].copy_from_slice(b"\x13\xc6\x01\x00\x04\x00\x00\x00");
    assert_eq!(recognize_raw_image_type_lite(&erf_dng), RawImageType::Erf);

    // 3. Sony ARW with DNG tags present (DNG has smaller requested_size 1024 < 10000)
    let mut arw_dng = make_arw_header(true);
    arw_dng[10..18].copy_from_slice(b"\x12\xc6\x01\x00\x04\x00\x00\x00");
    arw_dng[30..38].copy_from_slice(b"\x13\xc6\x01\x00\x04\x00\x00\x00");
    assert_eq!(recognize_raw_image_type_lite(&arw_dng), RawImageType::Dng);

    // 4. Nikon NEF with DNG tags present (DNG has smaller requested_size 1024 < 4000)
    let mut nef_dng = make_nef_header(true);
    nef_dng[10..18].copy_from_slice(b"\x12\xc6\x01\x00\x04\x00\x00\x00");
    nef_dng[30..38].copy_from_slice(b"\x13\xc6\x01\x00\x04\x00\x00\x00");
    // Ensure "NIKON" is still present elsewhere in range
    nef_dng[50..55].copy_from_slice(b"NIKON");
    assert_eq!(recognize_raw_image_type_lite(&nef_dng), RawImageType::Dng);

    // 5. Pentax PEF with DNG tags present (DNG has smaller requested_size 1024 < 1280)
    let mut pef_dng = make_pef_header(true, true);
    pef_dng[10..18].copy_from_slice(b"\x12\xc6\x01\x00\x04\x00\x00\x00");
    pef_dng[30..38].copy_from_slice(b"\x13\xc6\x01\x00\x04\x00\x00\x00");
    assert_eq!(recognize_raw_image_type_lite(&pef_dng), RawImageType::Dng);

    // 6. Canon CR2 with DNG tags present (CR2 has smaller requested_size 16 < 1024)
    let mut cr2_dng = make_cr2_header(true);
    cr2_dng.resize(1024, 0);
    cr2_dng[50..58].copy_from_slice(b"\x12\xc6\x01\x00\x04\x00\x00\x00");
    cr2_dng[70..78].copy_from_slice(b"\x13\xc6\x01\x00\x04\x00\x00\x00");
    assert_eq!(recognize_raw_image_type_lite(&cr2_dng), RawImageType::Cr2);

    // 7. Nikon NRW vs NEF precedence
    let nrw_header = make_nrw_header(true);
    assert_eq!(recognize_raw_image_type_lite(&nrw_header), RawImageType::Nrw);
}

#[test]
fn test_boundary_and_slice_safety() {
    // Empty slice
    assert_eq!(recognize_raw_image_type_lite(&[]), RawImageType::NonRaw);
    assert!(!is_raw_lite(&[]));

    // Single-byte slices
    for byte in [0x00, 0x49, 0x4D, 0x2A, 0x55, 0x89, 0xFF] {
        assert_eq!(recognize_raw_image_type_lite(&[byte]), RawImageType::NonRaw);
        assert!(!is_raw_lite(&[byte]));
    }

    // 2-byte slices
    for prefix in [b"II".as_slice(), b"MM", b"qk", b"FU", b"BM", b"CR", b"PX"] {
        assert_eq!(recognize_raw_image_type_lite(prefix), RawImageType::NonRaw);
        assert!(!is_raw_lite(prefix));
    }

    // 3-byte slices
    for prefix in [b"MRM".as_slice(), b"\x00MR", b"FOV", b"CR2", b"fty", b"PXN", b"DSC"] {
        assert_eq!(recognize_raw_image_type_lite(prefix), RawImageType::NonRaw);
        assert!(!is_raw_lite(prefix));
    }

    // Truncated headers for every RAW format (1 byte shorter than min requested)
    let arw_trunc = &make_arw_header(true)[..199];
    assert_eq!(recognize_raw_image_type_lite(arw_trunc), RawImageType::NonRaw);

    let cr2_trunc = &make_cr2_header(true)[..11];
    assert_eq!(recognize_raw_image_type_lite(cr2_trunc), RawImageType::NonRaw);

    let cr3_trunc = &make_cr3_header()[..11];
    assert_eq!(recognize_raw_image_type_lite(cr3_trunc), RawImageType::NonRaw);

    let crw_trunc = &make_crw_header(true)[..13];
    assert_eq!(recognize_raw_image_type_lite(crw_trunc), RawImageType::NonRaw);

    let erf_trunc = &make_erf_header(true)[..104];
    assert_eq!(recognize_raw_image_type_lite(erf_trunc), RawImageType::NonRaw);

    let mrw_trunc = &make_mrw_header()[..3];
    assert_eq!(recognize_raw_image_type_lite(mrw_trunc), RawImageType::NonRaw);

    let phase_one_trunc = &make_phase_one_header(0)[..7];
    assert_eq!(recognize_raw_image_type_lite(phase_one_trunc), RawImageType::NonRaw);

    let pxn_trunc = &make_pxn_header()[..3];
    assert_eq!(recognize_raw_image_type_lite(pxn_trunc), RawImageType::NonRaw);

    let qtk_trunc = &make_qtk_header(false)[..7];
    assert_eq!(recognize_raw_image_type_lite(qtk_trunc), RawImageType::NonRaw);

    let raf_trunc = &make_raf_header()[..7];
    assert_eq!(recognize_raw_image_type_lite(raf_trunc), RawImageType::NonRaw);

    let rdc_trunc = &make_rdc_header()[..15];
    assert_eq!(recognize_raw_image_type_lite(rdc_trunc), RawImageType::NonRaw);

    let rw2_trunc = &make_rw2_header(true)[..3];
    assert_eq!(recognize_raw_image_type_lite(rw2_trunc), RawImageType::NonRaw);

    let x3f_trunc = &make_x3f_header()[..3];
    assert_eq!(recognize_raw_image_type_lite(x3f_trunc), RawImageType::NonRaw);

    // Oversized buffer (> MAX_HEADER_INSPECTION_SIZE)
    let mut huge_arw = make_arw_header(true);
    huge_arw.resize(25_000, 0xEE);
    assert_eq!(recognize_raw_image_type_lite(&huge_arw), RawImageType::Arw);
    assert!(is_raw_lite(&huge_arw));
}

#[test]
fn test_header_size_functions() {
    assert_eq!(get_number_of_bytes_for_is_raw_lite(), MAX_HEADER_INSPECTION_SIZE);

    let expected_sizes: [(RawImageType, usize); 25] = [
        (RawImageType::Mrw, 4),
        (RawImageType::Pxn, 4),
        (RawImageType::Rw2, 4),
        (RawImageType::X3f, 4),
        (RawImageType::Qtk, 8),
        (RawImageType::Raf, 8),
        (RawImageType::Cr3, 12),
        (RawImageType::Crw, 14),
        (RawImageType::Cr2, 16),
        (RawImageType::PhaseOne, 16),
        (RawImageType::Rdc, 16),
        (RawImageType::RawContaxN, 36),
        (RawImageType::Srw, 256),
        (RawImageType::Erf, 1024),
        (RawImageType::Raw3fr, 1024),
        (RawImageType::Dng, 1024),
        (RawImageType::Pef, 1280),
        (RawImageType::Orf, 3000),
        (RawImageType::Nef, 4000),
        (RawImageType::Nrw, 4000),
        (RawImageType::Dcr, 5000),
        (RawImageType::Kdc, 5000),
        (RawImageType::Mos, 5000),
        (RawImageType::Arw, 10_000),
        (RawImageType::NonRaw, 0),
    ];

    for (raw_type, expected_bytes) in expected_sizes {
        let size = get_number_of_bytes_for_is_of_type(raw_type);
        assert_eq!(size, expected_bytes, "Mismatch requested byte size for {:?}", raw_type);
    }
}

#[test]
fn test_format_recognizer_extreme_offsets() {
    let sample = b"II\x2a\x00\x08\x00\x00\x00";

    // Valid reads
    assert_eq!(get_u16(sample, 0, Endian::Little), Some(0x4949));
    assert_eq!(get_u32(sample, 4, Endian::Little), Some(8));

    // Boundary and overflow reads
    assert_eq!(get_u16(sample, usize::MAX, Endian::Little), None);
    assert_eq!(get_u16(sample, usize::MAX - 1, Endian::Little), None);
    assert_eq!(get_u32(sample, usize::MAX, Endian::Little), None);
    assert_eq!(get_u32(sample, usize::MAX - 3, Endian::Little), None);

    // Signature matching with overflow offsets
    assert!(!is_signature_matched(sample, usize::MAX, b"II"));
    assert!(!is_signature_matched(sample, usize::MAX - 1, b"II\x2a"));
    assert!(!is_signature_found(sample, usize::MAX, 10, b"II"));
    assert!(!is_signature_found(sample, 0, usize::MAX, b"II"));
}

#[test]
fn test_regression_crash_01_is_raw_short_buffers() {
    // Truncated buffer with magic signature bytes (5 bytes long)
    let short_x3f = b"\x03FOVb";
    let mut cursor = Cursor::new(short_x3f.as_slice());

    // Should not panic and should consistently return false for truncated header
    let is_raw_res = piex::extractor::is_raw(&mut cursor);
    let is_raw_lite_res = is_raw_lite(short_x3f);
    let type_res = recognize_raw_image_type_lite(short_x3f);

    assert_eq!(type_res, RawImageType::NonRaw);
    assert!(!is_raw_res);
    assert!(!is_raw_lite_res);
}

#[test]
fn test_checker_precedence_ordering() {
    // Buffer with Big-Endian TIFF header, PENTAX make string, and 2 DNG tag signatures.
    // In C++, PefTypeChecker (1280) is evaluated AFTER Erf/Raw3fr/Dng (1024), but
    // for identical requested size formats (like PEF vs SRW vs DNG), specific camera
    // formats should take precedence.
    let mut pef_dng_buf = vec![0u8; 1280];
    pef_dng_buf[0..2].copy_from_slice(b"MM");
    pef_dng_buf[2..4].copy_from_slice(&0x2Au16.to_be_bytes());
    pef_dng_buf[4..8].copy_from_slice(&8u32.to_be_bytes());
    pef_dng_buf[16..24].copy_from_slice(b"PENTAX \x00");
    // Include 2 DNG tags (50706 and 50707)
    pef_dng_buf[100..108].copy_from_slice(b"\xc6\x12\x00\x01\x00\x00\x00\x04");
    pef_dng_buf[200..208].copy_from_slice(b"\xc6\x13\x00\x01\x00\x00\x00\x04");

    // Buffer with Big-Endian TIFF header, PENTAX make string, and Samsung SubIFD tag.
    // In C++, PefTypeChecker is evaluated BEFORE SrwTypeChecker when sizes match.
    let mut erf_dng_buf = vec![0u8; 1024];
    erf_dng_buf[0..2].copy_from_slice(b"II");
    erf_dng_buf[2..4].copy_from_slice(&0x2Au16.to_le_bytes());
    erf_dng_buf[4..8].copy_from_slice(&8u32.to_le_bytes());
    erf_dng_buf[16..22].copy_from_slice(b"EPSON\x00");
    erf_dng_buf[100..108].copy_from_slice(b"\x12\xc6\x01\x00\x04\x00\x00\x00");
    erf_dng_buf[200..208].copy_from_slice(b"\x13\xc6\x01\x00\x04\x00\x00\x00");

    assert_eq!(
        recognize_raw_image_type_lite(&erf_dng_buf),
        RawImageType::Erf,
        "ERF checker (1024) must take precedence over DNG checker (1024) per C++ TypeCheckerList"
    );
}

#[test]
fn test_type_check_window_ignores_bytes_past_requested_size() {
    // The Epson checker inspects the first 1024 bytes and accepts "EPSON"
    // starting at offsets 0..=1018.
    for &(file_len, signature_pos, expected) in &[
        (1040usize, 1018usize, RawImageType::Erf),
        (1040, 1019, RawImageType::NonRaw),
        (1040, 1022, RawImageType::NonRaw),
        (1040, 1024, RawImageType::NonRaw),
        (1024, 1018, RawImageType::Erf),
        (1023, 100, RawImageType::NonRaw),
        (1056, 1018, RawImageType::Erf),
        (1056, 1019, RawImageType::NonRaw),
        (1056, 1040, RawImageType::NonRaw),
    ] {
        let mut file = vec![0u8; file_len];
        file[..4].copy_from_slice(b"II*\0");
        file[signature_pos..signature_pos + 5].copy_from_slice(b"EPSON");
        let case = format!("file_len={file_len} signature_pos={signature_pos}");
        assert_eq!(recognize_raw_image_type_lite(&file), expected, "{case}");
        assert_eq!(is_of_type(&file, RawImageType::Erf), expected == RawImageType::Erf, "{case}");
    }
}
