#![forbid(unsafe_code)]

//! TIFF container parser implementation.
//!
//! Provides the `TiffParser` struct for traversing TIFF Image File Directories (IFDs),
//! parsing JPEG preview metadata, scanning JPEG stream markers for image dimensions,
//! and constructing `PreviewImageData` structures.

use crate::stream::{Endian, StreamInterface};
use crate::tiff::directory::{DirectoryEntry, TiffContent, TiffDirectory, TiffError, TiffType};
use crate::tiff::tags::gps::*;
use crate::tiff::tags::panasonic::*;
use crate::tiff::tags::tiff::*;
use crate::types::{ColorSpace, Image, ImageFormat, PreviewImageData, Rational};

/// Maximum pixel length along any dimension for an embedded image to be classified as a thumbnail.
pub const THUMBNAIL_MAX_DIMENSION: u32 = 512;

/// Enables parsing through data that complies with the TIFF/EP specification.
pub struct TiffParser<'a, S: StreamInterface + ?Sized> {
    stream: &'a S,
    tiff_offset: u32,
    endian: Endian,
}

impl<'a, S: StreamInterface + ?Sized> TiffParser<'a, S> {
    /// Creates a new `TiffParser` bound to `stream` with header starting at offset 0.
    ///
    /// Determines byte endianness and validates TIFF magic number from header bytes.
    pub fn new(stream: &'a S) -> Result<Self, TiffError> {
        Self::with_offset(stream, 0)
    }

    /// Creates a new `TiffParser` bound to `stream` with TIFF header at specified `tiff_offset`.
    ///
    /// Determines byte endianness and validates TIFF magic number at `tiff_offset`.
    pub fn with_offset(stream: &'a S, tiff_offset: u32) -> Result<Self, TiffError> {
        let endian = Self::get_endianness(stream, tiff_offset)?;
        Ok(Self { stream, tiff_offset, endian })
    }

    /// Creates a new `TiffParser` with specified `tiff_offset` and `endian`
    /// without re-checking magic bytes.
    pub fn with_offset_and_endian(stream: &'a S, tiff_offset: u32, endian: Endian) -> Self {
        Self { stream, tiff_offset, endian }
    }

    /// Returns the detected endianness of this parser instance.
    pub fn endian(&self) -> Endian {
        self.endian
    }

    /// Returns the base TIFF offset in bytes.
    pub fn tiff_offset(&self) -> u32 {
        self.tiff_offset
    }

    /// Retrieves byte order and validates magic number at `tiff_offset` in `stream`.
    ///
    /// Reads 4 bytes from `tiff_offset`:
    /// - Bytes 0..2: Byte order marker (`"II"` for Little-Endian, `"MM"` for Big-Endian).
    /// - Bytes 2..4: TIFF magic number (expected `42` / `0x002A`).
    pub fn get_endianness(stream: &S, tiff_offset: u32) -> Result<Endian, TiffError> {
        let mut header = [0u8; 4];
        stream.get_data(tiff_offset as usize, &mut header).map_err(TiffError::StreamError)?;

        let endian = match &header[0..2] {
            b"II" => Endian::Little,
            b"MM" => Endian::Big,
            _ => return Err(TiffError::InvalidByteOrder),
        };

        Ok(endian)
    }

    /// Parses an Image File Directory (IFD) at `ifd_offset`.
    ///
    /// Extracts directory entries matching `desired_tags` and populates `next_ifd_offset`.
    pub fn parse_directory(
        &self,
        ifd_offset: u32,
        desired_tags: &[u32],
        next_ifd_offset: &mut u32,
    ) -> Result<TiffDirectory, TiffError> {
        let number_of_entries = self.read_u16(ifd_offset)?;
        let mut directory = TiffDirectory::new(self.endian);

        let entries_base_offset =
            ifd_offset.checked_add(2).ok_or(TiffError::CorruptedEntry { offset: ifd_offset })?;

        for i in 0..number_of_entries as u32 {
            let entry_offset = entries_base_offset
                .checked_add(i.saturating_mul(12))
                .ok_or(TiffError::CorruptedEntry { offset: ifd_offset })?;

            let tag = self.read_u16(entry_offset)?;
            let type_val_offset = entry_offset
                .checked_add(2)
                .ok_or(TiffError::CorruptedEntry { offset: entry_offset })?;
            let type_val = self.read_u16(type_val_offset)?;

            let count_offset = entry_offset
                .checked_add(4)
                .ok_or(TiffError::CorruptedEntry { offset: entry_offset })?;
            let number_of_elements = self.read_u32(count_offset)?;

            let tag_u32 = tag as u32;
            if !desired_tags.is_empty() && !desired_tags.contains(&tag_u32) {
                continue;
            }

            let type_ = match TiffType::from_u16(type_val) {
                Some(t) => t,
                None => continue,
            };
            let type_size = type_.element_size();
            if type_size == 0 {
                continue;
            }

            if type_size != 0 && number_of_elements > (u32::MAX / type_size as u32) {
                return Err(TiffError::CorruptedEntry { offset: entry_offset });
            }

            let byte_count = (type_size as usize) * (number_of_elements as usize);

            let val_ptr_offset = entry_offset
                .checked_add(8)
                .ok_or(TiffError::CorruptedEntry { offset: entry_offset })?;

            let value_offset = if byte_count > 4 {
                let val_offset = self.read_u32(val_ptr_offset)?;
                self.tiff_offset
                    .checked_add(val_offset)
                    .ok_or(TiffError::CorruptedEntry { offset: entry_offset })?
            } else if byte_count != 0 {
                val_ptr_offset
            } else {
                continue;
            };

            let remaining = self.stream.len().saturating_sub(value_offset as usize);
            if byte_count > remaining {
                return Err(TiffError::CorruptedEntry { offset: value_offset });
            }

            if byte_count <= 8 {
                let mut buf = [0u8; 8];
                self.stream
                    .get_data(value_offset as usize, &mut buf[..byte_count])
                    .map_err(TiffError::StreamError)?;
                directory.add_directory_entry(
                    tag_u32,
                    DirectoryEntry::new_slice(
                        type_,
                        number_of_elements,
                        value_offset,
                        &buf[..byte_count],
                    ),
                );
            } else {
                let mut data = vec![0u8; byte_count];
                self.stream
                    .get_data(value_offset as usize, &mut data)
                    .map_err(TiffError::StreamError)?;
                directory.add_directory_entry(
                    tag_u32,
                    DirectoryEntry::new(type_, number_of_elements, value_offset, data),
                );
            }
        }

        let next_offset_pos = entries_base_offset
            .checked_add((number_of_entries as u32).saturating_mul(12))
            .ok_or(TiffError::CorruptedEntry { offset: ifd_offset })?;

        *next_ifd_offset = self.read_u32(next_offset_pos)?;

        Ok(directory)
    }

    /// Parses the primary TIFF IFD chain, EXIF sub-directory, and GPS sub-directory.
    ///
    /// Reads initial IFD offset at offset 4 of the TIFF header and traverses up to
    /// `max_number_ifds` directories.
    pub fn parse(
        &self,
        desired_tags: &[u32],
        max_number_ifds: u16,
    ) -> Result<TiffContent, TiffError> {
        let initial_offset_pos = self
            .tiff_offset
            .checked_add(4)
            .ok_or(TiffError::CorruptedEntry { offset: self.tiff_offset })?;
        let initial_ifd_offset = self.read_u32(initial_offset_pos)?;
        let mut content = TiffContent::new();
        let mut visited_offsets = Vec::new();

        let first_ifd_offset = self
            .tiff_offset
            .checked_add(initial_ifd_offset)
            .ok_or(TiffError::CorruptedEntry { offset: self.tiff_offset })?;

        self.parse_ifd_chain(
            first_ifd_offset,
            desired_tags,
            max_number_ifds,
            &mut content.tiff_directories,
            &mut visited_offsets,
        )?;

        // Get the Exif data.
        if let Some(tiff_ifd) =
            self.find_first_tag_in_ifds(EXIF_TAG_EXIF_IFD, &content.tiff_directories)
        {
            if let Ok(offset) = tiff_ifd.get_u32(EXIF_TAG_EXIF_IFD) {
                let gps_offset = tiff_ifd.get_u32(EXIF_TAG_GPS).ok();
                let exif_offset = self
                    .tiff_offset
                    .checked_add(offset)
                    .ok_or(TiffError::CorruptedEntry { offset: self.tiff_offset })?;
                let mut next_offset = 0u32;
                let exif_dir = self.parse_directory(exif_offset, desired_tags, &mut next_offset)?;
                content.exif_directory = Some(exif_dir);

                if let Some(gps_off) = gps_offset {
                    self.parse_gps_data(gps_off, &mut content)?;
                }
                return Ok(content);
            }
        }

        // Get the GPS data from the tiff ifd.
        if let Some(tiff_ifd) = self.find_first_tag_in_ifds(EXIF_TAG_GPS, &content.tiff_directories)
        {
            if let Ok(gps_off) = tiff_ifd.get_u32(EXIF_TAG_GPS) {
                self.parse_gps_data(gps_off, &mut content)?;
            }
            return Ok(content);
        }

        Ok(content)
    }

    /// Retrieves preview image location and dimension metadata from `dir`.
    pub fn get_image_data(&self, dir: &TiffDirectory) -> Result<Image, TiffError> {
        let length: u32;
        let offset: u32;
        let format;

        if dir.has(TIFF_TAG_JPEG_OFFSET) && dir.has(TIFF_TAG_JPEG_BYTE_COUNT) {
            offset = dir.get_u32(TIFF_TAG_JPEG_OFFSET)?;
            length = dir.get_u32(TIFF_TAG_JPEG_BYTE_COUNT)?;
            format = ImageFormat::JpegCompressed;
        } else if dir.has(TIFF_TAG_STRIP_OFFSETS) && dir.has(TIFF_TAG_STRIP_BYTE_COUNTS) {
            let strip_offsets = dir.get_u32_vec(TIFF_TAG_STRIP_OFFSETS)?;
            let strip_byte_counts = dir.get_u32_vec(TIFF_TAG_STRIP_BYTE_COUNTS)?;

            if !self.offsets_are_consecutive(&strip_offsets, &strip_byte_counts) {
                return Err(TiffError::TagNotFound { tag: TIFF_TAG_STRIP_OFFSETS });
            }

            let compression = dir.get_u32(TIFF_TAG_COMPRESSION)?;
            if dir.has(TIFF_TAG_PHOTOMETRIC) {
                let photometric = dir.get_u32(TIFF_TAG_PHOTOMETRIC)?;
                if photometric != 2 && photometric != 6 {
                    return Err(TiffError::TagNotFound { tag: TIFF_TAG_PHOTOMETRIC });
                }
            }

            format = match compression {
                1 => ImageFormat::UncompressedRgb,
                6 | 7 => ImageFormat::JpegCompressed,
                _ => return Err(TiffError::TagNotFound { tag: TIFF_TAG_COMPRESSION }),
            };

            length = strip_byte_counts
                .iter()
                .try_fold(0u32, |acc, &count| acc.checked_add(count))
                .ok_or(TiffError::CorruptedEntry {
                    offset: *strip_offsets.first().unwrap_or(&0),
                })?;
            offset = *strip_offsets.first().unwrap_or(&0);
        } else if dir.has(PANASONIC_TAG_JPEG_IMAGE) {
            let (off, len) =
                dir.get_offset_and_length(PANASONIC_TAG_JPEG_IMAGE, TiffType::Undefined)?;
            offset = off;
            length = len;
            format = ImageFormat::JpegCompressed;
        } else {
            return Err(TiffError::TagNotFound { tag: TIFF_TAG_JPEG_OFFSET });
        }

        let mut image = Image::new(0, 0, length, offset, format);
        self.fill_image_size(dir, &mut image)?;
        Ok(image)
    }

    /// Scans JPEG segment markers from `jpeg_offset` to extract width and height.
    pub fn get_jpeg_dimensions(&self, jpeg_offset: u32) -> Result<(u16, u16), TiffError> {
        let mut offset = jpeg_offset;

        loop {
            let segment = self.read_u16_be(offset)?;
            offset = offset.checked_add(2).ok_or(TiffError::CorruptedEntry { offset })?;

            match segment {
                0xFFD8 => continue, // Start of Image (SOI)
                0xFFC0 => {
                    // Start of Frame (SOF0: Baseline DCT only)
                    let height_offset =
                        offset.checked_add(3).ok_or(TiffError::CorruptedEntry { offset })?;
                    let height = self.read_u16_be(height_offset)?;

                    let width_offset =
                        offset.checked_add(5).ok_or(TiffError::CorruptedEntry { offset })?;
                    let width = self.read_u16_be(width_offset)?;
                    return Ok((width, height));
                }
                0xFFDA => break, // Start of Scan (SOS) - image payload begins
                _ => {
                    let length = self.read_u16_be(offset)?;
                    if length < 2 {
                        return Err(TiffError::CorruptedEntry { offset });
                    }
                    offset = offset
                        .checked_add(length as u32)
                        .ok_or(TiffError::CorruptedEntry { offset })?;
                }
            }
        }

        Err(TiffError::TagNotFound { tag: 0xFFC0 })
    }

    /// Checks if `image` dimensions qualify as a thumbnail.
    pub fn is_thumbnail(&self, image: &Image) -> bool {
        image.width <= THUMBNAIL_MAX_DIMENSION as u16
            && image.height <= THUMBNAIL_MAX_DIMENSION as u16
    }

    /// Extracts preview images and metadata from `content` into `preview_data`.
    pub fn get_preview_image_data(
        &self,
        content: &TiffContent,
        preview_data: &mut PreviewImageData,
    ) -> Result<(), TiffError> {
        let mut success = true;
        for tiff_directory in &content.tiff_directories {
            success = self.fill_preview_image_data(tiff_directory, preview_data);

            if success && tiff_directory.has(EXIF_TAG_EXIF_IFD) {
                if let Some(exif_dir) = &content.exif_directory {
                    success = self.fill_preview_image_data(exif_dir, preview_data);
                }
            }

            if success && tiff_directory.has(EXIF_TAG_GPS) {
                if let Some(gps_dir) = &content.gps_directory {
                    self.extract_gps_metadata(gps_dir, preview_data);
                }
            }

            for sub_directory in tiff_directory.sub_directories() {
                if success {
                    success = self.fill_preview_image_data(sub_directory, preview_data);
                }
            }
        }

        if success {
            Ok(())
        } else {
            Err(TiffError::CorruptedEntry { offset: self.tiff_offset })
        }
    }

    /// Helper method extracting EXIF metadata tags from an IFD directory into `PreviewImageData`.
    pub fn extract_exif_metadata(
        &self,
        dir: &TiffDirectory,
        preview_data: &mut PreviewImageData,
    ) -> bool {
        let mut success = true;
        if dir.has(EXIF_TAG_DATE_TIME_ORIGINAL) {
            if let Ok(dt) = dir.get_ascii_bytes(EXIF_TAG_DATE_TIME_ORIGINAL) {
                preview_data.date_time = dt.to_vec();
            } else {
                success = false;
            }
        }

        if dir.has(EXIF_TAG_ISO_SPEED) {
            if let Ok(iso) = dir.get_u32(EXIF_TAG_ISO_SPEED) {
                preview_data.iso = iso;
            } else {
                success = false;
            }
        }

        if dir.has(EXIF_TAG_EXPOSURE_TIME) {
            if let Ok(exp) = dir.get_rational(EXIF_TAG_EXPOSURE_TIME) {
                preview_data.exposure_time = exp;
            } else {
                success = false;
            }
        }

        if dir.has(EXIF_TAG_FNUMBER) {
            if let Ok(fn_val) = dir.get_rational(EXIF_TAG_FNUMBER) {
                preview_data.fnumber = fn_val;
            } else {
                success = false;
            }
        }

        if dir.has(EXIF_TAG_FOCAL_LENGTH) {
            if let Ok(fl) = dir.get_rational(EXIF_TAG_FOCAL_LENGTH) {
                preview_data.focal_length = fl;
            } else {
                success = false;
            }
        }

        success
    }

    /// Helper method extracting GPS tags from `gps_dir` into `preview_data.gps`.
    ///
    /// Mirrors C++ `FillGpsPreviewImageData`: fields are written in place one at
    /// a time, so a failure part-way leaves the earlier fields written and
    /// `is_valid == false`.
    pub fn extract_gps_metadata(
        &self,
        gps_dir: &TiffDirectory,
        preview_data: &mut PreviewImageData,
    ) {
        if !(gps_dir.has(GPS_TAG_LATITUDE_REF)
            && gps_dir.has(GPS_TAG_LATITUDE)
            && gps_dir.has(GPS_TAG_LONGITUDE_REF)
            && gps_dir.has(GPS_TAG_LONGITUDE)
            && gps_dir.has(GPS_TAG_TIME_STAMP)
            && gps_dir.has(GPS_TAG_DATE_STAMP))
        {
            return;
        }
        let gps = &mut preview_data.gps;
        gps.is_valid = false;
        let get_rationals = |tag: u32, out: &mut [Rational; 3]| -> bool {
            match gps_dir.get_rational_vec(tag) {
                Ok(v) if v.len() == 3 => {
                    out.copy_from_slice(&v);
                    true
                }
                _ => false,
            }
        };

        let lat_ref = match gps_dir.get_ascii_bytes(GPS_TAG_LATITUDE_REF) {
            Ok([first @ (b'N' | b'S'), ..]) => *first,
            _ => return,
        };
        if !get_rationals(GPS_TAG_LATITUDE, &mut gps.latitude) {
            return;
        }
        gps.latitude_ref = lat_ref;

        let lon_ref = match gps_dir.get_ascii_bytes(GPS_TAG_LONGITUDE_REF) {
            Ok([first @ (b'E' | b'W'), ..]) => *first,
            _ => return,
        };
        if !get_rationals(GPS_TAG_LONGITUDE, &mut gps.longitude) {
            return;
        }
        gps.longitude_ref = lon_ref;

        if !get_rationals(GPS_TAG_TIME_STAMP, &mut gps.time_stamp) {
            return;
        }

        // C++ assigns the full value first and only then checks its size.
        const GPS_DATE_STAMP_SIZE: usize = 11;
        let Ok(date_stamp) = gps_dir.get_ascii_bytes(GPS_TAG_DATE_STAMP) else {
            return;
        };
        if date_stamp.len() != GPS_DATE_STAMP_SIZE {
            gps.date_stamp = date_stamp.to_vec();
            return;
        }
        // Drop the trailing NUL, like C++ `resize(kGpsDateStampSize - 1)`.
        gps.date_stamp = date_stamp[..GPS_DATE_STAMP_SIZE - 1].to_vec();

        if gps_dir.has(GPS_TAG_ALTITUDE_REF) && gps_dir.has(GPS_TAG_ALTITUDE) {
            let altitude_ref = match gps_dir.get_bytes(GPS_TAG_ALTITUDE_REF) {
                Ok(bytes) if !bytes.is_empty() => bytes[0] != 0,
                _ => return,
            };
            match gps_dir.get_rational_vec(GPS_TAG_ALTITUDE) {
                Ok(v) if v.len() == 1 => gps.altitude = v[0],
                _ => return,
            }
            gps.altitude_ref = altitude_ref;
        }

        gps.is_valid = true;
    }

    // --- Internal helper functions ---

    fn parse_ifd_chain(
        &self,
        ifd_offset: u32,
        desired_tags: &[u32],
        max_number_ifds: u16,
        tiff_directories: &mut Vec<TiffDirectory>,
        visited_offsets: &mut Vec<u32>,
    ) -> Result<(), TiffError> {
        if visited_offsets.contains(&ifd_offset) {
            return Ok(());
        }

        if tiff_directories.len() >= max_number_ifds as usize {
            return Ok(());
        }

        visited_offsets.push(ifd_offset);

        let mut next_ifd_offset = 0u32;
        let mut tiff_ifd = self.parse_directory(ifd_offset, desired_tags, &mut next_ifd_offset)?;

        self.parse_sub_ifds(&mut tiff_ifd, desired_tags, max_number_ifds)?;
        tiff_directories.push(tiff_ifd);

        if next_ifd_offset != 0 && tiff_directories.len() < max_number_ifds as usize {
            let next_offset = self
                .tiff_offset
                .checked_add(next_ifd_offset)
                .ok_or(TiffError::CorruptedEntry { offset: ifd_offset })?;
            self.parse_ifd_chain(
                next_offset,
                desired_tags,
                max_number_ifds,
                tiff_directories,
                visited_offsets,
            )?;
        }

        Ok(())
    }

    fn parse_sub_ifds(
        &self,
        tiff_ifd: &mut TiffDirectory,
        desired_tags: &[u32],
        max_number_ifds: u16,
    ) -> Result<(), TiffError> {
        if tiff_ifd.get_offset_and_length(EXIF_TAG_SUB_IFDS, TiffType::Long).is_ok() {
            if let Ok(offsets) = tiff_ifd.get_u32_vec(EXIF_TAG_SUB_IFDS) {
                for &sub_offset in offsets.iter().take(max_number_ifds as usize) {
                    let mut next_offset = 0u32;
                    let sub_dir_offset = self
                        .tiff_offset
                        .checked_add(sub_offset)
                        .ok_or(TiffError::CorruptedEntry { offset: self.tiff_offset })?;
                    let sub_ifd =
                        self.parse_directory(sub_dir_offset, desired_tags, &mut next_offset)?;
                    tiff_ifd.add_sub_directory(sub_ifd);
                }
            }
        }
        Ok(())
    }

    fn parse_gps_data(&self, gps_offset: u32, content: &mut TiffContent) -> Result<(), TiffError> {
        let gps_tags = all_gps_tags();
        let mut next_offset = 0u32;
        let target_offset = self
            .tiff_offset
            .checked_add(gps_offset)
            .ok_or(TiffError::CorruptedEntry { offset: self.tiff_offset })?;
        let gps_dir = self.parse_directory(target_offset, gps_tags, &mut next_offset)?;
        content.gps_directory = Some(gps_dir);
        Ok(())
    }

    fn fill_image_size(&self, dir: &TiffDirectory, image: &mut Image) -> Result<(), TiffError> {
        match image.format {
            ImageFormat::UncompressedRgb => {
                if let Ok((w, h)) = self.get_full_dimension16(dir) {
                    image.width = w;
                    image.height = h;
                }
            }
            ImageFormat::JpegCompressed => {
                if let Ok((w, h)) = self.get_jpeg_dimensions(image.offset) {
                    image.width = w;
                    image.height = h;
                }
            }
            _ => {}
        }
        Ok(())
    }

    fn fill_preview_image_data(
        &self,
        dir: &TiffDirectory,
        preview_data: &mut PreviewImageData,
    ) -> bool {
        let mut success = true;

        if let Ok(image) = self.get_image_data(dir) {
            if self.is_thumbnail(&image) {
                preview_data.thumbnail = image;
            } else if image.format == ImageFormat::JpegCompressed {
                preview_data.preview = image;
            }
        }

        if dir.has(TIFF_TAG_ORIENTATION) && preview_data.exif_orientation == 1 {
            if let Ok(orient) = dir.get_u32(TIFF_TAG_ORIENTATION) {
                preview_data.exif_orientation = orient;
            } else {
                success = false;
            }
        }

        if dir.has(EXIF_TAG_COLOR_SPACE) {
            if let Ok(cs) = dir.get_u32(EXIF_TAG_COLOR_SPACE) {
                if cs == 1 {
                    preview_data.color_space = ColorSpace::Srgb;
                } else if cs == 2 || cs == 65535 {
                    preview_data.color_space = ColorSpace::AdobeRgb;
                }
            } else {
                success = false;
            }
        }

        let is_secondary = dir.has(TIFF_TAG_SUB_FILE_TYPE)
            && dir.get_u32(TIFF_TAG_SUB_FILE_TYPE).map(|v| v != 0).unwrap_or(false);

        if !is_secondary {
            if let Ok((fw, fh)) = self.get_full_dimension32(dir) {
                if fw > preview_data.full_width && fh > preview_data.full_height {
                    preview_data.full_width = fw;
                    preview_data.full_height = fh;
                }
            }
        }

        if dir.has(TIFF_TAG_MAKE) {
            if let Ok(maker) = dir.get_ascii_bytes(TIFF_TAG_MAKE) {
                preview_data.maker = maker.to_vec();
            } else {
                success = false;
            }
        }

        if dir.has(TIFF_TAG_MODEL) {
            if let Ok(model) = dir.get_ascii_bytes(TIFF_TAG_MODEL) {
                preview_data.model = model.to_vec();
            } else {
                success = false;
            }
        }

        if dir.has(TIFF_TAG_CFA_PATTERN_DIM) {
            if let Ok(cfa) = dir.get_u32_vec(TIFF_TAG_CFA_PATTERN_DIM) {
                if cfa.len() == 2 {
                    preview_data.cfa_pattern_dim = [cfa[0], cfa[1]];
                }
            }
        }

        success = self.extract_exif_metadata(dir, preview_data) && success;

        if dir.has(PANASONIC_TAG_ISO) && preview_data.iso == 0 {
            if let Ok(iso) = dir.get_u32(PANASONIC_TAG_ISO) {
                preview_data.iso = iso;
            } else {
                success = false;
            }
        }

        success
    }

    /// Extracts full image dimensions (width, height) in pixels from `dir`.
    pub fn get_full_dimension32(&self, dir: &TiffDirectory) -> Result<(u32, u32), TiffError> {
        if dir.has(EXIF_TAG_DEFAULT_CROP_SIZE) {
            if let Ok(crop) = dir.get_u32_vec(EXIF_TAG_DEFAULT_CROP_SIZE) {
                if crop.len() == 2 && crop[0] > 0 && crop[1] > 0 {
                    return Ok((crop[0], crop[1]));
                }
            }
            if let Ok(crop_rat) = dir.get_rational_vec(EXIF_TAG_DEFAULT_CROP_SIZE) {
                if crop_rat.len() == 2
                    && crop_rat[0].numerator > 0
                    && crop_rat[0].denominator > 0
                    && crop_rat[1].numerator > 0
                    && crop_rat[1].denominator > 0
                {
                    return Ok((
                        crop_rat[0].numerator / crop_rat[0].denominator,
                        crop_rat[1].numerator / crop_rat[1].denominator,
                    ));
                }
            }
            return Err(TiffError::InvalidCount {
                tag: EXIF_TAG_DEFAULT_CROP_SIZE,
                expected: 2,
                actual: 0,
            });
        } else if dir.has(EXIF_TAG_WIDTH) && dir.has(EXIF_TAG_HEIGHT) {
            let w = dir.get_u32(EXIF_TAG_WIDTH)?;
            let h = dir.get_u32(EXIF_TAG_HEIGHT)?;
            return Ok((w, h));
        } else if dir.has(TIFF_TAG_IMAGE_WIDTH) && dir.has(TIFF_TAG_IMAGE_LENGTH) {
            let w = dir.get_u32(TIFF_TAG_IMAGE_WIDTH)?;
            let h = dir.get_u32(TIFF_TAG_IMAGE_LENGTH)?;
            return Ok((w, h));
        } else if dir.has(PANASONIC_TAG_LEFT_BORDER)
            && dir.has(PANASONIC_TAG_RIGHT_BORDER)
            && dir.has(PANASONIC_TAG_TOP_BORDER)
            && dir.has(PANASONIC_TAG_BOTTOM_BORDER)
        {
            let left = dir.get_u32(PANASONIC_TAG_LEFT_BORDER)?;
            let right = dir.get_u32(PANASONIC_TAG_RIGHT_BORDER)?;
            let top = dir.get_u32(PANASONIC_TAG_TOP_BORDER)?;
            let bottom = dir.get_u32(PANASONIC_TAG_BOTTOM_BORDER)?;

            if right > left && bottom > top {
                return Ok((right - left, bottom - top));
            } else {
                return Err(TiffError::CorruptedEntry { offset: 0 });
            }
        }

        Err(TiffError::TagNotFound { tag: TIFF_TAG_IMAGE_WIDTH })
    }

    fn get_full_dimension16(&self, dir: &TiffDirectory) -> Result<(u16, u16), TiffError> {
        let (w, h) = self.get_full_dimension32(dir)?;
        if w <= u16::MAX as u32 && h <= u16::MAX as u32 {
            Ok((w as u16, h as u16))
        } else {
            Err(TiffError::InvalidCount { tag: TIFF_TAG_IMAGE_WIDTH, expected: 2, actual: 0 })
        }
    }

    fn find_first_tag_in_ifds<'b>(
        &self,
        tag: u32,
        ifds: &'b [TiffDirectory],
    ) -> Option<&'b TiffDirectory> {
        for dir in ifds {
            if dir.has(tag) {
                return Some(dir);
            }
            if let Some(sub) = self.find_first_tag_in_ifds(tag, dir.sub_directories()) {
                return Some(sub);
            }
        }
        None
    }

    fn offsets_are_consecutive(&self, strip_offsets: &[u32], strip_byte_counts: &[u32]) -> bool {
        if strip_offsets.is_empty() || strip_offsets.len() != strip_byte_counts.len() {
            return false;
        }
        strip_offsets
            .windows(2)
            .zip(strip_byte_counts.iter())
            .all(|(w, &count)| w[0].checked_add(count) == Some(w[1]))
    }

    fn read_u16(&self, offset: u32) -> Result<u16, TiffError> {
        let mut buf = [0u8; 2];
        self.stream.get_data(offset as usize, &mut buf).map_err(TiffError::StreamError)?;
        Ok(match self.endian {
            Endian::Little => u16::from_le_bytes(buf),
            Endian::Big => u16::from_be_bytes(buf),
        })
    }

    fn read_u32(&self, offset: u32) -> Result<u32, TiffError> {
        let mut buf = [0u8; 4];
        self.stream.get_data(offset as usize, &mut buf).map_err(TiffError::StreamError)?;
        Ok(match self.endian {
            Endian::Little => u32::from_le_bytes(buf),
            Endian::Big => u32::from_be_bytes(buf),
        })
    }

    fn read_u16_be(&self, offset: u32) -> Result<u16, TiffError> {
        let mut buf = [0u8; 2];
        self.stream.get_data(offset as usize, &mut buf).map_err(TiffError::StreamError)?;
        Ok(u16::from_be_bytes(buf))
    }
}
