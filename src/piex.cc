// Copyright 2015 Google Inc.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
////////////////////////////////////////////////////////////////////////////////

#include "src/piex.h"

#include <cstdint>
#include <limits>
#include <set>
#include <vector>

#include "src/binary_parse/range_checked_byte_ptr.h"
#include "src/image_type_recognition/image_type_recognition_lite.h"
#include "src/tiff_parser.h"

namespace piex {
namespace {

using binary_parse::RangeCheckedBytePtr;
using image_type_recognition::RawImageTypes;
using image_type_recognition::RecognizeRawImageTypeLite;
using tiff_directory::Endian;
using tiff_directory::TiffDirectory;

Error GetPreviewData(const TagSet& extended_tags,
                     const std::uint32_t tiff_offset,
                     const std::uint32_t number_of_ifds,
                     StreamInterface* stream, TiffContent* tiff_content,
                     PreviewImageData* preview_image_data) {
  TagSet desired_tags = {kExifTagColorSpace,   kExifTagDateTimeOriginal,
                         kExifTagExposureTime, kExifTagFnumber,
                         kExifTagFocalLength,  kExifTagGps,
                         kExifTagIsoSpeed,     kTiffTagDateTime,
                         kTiffTagExifIfd,      kTiffTagCfaPatternDim,
                         kTiffTagMake,         kTiffTagModel,
                         kTiffTagOrientation};
  desired_tags.insert(extended_tags.cbegin(), extended_tags.cend());

  TiffParser tiff_parser(stream, tiff_offset);
  Error error = tiff_parser.Parse(desired_tags, number_of_ifds, tiff_content);
  if (error != kOk) {
    return error;
  }
  if (tiff_content->tiff_directory.empty()) {
    // Returns kFail if the stream does not contain any TIFF structure.
    return kFail;
  }
  return tiff_parser.GetPreviewImageData(*tiff_content, preview_image_data);
}

Error GetPreviewData(const TagSet& extended_tags,
                     const std::uint32_t number_of_ifds,
                     StreamInterface* stream,
                     PreviewImageData* preview_image_data) {
  const std::uint32_t kTiffOffset = 0;
  TiffContent tiff_content;
  return GetPreviewData(extended_tags, kTiffOffset, number_of_ifds, stream,
                        &tiff_content, preview_image_data);
}

Error GetExifData(const std::uint32_t exif_offset, StreamInterface* stream,
                  PreviewImageData* preview_image_data) {
  const TagSet kExtendedTags = {kTiffTagJpegByteCount, kTiffTagJpegOffset};
  const std::uint32_t kNumberOfIfds = 2;
  TiffContent tiff_content;
  return GetPreviewData(kExtendedTags, exif_offset, kNumberOfIfds, stream,
                        &tiff_content, preview_image_data);
}

// Reads the jpeg compressed thumbnail information.
void GetThumbnailOffsetAndLength(const TagSet& extended_tags,
                                 StreamInterface* stream,
                                 PreviewImageData* preview_image_data) {
  TagSet desired_tags = {kTiffTagJpegByteCount, kTiffTagJpegOffset};
  desired_tags.insert(extended_tags.cbegin(), extended_tags.cend());

  const std::uint32_t kNumberOfIfds = 2;
  PreviewImageData thumbnail_data;
  if (GetPreviewData(desired_tags, kNumberOfIfds, stream, &thumbnail_data) ==
      kOk) {
    preview_image_data->thumbnail_offset = thumbnail_data.preview_offset;
    preview_image_data->thumbnail_length = thumbnail_data.preview_length;
  }
}

Error GetExifIfd(const Endian endian, StreamInterface* stream,
                 TiffDirectory* exif_ifd) {
  const std::uint32_t kTiffOffset = 0;
  std::uint32_t offset_to_ifd;
  if (!Get32u(stream, sizeof(offset_to_ifd), endian, &offset_to_ifd)) {
    return kFail;
  }

  std::uint32_t next_ifd_offset;
  TiffDirectory tiff_ifd(endian);
  Error error =
      ParseDirectory(kTiffOffset, offset_to_ifd, endian, {kTiffTagExifIfd},
                     stream, &tiff_ifd, &next_ifd_offset);
  if (error != kOk) {
    return error;
  }

  std::uint32_t exif_offset;
  if (!tiff_ifd.Get(kTiffTagExifIfd, &exif_offset)) {
    return kUnsupported;
  }

  return ParseDirectory(kTiffOffset, exif_offset, endian, {kExifTagMakernotes},
                        stream, exif_ifd, &next_ifd_offset);
}

struct Image {
  std::uint16_t width = 0;
  std::uint16_t height = 0;
  std::uint32_t length = 0;
  std::uint32_t offset = 0;

  bool operator>(const Image& rhs) const {
    return width > rhs.width && height > rhs.height;
  }
};

bool IsThumbnail(const Image& image) {
  // According to Tiff/EP a thumbnail has max 256 pixels per dimension.
  // http://standardsproposals.bsigroup.com/Home/getPDF/567
  const std::uint16_t kThumbnailAxis = 256;
  return image.width <= kThumbnailAxis && image.height <= kThumbnailAxis;
}

bool GetImageFromIfd(const TiffDirectory& ifd, StreamInterface* stream,
                     Image* image) {
  std::uint32_t compression;
  std::uint32_t photometric_interpretation;
  if (ifd.Get(kTiffTagPhotometric, &photometric_interpretation) &&
      ifd.Get(kTiffTagCompression, &compression)) {
    if (photometric_interpretation == 6 /* YCbCr */ &&
        (compression == 6 /* JPEG(old) */ || compression == 7 /* JPEG */)) {
      std::vector<std::uint32_t> strip_offsets;
      std::vector<std::uint32_t> byte_counts;
      if (ifd.Get(kTiffTagStripOffsets, &strip_offsets) &&
          ifd.Get(kTiffTagStripByteCounts, &byte_counts) &&
          strip_offsets.size() == 1 && byte_counts.size() == 1) {
        image->length = byte_counts[0];
        image->offset = strip_offsets[0];
        return GetPreviewDimensions(image->offset, stream, &image->width,
                                    &image->height);
      }
    }
  }
  return false;
}

Error GetMakernoteIfd(const TiffDirectory& exif_ifd, const Endian endian,
                      const std::uint32_t skip_offset, StreamInterface* stream,
                      std::uint32_t* makernote_offset,
                      TiffDirectory* makernote_ifd) {
  std::uint32_t makernote_length;
  if (!exif_ifd.GetOffsetAndLength(kExifTagMakernotes,
                                   tiff_directory::TIFF_TYPE_UNDEFINED,
                                   makernote_offset, &makernote_length)) {
    return kUnsupported;
  }

  std::uint32_t next_ifd_offset;
  return ParseDirectory(*makernote_offset, *makernote_offset + skip_offset,
                        endian, {kTiffTagImageWidth, kOlymTagCameraSettings,
                                 kOlymTagRawProcessing, kPentaxTagColorSpace},
                        stream, makernote_ifd, &next_ifd_offset);
}

Error GetCameraSettingsIfd(const TiffDirectory& makernote_ifd,
                           const std::uint32_t makernote_offset,
                           const Endian endian, StreamInterface* stream,
                           TiffDirectory* camera_settings_ifd) {
  std::uint32_t camera_settings_offset;
  std::uint32_t camera_settings_length;
  if (!makernote_ifd.GetOffsetAndLength(
          kOlymTagCameraSettings, tiff_directory::TIFF_IFD,
          &camera_settings_offset, &camera_settings_length)) {
    return kUnsupported;
  }

  std::uint32_t next_ifd_offset;
  if (!Get32u(stream, camera_settings_offset, endian,
              &camera_settings_offset)) {
    return kFail;
  }
  return ParseDirectory(makernote_offset,
                        makernote_offset + camera_settings_offset, endian,
                        {kTiffTagBitsPerSample, kTiffTagImageLength}, stream,
                        camera_settings_ifd, &next_ifd_offset);
}

Error GetRawProcessingIfd(const TagSet& desired_tags,
                          const TiffDirectory& makernote_ifd,
                          const std::uint32_t makernote_offset,
                          const Endian endian, StreamInterface* stream,
                          TiffDirectory* raw_processing_ifd) {
  std::uint32_t raw_processing_offset;
  std::uint32_t raw_processing_length;
  if (!makernote_ifd.GetOffsetAndLength(
          kOlymTagRawProcessing, tiff_directory::TIFF_IFD,
          &raw_processing_offset, &raw_processing_length)) {
    return kUnsupported;
  }

  std::uint32_t next_ifd_offset;
  if (!Get32u(stream, raw_processing_offset, endian, &raw_processing_offset)) {
    return kFail;
  }

  return ParseDirectory(
      makernote_offset, makernote_offset + raw_processing_offset, endian,
      desired_tags, stream, raw_processing_ifd, &next_ifd_offset);
}

// Retrieves the preview image offset and length from the camera settings and
// the 'full_width' and 'full_height' from the raw processing ifd in 'stream'.
// Returns kUnsupported if the camera settings are missing, since it is not able
// to get the preview data.
Error GetOlympusPreviewImage(StreamInterface* stream,
                             PreviewImageData* preview_image_data) {
  Endian endian;
  if (!GetEndianness(0 /* tiff offset */, stream, &endian)) {
    return kFail;
  }

  TiffDirectory exif_ifd(endian);
  Error error = GetExifIfd(endian, stream, &exif_ifd);
  if (error != kOk) {
    return error;
  }

  std::uint32_t makernote_offset;
  TiffDirectory makernote_ifd(endian);
  const std::uint32_t kSkipMakernoteStart = 12;
  error = GetMakernoteIfd(exif_ifd, endian, kSkipMakernoteStart, stream,
                          &makernote_offset, &makernote_ifd);
  if (error != kOk) {
    return error;
  }

  const std::uint32_t kThumbnailTag = 0x0100;
  if (makernote_ifd.Has(kThumbnailTag)) {
    if (!makernote_ifd.GetOffsetAndLength(
            kThumbnailTag, tiff_directory::TIFF_TYPE_UNDEFINED,
            &preview_image_data->thumbnail_offset,
            &preview_image_data->thumbnail_length)) {
      return kFail;
    }
  }

  TiffDirectory camera_settings_ifd(endian);
  error = GetCameraSettingsIfd(makernote_ifd, makernote_offset, endian, stream,
                               &camera_settings_ifd);
  if (error != kOk) {
    return error;
  }

  const std::uint32_t kPreviewOffset = 0x0101;
  const std::uint32_t kPreviewLength = 0x0102;
  if (!camera_settings_ifd.Has(kPreviewOffset) ||
      !camera_settings_ifd.Has(kPreviewLength)) {
    return kUnsupported;
  }

  camera_settings_ifd.Get(kPreviewOffset, &preview_image_data->preview_offset);
  preview_image_data->preview_offset += makernote_offset;
  camera_settings_ifd.Get(kPreviewLength, &preview_image_data->preview_length);

  // Get the crop size from the raw processing ifd.
  TiffDirectory raw_processing_ifd(endian);
  error = GetRawProcessingIfd({kOlymTagAspectFrame}, makernote_ifd,
                              makernote_offset, endian, stream,
                              &raw_processing_ifd);
  if (error != kOk) {
    return error;
  }

  if (raw_processing_ifd.Has(kOlymTagAspectFrame)) {
    std::vector<std::uint32_t> aspect_frame(4);
    if (raw_processing_ifd.Get(kOlymTagAspectFrame, &aspect_frame) &&
        aspect_frame[2] > aspect_frame[0] &&
        aspect_frame[3] > aspect_frame[1]) {
      preview_image_data->full_width = aspect_frame[2] - aspect_frame[0] + 1;
      preview_image_data->full_height = aspect_frame[3] - aspect_frame[1] + 1;
      if (preview_image_data->full_width < preview_image_data->full_height) {
        std::swap(preview_image_data->full_width,
                  preview_image_data->full_height);
      }
    }
  }

  return kOk;
}

Error PefGetColorSpace(StreamInterface* stream,
                       PreviewImageData* preview_image_data) {
  Endian endian;
  if (!GetEndianness(0 /* tiff offset */, stream, &endian)) {
    return kFail;
  }

  TiffDirectory exif_ifd(endian);
  Error error = GetExifIfd(endian, stream, &exif_ifd);
  if (error != kOk) {
    return error;
  }

  std::uint32_t makernote_offset;
  TiffDirectory makernote_ifd(endian);
  const std::uint32_t kSkipMakernoteStart = 6;
  error = GetMakernoteIfd(exif_ifd, endian, kSkipMakernoteStart, stream,
                          &makernote_offset, &makernote_ifd);
  if (error != kOk) {
    return error;
  }
  if (makernote_ifd.Has(kPentaxTagColorSpace)) {
    std::uint32_t color_space;
    if (!makernote_ifd.Get(kPentaxTagColorSpace, &color_space)) {
      return kFail;
    }
    preview_image_data->color_space = color_space == 0
                                          ? PreviewImageData::kSrgb
                                          : PreviewImageData::kAdobeRgb;
  }
  return kOk;
}

// Parses the Fuji Cfa header for the image width and height.
bool RafGetDimension(StreamInterface* stream, std::uint32_t* width,
                     std::uint32_t* height) {
  const Endian endian = tiff_directory::kBigEndian;
  std::uint32_t cfa_header_index = 0;  // actual position in the cfa header.
  std::uint32_t cfa_header_entries = 0;
  if (!Get32u(stream, 92 /* cfa header offset */, endian, &cfa_header_index) ||
      !Get32u(stream, cfa_header_index, endian, &cfa_header_entries)) {
    return false;
  }

  // Add 4 to point to the actual read position in the cfa header.
  cfa_header_index += 4;

  for (std::uint32_t i = 0; i < cfa_header_entries; ++i) {
    std::uint16_t id = 0;
    std::uint16_t length = 0;
    if (!Get16u(stream, cfa_header_index, endian, &id) ||
        !Get16u(stream, cfa_header_index + 2, endian, &length)) {
      return false;
    }

    std::uint16_t tmp_width = 0;
    std::uint16_t tmp_height = 0;
    if (id == 0x0111 /* tags the crop dimensions */ &&
        Get16u(stream, cfa_header_index + 4, endian, &tmp_height) &&
        Get16u(stream, cfa_header_index + 6, endian, &tmp_width)) {
      *width = tmp_width;
      *height = tmp_height;
      return true;
    }
    cfa_header_index += 4u + length;
  }
  return false;
}

Error ArwGetPreviewData(StreamInterface* stream,
                        PreviewImageData* preview_image_data) {
  const TagSet extended_tags = {kExifTagHeight, kExifTagWidth,
                                kTiffTagJpegByteCount, kTiffTagJpegOffset,
                                kTiffTagSubIfd};

  GetThumbnailOffsetAndLength(TagSet(), stream, preview_image_data);

  const std::uint32_t kNumberOfIfds = 1;
  return GetPreviewData(extended_tags, kNumberOfIfds, stream,
                        preview_image_data);
}

Error Cr2GetPreviewData(StreamInterface* stream,
                        PreviewImageData* preview_image_data) {
  const TagSet extended_tags = {kExifTagHeight, kExifTagWidth,
                                kTiffTagStripByteCounts, kTiffTagStripOffsets};

  GetThumbnailOffsetAndLength(TagSet(), stream, preview_image_data);

  const std::uint32_t kNumberOfIfds = 1;
  return GetPreviewData(extended_tags, kNumberOfIfds, stream,
                        preview_image_data);
}

Error DngGetPreviewData(StreamInterface* stream,
                        PreviewImageData* preview_image_data) {
  const TagSet extended_tags = {
      kExifTagDefaultCropSize, kTiffTagCompression,  kTiffTagPhotometric,
      kTiffTagStripByteCounts, kTiffTagStripOffsets, kTiffTagSubIfd};

  TiffContent tiff_content;
  const std::uint32_t kNumberOfIfds = 4;
  Error error = GetPreviewData(extended_tags, 0, kNumberOfIfds, stream,
                               &tiff_content, preview_image_data);
  if (error != kOk) {
    return error;
  }

  // Find the jpeg compressed thumbnail and preview image.
  Image preview;
  Image thumbnail;

  // Search for images in IFD0
  Image temp_image;
  if (GetImageFromIfd(tiff_content.tiff_directory[0], stream, &temp_image)) {
    if (IsThumbnail(temp_image)) {
      thumbnail = temp_image;
    } else {
      preview = temp_image;
    }
  }

  // Search for images in other IFDs
  for (const auto& ifd : tiff_content.tiff_directory[0].GetSubDirectories()) {
    if (GetImageFromIfd(ifd, stream, &temp_image)) {
      // Try to find the largest thumbnail/preview.
      if (IsThumbnail(temp_image)) {
        if (temp_image > thumbnail) {
          thumbnail = temp_image;
        }
      } else {
        if (temp_image > preview) {
          preview = temp_image;
        }
      }
    }
  }

  preview_image_data->preview_length = preview.length;
  preview_image_data->preview_offset = preview.offset;
  preview_image_data->thumbnail_length = thumbnail.length;
  preview_image_data->thumbnail_offset = thumbnail.offset;

  return kOk;
}

Error NefGetPreviewData(StreamInterface* stream,
                        PreviewImageData* preview_image_data) {
  const TagSet extended_tags = {kTiffTagImageWidth, kTiffTagImageLength,
                                kTiffTagJpegByteCount, kTiffTagJpegOffset,
                                kTiffTagSubIfd};
  const std::uint32_t kNumberOfIfds = 2;
  Error error =
      GetPreviewData(extended_tags, kNumberOfIfds, stream, preview_image_data);
  if (error != kOk) {
    return error;
  }

  PreviewImageData thumbnail_data;
  GetThumbnailOffsetAndLength(TagSet(), stream, &thumbnail_data);
  preview_image_data->thumbnail_offset = thumbnail_data.thumbnail_offset;
  preview_image_data->thumbnail_length = thumbnail_data.thumbnail_length;

  // The Nikon RAW data provides the dimensions of the sensor image, which are
  // slightly larger than the dimensions of the preview image. In order to
  // determine the correct full width and height of the image, the preview image
  // size needs to be taken into account. Based on experiments the preview image
  // dimensions must be at least 90% of the sensor image dimensions to let it be
  // a full size preview image.
  if (preview_image_data->preview_length > 0) {  // when preview image exists
    const float kEpsilon = 0.9f;

    std::uint16_t width;
    std::uint16_t height;
    if (!GetPreviewDimensions(preview_image_data->preview_offset, stream,
                              &width, &height) ||
        preview_image_data->full_width == 0 ||
        preview_image_data->full_height == 0) {
      return kUnsupported;
    }

    if (static_cast<float>(width) /
                static_cast<float>(preview_image_data->full_width) >
            kEpsilon ||
        static_cast<float>(height) /
                static_cast<float>(preview_image_data->full_height) >
            kEpsilon) {
      preview_image_data->full_width = width;
      preview_image_data->full_height = height;
    }
  }
  return kOk;
}

Error OrfGetPreviewData(StreamInterface* stream,
                        PreviewImageData* preview_image_data) {
  // Omit kUnsupported, because the exif data does not contain any preview
  // image.
  if (GetExifData(0, stream, preview_image_data) == kFail) {
    return kFail;
  }

  return GetOlympusPreviewImage(stream, preview_image_data);
}

Error PefGetPreviewData(StreamInterface* stream,
                        PreviewImageData* preview_image_data) {
  const TagSet extended_tags = {kTiffTagImageWidth, kTiffTagImageLength,
                                kTiffTagJpegByteCount, kTiffTagJpegOffset,
                                kTiffTagSubIfd};
  const std::uint32_t kNumberOfIfds = 3;
  Error error =
      GetPreviewData(extended_tags, kNumberOfIfds, stream, preview_image_data);
  if (error != kOk) {
    return error;
  }

  error = PefGetColorSpace(stream, preview_image_data);
  if (error != kOk) {
    return error;
  }

  PreviewImageData thumbnail_data;
  GetThumbnailOffsetAndLength(TagSet(), stream, &thumbnail_data);
  preview_image_data->thumbnail_offset = thumbnail_data.thumbnail_offset;
  preview_image_data->thumbnail_length = thumbnail_data.thumbnail_length;

  return kOk;
}

Error RafGetPreviewData(StreamInterface* stream,
                        PreviewImageData* preview_image_data) {
  // Parse the Fuji RAW header to get the offset and length of the preview
  // image, which contains the Exif information.
  const Endian endian = tiff_directory::kBigEndian;
  std::uint32_t preview_offset = 0;
  std::uint32_t preview_length = 0;
  if (!Get32u(stream, 84 /* preview offset */, endian, &preview_offset) ||
      !Get32u(stream, 88 /* preview length */, endian, &preview_length)) {
    return kFail;
  }

  if (!RafGetDimension(stream, &preview_image_data->full_width,
                       &preview_image_data->full_height)) {
    return kFail;
  }

  if (preview_length > 0) {  // when preview image exists
    // Parse the Exif information from the preview image. Omit kUnsupported,
    // because the exif data does not contain any preview image.
    const std::uint32_t exif_offset = preview_offset + 12;
    if (GetExifData(exif_offset, stream, preview_image_data) == kFail) {
      return kFail;
    }
  }

  // Merge the Exif data with the RAW data to form the preview_image_data.
  // The preview offset and length extracted from the Exif data are actually
  // the thumbnail offset and length.
  preview_image_data->thumbnail_offset = preview_image_data->preview_offset;
  preview_image_data->thumbnail_offset += 160;  // Skip the cfa header.
  preview_image_data->thumbnail_length = preview_image_data->preview_length;
  preview_image_data->preview_offset = preview_offset;
  preview_image_data->preview_length = preview_length;
  return kOk;
}

Error Rw2GetPreviewData(StreamInterface* stream,
                        PreviewImageData* preview_image_data) {
  const TagSet extended_tags = {kPanaTagTopBorder,     kPanaTagLeftBorder,
                                kPanaTagBottomBorder,  kPanaTagRightBorder,
                                kPanaTagIso,           kPanaTagJpegImage,
                                kTiffTagJpegByteCount, kTiffTagJpegOffset};
  // Parse the RAW data to get the ISO, offset and length of the preview image,
  // which contains the Exif information.
  const std::uint32_t kNumberOfIfds = 1;
  PreviewImageData preview_data;
  Error error =
      GetPreviewData(extended_tags, kNumberOfIfds, stream, &preview_data);
  if (error != kOk) {
    return error;
  }

  if (preview_data.preview_length > 0) {  // when preview image exists
    // Parse the Exif information from the preview image. Omit kUnsupported,
    // because the exif data does not contain any preview image.
    const std::uint32_t exif_offset = preview_data.preview_offset + 12;
    if (GetExifData(exif_offset, stream, preview_image_data) == kFail) {
      return kFail;
    }
    // The preview offset and length extracted from the Exif data are actually
    // the thumbnail offset and length.
    preview_image_data->thumbnail_offset =
        exif_offset + preview_image_data->preview_offset;
    preview_image_data->thumbnail_length = preview_image_data->preview_length;
  }

  // Merge the Exif data with the RAW data to form the preview_image_data.
  preview_image_data->preview_offset = preview_data.preview_offset;
  preview_image_data->preview_length = preview_data.preview_length;
  preview_image_data->iso = preview_data.iso;
  preview_image_data->full_width = preview_data.full_width;
  preview_image_data->full_height = preview_data.full_height;

  return kOk;
}

Error SrwGetPreviewData(StreamInterface* stream,
                        PreviewImageData* preview_image_data) {
  GetThumbnailOffsetAndLength({kTiffTagSubIfd}, stream, preview_image_data);

  const TagSet extended_tags = {kExifTagWidth, kExifTagHeight,
                                kTiffTagJpegByteCount, kTiffTagJpegOffset,
                                kTiffTagSubIfd};
  const std::uint32_t kNumberOfIfds = 1;
  return GetPreviewData(extended_tags, kNumberOfIfds, stream,
                        preview_image_data);
}

}  // namespace

size_t BytesRequiredForIsRaw() {
  return image_type_recognition::GetNumberOfBytesForIsRawLite();
}

bool IsRaw(StreamInterface* data) {
  const size_t bytes = BytesRequiredForIsRaw();
  if (data == nullptr) {
    return false;
  }

  // Read required number of bytes into a vector.
  std::vector<std::uint8_t> file_header(bytes);
  if (data->GetData(0, file_header.size(), file_header.data()) != kOk) {
    return false;
  }

  RangeCheckedBytePtr data_buffer(file_header.data(), file_header.size());

  return image_type_recognition::IsRawLite(data_buffer);
}

Error GetPreviewImageData(StreamInterface* data,
                          PreviewImageData* preview_image_data) {
  const size_t bytes = BytesRequiredForIsRaw();
  if (data == nullptr || bytes == 0) {
    return kFail;
  }

  std::vector<std::uint8_t> file_header(bytes);
  Error error = data->GetData(0, file_header.size(), file_header.data());
  if (error != kOk) {
    return error;
  }
  RangeCheckedBytePtr header_buffer(file_header.data(), file_header.size());

  switch (RecognizeRawImageTypeLite(header_buffer)) {
    case image_type_recognition::kArwImage:
      return ArwGetPreviewData(data, preview_image_data);
    case image_type_recognition::kCr2Image:
      return Cr2GetPreviewData(data, preview_image_data);
    case image_type_recognition::kDngImage:
      return DngGetPreviewData(data, preview_image_data);
    case image_type_recognition::kNefImage:
    case image_type_recognition::kNrwImage:
      return NefGetPreviewData(data, preview_image_data);
    case image_type_recognition::kOrfImage:
      return OrfGetPreviewData(data, preview_image_data);
    case image_type_recognition::kPefImage:
      return PefGetPreviewData(data, preview_image_data);
    case image_type_recognition::kRafImage:
      return RafGetPreviewData(data, preview_image_data);
    case image_type_recognition::kRw2Image:
      return Rw2GetPreviewData(data, preview_image_data);
    case image_type_recognition::kSrwImage:
      return SrwGetPreviewData(data, preview_image_data);
    default:
      return kUnsupported;
  }
}

std::vector<std::string> SupportedExtensions() {
  return {"ARW", "CR2", "DNG", "NEF", "NRW", "ORF", "PEF", "RAF", "RW2", "SRW"};
}

}  // namespace piex
