// Copyright 2020 Google Inc.
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

#include "src/piex_cr3.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <unordered_set>

#include "src/binary_parse/range_checked_byte_ptr.h"
#include "src/piex_types.h"
#include "src/tiff_directory/tiff_directory.h"
#include "src/tiff_parser.h"

namespace piex {
namespace {

constexpr size_t kUuidSize = 16;
using Uuid = std::array<std::uint8_t, kUuidSize>;
// Uuid of uuid box under the moov box.
constexpr Uuid kUuidMoov = {0x85, 0xc0, 0xb6, 0x87, 0x82, 0x0f, 0x11, 0xe0,
                            0x81, 0x11, 0xf4, 0xce, 0x46, 0x2b, 0x6a, 0x48};

// Uuid of uuid box containing PRVW box.
constexpr Uuid kUuidPrvw = {0xea, 0xf4, 0x2b, 0x5e, 0x1c, 0x98, 0x4b, 0x88,
                            0xb9, 0xfb, 0xb7, 0xdc, 0x40, 0x6e, 0x4d, 0x16};

constexpr size_t kTagSize = 4;
using BoxTag = std::array<std::uint8_t, kTagSize>;

constexpr BoxTag NewTag(const char s[kTagSize + 1]) {
  return BoxTag{s[0], s[1], s[2], s[3]};
}

constexpr BoxTag kUuidTag = NewTag("uuid");
constexpr BoxTag kPrvwTag = NewTag("PRVW");
constexpr BoxTag kThmbTag = NewTag("THMB");
constexpr BoxTag kCmt1Tag = NewTag("CMT1");
constexpr BoxTag kCmt2Tag = NewTag("CMT2");
constexpr BoxTag kStblTag = NewTag("stbl");
constexpr BoxTag kStsdTag = NewTag("stsd");
constexpr BoxTag kCrawTag = NewTag("CRAW");
constexpr BoxTag kStszTag = NewTag("stsz");
constexpr BoxTag kCo64Tag = NewTag("co64");
constexpr BoxTag kMdatTag = NewTag("mdat");

// Convenience class for a box.
class Box {
 public:
  Box()
      : is_valid_(false),
        tag_(BoxTag()),
        offset_(0),
        header_offset_(0),
        next_box_offset_(0) {}
  Box(const BoxTag& tag, size_t offset, size_t header_length, size_t length)
      : is_valid_(true),
        tag_(tag),
        offset_(offset),
        header_offset_(offset + header_length),
        next_box_offset_(offset + length) {}

  bool IsValid() const { return is_valid_ && next_box_offset_ > offset_; }
  const BoxTag& tag() const { return tag_; }

  // Returns offset from start of file.
  size_t offset() const { return offset_; }
  // Returns offset from start of file, including box's header.
  size_t header_offset() const { return header_offset_; }
  // Returns offset from start of file of the next box, accounting for size of
  // this box.
  size_t next_box_offset() const { return next_box_offset_; }

 private:
  bool is_valid_;
  BoxTag tag_;
  size_t offset_;
  size_t header_offset_;
  size_t next_box_offset_;
};

struct ProcessData {
  PreviewImageData* preview_image_data = nullptr;
  Image mdat_image;
  Image prvw_image;
};

// Wraps Get16u w/ assumption that CR3 is always big endian, based on
// ISO/IEC 14496-12 specification that all box fields are big endian.
bool Get16u(StreamInterface* stream, size_t offset, std::uint16_t* value) {
  return Get16u(stream, offset, tiff_directory::kBigEndian, value);
}

// Wraps Get32u w/ assumption that CR3 is always big endian, based on
// ISO/IEC 14496-12 specification that all box fields are big endian.
bool Get32u(StreamInterface* stream, size_t offset, std::uint32_t* value) {
  return Get32u(stream, offset, tiff_directory::kBigEndian, value);
}

// Always big endian, based on ISO/IEC 14496-12 specification that all box
// fields are big endian.
bool Get64u(StreamInterface* stream, size_t offset, std::uint64_t* value) {
  std::uint8_t data[8];
  if (stream->GetData(offset, 8, data) == kOk) {
    *value = (data[0] * 0x1000000u) | (data[1] * 0x10000u) |
             (data[2] * 0x100u) | data[3];
    *value <<= 32;
    *value = (data[4] * 0x1000000u) | (data[5] * 0x10000u) |
             (data[6] * 0x100u) | data[7];
    return true;
  } else {
    return false;
  }
}

// Jpeg box offsets based on the box tag. The expected layout is as follows:
//        Byte Offset Type     Meaning
//                  0 [long]   size of box
//                  4 [char[]] box tag
//       offset.width [short]  width of jpeg
//      offset.height [short]  height of jpeg
//   offset.jpeg_size [long]   number of bytes in jpeg
//   offset.jpeg_data [byte[]] start of jpeg data
struct JpegBoxOffset {
  size_t width = 0;
  size_t height = 0;
  size_t jpeg_size = 0;
  size_t jpeg_data = 0;
};

// Processes box w/ JPEG data. Box must be PRVW and THMB boxes.
bool ProcessJpegBox(StreamInterface* stream, const Box& box, Image* image) {
  static constexpr JpegBoxOffset kPrvwJpegOffsets{14, 16, 20, 24};
  static constexpr JpegBoxOffset kThmbJpegOffsets{12, 14, 16, 24};
  if (box.tag() != kPrvwTag && box.tag() != kThmbTag) {
    return false;
  }
  const JpegBoxOffset& offsets =
      box.tag() == kPrvwTag ? kPrvwJpegOffsets : kThmbJpegOffsets;
  uint16_t width, height;
  uint32_t jpeg_size;
  if (!Get16u(stream, box.offset() + offsets.width, &width)) {
    return false;
  }
  if (!Get16u(stream, box.offset() + offsets.height, &height)) {
    return false;
  }
  if (!Get32u(stream, box.offset() + offsets.jpeg_size, &jpeg_size)) {
    return false;
  }
  image->format = Image::kJpegCompressed;
  image->width = width;
  image->height = height;
  image->offset = box.offset() + offsets.jpeg_data;
  image->length = jpeg_size;
  return true;
}

// Parses the Exif IFD0 tags at tiff_offset.
bool ParseExifIfd0(StreamInterface* stream, size_t tiff_offset,
                   PreviewImageData* preview_image_data) {
  static const TagSet kIfd0TagSet = {kTiffTagModel, kTiffTagMake,
                                     kTiffTagOrientation, kTiffTagImageWidth,
                                     kTiffTagImageLength};
  TiffContent content;
  TiffParser(stream, tiff_offset).Parse(kIfd0TagSet, 1, &content);
  if (content.tiff_directory.size() != 1) {
    return false;
  }

  content.tiff_directory[0].Get(kTiffTagModel, &preview_image_data->model);
  content.tiff_directory[0].Get(kTiffTagMake, &preview_image_data->maker);
  content.tiff_directory[0].Get(kTiffTagOrientation,
                                &preview_image_data->exif_orientation);
  content.tiff_directory[0].Get(kTiffTagImageWidth,
                                &preview_image_data->full_width);
  content.tiff_directory[0].Get(kTiffTagImageLength,
                                &preview_image_data->full_height);
  return true;
}

// Parses the Exif Exif IFD tags at tiff_offset.
bool ParseExifExifIfd(StreamInterface* stream, size_t tiff_offset,
                      PreviewImageData* preview_image_data) {
  static const TagSet kExifIfdTagSet = {kExifTagDateTimeOriginal,
                                        kExifTagExposureTime, kExifTagFnumber,
                                        kExifTagFocalLength, kExifTagIsoSpeed};
  TiffContent content;
  TiffParser(stream, tiff_offset).Parse(kExifIfdTagSet, 1, &content);
  if (content.tiff_directory.size() != 1) {
    return false;
  }

  content.tiff_directory[0].Get(kExifTagDateTimeOriginal,
                                &preview_image_data->date_time);
  GetRational(kExifTagExposureTime, content.tiff_directory[0], 1,
              &preview_image_data->exposure_time);
  GetRational(kExifTagFnumber, content.tiff_directory[0], 1,
              &preview_image_data->fnumber);
  GetRational(kExifTagFocalLength, content.tiff_directory[0], 1,
              &preview_image_data->focal_length);
  content.tiff_directory[0].Get(kExifTagIsoSpeed, &preview_image_data->iso);
  return true;
}

// Returns the next box or an invalid box.
//
// Based on ISO/IEC 14496-12: boxes start with a header: size and type. The size
// can be compact (32-bits) or extended (64-bit, e.g. mdat box).
// The type can be compact (32 bits) or extended (full UUID, e.g. uuid boxes).
// values are stored after the compact size/type.
//
// Fields in a box are big-endian.
Box GetNextBox(StreamInterface* stream, size_t offset) {
  uint32_t length_32;
  if (!Get32u(stream, offset, &length_32)) {
    return Box();
  }
  BoxTag tag;
  Error status =
      stream->GetData(offset + sizeof(length_32), kTagSize, tag.data());
  if (status != kOk) {
    return Box();
  }
  size_t length;
  size_t header_offset = sizeof(length_32) + sizeof(tag);
  if (length_32 == 1) {
    // Magic number of 1 implies extended size.
    uint64_t length_64 = 0;
    if (!Get64u(stream, offset + header_offset, &length_64)) {
      return Box();
    }
    length = length_64;
    header_offset += sizeof(length_64);
  } else {
    // Compact size.
    length = length_32;
  }
  return Box(tag, offset, header_offset, length);
}

// Searches for the next box with the given tag.
Box GetNextBoxWithTag(StreamInterface* stream, size_t offset,
                      const BoxTag& expected_tag) {
  while (true) {
    Box box = GetNextBox(stream, offset);
    if (!box.IsValid() || box.tag() == expected_tag) {
      return box;
    }
    offset = box.next_box_offset();
  }
}

// Returns the width, height, and content type from the CRAW box.
bool ProcessCrawBox(StreamInterface* stream, const Box& craw_box,
                    uint16_t* width, uint16_t* height, uint16_t* content_type) {
  constexpr size_t kWidthOffset = 32;
  if (!Get16u(stream, craw_box.offset() + kWidthOffset, width)) {
    return false;
  }

  constexpr size_t kHeightOffset = 34;
  if (!Get16u(stream, craw_box.offset() + kHeightOffset, height)) {
    return false;
  }

  constexpr size_t kTypeOffset = 86;
  if (!Get16u(stream, craw_box.offset() + kTypeOffset, content_type)) {
    return false;
  }
  return true;
}

// stsz box offset:
//        Byte Offset Type     Meaning
//                  0 [long]   size of box
//                  4 [char[]] box tag
//                  8 [long]   version/flags
//                 12 [long]   sample size
//                 16 [long]   number of entries in sample table
//                 20 [long[]] sample table if samples size is 0
bool ProcessStszBox(StreamInterface* stream, const Box& stsz_box,
                    uint32_t* image_size) {
  uint32_t sample_size;
  if (!Get32u(stream, stsz_box.offset() + 12, &sample_size)) {
    return false;
  }
  if (sample_size > 0) {
    *image_size = sample_size;
    return true;
  }
  // sample_size of 0 implies the data is in the sample table. We expect only
  // one entry. This is true of Canon EOS RP Cr3 files.
  uint32_t count;
  if (!Get32u(stream, stsz_box.offset() + 16, &count)) {
    return false;
  }
  if (count != 1) {
    // Expect at most one entry in the table.
    return false;
  }
  return Get32u(stream, stsz_box.offset() + 20, image_size);
}

// co64 box offsets:
//        Byte Offset Type     Meaning
//                  0 [long]   size of box
//                  4 [char[]] box tag
//                  8 [long]   version
//                 12 [long]   count (expect to be value 1)
//                 16 [long]   offset of image data in mdat
bool ProcessCo64(StreamInterface* stream, const Box& co64_box,
                 uint32_t* image_offset) {
  uint32_t count = 0;
  if (!Get32u(stream, co64_box.header_offset() + 4, &count)) {
    return false;
  }
  if (count != 1) {
    return false;
  }
  return Get32u(stream, co64_box.header_offset() + 8, image_offset);
}

// Process the stbl box. Expected box layout:
// stbl
//   stsd
//     CRAW  (embedded image (JPEG) information)
//   (0 or more skipped boxes)
//   stsz (embedded image byte size)
//   (0 or more skipped boxes)
//   co64 (offset of embedded image, relative to mdat box)
bool ProcessStblBox(StreamInterface* stream, const Box& stbl_box,
                    ProcessData* data) {
  Box stsd_box = GetNextBoxWithTag(stream, stbl_box.header_offset(), kStsdTag);
  if (!stsd_box.IsValid()) {
    return false;
  }
  // This is either CRAW or CTMD. Skip when CTMD.
  Box craw_box = GetNextBox(stream, stsd_box.header_offset() + 8);
  if (!craw_box.IsValid()) {
    return false;
  }
  if (craw_box.tag() != kCrawTag) {
    return true;
  }
  // CRAW contains info about the full-size image embedded in the mdat box.
  // The image is either JPEG or HEVC.
  uint16_t image_width = 0;
  uint16_t image_height = 0;
  uint16_t content_type = 0;
  if (!ProcessCrawBox(stream, craw_box, &image_width, &image_height,
                      &content_type)) {
    return false;
  }
  // Only continue if JPEG or HEVC content.
  constexpr uint16_t kJpegContentType = 3;
  constexpr uint16_t kHevcContentType = 4;
  if (content_type != kJpegContentType && content_type != kHevcContentType) {
    return true;
  }

  // Skip until we find stsz, contains the size (# of bytes) of image data.
  Box stsz_box =
      GetNextBoxWithTag(stream, stsd_box.next_box_offset(), kStszTag);
  if (!stsz_box.IsValid()) {
    return false;
  }
  uint32_t image_size;
  if (!ProcessStszBox(stream, stsz_box, &image_size)) {
    return false;
  }

  // Skip until we find co64, contains the offset of image data.
  Box co64_box =
      GetNextBoxWithTag(stream, stsz_box.next_box_offset(), kCo64Tag);
  if (!co64_box.IsValid()) {
    return false;
  }

  uint32_t image_offset = 0;
  if (!ProcessCo64(stream, co64_box, &image_offset)) {
    return false;
  }

  data->mdat_image.format = content_type == kJpegContentType
                                ? Image::kJpegCompressed
                                : Image::kHevcCompressed;
  data->mdat_image.width = image_width;
  data->mdat_image.height = image_height;
  data->mdat_image.length = image_size;
  // This offset is relative to the position of the mdat box. The value will
  // be updated once mdat's offset is found.
  data->mdat_image.offset = image_offset;
  return true;
}

// Returns true if we should parse the children of the box.
bool DoProcessChildren(const BoxTag& tag) {
  static const std::set<BoxTag> kTags = {NewTag("trak"), NewTag("moov"),
                                         NewTag("mdia"), NewTag("minf")};
  return kTags.find(tag) != kTags.end();
}

// Processes box and returns offset of the next box to process.
// A return value of 0 indicates an error.
//
// Outline of hierarchy and important boxes:
// ftyp
// moov
//   uuid (id is kUuidMoov)
//     ... boxes we skip ...
//     CMT1 (EXIF data)
//     CMT2 (EXIF data)
//     ... boxes we skip ...
//     THMB (160x120 JPEG thumbnail, embedded in this box)
//   trak
//     tkhd
//     mdia
//     ... boxes we skip ...
//     minf
//       ... boxes we skip ...
//       stbl
//         stsd
//           CRAW (Full image preview, type (JPEG or HEVC), width, height. The
//                 image data is found in mdat box, below.)
//       ... boxes we skip ...
//       stsz (Size of preview, in bytes)
//       ... boxes we skip ...
//       co64 (Location/offset of full preview data in mdat)
//   .. boxes we skip ...
// uuid (id is kUuidPrvw)
//   PRVW (1620x1080 JPEG preview, embedded in this box)
// mdat
//   Full image preview (JPEG or HEVC)
//   ... RAW image data ...
size_t ProcessBox(StreamInterface* stream, const Box& box, ProcessData* data) {
  // Parse child boxes.
  if (box.tag() == kUuidTag) {
    // Uuid box have extended box types.
    Uuid uuid;
    if (stream->GetData(box.header_offset(), uuid.size(), uuid.data()) != kOk) {
      return 0;
    }
    if (uuid == kUuidPrvw) {
      return box.header_offset() + uuid.size() + 8;
    } else if (uuid == kUuidMoov) {
      return box.header_offset() + uuid.size();
    }  // else skip the box, below.
  } else if (DoProcessChildren(box.tag())) {
    return box.header_offset();
  }

  // Potentially process the data contained in the box.
  bool success;
  if (box.tag() == kMdatTag) {
    // mdat_image.offset is relative to mdat's header, update it to be absolute
    // offset to the image data.
    data->mdat_image.offset += box.header_offset();
    success = true;
  } else if (box.tag() == kStblTag) {
    success = ProcessStblBox(stream, box, data);
  } else if (box.tag() == kPrvwTag) {
    // Preview jpeg. 1620x1080 for EOS R.
    success = ProcessJpegBox(stream, box, &data->prvw_image);
  } else if (box.tag() == kThmbTag) {
    // Thumbnail jpeg. 160x120 for EOS R.
    success = ProcessJpegBox(stream, box, &data->preview_image_data->thumbnail);
  } else if (box.tag() == kCmt1Tag) {
    success =
        ParseExifIfd0(stream, box.header_offset(), data->preview_image_data);
  } else if (box.tag() == kCmt2Tag) {
    success =
        ParseExifExifIfd(stream, box.header_offset(), data->preview_image_data);
  } else {
    // This box isn't interesting, skip it.
    success = true;
  }
  return success ? box.next_box_offset() : 0;
}

bool ProcessStream(StreamInterface* stream, const BoxTag& last_chunk,
                   ProcessData* data) {
  size_t offset = 0;
  while (true) {
    Box box = GetNextBox(stream, offset);
    if (!box.IsValid()) {
      return false;
    }
    size_t new_offset = ProcessBox(stream, box, data);
    if (new_offset <= offset) {
      return false;
    }
    if (box.tag() == last_chunk) {
      return true;
    }
    offset = new_offset;
  }
}

bool IsImage(StreamInterface* stream, const Image& image) {
  if (image.format != Image::kJpegCompressed) {
    // Pass responsibility to the caller.
    return true;
  }
  // Check for JPEG magic number at start. This could be HEVC data.
  constexpr std::array<uint8_t, 3> kJpegMagicNumber = {0xFF, 0xD8, 0xFF};
  std::array<uint8_t, 3> magic_number;
  if (stream->GetData(image.offset, magic_number.size(), magic_number.data()) !=
      kOk) {
    return false;
  }
  return magic_number == kJpegMagicNumber;
}

}  // namespace

Error Cr3GetPreviewData(StreamInterface* stream,
                        PreviewImageData* preview_image_data) {
  ProcessData data{preview_image_data};
  if (!ProcessStream(stream, kMdatTag, &data)) {
    return kFail;
  }
  // Prefer image in mdata box, as spec ensures it is the largest image.
  if (data.mdat_image.length > 0 && IsImage(stream, data.mdat_image)) {
    preview_image_data->preview = data.mdat_image;
  } else if (data.prvw_image.length > 0 && IsImage(stream, data.prvw_image)) {
    preview_image_data->preview = data.prvw_image;
  } else {
    return kFail;
  }
  return kOk;
}

bool Cr3GetOrientation(StreamInterface* stream, std::uint32_t* orientation) {
  PreviewImageData preview_image_data;
  ProcessData data{&preview_image_data};
  if (ProcessStream(stream, kCmt1Tag, &data)) {
    *orientation = preview_image_data.exif_orientation;
    return true;
  }
  return false;
}

}  // namespace piex
