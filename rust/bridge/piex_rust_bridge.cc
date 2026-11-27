// Copyright 2026 Google LLC
//
// Drop-in C++ replacement implementation of PIEX public API
// wrapping pure Rust implementation via CXX bridge.

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "rust/include/piex_rust.h"
#include "rust/src/cxx_bridge.rs.h"
#include "src/image_type_recognition/image_type_recognition_lite.h"
#include "src/piex_types.h"

namespace piex_rust {

using ::piex::Error;
using ::piex::PreviewImageData;
using ::piex::StreamInterface;
using ::piex::image_type_recognition::RawImageTypes;
namespace {

namespace bridge = ::piex_bridge;

std::string ToString(const rust::Vec<uint8_t>& bytes) {
  return std::string(reinterpret_cast<const char*>(bytes.data()), bytes.size());
}

rust::Vec<uint8_t> ToBytes(const std::string& str) {
  // `rust::Vec` has no public `resize`/`set_len`, so C++ cannot memcpy into
  // it without writing past its length. Let Rust copy the bytes with a single
  // memcpy (`<[u8]>::to_vec`) instead.
  return bridge::rust_bytes_to_vec(rust::Slice<const uint8_t>(
      reinterpret_cast<const uint8_t*>(str.data()), str.size()));
}

void CopyPreviewData(const bridge::PreviewImageData& src,
                     PreviewImageData* dst) {
  if (!dst) {
    return;
  }
  dst->preview = src.preview;
  dst->thumbnail = src.thumbnail;
  dst->exif_orientation = src.exif_orientation;
  dst->color_space = src.color_space;
  dst->full_width = src.full_width;
  dst->full_height = src.full_height;
  dst->maker = ToString(src.maker);
  dst->model = ToString(src.model);
  dst->date_time = ToString(src.date_time);
  dst->iso = src.iso;
  dst->exposure_time = src.exposure_time;
  dst->fnumber = src.fnumber;
  dst->focal_length = src.focal_length;

  dst->gps.is_valid = src.gps.is_valid;
  dst->gps.latitude_ref = static_cast<char>(src.gps.latitude_ref);
  std::copy(src.gps.latitude.begin(), src.gps.latitude.end(),
            dst->gps.latitude);
  dst->gps.longitude_ref = static_cast<char>(src.gps.longitude_ref);
  std::copy(src.gps.longitude.begin(), src.gps.longitude.end(),
            dst->gps.longitude);
  dst->gps.altitude_ref = src.gps.altitude_ref;
  dst->gps.altitude = src.gps.altitude;
  std::copy(src.gps.time_stamp.begin(), src.gps.time_stamp.end(),
            dst->gps.time_stamp);
  dst->gps.date_stamp = ToString(src.gps.date_stamp);

  // Like C++ PIEX, write the two CFA dimensions in place.
  if (dst->cfa_pattern_dim.size() >= 2 || src.cfa_pattern_dim[0] != 0 ||
      src.cfa_pattern_dim[1] != 0) {
    dst->cfa_pattern_dim[0] = src.cfa_pattern_dim[0];
    dst->cfa_pattern_dim[1] = src.cfa_pattern_dim[1];
  }
}

// Converts caller state into bridge state so that the Rust implementation
// updates it in place, matching C++ PIEX semantics (fields it does not parse
// stay untouched, partial writes on failure are preserved).
bridge::PreviewImageData ToBridge(const PreviewImageData& src) {
  bridge::PreviewImageData dst{};
  dst.preview = src.preview;
  dst.thumbnail = src.thumbnail;
  dst.exif_orientation = src.exif_orientation;
  dst.color_space = src.color_space;
  dst.full_width = src.full_width;
  dst.full_height = src.full_height;
  dst.maker = ToBytes(src.maker);
  dst.model = ToBytes(src.model);
  dst.date_time = ToBytes(src.date_time);
  dst.iso = src.iso;
  dst.exposure_time = src.exposure_time;
  dst.fnumber = src.fnumber;
  dst.focal_length = src.focal_length;
  dst.gps.is_valid = src.gps.is_valid;
  dst.gps.latitude_ref = static_cast<uint8_t>(src.gps.latitude_ref);
  std::copy(std::begin(src.gps.latitude), std::end(src.gps.latitude),
            dst.gps.latitude.begin());
  dst.gps.longitude_ref = static_cast<uint8_t>(src.gps.longitude_ref);
  std::copy(std::begin(src.gps.longitude), std::end(src.gps.longitude),
            dst.gps.longitude.begin());
  dst.gps.altitude_ref = src.gps.altitude_ref;
  dst.gps.altitude = src.gps.altitude;
  std::copy(std::begin(src.gps.time_stamp), std::end(src.gps.time_stamp),
            dst.gps.time_stamp.begin());
  dst.gps.date_stamp = ToBytes(src.gps.date_stamp);
  if (src.cfa_pattern_dim.size() >= 2) {
    dst.cfa_pattern_dim[0] = src.cfa_pattern_dim[0];
    dst.cfa_pattern_dim[1] = src.cfa_pattern_dim[1];
  }
  return dst;
}

Error MapErrorCode(int code) {
  switch (code) {
    case 0:
      return Error::kOk;
    case 2:
      return Error::kUnsupported;
    default:
      return Error::kFail;
  }
}

}  // namespace

size_t BytesRequiredForIsRaw() {
  return bridge::rust_bytes_required_for_is_raw();
}

bool IsRaw(StreamInterface* data) {
  if (!data) {
    return false;
  }
  return bridge::rust_is_raw(*data);
}

Error GetPreviewImageData(StreamInterface* data,
                          PreviewImageData* preview_image_data,
                          RawImageTypes* output_type) {
  if (!data || !preview_image_data) {
    return Error::kFail;
  }

  // The Rust side only writes the fields (and type) that C++ PIEX would write.
  // `*output_type` may be uninitialized, so it is never read: a sentinel tells
  // whether the Rust side determined a type.
  bridge::PreviewImageData bridge_data = ToBridge(*preview_image_data);
  constexpr uint32_t kTypeNotWritten = UINT32_MAX;
  uint32_t type_val = kTypeNotWritten;

  int status =
      bridge::rust_get_preview_image_data(*data, bridge_data, type_val);
  if (output_type && type_val != kTypeNotWritten) {
    *output_type = static_cast<RawImageTypes>(type_val);
  }
  CopyPreviewData(bridge_data, preview_image_data);
  return MapErrorCode(status);
}

bool GetDngInformation(StreamInterface* data, std::uint32_t* width,
                       std::uint32_t* height,
                       std::vector<std::uint32_t>* cfa_pattern_dim) {
  if (!data || !width || !height || !cfa_pattern_dim) {
    return false;
  }
  std::array<uint32_t, 2> cfa{0, 0};
  if (!bridge::rust_get_dng_information(*data, *width, *height, cfa)) {
    return false;
  }
  *cfa_pattern_dim = {cfa[0], cfa[1]};
  return true;
}

bool GetOrientation(StreamInterface* data, std::uint32_t* orientation) {
  if (!data || !orientation) {
    return false;
  }
  return bridge::rust_get_orientation(*data, *orientation);
}

std::vector<std::string> SupportedExtensions() {
  return {"ARW", "CR2", "CR3", "DNG", "NEF", "NRW",
          "ORF", "PEF", "RAF", "RW2", "SRW"};
}

// CR3 APIs
bool Cr3GetOrientation(StreamInterface* stream, std::uint32_t* orientation) {
  if (!stream || !orientation) {
    return false;
  }
  return bridge::rust_cr3_get_orientation(*stream, *orientation);
}

Error Cr3GetPreviewData(StreamInterface* stream,
                        PreviewImageData* preview_image_data) {
  if (!stream || !preview_image_data) {
    return Error::kFail;
  }
  bridge::PreviewImageData bridge_data = ToBridge(*preview_image_data);
  int status = bridge::rust_cr3_get_preview_data(*stream, bridge_data);
  CopyPreviewData(bridge_data, preview_image_data);
  return MapErrorCode(status);
}

}  // namespace piex_rust
