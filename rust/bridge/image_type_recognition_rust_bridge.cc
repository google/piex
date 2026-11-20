// Copyright 2026 Google LLC
//
// Drop-in C++ replacement implementation of image_type_recognition_lite
// wrapping pure Rust implementation via CXX bridge.

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "rust/include/piex_rust.h"
#include "rust/src/cxx_bridge.rs.h"
#include "src/binary_parse/range_checked_byte_ptr.h"
#include "src/image_type_recognition/image_type_recognition_lite.h"

namespace piex_rust {
namespace image_type_recognition {

using ::piex::binary_parse::RangeCheckedBytePtr;
using ::piex::image_type_recognition::kNonRawImage;
using ::piex::image_type_recognition::RawImageTypes;

namespace bridge = ::piex_bridge;

bool IsRaw(const RawImageTypes type) {
  return bridge::rust_is_raw_type(static_cast<uint32_t>(type));
}

bool IsOfType(const RangeCheckedBytePtr& source, const RawImageTypes type) {
  const size_t req = bridge::rust_get_number_of_bytes_for_is_of_type(
      static_cast<uint32_t>(type));
  if (req == 0) {
    return false;
  }
  const RangeCheckedBytePtr limited_source = source.pointerToSubArray(0, req);
  if (limited_source.errorOccurred()) {
    return false;
  }
  const std::vector<uint8_t> buf = limited_source.extractBytes(0, req);
  const rust::Slice<const uint8_t> slice(buf.data(), buf.size());
  return bridge::rust_is_of_type(slice, static_cast<uint32_t>(type));
}

RawImageTypes RecognizeRawImageTypeLite(const RangeCheckedBytePtr& source) {
  // Step through the ascending checker window sizes reported by Rust so each
  // tier narrows `source` via pointerToSubArray(0, tier) identically to C++
  // TypeChecker::LimitSource: a match at a smaller tier never reads or
  // requires pages from larger tiers, and a shorter buffer sets
  // source.errorOccurred() at the first unfulfilled tier. Bytes are appended
  // incrementally into one buffer and only the checkers for `tier` are run.
  std::vector<uint8_t> buf;
  size_t prev_tier = 0;
  for (size_t tier = bridge::rust_next_recognition_window_size(0); tier != 0;
       tier = bridge::rust_next_recognition_window_size(tier)) {
    const RangeCheckedBytePtr limited_source =
        source.pointerToSubArray(0, tier);
    if (limited_source.errorOccurred()) {
      return kNonRawImage;
    }
    if (tier > buf.capacity()) {
      buf.reserve(std::max(tier, buf.capacity() * 2));
    }
    for (size_t i = prev_tier; i < tier; ++i) {
      buf.push_back(limited_source[i]);
    }
    prev_tier = tier;
    const rust::Slice<const uint8_t> slice(buf.data(), buf.size());
    const auto matched = static_cast<RawImageTypes>(
        bridge::rust_recognize_raw_image_type_for_window(slice, tier));
    if (matched != kNonRawImage) {
      return matched;
    }
  }
  return kNonRawImage;
}

size_t GetNumberOfBytesForIsRawLite() {
  return bridge::rust_get_number_of_bytes_for_is_raw_lite();
}

size_t GetNumberOfBytesForIsOfType(const RawImageTypes type) {
  return bridge::rust_get_number_of_bytes_for_is_of_type(
      static_cast<uint32_t>(type));
}

bool IsRawLite(const RangeCheckedBytePtr& source) {
  return bridge::rust_is_raw_type(
      static_cast<uint32_t>(RecognizeRawImageTypeLite(source)));
}

}  // namespace image_type_recognition
}  // namespace piex_rust
