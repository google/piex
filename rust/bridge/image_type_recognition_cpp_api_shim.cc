// Copyright 2026 Google LLC
//
// Drop-in C++ replacement entry points for piex::image_type_recognition,
// forwarding to the Rust implementation in piex_rust::image_type_recognition.

#include <cstddef>

#include "rust/include/piex_rust.h"
#include "src/binary_parse/range_checked_byte_ptr.h"
#include "src/image_type_recognition/image_type_recognition_lite.h"

namespace piex {
namespace image_type_recognition {

bool IsRaw(RawImageTypes type) {
  return ::piex_rust::image_type_recognition::IsRaw(type);
}

bool IsOfType(const binary_parse::RangeCheckedBytePtr& source,
              RawImageTypes type) {
  return ::piex_rust::image_type_recognition::IsOfType(source, type);
}

RawImageTypes RecognizeRawImageTypeLite(
    const binary_parse::RangeCheckedBytePtr& source) {
  return ::piex_rust::image_type_recognition::RecognizeRawImageTypeLite(source);
}

size_t GetNumberOfBytesForIsRawLite() {
  return ::piex_rust::image_type_recognition::GetNumberOfBytesForIsRawLite();
}

size_t GetNumberOfBytesForIsOfType(RawImageTypes type) {
  return ::piex_rust::image_type_recognition::GetNumberOfBytesForIsOfType(type);
}

bool IsRawLite(const binary_parse::RangeCheckedBytePtr& source) {
  return ::piex_rust::image_type_recognition::IsRawLite(source);
}

}  // namespace image_type_recognition
}  // namespace piex
