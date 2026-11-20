// Copyright 2026 Google LLC
//
// Drop-in C++ replacement entry points for namespace piex (piex.h and
// piex_cr3.h), forwarding to the Rust implementation in namespace piex_rust.

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "rust/include/piex_rust.h"
#include "src/image_type_recognition/image_type_recognition_lite.h"
#include "src/piex.h"
#include "src/piex_cr3.h"
#include "src/piex_types.h"

namespace piex {

size_t BytesRequiredForIsRaw() { return ::piex_rust::BytesRequiredForIsRaw(); }

bool IsRaw(StreamInterface* data) { return ::piex_rust::IsRaw(data); }

Error GetPreviewImageData(StreamInterface* data,
                          PreviewImageData* preview_image_data,
                          image_type_recognition::RawImageTypes* output_type) {
  return ::piex_rust::GetPreviewImageData(data, preview_image_data,
                                          output_type);
}

bool GetDngInformation(StreamInterface* data, std::uint32_t* width,
                       std::uint32_t* height,
                       std::vector<std::uint32_t>* cfa_pattern_dim) {
  return ::piex_rust::GetDngInformation(data, width, height, cfa_pattern_dim);
}

bool GetOrientation(StreamInterface* data, std::uint32_t* orientation) {
  return ::piex_rust::GetOrientation(data, orientation);
}

std::vector<std::string> SupportedExtensions() {
  return ::piex_rust::SupportedExtensions();
}

bool Cr3GetOrientation(StreamInterface* stream, std::uint32_t* orientation) {
  return ::piex_rust::Cr3GetOrientation(stream, orientation);
}

Error Cr3GetPreviewData(StreamInterface* stream,
                        PreviewImageData* preview_image_data) {
  return ::piex_rust::Cr3GetPreviewData(stream, preview_image_data);
}

}  // namespace piex
