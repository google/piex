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

#include "src/piex_types.h"

#ifndef PIEX_PIEX_CR3_H_
#define PIEX_PIEX_CR3_H_

namespace piex {

// Gets the EXIF orientation of a CR3 stream, returning true on success.
bool Cr3GetOrientation(StreamInterface* stream, std::uint32_t* orientation);

// Gets preview images of a CR3 stream, returning kOk on success. Assumes the
// stream is a CR3 stream.
//
// Canon's CR3 is based on ISO/IEC 14496-12: ISO base media file format. (CR2 is
// TIFF based.) A Canon CR3 contains multiple embedded images. Most cameras
// output CR3 files that contain a full-size JPEG, a 1620x1080 preview JPEG, and
// a 160x120 thumbnail JPEG.
// The Canon EOS 1D X Mark III, though, contains a full-size HEVC image, a
// 1620x1080 preview JPEG, and a 160x120 thumbnail JPEG.
// Until support for HEVC is added, this method returns the largest embedded
// JPEG in preview_image_data->preview.
//
Error Cr3GetPreviewData(StreamInterface* stream,
                        PreviewImageData* preview_image_data);
}  // namespace piex

#endif  // PIEX_PIEX_CR3_H_
