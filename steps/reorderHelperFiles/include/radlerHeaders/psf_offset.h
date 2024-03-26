// SPDX-License-Identifier: LGPL-3.0-only

#ifndef RADLER_PSF_OFFSET_H_
#define RADLER_PSF_OFFSET_H_

#include <iomanip>
#include <ostream>

namespace radler {

class PsfOffset {
 public:
  explicit PsfOffset(size_t x, size_t y) : x(x), y(y) {}
  PsfOffset() = default;

  PsfOffset(const PsfOffset&) = default;
  PsfOffset& operator=(const PsfOffset&) = default;

  PsfOffset(PsfOffset&&) = default;
  PsfOffset& operator=(PsfOffset&&) = default;

  /** The x-offset in pixels from the corner position. */
  size_t x{0};

  /** The y-offset in pixels from the corner position. */
  size_t y{0};

  friend std::ostream& operator<<(std::ostream& output, const PsfOffset& psf) {
    return output << "[x: " << psf.x << ", y: " << psf.y << ']';
  }
};

}  // namespace radler

#endif  // RADLER_PSF_OFFSET_H_
