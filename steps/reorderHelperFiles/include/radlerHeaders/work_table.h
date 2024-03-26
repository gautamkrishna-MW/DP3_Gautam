// SPDX-License-Identifier: LGPL-3.0-only

#ifndef RADLER_WORK_TABLE_H_
#define RADLER_WORK_TABLE_H_

#include "psf_offset.h"
#include "work_table_entry.h"

#include <functional>
#include <memory>
#include <ostream>
#include <vector>

namespace radler {
/**
 * The WorkTable contains WorkTableEntry's and groups entries
 * that have the same squaredDeconvolutionIndex.
 */
class WorkTable {
 public:
  using Entries = std::vector<std::unique_ptr<WorkTableEntry>>;
  using Group = std::vector<const WorkTableEntry*>;

  /**
   * Iterator-like class which (only) supports a range-based loop over entries.
   *
   * Dereferencing this class yields a reference to the actual object instead
   * of a reference to the pointer for the object.
   */
  class EntryIteratorLite {
    using BaseIterator = Entries::const_iterator;

   public:
    explicit EntryIteratorLite(BaseIterator base_iterator)
        : base_iterator_(base_iterator) {}

    EntryIteratorLite(const EntryIteratorLite&) = default;
    EntryIteratorLite(EntryIteratorLite&&) = default;
    EntryIteratorLite& operator=(const EntryIteratorLite&) = default;
    EntryIteratorLite& operator=(EntryIteratorLite&&) = default;

    const WorkTableEntry& operator*() const { return **base_iterator_; }
    EntryIteratorLite& operator++() {
      ++base_iterator_;
      return *this;
    }
    bool operator!=(const EntryIteratorLite& other) const {
      return base_iterator_ != other.base_iterator_;
    }
    bool operator==(const EntryIteratorLite& other) const {
      return base_iterator_ == other.base_iterator_;
    }

   private:
    BaseIterator base_iterator_;
  };

  /**
   * @brief Constructs a new WorkTable object.
   *
   * @param n_original_groups The number of original channel groups. When adding
   * entries, their original channel index must be less than the number of
   * original groups. If the value is zero, one group is used.
   * @param n_deconvolution_groups The number of deconvolution groups, which is
   * the number of channels used during deconvolution.
   * A deconvolution group consist of one or more original channel groups, which
   * are then joinedly deconvolved by averaging them before deconvolution and
   * interpolating them after deconvolution.
   * If the value is zero, or larger than the number of original groups,
   * all channels are deconvolved separately.
   * @param channel_index_offset The index of the first channel in the caller.
   */
  explicit WorkTable(std::vector<PsfOffset> psf_offsets,
                     std::size_t n_original_groups,
                     std::size_t n_deconvolution_groups,
                     std::size_t channel_index_offset = 0);

  // TODO(AST-912) Make copy/move operations Google Style compliant.
  WorkTable(const WorkTable&) = default;
  WorkTable(WorkTable&&) = default;
  WorkTable& operator=(const WorkTable&) = delete;
  WorkTable& operator=(WorkTable&&) = delete;

  /**
   * @return The table entries, grouped by their original channel index.
   * @see AddEntry()
   */
  const std::vector<Group>& OriginalGroups() const { return original_groups_; }

  /**
   * @return The original group indices for each deconvolution group.
   */
  const std::vector<std::vector<std::size_t>>& DeconvolutionGroups() const {
    return deconvolution_groups_;
  }

  /**
   * Find the first group of original channels, given a deconvolution group
   * index.
   *
   * @param deconvolution_index Index for a deconvolution group. Must be less
   * than the number of deconvolution groups.
   * @return A reference to the first original group for the deconvolution
   * group.
   */
  const Group& FirstOriginalGroup(size_t deconvolution_index) const {
    return original_groups_[deconvolution_groups_[deconvolution_index].front()];
  }

  EntryIteratorLite Begin() const {
    return EntryIteratorLite(entries_.begin());
  }
  EntryIteratorLite End() const { return EntryIteratorLite(entries_.end()); }

  /**
   * @brief Adds an entry to the table.
   *
   * The original channel index of the entry determines the original group for
   * the entry. It must be less than the number of original channel groups, as
   * given in the constructor.
   *
   * @param entry A new entry.
   */
  void AddEntry(std::unique_ptr<WorkTableEntry> entry);

  /**
   * @return A reference to the first entry.
   */
  const WorkTableEntry& Front() const { return *entries_.front(); }

  /**
   * @return The number of entries in the table.
   */
  size_t Size() const { return entries_.size(); }

  /**
   * @return The channel index offset, which was set in the constructor.
   */
  size_t GetChannelIndexOffset() const { return channel_index_offset_; }

  const std::vector<PsfOffset>& PsfOffsets() const noexcept {
    return psf_offsets_;
  }

  /**
   * Validates the invariant for @ref psf_offsets_ and the
   * @ref WorkTableEntry::psf_accessors in @ref entries_.
   * @throw std::runtime_error If the worktable is incorrect.
   */
  void ValidatePsfs() const;

 private:
  Entries entries_;

  /**
   * The direction-dependent PSF offsets.
   *
   * When no direction-dependent PSF is used the vector is empty.
   *
   * All @ref Worktable::entries_ use the same PSF offsets. The number of PSF
   * accessors shall be equal to the number of elements of this vector, or
   * shall be 1 when no direction-dependant PSF is used.  When the
   * deconvolution is executed this invariant shall be valid. This is validated
   * by @ref ValidatePsfs when calling @ref Radler::Perform.
   */
  std::vector<PsfOffset> psf_offsets_;

  /**
   * A user of the WorkTable may use different channel indices than
   * the WorkTable. This offset is the difference between those
   * indices.
   * For example, with three channels, the WorkTable indices are always
   * 0, 1, and 2. When the user indices are 4, 5, and 6, this offset will be 4.
   */
  const std::size_t channel_index_offset_;

  /**
   * An original group has entries with equal original channel indices.
   */
  std::vector<Group> original_groups_;

  /**
   * A deconvolution group consists of one or more original groups, which
   * are deconvolved together. Each entry contains the indices of the original
   * groups that are part of the deconvolution group.
   */
  std::vector<std::vector<std::size_t>> deconvolution_groups_;

  /**
   * begin() and end() allow writing range-based loops over all entries.
   * @{
   */
  friend EntryIteratorLite begin(const WorkTable& table) {
    return table.Begin();
  }
  friend EntryIteratorLite end(const WorkTable& table) { return table.End(); }
  /** @} */

  friend std::ostream& operator<<(std::ostream& output,
                                  const WorkTable& work_table) {
    output << "=== IMAGING TABLE ==="
           << "\nOriginal groups       " << work_table.original_groups_.size()
           << "\nDeconvolution groups  "
           << work_table.deconvolution_groups_.size()
           << "\nChannel index         " << work_table.channel_index_offset_
           << '\n';

    if (!work_table.entries_.empty()) {
      output << "   # Pol Ch Mask Interval Weight Freq(MHz)\n";

      for (const auto& entry : work_table.entries_) {
        output << *entry;
      }
    }

    if (!work_table.psf_offsets_.empty()) {
      output << "=== PSFs ===\n";
      for (const auto& psf : work_table.psf_offsets_) {
        output << psf << '\n';
      }
    }
    return output;
  }
};
}  // namespace radler

#endif
