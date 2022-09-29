// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_STEPS_TEST_UNIT_THROWSTEP_H_
#define DP3_STEPS_TEST_UNIT_THROWSTEP_H_

#include "../../../Step.h"

namespace dp3 {
namespace steps {
namespace test {
/**
 * Common Step class for use in tests.
 * All methods raise "Unexpected call" errors when called.
 */
class ThrowStep : public Step {
 public:
  ThrowStep();
  ~ThrowStep() override;

  common::Fields getRequiredFields() const override;
  common::Fields getProvidedFields() const override;

  void updateInfo(const base::DPInfo&) override;
  bool process(const base::DPBuffer&) override;
  bool process(std::unique_ptr<base::BDABuffer>) override;
  void finish() override;
  void show(std::ostream&) const override;
};
}  // namespace test
}  // namespace steps
}  // namespace dp3

#endif
