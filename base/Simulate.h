// Simulate.h: Simulate visibilities for a patch of sources.
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// \file
/// Simulate visibilities for a patch of sources.

#ifndef DPPP_SIMULATE_H
#define DPPP_SIMULATE_H

#include <casacore/casa/Arrays/Matrix.h>
#include <casacore/casa/Arrays/Vector.h>

#include "Baseline.h"
#include "Cursor.h"
#include "Patch.h"
#include "Direction.h"

namespace dp3 {
namespace base {

/// Setup the splitting of the baseline UVWs into station UVWs. It returns
/// the indices of the baselines needed to split the baseline UVWs into
/// station UVWs. They are in such an order that the UVW of a station is known
/// before used in another baseline to derive the UVW of the other station.
/// It can handle cases where baselines occur in disjoint station groups
/// like 0-1, 0-2, 1-2 and 3-4, 4-5, 5-6.
/// Note that the first station of a group gets UVW=0. All other station UVWs
/// are relative to it using the baseline UVWs.
/// Also note that nr of groups can be derived from the size of the returned
/// vector (because it contains no entry for the first antenna in a group).
std::vector<int> nsetupSplitUVW(unsigned int nant,
                                const casacore::Vector<int>& ant1,
                                const casacore::Vector<int>& ant2);

/// Do the actual splitting of baseline UVWs into station UVWs using
/// the index vector generated by setupSplitUVW.
void nsplitUVW(const std::vector<int>& blindex,
               const std::vector<Baseline>& baselines,
               const casacore::Matrix<double>& uvwbl,
               casacore::Matrix<double>& uvwant);

/// Split baseline UVW coordinates into station UVW coordinates (by assuming the
/// station with index 0 has UVW coordinates of (0, 0, 0)).
///
/// \param[in]   nStation
/// The number of stations.
/// \param[in]   nBaseline
/// The number of baselines.
/// \param[in]   baselines
/// A cursor for a 1-D buffer of baselines of shape (\p nBaseline).
/// \param[in]   uvw
/// A cursor for a 2-D buffer of UVW coordinates of shape (\p nBaseline, 3).
/// \param[in]   split
/// A cursor for a 2-D buffer of station UVW coordinates of shape
/// (\p nStation, 3).
void splitUVW(size_t nStation, size_t nBaseline,
              const_cursor<Baseline> baselines, const_cursor<double> uvw,
              cursor<double> split);

/// Transform UVW coordinates from phase reference position \p from to phase
/// reference position \p to. The transformation is performed in place.
///
/// \param[in]   from
/// Current phase reference position for the UVW coordinates.
/// \param[in]   to
/// New phase reference position for the UVW coordinates.
/// \param[in]   nUVW
/// The number of UVW coordinates to transform.
/// \param[in]   uvw
/// A 2-D buffer of UVW coordinates of shape (\p UVW, 3).
void rotateUVW(const Direction& from, const Direction& to, size_t nUVW,
               double* uvw);

/// Simulate visibilities for a patch of sources. The computed visibilities are
/// added to \p vis.
///
/// \param[in]   reference
/// Phase reference position.
/// \param[in]   patch
/// Patch of sources to simulate visibilities for.
/// \param[in]   nStation
/// The number of stations.
/// \param[in]   nBaseline
/// The number of baselines.
/// \param[in]   nChannel
/// The number of frequency channels.
/// \param[in]   baselines
/// A cursor for a 1-D buffer of baselines of shape (\p nBaseline).
/// \param[in]   freq
/// A cursor for a 1-D buffer of channel frequencies of shape (\p nChannel).
/// \param[in]   uvw
/// A cursor for a 2-D buffer of station UVW coordinates of shape
/// (\p nStation, 3).
/// \param[in]   buffer
/// A cursor for a 3-D buffer of shape (\p nBaseline, \p nChannel, 4) into which
/// the simulated visibilities will be written.
void simulate(const Direction& reference, const Patch::ConstPtr& patch,
              size_t nStation, size_t nBaseline, size_t nChannel,
              const_cursor<Baseline> baselines, const_cursor<double> freq,
              const_cursor<double> uvw, cursor<dcomplex> buffer);

}  // namespace base
}  // namespace dp3

#endif
