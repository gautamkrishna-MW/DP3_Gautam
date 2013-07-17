//# ApplyCal.h: DPPP step class to apply a calibration correction to the data
//# Copyright (C) 2013
//# ASTRON (Netherlands Institute for Radio Astronomy)
//# P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
//#
//# This file is part of the LOFAR software suite.
//# The LOFAR software suite is free software: you can redistribute it and/or
//# modify it under the terms of the GNU General Public License as published
//# by the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The LOFAR software suite is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# GNU General Public License for more details.
//#
//# You should have received a copy of the GNU General Public License along
//# with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.
//#
//# $Id: ApplyCal.h 21598 2012-07-16 08:07:34Z diepen $
//#
//# @author Ger van Diepen

#ifndef DPPP_APPLYCAL_H
#define DPPP_APPLYCAL_H

// @file
// @brief DPPP step class to apply a calibration correction to the data

#include <DPPP/DPInput.h>
#include <DPPP/DPBuffer.h>
#include <ParmDB/ParmFacade.h>
#include <ParmDB/ParmSet.h>
#include <ParmDB/Parm.h>
#include <casa/Arrays/Cube.h>
#include <casa/Arrays/ArrayMath.h>

namespace LOFAR {

  class ParameterSet;

  namespace DPPP {
    // @ingroup NDPPP

    // This class is a DPStep class applying calibration parameters to the data.

    class ApplyCal: public DPStep
    {
    public:
      // Construct the object.
      // Parameters are obtained from the parset using the given prefix.
      ApplyCal (DPInput*, const ParameterSet&, const string& prefix);

      virtual ~ApplyCal();

      // Process the data.
      // It keeps the data.
      // When processed, it invokes the process function of the next step.
      virtual bool process (const DPBuffer&);

      // Finish the processing of this step and subsequent steps.
      virtual void finish();

      // Update the general info.
      virtual void updateInfo (const DPInfo&);

      // Show the step parameters.
      virtual void show (std::ostream&) const;

      // Show the timings.
      virtual void showTimings (std::ostream&, double duration) const;


    private:
      void applyTEC (casa::Complex* vis, const double tec, const double freq);
      void applyClock (casa::Complex* vis, const double clockA,
                                     const double clockB);
      void applyBandpass (casa::Complex* vis, const casa::DComplex* lhs,
                                        const casa::DComplex* rhs);
      void applyRM (casa::Complex* vis, const casa::DComplex* lhs,
                                  const casa::DComplex* rhs);
      void applyJones (casa::Complex* vis, const casa::DComplex* lhs,
                                     const casa::DComplex* rhs);
      void applyGain (casa::Complex* vis, int ant1, int ant2, int chan,
                       int time);

      void updateParms (const double bufStartTime);

      void invert (casa::DComplex* v) const;
      void invert (casa::Complex* v) const;

      //# Data members.
      DPInput*         itsInput;
      string           itsName;
      string           itsParmDBName;
      boost::shared_ptr<BBS::ParmFacade> itsParmDB;
      string           itsCorrectType;
      // Expressions to search for in itsParmDB

      vector<casa::String>   itsParmExprs;

      // itsParms contains the parameters to a grid, first for all parameters
      // (e.g. Gain:0:0 and Gain:1:1), next all antennas, next over frec * time
      // as returned by ParmDB
      vector<vector<vector<double> > > itsParms;
      int             itsBufStep;
      int             itsNCorr;
      double          itsSigma;
      double          itsTimeInterval;
      double          itsLastTime;
      bool            itsUseAP;      //# use ampl/phase or real/imag
      bool            itsHasCrossGain;  //# use terms Gain:0:1:* and Gain:1:0:*
      NSTimer          itsTimer;
    };

  } //# end namespace
}

#endif
