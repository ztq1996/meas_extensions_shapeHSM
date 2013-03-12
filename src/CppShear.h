// -*- lsst-c++ -*-
#ifndef LSST_MEAS_EXTENSIONS_SHAPEHSM_CppShear_h_INCLUDED
#define LSST_MEAS_EXTENSIONS_SHAPEHSM_CppShear_h_INCLUDED

/*
 * This is an adapter file that provides a minimal reimplementation
 * of the GalSim CppShear object using LSST objects, defined here
 * so that we can use PSFCorr.[h,cc] with minimal modification.
 */

#include "lsst/afw/geom/ellipses/Distortion.h"

namespace galsim {

struct CppShear {

    lsst::afw::geom::ellipses::Distortion distortion;

    void setE1E2(double e1, double e2) {
        distortion.setE1(e1);
        distortion.setE2(e2);
    }

};

} // namespace galsim

#endif // !LSST_MEAS_EXTENSIONS_SHAPEHSM_CppShear_h_INCLUDED
