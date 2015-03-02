// -*- LSST-C++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2015 AURA/LSST
 *
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the LSST License Statement and
 * the GNU General Public License along with this program.  If not,
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */

#ifndef LSST_MEAS_EXTENSIONS_SHAPEHSM_HsmShapeControl_h_INCLUDED
#define LSST_MEAS_EXTENSIONS_SHAPEHSM_HsmShapeControl_h_INCLUDED

#include "lsst/pex/config.h"
#include "lsst/meas/base/Algorithm.h"
#include "lsst/meas/base/ShapeUtilities.h"
#include "lsst/meas/base/CentroidUtilities.h"
#include "lsst/meas/base/FlagHandler.h"
#include "lsst/meas/base/InputUtilities.h"

namespace lsst { namespace meas { namespace extensions { namespace shapeHSM {

enum MeasType {
    ELLIPTICITY,                        // We will measure e1,e2
    SHEAR,                              // We will measure g1,g2
};

/**
 *  @brief A C++ control class to handle HsmShapeAlgorithm's configuration
 *
 *  In C++, we define Control objects to handle configuration information.  Using the LSST_CONTROL_FIELD
 *  macro and lsst.pex.config.wrap.makeConfigClass, we can turn these into more full-featured Config classes
 *  in Python.  While the user will usually interact with the Config class, the plugin wrapper system will
 *  turn Config instances into Control instances when passing them to C++.
 *
 *  This should logically be an inner class, but Swig doesn't know how to parse those.
 */
class HsmShapeControl {
public:
    LSST_CONTROL_FIELD(deblendNChild, std::string, "Field name for number of deblend children");
    LSST_CONTROL_FIELD(badMaskPlanes, std::vector<std::string>, "Mask planes that indicate pixels that should be excluded from the fit");

    /**
     *  @brief Default constructor
     *
     *  All control classes should define a default constructor that sets all fields to their default values.
     */
    HsmShapeControl() : _name("") {
        badMaskPlanes.push_back("BAD");
        badMaskPlanes.push_back("SAT");
        badMaskPlanes.push_back("INTRP");
    }

    HsmShapeControl(std::string const & name) : _name(name) {
        badMaskPlanes.push_back("BAD");
        badMaskPlanes.push_back("SAT");
        badMaskPlanes.push_back("INTRP");
    }

private:
    std::string _name;
};

class HsmShapeBjControl : public HsmShapeControl {
public:
    HsmShapeBjControl() : HsmShapeControl("ext_shapeHSM_HsmShapeBj") {}
};

class HsmShapeLinearControl : public HsmShapeControl {
public:
    HsmShapeLinearControl() : HsmShapeControl("ext_shapeHSM_HsmShapeLinear") {}
};

class HsmShapeKsbControl : public HsmShapeControl {
public:
    HsmShapeKsbControl() : HsmShapeControl("ext_shapeHSM_HsmShapeKsb") {}
};

class HsmShapeRegaussControl : public HsmShapeControl {
public:
    HsmShapeRegaussControl() : HsmShapeControl("ext_shapeHSM_HsmShapeRegauss") {}
};

/*
 * @brief A measurement algorithm for measuring HSM shape 
 * HSM shape algorithm class - we use one class for all algorithms; all the specialized
 * work is done by the control class _makeAlgorithm implementations.
 *
 * We don't inherit from ShapeAlgorithm because none of these compute any errors,
 * and we don't want to add fields that are always NaNs.  We probably need to think
 * about how to handle the "is-a" ShapeAlgorithm test a little better, but we don't
 * rely on it anywhere presently.
 */
class HsmShapeAlgorithm : public base::SimpleAlgorithm {
public:
    enum {
        FAILURE=base::FlagHandler::FAILURE,
        NO_PIXELS,
        NOT_CONTAINED,
        HAS_DEBLEND,
        IGNORE_PARENT_SOURCE,
        N_FLAGS
    };

    /// A typedef to the Control object for this algorithm, defined above.
    /// The control object contains the configuration parameters for this algorithm.
    typedef HsmShapeControl Control;

    HsmShapeAlgorithm(Control const & ctrl, std::string const & name, std::string const & shearType, MeasType measType,
                      std::string const & doc, afw::table::Schema & schema);

protected:
    void measure(
        afw::table::SourceRecord & measRecord,
        afw::image::Exposure<float> const & exposure
    ) const;

    void fail(
        afw::table::SourceRecord & measRecord,
        meas::base::MeasurementError * error=NULL
    ) const;

private:
    Control _ctrl;
    std::string _shearType;
    MeasType _measType;
    std::string _doc;

    base::CentroidResultKey _centroidResultKey;
    base::ShapeResultKey _momentsKey;
    base::ShapeResultKey _psfMomentsKey;
    afw::table::Key<double> _e1Key;
    afw::table::Key<double> _e2Key;
    afw::table::Key<double> _sigmaKey;
    afw::table::Key<double> _resolutionKey;
    afw::table::Key<int> _deblendKey;
    base::FlagHandler _flagHandler;
    base::SafeCentroidExtractor _centroidExtractor;
    bool _hasDeblendKey;
};

class HsmShapeBjAlgorithm : public HsmShapeAlgorithm {
public:
    typedef HsmShapeBjControl Control;
    HsmShapeBjAlgorithm(Control const & ctrl, std::string const & name, afw::table::Schema & schema) :
        HsmShapeAlgorithm( ctrl, name, "BJ", ELLIPTICITY, "PSF-corrected shear using Bernstein & Jarvis (2002) method", schema) {}
};

class HsmShapeLinearAlgorithm : public HsmShapeAlgorithm {
public:
    typedef HsmShapeLinearControl Control;
    HsmShapeLinearAlgorithm(Control const & ctrl, std::string const & name, afw::table::Schema & schema) :
        HsmShapeAlgorithm(ctrl, name, "LINEAR", ELLIPTICITY, "PSF-corrected shear using Hirata & Seljak (2003) 'linear' method", schema) {}
};

class HsmShapeKsbAlgorithm : public HsmShapeAlgorithm {
public:
    typedef HsmShapeKsbControl Control;
    HsmShapeKsbAlgorithm(Control const & ctrl, std::string const & name, afw::table::Schema & schema) :
        HsmShapeAlgorithm(ctrl, name, "KSB", SHEAR, "PSF-corrected shear using KSB method",schema) {}
};

class HsmShapeRegaussAlgorithm : public HsmShapeAlgorithm {
public:
    typedef HsmShapeRegaussControl Control;
    explicit HsmShapeRegaussAlgorithm(Control const & ctrl, std::string const & name, afw::table::Schema & schema) :
        HsmShapeAlgorithm(ctrl, name, "REGAUSS", ELLIPTICITY, "PSF-corrected shear using Hirata & Seljak (2003) 'regaussianization' method", schema) {}
};
}}}} // namespace lsst::meas::extensions::shapeHSM

#endif // !LSST_MEAS_EXTENSIONS_SHAPEHSM_HsmShapeControl_h_INCLUDED
