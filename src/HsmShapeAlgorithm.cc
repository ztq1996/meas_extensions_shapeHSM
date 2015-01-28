// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2013 LSST Corporation.
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

#include "lsst/pex/exceptions.h"
#include "lsst/pex/logging/Trace.h"
#include "lsst/afw/image.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/meas/algorithms/Measure.h"

#include "lsst/meas/extensions/shapeHSM/HsmShapeAdapter.h"

#include "lsst/afw/table/Source.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/detection/FootprintArray.h"
#include "lsst/afw/detection/FootprintArray.cc"
#include "lsst/meas/extensions/shapeHSM/HsmShapeAlgorithm.h"

namespace lsst { namespace meas { namespace extensions { namespace shapeHSM {

HsmShapeShapeletAlgorithm::HsmShapeShapeletAlgorithm(
    Control const & ctrl,
    std::string const & name,
    afw::table::Schema & schema
    ) : HsmShapeAlgorithm(ctrl, name + "Shapelet", "SHAPELET", 'g', "PSF-corrected shear using KSB method", schema) {}

HsmShapeAlgorithm::HsmShapeAlgorithm(
    Control const & ctrl,
    std::string const & name,
    std::string const & shearType,
    char measType,
    std::string const & doc,
    afw::table::Schema & schema
) : _ctrl(ctrl),
    _shearType(shearType),
    _measType(measType),
    _doc(doc), 
    _centroidExtractor(schema, name)
{
    static boost::array<base::FlagDefinition,N_FLAGS> const flagDefs = {{
        {"flag", "general failure flag"},
    }};
    _flagHandler = base::FlagHandler::addFields(schema, name, flagDefs.begin(), flagDefs.end());
}

void HsmShapeAlgorithm::measure(
    afw::table::SourceRecord & measRecord,
    afw::image::Exposure<float> const & exposure
) const {
    PTR(afw::detection::Psf const) psf = exposure.getPsf();
    if (!psf) {
        throw LSST_EXCEPT(
            base::FatalAlgorithmError,
            "HsmShape algorithm requires a Psf with every exposure"
        );
    }
    afw::geom::Point2D position = _centroidExtractor(measRecord, _flagHandler);
    PTR(afw::detection::Psf::Image) psfImage = psf->computeImage(position);
    afw::geom::Box2I fitBBox = psfImage->getBBox();
    fitBBox.clip(exposure.getBBox());

    std::vector<std::string> const & badMaskPlanes 
        = _ctrl.badMaskPlanes;

    afw::image::MaskPixel badPixelMask 
        = exposure.getMaskedImage().getMask()->getPlaneBitMask(badMaskPlanes);
    HsmShapeAdapter< afw::image::Exposure<float> > shearEst(
        exposure, position, *measRecord.getFootprint(), badPixelMask
    );
    
    short status = shearEst.measure(_shearType);
    base::CentroidResult cent;
    cent.setCentroid(shearEst.getCentroid());
    measRecord.set(_centroidResultKey, cent);
    base::ShapeResult shape;
    shape.setShape(shearEst.getMoments());
    measRecord.set(_momentsKey, shape);
    base::ShapeResult psfShape;
    psfShape.setShape(shearEst.getPsfMoments());

    measRecord.set(_sigmaKey, shearEst.getSigma());
    measRecord.set(_e1Key, shearEst.getE1());
    measRecord.set(_e2Key, shearEst.getE2());
    measRecord.set(_resolutionKey, shearEst.getResolution());


    char meas_type = shearEst.getMeasType();
    if (meas_type == 'e') {
        measRecord.set(_errKey, 0.5 * shearEst.getShearSig());
    } else if (meas_type == 'g') {
        measRecord.set(_errKey, shearEst.getShearSig());
    }

    // FIXME measRecord.set(_flagKey, status); 
}

void HsmShapeAlgorithm::fail(afw::table::SourceRecord & measRecord, base::MeasurementError * error) const {
    _flagHandler.handleFailure(measRecord, error);
}


}}}} // namespace lsst::meas::extensions::shapeHSM
