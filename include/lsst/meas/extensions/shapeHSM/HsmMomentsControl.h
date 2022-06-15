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

#ifndef LSST_MEAS_EXTENSIONS_SHAPEHSM_HsmMomentsControl_h_INCLUDED
#define LSST_MEAS_EXTENSIONS_SHAPEHSM_HsmMomentsControl_h_INCLUDED

#include "lsst/pex/config.h"
#include "lsst/meas/base/Algorithm.h"
#include "lsst/meas/base/ShapeUtilities.h"
#include "lsst/meas/base/CentroidUtilities.h"
#include "lsst/meas/base/FlagHandler.h"
#include "lsst/meas/base/InputUtilities.h"
#include "lsst/afw/detection/Psf.h"

namespace lsst { namespace meas { namespace extensions { namespace shapeHSM {

// This is here just as a parent class for the following
class HsmMomentsControl {
};

class HsmSourceMomentsControl {
public:
    LSST_CONTROL_FIELD(badMaskPlanes, std::vector<std::string>, "Mask planes used to reject bad pixels.");
    LSST_CONTROL_FIELD(roundMoments, bool, "Use round weight function?");
    LSST_CONTROL_FIELD(addFlux, bool, "Store measured flux?");

    HsmSourceMomentsControl(HsmSourceMomentsControl const & other) :
        badMaskPlanes(other.badMaskPlanes), roundMoments(other.roundMoments), addFlux(other.addFlux) {}

    HsmSourceMomentsControl() :
        roundMoments(false),
        addFlux(false)
    {
        badMaskPlanes.push_back("BAD");
        badMaskPlanes.push_back("SAT");
    }
};


class HsmSourceMomentsRoundControl : public HsmSourceMomentsControl {
public:
    HsmSourceMomentsRoundControl(HsmSourceMomentsRoundControl const & other) : HsmSourceMomentsControl(other) {}

    HsmSourceMomentsRoundControl() : HsmSourceMomentsControl() {
        roundMoments = true;
        addFlux = true;
    }
};


class HsmPsfMomentsControl {
public:
    LSST_CONTROL_FIELD(
        useSourceCentroidOffset,
        bool,
        "Include subpixel source centroid offset when drawing PSF to be measured?"
    );

    HsmPsfMomentsControl() : useSourceCentroidOffset(false) {}

    HsmPsfMomentsControl(HsmPsfMomentsControl const & other) :
        useSourceCentroidOffset(other.useSourceCentroidOffset) {}
};


class HsmPsfMomentsDebiasedControl : public HsmPsfMomentsControl {
public:
    LSST_CONTROL_FIELD(
        noiseSource,
        std::string,
        "Noise source.  How to choose variance of the zero-mean Gaussian noise added to image.\n"
        "Allowed values:\n"
        "  'meta': variance = the 'BGMEAN' metadata entry\n"
        "  'variance': variance = the image's variance plane'\n"
    );
    LSST_CONTROL_FIELD(seedOffset, int, "Random seed offset");
    LSST_CONTROL_FIELD(badMaskPlanes, std::vector<std::string>, "Mask planes used to reject bad pixels.");

    HsmPsfMomentsDebiasedControl() :
        HsmPsfMomentsControl(),
        noiseSource("variance"),
        seedOffset(0)
    {
        useSourceCentroidOffset=true;
        badMaskPlanes.push_back("BAD");
        badMaskPlanes.push_back("SAT");
    }

    HsmPsfMomentsDebiasedControl(HsmPsfMomentsDebiasedControl const & other) :
        HsmPsfMomentsControl(other),
        noiseSource(other.noiseSource),
        seedOffset(other.seedOffset),
        badMaskPlanes(other.badMaskPlanes)
    {}
};


/// Base class to measure HSM adaptive moments
///
/// Use this to consolidate common code for HsmSourceMoments and HsmPsfMoments
class HsmMomentsAlgorithm : public base::SimpleAlgorithm {
public:

    // Structures and routines to manage flaghandler
    static base::FlagDefinitionList const & getFlagDefinitions();
    static base::FlagDefinition const FAILURE;
    static base::FlagDefinition const NO_PIXELS;
    static base::FlagDefinition const NOT_CONTAINED;
    static base::FlagDefinition const PARENT_SOURCE;
    static base::FlagDefinition const GALSIM;

    typedef HsmMomentsControl Control;

protected:
    HsmMomentsAlgorithm(std::string const & name, afw::table::Schema & schema, char const* doc) :
    _doc(doc),
    _centroidResultKey(
        base::CentroidResultKey::addFields(schema, name, doc, base::NO_UNCERTAINTY)),
    _momentsKey(
        base::ShapeResultKey::addFields(schema, name, doc, base::NO_UNCERTAINTY)
    ),
    _centroidExtractor(schema, name)
    {
        _flagHandler = base::FlagHandler::addFields(schema, name, getFlagDefinitions());
    }

    /// Calculate moments
    template<typename PixelT>
    void calculate(
        afw::table::SourceRecord& source, // Source for recording moments
        std::shared_ptr<afw::image::Image<PixelT>> const& afwImage, // Image on which to measure moments
        std::shared_ptr<afw::image::Mask<afw::image::MaskPixel>> const& afwMask, // Mask for image
        geom::Box2I const& bbox,     // Bounding box
        geom::Point2D const& center, // Starting center for measuring moments
        afw::image::MaskPixel const badPixelMask, // Bitmask for bad pixels
        float const width,            // PSF width estimate, for starting moments
        bool roundMoments=false, // Use round weight function
        bool addFlux=false, // add Flux to output
        bool subtractCenter=false // subtract starting center from x/y outputs?
    ) const;

public:
    void fail(
        afw::table::SourceRecord & measRecord,
        meas::base::MeasurementError * error=NULL
    ) const;

protected:
    std::string _doc;
    base::CentroidResultKey _centroidResultKey;
    base::ShapeResultKey _momentsKey;
    base::FlagHandler _flagHandler;
    base::SafeCentroidExtractor _centroidExtractor;
    afw::table::Key<float> _fluxKey;
};


/// Class to measure HSM adaptive moments of source
class HsmSourceMomentsAlgorithm : public HsmMomentsAlgorithm {
public:
    /// A typedef to the Control object for this algorithm, defined above.
    ///  The control object contains the configuration parameters for this algorithm.
    typedef HsmSourceMomentsControl Control;

    /// @brief Initialize with standard field names and customized documentation.
    HsmSourceMomentsAlgorithm(Control const & ctrl, std::string const & name, afw::table::Schema & schema) :
        HsmMomentsAlgorithm(name, schema, "Source adaptive moments algorithm from HSM"), _ctrl(ctrl)
    {
        if (ctrl.addFlux) {
            _fluxKey = schema.addField<float>(name + "_Flux", "HSM flux");
        }
    }

    void measure(
        afw::table::SourceRecord & measRecord,
        afw::image::Exposure<float> const & exposure
    ) const;

private:
    Control _ctrl;
};


/// Class to measure HSM adaptive moments of PSF
class HsmPsfMomentsAlgorithm : public HsmMomentsAlgorithm {
public:
    /// A typedef to the Control object for this algorithm, defined above.
    ///  The control object contains the configuration parameters for this algorithm.
    typedef std::shared_ptr<const HsmPsfMomentsControl> Control;

    /// @brief Initialize with standard field names and customized documentation.
    HsmPsfMomentsAlgorithm(Control const & ctrl, std::string const & name, afw::table::Schema & schema) :
        HsmMomentsAlgorithm(name, schema, "Psf adaptive moments algorithm from HSM"), _ctrl(ctrl) {}

    void measure(
        afw::table::SourceRecord & measRecord,
        afw::image::Exposure<float> const & exposure
    ) const;

protected:
    virtual std::shared_ptr<afw::detection::Psf::Image> getPsfImage(
        geom::Point2D center,
        afw::table::SourceRecord & source,
        afw::image::Exposure<float> const & exposure
    ) const;

    virtual std::shared_ptr<afw::image::Mask<afw::image::MaskPixel>> getPsfMask(
        geom::Point2D center,
        afw::table::SourceRecord & source,
        afw::image::Exposure<float> const & exposure
    ) const;

    virtual afw::image::MaskPixel const getBadPixelMask(
        afw::image::Exposure<float> const & exposure
    ) const { return 0; }

    Control _ctrl;
};


/// Measure PSF moments while making image signal-to-noise look similar to underlying star,
/// thereby removing any signal-to-noise related biases.
class HsmPsfMomentsDebiasedAlgorithm : public HsmPsfMomentsAlgorithm {
public:
    static base::FlagDefinition const EDGE;

    typedef std::shared_ptr<const HsmPsfMomentsDebiasedControl> Control;

    HsmPsfMomentsDebiasedAlgorithm(Control const & ctrl, std::string const & name, afw::table::Schema & schema) :
        HsmPsfMomentsAlgorithm(ctrl, name, schema) {
            auto thisCtrl = std::static_pointer_cast<Control::element_type>(_ctrl);
            if ((thisCtrl->noiseSource != "meta") && (thisCtrl->noiseSource != "variance")) {
                throw LSST_EXCEPT(base::MeasurementError, "invalid noiseSorce", FAILURE.number);
            }
        }

private:
    std::shared_ptr<afw::detection::Psf::Image> getPsfImage(
        geom::Point2D center,
        afw::table::SourceRecord & source,
        afw::image::Exposure<float> const & exposure
    ) const override;

    std::shared_ptr<afw::image::Mask<afw::image::MaskPixel>> getPsfMask(
        geom::Point2D center,
        afw::table::SourceRecord & source,
        afw::image::Exposure<float> const & exposure
    ) const override;

    afw::image::MaskPixel const getBadPixelMask(
        afw::image::Exposure<float> const & exposure
    ) const override;
};
}}}} // namespace lsst::meas::extensions::shapeHSM

#endif // !LSST_MEAS_EXTENSIONS_SHAPEHSM_HsmMomentsControl_h_INCLUDED
