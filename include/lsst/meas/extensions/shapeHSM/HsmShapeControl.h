// -*- lsst-c++ -*-
#ifndef LSST_MEAS_EXTENSIONS_SHAPEHSM_HsmShapeControl_h_INCLUDED
#define LSST_MEAS_EXTENSIONS_SHAPEHSM_HsmShapeControl_h_INCLUDED

#include "lsst/meas/algorithms/Algorithm.h"

namespace lsst { namespace meas { namespace extensions { namespace shapeHSM {

class HsmShapeControl : public algorithms::AlgorithmControl {
public:
    LSST_CONTROL_FIELD(badMaskPlanes, std::vector<std::string>, "Mask planes used to reject bad pixels.");
    LSST_CONTROL_FIELD(deblendNChild, std::string, "Field name for number of deblend children");
protected:

    HsmShapeControl(HsmShapeControl const & other) :
        AlgorithmControl(other), badMaskPlanes(other.badMaskPlanes), deblendNChild(other.deblendNChild) {}

    HsmShapeControl(std::string const & name) : AlgorithmControl(name, 1.0), deblendNChild("") {
        badMaskPlanes.push_back("BAD");
        badMaskPlanes.push_back("SAT");
        badMaskPlanes.push_back("INTRP");
    }
};

class HsmShapeBjControl : public HsmShapeControl {
public:
    HsmShapeBjControl() : HsmShapeControl("shape.hsm.bj") {}
private:
    virtual PTR(algorithms::AlgorithmControl) _clone() const;
    virtual PTR(algorithms::Algorithm) _makeAlgorithm(
        afw::table::Schema & schema, PTR(daf::base::PropertyList) const & metadata
    ) const;
};

class HsmShapeLinearControl : public HsmShapeControl {
public:
    HsmShapeLinearControl() : HsmShapeControl("shape.hsm.linear") {}
private:
    virtual PTR(algorithms::AlgorithmControl) _clone() const;
    virtual PTR(algorithms::Algorithm) _makeAlgorithm(
        afw::table::Schema & schema, PTR(daf::base::PropertyList) const & metadata
    ) const;
};

class HsmShapeKsbControl : public HsmShapeControl {
public:
    HsmShapeKsbControl() : HsmShapeControl("shape.hsm.ksb") {}
private:
    virtual PTR(algorithms::AlgorithmControl) _clone() const;
    virtual PTR(algorithms::Algorithm) _makeAlgorithm(
        afw::table::Schema & schema, PTR(daf::base::PropertyList) const & metadata
    ) const;
};

class HsmShapeRegaussControl : public HsmShapeControl {
public:
    HsmShapeRegaussControl() : HsmShapeControl("shape.hsm.regauss") {}
private:
    virtual PTR(algorithms::AlgorithmControl) _clone() const;
    virtual PTR(algorithms::Algorithm) _makeAlgorithm(
        afw::table::Schema & schema, PTR(daf::base::PropertyList) const & metadata
    ) const;
};

class HsmShapeShapeletControl : public HsmShapeControl {
public:
    LSST_CONTROL_FIELD(maxOrderPsf, int, "Maximum shapelet order for PSF");
    LSST_CONTROL_FIELD(maxOrderGalaxy, int, "Maximum shapelet order for galaxy");

    HsmShapeShapeletControl() : HsmShapeControl("shape.hsm.shapelet"), maxOrderPsf(8), maxOrderGalaxy(8) {}
private:
    virtual PTR(algorithms::AlgorithmControl) _clone() const;
    virtual PTR(algorithms::Algorithm) _makeAlgorithm(
        afw::table::Schema & schema, PTR(daf::base::PropertyList) const & metadata
    ) const;
};

}}}} // namespace lsst::meas::extensions::shapeHSM

#endif // !LSST_MEAS_EXTENSIONS_SHAPEHSM_HsmShapeControl_h_INCLUDED
