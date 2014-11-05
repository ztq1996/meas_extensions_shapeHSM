// -*- lsst-c++ -*-
#ifndef LSST_MEAS_EXTENSIONS_SHAPEHSM_HsmMomentsControl_h_INCLUDED
#define LSST_MEAS_EXTENSIONS_SHAPEHSM_HsmMomentsControl_h_INCLUDED

#include "lsst/meas/algorithms/ShapeControl.h"

namespace lsst { namespace meas { namespace extensions { namespace shapeHSM {

class HsmSourceMomentsControl : public algorithms::ShapeControl {
public:
    LSST_CONTROL_FIELD(badMaskPlanes, std::vector<std::string>, "Mask planes used to reject bad pixels.");

    HsmSourceMomentsControl() : ShapeControl("shape.hsm.moments") {}

    HsmSourceMomentsControl(HsmSourceMomentsControl const & other) :
        algorithms::ShapeControl(other), badMaskPlanes(other.badMaskPlanes) {}

    HsmSourceMomentsControl(std::string const & name) : algorithms::ShapeControl(name) {
        badMaskPlanes.push_back("BAD");
        badMaskPlanes.push_back("SAT");
        badMaskPlanes.push_back("INTRP");
    }
private:
    virtual PTR(algorithms::AlgorithmControl) _clone() const;
    virtual PTR(algorithms::Algorithm) _makeAlgorithm(
        afw::table::Schema & schema, PTR(daf::base::PropertyList) const & metadata
    ) const;
};


class HsmPsfMomentsControl : public algorithms::ShapeControl {
public:
    HsmPsfMomentsControl() : ShapeControl("shape.hsm.psfMoments") {}
    HsmPsfMomentsControl(std::string const & name) : algorithms::ShapeControl(name) {}
    HsmPsfMomentsControl(HsmPsfMomentsControl const & other) : algorithms::ShapeControl(other) {}
private:
    virtual PTR(algorithms::AlgorithmControl) _clone() const;
    virtual PTR(algorithms::Algorithm) _makeAlgorithm(
        afw::table::Schema & schema, PTR(daf::base::PropertyList) const & metadata
    ) const;
};


}}}} // namespace lsst::meas::extensions::shapeHSM

#endif // !LSST_MEAS_EXTENSIONS_SHAPEHSM_HsmMomentsControl_h_INCLUDED
