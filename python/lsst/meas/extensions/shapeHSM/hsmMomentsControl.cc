/*
 * LSST Data Management System
 *
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 * See the COPYRIGHT file
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
 * see <https://www.lsstcorp.org/LegalNotices/>.
 */
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "lsst/cpputils/python.h"

#include "lsst/meas/extensions/shapeHSM/HsmMomentsControl.h"
#include "lsst/pex/config/python.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace extensions {
namespace shapeHSM {

void wrapHsmMomentsControl(lsst::cpputils::python::WrapperCollection &wrappers) {

    using PyHsmMomentsAlgorithm = py::class_<HsmMomentsAlgorithm, std::shared_ptr<HsmMomentsAlgorithm>, lsst::meas::base::SimpleAlgorithm>;
    wrappers.wrapType(PyHsmMomentsAlgorithm(wrappers.module, "HsmMomentsAlgorithm"), [](auto &mod, auto &cls) {
        cls.def("fail", &HsmMomentsAlgorithm::fail);
    });

    using PyHsmSourceMomentsAlgorithm = py::class_<HsmSourceMomentsAlgorithm, std::shared_ptr<HsmSourceMomentsAlgorithm>, HsmMomentsAlgorithm>;
    wrappers.wrapType(PyHsmSourceMomentsAlgorithm(wrappers.module, "HsmSourceMomentsAlgorithm"), [](auto &mod, auto &cls) {
        cls.def(py::init<HsmSourceMomentsAlgorithm::Control const &, std::string const &, afw::table::Schema &>(),
                "ctrl"_a, "name"_a, "schema"_a);
    });

    using PyHsmSourceMomentsControl = py::class_<HsmSourceMomentsControl>;
    wrappers.wrapType(PyHsmSourceMomentsControl(wrappers.module, "HsmSourceMomentsControl"), [](auto &mod, auto &cls) {
        cls.def(py::init<>());
        LSST_DECLARE_CONTROL_FIELD(cls, HsmSourceMomentsControl, badMaskPlanes);
        LSST_DECLARE_CONTROL_FIELD(cls, HsmSourceMomentsControl, roundMoments);
        LSST_DECLARE_CONTROL_FIELD(cls, HsmSourceMomentsControl, addFlux);
    });

    using PyHsmSourceMomentsRoundAlgorithm = py::class_<HsmSourceMomentsRoundAlgorithm, std::shared_ptr<HsmSourceMomentsRoundAlgorithm>, HsmSourceMomentsAlgorithm, HsmMomentsAlgorithm>;
    wrappers.wrapType(PyHsmSourceMomentsRoundAlgorithm(wrappers.module, "HsmSourceMomentsRoundAlgorithm"), [](auto &mod, auto &cls) {
        cls.def(py::init<HsmSourceMomentsRoundAlgorithm::Control const &, std::string const &, afw::table::Schema &>(),
                "ctrl"_a, "name"_a, "schema"_a);
    });

    using PyHsmSourceMomentsRoundControl = py::class_<HsmSourceMomentsRoundControl, HsmSourceMomentsControl>;
    wrappers.wrapType(PyHsmSourceMomentsRoundControl(wrappers.module, "HsmSourceMomentsRoundControl"), [](auto &mod, auto &cls) {
        cls.def(py::init<>());
    });

    using PyHsmPsfMomentsAlgorithm = py::class_<HsmPsfMomentsAlgorithm, std::shared_ptr<HsmPsfMomentsAlgorithm>, HsmMomentsAlgorithm>;
    wrappers.wrapType(PyHsmPsfMomentsAlgorithm(wrappers.module, "HsmPsfMomentsAlgorithm"), [](auto &mod, auto &cls) {
        cls.def(py::init<HsmPsfMomentsAlgorithm::Control const &, std::string const &, afw::table::Schema &>(),
                "ctrl"_a, "name"_a, "schema"_a);
    });

    using PyHsmPsfMomentsControl = py::class_<HsmPsfMomentsControl, std::shared_ptr<HsmPsfMomentsControl>>;
    wrappers.wrapType(PyHsmPsfMomentsControl(wrappers.module, "HsmPsfMomentsControl"), [](auto &mod, auto &cls) {
        cls.def(py::init<>());
        LSST_DECLARE_CONTROL_FIELD(cls, HsmPsfMomentsControl, useSourceCentroidOffset);
    });

    using PyHsmPsfMomentsDebiasedAlgorithm = py::class_<HsmPsfMomentsDebiasedAlgorithm, std::shared_ptr<HsmPsfMomentsDebiasedAlgorithm>,
            HsmPsfMomentsAlgorithm, HsmMomentsAlgorithm>;
    wrappers.wrapType(PyHsmPsfMomentsDebiasedAlgorithm(wrappers.module, "HsmPsfMomentsDebiasedAlgorithm"), [](auto &mod, auto &cls) {
        cls.def(py::init<HsmPsfMomentsDebiasedAlgorithm::Control const &,
                        std::string const &,
                        afw::table::Schema &>(),
                "ctrl"_a, "name"_a, "schema"_a);
    });

    using PyHsmPsfMomentsDebiasedControl =  py::class_<HsmPsfMomentsDebiasedControl, std::shared_ptr<HsmPsfMomentsDebiasedControl>>;
    wrappers.wrapType(PyHsmPsfMomentsDebiasedControl(wrappers.module, "HsmPsfMomentsDebiasedControl"), [](auto &mod, auto &cls) {
        cls.def(py::init<>());
        LSST_DECLARE_CONTROL_FIELD(
                cls, HsmPsfMomentsDebiasedControl, useSourceCentroidOffset);
        LSST_DECLARE_CONTROL_FIELD(cls, HsmPsfMomentsDebiasedControl, noiseSource);
        LSST_DECLARE_CONTROL_FIELD(cls, HsmPsfMomentsDebiasedControl, seedOffset);
        LSST_DECLARE_CONTROL_FIELD(cls, HsmPsfMomentsDebiasedControl, badMaskPlanes);
    });

}

}  // shapeHSM
}  // extensions
}  // meas
}  // lsst
