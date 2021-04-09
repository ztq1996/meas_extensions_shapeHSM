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

#include "lsst/meas/extensions/shapeHSM/HsmMomentsControl.h"
#include "lsst/pex/config/python.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace extensions {
namespace shapeHSM {

PYBIND11_MODULE(hsmMomentsControl, mod) {
    py::module::import("lsst.afw.table");
    py::module::import("lsst.meas.base");

    /* Module level */
    py::class_<HsmMomentsAlgorithm, std::shared_ptr<HsmMomentsAlgorithm>, base::SimpleAlgorithm>
            clsHsmMomentsAlgorithm(mod, "HsmMomentsAlgorithm");

    py::class_<HsmSourceMomentsAlgorithm, std::shared_ptr<HsmSourceMomentsAlgorithm>, HsmMomentsAlgorithm>
            clsHsmSourceMomentsAlgorithm(mod, "HsmSourceMomentsAlgorithm");
    py::class_<HsmSourceMomentsControl> clsHsmSourceMomentsControl(mod, "HsmSourceMomentsControl");
    py::class_<HsmSourceMomentsRoundControl, HsmSourceMomentsControl>
            clsHsmSourceMomentsRoundControl(mod, "HsmSourceMomentsRoundControl");

    py::class_<HsmPsfMomentsAlgorithm, std::shared_ptr<HsmPsfMomentsAlgorithm>, HsmMomentsAlgorithm>
            clsHsmPsfMomentsAlgorithm(mod, "HsmPsfMomentsAlgorithm");
    py::class_<HsmPsfMomentsControl, std::shared_ptr<HsmPsfMomentsControl>>
            clsHsmPsfMomentsControl(mod, "HsmPsfMomentsControl");

    py::class_<HsmPsfMomentsDebiasedAlgorithm, std::shared_ptr<HsmPsfMomentsDebiasedAlgorithm>,
               HsmPsfMomentsAlgorithm, HsmMomentsAlgorithm>
            clsHsmPsfMomentsDebiasedAlgorithm(mod, "HsmPsfMomentsDebiasedAlgorithm");
    py::class_<HsmPsfMomentsDebiasedControl, std::shared_ptr<HsmPsfMomentsDebiasedControl>>
            clsHsmPsfMomentsDebiasedControl(mod, "HsmPsfMomentsDebiasedControl");

    /* Constructors */
    clsHsmSourceMomentsAlgorithm.def(
            py::init<HsmSourceMomentsAlgorithm::Control const &, std::string const &, afw::table::Schema &>(),
            "ctrl"_a, "name"_a, "schema"_a);
    clsHsmSourceMomentsControl.def(py::init<>());
    clsHsmSourceMomentsRoundControl.def(py::init<>());

    clsHsmPsfMomentsAlgorithm.def(
            py::init<HsmPsfMomentsAlgorithm::Control const &, std::string const &, afw::table::Schema &>(),
            "ctrl"_a, "name"_a, "schema"_a);
    clsHsmPsfMomentsControl.def(py::init<>());

    clsHsmPsfMomentsDebiasedAlgorithm.def(
            py::init<HsmPsfMomentsDebiasedAlgorithm::Control const &,
                     std::string const &,
                     afw::table::Schema &>(),
            "ctrl"_a, "name"_a, "schema"_a);
    clsHsmPsfMomentsDebiasedControl.def(py::init<>());

    /* Members */
    LSST_DECLARE_CONTROL_FIELD(clsHsmSourceMomentsControl, HsmSourceMomentsControl, badMaskPlanes);
    LSST_DECLARE_CONTROL_FIELD(clsHsmSourceMomentsControl, HsmSourceMomentsControl, roundMoments);
    LSST_DECLARE_CONTROL_FIELD(clsHsmSourceMomentsControl, HsmSourceMomentsControl, addFlux);

    LSST_DECLARE_CONTROL_FIELD(clsHsmPsfMomentsControl, HsmPsfMomentsControl, useSourceCentroidOffset);

    LSST_DECLARE_CONTROL_FIELD(
            clsHsmPsfMomentsDebiasedControl, HsmPsfMomentsDebiasedControl, useSourceCentroidOffset);
    LSST_DECLARE_CONTROL_FIELD(clsHsmPsfMomentsDebiasedControl, HsmPsfMomentsDebiasedControl, noiseSource);
    LSST_DECLARE_CONTROL_FIELD(clsHsmPsfMomentsDebiasedControl, HsmPsfMomentsDebiasedControl, seedOffset);
    LSST_DECLARE_CONTROL_FIELD(clsHsmPsfMomentsDebiasedControl, HsmPsfMomentsDebiasedControl, badMaskPlanes);

    clsHsmMomentsAlgorithm.def("fail", &HsmMomentsAlgorithm::fail);
}

}  // shapeHSM
}  // extensions
}  // meas
}  // lsst
