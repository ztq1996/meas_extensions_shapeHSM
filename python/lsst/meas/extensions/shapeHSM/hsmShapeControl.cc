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

#include "lsst/meas/extensions/shapeHSM/HsmShapeControl.h"
#include "lsst/pex/config/python.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace extensions {
namespace shapeHSM {

PYBIND11_MODULE(hsmShapeControl, mod) {
    py::module::import("lsst.afw.table");
    py::module::import("lsst.meas.base");

    /* Module level */
    py::class_<HsmShapeAlgorithm, base::SimpleAlgorithm>
            clsHsmShapeAlgorithm(mod, "HsmShapeAlgorithm");
    py::class_<HsmShapeControl> clsHsmShapeControl(mod, "HsmShapeControl");

    py::class_<HsmShapeBjAlgorithm, HsmShapeAlgorithm>
            clsHsmShapeBjAlgorithm(mod, "HsmShapeBjAlgorithm");
    py::class_<HsmShapeBjControl, HsmShapeControl> clsHsmShapeBjControl(mod, "HsmShapeBjControl");

    py::class_<HsmShapeLinearAlgorithm, HsmShapeAlgorithm>
            clsHsmShapeLinearAlgorithm(mod, "HsmShapeLinearAlgorithm");
    py::class_<HsmShapeLinearControl, HsmShapeControl> clsHsmShapeLinearControl(mod, "HsmShapeLinearControl");

    py::class_<HsmShapeKsbAlgorithm, HsmShapeAlgorithm>
            clsHsmShapeKsbAlgorithm(mod, "HsmShapeKsbAlgorithm");
    py::class_<HsmShapeKsbControl, HsmShapeControl> clsHsmShapeKsbControl(mod, "HsmShapeKsbControl");

    py::class_<HsmShapeRegaussAlgorithm, HsmShapeAlgorithm>
            clsHsmShapeRegaussAlgorithm(mod, "HsmShapeRegaussAlgorithm");
    py::class_<HsmShapeRegaussControl, HsmShapeControl> clsHsmShapeRegaussControl(mod,
                                                                                  "HsmShapeRegaussControl");
    /* Constructors */
    clsHsmShapeBjAlgorithm.def(
            py::init<HsmShapeBjAlgorithm::Control const &, std::string const &, afw::table::Schema &>(),
            "ctrl"_a, "name"_a, "schema"_a);
    clsHsmShapeBjControl.def(py::init<>());

    clsHsmShapeLinearAlgorithm.def(
            py::init<HsmShapeLinearAlgorithm::Control const &, std::string const &, afw::table::Schema &>(),
            "ctrl"_a, "name"_a, "schema"_a);
    clsHsmShapeLinearControl.def(py::init<>());

    clsHsmShapeKsbAlgorithm.def(
            py::init<HsmShapeKsbAlgorithm::Control const &, std::string const &, afw::table::Schema &>(),
            "ctrl"_a, "name"_a, "schema"_a);
    clsHsmShapeKsbControl.def(py::init<>());

    clsHsmShapeRegaussAlgorithm.def(
            py::init<HsmShapeRegaussAlgorithm::Control const &, std::string const &, afw::table::Schema &>(),
            "ctrl"_a, "name"_a, "schema"_a);
    clsHsmShapeRegaussControl.def(py::init<>());

    /* Members */
    LSST_DECLARE_CONTROL_FIELD(clsHsmShapeControl, HsmShapeControl, badMaskPlanes);
    LSST_DECLARE_CONTROL_FIELD(clsHsmShapeControl, HsmShapeControl, deblendNChild);

    clsHsmShapeAlgorithm.def("measure", &HsmShapeAlgorithm::measure);
    clsHsmShapeAlgorithm.def("fail", &HsmShapeAlgorithm::fail);
}

}  // shapeHSM
}  // extensions
}  // meas
}  // lsst
