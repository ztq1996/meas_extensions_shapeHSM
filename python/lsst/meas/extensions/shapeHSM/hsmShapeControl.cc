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
#include "lsst/cpputils/python.h"

#include "lsst/meas/extensions/shapeHSM/HsmShapeControl.h"
#include "lsst/pex/config/python.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace extensions {
namespace shapeHSM {
void wrapHsmShapeControl(lsst::cpputils::python::WrapperCollection &wrappers) {
    wrappers.addInheritanceDependency("lsst.meas.base");
    wrappers.addSignatureDependency("lsst.afw.table");

    /* Module level */
    using PyHsmShapeAlgorithm =  py::class_<HsmShapeAlgorithm, std::shared_ptr<HsmShapeAlgorithm>, base::SimpleAlgorithm>;
    wrappers.wrapType(PyHsmShapeAlgorithm(wrappers.module, "HsmShapeAlgorithm"), [](auto &mod, auto &cls) {
        cls.def("measure", &HsmShapeAlgorithm::measure);
        cls.def("fail", &HsmShapeAlgorithm::fail);
    });

    using PyHsmShapeControl =  py::class_<HsmShapeControl>;
    wrappers.wrapType(PyHsmShapeControl(wrappers.module, "HsmShapeControl"), [](auto &mod, auto &cls) {
        LSST_DECLARE_CONTROL_FIELD(cls, HsmShapeControl, badMaskPlanes);
        LSST_DECLARE_CONTROL_FIELD(cls, HsmShapeControl, deblendNChild);
    });

    using PyHsmShapeBjAlgorithm = py::class_<HsmShapeBjAlgorithm, std::shared_ptr<HsmShapeBjAlgorithm>, HsmShapeAlgorithm>;
    wrappers.wrapType(PyHsmShapeBjAlgorithm(wrappers.module, "HsmShapeBjAlgorithm"), [](auto &mod, auto &cls) {
        cls.def(py::init<HsmShapeBjAlgorithm::Control const &, std::string const &, afw::table::Schema &>(),
                "ctrl"_a, "name"_a, "schema"_a);
    });

    using PyHsmShapeBjControl = py::class_<HsmShapeBjControl, HsmShapeControl>;
    wrappers.wrapType(PyHsmShapeBjControl(wrappers.module, "HsmShapeBjControl"), [](auto &mod, auto &cls) {
        cls.def(py::init<>());
    });

    using PyHsmShapeLinearAlgorithm =  py::class_<HsmShapeLinearAlgorithm, std::shared_ptr<HsmShapeLinearAlgorithm>, HsmShapeAlgorithm>;
    wrappers.wrapType(PyHsmShapeLinearAlgorithm(wrappers.module, "HsmShapeLinearAlgorithm"), [](auto &mod, auto &cls) {
        cls.def(py::init<HsmShapeLinearAlgorithm::Control const &, std::string const &, afw::table::Schema &>(),
                "ctrl"_a, "name"_a, "schema"_a);
    });

    using PyHsmShapeLinearControl = py::class_<HsmShapeLinearControl, HsmShapeControl>;
    wrappers.wrapType(PyHsmShapeLinearControl(wrappers.module, "HsmShapeLinearControl"), [](auto &mod, auto &cls) {
        cls.def(py::init<>());
    });

    using PyHsmShapeKsbAlgorithm =
            py::class_<HsmShapeKsbAlgorithm, std::shared_ptr<HsmShapeKsbAlgorithm>, HsmShapeAlgorithm>;
    wrappers.wrapType(PyHsmShapeKsbAlgorithm(wrappers.module, "HsmShapeKsbAlgorithm"), [](auto &mod, auto &cls) {
        cls.def(
                py::init<HsmShapeKsbAlgorithm::Control const &, std::string const &, afw::table::Schema &>(),
                "ctrl"_a, "name"_a, "schema"_a);
    });

    using PyHsmShapeKsbControl = py::class_<HsmShapeKsbControl, HsmShapeControl>;
    wrappers.wrapType(PyHsmShapeKsbControl(wrappers.module, "HsmShapeKsbControl"), [](auto &mod, auto &cls) {
        cls.def(py::init<>());
    });

    using PyHsmShapeRegaussAlgorithm = py::class_<HsmShapeRegaussAlgorithm, std::shared_ptr<HsmShapeRegaussAlgorithm>, HsmShapeAlgorithm>;
    wrappers.wrapType(PyHsmShapeRegaussAlgorithm(wrappers.module, "HsmShapeRegaussAlgorithm"), [](auto &mod, auto &cls) {
        cls.def(
                py::init<HsmShapeRegaussAlgorithm::Control const &, std::string const &, afw::table::Schema &>(),
                "ctrl"_a, "name"_a, "schema"_a);
    });

    using PyHsmShapeRegaussControl = py::class_<HsmShapeRegaussControl, HsmShapeControl>;
    wrappers.wrapType(PyHsmShapeRegaussControl(wrappers.module, "HsmShapeRegaussControl"), [](auto &mod, auto &cls) {
        cls.def(py::init<>());
    });

}

}  // shapeHSM
}  // extensions
}  // meas
}  // lsst
