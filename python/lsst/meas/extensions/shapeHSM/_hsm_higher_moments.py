# This file is part of meas_extensions_shapeHSM.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

__all__ = []

import logging
import lsst.meas.base as measBase
from lsst.pex.confg import Field


class HigherOrderMomentsConfig(measBase.SingleFramePluginConfig):
    psfSigma = Field[float](
        doc="PSF sigma in pixels",
        default=5.0,
    )

    order = Field[int](
        doc="Order of moments to compute",
        default=4,
    )


@measBase.register("ext_shapeHSM_HigherOrderMoments")
class HigherOrderMomentsPlugin(measBase.SingleFramePlugin):
    ConfigClass = HigherOrderMomentsConfig

    @classmethod
    def getExecutionOrder(cls):
        return cls.SHAPE_ORDER

    def __init__(self, config, name, schema, metadata, logName=None):
        super().__init__(config, name, schema, metadata, logName=logName)

        # Add the output columns to the schema
        # Attn: TQ
        schema.addField(schema.join(name, "blah"), type=float, doc="blah")

        # Define flags for possible issues that might arise during measurement.
        flagDefs = measBase.FlagDefinitionList()
        self.FAILURE = flagDefs.addFailureFlag("General failure flag, set if anything went wrong")

        # Embed the flag definitions in the schema using a flag handler.
        self.flagHandler = measBase.FlagHandler.addFields(schema, name, flagDefs)

        # Utilize a safe centroid extractor that uses the detection footprint
        # as a fallback if necessary.
        self.centroidExtractor = measBase.SafeCentroidExtractor(schema, name)
        self.log = logging.getLogger(self.logName)

    def fail(self, record, error=None):
        # Docstring inherited.
        self.flagHandler.handleFailure(record)

    def measure(self, record, exposure):
        center = self.centroidExtractor(record, self.flagHandler)
