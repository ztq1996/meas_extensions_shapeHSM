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

import galsim
import lsst.afw.geom as afwGeom
import lsst.afw.table as afwTable
import lsst.meas.base as measBase
import lsst.pex.config as pexConfig
from lsst.meas.base import (
    FlagDefinitionList,
    FlagHandler,
    MeasurementError,
    SingleFramePlugin,
    SingleFramePluginConfig,
)
from lsst.pex.exceptions import InvalidParameterError
from lsst.geom import Point2I, Point2D


class HsmMomentsConfig(SingleFramePluginConfig):
    roundMoments = pexConfig.Field[bool](doc="Use round weight function?", default=False)
    addFlux = pexConfig.Field[bool](doc="Store measured flux?", default=False)
    subtractCenter = pexConfig.Field[bool](doc="Subtract starting center from x/y outputs?", default=False)


class HsmMomentsPlugin(SingleFramePlugin):
    ConfigClass = HsmMomentsConfig

    # def __init__(self, config, name, schema, metadata, logName=None) -> None:
    def __init__(self, config, name, schema, metadata, logName=None):
        super().__init__(config, name, schema, logName)

        self.centroidResultKey = afwTable.Point2DKey.addFields(
            schema, name, "The centroid measured by HSM", "pixel"
        )
        self.shapeKey = afwTable.QuadrupoleKey.addFields(schema, name, doc="HSM source moments.")
        if config.addFlux:
            self.fluxKey = measBase.FluxResultKey.addFields(schema, name + "_Flux", "HSM flux", "dn")

        # Flag definitions.
        flagDefs = FlagDefinitionList()
        self.FAILURE = flagDefs.addFailureFlag("General failure flag, set if anything went wrong")
        self.NO_PIXELS = flagDefs.add("flag_no_pixels", "No pixels to measure")
        self.NOT_CONTAINED = flagDefs.add(
            "flag_not_contained", "Center not contained in footprint bounding box"
        )
        self.PARENT_SOURCE = flagDefs.add("flag_parent_source", "Parent source, ignored")
        self.GALSIM = flagDefs.add("flag_galsim", "GalSim failure")
        self.EDGE = flagDefs.add("flag_edge", "Variance undefined outside image edge")

        # Flag handler.
        self.flagHandler = FlagHandler.addFields(schema, name, flagDefs)

        # Safe centroid extractor.
        self.centroidExtractor = measBase.SafeCentroidExtractor(schema, name)

    @classmethod
    def getExecutionOrder(cls):
        return cls.SHAPE_ORDER

    def fail(self, record, error=None):
        # docstring inherited.
        self.flagHandler.handleFailure(record)
        if error:
            center = record.getCentroid()
            self.log.debug(
                "Failed to measure shape for %d at (%f, %f): %s",
                record.getId(),
                center.getX(),
                center.getY(),
                error,
            )


class HsmSourceMomentsConfig(HsmMomentsConfig):
    # later we can put some doscstring
    badMaskPlanes = pexConfig.ListField[str](
        doc="Mask planes used to reject bad pixels.", default=["BAD", "SAT"]
    )


@measBase.register("ext_shapeHSM_HsmSourceMoments")
class HsmSourceMomentsPlugin(HsmMomentsPlugin):
    ConfigClass = HsmSourceMomentsConfig

    # def __init__(self, config, name, schema, metadata, logName=None):
    #     super().__init__(config, name, schema, logName)
    #     self.centroidResultKey = afwTable.Point2DKey.addFields(
    #         schema, name, "The centroid measured by HSM", "pixel"
    #     )
    #     self.shapeKey = afwTable.QuadrupoleKey.addFields(schema, name, doc="HSM source moments.")
    #     if config.addFlux:
    #         self.fluxKey = measBase.FluxResultKey.addFields(schema, name + "_Flux", "HSM flux", "dn")

    def measure(self, record, exposure):
        """... in place"""
        # source in c++ is record in py. exposure the same

        # center = record.getCentroid()
        center = self.centroidExtractor(record, self.flagHandler)

        bbox = record.getFootprint().getBBox()
        if bbox.getArea() == 0:
            raise MeasurementError(self.NO_PIXELS.doc, self.NO_PIXELS.number)

        if not bbox.contains(Point2I(center)):
            raise MeasurementError(self.NOT_CONTAINED.doc, self.NOT_CONTAINED.number)

        psfSigma = exposure.getPsf().computeShape(center).getTraceRadius()

        # galsim
        image_array = exposure[bbox].getImage().array
        xmin, xmax = bbox.getMinX() + 1, bbox.getMaxX() + 1
        ymin, ymax = bbox.getMinY() + 1, bbox.getMaxY() + 1
        bounds = galsim.bounds.BoundsI(xmin, xmax, ymin, ymax)
        image = galsim.Image(image_array, bounds=bounds, copy=False)
        # GalSim's HSM uses the FITS convention of 1,1 for the lower-leftcorner
        # account for difference in origin while converting center to galsim
        # centerid (guessed)
        # make sure bbox is centered properly

        bounds = galsim.bounds.BoundsI(xmin, xmax, ymin, ymax)

        # Get the `lsst.meas.base` mask for bad pixels.
        badpix = exposure[bbox].mask.array.copy()
        bitValue = exposure.mask.getPlaneBitMask(self.config.badMaskPlanes)
        badpix &= bitValue

        # Convert to `galsim` mask (Note that galsim.Image will match whatever
        # dtype the input array is (here int32).
        badpix = galsim.Image(badpix, bounds=bounds)

        # Convert centroid to galsim.PositionD
        guess_centroid = galsim.PositionD(center.x, center.y)
        breakpoint()

        try:
            shape = galsim.hsm.FindAdaptiveMom(
                image,
                weight=None,
                badpix=badpix,
                guess_sig=2.5 * psfSigma,
                precision=1.0e-6,
                guess_centroid=guess_centroid,
                strict=True,
                round_moments=self.config.roundMoments,
                hsmparams=None,
            )
        except galsim.hsm.GalSimHSMError as e:
            raise MeasurementError(e, self.GALSIM.number)

        # now that we have shape, we can get the values below:
        determinantRadius = shape.moments_sigma
        centroidResult = shape.moments_centroid

        # subtract center should be implemented as well
        if self.config.subtractCenter:
            centroidResult.x -= center.getX()
            centroidResult.y -= center.getY()

        # Convert the centroid to a `lsst.geom.Point2D`.
        centroidResult = Point2D(centroidResult.x, centroidResult.y)

        # centroidResult
        record.set(self.centroidResultKey, centroidResult)

        # record["ext_shapeHSM_HsmSourceMoments_x"]

        # Convert galsim measurements to lsst measurements.
        try:
            ellipse = afwGeom.ellipses.SeparableDistortionDeterminantRadius(
                e1=shape.observed_shape.e1,
                e2=shape.observed_shape.e2,
                radius=determinantRadius,
                normalize=True,  # Fail if |e|>1.
            )
            quad = afwGeom.ellipses.Quadrupole(
                ellipse,
            )
        except InvalidParameterError as e:
            raise MeasurementError(e)

        record.set(self.shapeKey, quad)

        if self.config.addFlux:
            record.set(self.fluxKey, shape.moments_amp)

        # TODO: calculate errors in shape, centroid?
        # galsim does not give us any errors!
        # centroidErr = record.getCentroidErr()


class HsmSourceMomentsRoundConfig(HsmSourceMomentsConfig):
    def setDefaults(self):
        super().setDefaults()
        self.roundMoments = True

    def validate(self):
        if not self.roundMoments:
            raise pexConfig.FieldValidationError(
                self.roundMoments, self, "roundMoments should be set to `True`."
            )
        super().validate()


class HsmSourceMomentsRoundPlugin(HsmSourceMomentsPlugin):
    ConfigClass = HsmSourceMomentsRoundConfig


class HsmPsfMomentsConfig(HsmMomentsConfig):
    useSourceCentroidOffset = pexConfig.Field[bool](doc="Use source centroid offset?", default=False)


class HsmPsfMomentsDebiasedConfig(HsmPsfMomentsConfig):
    # HsmPsfMomentsConfig does not have badMaskPlanes so we add it here
    badMaskPlanes = pexConfig.ListField[str](
        doc="Mask planes used to reject bad pixels.", default=["BAD", "SAT"]
    )
    noiseSource = pexConfig.ChoiceField[str](
        doc="...",
        allowed={
            "meta": "variance = the 'BGMEAN' metadata entry",
            "variance": "variance = the image's variance plane",
        },
        default="variance",
    )
    seedOffset = pexConfig.Field[int](doc="Seed offset for random number generator.", default=0)

    def setDefaults(self):
        super().setDefaults()
        self.useSourceCentroidOffset = True
