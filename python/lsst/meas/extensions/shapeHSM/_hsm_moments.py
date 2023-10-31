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
import lsst.afw.image as afwImage
import lsst.afw.table as afwTable
import lsst.meas.base as measBase
import lsst.pex.config as pexConfig
from lsst.geom import Point2D, Point2I
from lsst.pex.exceptions import InvalidParameterError


class HsmMomentsConfig(measBase.SingleFramePluginConfig):
    roundMoments = pexConfig.Field[bool](doc="Use round weight function?", default=False)
    addFlux = pexConfig.Field[bool](doc="Store measured flux?", default=False)
    subtractCenter = pexConfig.Field[bool](doc="Subtract starting center from x/y outputs?", default=False)


class HsmMomentsPlugin(measBase.SingleFramePlugin):
    ConfigClass = HsmMomentsConfig

    # def __init__(self, config, name, schema, metadata, logName=None) -> None:
    def __init__(self, config, name, schema, metadata, logName=None):
        super().__init__(config, name, schema, logName)
        # Flag definitions.
        flagDefs = measBase.FlagDefinitionList()
        self.FAILURE = flagDefs.addFailureFlag("General failure flag, set if anything went wrong")
        self.NO_PIXELS = flagDefs.add("flag_no_pixels", "No pixels to measure")
        self.NOT_CONTAINED = flagDefs.add(
            "flag_not_contained", "Center not contained in footprint bounding box"
        )
        self.PARENT_SOURCE = flagDefs.add("flag_parent_source", "Parent source, ignored")
        self.GALSIM = flagDefs.add("flag_galsim", "GalSim failure")
        self.INVALID_PARAM = flagDefs.add("flag_invalid_param", "Invalid parameters")
        self.EDGE = flagDefs.add("flag_edge", "Variance undefined outside image edge")
        self.NO_PSF = flagDefs.add("flag_no_psf", "Exposure lacks PSF")

        # Flag handler.
        self.flagHandler = measBase.FlagHandler.addFields(schema, name, flagDefs)

        # Safe centroid extractor.
        self.centroidExtractor = measBase.SafeCentroidExtractor(schema, name)

    @classmethod
    def getExecutionOrder(cls):
        return cls.SHAPE_ORDER

    def calculate(
        self,
        record,
        *,
        image: galsim.Image,
        sigma: float,
        badpix: galsim.Image | None,
        centeroid: Point2D,
    ):
        # Convert centroid to galsim.PositionD
        guessCentroid = galsim.PositionD(centeroid.x, centeroid.y)

        try:
            shape = galsim.hsm.FindAdaptiveMom(
                image,
                weight=None,
                badpix=badpix,
                guess_sig=sigma,
                precision=1.0e-6,
                guess_centroid=guessCentroid,
                strict=True,
                round_moments=self.config.roundMoments,
                hsmparams=None,
            )
        except galsim.hsm.GalSimHSMError as error:
            raise measBase.MeasurementError(error, self.GALSIM.number)

        # now that we have shape, we can get the values below:
        determinantRadius = shape.moments_sigma
        centroidResult = shape.moments_centroid

        # subtract center should be implemented as well
        if self.config.subtractCenter:
            centroidResult.x -= centeroid.getX()
            centroidResult.y -= centeroid.getY()

        # Make a `lsst.geom.Point2D` object for the centroid.
        centroidResult = Point2D(centroidResult.x, centroidResult.y)

        # Populate the record with the results.
        record.set(self.centroidResultKey, centroidResult)

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
        except InvalidParameterError as error:
            raise measBase.MeasurementError(error, self.INVALID_PARAM.number)

        record.set(self.shapeKey, quad)

        if self.config.addFlux:
            record.set(self.fluxKey, shape.moments_amp)

    def fail(self, record, error=None):
        # docstring inherited.
        self.flagHandler.handleFailure(record)
        if error:
            centeroid = record.getCentroid()
            self.log.debug(
                "Failed to measure shape for %d at (%f, %f): %s",
                record.getId(),
                centeroid.getX(),
                centeroid.getY(),
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

    def __init__(self, config, name, schema, metadata, logName=None):
        super().__init__(config, name, schema, logName)
        self.centroidResultKey = afwTable.Point2DKey.addFields(schema, name, "HSM source centroid", "pixel")
        self.shapeKey = afwTable.QuadrupoleKey.addFields(
            schema, name, "HSM source moments", afwTable.CoordinateType.PIXEL
        )
        if config.addFlux:
            self.fluxKey = measBase.FluxResultKey.addFields(schema, name + "_Flux", "HSM source flux", "dn")

    def measure(self, record, exposure):
        """... in place"""
        # source in c++ is record in py. exposure the same

        # center = record.getCentroid()
        center = self.centroidExtractor(record, self.flagHandler)

        bbox = record.getFootprint().getBBox()
        if bbox.getArea() == 0:
            raise measBase.MeasurementError(self.NO_PIXELS.doc, self.NO_PIXELS.number)

        if not bbox.contains(Point2I(center)):
            raise measBase.MeasurementError(self.NOT_CONTAINED.doc, self.NOT_CONTAINED.number)

        psfSigma = exposure.getPsf().computeShape(center).getTraceRadius()

        # galsim
        imageArray = exposure[bbox].getImage().array
        xmin, xmax = bbox.getMinX(), bbox.getMaxX()
        ymin, ymax = bbox.getMinY(), bbox.getMaxY()
        bounds = galsim.bounds.BoundsI(xmin, xmax, ymin, ymax)
        image = galsim.Image(imageArray, bounds=bounds, copy=False)
        # GalSim's HSM uses the FITS convention of 1,1 for the lower-leftcorner
        # account for difference in origin while converting center to galsim
        # centerid (guessed)
        # make sure bbox is centered properly

        # Get the `lsst.meas.base` mask for bad pixels.
        badpix = exposure[bbox].mask.array.copy()
        bitValue = exposure.mask.getPlaneBitMask(self.config.badMaskPlanes)
        badpix &= 2**bitValue  # let's go back to this!

        # Convert to `galsim` mask (Note that galsim.Image will match whatever
        # dtype the input array is (here int32).
        badpix = galsim.Image(badpix, bounds=bounds)

        self.calculate(
            record,
            image=image,
            sigma=2.5 * psfSigma,
            badpix=badpix,
            centeroid=center,
        )

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


@measBase.register("ext_shapeHSM_HsmSourceMomentsRound")
class HsmSourceMomentsRoundPlugin(HsmSourceMomentsPlugin):
    ConfigClass = HsmSourceMomentsRoundConfig


class HsmPsfMomentsConfig(HsmMomentsConfig):
    useSourceCentroidOffset = pexConfig.Field[bool](doc="Use source centroid offset?", default=False)


@measBase.register("ext_shapeHSM_HsmPsfMoments")
class HsmPsfMomentsPlugin(HsmMomentsPlugin):
    ConfigClass = HsmPsfMomentsConfig

    def __init__(self, config, name, schema, metadata, logName=None):
        super().__init__(config, name, schema, logName)
        self.centroidResultKey = afwTable.Point2DKey.addFields(
            schema, name, "HSM PSF centroid", "pixel"
        )  # not sure if we report this but let's keep it here for now
        self.shapeKey = afwTable.QuadrupoleKey.addFields(
            schema, name, "HSM PSF moments", afwTable.CoordinateType.PIXEL
        )
        if config.addFlux:
            self.fluxKey = measBase.FluxResultKey.addFields(schema, name + "_Flux", "HSM PSF flux", "dn")

    def measure(self, record, exposure):
        """... in place PSF"""
        # source in c++ is record in py. exposure the same

        # center = record.getCentroid()
        center = self.centroidExtractor(record, self.flagHandler)

        psf = exposure.getPsf()
        if not psf:
            raise measBase.MeasurementError(self.NO_PSF.doc, self.NO_PSF.number)

        if self.config.useSourceCentroidOffset:
            psfImage = psf.computeImage(center)
        else:
            psfImage = psf.computeKernelImage(center)
            psfImage.setXY0(psf.computeImageBBox(center).getMin())

        # psfMask =  ... don't really need it here

        psfSigma = psf.computeShape(center).getTraceRadius()

        # galsim
        imageArray = psfImage.array

        # Get the bounding box of the PSF image in the parent image.
        bbox = psfImage.getBBox(afwImage.PARENT)

        xmin, xmax = bbox.getMinX(), bbox.getMaxX()
        ymin, ymax = bbox.getMinY(), bbox.getMaxY()
        bounds = galsim.bounds.BoundsI(xmin, xmax, ymin, ymax)
        image = galsim.Image(imageArray, bounds=bounds, copy=False)

        # No psfMask, no badpix calculation for PSF.
        badpix = None

        if self.config.useSourceCentroidOffset:
            centroid = center
        else:
            centroid = Point2D(bbox.getMin() + bbox.getDimensions() / 2)

        self.calculate(
            record,
            image=image,
            sigma=psfSigma,
            badpix=badpix,
            centeroid=centroid,
        )


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


@measBase.register("ext_shapeHSM_HsmPsfMomentsDebiased")
class HsmPsfMomentsDebiasedPlugin(HsmPsfMomentsPlugin):
    ConfigClass = HsmPsfMomentsDebiasedConfig

    @classmethod
    def getExecutionOrder(cls):
        return cls.FLUX_ORDER + 1
