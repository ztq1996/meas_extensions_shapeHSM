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

import logging

import galsim
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.table as afwTable
import lsst.meas.base as measBase
import lsst.pex.config as pexConfig
import numpy as np
from lsst.geom import Box2I, Point2D, Point2I
from lsst.pex.exceptions import InvalidParameterError, NotFoundError

__all__ = [
    "HsmSourceMomentsConfig",
    "HsmSourceMomentsPlugin",
    "HsmSourceMomentsRoundConfig",
    "HsmSourceMomentsRoundPlugin",
    "HsmPsfMomentsConfig",
    "HsmPsfMomentsPlugin",
    "HsmPsfMomentsDebiasedConfig",
    "HsmPsfMomentsDebiasedPlugin",
]


def _quickGalsimImageView(array, bounds):
    match array.dtype:
        case np.float64:
            gcls = galsim._galsim.ImageViewD
        case np.float32:
            gcls = galsim._galsim.ImageViewF
        case np.int32:
            gcls = galsim._galsim.ImageViewI

    return gcls(
        array.__array_interface__["data"][0],
        array.strides[1]//array.itemsize,
        array.strides[0]//array.itemsize,
        bounds._b,
    )


class HsmMomentsConfig(measBase.SingleFramePluginConfig):
    """Base configuration for HSM adaptive moments measurement."""

    roundMoments = pexConfig.Field[bool](doc="Use round weight function?", default=False)
    addFlux = pexConfig.Field[bool](doc="Store measured flux?", default=False)
    subtractCenter = pexConfig.Field[bool](doc="Subtract starting center from x/y outputs?", default=False)


class HsmMomentsPlugin(measBase.SingleFramePlugin):
    """Base plugin for HSM adaptive moments measurement."""

    ConfigClass = HsmMomentsConfig

    def __init__(self, config, name, schema, metadata, logName=None):
        if logName is None:
            logName = __name__
        super().__init__(config, name, schema, metadata, logName=logName)

        # Define flags for possible issues that might arise during measurement.
        flagDefs = measBase.FlagDefinitionList()
        self.FAILURE = flagDefs.addFailureFlag("General failure flag, set if anything went wrong")
        self.NO_PIXELS = flagDefs.add("flag_no_pixels", "No pixels to measure")
        self.NOT_CONTAINED = flagDefs.add(
            "flag_not_contained", "Center not contained in footprint bounding box"
        )
        self.PARENT_SOURCE = flagDefs.add("flag_parent_source", "Parent source, ignored")
        self.GALSIM = flagDefs.add("flag_galsim", "GalSim failure")
        self.INVALID_PARAM = flagDefs.add("flag_invalid_param", "Invalid combination of moments")
        self.EDGE = flagDefs.add("flag_edge", "Variance undefined outside image edge")
        self.NO_PSF = flagDefs.add("flag_no_psf", "Exposure lacks PSF")

        # Embed the flag definitions in the schema using a flag handler.
        self.flagHandler = measBase.FlagHandler.addFields(schema, name, flagDefs)

        # Utilize a safe centroid extractor that uses the detection footprint
        # as a fallback if necessary.
        self.centroidExtractor = measBase.SafeCentroidExtractor(schema, name)
        self.log = logging.getLogger(self.logName)

    @classmethod
    def getExecutionOrder(cls):
        return cls.SHAPE_ORDER

    def _calculate(
        self,
        record: afwTable.SourceRecord,
        *,
        image: galsim._galsim.ImageViewD | galsim._galsim.ImageViewF,
        weight_image: galsim._galsim.ImageViewI,
        centroid: Point2D,
        sigma: float = 5.0,
        precision: float = 1.0e-6,
    ) -> None:
        """
        Calculate adaptive moments using GalSim's HSM and modify the record in
        place.

        Parameters
        ----------
        record : `~lsst.afw.table.SourceRecord`
            Record to store measurements.
        image : `~galsim._galsim.ImageViewF` or `~galsim._galsim.ImageViewD`
            Image on which to perform measurements.
        weight_image : `~galsim._galsim.ImageViewI`
            The combined badpix/weight image for input to galsim HSM code.
        centroid : `~lsst.geom.Point2D`
            Centroid guess for HSM adaptive moments.
        sigma : `float`, optional
            Estimate of object's Gaussian sigma in pixels. Default is 5.0.
        precision : `float`, optional
            Precision for HSM adaptive moments. Default is 1.0e-6.

        Raises
        ------
        MeasurementError
            Raised for errors in measurement.
        """
        # Convert centroid to GalSim's PositionD type.
        guessCentroid = galsim.PositionD(centroid.x, centroid.y)
        try:
            # Attempt to compute HSM moments.

            # Use galsim c++/python interface directly.
            shape = galsim.hsm.ShapeData()
            hsmparams = galsim.hsm.HSMParams.default

            galsim._galsim.FindAdaptiveMomView(
                shape._data,
                image,
                weight_image,
                float(sigma),
                float(precision),
                guessCentroid._p,
                bool(self.config.roundMoments),
                hsmparams._hsmp,
            )

        except RuntimeError as error:
            raise measBase.MeasurementError(str(error), self.GALSIM.number)

        # Retrieve computed moments sigma and centroid.
        determinantRadius = shape.moments_sigma
        centroidResult = shape.moments_centroid

        # Subtract center if required by configuration.
        if self.config.subtractCenter:
            centroidResult.x -= centroid.getX()
            centroidResult.y -= centroid.getY()

        # Convert GalSim's `galsim.PositionD` to `lsst.geom.Point2D`.
        centroidResult = Point2D(centroidResult.x, centroidResult.y)

        # Populate the record with the centroid results.
        record.set(self.centroidResultKey, centroidResult)

        # Convert GalSim measurements to lsst measurements.
        try:
            # Create an ellipse for the shape.
            ellipse = afwGeom.ellipses.SeparableDistortionDeterminantRadius(
                e1=shape.observed_shape.e1,
                e2=shape.observed_shape.e2,
                radius=determinantRadius,
                normalize=True,  # Fail if |e|>1.
            )
            # Get the quadrupole moments from the ellipse.
            quad = afwGeom.ellipses.Quadrupole(ellipse)
        except InvalidParameterError as error:
            raise measBase.MeasurementError(error, self.INVALID_PARAM.number)

        # Store the quadrupole moments in the record.
        record.set(self.shapeKey, quad)

        # Store the flux if required by configuration.
        if self.config.addFlux:
            record.set(self.fluxKey, shape.moments_amp)

    def fail(self, record, error=None):
        # Docstring inherited.
        self.flagHandler.handleFailure(record)
        if error:
            centroid = self.centroidExtractor(record, self.flagHandler)
            self.log.debug(
                "Failed to measure shape for %d at (%f, %f): %s",
                record.getId(),
                centroid.getX(),
                centroid.getY(),
                error,
            )


class HsmSourceMomentsConfig(HsmMomentsConfig):
    """Configuration for HSM adaptive moments measurement for sources."""

    badMaskPlanes = pexConfig.ListField[str](
        doc="Mask planes used to reject bad pixels.", default=["BAD", "SAT"]
    )


@measBase.register("ext_shapeHSM_HsmSourceMoments")
class HsmSourceMomentsPlugin(HsmMomentsPlugin):
    """Plugin for HSM adaptive moments measurement for sources."""

    ConfigClass = HsmSourceMomentsConfig

    def __init__(self, config, name, schema, metadata, logName=None):
        super().__init__(config, name, schema, metadata, logName=logName)
        self.centroidResultKey = afwTable.Point2DKey.addFields(
            schema, name, "Centroid of the source via the HSM shape algorithm", "pixel"
        )
        self.shapeKey = afwTable.QuadrupoleKey.addFields(
            schema,
            name,
            "Adaptive moments of the source via the HSM shape algorithm",
            afwTable.CoordinateType.PIXEL,
        )
        if config.addFlux:
            self.fluxKey = schema.addField(
                schema.join(name, "Flux"), type=float, doc="Flux of the source via the HSM shape algorithm"
            )

    def measure(self, record, exposure):
        """
        Measure adaptive moments of sources given an exposure and set the
        results in the record in place.

        Parameters
        ----------
        record : `~lsst.afw.table.SourceRecord`
            The record where measurement outputs will be stored.
        exposure : `~lsst.afw.image.Exposure`
            The exposure containing the source which needs measurement.

        Raises
        ------
        MeasurementError
            Raised for errors in measurement.
        """
        # Extract the centroid from the record.
        center = self.centroidExtractor(record, self.flagHandler)

        # Get the bounding box of the source's footprint.
        bbox = record.getFootprint().getBBox()

        # Check that the bounding box has non-zero area.
        if bbox.getArea() == 0:
            raise measBase.MeasurementError(self.NO_PIXELS.doc, self.NO_PIXELS.number)

        # Ensure that the centroid is within the bounding box.
        if not bbox.contains(Point2I(center)):
            raise measBase.MeasurementError(self.NOT_CONTAINED.doc, self.NOT_CONTAINED.number)

        # Get the trace radius of the PSF.
        psfSigma = exposure.getPsf().computeShape(center).getTraceRadius()

        # Turn bounding box corners into GalSim bounds.
        xmin, xmax = bbox.getMinX(), bbox.getMaxX()
        ymin, ymax = bbox.getMinY(), bbox.getMaxY()
        bounds = galsim.bounds.BoundsI(xmin, xmax, ymin, ymax)

        # Get the `lsst.meas.base` mask for bad pixels.
        badpix = exposure.mask[bbox].array.copy()
        bitValue = exposure.mask.getPlaneBitMask(self.config.badMaskPlanes)
        badpix &= bitValue

        # Extract the numpy array underlying the image within the bounding box
        # of the source.
        imageArray = exposure.image[bbox].array

        # Create a GalSim image using the extracted array.
        # NOTE: GalSim's HSM uses the FITS convention of 1,1 for the
        # lower-left corner.
        _image = _quickGalsimImageView(imageArray, bounds)

        # Convert the mask of bad pixels to a format suitable for galsim.
        gd = (badpix == 0)
        badpix[gd] = 1
        badpix[~gd] = 0

        _weight_image = _quickGalsimImageView(badpix, bounds)

        # Call the internal method to calculate adaptive moments using GalSim.
        self._calculate(
            record,
            image=_image,
            weight_image=_weight_image,
            sigma=2.5 * psfSigma,
            precision=1.0e-6,
            centroid=center,
        )


class HsmSourceMomentsRoundConfig(HsmSourceMomentsConfig):
    """Configuration for HSM adaptive moments measurement for sources using
    round weight function.
    """

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
    """Plugin for HSM adaptive moments measurement for sources using round
    weight function.
    """

    ConfigClass = HsmSourceMomentsRoundConfig


class HsmPsfMomentsConfig(HsmMomentsConfig):
    """Configuration for HSM adaptive moments measurement for PSFs."""

    useSourceCentroidOffset = pexConfig.Field[bool](
        doc="If True, then draw the PSF to be measured in the coordinate "
        "system of the original image (the PSF model origin - which is "
        "commonly the PSF centroid - may end up near a pixel edge or corner). "
        "If False, then draw the PSF to be measured in a shifted coordinate "
        "system such that the PSF model origin lands precisely in the center "
        "of the central pixel of the PSF image.",
        default=False,
    )

    def setDefaults(self):
        super().setDefaults()
        self.subtractCenter = True


@measBase.register("ext_shapeHSM_HsmPsfMoments")
class HsmPsfMomentsPlugin(HsmMomentsPlugin):
    """Plugin for HSM adaptive moments measurement for PSFs."""

    ConfigClass = HsmPsfMomentsConfig
    _debiased = False

    def __init__(self, config, name, schema, metadata, logName=None):
        super().__init__(config, name, schema, metadata, logName=logName)
        docPrefix = "Debiased centroid" if self._debiased else "Centroid"
        self.centroidResultKey = afwTable.Point2DKey.addFields(
            schema, name, docPrefix + " of the PSF via the HSM shape algorithm", "pixel"
        )
        docPrefix = "Debiased adaptive" if self._debiased else "Adaptive"
        self.shapeKey = afwTable.QuadrupoleKey.addFields(
            schema,
            name,
            docPrefix + " moments of the PSF via the HSM shape algorithm",
            afwTable.CoordinateType.PIXEL,
        )
        if config.addFlux:
            self.fluxKey = schema.addField(
                schema.join(name, "Flux"),
                type=float,
                doc="Flux of the PSF via the HSM shape algorithm",
            )

    def measure(self, record, exposure):
        """
        Measure adaptive moments of the PSF given an exposure and set the
        results in the record in place.

        Parameters
        ----------
        record : `~lsst.afw.table.SourceRecord`
            The record where measurement outputs will be stored.
        exposure : `~lsst.afw.image.Exposure`
            The exposure containing the PSF which needs measurement.

        Raises
        ------
        MeasurementError
            Raised for errors in measurement.
        """
        # Extract the centroid from the record.
        center = self.centroidExtractor(record, self.flagHandler)

        # Retrieve the PSF from the exposure.
        psf = exposure.getPsf()

        # Check that the PSF is not None.
        if not psf:
            raise measBase.MeasurementError(self.NO_PSF.doc, self.NO_PSF.number)

        # Get the bounding box of the PSF.
        psfBBox = psf.computeImageBBox(center)

        # Two methods for getting PSF image evaluated at the source centroid:
        if self.config.useSourceCentroidOffset:
            # 1. Using `computeImage()` returns an image in the same coordinate
            # system as the pixelized image.
            psfImage = psf.computeImage(center)
        else:
            psfImage = psf.computeKernelImage(center)
            # 2. Using `computeKernelImage()` to return an image does not
            # retain any information about the original bounding box of the
            # PSF. We therefore reset the origin to be the same as the
            # pixelized image.
            psfImage.setXY0(psfBBox.getMin())

        # Get the trace radius of the PSF.
        psfSigma = psf.computeShape(center).getTraceRadius()

        # Get the bounding box in the parent coordinate system.
        bbox = psfImage.getBBox(afwImage.PARENT)

        # Turn bounding box corners into GalSim bounds.
        xmin, xmax = bbox.getMinX(), bbox.getMaxX()
        ymin, ymax = bbox.getMinY(), bbox.getMaxY()
        bounds = galsim.bounds.BoundsI(xmin, xmax, ymin, ymax)

        # Adjust the psfImage for noise as needed, and retrieve the mask of bad
        # pixels.
        badpix = self._adjustNoise(psfImage, psfBBox, exposure, record, bounds)

        # Extract the numpy array underlying the PSF image.
        imageArray = psfImage.array

        # Create a GalSim image using the PSF image array.
        _image = _quickGalsimImageView(imageArray, bounds)

        if badpix is not None:
            gd = (badpix == 0)
            badpix[gd] = 1
            badpix[~gd] = 0

            _weight_image = _quickGalsimImageView(badpix, bounds)
        else:
            arr = np.ones(imageArray.shape, dtype=np.int32)
            _weight_image = _quickGalsimImageView(arr, bounds)

        # Decide on the centroid position based on configuration.
        if self.config.useSourceCentroidOffset:
            # If the source centroid offset should be used, use the source
            # centroid.
            centroid = center
        else:
            # Otherwise, use the center of the bounding box of psfImage.
            centroid = Point2D(psfBBox.getMin() + psfBBox.getDimensions() // 2)

        # Call the internal method to calculate adaptive moments using GalSim.
        self._calculate(
            record,
            image=_image,
            weight_image=_weight_image,
            sigma=psfSigma,
            centroid=centroid,
        )

    def _adjustNoise(self, *args) -> None:
        """A noop in the base class, returning None for the bad pixel mask.
        This method is designed to be overridden in subclasses."""
        pass


class HsmPsfMomentsDebiasedConfig(HsmPsfMomentsConfig):
    """Configuration for debiased HSM adaptive moments measurement for PSFs."""

    noiseSource = pexConfig.ChoiceField[str](
        doc="Noise source. How to choose variance of the zero-mean Gaussian noise added to image.",
        allowed={
            "meta": "variance = the 'BGMEAN' metadata entry",
            "variance": "variance = the image's variance plane",
        },
        default="variance",
    )
    seedOffset = pexConfig.Field[int](doc="Seed offset for random number generator.", default=0)
    badMaskPlanes = pexConfig.ListField[str](
        doc="Mask planes used to reject bad pixels.", default=["BAD", "SAT"]
    )

    def setDefaults(self):
        super().setDefaults()
        self.useSourceCentroidOffset = True


@measBase.register("ext_shapeHSM_HsmPsfMomentsDebiased")
class HsmPsfMomentsDebiasedPlugin(HsmPsfMomentsPlugin):
    """Plugin for debiased HSM adaptive moments measurement for PSFs."""

    ConfigClass = HsmPsfMomentsDebiasedConfig
    _debiased = True

    @classmethod
    def getExecutionOrder(cls):
        # Since the standard execution order increases in steps of 1, it's
        # safer to keep the increase by hand to less than 1. The exact value
        # does not matter.
        return cls.FLUX_ORDER + 0.1

    def _adjustNoise(
        self,
        psfImage: afwImage.Image,
        psfBBox: Box2I,
        exposure: afwImage.Exposure,
        record: afwTable.SourceRecord,
        bounds: galsim.bounds.BoundsI,
    ) -> np.ndarray:
        """
        Adjusts noise in the PSF image and updates the bad pixel mask based on
        exposure data. This method modifies `psfImage` in place and returns a
        newly created `badpix` mask.

        Parameters
        ----------
        psfImage : `~lsst.afw.image.Image`
            The PSF image to be adjusted. This image is modified in place.
        psfBBox : `~lsst.geom.Box2I`
            The bounding box of the PSF.
        exposure : `~lsst.afw.image.Exposure`
            The exposure object containing relevant metadata and mask
            information.
        record : `~lsst.afw.table.SourceRecord`
            The source record where measurement outputs will be stored. May be
            modified in place to set flags.
        bounds : `~galsim.bounds.BoundsI`
            The bounding box of the PSF as a GalSim bounds object.

        Returns
        -------
        badpix : `~np.ndarray`
            Numpy image array (np.int32) representing bad pixels, where zero
            indicates good pixels and any nonzero value denotes a bad pixel.

        Raises
        ------
        MeasurementError
            If there's an issue during the noise adjustment process.
        FatalAlgorithmError
            If BGMEAN is not present in the metadata when using the meta noise
            source.
        """
        # Psf image crossing exposure edge is fine if we're getting the
        # variance from metadata, but not okay if we're getting the
        # variance from the variance plane. In both cases, set the EDGE
        # flag, but only fail hard if using variance plane.
        overlap = psfImage.getBBox()
        overlap.clip(exposure.getBBox())
        if overlap != psfImage.getBBox():
            self.flagHandler.setValue(record, self.EDGE.number, True)
            if self.config.noiseSource == "variance":
                self.flagHandler.setValue(record, self.FAILURE.number, True)
                raise measBase.MeasurementError(self.EDGE.doc, self.EDGE.number)

        # Match PSF flux to source.
        psfImage *= record.getPsfInstFlux()

        # Add Gaussian noise to image in 4 steps:
        # 1. Initialize the noise image and random number generator.
        noise = afwImage.Image(psfImage.getBBox(), dtype=psfImage.dtype, initialValue=0.0)
        seed = record.getId() + self.config.seedOffset
        rand = afwMath.Random("MT19937", seed)

        # 2. Generate Gaussian noise image.
        afwMath.randomGaussianImage(noise, rand)

        # 3. Determine the noise scaling based on the noise source.
        if self.config.noiseSource == "meta":
            # Retrieve BGMEAN from the exposure metadata.
            try:
                bgmean = exposure.getMetadata().getAsDouble("BGMEAN")
            except NotFoundError as error:
                raise measBase.FatalAlgorithmError(str(error))
            # Scale the noise by the square root of the background mean.
            noise *= np.sqrt(bgmean)
        elif self.config.noiseSource == "variance":
            # Get the variance image from the exposure and restrict to the
            # PSF bounding box.
            var = afwImage.Image(exposure.variance[psfImage.getBBox()], dtype=psfImage.dtype, deep=True)
            # Scale the noise by the square root of the variance.
            var.sqrt()  # In-place square root.
            noise *= var

        # 4. Add the scaled noise to the PSF image.
        psfImage += noise

        # Masking is needed for debiased PSF moments.
        badpix = afwImage.Mask(psfBBox)
        # NOTE: We repeat the `overlap` calculation in the two lines below to
        # align with the old C++ version and minimize potential discrepancies.
        # There's zero chance this will be a time sink, and using the bbox from
        # the image that's about to be cropped seems safer than using the bbox
        # from a different image, even if they're nominally supposed to have
        # the same bounds.
        overlap = badpix.getBBox()
        overlap.clip(exposure.getBBox())
        badpix[overlap] = exposure.mask[overlap]
        # Pull out the numpy view of the badpix mask image.
        badpix = badpix.array

        bitValue = exposure.mask.getPlaneBitMask(self.config.badMaskPlanes)
        badpix &= bitValue

        return badpix
