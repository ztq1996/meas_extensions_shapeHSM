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
import lsst.afw.math as afwMath
import lsst.afw.table as afwTable
import lsst.meas.base as measBase
import lsst.pex.config as pexConfig
import numpy as np
from lsst.geom import Point2D, Point2I
from lsst.pex.exceptions import InvalidParameterError, NotFoundError


class HsmMomentsConfig(measBase.SingleFramePluginConfig):
    """Base configuration for HSM adaptive moments measurement."""

    roundMoments = pexConfig.Field[bool](doc="Use round weight function?", default=False)
    addFlux = pexConfig.Field[bool](doc="Store measured flux?", default=False)
    subtractCenter = pexConfig.Field[bool](doc="Subtract starting center from x/y outputs?", default=False)


class HsmMomentsPlugin(measBase.SingleFramePlugin):
    """Base plugin for HSM adaptive moments measurement."""

    ConfigClass = HsmMomentsConfig

    def __init__(self, config, name, schema, metadata, logName=None):
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

    @classmethod
    def getExecutionOrder(cls):
        return cls.SHAPE_ORDER

    def _calculate(
        self,
        record: afwTable.SourceRecord,
        *,
        image: galsim.Image,
        weight: galsim.Image | None,
        badpix: galsim.Image | None,
        sigma: float,
        precision: float = 1.0e-6,
        centroid: Point2D,
    ):
        """
        Calculate adaptive moments using GalSim's HSM and modify the record in
        place.

        Parameters
        ----------
        record : `~lsst.afw.table.SourceRecord`
            Record to store measurements.
        image : `~galsim.Image`
            Image on which to perform measurements.
        weight : `~galsim.Image`
            The weight image for the galaxy being measured. Can be an int or a
            float array.
        badpix : `~galsim.Image`
            Image representing bad pixels.
        sigma : `float`
            Estimate of object's Gaussian sigma.
        precision : `float`, optional
            Precision for HSM adaptive moments. Default is 1.0e-6.
        centroid : `~lsst.geom.Point2D`
            Centroid guess for HSM adaptive moments.

        Raises
        ------
        MeasurementError
            Raised for errors in measurement.
        """
        # Convert centroid to GalSim's PositionD type.
        guessCentroid = galsim.PositionD(centroid.x, centroid.y)

        try:
            # Attempt to compute HSM moments.
            shape = galsim.hsm.FindAdaptiveMom(
                image,
                weight=weight,
                badpix=badpix,
                guess_sig=sigma,
                precision=precision,
                guess_centroid=guessCentroid,
                strict=True,
                round_moments=self.config.roundMoments,
                hsmparams=None,
            )
        except galsim.hsm.GalSimHSMError as error:
            raise measBase.MeasurementError(error, self.GALSIM.number)

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

        # TODO: calculate errors in shape, centroid?

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
        self.centroidResultKey = afwTable.Point2DKey.addFields(schema, name, "HSM source centroid", "pixel")
        self.shapeKey = afwTable.QuadrupoleKey.addFields(
            schema, name, "HSM source moments", afwTable.CoordinateType.PIXEL
        )
        if config.addFlux:
            self.fluxKey = measBase.FluxResultKey.addFields(schema, name + "_Flux", "HSM source flux")

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
        badpix = exposure[bbox].mask.array.copy()
        bitValue = exposure.mask.getPlaneBitMask(self.config.badMaskPlanes)
        badpix &= bitValue

        # Extract the numpy array underlying the image within the bounding box
        # of the source.
        imageArray = exposure[bbox].getImage().array

        # Create a GalSim image using the extracted array.
        # NOTE: GalSim's HSM uses the FITS convention of 1,1 for the
        # lower-left corner.
        image = galsim.Image(imageArray, bounds=bounds, copy=False)

        # Convert the mask of bad pixels to a format suitable for galsim.
        # NOTE: galsim.Image will match whatever dtype the input array is
        # (here int32).
        badpix = galsim.Image(badpix, bounds=bounds)

        # Call the internal method to calculate adaptive moments using GalSim.
        self._calculate(
            record,
            image=image,
            weight=None,
            badpix=badpix,
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

    useSourceCentroidOffset = pexConfig.Field[bool](doc="Use source centroid offset?", default=False)
    debiasedPsfMoments = pexConfig.Field[bool](doc="Debias PSF moments?", default=False)
    # The rest of the config is only relevant if debias is True.
    noiseSource = pexConfig.ChoiceField[str](
        doc="Noise source. How to choose variance of the zero-mean Gaussian noise added to image.",
        allowed={
            "meta": "variance = the 'BGMEAN' metadata entry",
            "variance": "variance = the image's variance plane",
        },
        default="variance",
    )
    seedOffset = pexConfig.Field[int](doc="Seed offset for random number generator.", default=0)


@measBase.register("ext_shapeHSM_HsmPsfMoments")
class HsmPsfMomentsPlugin(HsmMomentsPlugin):
    """Plugin for HSM adaptive moments measurement for PSFs."""

    ConfigClass = HsmPsfMomentsConfig

    def __init__(self, config, name, schema, metadata, logName=None):
        super().__init__(config, name, schema, metadata, logName=logName)
        self.centroidResultKey = afwTable.Point2DKey.addFields(
            schema, name, "HSM PSF centroid", "pixel"
        )  # not sure if we report this but let's keep it here for now!!!
        self.shapeKey = afwTable.QuadrupoleKey.addFields(
            schema, name, "HSM PSF moments", afwTable.CoordinateType.PIXEL
        )
        if config.addFlux:
            self.fluxKey = measBase.FluxResultKey.addFields(schema, name + "_Flux", "HSM PSF flux")

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

        # Two methods for getting PSF image evaluated at the source centroid:
        if self.config.useSourceCentroidOffset:
            # 1. Using `computeImage()` returns an image in the same coordinate
            # system as the pixelized image.
            psfImage = psf.computeImage(center)
        else:
            psfImage = psf.computeKernelImage(center)
            # 2. Using `computeKernelImage()` to return an image does not
            # retain any information about the origial bounding box of the
            # PSF. We therefore reset the origin to be the same as the
            # pixelized image.
            psfImage.setXY0(psf.computeImageBBox(center).getMin())

        # Get the trace radius of the PSF.
        psfSigma = psf.computeShape(center).getTraceRadius()

        # Get the bounding box in the parent coordinate system.
        bbox = psfImage.getBBox(afwImage.PARENT)

        # Turn bounding box corners into GalSim bounds.
        xmin, xmax = bbox.getMinX(), bbox.getMaxX()
        ymin, ymax = bbox.getMinY(), bbox.getMaxY()
        bounds = galsim.bounds.BoundsI(xmin, xmax, ymin, ymax)

        if self.config.debiasedPsfMoments:
            # Match PSF flux to source.
            psfImage *= record.getPsfInstFlux()

            # Add Gaussian noise to image in 4 steps:
            # 1. Initialize the noise image and random number generator.
            # if isinstance(psfImage, afwImage.ImageF):
            #     noise = afwImage.ImageF(psfImage.getBBox())
            # elif isinstance(psfImage, afwImage.ImageD):
            #     noise = afwImage.ImageD(psfImage.getBBox())
            noise = afwImage.Image(psfImage.getBBox(), dtype=psfImage.dtype)
            seed = record.getId() + self.config.seedOffset
            rand = afwMath.Random("MT19937", seed)

            # 2. Generate Gaussian noise image.
            afwMath.randomGaussianImage(noise, rand)

            # 3. Determine the noise scaling based on the noise source.
            if self.config.noiseSource == "meta":
                # Retrieve the BGMEAN from the exposure metadata.
                try:
                    bgmean = exposure.getMetadata().getAsDouble("BGMEAN")
                except NotFoundError as error:
                    raise measBase.FatalAlgorithmError(error)
                # Scale the noise by the square root of the background mean.
                noise *= np.sqrt(bgmean)
            elif self.config.noiseSource == "variance":
                # Get the variance image from the exposure and restrict to the
                # PSF bounding box.
                var = afwImage.ImageF(
                    exposure.getMaskedImage().getVariance(),
                    psfImage.getBBox(),
                    afwImage.PARENT,
                    True,
                )
                # Scale the noise by the square root of the variance.
                noise *= np.sqrt(var)

            # 4. Add the scaled noise to the PSF image
            psfImage += noise

            # Masking is needed for debiased PSF moments.
            # FIXME: setting badpix leads to nans!!!
            badpix = exposure[bbox].mask.array.copy()
            bitValue = exposure.mask.getPlaneBitMask(self.config.badMaskPlanes)
            badpix &= bitValue
        else:
            badpix = None

        # Extract the numpy array underlying the PSF image.
        imageArray = psfImage.array

        # Create a GalSim image using the PSF image array.
        image = galsim.Image(imageArray, bounds=bounds, copy=False)

        # Decide on the centroid position based on configuration.
        if self.config.useSourceCentroidOffset:
            # If the source centroid offset should be used, use the source
            # centroid.
            centroid = center
        else:
            # Otherwise, use the center of the bounding box.
            centroid = Point2D(bbox.getMin() + bbox.getDimensions() / 2)

        # Call the internal method to calculate adaptive moments using GalSim.
        self._calculate(
            record,
            image=image,
            weight=None,
            badpix=badpix,
            sigma=psfSigma,
            centroid=centroid,
        )


class HsmPsfMomentsDebiasedConfig(HsmPsfMomentsConfig):
    """Configuration for debiased HSM adaptive moments measurement for PSFs."""

    def setDefaults(self):
        super().setDefaults()
        self.useSourceCentroidOffset = True
        self.debiasedPsfMoments = True


@measBase.register("ext_shapeHSM_HsmPsfMomentsDebiased")
class HsmPsfMomentsDebiasedPlugin(HsmPsfMomentsPlugin):
    """Plugin for debiased HSM adaptive moments measurement for PSFs."""

    ConfigClass = HsmPsfMomentsDebiasedConfig

    @classmethod
    def getExecutionOrder(cls):
        return cls.FLUX_ORDER + 1
