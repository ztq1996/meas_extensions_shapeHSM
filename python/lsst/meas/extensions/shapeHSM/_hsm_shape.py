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
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.meas.base as measBase
import lsst.pex.config as pexConfig
from lsst.geom import Point2I

__all__ = [
    "HsmShapeBjConfig",
    "HsmShapeBjPlugin",
    "HsmShapeLinearConfig",
    "HsmShapeLinearPlugin",
    "HsmShapeKsbConfig",
    "HsmShapeKsbPlugin",
    "HsmShapeRegaussConfig",
    "HsmShapeRegaussPlugin",
]


class HsmShapeConfig(measBase.SingleFramePluginConfig):
    """Base configuration for HSM shape measurement."""

    shearType = pexConfig.ChoiceField[str](
        doc="The desired method of PSF correction using GalSim. The first three options return an e-type "
            "distortion, whereas the last option returns a g-type shear.",
        allowed={
            "REGAUSS": "Regaussianization method from Hirata & Seljak (2003)",
            "LINEAR": "A modification by Hirata & Seljak (2003) of methods in Bernstein & Jarvis (2002)",
            "BJ": "From Bernstein & Jarvis (2002)",
            "KSB": "From Kaiser, Squires, & Broadhurst (1995)",
        },
        default="REGAUSS",
    )

    deblendNChild = pexConfig.Field[str](
        doc="Field name for number of deblend children.",
        default="",
    )

    badMaskPlanes = pexConfig.ListField[str](
        doc="Mask planes that indicate pixels that should be excluded from the fit.",
        default=["BAD", "SAT"],
    )


class HsmShapePlugin(measBase.SingleFramePlugin):
    """Base plugin for HSM shape measurement."""

    ConfigClass = HsmShapeConfig
    doc = ""

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

        # Embed the flag definitions in the schema using a flag handler.
        self.flagHandler = measBase.FlagHandler.addFields(schema, name, flagDefs)

        # Utilize a safe centroid extractor that uses the detection footprint
        # as a fallback if necessary.
        self.centroidExtractor = measBase.SafeCentroidExtractor(schema, name)

        self.e1Key = self._addEllipticityField(name, 1, schema, self.doc)
        self.e2Key = self._addEllipticityField(name, 2, schema, self.doc)
        self.sigmaKey = schema.addField(
            schema.join(name, "sigma"),
            type=float,
            doc=f"{self.doc} (shape measurement uncertainty per component)",
        )
        self.resolutionKey = schema.addField(
            schema.join(name, "resolution"), type=float, doc="Resolution factor (0=unresolved, 1=resolved)"
        )
        self.hasDeblendKey = len(config.deblendNChild) > 0

        if self.hasDeblendKey:
            self.deblendKey = schema[config.deblendNChild]

        self.log = logging.getLogger(self.logName)

    @classmethod
    def getExecutionOrder(cls):
        return cls.SHAPE_ORDER

    @staticmethod
    def bboxToGalSimBounds(bbox):
        xmin, xmax = bbox.getMinX(), bbox.getMaxX()
        ymin, ymax = bbox.getMinY(), bbox.getMaxY()
        return galsim._BoundsI(xmin, xmax, ymin, ymax)

    def _addEllipticityField(self, name, n, schema, doc):
        """
        Helper function to add an ellipticity field to a measurement schema.

        Parameters
        ----------
        name : `str`
            Base name of the field.
        n : `int`
            Specifies whether the field is for the first (1) or second (2)
            component.
        schema : `~lsst.afw.table.Schema`
            The schema to which the field is added.
        doc : `str`
            The documentation string that needs to be updated to reflect the
            type and component of the measurement.

        Returns
        -------
        `~lsst.afw.table.KeyD`
            The key associated with the added field in the schema.
        """
        componentLookup = {1: "+ component", 2: "x component"}
        typeLookup = {"e": " of ellipticity", "g": " of estimated shear"}
        name = f"{name}_{self.measTypeSymbol}{n}"
        updatedDoc = f"{doc} ({componentLookup[n]}{typeLookup[self.measTypeSymbol]})"
        return schema.addField(name, type=float, doc=updatedDoc)

    def measure(self, record, exposure):
        """
        Measure the shape of sources given an exposure and set the results in
        the record in place.

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

        if self.hasDeblendKey and record.get(self.deblendKey) > 0:
            raise measBase.MeasurementError(self.PARENT_SOURCE.doc, self.PARENT_SOURCE.number)

        # Get the bounding box of the source's footprint.
        bbox = record.getFootprint().getBBox()

        # Check that the bounding box has non-zero area.
        if bbox.getArea() == 0:
            raise measBase.MeasurementError(self.NO_PIXELS.doc, self.NO_PIXELS.number)

        # Ensure that the centroid is within the bounding box.
        if not bbox.contains(Point2I(center)):
            raise measBase.MeasurementError(self.NOT_CONTAINED.doc, self.NOT_CONTAINED.number)

        # Get the PSF image evaluated at the source centroid.
        psfImage = exposure.getPsf().computeImage(center)
        psfImage.setXY0(0, 0)

        # Get the trace radius of the PSF.
        psfSigma = exposure.getPsf().computeShape(center).getTraceRadius()

        # Turn bounding box corners into GalSim bounds.
        bounds = self.bboxToGalSimBounds(bbox)

        # Get the bounding box of the PSF in the parent coordinate system.
        psfBBox = psfImage.getBBox(afwImage.PARENT)

        # Turn the PSF bounding box corners into GalSim bounds.
        psfBounds = self.bboxToGalSimBounds(psfBBox)

        # Each GalSim image below will match whatever dtype the input array is.
        # NOTE: PSF is already restricted to a small image, so no bbox for the
        # PSF is expected.
        image = galsim._Image(exposure.image[bbox].array, bounds, wcs=None)
        psf = galsim._Image(psfImage.array, psfBounds, wcs=None)

        # Get the `lsst.meas.base` mask for bad pixels.
        subMask = exposure.mask[bbox]
        badpix = subMask.array.copy()  # Copy it since badpix gets modified.
        bitValue = exposure.mask.getPlaneBitMask(self.config.badMaskPlanes)
        badpix &= bitValue

        # Turn badpix to weight where elements set to 1 indicate 'use pixel'
        # and those set to 0 mean 'do not use pixel'. Now, weight will assume
        # the role of badpix, and we will no longer use badpix in our call to
        # EstimateShear().
        gd = badpix == 0
        badpix[gd] = 1
        badpix[~gd] = 0
        weight = galsim._Image(badpix, bounds, wcs=None)

        # Get the statistics control object for sky variance estimation.
        sctrl = afwMath.StatisticsControl()
        sctrl.setAndMask(bitValue)

        # Create a variance image from the exposure.
        # NOTE: Origin defaults to PARENT in all cases accessible from Python.
        variance = afwImage.Image(
            exposure.variance[bbox],
            dtype=exposure.variance.dtype,
            deep=False,
        )

        # Calculate median sky variance for use in shear estimation.
        stat = afwMath.makeStatistics(variance, subMask, afwMath.MEDIAN, sctrl)
        skyvar = stat.getValue(afwMath.MEDIAN)

        # Directly use GalSim's C++/Python interface for shear estimation.
        try:
            # Initialize an instance of ShapeData to store the results.
            shape = galsim.hsm.ShapeData(
                image_bounds=galsim._BoundsI(0, 0, 1, 1),
                observed_shape=galsim._Shear(0j),
                psf_shape=galsim._Shear(0j),
                moments_centroid=galsim._PositionD(0, 0),
            )

            # Prepare various values for GalSim's EstimateShearView.
            recomputeFlux = "FIT"
            precision = 1.0e-6
            guessCentroid = galsim._PositionD(center.getX(), center.getY())
            hsmparams = galsim.hsm.HSMParams.default

            # Estimate shear using GalSim. Arguments are passed positionally
            # to the C++ function. Inline comments specify the Python layer
            # equivalent of each argument for clarity.
            # TODO: [DM-42047] Change to public API when an optimized version
            # is available.
            galsim._galsim.EstimateShearView(
                shape._data,  # shape data buffer (not passed in pure Python)
                image._image,  # gal_image
                psf._image,  # PSF_image
                weight._image,  # weight
                float(skyvar),  # sky_var
                self.config.shearType.upper(),  # shear_est
                recomputeFlux.upper(),  # recompute_flux
                float(2.5 * psfSigma),  # guess_sig_gal
                float(psfSigma),  # guess_sig_PSF
                float(precision),  # precision
                guessCentroid._p,  # guess_centroid
                hsmparams._hsmp,  # hsmparams
            )
        except galsim.hsm.GalSimHSMError as error:
            raise measBase.MeasurementError(str(error), self.GALSIM.number)

        # Set ellipticity and error values based on measurement type.
        if shape.meas_type == "e":
            record.set(self.e1Key, shape.corrected_e1)
            record.set(self.e2Key, shape.corrected_e2)
            record.set(self.sigmaKey, 2.0 * shape.corrected_shape_err)
        else:
            record.set(self.e1Key, shape.corrected_g1)
            record.set(self.e2Key, shape.corrected_g2)
            record.set(self.sigmaKey, shape.corrected_shape_err)

        record.set(self.resolutionKey, shape.resolution_factor)
        self.flagHandler.setValue(record, self.FAILURE.number, shape.correction_status != 0)

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


class HsmShapeBjConfig(HsmShapeConfig):
    """Configuration for HSM shape measurement for the BJ estimator."""

    def setDefaults(self):
        super().setDefaults()
        self.shearType = "BJ"

    def validate(self):
        if self.shearType != "BJ":
            raise pexConfig.FieldValidationError(self.shearType, self, "shearType should be set to 'BJ'.")
        super().validate()


@measBase.register("ext_shapeHSM_HsmShapeBj")
class HsmShapeBjPlugin(HsmShapePlugin):
    """Plugin for HSM shape measurement for the BJ estimator."""

    ConfigClass = HsmShapeBjConfig
    measTypeSymbol = "e"
    doc = "PSF-corrected shear using Bernstein & Jarvis (2002) method"


class HsmShapeLinearConfig(HsmShapeConfig):
    """Configuration for HSM shape measurement for the LINEAR estimator."""

    def setDefaults(self):
        super().setDefaults()
        self.shearType = "LINEAR"

    def validate(self):
        if self.shearType != "LINEAR":
            raise pexConfig.FieldValidationError(self.shearType, self, "shearType should be set to 'LINEAR'.")
        super().validate()


@measBase.register("ext_shapeHSM_HsmShapeLinear")
class HsmShapeLinearPlugin(HsmShapePlugin):
    """Plugin for HSM shape measurement for the LINEAR estimator."""

    ConfigClass = HsmShapeLinearConfig
    measTypeSymbol = "e"
    doc = "PSF-corrected shear using Hirata & Seljak (2003) 'linear' method"


class HsmShapeKsbConfig(HsmShapeConfig):
    """Configuration for HSM shape measurement for the KSB estimator."""

    def setDefaults(self):
        super().setDefaults()
        self.shearType = "KSB"

    def validate(self):
        if self.shearType != "KSB":
            raise pexConfig.FieldValidationError(self.shearType, self, "shearType should be set to 'KSB'.")
        super().validate()


@measBase.register("ext_shapeHSM_HsmShapeKsb")
class HsmShapeKsbPlugin(HsmShapePlugin):
    """Plugin for HSM shape measurement for the KSB estimator."""

    ConfigClass = HsmShapeKsbConfig
    measTypeSymbol = "g"
    doc = "PSF-corrected shear using Kaiser, Squires, & Broadhurst (1995) method"


class HsmShapeRegaussConfig(HsmShapeConfig):
    """Configuration for HSM shape measurement for the REGAUSS estimator."""

    def setDefaults(self):
        super().setDefaults()
        self.shearType = "REGAUSS"

    def validate(self):
        if self.shearType != "REGAUSS":
            raise pexConfig.FieldValidationError(
                self.shearType, self, "shearType should be set to 'REGAUSS'."
            )
        super().validate()


@measBase.register("ext_shapeHSM_HsmShapeRegauss")
class HsmShapeRegaussPlugin(HsmShapePlugin):
    """Plugin for HSM shape measurement for the REGAUSS estimator."""

    ConfigClass = HsmShapeRegaussConfig
    measTypeSymbol = "e"
    doc = "PSF-corrected shear using Hirata & Seljak (2003) 'regaussianization' method"
