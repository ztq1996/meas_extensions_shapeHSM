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

__all__ = (
    "HigherOrderMomentsConfig",
    "HigherOrderMomentsPlugin",
    "HigherOrderMomentsPSFConfig",
    "HigherOrderMomentsPSFPlugin",
    "HigherOrderMomentsSourceConfig",
    "HigherOrderMomentsSourcePlugin",
)

import lsst.geom as geom
import lsst.meas.base as measBase
import numpy as np
from lsst.pex.config import Field, FieldValidationError, ListField


class HigherOrderMomentsConfig(measBase.SingleFramePluginConfig):
    min_order = Field[int](
        doc="Minimum order of the higher order moments to compute",
        default=3,
    )

    max_order = Field[int](
        doc="Maximum order of the higher order moments to compute",
        default=4,
    )

    def validate(self):
        if self.min_order > self.max_order:
            raise FieldValidationError(
                self.__class__.min_order, self, "min_order must be less than or equal to max_order"
            )
        super().validate()


class HigherOrderMomentsPlugin(measBase.SingleFramePlugin):
    """Base plugin for higher moments measurement"""

    ConfigClass = HigherOrderMomentsConfig

    def __init__(self, config, name, schema, metadata, logName=None):
        super().__init__(config, name, schema, metadata, logName=logName)

        # Define flags for possible issues that might arise during measurement.
        flagDefs = measBase.FlagDefinitionList()
        self.FAILURE = flagDefs.addFailureFlag("General failure flag, set if anything went wrong")

        # Embed the flag definitions in the schema using a flag handler.
        self.flagHandler = measBase.FlagHandler.addFields(schema, name, flagDefs)

        self.pqlist = self._get_pq_full()

    @classmethod
    def getExecutionOrder(cls):
        return cls.FLUX_ORDER

    def fail(self, record, error=None):
        # Docstring inherited.
        self.flagHandler.handleFailure(record)

    def _get_pq_full(self):
        """Get a list of the orders to measure as a tuple.

        Returns
        -------
        pqlist: `list` [`tuples`]
            A list of tuples of the form (p, q) where p and q denote the order
            in x and y direction.
        """
        pq_list = []

        for n in range(self.config.min_order, self.config.max_order + 1):
            p = 0
            q = n

            pq_list.append((p, q))

            while p < n:
                p += 1
                q -= 1
                pq_list.append((p, q))

        return pq_list

    def _generate_suffixes(self):
        """Generator of suffixes 'pq'."""
        for p, q in self.pqlist:
            yield f"{p}{q}"

    def _generate_powers_of_standard_positions(self, std_x, std_y):
        std_x_powers, std_y_powers = {0: 1.0, 1: std_x}, {0: 1.0, 1: std_y}

        for p in range(2, self.config.max_order + 1):
            std_x_powers[p] = std_x_powers[p - 1] * std_x

        for q in range(2, self.config.max_order + 1):
            std_y_powers[q] = std_y_powers[q - 1] * std_y

        return std_x_powers, std_y_powers

    def _calculate_higher_order_moments(
        self,
        image,
        center,
        M,
        badpix=None,
        set_masked_pixels_to_zero=False,
        use_linear_algebra=False,
    ):
        """
        Calculate the higher order moments of an image.

        Parameters
        ----------
        image : `~lsst.afw.image.Image`
            Image from which the moments need to be measured (source or PSF).
        center: `~lsst.geom.Point2D`
            First order moments of ``image``. This is used as the peak of the
            Gaussian weight image.
        M : `~numpy.ndarray`
            A 2x2 numpy array representing the second order moments of
            ``image``. This is used to generate the Gaussian weight image.
        badpix : `~numpy.ndarray` or None
            A 2D array having the same shape and orientation as ``image.array``
            that denotes which pixels are bad and should not be accounted for
            when computing the moments.
        set_masked_pixels_to_zero: `bool`
            Whether to treat pixels corresponding to ``badpix`` should be set
            to zero, or replaced by a scaled version of the weight image.
            This is ignored if ``badpix`` is None.
        use_linear_algebra: `bool`
            Use linear algebra operations (eigen decomposition and inverse) to
            calculate the moments? If False, use the specialized formulae for
            2x2 matrix.

        Returns
        -------
        results : `dict`
            A dictionary mapping the order of the moments expressed as tuples
            to the corresponding higher order moments.
        """

        bbox = image.getBBox()
        image_array = image.array

        y, x = np.mgrid[: image_array.shape[0], : image_array.shape[1]]

        if use_linear_algebra:
            inv_M = np.linalg.inv(M)

            evalues, evectors = np.linalg.eig(inv_M)

            sqrt_inv_M = evectors * np.sqrt(evalues) @ np.linalg.inv(evectors)
        else:
            # This is the implementation of Eq. 6 in Hirata & Seljak (2003):
            # https://arxiv.org/pdf/astro-ph/0301054.pdf
            D = M[0, 0] * M[1, 1] - M[0, 1] * M[1, 0]
            sqrt_D = D**0.5
            sqrt_eta = (D * (M[0, 0] + M[1, 1] + 2 * sqrt_D)) ** 0.5
            sqrt_inv_M = (1 / sqrt_eta) * np.array(
                [[M[1, 1] + sqrt_D, -M[0, 1]], [-M[1, 0], M[0, 0] + sqrt_D]]
            )

        pos = np.array([x - (center.getX() - bbox.getMinX()), y - (center.getY() - bbox.getMinY())])

        std_pos = np.einsum("ij,jqp->iqp", sqrt_inv_M, pos)
        weight = np.exp(-0.5 * np.einsum("ijk,ijk->jk", std_pos, std_pos))

        image_weight = weight * image_array

        # Modify only the weight, not the image_array, since it will change the
        # pixel values forever!!!
        if badpix is not None and badpix.any():
            if set_masked_pixels_to_zero:
                # This is how HSM treats bad pixels to compute the quadrupole
                # moments.
                image_weight[badpix] = 0.0
            else:
                # This is how Piff treats bad pixels to compute the
                # higher-order moments.
                scale = image_array[~badpix].sum() / weight[~badpix].sum()
                image_weight[badpix] = (weight[badpix] ** 2) * scale

        normalization = np.sum(image_weight)

        std_x, std_y = std_pos
        std_x_powers, std_y_powers = self._generate_powers_of_standard_positions(std_x, std_y)

        results = {}
        for p, q in self.pqlist:
            results[(p, q)] = np.sum(std_x_powers[p] * std_y_powers[q] * image_weight) / normalization

        return results


class HigherOrderMomentsSourceConfig(HigherOrderMomentsConfig):
    """Configuration for the measurement of higher order moments of objects."""

    badMaskPlanes = ListField[str](
        doc="Mask planes used to reject bad pixels.",
        default=["BAD", "SAT"],
    )

    setMaskedPixelsToZero = Field[bool](
        doc="Set masked pixels to zero? If False, they are replaced by the "
        "scaled version of the adaptive weights.",
        default=False,
    )


@measBase.register("ext_shapeHSM_HigherOrderMomentsSource")
class HigherOrderMomentsSourcePlugin(HigherOrderMomentsPlugin):
    """Plugin for Higher Order Moments measurement of objects.

    The moments are measured in normalized coordinates, where the normalized x
    axis is along the major axis and the normalized y axis along the minor.
    The moments are dependent only on the light profile, and does not scale
    with the size or orientation of the object.

    For any well-sampled image, the zeroth order moment is 1,
    the first order moments are 0, and the second order moments are 0.5 for xx
    and yy and 0 for xy. For a symmetric profile, the moments are zeros if
    either of the indices is odd.

    Notes
    -----
    This plugin requires the `ext_shapeHSM_HsmSourceMoments` plugin to be
    enabled in order to measure the higher order moments, and raises a
    FatalAlgorithmError otherwise. For accurate results, the weight function
    used must match those used for first and second order moments. Hence, this
    plugin does not use slots for centroids and shapes, but instead uses those
    measured by the `ext_shapeHSM_HsmSourceMoments` explicitly.

    The only known failure mode of this plugin is if
    `ext_shapeHSM_HsmSourceMoments` measurement failed. The flags of that
    plugin are informative here as well and should be used to filter out
    unreliable measurements.
    """

    ConfigClass = HigherOrderMomentsSourceConfig

    def __init__(self, config, name, schema, metadata, logName=None):
        super().__init__(config, name, schema, metadata, logName=logName)

        for suffix in self._generate_suffixes():
            schema.addField(
                schema.join(name, suffix),
                type=float,
                doc=f"Higher order moments M_{suffix} for source",
            )

    def measure(self, record, exposure):
        # Docstring inherited.
        M = np.zeros((2, 2))
        try:
            center = geom.Point2D(
                record["ext_shapeHSM_HsmSourceMoments_x"],
                record["ext_shapeHSM_HsmSourceMoments_y"],
            )
            M[0, 0] = record["ext_shapeHSM_HsmSourceMoments_xx"]
            M[1, 1] = record["ext_shapeHSM_HsmSourceMoments_yy"]
            M[0, 1] = M[1, 0] = record["ext_shapeHSM_HsmSourceMoments_xy"]
        except KeyError:
            raise measBase.FatalAlgorithmError("'ext_shapeHSM_HsmSourceMoments' plugin must be enabled.")

        # Obtain the bounding box of the source footprint
        bbox = record.getFootprint().getBBox()

        bitValue = exposure.mask.getPlaneBitMask(self.config.badMaskPlanes)
        badpix = (exposure.mask[bbox].array & bitValue) != 0

        # Measure all the moments together to save time
        try:
            hm_measurement = self._calculate_higher_order_moments(
                exposure.image[bbox],
                center,
                M,
                badpix,
                set_masked_pixels_to_zero=self.config.setMaskedPixelsToZero,
            )
        except Exception as e:
            raise measBase.MeasurementError(e)

        # Record the moments
        for (p, q), M_pq in hm_measurement.items():
            column_key = self.name + f"_{p}{q}"
            record.set(column_key, M_pq)


class HigherOrderMomentsPSFConfig(HigherOrderMomentsConfig):
    """Configuration for the higher order moments of the PSF."""

    useSourceCentroidOffset = Field[bool](
        doc="Use source centroid offset?",
        default=False,
    )


@measBase.register("ext_shapeHSM_HigherOrderMomentsPSF")
class HigherOrderMomentsPSFPlugin(HigherOrderMomentsPlugin):
    """Plugin for Higher Order Moments measurement of PSF models.

    The moments are measured in normalized coordinates, where the normalized x
    axis is along the major axis and the normalized y axis along the minor.
    The moments are dependent only on the light profile, and does not scale
    with the size or orientation of the object.

    For any well-sampled image, the zeroth order moment is 1,
    the first order moments are 0, and the second order moments are 0.5 for xx
    and yy and 0 for xy. For a symmetric profile, the moments are zeros if
    either of the indices is odd.

    Notes
    -----
    This plugin requires the `ext_shapeHSM_HsmPsfMoments` plugin to be
    enabled in order to measure the higher order moments, and raises a
    FatalAlgorithmError otherwise. The weight function is parametrized by the
    shape measured from `ext_shapeHSM_HsmPsfMoments` but for efficiency
    reasons, uses the slot centroid to evaluate the PSF model.

    The only known failure mode of this plugin is if
    `ext_shapeHSM_HsmPsfMoments` measurement failed. The flags of that
    plugin are informative here as well and should be used to filter out
    unreliable measurements.
    """

    ConfigClass = HigherOrderMomentsPSFConfig

    def __init__(self, config, name, schema, metadata, logName=None):
        super().__init__(config, name, schema, metadata, logName=logName)
        # Use the standard slot centroid to use the shared PSF model with
        # other plugins.
        self.centroidExtractor = measBase.SafeCentroidExtractor(schema, name)

        for suffix in self._generate_suffixes():
            schema.addField(
                schema.join(name, suffix),
                type=float,
                doc=f"Higher order moments M_{suffix} for PSF",
            )

    def measure(self, record, exposure):
        # Docstring inherited.
        M = np.zeros((2, 2))
        try:
            M[0, 0] = record["ext_shapeHSM_HsmPsfMoments_xx"]
            M[1, 1] = record["ext_shapeHSM_HsmPsfMoments_yy"]
            M[0, 1] = M[1, 0] = record["ext_shapeHSM_HsmPsfMoments_xy"]
        except KeyError:
            raise measBase.FatalAlgorithmError("'ext_shapeHSM_HsmPsfMoments' plugin must be enabled.")

        psf = exposure.getPsf()

        centroid = self.centroidExtractor(record, self.flagHandler)
        if self.config.useSourceCentroidOffset:
            psfImage = psf.computeImage(centroid)
            psfCenter = centroid
            # Undo what subtractCenter config did.
            # This operation assumes subtractCenter was set to True (default)
            # in the ext_shapeHSM_HsmPsfMomentsConfig and does not have
            # access to it.
            psfCenter.x += record["ext_shapeHSM_HsmPsfMoments_x"]
            psfCenter.y += record["ext_shapeHSM_HsmPsfMoments_y"]
        else:
            psfImage = psf.computeKernelImage(centroid)
            center0 = geom.Point2I(centroid)
            xy0 = geom.Point2I(center0.x + psfImage.getX0(), center0.y + psfImage.getY0())
            psfImage.setXY0(xy0)
            psfBBox = psfImage.getBBox()
            psfCenter = geom.Point2D(psfBBox.getMin() + psfBBox.getDimensions() // 2)

        # Measure all the moments together to save time
        try:
            hm_measurement = self._calculate_higher_order_moments(psfImage, psfCenter, M)
        except Exception as e:
            raise measBase.MeasurementError(e)

        # Record the moments
        for (p, q), M_pq in hm_measurement.items():
            column_key = self.name + f"_{p}{q}"
            record.set(column_key, M_pq)
