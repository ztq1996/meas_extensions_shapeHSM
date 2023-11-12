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
import lsst.afw.image as afwImage
import lsst.pex.config as pexConfig
from numpy import mgrid, sum
import galsim




PLUGIN_NAME = "ext_shapeHSM_HigherOrderMoments"

class HigherOrderMomentsConfig(measBase.SingleFramePluginConfig):
    orderMax = Field[int](
        doc="Maximum order of moments to compute, must be >2",
        default=4,
    )


class HigherOrderMomentsPlugin(measBase.SingleFramePlugin):
    '''Base plugin for higher moments measurement '''

    ConfigClass = HigherOrderMomentsConfig


    def __init__(self, config, name, schema, metadata, logName=None):
        super().__init__(config, name, schema, metadata, logName=logName)

        # Add the output columns to the schema
        # Attn: TQ

        # Define flags for possible issues that might arise during measurement.
        flagDefs = measBase.FlagDefinitionList()
        self.FAILURE = flagDefs.addFailureFlag("General failure flag, set if anything went wrong")

        # Embed the flag definitions in the schema using a flag handler.
        self.flagHandler = measBase.FlagHandler.addFields(schema, name, flagDefs)

        # Utilize a safe centroid extractor that uses the detection footprint
        # as a fallback if necessary.
        self.centroidExtractor = measBase.SafeCentroidExtractor(schema, name)
        self.log = logging.getLogger(self.logName)

    @classmethod
    def getExecutionOrder(cls):
        return cls.SHAPE_ORDER


    def fail(self, record, error=None):
        # Docstring inherited.
        self.flagHandler.handleFailure(record)



    def getAllNames(self, ispsf = False):
        for (p, q) in self._get_pq_full(self.config.orderMax):
            if p + q > 2:
                yield f"{p}{q}"


    @staticmethod
    def _get_pq_full(orderMax):

        pq_list = []

        for n in range(3, orderMax+1):
            p = 0
            q = n

            pq_list.append((p,q))

            while p<n:
                p+=1
                q-=1
                pq_list.append((p,q))
        return pq_list

    @staticmethod
    def _get_suffix(p, q):
        return f"{p}{q}"


    def calc_higher_moments(image,pqlist):
        """
        image : galsim.Image
            Image of the PSF
        pqlist : list of tuples
            Order of moments to measure
        """

        results_list = []
        
        image_array = image.array
        
        y, x = mgrid[:image_array.shape[0],:image_array.shape[1]]+1
        psfresults = galsim.hsm.FindAdaptiveMom(image)
        M = np.zeros((2,2))
        e1 = psfresults.observed_shape.e1
        e2 = psfresults.observed_shape.e2
        sigma4 = psfresults.moments_sigma**4
        c = (1+e1)/(1-e1)
        M[1][1] = np.sqrt(sigma4/(c-0.25*e2**2*(1+c)**2))
        M[0][0] = c*M[1][1]
        M[0][1] = 0.5*e2*(M[1][1]+M[0][0])
        M[1][0] = M[0][1]

        pos = np.array([x-psfresults.moments_centroid.x, y-psfresults.moments_centroid.y])
        inv_M = np.linalg.inv(M)
        sqrt_inv_M = alg.sqrtm(inv_M)
        
        std_pos = np.einsum('ij,jqp->iqp',sqrt_inv_M,pos)
        weight = np.exp(-0.5* np.einsum('ijk,ijk->jk',std_pos,std_pos ))

        std_x, std_y = std_pos[0],std_pos[1]
        
        normalization = sum(image_array*weight)
        image_weight = weight*image_array
        for tup in pqlist:
            p = tup[0]
            q = tup[1]
            
            if p+q<=2:
                raise ValueError('Standardized Moments Work with Order>2')
            
            this_moment = sum(std_x**p*std_y**q*image_weight)/normalization
            results_list.append(this_moment)
            
        return np.array(results_list)



class HigherOrderMomentsSourceConfig(HigherOrderMomentsConfig):

    """
    Configuration for the higher order moments of the source
    """

    badMaskPlanes = pexConfig.ListField[str](
        doc="Mask planes used to reject bad pixels.", default=["BAD", "SAT"]
    )


@measBase.register("ext_shapeHSM_HigherOrderMomentsSource")
class HigherOrderMomentsSourcePlugin(HigherOrderMomentsPlugin)
    """
    Plugin class for Higher Order Moments measurement of the source
    """

    ConfigClass = HigherOrderMomentsSourceConfig

    def __init__(self, config, name, schema, metadata, logName=None):
        super().__init__(config, name, schema, metadata, logName=logName)

        # Add column names for the sources

        for suffix in self.getAllNames():
            schema.addField(schema.join(name, suffix), type=float, doc=f"Higher order moments M_{suffix} for source")


    def measure(self, record, exposure):
        center = self.centroidExtractor(record, self.flagHandler)

        # get the bounding box of the source footprint
        bbox = record.getFootprint().getBBox()

        # get Galsim bounds from bbox

        xmin, xmax = bbox.getMinX(), bbox.getMaxX()
        ymin, ymax = bbox.getMinY(), bbox.getMaxY()
        bounds = galsim.bounds.BoundsI(xmin, xmax, ymin, ymax)

        # get the image array from exposure

        imageArray = exposure[bbox].getImage().array

        # get Galsim image from the image array

        image = galsim.Image(imageArray, bounds=bounds, copy=False)

        # Measure all the moments together to save time
        pqlist = self._get_pq_full(self.config.orderMax)
        try:
            hm_measurement = self.calc_higher_moments(image,pqlist)
        except Exception as e:
            raise measBase.MeasurementError(e)

        # Record the moments 
        for (p,q) in pqlist:

            suffix = self._get_suffix(p,q)
            record.set(self.schema.join(self.name, suffix), M_pq)




class HigherOrderMomentsPSFConfig(HigherOrderMomentsConfig):

    """
    Configuration for the higher order moments of the PSF
    """

    useSourceCentroidOffset = pexConfig.Field[bool](doc="Use source centroid offset?", default=False)


@measBase.register("ext_shapeHSM_HigherOrderMomentsPSF")
class HigherOrderMomentsPSFPlugin(HigherOrderMomentsPlugin)
    """
    Plugin class for Higher Order Moments measurement of the source
    """

    ConfigClass = HigherOrderMomentsPSFConfig

    def __init__(self, config, name, schema, metadata, logName=None):
        super().__init__(config, name, schema, metadata, logName=logName)

        # Add column names for the sources

        for suffix in self.getAllNames():
            schema.addField(schema.join(schema.join(name,'Psf'), suffix), type=float, doc=f"Higher order moments M_{suffix} for PSF")


    def measure(self, record, exposure):
        center = self.centroidExtractor(record, self.flagHandler)


        # get the psf and psf image from the exposure

        psf = exposure.getPsf()

        psfImage = psf.computeImage(center)

        # get the bounding box of the source footprint
        bbox = psfImage.getBBox(afwImage.PARENT)

        # get Galsim bounds from bbox

        xmin, xmax = bbox.getMinX(), bbox.getMaxX()
        ymin, ymax = bbox.getMinY(), bbox.getMaxY()
        bounds = galsim.bounds.BoundsI(xmin, xmax, ymin, ymax)

        # get the image array from exposure

        imageArray = psfImage.array

        # get Galsim image from the image array

        image = galsim.Image(imageArray, bounds=bounds, copy=False)

        # Measure all the moments together to save time
        pqlist = self._get_pq_full(self.config.orderMax)
        try:
            hm_measurement = self.calc_higher_moments(image,pqlist)
        except Exception as e:
            raise measBase.MeasurementError(e)

        # Record the moments 
        for (p,q) in pqlist:

            suffix = self._get_suffix(p,q)
            record.set(schema.join(schema.join(name,'Psf'), suffix), M_pq)



