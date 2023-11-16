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

__all__ = ["HigherOrderMomentsPlugin","HigherOrderMomentsSourcePlugin", "HigherOrderMomentsPSFPlugin" ]

import logging
import lsst.meas.base as measBase
from lsst.pex.config import Field
import lsst.afw.image as afwImage
import lsst.pex.config as pexConfig
import galsim
import numpy as np

import logging
import sys
import scipy




PLUGIN_NAME = "ext_shapeHSM_HigherOrderMoments"

class HigherOrderMomentsConfig(measBase.SingleFramePluginConfig):
    orderMax = Field(
        doc="Maximum order of moments to compute, must be >2",
        default=4, dtype = int
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

        self.NO_PIXELS = flagDefs.add("flag_no_pixels", "No pixels to measure")
        self.NOT_CONTAINED = flagDefs.add(
            "flag_not_contained", "Center not contained in footprint bounding box"
        )

        self.GALSIM = flagDefs.add("flag_galsim", "GalSim failure")

        self.NO_PSF = flagDefs.add("flag_no_psf", "Exposure lacks PSF")

        # Embed the flag definitions in the schema using a flag handler.
        self.flagHandler = measBase.FlagHandler.addFields(schema, name, flagDefs)

        # Utilize a safe centroid extractor that uses the detection footprint
        # as a fallback if necessary.
        self.centroidExtractor = measBase.SafeCentroidExtractor(schema, name)
        self.log = logging.getLogger(self.logName)

        # initialize a logger
        
        logger = logging.getLogger()
        logger.setLevel(logging.INFO)
        formatter = logging.Formatter('%(asctime)s | %(levelname)s | %(message)s')

        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.DEBUG)
        stdout_handler.setFormatter(formatter)

        file_handler = logging.FileHandler('logs_source.log')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)


        logger.addHandler(file_handler)
        logger.addHandler(stdout_handler)
        self.logger = logger

    @classmethod
    def getExecutionOrder(cls):
        return cls.SHAPE_ORDER


    def fail(self, record, error=None):
        # Docstring inherited.
        self.flagHandler.handleFailure(record)



    def getAllNames(self):
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


    def calc_higher_moments(self,image,pqlist):
        """
        image : galsim.Image
            Image of the PSF
        pqlist : list of tuples
            Order of moments to measure
        """

        results_list = []
        
        image_array = image.array
        
        # self.logger.info(image.bounds)
        
        y, x = np.mgrid[:image_array.shape[0],:image_array.shape[1]]
        
        # self.logger.info('try hsm findadaptivemom')
        try:
            
            psfresults = galsim.hsm.FindAdaptiveMom(image)
        except galsim.hsm.GalSimHSMError as error:
            
            raise measBase.MeasurementError(str(error), self.GALSIM.number)
        # self.logger.info(psfresults)
        
        
        M = np.zeros((2,2))
        e1 = psfresults.observed_shape.e1
        e2 = psfresults.observed_shape.e2
        sigma4 = psfresults.moments_sigma**4
        c = (1+e1)/(1-e1)
        M[1][1] = np.sqrt(sigma4/(c-0.25*e2**2*(1+c)**2))
        M[0][0] = c*M[1][1]
        M[0][1] = 0.5*e2*(M[1][1]+M[0][0])
        M[1][0] = M[0][1]
        
        inv_M = np.linalg.inv(M)
                
        #find the sqrt of inv_M
        evalues, evectors = np.linalg.eig(inv_M)
        # Ensuring square root matrix exists
        
        sqrt_inv_M = evectors * np.sqrt(evalues) @ np.linalg.inv(evectors)

#         self.logger.info("Finish sqrt_inv_M")
#         self.logger.info(image.bounds.getXMin())
        
        
        pos = np.array([x-(psfresults.moments_centroid.x - image.bounds.getXMin()), y-(psfresults.moments_centroid.y - image.bounds.getYMin())])

        
        std_pos = np.einsum('ij,jqp->iqp',sqrt_inv_M,pos)
        weight = np.exp(-0.5* np.einsum('ijk,ijk->jk',std_pos,std_pos ))

        std_x, std_y = std_pos[0],std_pos[1]
        
        normalization = np.sum(image_array*weight)
        image_weight = weight*image_array
        
        
        for tup in pqlist:
            p = tup[0]
            q = tup[1]
            
            if p+q<=2:
                raise ValueError('Standardized Moments Work with Order>2')
            
            this_moment = np.sum(std_x**p*std_y**q*image_weight)/normalization
            results_list.append(this_moment)
            
        return np.array(results_list)


class HigherOrderMomentsSourceConfig(HigherOrderMomentsConfig):

    """
    Configuration for the higher order moments of the source
    """

    badMaskPlanes = pexConfig.ListField(
        doc="Mask planes used to reject bad pixels.", default=["BAD", "SAT"], dtype = str
    )


@measBase.register("ext_shapeHSM_HigherOrderMomentsSource")
class HigherOrderMomentsSourcePlugin(HigherOrderMomentsPlugin):
    """
    Plugin class for Higher Order Moments measurement of the source
    """

    ConfigClass = HigherOrderMomentsSourceConfig

    def __init__(self, config, name, schema, metadata, logName=None):
        super().__init__(config, name, schema, metadata, logName=logName)

        # Add column names for the sources

        for suffix in self.getAllNames():
            schema.addField(schema.join(name, suffix), type=float, doc=f"Higher order moments M_{suffix} for source")
            
        # self.logger.info('finish initialization')


    def measure(self, record, exposure):
        
        center = self.centroidExtractor(record, self.flagHandler)
        
        # get the bounding box of the source footprint
        bbox = record.getFootprint().getBBox()
        
        # self.logger.info(bbox)

        # Check that the bounding box has non-zero area.
        if bbox.getArea() == 0:
            raise measBase.MeasurementError(self.NO_PIXELS.doc, self.NO_PIXELS.number)

        # # Ensure that the centroid is within the bounding box.
        # if not bbox.contains(Point2I(center)):
        #     raise measBase.MeasurementError(self.NOT_CONTAINED.doc, self.NOT_CONTAINED.number)

        # get Galsim bounds from bbox

        xmin, xmax = bbox.getMinX(), bbox.getMaxX()
        ymin, ymax = bbox.getMinY(), bbox.getMaxY()
        bounds = galsim.bounds.BoundsI(xmin, xmax, ymin, ymax)

        # get the image array from exposure

        imageArray = exposure[bbox].getImage().array
        # np.savetxt('savefig.txt', imageArray)
        # self.logger.info('array' + str(imageArray))

        # get Galsim image from the image array

        image = galsim.Image(imageArray, bounds=bounds, copy=False)
        

        # Measure all the moments together to save time
        pqlist = self._get_pq_full(self.config.orderMax)
        try:
            hm_measurement = self.calc_higher_moments(image,pqlist)
        except Exception as e:
            raise measBase.MeasurementError(e)
        
        
        # Record the moments 
        for i in range(len(hm_measurement)):

            (p,q) = pqlist[i]
            M_pq = hm_measurement[i]
            suffix = self._get_suffix(p,q)
            this_column_name = self.name + f"_{p}{q}"
            
            record.set(this_column_name, M_pq)





class HigherOrderMomentsPSFConfig(HigherOrderMomentsConfig):

    """
    Configuration for the higher order moments of the PSF
    """

    useSourceCentroidOffset = pexConfig.Field(doc="Use source centroid offset?", default=False, dtype = bool)


@measBase.register("ext_shapeHSM_HigherOrderMomentsPSF")
class HigherOrderMomentsPSFPlugin(HigherOrderMomentsPlugin):
    """
    Plugin class for Higher Order Moments measurement of the source
    """

    ConfigClass = HigherOrderMomentsPSFConfig

    def __init__(self, config, name, schema, metadata, logName=None):
        super().__init__(config, name, schema, metadata, logName=logName)

        # Add column names for the PSFs

        for suffix in self.getAllNames():
            schema.addField(schema.join(name, suffix), type=float, doc=f"Higher order moments M_{suffix} for PSF")


    def measure(self, record, exposure):
        center = self.centroidExtractor(record, self.flagHandler)


        # get the psf and psf image from the exposure

        psf = exposure.getPsf()

        # check if the psf is none:

        if not psf:
            raise measBase.MeasurementError(self.NO_PSF.doc, self.NO_PSF.number)

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
            
        # self.logger.info(hm_measurement)

        # Record the moments 
        for i in range(len(hm_measurement)):

            (p,q) = pqlist[i]
            M_pq = hm_measurement[i]

            suffix = self._get_suffix(p,q)
            this_column_name = self.name + f"_{p}{q}"
            record.set(this_column_name, M_pq)



