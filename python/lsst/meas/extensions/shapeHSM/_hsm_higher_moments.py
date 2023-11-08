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
from numpy import mgrid, sum


PLUGIN_NAME = "ext_shapeHSM_HigherOrderMoments"

class HigherOrderMomentsConfig(measBase.SingleFramePluginConfig):
    order = Field[int](
        doc="Order of moments to compute",
        default=4,
    )


@measBase.register(PLUGIN_NAME)
class HigherOrderMomentsPlugin(measBase.SingleFramePlugin):
    ConfigClass = HigherOrderMomentsConfig

    @classmethod
    def getExecutionOrder(cls):
        return cls.SHAPE_ORDER

    def __init__(self, config, name, schema, metadata, logName=None):
        super().__init__(config, name, schema, metadata, logName=logName)

        # Add the output columns to the schema
        # Attn: TQ
        for suffix in self.getAllNames():
            schema.addField(schema.join(name, suffix), type=float, doc="Higher order moments")

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

        ## Attn: TQ Do or call the main computation here

        ## Attn: Arun - store results below this.
        for (p,q) in self._get_pq_full(self.config.order):
            try:
                # Compute the pq moment as M_pq
                M_pq = -99  # dummy value
            except Exception as e:
                raise measBase.MeasurementError(e)

            suffix = self._get_suffix(p,q)
            record.set(self.schema.join(self.name, suffix), M_pq)

    def getAllNames(self):
        for (p, q) in self._get_pq_full(self.config.order):
            if p + q > 2:
                yield f"{p}{q}"


    @staticmethod
    def _get_pq_full(nmax):

        pq_list = []

        for n in range(2, nmax+1):
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

def get_all_moments_fast(image,pqlist):
    """
    image : galsim.Image
        Image of the PSF
    pqlist : list of tuples
        Order of moments to measure
    """

    results_list = []

    image_results = galsim.hsm.FindAdaptiveMom(image)

    image_array = image.array
​
    y, x = mgrid[:image_array.shape[0],:image_array.shape[1]]+1
​
    image.scale = 1.0

    ### Build the second moment metrix of the image (maybe use the standard way to convert e1, e2, sigma -> Mxx, Mxy, Myy)
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
​
    pos = np.array([x-psfresults.moments_centroid.x, y-psfresults.moments_centroid.y])
    pos = np.swapaxes(pos,0,1)
    pos = np.swapaxes(pos,1,2)
​
    inv_M = np.linalg.inv(M)
    sqrt_inv_M = alg.sqrtm(inv_M)
    std_pos = np.zeros(pos.shape)
    weight = np.zeros(pos.shape[0:2])
    for i in range(pos.shape[0]):
        for j in range(pos.shape[1]):
            this_pos = pos[i][j]
            this_standard_pos = np.matmul(sqrt_inv_M, this_pos)
            std_pos[i][j] = this_standard_pos
            weight[i][j] = np.exp(-0.5* this_standard_pos.dot(this_standard_pos))
​
    std_x, std_y = std_pos[:,:,0],std_pos[:,:,1]

    for tup in pqlist:
        p = tup[0]
        q = tup[1]

        if q+p==2:
            if p==2:
                this_moment = image_results.observed_shape.e1
            elif q==2:
                this_moment = image_results.observed_shape.e2
            else:
                this_moment = image_results.moments_sigma
        else:
            this_moment = sum(std_x**p*std_y**q*weight*image)/sum(image*weight)

        results_list.append(this_moment)

    return np.array(results_list)
