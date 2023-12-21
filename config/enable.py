# Enable HSM measurements
#
# The 'config' should be a SourceMeasurementConfig.
#
# We activate the REGAUSS PSF-corrected shape measurement, the adaptive moments of the source and PSF and the
# higher-order moments measurements of the same.

import lsst.meas.extensions.shapeHSM

config.plugins.names |= [
    "ext_shapeHSM_HsmShapeRegauss",
    "ext_shapeHSM_HsmSourceMoments",
    "ext_shapeHSM_HsmPsfMoments",
    "ext_shapeHSM_HsmSourceMomentsRound",
    "ext_shapeHSM_HigherOrderMomentsSource",
    "ext_shapeHSM_HigherOrderMomentsPSF",
]
config.slots.shape = "ext_shapeHSM_HsmSourceMoments"
config.slots.psfShape = "ext_shapeHSM_HsmPsfMoments"
