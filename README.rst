########################
meas_extensions_shapeHSM
########################

``meas_extensions_shapeHSM`` is a package in the `LSST Science Pipelines <https://pipelines.lsst.io>`_.

The ``lsst.meas.extensions.shapeHSM`` Python module provides algorithms for HSM shape measurement.
The algorithm was initially described in https://ui.adsabs.harvard.edu/abs/2003MNRAS.343..459H/abstract,
and was modified later in https://ui.adsabs.harvard.edu/abs/2005MNRAS.361.1287M/abstract.
HSM is named after the primary authors: i) Hirata, Christopher, ii) Seljak, Uros and iii) Mandelbaum, Rachel.
Their implementation of this algorithm lives within `GalSim <https://github.com/GalSim-developers/GalSim>`_,
and this package interacts with the Python layer of GalSim to make the measurements.
We use several GalSim primitives to reduce overhead by bypassing extensive sanity checks.

Documentation is available at https://pipelines.lsst.io/modules/lsst.meas.extensions.shapeHSM/index.html.
