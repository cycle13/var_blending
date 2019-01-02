methods
=======


Cost function
^^^^^^^^^^^^^
At each level, the Fast Fourier Transform (FFT) is carrid out before the minimization of
the following cost function (here we take the blending of GFS and WRF as an example):

:math:`\mathbf{J} = \mathbf{J_{GFS}} + \mathbf{J_{WRF}} + \mathbf{J_{P}}`

Where :math:`\mathbf{J_{GFS}}` and :math:`\mathbf{J_{WRF}}` represent the GFS and WRF term, respectively, and :math:`J_{P}`
is the power constrain term.

GFS and WRF terms
"""""""""""""""""
The GFS and WRF terms determine the ratio of contributions from GFS and WRF, we have:

* :math:`\mathbf{J_{GFS}} = (x-x_{GFS})\mathbf{B_{GFS}}^{-1}(x-x_{GFS})^{T}`
* :math:`\mathbf{J_{WRF}} = (x-x_{WRF})\mathbf{B_{WRF}}^{-1}(x-x_{WRF})^{T}`

Where :math:`x`, :math:`x_{GFS}` and :math:`x_{WRF}` represent the fields of the blended analysis,
GFS and WRF forecasts, respectively, in the spectrum space. :math:`\mathbf{B_{GFS}}` and :math:`\mathbf{B_{WRF}}`
are GFS and WRF errors, which are estimated from:

* :math:`\mathbf{B_{GFS}} = \sum_{i=0}^{N}\frac{\sqrt{[f(x_{GFS})_{i} - f(x_{GA})_{i}]^{2}}}{N}`
* :math:`\mathbf{B_{WRF}} = \sum_{i=0}^{N}\frac{\sqrt{[f(x_{WRF})_{i} - f(x_{GA})_{i}]^{2}}}{N}`

Where :math:`f` represents the inverse FFT (iFFT), :math:`x_{GA}` is the global analysis and :math:`N` represents
the total number of cases used for estimating the model errors. The error matrix is calculated for the range between
the wavenumber 0 and wavenumber :math:`v`, which is determined by the distance between two radiosondes.

Power constrain term
""""""""""""""""""""
Analysis power constrain term is used to make sure that the total analyzed power to be consistent with the power of
WRF forecasts:

* :math:`\mathbf{J_{P}} = [P(x)-P(x_{WRF})]\mathbf{B_{P}}^{-1}[P(x)-P(x_{WRF})]^{T}`

Where :math:`P` represents the Welch's power calculated by the analysis and WRF. :math:`\mathbf{B_{P}}` is the
predefined error for the power, which can be considered as the most tolerable power difference between the analysis
and WRF.
