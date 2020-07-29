.. highlight:: matlab
.. _example_extracting_gauss_constraints:

***************************************************************************************
Extracting Gaussian constraints from a parameter-free distribution fit
***************************************************************************************

**Script:**


.. literalinclude:: ../../../examples/DL_extracting_gauss_constraints.m


**Output:**

.. code-block:: none

        Goodness of fit
          V{1}: χ2 = 1.0040  RMSD  = 0.0103 
        Fitted parameters and 95%-confidence intervals
          bg{1}(1):   0.0765482  (0.0094946, 0.1436018)  Decay Rate (us^-1)
          ex{1}(1):   0.3046161  (0.2664917, 0.3427405)  Modulation depth (  )

        Gaussian constraints:
          parfit(1) = 3.03 (2.95, 3.10) Center of 1st Gaussian
          parfit(2) = 0.40 (0.28, 0.52) FWHM of 1st Gaussian
          parfit(3) = 0.29 (0.19, 0.39) Amplitude of 1st Gaussian
          parfit(4) = 4.25 (4.17, 4.33) Center of 2nd Gaussian
          parfit(5) = 0.97 (0.75, 1.19) FWHM of 2nd Gaussian
          parfit(6) = 0.79 (0.75, 0.82) Amplitude of 2nd Gaussian


.. image:: ../images/example_extracting_gauss_constraints.png
   :width: 30%