Distance distribution models
=======================================

In DeerLab, distance distributions can be represented either as parameter-free models or as parametric models.

Parameter-free distance distributions
------------------------------------------------

A parameter-free distance distribution is a vector of population densities ``P`` defined over a vector of distances ``r``. Outside the range of ``r``, the distance distribution is considered zero, i.e. it is truncated to the range of ``r``. Such parameter-free distance distributions are returned by the fit functions :doc:`fitregmodel<api/fitregmodel>` and :doc:`obir<api/obir>`. A ``P`` can be converted to a time-domain signal using :doc:`dipolarsignal<api/dipolarsignal>` or :doc:`dipolarkernel<api/dipolarkernel>`.

Parameter-free distance distributions are preferred over parametric distance distributions, since they make fewer assumptions about the distribution and are more flexible. They introduce less bias.

Parametric distance distributions
------------------------------------------------

DeerLab provide a series of parametric distance distributions models. These can be used with several functions, including :doc:`fitparamodel<api/fitparamodel>` and :doc:`selectmodel<api/selectmodel>`. All parametric distance distribution models are functions start with the prefix ``dd_``. They take a distance vector ``r`` and a parameter vector ``param`` as inputs and return the distance distribution as a vector ``P``. Here is an example:

.. code-block:: matlab

        P = dd_onerice(r,params)


There are several classes of parametric distance distribution models. The first class contains models that are based on basis functions and linear combinations thereof.

=================================================  ===================================================
  Model function                                       Description
=================================================  ===================================================
:doc:`dd_onerice<models/dd_onerice>`               Single 3D-Rice distribution
:doc:`dd_tworice<models/dd_tworice>`               Two 3D-Rice distributions
:doc:`dd_threerice<models/dd_threerice>`           Three 3D-Rice distributions
:doc:`dd_onegauss<models/dd_onegauss>`             Single Gaussian distribution
:doc:`dd_twogauss<models/dd_twogauss>`             Two Gaussians distributions
:doc:`dd_threegauss<models/dd_threegauss>`         Three Gaussians distributions
:doc:`dd_fourgauss<models/dd_threegauss>`          Four Gaussians distributions
:doc:`dd_fivegauss<models/dd_threegauss>`          Five Gaussians distributions
:doc:`dd_gengauss<models/dd_gengauss>`             Single generalized Gaussian distribution
:doc:`dd_skewgauss<models/dd_skewgauss>`           Single skew Gaussian distribution
=================================================  ===================================================

The second group contains distance distribution models that represent three-dimensional disordered segmented objects such as proteins and other polymers.

=================================================  ===================================================
  Model function                                       Description
=================================================  ===================================================
:doc:`dd_wormchain<models/dd_wormchain>`           Worm-like chain
:doc:`dd_wormgauss<models/dd_wormgauss>`           Worm-like chain with Gaussian convolution
:doc:`dd_randcoil<models/dd_randcoil>`             Random coil
=================================================  ===================================================

Finally, DeerLab provides a series of models for distributions of labels in simple confined spaces such as spheres and spherical shells.

====================================================  =============================================================================
  Model function                                       Description
====================================================  =============================================================================
:doc:`dd_sphere<models/dd_sphere>`                    Labels distributed in a sphere
:doc:`dd_spheresurf<models/dd_spheresurf>`            Labels distributed on a sphere's surface
:doc:`dd_spherepoint<models/dd_spherepoint>`          One label fixed, other label distributed in a sphere
:doc:`dd_shell<models/dd_shell>`                      Labels distributed on a spherical shell
:doc:`dd_shellsphere<models/dd_shellsphere>`          Labels distributed on a sphere inside a spherical shell
:doc:`dd_shellshell<models/dd_shellshell>`            Labels distributed on a spherical shell inside another spherical shell
:doc:`dd_shellvoidshell<models/dd_shellvoidshell>`    Labels distributed on two spherical shells separated by a void
:doc:`dd_shellvoidsphere<models/dd_shellvoidsphere>`  Labels distributed on a sphere inside a spherical shell separated by a void 
====================================================  =============================================================================


.. toctree::
    :maxdepth: 1
    :hidden:
    :caption: Distance distribution models

    ./models/dd_onegauss
    ./models/dd_twogauss
    ./models/dd_threegauss
    ./models/dd_fourgauss
    ./models/dd_fivegauss
    ./models/dd_gengauss
    ./models/dd_skewgauss
    ./models/dd_onerice
    ./models/dd_tworice
    ./models/dd_threerice
    ./models/dd_wormchain
    ./models/dd_wormgauss
    ./models/dd_randcoil
    ./models/dd_shellshell
    ./models/dd_shell
    ./models/dd_shellsphere
    ./models/dd_sphere
    ./models/dd_spheresurf
    ./models/dd_spherepoint
    ./models/dd_shellvoidsphere
    ./models/dd_shellvoidshell

..
  .. raw:: html

   <br>
   <div class="tutorialSectionBox ", style="border: #5c87d6 1.5px solid; padding-bottom:15px;">
     <div class="sectionLogo"><img class="avatar" src="./_static/img/models_section_distance.png"></div>
     <div class="sectionTitle"><h3>Parametric distance distribution models</h3></div>
         <ul class="simple" style="width: 90%;margin: auto;">
         <li><a class="reference internal" href="models/dd_onegauss.html#dd_onegauss"><span class="std std-ref" style="font-family: monospace">dd_onegauss</span></a>
         &nbsp - &nbsp Unimodal Gaussian distribution </li>
         <li><a class="reference internal" href="models/dd_twogauss.html#dd_twogauss"><span class="std std-ref" style="font-family: monospace">dd_twogauss</span></a>
         &nbsp - &nbsp Bimodal Gaussian distribution </li>
         <li><a class="reference internal" href="models/dd_threegauss.html#dd_threegauss"><span class="std std-ref" style="font-family: monospace">dd_threegauss</span></a>
         &nbsp - &nbsp Trimodal Gaussian distribution </li>
         <li><a class="reference internal" href="models/dd_fourgauss.html#dd_fourgauss"><span class="std std-ref" style="font-family: monospace">dd_fourgauss</span></a>
         &nbsp - &nbsp Tetramodal Gaussian distribution </li>
         <li><a class="reference internal" href="models/dd_fivegauss.html#dd_fivegauss"><span class="std std-ref" style="font-family: monospace">dd_fivegauss</span></a>
         &nbsp - &nbsp Pentamodal Gaussian distribution </li>
         <li><a class="reference internal" href="models/dd_gengauss.html#dd_gengauss"><span class="std std-ref" style="font-family: monospace">dd_gengauss</span></a>
         &nbsp - &nbsp Generalized Gaussian distribution </li>
         <li><a class="reference internal" href="models/dd_skewgauss.html#dd_skewgauss"><span class="std std-ref" style="font-family: monospace">dd_skewgauss</span></a>
         &nbsp - &nbsp Skew Gaussian distribution </li>
         <li><a class="reference internal" href="models/dd_onerice.html#dd_onerice"><span class="std std-ref" style="font-family: monospace">dd_onerice</span></a>
         &nbsp - &nbsp Unimodal Rician distribution </li>
         <li><a class="reference internal" href="models/dd_tworice.html#dd_tworice"><span class="std std-ref" style="font-family: monospace">dd_tworice</span></a>
         &nbsp - &nbsp Bimodal Rician distribution </li>
         <li><a class="reference internal" href="models/dd_threerice.html#dd_threerice"><span class="std std-ref" style="font-family: monospace">dd_threerice</span></a>
         &nbsp - &nbsp Trimodal Rician distribution </li>
         <li><a class="reference internal" href="models/dd_wormchain.html#dd_wormchain"><span class="std std-ref" style="font-family: monospace">dd_wormchain</span></a>
         &nbsp - &nbsp Worm-like chain distribution </li>
         <li><a class="reference internal" href="models/dd_wormgauss.html#dd_wormgauss"><span class="std std-ref" style="font-family: monospace">dd_wormgauss</span></a>
         &nbsp - &nbsp Worm-like chain with Gaussian convolution distribution </li>
         <li><a class="reference internal" href="models/dd_randcoil.html#dd_randcoil"><span class="std std-ref" style="font-family: monospace">dd_randcoil</span></a>
         &nbsp - &nbsp Random coil distribution </li>
         <li><a class="reference internal" href="models/dd_sphere.html#dd_sphere"><span class="std std-ref" style="font-family: monospace">dd_sphere</span></a>
         &nbsp - &nbsp Particles distributed on a sphere</li>
         <li><a class="reference internal" href="models/dd_spheresurf.html#dd_spheresurf"><span class="std std-ref" style="font-family: monospace">dd_spheresurf</span></a>
         &nbsp - &nbsp Particles distributed on a sphere's surface</li>
         <li><a class="reference internal" href="models/dd_shell.html#dd_shell"><span class="std std-ref" style="font-family: monospace">dd_shell</span></a>
         &nbsp - &nbsp Particles distributed on a spherical shell</li>
         <li><a class="reference internal" href="models/dd_shellsphere.html#dd_shellsphere"><span class="std std-ref" style="font-family: monospace">dd_shellsphere</span></a>
         &nbsp - &nbsp Particles distributed on a sphere inside a spherical shell </li>
         <li><a class="reference internal" href="models/dd_shellshell.html#dd_shellshell"><span class="std std-ref" style="font-family: monospace">dd_shellshell</span></a>
         &nbsp - &nbsp Particles distributed on a spherical shell inside another spherical shell</li>
         <li><a class="reference internal" href="models/dd_spherepoint.html#dd_spherepoint"><span class="std std-ref" style="font-family: monospace">dd_spherepoint</span></a>
         &nbsp - &nbsp One particle distanced from particles distributed on a sphere </li>
         <li><a class="reference internal" href="models/dd_shellvoidsphere.html#dd_shellvoidsphere"><span class="std std-ref" style="font-family: monospace">dd_shellvoidsphere</span></a>
         &nbsp - &nbsp  Particles distributed on a sphere inside a spherical shell separated by a void </li>
         <li><a class="reference internal" href="models/dd_shellvoidshell.html#dd_shellvoidshell"><span class="std std-ref" style="font-family: monospace">dd_shellvoidshell</span></a>
         &nbsp - &nbsp  Particles distributed on two spherical shells separated by a void </li>
         </ul>
     </div>
   <br>