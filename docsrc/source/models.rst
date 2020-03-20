.. _parametric_models:

Parametric Models
======================

DeerLab includes the following collection of parametric models. The model names are categorized depending on whether the model is defined in time-domain (model name starts with the prefix ``bg_``) or in distance model (model name starts with the prefix ``dd_``). 

.. toctree::
    :maxdepth: 1
    :hidden:
    :caption: Distance-distribution models

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


.. toctree::
    :maxdepth: 1
    :hidden:
    :caption: Background models

    ./models/bg_exp
    ./models/bg_strexp
    ./models/bg_sumstrexp
    ./models/bg_prodstrexp
    ./models/bg_poly1
    ./models/bg_poly2
    ./models/bg_poly3


.. raw:: html

   <br>
   <div class="tutorialSectionBox ", style="border: #5c87d6 1.5px solid; padding-bottom:15px;">
     <div class="sectionLogo"><img class="avatar" src="./_static/img/models_section_time.png"></div>
     <div class="sectionTitle"><h3>Background models</h3></div>
         <ul class="simple" style="width: 90%;margin: auto;">
         <li><a class="reference internal" href="models/bg_exp.html#bg_exp"><span class="std std-ref" style="font-family: monospace">bg_exp</span></a>
         &nbsp - &nbsp Exponential function </li>
         <li><a class="reference internal" href="models/bg_strexp.html#bg_strexp"><span class="std std-ref" style="font-family: monospace">bg_strexp</span></a>
         &nbsp - &nbsp Stretched exponential function </li>
         <li><a class="reference internal" href="models/bg_sumstrexp.html#bg_sumstrexp"><span class="std std-ref" style="font-family: monospace">bg_sumstrexp</span></a>
         &nbsp - &nbsp Sum of two stretched exponentials </li>
         <li><a class="reference internal" href="models/bg_prodstrexp.html#bg_prodstrexp"><span class="std std-ref" style="font-family: monospace">bg_prodstrexp</span></a>
         &nbsp - &nbsp Product of two stretched exponentials </li>
         <li><a class="reference internal" href="models/bg_poly1.html#bg_poly1"><span class="std std-ref" style="font-family: monospace">bg_poly1</span></a>
         &nbsp - &nbsp Polynomial 1st order </li>
         <li><a class="reference internal" href="models/bg_poly2.html#bg_poly2"><span class="std std-ref" style="font-family: monospace">bg_poly2</span></a>
         &nbsp - &nbsp Polynomial 2nd order </li>
         <li><a class="reference internal" href="models/bg_poly3.html#bg_poly3"><span class="std std-ref" style="font-family: monospace">bg_poly3</span></a>
         &nbsp - &nbsp Polynomial 3rd order </li>
         </ul>
     </div>

.. raw:: html

   <br>
   <div class="tutorialSectionBox ", style="border: #5c87d6 1.5px solid; padding-bottom:15px;">
     <div class="sectionLogo"><img class="avatar" src="./_static/img/models_section_distance.png"></div>
     <div class="sectionTitle"><h3>Distance distribution models</h3></div>
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
         &nbsp - &nbsp Worm-like chain (WLC) distribution </li>
         <li><a class="reference internal" href="models/dd_wormgauss.html#dd_wormgauss"><span class="std std-ref" style="font-family: monospace">dd_wormgauss</span></a>
         &nbsp - &nbsp Worm-like chain (WLC) with Gaussian convolution distribution </li>
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