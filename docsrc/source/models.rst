.. _parametric_models:

Parametric Models
======================

DeerLab includes the following collection of parametric models. The model names are categorized depending on whether the model is defined in time-domain (model name starts with the prefix ``td_``) or in distance model (model name starts with the prefix ``rd_``). 

.. toctree::
    :maxdepth: 1
    :hidden:
    :caption: Distance-domain

    ./models/rd_onegaussian
    ./models/rd_twogaussian
    ./models/rd_threegaussian
    ./models/rd_fourgaussian
    ./models/rd_fivegaussian
    ./models/rd_gengaussian
    ./models/rd_skewgaussian
    ./models/rd_onerice
    ./models/rd_tworice
    ./models/rd_threerice
    ./models/rd_wormchain
    ./models/rd_wormgauss
    ./models/rd_randcoil
    ./models/rd_shellshell
    ./models/rd_shell
    ./models/rd_shellsphere
    ./models/rd_sphere
    ./models/rd_spheresurf
    ./models/rd_spherepoint
    ./models/rd_shellvoidsphere
    ./models/rd_shellvoidshell


.. toctree::
    :maxdepth: 1
    :hidden:
    :caption: Time-domain

    ./models/td_exp
    ./models/td_strexp
    ./models/td_sumstrexp
    ./models/td_prodstrexp
    ./models/td_poly1
    ./models/td_poly2
    ./models/td_poly3


.. raw:: html

   <br>
   <div class="tutorialSectionBox ", style="border: #5c87d6 1.5px solid; padding-bottom:15px;">
     <div class="sectionLogo"><img class="avatar" src="./_static/img/models_section_time.png"></div>
     <div class="sectionTitle"><h3>Time-domain models</h3></div>
         <ul class="simple" style="width: 90%;margin: auto;">
         <li><a class="reference internal" href="models/td_exp.html#td_exp"><span class="std std-ref" style="font-family: monospace">td_exp</span></a>
         &nbsp - &nbsp Exponential function </li>
         <li><a class="reference internal" href="models/td_strexp.html#td_strexp"><span class="std std-ref" style="font-family: monospace">td_strexp</span></a>
         &nbsp - &nbsp Stretched exponential function </li>
         <li><a class="reference internal" href="models/td_sumstrexp.html#td_sumstrexp"><span class="std std-ref" style="font-family: monospace">td_sumstrexp</span></a>
         &nbsp - &nbsp Sum of two stretched exponentials </li>
         <li><a class="reference internal" href="models/td_prodstrexp.html#td_prodstrexp"><span class="std std-ref" style="font-family: monospace">td_prodstrexp</span></a>
         &nbsp - &nbsp Product of two stretched exponentials </li>
         <li><a class="reference internal" href="models/td_poly1.html#td_poly1"><span class="std std-ref" style="font-family: monospace">td_poly1</span></a>
         &nbsp - &nbsp Polynomial 1st order </li>
         <li><a class="reference internal" href="models/td_poly2.html#td_poly2"><span class="std std-ref" style="font-family: monospace">td_poly2</span></a>
         &nbsp - &nbsp Polynomial 2nd order </li>
         <li><a class="reference internal" href="models/td_poly3.html#td_poly3"><span class="std std-ref" style="font-family: monospace">td_poly3</span></a>
         &nbsp - &nbsp Polynomial 3rd order </li>
         </ul>
     </div>

.. raw:: html

   <br>
   <div class="tutorialSectionBox ", style="border: #5c87d6 1.5px solid; padding-bottom:15px;">
     <div class="sectionLogo"><img class="avatar" src="./_static/img/models_section_distance.png"></div>
     <div class="sectionTitle"><h3>Distance-domain models</h3></div>
         <ul class="simple" style="width: 90%;margin: auto;">
         <li><a class="reference internal" href="models/rd_onegaussian.html#rd_onegaussian"><span class="std std-ref" style="font-family: monospace">rd_onegaussian</span></a>
         &nbsp - &nbsp Unimodal Gaussian distribution </li>
         <li><a class="reference internal" href="models/rd_twogaussian.html#rd_twogaussian"><span class="std std-ref" style="font-family: monospace">rd_twogaussian</span></a>
         &nbsp - &nbsp Bimodal Gaussian distribution </li>
         <li><a class="reference internal" href="models/rd_threegaussian.html#rd_threegaussian"><span class="std std-ref" style="font-family: monospace">rd_threegaussian</span></a>
         &nbsp - &nbsp Trimodal Gaussian distribution </li>
         <li><a class="reference internal" href="models/rd_fourgaussian.html#rd_fourgaussian"><span class="std std-ref" style="font-family: monospace">rd_fourgaussian</span></a>
         &nbsp - &nbsp Tetramodal Gaussian distribution </li>
         <li><a class="reference internal" href="models/rd_fivegaussian.html#rd_fivegaussian"><span class="std std-ref" style="font-family: monospace">rd_fivegaussian</span></a>
         &nbsp - &nbsp Tetramodal Gaussian distribution </li>
         <li><a class="reference internal" href="models/rd_gengaussian.html#rd_gengaussian"><span class="std std-ref" style="font-family: monospace">rd_gengaussian</span></a>
         &nbsp - &nbsp Generalized Gaussian distribution </li>
         <li><a class="reference internal" href="models/rd_skewgaussian.html#rd_skewgaussian"><span class="std std-ref" style="font-family: monospace">rd_skewgaussian</span></a>
         &nbsp - &nbsp Skew Gaussian distribution </li>
         <li><a class="reference internal" href="models/rd_onerice.html#rd_onerice"><span class="std std-ref" style="font-family: monospace">rd_onerice</span></a>
         &nbsp - &nbsp Unimodal Rician distribution </li>
         <li><a class="reference internal" href="models/rd_tworice.html#rd_tworice"><span class="std std-ref" style="font-family: monospace">rd_tworice</span></a>
         &nbsp - &nbsp Bimodal Rician distribution </li>
         <li><a class="reference internal" href="models/rd_threerice.html#rd_threerice"><span class="std std-ref" style="font-family: monospace">rd_threerice</span></a>
         &nbsp - &nbsp Trimodal Rician distribution </li>
         <li><a class="reference internal" href="models/rd_wormchain.html#rd_wormchain"><span class="std std-ref" style="font-family: monospace">rd_wormchain</span></a>
         &nbsp - &nbsp Worm-like chain (WLC) distribution </li>
         <li><a class="reference internal" href="models/rd_wormgauss.html#rd_wormgauss"><span class="std std-ref" style="font-family: monospace">rd_wormgauss</span></a>
         &nbsp - &nbsp Worm-like chain (WLC) with Gaussian convolution distribution </li>
         <li><a class="reference internal" href="models/rd_randcoil.html#rd_randcoil"><span class="std std-ref" style="font-family: monospace">rd_randcoil</span></a>
         &nbsp - &nbsp Random coil distribution </li>
         <li><a class="reference internal" href="models/rd_sphere.html#rd_sphere"><span class="std std-ref" style="font-family: monospace">rd_sphere</span></a>
         &nbsp - &nbsp Particles distributed on a sphere</li>
         <li><a class="reference internal" href="models/rd_spheresurf.html#rd_spheresurf"><span class="std std-ref" style="font-family: monospace">rd_spheresurf</span></a>
         &nbsp - &nbsp Particles distributed on a sphere's surface</li>
         <li><a class="reference internal" href="models/rd_shell.html#rd_shell"><span class="std std-ref" style="font-family: monospace">rd_shell</span></a>
         &nbsp - &nbsp Particles distributed on a spherical shell</li>
         <li><a class="reference internal" href="models/rd_shellsphere.html#rd_shellsphere"><span class="std std-ref" style="font-family: monospace">rd_shellsphere</span></a>
         &nbsp - &nbsp Particles distributed on a sphere inside a spherical shell </li>
         <li><a class="reference internal" href="models/rd_shellshell.html#rd_shellshell"><span class="std std-ref" style="font-family: monospace">rd_shellshell</span></a>
         &nbsp - &nbsp Particles distributed on a spherical shell inside another spherical shell</li>
         <li><a class="reference internal" href="models/rd_spherepoint.html#rd_spherepoint"><span class="std std-ref" style="font-family: monospace">rd_spherepoint</span></a>
         &nbsp - &nbsp One particle distanced from particles distributed on a sphere </li>
         <li><a class="reference internal" href="models/rd_shellvoidsphere.html#rd_shellvoidsphere"><span class="std std-ref" style="font-family: monospace">rd_shellvoidsphere</span></a>
         &nbsp - &nbsp  Particles distributed on a sphere inside a spherical shell separated by a void </li>
         <li><a class="reference internal" href="models/rd_shellvoidshell.html#rd_shellvoidshell"><span class="std std-ref" style="font-family: monospace">rd_shellvoidshell</span></a>
         &nbsp - &nbsp  Particles distributed on two spherical shells separated by a void </li>
         </ul>
     </div>
   <br>