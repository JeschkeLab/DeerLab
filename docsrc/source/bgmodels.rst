Background models
======================

DeerLab includes a collection of parametric models that can be used to model the background signal, i.e. the signal due to randomly distributed spins in the sample that are not part of the spin-labelled protein or object.

All background model functions start with the prefix ``bg_``. They take a time axis vector ``t`` (in microseconds) and a paramete vector ``params`` as inputs and return a background function ``B``, defined over ``t``, as output. Here is an example:

.. code-block:: matlab

        B = bg_strexp(t,params)


The first class of models involves exponential decays of various forms. The exponential and stretched exponential decays are physical models that represent homogeneous uniform distributions of spin labels.

=============================================  ============================================================
Model function                                 Description
=============================================  ============================================================
:doc:`bg_exp<models/bg_exp>`                   Exponential decay
:doc:`bg_strexp<models/bg_strexp>`             Stretched exponential decay
:doc:`bg_sumstrexp<models/bg_sumstrexp>`       Sum of two stretched exponential decays
:doc:`bg_prodstrexp<models/bg_prodstrexp>`     Product of two stretched exponential decays
=============================================  ============================================================

The second class contains polynomial functions. These are purely phenomenological and lack a physical basis.

=============================================  ============================================================
Model function                                 Description
=============================================  ============================================================
:doc:`bg_poly1<models/bg_poly1>`               Linear polynomial
:doc:`bg_poly2<models/bg_poly2>`               Quadratic polynomial
:doc:`bg_poly3<models/bg_poly3>`               Cubic polynomial
=============================================  ============================================================

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


..
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
