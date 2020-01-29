Functions
======================

This is the official documentation for the DeerAnalysis toolbox functions. The following list contains the names of the different function and a brief description of their functionality. The parametric model functions are listed in a separate section (see :ref:`parametric_models`).

.. toctree::
    :hidden:
    :maxdepth: 0

    ./api/apt
    ./api/aptkernel
    ./api/backgroundstart
    ./api/correctphase
    ./api/correctscale
    ./api/correctzerotime
    ./api/deerload
    ./api/dipolarkernel
    ./api/dipolarsignal
    ./api/fftspec
    ./api/fitbackground
    ./api/fitmultigauss
    ./api/fitparamodel
    ./api/fitregmodel
    ./api/longpass
    ./api/mixmodels
    ./api/noiselevel
    ./api/obir
    ./api/overtones
    ./api/paramodel
    ./api/prepvalidation
    ./api/regoperator
    ./api/regparamrange
    ./api/selectmodel
    ./api/selregparam
    ./api/sensitivan
    ./api/suppressghost
    ./api/time2dist
    ./api/time2freq
    ./api/whitenoise
    ./api/winlowpass
    ./api/docs

.. raw:: html

   <br>
   <div class="tutorialSectionBox ", style="border: #5c87d6 1.5px solid; padding-bottom:15px;">
     <div class="sectionLogo"><img class="avatar" src="./_static/img/function_section_model.png"></div>
     <div class="sectionTitle"><h3>Modelling</h3></div>
         <ul class="functionsList" style="width: 90%;margin: auto;">
         <li><a class="reference internal" href="api/dipolarkernel.html#dipolarkernel"><span class="std std-ref" style="font-family: monospace;">dipolarkernel</span></a>
         &nbsp - &nbsp Dipolar kernel contructor </li>
         <li><a class="reference internal" href="api/aptkernel.html#aptkernel"><span class="std std-ref" style="font-family: monospace">aptkernel</span></a>
         &nbsp - &nbsp APT kernel contructor </li>
         <li><a class="reference internal" href="api/dipolarsignal.html#dipolarsignal"><span class="std std-ref" style="font-family: monospace">dipolarsignal</span></a>
         &nbsp - &nbsp Dipolar signal simulator </li>
         <li><a class="reference internal" href="api/whitenoise.html#whitenoise"><span class="std std-ref" style="font-family: monospace">whitegaussnoise</span></a>
         &nbsp - &nbsp Gaussian white noise generator </li>
         <li><a class="reference internal" href="api/mixmodels.html#mixmodels"><span class="std std-ref" style="font-family: monospace">mixmodels</span></a>
         &nbsp - &nbsp Parametric model mixer </li>
         </ul>
   </div>
   <br>


.. raw:: html

   <br>
   <div class="tutorialSectionBox ", style="border: #5c87d6 1.5px solid; padding-bottom:15px;">
     <div class="sectionLogo"><img class="avatar" src="./_static/img/function_section_analysis.png"></div>
     <div class="sectionTitle"><h3>Analysis</h3></div>
         <ul class="functionsList" style="width: 90%;margin: auto;">
         <li><a class="reference internal" href="api/backgroundstart.html#backgroundstart"><span class="std std-ref" style="font-family: monospace">backgroundstart</span></a>
         &nbsp - &nbsp Background fit start point optimizer </li>
         <li><a class="reference internal" href="api/fitbackground.html#fitbackground"><span class="std std-ref" style="font-family: monospace">fitbackground</span></a>
         &nbsp - &nbsp Background fitting engine </li>
         <li><a class="reference internal" href="api/fitmultigauss.html#fitmultigauss"><span class="std std-ref" style="font-family: monospace">fitmultigauss</span></a>
         &nbsp - &nbsp Multi-Gauss fitting engine </li>
         <li><a class="reference internal" href="api/fitparamodel.html#fitparamodel"><span class="std std-ref" style="font-family: monospace">fitparamodel</span></a>
         &nbsp - &nbsp Parametric model fitting engine </li>
         <li><a class="reference internal" href="api/fitregmodel.html#fitregmodel"><span class="std std-ref" style="font-family: monospace">fitregmodel</span></a>
         &nbsp - &nbsp Regularization fitting engine </li>
         <li><a class="reference internal" href="api/obir.html#obir"><span class="std std-ref" style="font-family: monospace">obir</span></a>
         &nbsp - &nbsp Osher-Bregman iterative regularization </li>
         <li><a class="reference internal" href="api/apt.html#apt"><span class="std std-ref" style="font-family: monospace">apt</span></a>
         &nbsp - &nbsp Approximate Pake transformation </li>
         <li><a class="reference internal" href="api/sensitivan.html#sensitivan"><span class="std std-ref" style="font-family: monospace">sensitivan</span></a>
         &nbsp - &nbsp Sensitivity analysis engine </li>
         </ul>
     </div>

.. raw:: html

   <br>
   <div class="tutorialSectionBox ", style="border: #5c87d6 1.5px solid; padding-bottom:15px;">
     <div class="sectionLogo"><img class="avatar" src="./_static/img/function_section_preprocessing.png"></div>
     <div class="sectionTitle"><h3>Pre-Processing</h3></div>
         <ul class="functionsList" style="width: 90%;margin: auto;">
         <li><a class="reference internal" href="api/correctphase.html#correctphase"><span class="std std-ref" style="font-family: monospace">correctphase</span></a>
         &nbsp - &nbsp IQ Phase correction </li>
         <li><a class="reference internal" href="api/correctzerotime.html#correctzerotime"><span class="std std-ref" style="font-family: monospace">correctzerotime</span></a>
         &nbsp - &nbsp Dipolar zero-time correction </li>
         <li><a class="reference internal" href="api/correctscale.html#correctscale" style="font-family: monospace">correctscale</a>
         &nbsp - &nbsp Dipolar signal amplitude rescaling </li>
         <li><a class="reference internal" href="api/suppressghost.html#suppressghost"><span class="std std-ref" style="font-family: monospace">suppressghost</span></a>
         &nbsp - &nbsp Ghost-distance suppression </li>
         <li><a class="reference internal" href="api/longpass.html#longpass"><span class="std std-ref" style="font-family: monospace">longpass</span></a>
         &nbsp - &nbsp Longpass filtering </li>
         <li><a class="reference internal" href="api/winlowpass.html#winlowpass"><span class="std std-ref" style="font-family: monospace">winlowpass</span></a>
         &nbsp - &nbsp Windowed-lowpass filtering </li>
         </ul>
     </div>
   

.. raw:: html

   <br>
   <div class="tutorialSectionBox ", style="border: #5c87d6 1.5px solid; padding-bottom:15px;">
     <div class="sectionLogo"><img class="avatar" src="./_static/img/function_section_selection.png"></div>
     <div class="sectionTitle"><h3>Model Selection</h3></div>
         <ul class="functionsList" style="width: 90%;margin: auto;">
         <li><a class="reference internal" href="api/selectmodel.html#selectmodel"><span class="std std-ref" style="font-family: monospace">selectmodel</span></a>
         &nbsp - &nbsp Parametric model selector </li>
         <li><a class="reference internal" href="api/selregparam.html#selregparam"><span class="std std-ref" style="font-family: monospace">selregparam</span></a>
         &nbsp - &nbsp Regularization parameter selector </li>
         <li><a class="reference internal" href="api/regparamrange.html#regparamrange"><span class="std std-ref" style="font-family: monospace">regparamrange</span></a>
         &nbsp - &nbsp Regularization parameter range selector </li>
         </ul>
     </div>


.. raw:: html

   <br>
   <div class="tutorialSectionBox ", style="border: #5c87d6 1.5px solid; padding-bottom:15px;">
     <div class="sectionLogo"><img class="avatar" src="./_static/img/function_section_utilities.png"></div>
     <div class="sectionTitle"><h3>Utilities</h3></div>
         <ul class="functionsList" style="width: 90%;margin: auto;">
         <li><a class="reference internal" href="api/deerload.html#deerload"><span class="std std-ref" style="font-family: monospace">deerload</span></a>
         &nbsp - &nbsp Spectrometer data loader </li>
         <li><a class="reference internal" href="api/time2freq.html#time2freq"><span class="std std-ref" style="font-family: monospace">time2freq</span></a>
         &nbsp - &nbsp Time to frequency axis convertor </li>
         <li><a class="reference internal" href="api/time2dist.html#time2dist"><span class="std std-ref" style="font-family: monospace">time2dist</span></a>
         &nbsp - &nbsp Time to distance axis convertor </li>
         <li><a class="reference internal" href="api/noiselevel.html#noiselevel"><span class="std std-ref" style="font-family: monospace">noiselevel</span></a>
         &nbsp - &nbsp Noise level estimator </li>
         <li><a class="reference internal" href="api/fftspec.html#fftspec"><span class="std std-ref" style="font-family: monospace">fftspec</span></a>
         &nbsp - &nbsp Fast-Fourier trasnform spectrum </li>
         <li><a class="reference internal" href="api/docs.html#docs"><span class="std std-ref" style="font-family: monospace">docs</span></a>
         &nbsp - &nbsp Documentation caller </li>
         </ul>
   </div>
   <br>
   
