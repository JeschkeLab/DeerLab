.. _functions:


Functions
---------------------

This is the official documentation for the DeerLab toolbox functions. The following list contains the names of the different function and a brief description of their functionality. The parametric model functions are listed in a separate section (see the :ref:`parametric_models` section).


---------------------

Modelling
=========================================

.. toctree::
    :hidden:
    :glob:
    :maxdepth: 0

    ./api/dipolarkernel
    ./api/aptkernel
    ./api/dipolarsignal
    ./api/whitenoise
    ./api/paramodel
    ./api/mixmodels

.. raw:: html


	<div class="limiter">
	<div class="container-table100">
	<div class="wrap-table100">
	<div class="table">


		<div class="row" onclick="window.location='api/dipolarkernel.html#dipolarkernel'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">dipolarkernel</span>			</div>
			<div class="cell">
				Dipolar kernel contructor
			</div>
		</div>


		<div class="row" onclick="window.location='api/aptkernel.html#aptkernel'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">aptkernel</span>			</div>
			<div class="cell">
				APT kernel contructor
			</div>
		</div>


		<div class="row" onclick="window.location='api/dipolarsignal.html#dipolarsignal'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">dipolarsignal</span>			</div>
			<div class="cell">
				Dipolar signal simulator
			</div>
		</div>


		<div class="row" onclick="window.location='api/whitenoise.html#whitenoise'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">whitegaussnoise</span>			</div>
			<div class="cell">
				Gaussian white noise generator
			</div>
		</div>


		<div class="row" onclick="window.location='api/paramodel.html#paramodel'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">paramodel</span>			</div>
			<div class="cell">
				Parametric model builder
			</div>
		</div>


		<div class="row" onclick="window.location='api/mixmodels.html#mixmodels'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">mixmodels</span>			</div>
			<div class="cell">
				Parametric model mixer
			</div>
		</div>
	</div>
	</div>
	</div>
	</div>


---------------------

Analysis
=========================================


.. toctree::
    :hidden:
    :glob:
    :maxdepth: 0

    ./api/apt
    ./api/backgroundstart
    ./api/fitbackground
    ./api/fitmultigauss
    ./api/fitparamodel
    ./api/fitregmodel
    ./api/obir
    ./api/regoperator
    ./api/sensitivan

.. raw:: html

	<div class="limiter">
	<div class="container-table100">
	<div class="wrap-table100">
	<div class="table">


		<div class="row" onclick="window.location='api/backgroundstart.html#backgroundstart'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">backgroundstart</span>			</div>
			<div class="cell">
				Background fit start point optimizer
			</div>
		</div>


		<div class="row" onclick="window.location='api/fitbackground.html#fitbackground'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">fitbackground</span>			</div>
			<div class="cell">
				Background fitting engine
			</div>
		</div>


		<div class="row" onclick="window.location='api/fitmultigauss.html#fitmultigauss'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">fitmultigauss</span>			</div>
			<div class="cell">
				Multi-Gauss fitting engine
			</div>
		</div>


		<div class="row" onclick="window.location='api/fitparamodel.html#fitparamodel'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">fitparamodel</span>			</div>
			<div class="cell">
				Parametric model fitting engine
			</div>
		</div>


		<div class="row" onclick="window.location='api/fitregmodel.html#fitregmodel'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">fitregmodel</span>			</div>
			<div class="cell">
				Regularization fitting engine
			</div>
		</div>


		<div class="row" onclick="window.location='api/obir.html#obir'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">obir</span>			</div>
			<div class="cell">
				Osher-Bregman iterative regularization
			</div>
		</div>

		<div class="row" onclick="window.location='api/regoperator.html#regoperator'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">regoperator</span>			</div>
			<div class="cell">
				Regularization operator constructor
			</div>
		</div>

		<div class="row" onclick="window.location='api/apt.html#apt'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">apt</span>			</div>
			<div class="cell">
				Approximate Pake transformation
			</div>
		</div>


		<div class="row" onclick="window.location='api/sensitivan.html#sensitivan'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">sensitivan</span>			</div>
			<div class="cell">
				Sensitivity analysis engine
			</div>
		</div>
	</div>
	</div>
	</div>
	</div>

----------------

Pre-Processing
=========================================


.. toctree::
    :hidden:
    :glob:
    :maxdepth: 0

    ./api/correctphase
    ./api/correctzerotime
    ./api/correctscale
    ./api/suppressghost
    ./api/longpass
    ./api/winlowpass


.. raw:: html

	<div class="limiter">
	<div class="container-table100">
	<div class="wrap-table100">
	<div class="table">


		<div class="row" onclick="window.location='api/correctphase.html#correctphase'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">correctphase</span>			</div>
			<div class="cell">
				IQ Phase correction
			</div>
		</div>


		<div class="row" onclick="window.location='api/correctzerotime.html#correctzerotime'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">correctzerotime</span>			</div>
			<div class="cell">
				Dipolar zero-time correction
			</div>
		</div>


		<div class="row" onclick="window.location='api/correctscale.html#correctscale'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">correctscale</span>			</div>
			<div class="cell">
				Dipolar signal amplitude rescaling
			</div>
		</div>


		<div class="row" onclick="window.location='api/suppressghost.html#suppressghost'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">suppressghost</span>			</div>
			<div class="cell">
				Ghost-distance suppression
			</div>
		</div>


		<div class="row" onclick="window.location='api/longpass.html#longpass'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">longpass</span>			</div>
			<div class="cell">
				Longpass filtering
			</div>
		</div>


		<div class="row" onclick="window.location='api/winlowpass.html#winlowpass'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">winlowpass</span>			</div>
			<div class="cell">
				Windowed-lowpass filtering
			</div>
		</div>
	</div>
	</div>
	</div>
	</div>

------------------

Model Selection
=========================================


.. toctree::
    :hidden:
    :glob:
    :maxdepth: 0

    ./api/selectmodel
    ./api/selregparam
    ./api/regparamrange


.. raw:: html

	<div class="limiter">
	<div class="container-table100">
	<div class="wrap-table100">
	<div class="table">


		<div class="row" onclick="window.location='api/selectmodel.html#selectmodel'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">selectmodel</span>			</div>
			<div class="cell">
				Parametric model selector
			</div>
		</div>


		<div class="row" onclick="window.location='api/selregparam.html#selregparam'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">selregparam</span>			</div>
			<div class="cell">
				Regularization parameter selector
			</div>
		</div>


		<div class="row" onclick="window.location='api/regparamrange.html#regparamrange'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">regparamrange</span>			</div>
			<div class="cell">
				Regularization parameter range selector
			</div>
		</div>


	</div>
	</div>
	</div>
	</div>


------------------


Utilities
=========================================


.. toctree::
    :hidden:
    :glob:
    :maxdepth: 0

    ./api/deerload
    ./api/time2freq
    ./api/time2dist
    ./api/noiselevel
    ./api/fftspec
    ./api/prepvalidation

.. raw:: html

	<div class="limiter">
	<div class="container-table100">
	<div class="wrap-table100">
	<div class="table">

		<div class="row" onclick="window.location='api/deerload.html#deerload'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">deerload</span>			</div>
			<div class="cell">
				Spectrometer data loader
			</div>
		</div>

		<div class="row" onclick="window.location='api/time2freq.html#time2freq'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">time2freq</span>			</div>
			<div class="cell">
				Time to frequency axis convertor
			</div>
		</div>

		<div class="row" onclick="window.location='api/time2dist.html#time2dist'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">time2dist</span>			</div>
			<div class="cell">
				Time to distance axis convertor
			</div>
		</div>

		<div class="row" onclick="window.location='api/noiselevel.html#noiselevel'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">noiselevel</span>			</div>
			<div class="cell">
				Noise level estimator
			</div>
		</div>

		<div class="row" onclick="window.location='api/fftspec.html#fftspec'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">fftspec</span>			</div>
			<div class="cell">
				Fast-Fourier trasnform spectrum
			</div>
		</div>


		<div class="row" onclick="window.location='api/prepvalidation.html#prepvalidation'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">prepvalidation</span>			</div>
			<div class="cell">
				Full-factorial analysis preparation
			</div>
		</div>
	</div>
	</div>
	</div>
	</div>

