Background models
======================

DeerLab includes a collection of parametric models that can be used to model the background signal, i.e. the signal due to randomly distributed spins in the sample that are not part of the spin-labelled protein or object.

All background model functions start with the prefix ``bg_``. They take a time axis vector ``t`` (in microseconds) and a paramete vector ``params`` as inputs and return a background function ``B``, defined over ``t``, as output. Here is an example:

.. code-block:: matlab

        B = bg_strexp(t,params)



The first class of models involves exponential decays of various forms. The exponential and stretched exponential decays are physical models that represent homogeneous uniform distributions of spin labels.


.. raw:: html

	<div class="limiter">
	<div class="container-table100">
	<div class="wrap-table100">
	<div class="table">

		<div class="row" onclick="window.location='models/bg_exp.html#bg_exp'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">bg_exp</span>			</div>
			<div class="cell">
				Exponential decay
			</div>
		</div>


		<div class="row" onclick="window.location='models/bg_strexp.html#bg_strexp'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">bg_strexp</span>			</div>
			<div class="cell">
				Stretched exponential decay
			</div>
		</div>


		<div class="row" onclick="window.location='models/bg_sumstrexp.html#bg_sumstrexp'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">bg_sumstrexp</span>			</div>
			<div class="cell">
				Sum of two stretched exponential decays
			</div>
		</div>


		<div class="row" onclick="window.location='models/bg_prodstrexp.html#bg_prodstrexp'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">bg_prodstrexp</span>			</div>
			<div class="cell">
				Product of two stretched exponential decays
			</div>
		</div>

	</div>
	</div>
	</div>
	</div>


The second class contains polynomial functions. These are purely phenomenological and lack a physical basis.


.. raw:: html

	<div class="limiter">
	<div class="container-table100">
	<div class="wrap-table100">
	<div class="table">


		<div class="row" onclick="window.location='models/bg_poly1.html#bg_poly1'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">bg_poly1</span>			</div>
			<div class="cell">
				Linear polynomial
			</div>
		</div>


		<div class="row" onclick="window.location='models/bg_poly2.html#bg_poly2'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">bg_poly2</span>			</div>
			<div class="cell">
				Quadratic polynomial
			</div>
		</div>


		<div class="row" onclick="window.location='models/bg_poly3.html#bg_poly3'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">bg_poly3</span>			</div>
			<div class="cell">
				Cubic polynomial
			</div>
		</div>

	</div>
	</div>
	</div>
	</div>


.. toctree::
    :maxdepth: 0
    :hidden:
    :includehidden:
    :glob:
    :caption: Background models

    ./models/bg_*