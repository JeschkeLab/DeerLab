Distribution models
=========================================

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


.. raw:: html

	<div class="limiter">
	<div class="container-table100">
	<div class="wrap-table100">
	<div class="table">


		<div class="row" onclick="window.location='models/dd_onerice.html#dd_onerice'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">dd_onerice</span>			</div>
			<div class="cell">
				Single 3D-Rice distribution
			</div>
		</div>
		

		<div class="row" onclick="window.location='models/dd_tworice.html#dd_tworice'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">dd_tworice</span>			</div>
			<div class="cell">
				Two 3D-Rice distributions
			</div>
		</div>

		<div class="row" onclick="window.location='models/dd_threerice.html#dd_threerice'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">dd_threerice</span>			</div>
			<div class="cell">
				Three 3D-Rice distributions
			</div>
		</div>

		<div class="row" onclick="window.location='models/dd_onegauss.html#dd_onegauss'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">dd_onegauss</span>			</div>
			<div class="cell">
				Single Gaussian distribution
			</div>
		</div>

		<div class="row" onclick="window.location='models/dd_twogauss.html#dd_twogauss'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">dd_twogauss</span>			</div>
			<div class="cell">
				Two Gaussians distributions
			</div>
		</div>

		<div class="row" onclick="window.location='models/dd_threegauss.html#dd_threegauss'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">dd_threegauss</span>			</div>
			<div class="cell">
				Three Gaussians distributions
			</div>
		</div>

		<div class="row" onclick="window.location='models/dd_fourgauss.html#dd_fourgauss'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">dd_fourgauss</span>			</div>
			<div class="cell">
				Four Gaussians distributions
			</div>
		</div>

		<div class="row" onclick="window.location='models/dd_fivegauss.html#dd_fivegauss'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">dd_fivegauss</span>			</div>
			<div class="cell">
				Five Gaussians distributions
			</div>
		</div>

		<div class="row" onclick="window.location='models/dd_gengauss.html#dd_gengauss'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">dd_gengauss</span>			</div>
			<div class="cell">
				Single generalized Gaussian distribution
			</div>
		</div>

		<div class="row" onclick="window.location='models/dd_skewgauss.html#dd_skewgauss'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">dd_skewgauss</span>			</div>
			<div class="cell">
				Single skew Gaussian distribution
			</div>
		</div>

	</div>
	</div>
	</div>
	</div>



The second group contains distance distribution models that represent three-dimensional disordered segmented objects such as proteins and other polymers.


.. raw:: html

	<div class="limiter">
	<div class="container-table100">
	<div class="wrap-table100">
	<div class="table">

		<div class="row" onclick="window.location='models/dd_wormchain.html#dd_wormchain'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">dd_wormchain</span>			</div>
			<div class="cell">
				Worm-like chain
			</div>
		</div>

		<div class="row" onclick="window.location='models/dd_wormgauss.html#dd_wormgauss'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">dd_wormgauss</span>			</div>
			<div class="cell">
				Worm-like chain with Gaussian convolution
			</div>
		</div>

		<div class="row" onclick="window.location='models/dd_randcoil.html#dd_randcoil'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">dd_randcoil</span>			</div>
			<div class="cell">
				Random coil
			</div>
		</div>

	</div>
	</div>
	</div>
	</div>


Finally, DeerLab provides a series of models for distributions of labels in simple confined spaces such as spheres and spherical shells.


.. raw:: html

	<div class="limiter">
	<div class="container-table100">
	<div class="wrap-table100">
	<div class="table">

		<div class="row" onclick="window.location='models/dd_sphere.html#dd_sphere'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">dd_sphere</span>			</div>
			<div class="cell">
				Labels distributed in a sphere
			</div>
		</div>

		<div class="row" onclick="window.location='models/dd_spheresurf.html#dd_spheresurf'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">dd_spheresurf</span>			</div>
			<div class="cell">
				Labels distributed on a sphere's surface
			</div>
		</div>

		<div class="row" onclick="window.location='models/dd_spherepoint.html#dd_spherepoint'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">dd_spherepoint</span>			</div>
			<div class="cell">
				One label fixed, other label distributed in a sphere
			</div>
		</div>

		<div class="row" onclick="window.location='models/dd_shell.html#dd_shell'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">dd_shell</span>			</div>
			<div class="cell">
				Labels distributed on a spherical shell
			</div>
		</div>

		<div class="row" onclick="window.location='models/dd_shellsphere.html#dd_shellsphere'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">dd_shellsphere</span>			</div>
			<div class="cell">
				Labels distributed on a sphere inside a spherical shell
			</div>
		</div>

		<div class="row" onclick="window.location='models/dd_shellshell.html#dd_shellshell'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">dd_shellshell</span>			</div>
			<div class="cell">
				Labels distributed on a spherical shell inside another spherical shell
			</div>
		</div>

		<div class="row" onclick="window.location='models/dd_shellvoidshell.html#dd_shellvoidshell'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">dd_shellvoidshell</span>			</div>
			<div class="cell">
				Labels distributed on two spherical shells separated by a void
			</div>
		</div>

		<div class="row" onclick="window.location='models/dd_shellvoidsphere.html#dd_shellvoidsphere'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">dd_shellvoidsphere</span>			</div>
			<div class="cell">
				Labels distributed on a sphere inside a spherical shell separated by a void
			</div>
		</div>


	</div>
	</div>
	</div>
	</div>


.. toctree::
    :maxdepth: 0
    :hidden:
    :glob:
    :caption: Background models

    ./models/dd_*
