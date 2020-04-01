.. _parametric_models:

Models
---------------------

DeerLab provides a wide colltection of parametric models for the simulation and fitting of different elements: distance distributions, background parameters and experiment properties. Each of these models is labeled with a different prefix according to the model family they belong to. 
Specific information on the model, their parameters, and mathematical expressions can be accessed by selecting the desired model from these lists.



------------------------


Distribution models
=========================================


.. toctree::
    :maxdepth: 0
    :hidden:
    :glob:
    :caption: Distribution models

    ./models/dd_*

Distance distribution models aim simulate the distirbutions of possible distances between a pair of unpaired electrons via analytical equations. 

.. raw:: html

	<div class="limiter">
	<div class="container-table100">
	<div class="wrap-table100">
	<div class="table">

		<div class="row" onclick="window.location='models/dd_onegauss.html#dd_onegauss'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">dd_onegauss</span>			</div>
			<div class="cell">
				Unimodal Gaussian distribution
			</div>
		</div>


		<div class="row" onclick="window.location='models/dd_twogauss.html#dd_twogauss'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">dd_twogauss</span>			</div>
			<div class="cell">
				Bimodal Gaussian distribution
			</div>
		</div>


		<div class="row" onclick="window.location='models/dd_threegauss.html#dd_threegauss'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">dd_threegauss</span>			</div>
			<div class="cell">
				Trimodal Gaussian distribution
			</div>
		</div>


		<div class="row" onclick="window.location='models/dd_fourgauss.html#dd_fourgauss'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">dd_fourgauss</span>			</div>
			<div class="cell">
				Tetramodal Gaussian distribution
			</div>
		</div>


		<div class="row" onclick="window.location='models/dd_fivegauss.html#dd_fivegauss'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">dd_fivegauss</span>			</div>
			<div class="cell">
				Pentamodal Gaussian distribution
			</div>
		</div>


		<div class="row" onclick="window.location='models/dd_gengauss.html#dd_gengauss'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">dd_gengauss</span>			</div>
			<div class="cell">
				Generalized Gaussian distribution
			</div>
		</div>


		<div class="row" onclick="window.location='models/dd_skewgauss.html#dd_skewgauss'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">dd_skewgauss</span>			</div>
			<div class="cell">
				Skew Gaussian distribution
			</div>
		</div>


		<div class="row" onclick="window.location='models/dd_onerice.html#dd_onerice'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">dd_onerice</span>			</div>
			<div class="cell">
				Unimodal Rician distribution
			</div>
		</div>
		

		<div class="row" onclick="window.location='models/dd_tworice.html#dd_tworice'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">dd_tworice</span>			</div>
			<div class="cell">
				Bimodal Rician distribution
			</div>
		</div>

		<div class="row" onclick="window.location='models/dd_threerice.html#dd_threerice'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">dd_threerice</span>			</div>
			<div class="cell">
				Trimodal Rician distribution
			</div>
		</div>


		<div class="row" onclick="window.location='models/dd_wormchain.html#dd_wormchain'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">dd_wormchain</span>			</div>
			<div class="cell">
				Worm-like chain (WLC) distribution
			</div>
		</div>


		<div class="row" onclick="window.location='models/dd_wormgauss.html#dd_wormgauss'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">dd_wormgauss</span>			</div>
			<div class="cell">
				Worm-like chain (WLC) with Gaussian convolution distribution
			</div>
		</div>


		<div class="row" onclick="window.location='models/dd_randcoil.html#dd_randcoil'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">dd_randcoil</span>			</div>
			<div class="cell">
				Random coil distribution
			</div>
		</div>


		<div class="row" onclick="window.location='models/dd_sphere.html#dd_sphere'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">dd_sphere</span>			</div>
			<div class="cell">
				Particles distributed on a sphere
			</div>
		</div>


		<div class="row" onclick="window.location='models/dd_spheresurf.html#dd_spheresurf'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">dd_spheresurf</span>			</div>
			<div class="cell">
				Particles distributed on a sphere's surface
			</div>
		</div>


		<div class="row" onclick="window.location='models/dd_shell.html#dd_shell'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">dd_shell</span>			</div>
			<div class="cell">
				Particles distributed on a spherical shell
			</div>
		</div>


		<div class="row" onclick="window.location='models/dd_shellsphere.html#dd_shellsphere'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">dd_shellsphere</span>			</div>
			<div class="cell">
				Particles distributed on a sphere inside a spherical shell
			</div>
		</div>


		<div class="row" onclick="window.location='models/dd_shellshell.html#dd_shellshell'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">dd_shellshell</span>			</div>
			<div class="cell">
				Particles distributed on a spherical shell inside another spherical shell
			</div>
		</div>


		<div class="row" onclick="window.location='models/dd_spherepoint.html#dd_spherepoint'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">dd_spherepoint</span>			</div>
			<div class="cell">
				One particle distanced from particles distributed on a sphere
			</div>
		</div>


		<div class="row" onclick="window.location='models/dd_shellvoidsphere.html#dd_shellvoidsphere'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">dd_shellvoidsphere</span>			</div>
			<div class="cell">
				Particles distributed on a sphere inside a spherical shell separated by a void
			</div>
		</div>


		<div class="row" onclick="window.location='models/dd_shellvoidshell.html#dd_shellvoidshell'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">dd_shellvoidshell</span>			</div>
			<div class="cell">
				Particles distributed on two spherical shells separated by a void
			</div>
		</div>


	</div>
	</div>
	</div>
	</div>



---------------------

Background models
=========================================

Backgound models aim simulate the inter-molecular interactions between multiple pairs of unpaired electrons, which lead to a decaying signal which is detected along the desired foreground diplar signal.

.. toctree::
    :maxdepth: 0
    :hidden:
    :includehidden:
    :glob:
    :caption: Background models

    ./models/bg_*

.. raw:: html

	<div class="limiter">
	<div class="container-table100">
	<div class="wrap-table100">
	<div class="table">


		<div class="row" onclick="window.location='models/bg_exp.html#bg_exp'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">bg_exp</span>			</div>
			<div class="cell">
				Exponential (homogeneous) background
			</div>
		</div>


		<div class="row" onclick="window.location='models/bg_strexp.html#bg_strexp'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">bg_strexp</span>			</div>
			<div class="cell">
				Stretched exponential (fractal) function
			</div>
		</div>


		<div class="row" onclick="window.location='models/bg_sumstrexp.html#bg_sumstrexp'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">bg_sumstrexp</span>			</div>
			<div class="cell">
				Sum of two stretched exponential (fractal) backgrounds
			</div>
		</div>


		<div class="row" onclick="window.location='models/bg_prodstrexp.html#bg_prodstrexp'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">bg_prodstrexp</span>			</div>
			<div class="cell">
				Product of two stretched exponential (fractal) backgrounds
			</div>
		</div>


		<div class="row" onclick="window.location='models/bg_poly1.html#bg_poly1'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">bg_poly1</span>			</div>
			<div class="cell">
				1st-order polynomial background
			</div>
		</div>


		<div class="row" onclick="window.location='models/bg_poly2.html#bg_poly2'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">bg_poly2</span>			</div>
			<div class="cell">
				2nd-order polynomial background
			</div>
		</div>


		<div class="row" onclick="window.location='models/bg_poly3.html#bg_poly3'">
			<div class="cell">
				<span class="std std-ref" style="font-family: monospace">bg_poly3</span>			</div>
			<div class="cell">
				3rd-order polynomial background
			</div>
		</div>

	</div>
	</div>
	</div>
	</div>
