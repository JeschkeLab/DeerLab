Models
......................................................................


.. _ddmodels:

Distribution models
=========================================

In DeerLab, distance distributions can be represented either as parameter-free models or as parametric models.

Parameter-free distance distributions
------------------------------------------------

.. rst-class:: coderef

A parameter-free distance distribution is a vector of population densities ``P`` defined over a vector of distances ``r``. Outside the range of ``r``, the distance distribution is considered zero, i.e. it is truncated to the range of ``r``. Such parameter-free distance distributions are returned by the fit functions :ref:`fitregmodel` and :ref:`obir`. A ``P`` can be converted to a time-domain signal using :ref:`dipolarsignal` or :ref:`dipolarkernel`.

Parameter-free distance distributions are preferred over parametric distance distributions, since they make fewer assumptions about the distribution and are more flexible. They introduce less bias.

Parametric distance distributions
------------------------------------------------

.. rst-class:: coderef

DeerLab provide a series of parametric distance distributions models. These can be used with several functions, including :ref:`fitparamodel` and :ref:`selectmodel`. All parametric distance distribution models are functions start with the prefix ``dd_``. They take a distance vector ``r`` and a parameter vector ``param`` as inputs and return the distance distribution as a vector ``P``. Here is an example:

.. code-block:: matlab

        P = dd_gauss(r,param)


Outside the range of ``r``, the distance distribution is set to zero, i.e. it is truncated to the range of ``r`` and normalized such that its integral over the range of ``r`` equals one.

DeerLab's parametric distance distribution models can be classified into several groups. There are models that are based on basis functions and linear combinations thereof. Gaussian and 3D-Rice functions are the most important ones. To use a single Gaussian distribution, use :ref:`dd_gauss`. To combine 2, 3, etc. Gaussians distributions, use :ref:`dd_gauss2`, :ref:`dd_gauss3`, etc. Similarly, for 3D-Rice distributions, use :ref:`dd_rice`, :ref:`dd_rice2`, etc.

A second group of distance distribution models represent three-dimensional disordered segmented objects such as proteins and other polymers: :ref:`dd_wormchain`, :ref:`dd_wormgauss`, :ref:`dd_randcoil`.

A third group provides models for distributions of labels in simple confined spaces such as spheres and spherical shells.

Finally, DeerLab provides a series of distance distribution models that represent simple geometric shapes. These models are worthless for practical purposes, but are useful as toy distributions for methods development: :ref:`dd_circle`, :ref:`dd_triangle`, and :ref:`dd_uniform`.


.. _bgmodels:

Background models
======================

DeerLab includes a collection of parametric models that can be used to model the background signal, i.e. the signal due to randomly distributed spins in the sample that are not part of the spin-labeled protein or object.

All background model functions start with the prefix ``bg_`` and have the same calling syntax. The inputs a time axis vector ``t`` (in microseconds), a parameter vector ``param``, and a modulation amplitude ``lambda`` (between 0 and 1) . The length of ``param``, and the meaning of the elements, depends on the particular model. If ``lambda`` is not provided, it is set to one. The output is a background decay vector ``B``, defined over ``t``.

In the following example, the parameter is a single number and represents the spin concentration:

.. code-block:: matlab

        t = linspace(-1,5,301); % us
        c = 100; % uM
        lambda = 0.5;
        B = bg_hom3d(t,c,lambda);

All models are symmetric with respect to ``t=0``, meaning that they depend only on the magnitude of ``t``, ``|t|``.

Background models fall into two categories, physical and phenomenological.

Physical models describe particular distributions of spin labels in space. These models depend on physical parameters such as spin concentration, exclusion distances, and dimensionality. The most common one is :ref:`bg_hom3d`, which describes the signal due to a homogeneous three-dimensional distribution of spins of a given concentration. A homogeneous distribution in a fractal dimensions is available with :ref:`bg_homfractal`, and excluded-volume effects can be modelled using :ref:`bg_hom3dex`.


Phenomenological models represent various mathematical functions that are intended to mimick the background decay, without reference to a particular spatial distribution of spins. The parameters of these models do no have direct physical meaning. In general, it is preferable to use the physical instead of phenomenological models.

For a complete list of all background models, refer to the :ref:`background model reference section <modelsref_bg>`.


.. _exmodels:

Experiment models
=========================================

DeerLab supports a wide range of dipolar EPR experiments. For each type of supported dipolar EPR experiment, there is a dedicated experiment model function starting with ``ex_``. These functions take as inputs the time axis ``t`` and an array of parameters ``param``. As output, they return an array containing information about the dipolar pathways of the experiment model.

For example, the model function representing the typical model for a 4-pulse DEER signal is ``ex_4pdeer``:

.. code-block:: matlab

        t = linspace(0,3,151);
        lambda = 0.3;
        pathways = ex_4pdeer(t,lambda)

The returned output is

.. code-block:: matlab

    pathways =
              0.7          NaN
              0.3            0

Each row of this array holds information about one pathway. The first column is modulation amplitude, and the second column is the refocusing point. In the above example, the first row shows a pathway with amplitude 0.7 and no refocusing time, indicating that it represents the unmodulated contribution. The pathway of the second row shows amplitude of 0.3 and refocusing time 0, i.e. this is the primary dipolar pathway.

Experiment model functions are needed as inputs to ``fitsignal`` and ``dipolarsignal``.

Experiments differ in the number of dipolar modulation components and their refocusing times.

For a complete list of all experiment models, refer to the :ref:`experiment model reference section <modelsref_ex>`.
