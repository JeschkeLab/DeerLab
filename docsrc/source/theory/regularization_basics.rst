.. highlight:: matlab
.. _regularization_basics:

***********************
Regularization - Basics
***********************

To solve the ill-posed problem

.. math:: \mathbf{K}\mathbf{P} = \mathbf{S}

the typical approach is to replace the original problem by a least-squares minimization problem. However, the solution to this problem is still highly sensitive to perturbation. In order to stabilize the solution, a penalty term :math:`\Gamma[\mathbf{L}\mathbf{P}]` is added. Therefore, the general form of a regularization method is given by

.. math:: \arg\min_{\mathbf{P}\geq 0}\left\{ \frac{1}{2}\Vert \mathbf{K}\mathbf{P} - \mathbf{S} \Vert_2^2 + \alpha^2 \Gamma[\mathbf{L}\mathbf{P}]\right\},

where :math:`\alpha` is the regularization parameter which controls the balance between the regularization penalty and the least-squares agreement of the data. Depending on the choice of :math:`\Gamma[\mathbf{L}\mathbf{P}]` the regularization penalty will react different to changes in :math:`\mathbf{P}`.

============ =========================================================
   Penalty               :math:`\Gamma[\mathbf{L}\mathbf{P}]`
============ =========================================================
Tikhonov     :math:`\Vert \mathbf{L}\mathbf{P} \Vert_2^2`
TV           :math:`\sum \sqrt{(\mathbf{L}\mathbf{P})^2 + \beta^2 }`
pseudo-Huber :math:`\sum \sqrt{(\mathbf{L}\mathbf{P}/\eta)^2 + 1 }-1`
============ =========================================================

Global fitting
""""""""""""""""""
If several signals :math:`S_1,S_2,...S_N` corresponding to the same distance distribution are available, they may be fitted to that single distance distribution simultaneously. This procedure has been often referred to as global fitting, which aims to solve the following regularization problem

.. math:: \arg\min_{\mathbf{P}\geq 0}\left\{ \sum_i^N w_i\frac{1}{2}\Vert \mathbf{K}_i\mathbf{P} - \mathbf{S}_i \Vert_2^2 + \alpha^2 \Gamma[\mathbf{L}\mathbf{P}]\right\},

where :math:`w_i` are to so-called global-fit weights. They distribute the influence each signal has on the quality of the least-squares fitting. In DeerAnalysis2 these weights are computed according to the contribution of each signal to the ill-posedness of the inverse problem.

.. math:: w_i = \frac{\sum_j^N M_j\sigma_j}{N_i\sigma_i}, \sum_i^N w_i = 1

where :math:`M_i` is the length and :math:`\sigma_i` the noise level of the signal :math:`S_i`.


References
=========================================