Parametric Models
======================

DeerAnalysis includes the following collection of parametric models. The model names are categorized depending on whether the model is defined in time-domain (model name starts with the prefix ``td_``) or in distance model (model name starts with the prefix ``rd_``). 

.. toctree::
    :maxdepth: 1
    :caption: Distance-domain

    ./models/rd_onegaussian
    ./models/rd_twogaussian
    ./models/rd_threegaussian
    ./models/rd_fourgaussian
    ./models/rd_fivegaussian
    ./models/rd_onerice
    ./models/rd_tworice
    ./models/rd_threerice
    ./models/rd_wormchain
    ./models/rd_randcoil

.. toctree::
    :maxdepth: 1
    :caption: Time-domain

    ./models/td_exp
    ./models/td_strexp
    ./models/td_sumstrexp
    ./models/td_prodstrexp
    ./models/td_poly1
    ./models/td_poly2
    ./models/td_poly3
