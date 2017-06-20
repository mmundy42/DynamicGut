Run a simulation
================

A DynamicGut simulation predicts the population densities of a microbial community
over time given initial diet conditions. See :doc:`prepare` for details on creating
the initial density and diet files.

A simulation runs over a range of time points. At each time point, DynamicGut
performs the following operations.

Optimize single species models
------------------------------

Each single species model is optimized using the current diet conditions. Why do
we do this -- need to have a good explanation.

The results are stored in a file named "single-rates-NNNN.csv" that has this
format::

    ID,RATE

Calculate growth rates of each pair community model
---------------------------------------------------

Each pair community model is optimized three times: (1) with both pairs, (2) with
first species knocked out, and (3) with second species knocked out.

The results are stored in a file named "pair-rates-NNNN.csv" that has this format::

    GROWTH_OF,X,Y,Z

Calculate the effects of growth rates
-------------------------------------

The effect is the percent change in growth of one species in the presence of
another species.

The results are stored in a file named "effects-matrix-NNN.csv" that has this
format::

    PERCENT_CHANGE,X,Y,Z

Calculate new population densities
----------------------------------

Run the Leslie-Gower algorithm to calculate updated population densities.

The results are stored in a file named "density-NNNN.csv" that has this format::

    ID,DENSITY

Calculate new diet conditions
-----------------------------

Update the diet conditions by adjusting the bounds on the exchange reactions
based on the single species optimization and new population densities.

The results are stored in a file named "diet-NNNN.json" that has this format::

    {
        "EX_ascb_L_LPAREN_e_RPAREN_": 0.745862875,
        "EX_gncore2_LPAREN_e_RPAREN_": 0.1,
        "EX_fol_LPAREN_e_RPAREN_": 0.1
        ...
    }
