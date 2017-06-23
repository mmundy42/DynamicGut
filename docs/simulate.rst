Run a simulation
================

A DynamicGut simulation predicts the population densities of a microbial community
over time given initial diet conditions. See :doc:`prepare` for details on creating
the initial density and diet files.

A simulation runs over a range of time points. In the file names below, "NNNN" is a
number identifying the time point. For example, "0002" is for the second time point.
At each time point, DynamicGut performs the following operations.

Calculate growth rates of every pair in community
-------------------------------------------------

Each two species community model is optimized three times: (1) with both species, (2) with
first species knocked out, and (3) with second species knocked out.

The results are stored in a file named "pair-rates-NNNN.csv" that has this format::

    A_ID,B_ID,A_TOGETHER,A_ALONE,A_CHANGE,B_TOGETHER,B_ALONE,B_CHANGE
    Btheta,Ecoli,10.1946336186,0.00199613435042,5107.18810908,0.0,4.0422673957,0.0
    Btheta,Erectale,0.00199613435042,0.00199613435042,1.0,3.97515044271,3.97592893938,0.999804197539
    Ecoli,Erectale,4.0422673957,4.0422673957,1.0,0.0,3.97592893938,0.0

Calculate an effects matrix from growth rates
---------------------------------------------

Each cell in an effects matrix is the effect on the growth of one species in
the presence of another species. A row gives the percent change in growth of
one species in the presence of all of the other species in the community. The
diagonal in the matrix is always 1.

The results are stored in a file named "effects-matrix-NNN.csv" that has this
format::

    PERCENT_CHANGE,Btheta,Ecoli,Erectale
    Btheta,1.0,5107.18810916,1.0
    Ecoli,0.0,1.0,1.0
    Erectale,0.999804197539,0.0,1.0

Note the column headers are the model IDs and are different based on the members
in the community.

Calculate new population densities
----------------------------------

Run the Leslie-Gower algorithm to calculate updated population densities.

The results are stored in a file named "density-NNNN.csv" that has this format::

    ID,DENSITY

Calculate new diet conditions
-----------------------------

At the end of the time step, the diet conditions are adjusted to reflect the
consumption and production of metabolites by the members of the community. The
new diet conditions are used for the calculations in the next time point.

Each single species model is optimized on the current diet conditions for
a good reason. For each species, the bound on each exchange reaction in the
current diet is updated by the flux of the exchange reaction adjusted by the
species population density. Remember that there is variability in the results
of a flux balance analysis so the values are not exactly the same given the
same initial diet and population densities.

The results of the single species optimizations are stored in a file named
"single-rates-NNNN.csv" that has this format::

    ID,STATUS,GROWTH_RATE
    Bacteroides_thetaiotaomicron_VPI_5482,optimal,0.00199613435042
    Escherichia_coli_str_K_12_substr_MG1655,optimal,4.0422673957
    Eubacterium_rectale_ATCC_33656,optimal,3.97592893938

The new diet conditions are stored in a file named "diet-NNNN.json" that has
this format::

    {
        "EX_ascb_L_LPAREN_e_RPAREN_": 0.745862875,
        "EX_gncore2_LPAREN_e_RPAREN_": 0.1,
        "EX_fol_LPAREN_e_RPAREN_": 0.1
        ...
    }

