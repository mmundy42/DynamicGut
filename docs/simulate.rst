Run a simulation
================

A DynamicGut simulation predicts the population sizes of a microbial community
over time given initial diet conditions. See :doc:`prepare` for details on creating
the input files.

A simulation runs over a range of time points. The default unit for a time point is
one hour. At each time point, DynamicGut performs the operations described below.

Note, in the file names below, "NNNN" is a number identifying the time point. For
example, "0002" is for the second time point.

Calculate growth rates of every pair in community
-------------------------------------------------

Each two species community model is optimized three times: (1) with both species,
(2) with first species knocked out, and (3) with second species knocked out.

The results are stored in a file named "pair-rates-NNNN.csv" that has this format::

    A_ID,B_ID,A_TOGETHER,A_ALONE,A_EFFECT,B_TOGETHER,B_ALONE,B_EFFECT
    Btheta,Ecoli,10.1946336186,0.00199613435042,5107.18810908,0.0,4.0422673957,0.0
    Btheta,Erectale,0.00199613435042,0.00199613435042,1.0,3.97515044271,3.97592893938,0.999804197539
    Ecoli,Erectale,4.0422673957,4.0422673957,1.0,0.0,3.97592893938,0.0

Each row shows the growth rates for a pair of species in the community on the
current diet conditions where:

* A_ID is the model ID of species A.
* B_ID is the model ID of species B.
* A_TOGETHER is the growth rate of species A together with species B.
* A_ALONE is the growth rate of species A with species B knocked out.
* A_EFFECT is the magnitude of the effect in the growth rate of species A in the
  presence of species B calculated as ``A_TOGETHER / A_ALONE``.
* B_TOGETHER is the growth rate of species B together with species A.
* B_ALONE is the growth rate of species B with species A knocked out.
* B_EFFECT is the magnitude of the effect in the growth rate of species B in the
  presence of species A calculated as ``B_TOGETHER / B_ALONE``.

Calculate an effects matrix from growth rates
---------------------------------------------

Each cell in an effects matrix is the effect on the growth of one species in
the presence of another species. A row gives the magnitude of the effect of the
growth of one species in the presence of all of the other species in the community.
The diagonal in the matrix is always 1.

The results are stored in a file named "effects-matrix-NNNN.csv" that has this
format::

    EFFECT,Btheta,Ecoli,Erectale
    Btheta,1.0,5107.18810916,1.0
    Ecoli,0.0,1.0,1.0
    Erectale,0.999804197539,0.0,1.0

For example, the effect of ``Erectale`` on the growth rate of ``Btheta`` is ``1.0``.
Note the column headers are the model IDs and are different based on the members
in the community.

Update population sizes
-----------------------

The Leslie equation for calculating the population sizes of two competing species
uses the following terms:

* *lambda* is a logistic parameter for the species when it is living alone and is
  defined to be ``1 + B(t) - D(t)`` where *B(t)* is the birth rate at time point t and
  *D(t)* is the death rate at time point t (Equation 2-3). When a model is optimized
  by DynamicGut, the output is the growth rate at time point t. DynamicGut makes
  the assumption that ``G(t) = B(t) - D(t)`` where *G(t)* is the growth rate at time
  point t. So DynamicGut calculates *lambda* as ``1 + G(t)`` .
* *K* is the upper asymptote of the population size or the carrying capacity (i.e.
  the maximum size of the population that the environment has the capacity to
  support). DynamicGut uses only one value of *K* for all species in the microbial
  community.
* *alpha* is a logistic parameter for the species and is defined to be
  ``lambda - 1 / K`` (Equation 3-3).
* *gamma* is the magnitude of the effect which each species has on the rate of
  increase of the other species.
* *N(t)* is the population size at time point t.
* *N(t+1)* is the population size at time point t+1 (or the next time point).

Using Equation 3-4, the population sizes of the two species at the next time point
is defined to be:

N\ :sub:`1`\ (t+1) = lambda\ :sub:`1` * N\ :sub:`1`\ (t) / 1 + alpha\ :sub:`1` * N\ :sub:`1`\ (t) + gamma\ :sub:`1` * N\ :sub:`2`\ (t)

N\ :sub:`2`\ (t+1) = lambda\ :sub:`2` * N\ :sub:`2`\ (t) / 1 + alpha\ :sub:`2` * N\ :sub:`2`\ (t) + gamma\ :sub:`2` * N\ :sub:`1`\ (t)

The results are stored in a file named "density-NNNN.csv" that has this format::

    ID,DENSITY
    Btheta,1.22723737851e-05
    Ecoli,9.12613613888e-05
    Erectale,8.92668652468e-05

where:

* ID is the ID of the single species model
* DENSITY is the population size of the species at the end of the time point

Update diet conditions
----------------------

At the end of the time step, the diet conditions are adjusted to reflect the
consumption and production of metabolites by the members of the community. The
new diet conditions are used for the calculations in the next time point.

Each single species model is optimized on the current diet conditions for
*a good reason*. For each species, the bound on each exchange reaction in the
current diet is updated by the flux of the exchange reaction adjusted by the
species population size. Remember that there is variability in the results
of a flux balance analysis so the values are not exactly the same given the
same initial diet and population sizes.

The results of the single species optimizations are stored in a file named
"single-rates-NNNN.csv" that has this format::

    ID,STATUS,GROWTH_RATE
    Bacteroides_thetaiotaomicron_VPI_5482,optimal,0.00199613435042
    Escherichia_coli_str_K_12_substr_MG1655,optimal,4.0422673957
    Eubacterium_rectale_ATCC_33656,optimal,3.97592893938

where:

* ID is the ID of the single species model.
* STATUS is the status of the optimization as returned by the solver.
* GROWTH_RATE is the growth rate of the species on the current diet conditions.

The new diet conditions are stored in a file named "diet-NNNN.json" that has
this format::

    {
        "EX_ascb_L_LPAREN_e_RPAREN_": 0.745862875,
        "EX_gncore2_LPAREN_e_RPAREN_": 0.1,
        "EX_fol_LPAREN_e_RPAREN_": 0.1
        ...
    }

Output files
------------

At the end of the simulation, the following files are stored in the folder specified
with the ``data_folder`` parameter.

* initial-diet.json contains the initial diet conditions used as input to the first
  time point in the simulation
* final-diet.json contains the final diet conditions that are output from the last
  time point in the simulation
* timepoint-NNNN is a folder with the files created at each time point in the
  simulation that are described above