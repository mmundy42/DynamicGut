Prepare a simulation
====================

Before running a simulation with DynamicGut, you must prepare two input files
and construct pair models for all combinations of two species in the community.

Create an initial diet file
---------------------------

A diet identifies the metabolites that are available for the members of the
community to consume. An initial diet is the starting point for the simulation.
By defining different diets you can simulate how the community responds to
different nutrient conditions.

A diet is expressed in terms of exchange reactions that allow metabolites to
move between the boundary and either the extracellular compartment in a single
species model or the shared compartment in a community model. The bound on the
exchange reaction controls how much of the metabolite is available to be consumed.
Exchange reactions must be identified with an ID that starts with ``EX_`` (which
is typical in many models).

A diet is defined in a diet file that lists all of the metabolite exchange
reactions and a bound for the reaction in JSON format. The bound must be a
positive value. For example::

    {
        "EX_ascb_L_LPAREN_e_RPAREN_": 0.745862875,
        "EX_gncore2_LPAREN_e_RPAREN_": 0.1,
        "EX_fol_LPAREN_e_RPAREN_": 0.1
        ...
    }

There must be an entry in the diet file for every exchange reaction that is
defined in the single species models. Any exchange reactions that are defined
in the diet file but are not in any single species model are ignored.

Create an initial density file
------------------------------

A population density is the size of an organism's population in the community.
An initial population density is the starting point for the simulation. By defining
different population densities you can simulate how the community responds to
different compositions of the community.

A population density is expressed as a value of something. The population density
for each organism in the community is defined in a density file that lists the
IDs of all of the single species models and the value in CSV format. For example::

    ID,DENSITY
    Btheta,1e-5
    Ecoli,1e-5
    Erectale,1e-5

Construct two species community models
--------------------------------------

When running a simulation, DynamicGut optimizes a two species community model
for every pair of organisms in the community under the current nutrient
conditions. As input, you provide the single species models for each member of
the community. Each single species model should produce growth using the medium
set in the input model. The ``prepare()`` function creates two species community
models from the single species models and stores the pair model files in a folder.
You can optionally have the ``prepare()`` function confirm that the single species
models produce growth. For example::

    import dynamicgut
    pair_file_names = dynamicgut.prepare()
