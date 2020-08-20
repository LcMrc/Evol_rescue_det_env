# Evol_rescue_det_env
Matlab implementation of the code described in "Adapt or perish: Evolutionary rescue in a gradually deteriorating environment", by Loïc Marrec and Anne-Florence Bitbol.

This code contains a Matlab implementation of the code described in "Adapt or perish: Evolutionary rescue in a gradually deteriorating environment", Loïc Marrec and Anne-Florence Bitbol.

Archived version: DOI: 10.5281/zenodo.3993272 

Briefly, we perform stochastic simulations of a microbial population evolving in a deteriorating environment where two types of microbes can exist, namely wild-type and mutant (either generalist or mutant). The wild-type fitness (or reproduction rate) decreases according to a Hill function, while the generalist mutant fitness remains constant over time and the specialist mutant fitness increases according to a Hill function.

In order to use the code, please run either "generate_data" in the subfolder "Generalist_mutant" or "Specialist_mutant" of the folders "Fixation_probability_pfix", "Rescue_probability_pr" and "Time_of_appearance_of_the_mutant_that_fixes_tauaf" under Matlab.

The source code is freely available under the GNU GPLv3 license.

If you find this code useful for your research, please cite the associated reference, "Adapt or perish: Evolutionary rescue in a gradually deteriorating environment", by Loïc Marrec and Anne-Florence Bitbol.
