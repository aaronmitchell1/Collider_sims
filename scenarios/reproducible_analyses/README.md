# Scenario 1

Comparing bias in P vs Pobs. 

The pobs_vs_p_bias.R script shows how to reproduce the analysis and generate performance metrics.

loop_methods_P.RData contains results for methods run on P (on all individuals) for 1000 iterations.

loop_methods_P_obs.RData contains results for methods run on Pobs (on only individuals with the disease) for 1000 iterations.

loop_var_interaction.RData contains the beta*var(U) needed to calculate the true value from this DGM (this does not depend on P/Pobs).
