# Scenario 1

Comparing bias in P vs Pobs. 

The pobs_vs_p_bias.R script shows how to reproduce the analysis and generate performance metrics with the baseline scenario used in our study, incidence.SNPs <- 1:90, progression.SNPs <- 91:100, binary I and continuous P (to avoid issues with non-collapsibility of ORs and to make calculating the true value of θ easier due to DAG path rules), no epistatic interactions.

loop_methods_P.RData contains results for methods run on P (on all individuals) for 1000 iterations.

loop_methods_Pobs.RData contains results for methods run on Pobs (on only individuals with the disease) for 1000 iterations.

loop_var_interaction.RData contains the β I*var(U) variable needed to calculate the true value of θ from this DGM (this does not depend on P/Pobs).
