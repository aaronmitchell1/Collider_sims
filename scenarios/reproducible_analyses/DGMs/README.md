This folder contains all DGMs that were used in our simulation study. Methods were run on Pobs in all scenarios but an additional GWAS was also simulated on P in a hypothetical complete case analysis. All scenarios used a binary incidence trait and a continuous progression trait.

Scripts in the format '401050' refers to n incidence SNPs, progression SNPs and SNPs affecting both in that order. 

All of the following scenarios have 90 SNPs affecting incidence and 10 affecting progression with no genetic overlap, a scenario highlighted in the Ganna 2023 paper.

Scripts in the format '0.2' refers to altering the P*U effect when generating progression probabilities. 

'incb0.1' refers to reducing the SD of incidence betas.

'inc_int' adds interactions for 5 SNPs when generating incidence betas with beta 0.1*U which are then fed into a GWAS which is naïve to the interactions.

'prog_int' adds interactions for 5 SNPs when generating progression betas with beta 0.1*U which are then fed into a GWAS which is naïve to the interactions.

'post-sim-analysis' shows the analysis pipeline for each method in rsimsum and subsequent P/Pobs calculations.
