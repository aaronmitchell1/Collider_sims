##Function to generate sims - first attempt at writing a function so very much under development!

#nSNPs = total number of SNPs in analysis
#nind = number of individuals in 'GWAS'
#incidenceSNPs = how many SNPs associated with incidence, this needs to be inputted as a ratio i.e., 1:90
#progressionSNPs = how many SNPs associated with progression, this also needs to be inputted as a ratio
#genetic correlation between I and P SNPs can then be induced by setting the ratio of SNPs between incidence and progression
#minallelef = minor allele frequency in 'GWAS' population , default is to generate MAF between 0.05%-5% from uniform distribution
#collider = induce collider bias through common cause of I and P, either yes or no so you can get the 'right' answer
#seed = random seed for data generation to improve reproducibility

generate_sim <- function(nSNPs, nind, incidenceSNPs, progressionSNPs, minallelef = c(0.05, 0.5), collider, seed=54545)

## Auxiliary function to simulate logistic data.
expit <- function (x) exp(x) / (1 + exp(x))

## Set up the simulation - set sample size and common cause.
                         
set.seed(seed)

if (collider == "yes") {
    U <- rnorm(n, 0, 1)
  } 

else if (collider == "no") {
    U <- NA
  } 

##Simulate one-sample genetic data from individual-level data - number of SNPs, simulate MAF, ratios for how many affect I and P.
maf <- runif(nSNPs, min =  minallelef[1], max = minallelef[2])
G <- matrix(0, nind, nSNPs)
for (j in 1:nSNPs) G[, j] <- rbinom(n, 2, maf[j])
incidence.betas <- rep(0, nSNPs)
incidence.betas[incidence.SNPs] <- rnorm(length(incidence.SNPs), 0, 0.2)
incidence.probs <- expit(-1 + as.vector(G %*% incidence.betas) + 1 * U)
I <- rbinom(n, 1, incidence.probs)
progression.betas <- rep(0, nSNPs)
progression.betas[progression.SNPs] <- rnorm(length(progression.SNPs), 0, 0.2)
progression.probs <- expit(-1 + as.vector(G %*% progression.betas) + 1 * U)
P <- rbinom(n, 1, progression.probs)

##Simulate GWAS of incidence.
#To include effects of progression SNPs on incidence and vice versa to conduct subsequent analyses,
#don't need to separate this into incidence and progression SNPs, can just be one loop for nSNPs.

incidence_GWAS <- matrix(0, nSNPs, 3); colnames(incidence_GWAS) <- c("Estimate", "StdErr", "Pval")

for (j in 1:nSNPs) {
  incidence_model <- glm(I ~ G[,j], family = binomial)
  incidence_GWAS[j, 1] <- summary(incidence_model)$coefficients["G[, j]", "Estimate"]
  incidence_GWAS[j, 2] <- summary(incidence_model)$coefficients["G[, j]", "Std. Error"]
  incidence_GWAS[j, 3] <- summary(incidence_model)$coefficients["G[, j]", "Pr(>|z|)"]
}

##Simulate GWAS of progression.

progression_GWAS <- matrix(0, nSNPs, 3); colnames(progression_GWAS) <- c("Estimate", "StdErr", "Pval")

for (j in 1:nSNPs) {
  progression_model <- glm(P ~ G[,j], family = binomial)
  progression_GWAS[j, 1] <- summary(progression_model)$coefficients["G[, j]", "Estimate"]
  progression_GWAS[j, 2] <- summary(progression_model)$coefficients["G[, j]", "Std. Error"]
  progression_GWAS[j, 3] <- summary(progression_model)$coefficients["G[, j]", "Pr(>|z|)"]

}

##For working with HPC, need to take into account potentially performing analyses on separate cores in parallel. Section 4.1.1 of Tim Morris' et al. ADEMP paper mentions using rstream to do this.
##Could potentially add a systemtime wrapper for each method as time/performance trade-off could be part of the paper
#Storing systemtime is simple to implement like so:
#t1 <- system.time (x)
