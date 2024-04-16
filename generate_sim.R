## Auxiliary function to simulate logistic data.
expit <- function (x) exp(x) / (1 + exp(x))

## Set up the simulation - set sample size and common cause.
set.seed(54545)
n <- 1e5
U <- rnorm(n, 0, 1)

##Simulate one-sample genetic data from individual-level data - number of SNPs, simulate MAF, ratios for how many affect I and P.
nSNPs <- 100
maf <- runif(nSNPs, 0.05, 0.5)
incidence.SNPs <- 1:90
progression.SNPs <- 91:100
G <- matrix(0, n, nSNPs)
for (j in 1:nSNPs) G[, j] <- rbinom(n, 2, maf[j])
incidence.betas <- rep(0, nSNPs)
incidence.betas[incidence.SNPs] <- rnorm(length(incidence.SNPs), 0, 0.2)
incidence.probs <- expit(-1 + as.vector(G %*% incidence.betas) + 1 * U)
I <- rbinom(n, 1, incidence.probs)
progression.betas <- rep(0, nSNPs)
progression.betas[progression.SNPs] <- rnorm(length(progression.SNPs), 0, 0.2)
progression.probs <- expit(-1 + as.vector(G %*% progression.betas) + 1 * U)
P <- rbinom(n, 1, progression.probs)

##This seems to work - it gives you a matrix with a seperate beta for 90 SNPs, but not sure if it is right as I cannot get it to work for progression.

incidence_GWAS <- matrix(0, length(incidence.SNPs))

for (j in incidence.SNPs) { model <- glm(I ~ G[,j], family = binomial) 
+ incidence_results[j] <- summary(model)$coefficients[2] }
