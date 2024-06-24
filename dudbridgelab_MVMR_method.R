#Siyang Cai and Frank Dudbridge method proposed at the MR conference 2024. Inspired by Apostolis' IVsel work.
This sim is based on individual-level data (exposure, disease and progression drawn from the same population).

#Potential MVMR methods to evaluate, looking at the effect of weak instrument bias on estimates at each stage:

#Egger/Grapple/Horse/Robust/Median/Lasso
  
#Traditional MVMR estimation - ideally in practice you would need individual-level data to estimate σi,j (pairwise covariances between G and X1,2,n..etc to be known across all SNPs). 
#How does estimation perform when 1. you have individual-level genetic data and σi,j is known 2. phenotypic correlations are estimated from other open-access data 3. you assume σi,j to be 0, Spiller et al. sims showed this 
#only provides correct estimation when global β on the outcome is null.
  
#Udp (incidence-progression confounders, i.e., traditional collider)
#Uep (exposure-progression confounders) - these are likely to be problematic as we are conducting a case-only study where the exposure that 'caused' individuals to get the disease is not 
#accounted for in the estimation, difficult to account for without longitudinal individual-level data...

#Under what circumstances do you get weak instrument bias, and is it a problem/worse than existing methods? 
#Is this method more useful for applied researchers than current GWAS correction methods (pleiotropy-robust IV methods seem to perform well under most circumstances).

n.ind <- 1e5
nSNPs <- 100
G <- matrix(0, n.ind, nSNPs)
Udy <- rnorm(n.ind, 0, 1)
Uxy <- rnorm(n.ind, 0, 1)

maf <- runif(nSNPs, 0.05, 0.5)
  for (j in 1:nSNPs) G[, j] <- rbinom(n.ind, 2, maf[j])
  X.betas <- rep(0, nSNPs)
  X.betas[X.SNPs] <- rnorm(length(X.SNPs), 0, 0.3)
  X.probs <- expit(-1 + as.vector(G %*% X.betas) + 1 * Uxy)
  X <- rbinom(n.ind, 1, X.probs)  D.betas <- rep(0, nSNPs)

  D.betas[D.SNPs] <- rnorm(length(D.SNPs), 0, 0.3)
  D.probs <- expit(-1 + as.vector(G %*% D.betas) + 1 * Udy)
  D <- rbinom(n.ind, 1, D.probs)

  Y.betas[Y.SNPs] <- rnorm(length(Y.SNPs), 0, 0.3)
  Y.probs <- expit(-1 + as.vector(G %*% Y.betas) + 1 * Udy)
  Y <- rbinom(n.ind, 1, Y.probs)
  Yobs <- Y; Yobs[Y == 0] <- NA

MVMR regression


