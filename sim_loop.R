expit <- function (x) exp(x) / (1 + exp(x))
set.seed(54545)
iter <- 1000
n <- 1e5
nSNPs <- 100
incidence.SNPs <- 1:90
progression.SNPs <- 91:100
library(MendelianRandomization)
library(mr.raps)
library(SlopeHunter)

loop_methods <- vector('list', iter)

for (i in 1:iter) {

G <- matrix(0, n, nSNPs)

##Simulate common cause (collider).
U <- rnorm(n, 0, 1)

##Simulate one-sample genetic data from individual-level data - number of SNPs, generate MAF, ratios for how many affect I and P.
maf <- runif(nSNPs, 0.05, 0.5)
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
  
for (j in 1:nSNPs) {
  
  incidence_model <- glm(I ~ G[,j], family = binomial)
  incidence_GWAS[i, j, 1] <- summary(incidence_model)$coefficients["G[, j]", "Estimate"]
  incidence_GWAS[i, j, 2] <- summary(incidence_model)$coefficients["G[, j]", "Std. Error"]
  incidence_GWAS[i, j, 3] <- summary(incidence_model)$coefficients["G[, j]", "Pr(>|z|)"]
  incidence_GWAS[i, j, 4] <- coef(incidence_model)[1]

}

##Simulate GWAS of progression.
  
for (j in 1:nSNPs) {
  progression_model <- glm(P ~ G[,j], family = binomial)
  progression_GWAS[i, j, 1] <- summary(progression_model)$coefficients["G[, j]", "Estimate"]
  progression_GWAS[i, j, 2] <- summary(progression_model)$coefficients["G[, j]", "Std. Error"]
  progression_GWAS[i, j, 3] <- summary(progression_model)$coefficients["G[, j]", "Pr(>|z|)"]
  progression_GWAS[i, j, 4] <- coef(progression_model)[1]

}

collider_bias_results <- data.frame()

##Run the methods

##Dudbridge method (updated CWLS from Cai et al. paper), based on April Hartley's code.

ivw <- mr_ivw(incidence_GWAS[, 1], progression_GWAS[, 1], 
                 incidence_GWAS[, 2], progression_GWAS[, 2])
dudbridgeweights <- 1/progression_GWAS[, 2]^2
weighting <- (sum(dudbridgeweights*incidence_GWAS[, 1]^2))/
            ((sum(dudbridgeweights*incidence_GWAS[, 1]^2))-(sum(dudbridgeweights*incidence_GWAS[, 2]^2)))
cf.db <- ivw$b*weighting
cf.se.db <- ivw$se*weighting

collider_bias_results <- rbind(collider_bias_results,
                               data.frame(Method = 'Dudbridge',
                               Correction_Beta = cf.db,
                               Correction_SE = cf.se.db
                                 ))
  
##Weighted median

Weighted_Median_Res <- mr_median((mr_input(bx = incidence_GWAS[, 1], by = progression_GWAS[, 1], 
                                           bxse = incidence_GWAS[, 2], byse = progression_GWAS[, 2])))

Weighted_Median_Beta <- Weighted_Median_Res$Estimate
Weighted_Median_SE <- Weighted_Median_Res$StdError
  
collider_bias_results <- rbind(collider_bias_results,
                                 data.frame(Method = 'Weighted_median',
                                 Correction_Beta = Weighted_Median_Beta,
                                 Correction_SE = Weighted_Median_SE
                                 ))
  
##MR-RAPS

MR_RAPS_Res <- mr.raps(incidence_GWAS[, 1], incidence_GWAS[, 2], 
                       progression_GWAS[, 1], progression_GWAS[, 2])

MR_RAPS_Beta <- MR_RAPS_Res$beta.hat
MR_RAPS_SE <- MR_RAPS_Res$beta.se
  
collider_bias_results <- rbind(collider_bias_results,
                                 data.frame(Method = 'MR_RAPS',
                                 Correction_Beta = MR_RAPS_Beta,
                                 Correction_SE = MR_RAPS_SE
                                 ))

##MR-Horse
##Reformat data and run the MR-Horse method.

MR_Horse_Data <- data.frame(betaX = incidence_GWAS[, 1], betaY = progression_GWAS[, 1], 
                            betaXse = incidence_GWAS[, 2], betaYse = progression_GWAS[, 2])
  
MR_Horse_Res <- mr_horse(MR_Horse_Data)
  
MR_Horse_Beta <- MR_Horse_Res$MR_Estimate$Estimate
MR_Horse_SE <- MR_Horse_Res$MR_Estimate$SD
  
collider_bias_results <- rbind(collider_bias_results,
                               data.frame(Method = 'MR_Horse',
                               Correction_Beta = MR_Horse_Beta,
                               Correction_SE = MR_Horse_SE
                                 )) 
  
##Slope-Hunter
##Reformat data and run the SlopeHunter method.

SlopeHunter_Data <- data.frame(xbeta_col = incidence_GWAS[, 1], xse_col = incidence_GWAS[, 2], 
                               ybeta_col = progression_GWAS[, 1], yse_col = progression_GWAS[, 2])
  
SlopeHunter_Res <- hunt(dat = SlopeHunter_Data, xbeta_col = 'xbeta_col', 
                        xse_col = 'xse_col', ybeta_col = 'ybeta_col', yse_col = 'yse_col')
  
Slopehunter_Beta <- SlopeHunter_Res$b
Slopehunter_SE <- SlopeHunter_Res$bse
  
collider_bias_results <- rbind(collider_bias_results,
                               data.frame(Method = 'Slopehunter',
                               Correction_Beta = Slopehunter_Beta,
                               Correction_SE = Slopehunter_SE
                                 ))

##Summarise methods results for each iteration
  
loop_methods[[i]] <- list(collider_bias_results = collider_bias_results)
  
}
