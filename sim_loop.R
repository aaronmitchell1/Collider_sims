expit <- function (x) exp(x) / (1 + exp(x))
set.seed(54545)
iter <- 1000
n <- 1e5
nSNPs <- 100
incidence.SNPs <- 1:90
progression.SNPs <- 91:100

loop_results <- list()
loop_methods <- list()

collider_bias_results <- data.frame(
  Method = character(),
  Correction_Beta = numeric(),
  Correction_SE = numeric()
)

collider_bias_type <- list(
  Slopehunter = "Slopehunter"
  Dudbridge = "Dudbridge",
  Weighted_median = "Weighted_median",
  MR_Horse = "MR_Horse",
  MR_RAPS = "MR_RAPS"
)

G <- matrix(0, n, nSNPs)

for (i in 1:iter) {

##Simulate common cause (collider).
U <- rnorm(n, 0, 1)

##Simulate one-sample genetic data from individual-level data - number of SNPs, simulate MAF, ratios for how many affect I and P.
maf <- runif(nSNPs, 0.05, 0.5)
G[, j] <- rbinom(n, 2, maf[j])
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
  incidence_GWAS[j, 1] <- summary(incidence_model)$coefficients["G[, j]", "Estimate"]
  incidence_GWAS[j, 2] <- summary(incidence_model)$coefficients["G[, j]", "Std. Error"]
  incidence_GWAS[j, 3] <- summary(incidence_model)$coefficients["G[, j]", "Pr(>|z|)"]
}

##Simulate GWAS of progression.
  
for (j in 1:nSNPs) {
  progression_model <- glm(P ~ G[,j], family = binomial)
  progression_GWAS[j, 1] <- summary(progression_model)$coefficients["G[, j]", "Estimate"]
  progression_GWAS[j, 2] <- summary(progression_model)$coefficients["G[, j]", "Std. Error"]
  progression_GWAS[j, 3] <- summary(progression_model)$coefficients["G[, j]", "Pr(>|z|)"]
}
  
##Store results 
loop_results[[i]] <- list(incidence_GWAS_b = incidence_GWAS[j, 1], 
                          incidence_GWAS_se = incidence_GWAS[j, 2],
                          incidence_GWAS_p = incidence_GWAS[j, 3],
                          progression_GWAS_b = progression_GWAS[j, 1],
                          progression_GWAS_se = progression_GWAS[j, 2],
                          progression_GWAS_p = progression_GWAS[j, 3]
  )

##Run the methods

##Dudbridge method, based on April Hartley's code.

ivw <- mr_ivw(incidence_GWAS[, 1], progression_GWAS[, 1], 
                 incidence_GWAS[, 2], progression_GWAS[, 2])
dudbridgeweights <- 1/progression_GWAS[, 2]^2
weighting <- (sum(dudbridgeweights*incidence_GWAS[, 1]^2))/
            ((sum(dudbridgeweights*incidence_GWAS[, 1]^2))-(sum(dudbridgeweights*incidence_GWAS[, 2]^2)))
cf.db[i,] <- ivw$b*weighting
cf.se.db[i,] <- ivw$se*weighting

##Weighted median method.

library(MendelianRandomization)
Weighted_Median_Res[i,] <- mr_median((mr_input(bx = incidence_GWAS$Estimate, bxse = incidence_GWAS$StdErr,
                                       by = progression_GWAS$Estimate, byse = progression_GWAS$StdErr)))

##MR-RAPS

library(mr.raps)
MR_RAPS_Res[i,] <- mr.raps(incidence_GWAS$Estimate, incidence_GWAS$StdErr, 
                       progression_GWAS$Estimate, progression_GWAS$StdErr)

##MR-Horse
##Reformat data and run the MR-Horse method.

MR_Horse_Data <- data.frame(incidence_GWAS$Estimate, progression_GWAS$Estimate, 
                            incidence_GWAS$StdErr, progression_GWAS$StdErr)
names(MR_Horse_Data)[1] <- "betaX"
names(MR_Horse_Data)[2] <- "betaY"
names(MR_Horse_Data)[3] <- "betaXse"
names(MR_Horse_Data)[4] <- "betaYse"
MR_Horse_Res[i,] <- mr_horse(MR_Horse_Data)

##Slope-Hunter
##Reformat data and run the SlopeHunter method.

SlopeHunter_Data <- data.frame(incidence_GWAS$Estimate, incidence_GWAS$StdErr, 
                               progression_GWAS$Estimate, progression_GWAS$StdEr)
names(SlopeHunter_Data)[1] <- "xbeta_col"
names(SlopeHunter_Data)[2] <- "ybeta_col"
names(SlopeHunter_Data)[3] <- "xse_col"
names(SlopeHunter_Data)[4] <- "yse_col"

SlopeHunter_Res[i,] <- hunt(dat = SlopeHunter_Data, xbeta_col = "xbeta_col", 
                        xse_col = "xse_col", ybeta_col = "ybeta_col", yse_col = "yse_col")

##Summarise results
##Add Slopehunter results

collider_bias_results[i,] <- dplyr::add_row(collider_bias_results,
                                                Method = collider_bias_type$Slopehunter,
                                                Correction_Beta = SlopeHunter_Res$b,
                                                Correction_SE = SlopeHunter_Res$bse
        )

##Add Dudbridge results

collider_bias_results[i,] <- dplyr::add_row(collider_bias_results,
                                                Method = collider_bias_type$Dudbridge,
                                                Correction_Beta = cf.db,
                                                Correction_SE = cf.se.db
        )

##Add Weighted Median results

collider_bias_results[i,] <- dplyr::add_row(collider_bias_results,
                                        Method = collider_bias_type$Weighted_median,
                                        Correction_Beta = Weighted_Median_Res$b,
                                        Correction_SE = Weighted_Median_Res$se
)

##Add MR-RAPS results

collider_bias_results[i,] <- dplyr::add_row(collider_bias_results,
                                        Method = collider_bias_type$MR_RAPS,
                                        Correction_Beta = MR_RAPS_Res$beta.hat,
                                        Correction_SE = MR_RAPS_Res$beta.se
)

##Add MR-Horse results

collider_bias_results[i,] <- dplyr::add_row(collider_bias_results,
                                        Method = collider_bias_type$MR_Horse,
                                        Correction_Beta = MR_Horse_Res$MR_Estimate$Estimate,
                                        Correction_SE = MR_Horse_Res$MR_Estimate$SD
)

##Store results
loop_methods[[i]] <- list(collider_bias_results = collider_bias_results)
}
