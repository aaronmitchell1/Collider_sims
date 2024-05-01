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
library(R2jags)

mr_horse = function(D, no_ini = 3, variable.names = "theta", n.iter = 10000, n.burnin = 10000){
  if("theta" %in% variable.names){
    variable.names = variable.names
  } else{
    variable.names = c("theta", variable.names)
  }
  jags_fit = jags(data = list(by = D$betaY, bx = D$betaX, sy = D$betaYse, sx = D$betaXse, N = length(D$betaY)),
                  parameters.to.save = variable.names,
                  n.chains = no_ini,
                  n.iter = n.burnin + n.iter,
                  n.burnin = n.burnin,
                  model.file = mr_horse_model)
  mr.coda = as.mcmc(jags_fit)
  mr_estimate = data.frame("Estimate" = round(unname(summary(mr.coda[, "theta"])$statistics[1]), 3),
                           "SD" = round(unname(summary(mr.coda[, "theta"])$statistics[2]), 3),
                           "2.5% quantile" = round(unname(summary(mr.coda[, "theta"])$quantiles[1]), 3),
                           "97.5% quantile" = round(unname(summary(mr.coda[, "theta"])$quantiles[5]), 3),
                           "Rhat" = round(unname(gelman.diag(mr.coda)$psrf[1]), 3))
  names(mr_estimate) = c("Estimate", "SD", "2.5% quantile", "97.5% quantile", "Rhat")
  return(list("MR_Estimate" = mr_estimate, "MR_Coda" = mr.coda))
}

mr_horse_model = function() {
  for (i in 1:N){
    by[i] ~ dnorm(mu[i], 1/(sy[i] * sy[i]))
    mu[i] = theta * bx0[i] + alpha[i]
    bx[i] ~ dnorm(bx0[i], 1 / (sx[i] * sx[i]))
    
    bx0[i] ~ dnorm(mx0 + (sqrt(vx0)/(tau * phi[i])) * rho[i] * alpha[i], 1 / ((1 - rho[i]^2) * vx0))
    r[i] ~ dbeta(10, 10);T(, 1)
    rho[i] = 2*r[i] -1
    
    alpha[i] ~ dnorm(0, 1 / (tau * tau * phi[i] * phi[i]))
    phi[i] = a[i] / sqrt(b[i])
    a[i] ~ dnorm(0, 1);T(0, )
    b[i] ~ dgamma(0.5, 0.5)
  }
  
  c ~ dnorm(0, 1);T(0, )
  d ~ dgamma(0.5, 0.5)
  tau = c / sqrt(d)
  
  vx0 ~ dnorm(0, 1);T(0, )
  mx0 ~ dnorm(0, 1)
  
  theta ~ dunif(-10, 10)
}

loop_methods <- vector('list', iter)
loop_incidence_GWAS <- vector('list', iter)
loop_progression_GWAS <- vector('list', iter)
incidence_rsq <- matrix(0, iter)
progression_rsq <- matrix(0, iter)
correlation_coef <- matrix(0, iter)
incidence_mean <- matrix(0, iter)
progression_mean <- matrix(0, iter)
ci_coverage <- matrix(0, iter)

for (i in 1:iter) {

G <- matrix(0, n, nSNPs)
incidence_GWAS <- matrix(0, nSNPs, 3)
progression_GWAS <- matrix(0, nSNPs, 3)

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
  
##Overwrite P as we generated it for all 'individuals' but obviously don't want it for those who don't have the disease!
P[I == 0] <- NA

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

collider_bias_results <- data.frame()

##Run the methods

##Dudbridge method (updated CWLS from Cai et al. paper), based on April Hartley's code.

ivw <- mr_ivw(incidence_GWAS[, 1], progression_GWAS[, 1], 
              incidence_GWAS[, 2], progression_GWAS[, 2])
dudbridgeweights <- 1/progression_GWAS[, 2]^2
weighting <- (sum(dudbridgeweights*incidence_GWAS[, 1]^2))/
            ((sum(dudbridgeweights*incidence_GWAS[, 1]^2))-(sum(dudbridgeweights*incidence_GWAS[, 2]^2)))
cf.db <- ivw$Estimate*weighting
cf.se.db <- ivw$StdError*weighting

collider_bias_results <- rbind(collider_bias_results,
                               data.frame(Method = 'Dudbridge',
                               Correction_Beta = cf.db,
                               Correction_SE = cf.se.db,
                               Correction_Beta_SD <- sd(cf.db),
                               MCSE <- Correction_Beta_SD / sqrt(iter) 
                               ))
  
##Weighted median

Weighted_Median_Res <- mr_median((mr_input(bx = incidence_GWAS[, 1], by = progression_GWAS[, 1], 
                                           bxse = incidence_GWAS[, 2], byse = progression_GWAS[, 2])))

Weighted_Median_Beta <- Weighted_Median_Res$Estimate
Weighted_Median_SE <- Weighted_Median_Res$StdError
  
collider_bias_results <- rbind(collider_bias_results,
                                 data.frame(Method = 'Weighted_median',
                                 Correction_Beta = Weighted_Median_Beta,
                                 Correction_SE = Weighted_Median_SE,
                                 Correction_Beta_SD <- sd(Weighted_Median_Beta),
                                 MCSE <- Correction_Beta_SD / sqrt(iter) 
                                 ))
  
##MR-RAPS

MR_RAPS_Res <- mr.raps(incidence_GWAS[, 1], incidence_GWAS[, 2], 
                       progression_GWAS[, 1], progression_GWAS[, 2])

MR_RAPS_Beta <- MR_RAPS_Res$beta.hat
MR_RAPS_SE <- MR_RAPS_Res$beta.se
  
collider_bias_results <- rbind(collider_bias_results,
                                 data.frame(Method = 'MR_RAPS',
                                 Correction_Beta = MR_RAPS_Beta,
                                 Correction_SE = MR_RAPS_SE,
                                 Correction_Beta_SD <- sd(MR_RAPS_Beta),
                                 MCSE <- Correction_Beta_SD / sqrt(iter)               
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
                               Correction_SE = MR_Horse_SE,
                               Correction_Beta_SD <- sd(MR_Horse_Beta),
                               MCSE <- Correction_Beta_SD / sqrt(iter)
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
                               Correction_SE = Slopehunter_SE,
                               Correction_Beta_SD <- sd(Slopehunter_Beta),
                               MCSE <- Correction_Beta_SD / sqrt(iter)
                               ))

##Summarise methods results for each iteration
  
loop_methods[[i]] <- list(collider_bias_results = collider_bias_results)
loop_incidence_GWAS[[i]] <- list(incidence_GWAS = incidence_GWAS)
loop_progression_GWAS[[i]] <- list(progression_GWAS = progression_GWAS)

##Diagnostics

incidence_rsq[, i] <- summary(lm(I ~ G))$r.squared
progression_rsq[, i] <- summary(lm(P ~ G))$r.squared
correlation_coef[, i] <- cor(I,P)
incidence_mean[, i] <- mean(I)
progression_mean[, i] <- mean(P)
  
}

save(incidence_GWAS, progression_GWAS, loop_methods, incidence_rsq, progression_rsq, correlation_coef, incidence_mean, progression_mean, file = 'sim_results.RData')
