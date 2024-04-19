##Dudbridge method, based on April Hartley's code.

ivw <- mr_ivw(incidence_GWAS[, 1], progression_GWAS[, 1], 
                 incidence_GWAS[, 2], progression_GWAS[, 2])
dudbridgeweights <- 1/progression_GWAS[, 2]^2
weighting <- (sum(dudbridgeweights*incidence_GWAS[, 1]^2))/
            ((sum(dudbridgeweights*incidence_GWAS[, 1]^2))-(sum(dudbridgeweights*incidence_GWAS[, 2]^2)))
cf.db <- ivw$b*weighting
cf.se.db <- ivw$se*weighting

##Weighted median method, simplest to implement in Stephen Burgess' MendelianRandomization package.

library(MendelianRandomization)
Weighted_Median_Res <- mr_median((mr_input(bx = incidence_GWAS$Estimate, bxse = incidence_GWAS$StdErr,
                                       by = progression_GWAS$Estimate, byse = progression_GWAS$StdErr)))

##MR-RAPS

library(mr.raps)
MR_RAPS_Res <- mr.raps(incidence_GWAS$Estimate, incidence_GWAS$StdErr, 
                       progression_GWAS$Estimate, progression_GWAS$StdErr)

##MR-Horse
##Functions to run the new method from Stephen Burgess and Andrew Grant, also requires JAGS to be installed (seperate to R).
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

##Reformat data and run the MR-Horse method.

MR_Horse_Data <- data.frame(incidence_GWAS$Estimate, progression_GWAS$Estimate, 
                            incidence_GWAS$StdErr, progression_GWAS$StdErr)
names(MR_Horse_Data)[1] <- "betaX"
names(MR_Horse_Data)[2] <- "betaY"
names(MR_Horse_Data)[3] <- "betaXse"
names(MR_Horse_Data)[4] <- "betaYse"
MR_Horse_Res <- mr_horse(MR_Horse_Data)

##Slope-Hunter
library(SlopeHunter)

##Reformat data and run the SlopeHunter method.
SlopeHunter_Data <- data.frame(incidence_GWAS$Estimate, incidence_GWAS$StdErr, 
                               progression_GWAS$Estimate, progression_GWAS$StdEr)
names(SlopeHunter_Data)[1] <- "xbeta_col"
names(SlopeHunter_Data)[2] <- "ybeta_col"
names(SlopeHunter_Data)[3] <- "xse_col"
names(SlopeHunter_Data)[4] <- "yse_col"

SlopeHunter_Res <- hunt(dat = SlopeHunter_Data, xbeta_col = "xbeta_col", 
                        xse_col = "xse_col", ybeta_col = "ybeta_col", yse_col = "yse_col")

##Summarise results
##Idea taken from Andrew Elmore.
library(dplyr)

collider_bias_results <- data.frame(
  Method = character(),
  Correction_Beta = numeric(),
  Correction_SE = numeric()
)

##Add more to this later

collider_bias_type <- list(
  Slopehunter = "Slopehunter"
  Dudbridge = "Dudbridge",
  Weighted_median = "Weighted_median",
  MR_Horse = "MR_Horse",
  MR_RAPS = "MR_RAPS"
)

##Add Slopehunter results

collider_bias_results <- dplyr::add_row(collider_bias_results,
                                                Method = collider_bias_type$Slopehunter,
                                                Correction_Beta = SlopeHunter_Res$b,
                                                Correction_SE = SlopeHunter_Res$bse
        )

##Add Dudbridge results

collider_bias_results <- dplyr::add_row(collider_bias_results,
                                                Method = collider_bias_type$Dudbridge,
                                                Correction_Beta = cf.db,
                                                Correction_SE = cf.se.db
        )

##Add Weighted Median results

collider_bias_results <- dplyr::add_row(collider_bias_results,
                                        Method = collider_bias_type$Weighted_median,
                                        Correction_Beta = Weighted_Median_Res$b,
                                        Correction_SE = Weighted_Median_Res$se
)

##Add MR-RAPS results

collider_bias_results <- dplyr::add_row(collider_bias_results,
                                        Method = collider_bias_type$MR_RAPS,
                                        Correction_Beta = MR_RAPS_Res$beta.hat,
                                        Correction_SE = MR_RAPS_Res$beta.se
)

##Add MR-Horse results
##Bayesian method so outputs SD instead of SE?

collider_bias_results <- dplyr::add_row(collider_bias_results,
                                        Method = collider_bias_type$MR_Horse,
                                        Correction_Beta = MR_Horse_Res$MR_Estimate$Estimate,
                                        Correction_SE = MR_Horse_Res$MR_Estimate$SD
)
