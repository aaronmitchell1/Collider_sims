##Dudbridge method, based on April Hartley's code.

incidence_GWAS <- data.frame(incidence_GWAS)
progression_GWAS <- data.frame(progresion_GWAS)
ivw <- mr_ivw(incidence_GWAS$Estimate, progression_GWAS$Estimate, incidence_GWAS$StdErr, progression_GWAS$StdErr)
dudbridgeweights <- 1/progression_GWAS$StdErr^2
weighting <- (sum(dudbridgeweights*incidence_GWAS$Estimate^2))/((sum(dudbridgeweights*incidence_GWAS$Estimate^2))-(sum(dudbridgeweights*incidence_GWAS$StdErr^2)))
cf.db <- ivw$b*weighting
cf.se.db <- ivw$se*weighting

##Weighted median method, simplest to implement in Stephen Burgess' MendelianRandomization package.

library(MendelianRandomization)
weighted_median <- mr_median((mr_input(bx = incidence_GWAS$Estimate, bxse = incidence_GWAS$StdErr,
                                       by = progression_GWAS$Estimate, byse = progression_GWAS$StdErr)))

##MR-RAPS

library(mr.raps)
MR_RAPS_res <- mr.raps.all(incidence_GWAS$Estimate, incidence_GWAS$StdErr, progression_GWAS$Estimate, progression_GWAS$StdErr)

##MR-Horse
##Function to run the new method taken from Stephen Burgess and Andrew Grant, requires JAGS to be installed.

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

##Reformat data and run the MR-Horse method.

D <- data.frame(incidence_GWAS$Estimate, progression_GWAS$Estimate, incidence_GWAS$StdErr, progression_GWAS$StdErr)
names(D)[1] <- "betaX"
names(D)[2] <- "betaY"
names(D)[3] <- "betaXse"
names(D)[4] <- "betaYse"
MR_Horse_res <- mr_horse(D)

##Slope-Hunter

