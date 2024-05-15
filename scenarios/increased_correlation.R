rm(list=ls())
setwd("/user/home/vc23656/Collider_sims")
.libPaths("/user/work/vc23656/")
expit <- function (x) exp(x) / (1 + exp(x))
set.seed(54545)
reps <- 1000
n.ind <- 1e5
nSNPs <- 100
incidence.SNPs <- 1:80
progression.SNPs <- 71:100
library(MendelianRandomization)
library(mr.raps)
library(mclust)
library(R2jags)
library(dplyr)

shclust <- function(gwas, pi0, sxy1){
  # binding variable locally to the function:
  ## To avoid Notes: e.g. "shclust: no visible binding for global variable ‘xbeta’"
  xbeta <- ybeta <- clusters <- NULL
  
  sx0 = sx1 = gwas %>% summarise(sd(xbeta)) %>% pull
  sy0 = sy1 = gwas %>% summarise(sd(ybeta)) %>% pull
  dir0 = gwas %>% summarise(cov(xbeta, ybeta)) %>% pull %>% sign()
  if (dir0==0) stop("All associations with at least either x or y are constant")
  
  # convergence criterion
  loglkl_ck = 0
  
  ### EM algorithm
  for(iter in 1:50000){
    #### The E step:
    # covariance matrix for the target component (f0)
    sxy0 = sx0*sy0*dir0*0.95       # the x & y perfectly correlated under 1st component  #===========
    sigma0 = matrix(c(sx0^2,sxy0,sxy0,sy0^2), 2, 2)
    
    # 1st component
    f0 = gwas %>%
      dplyr::select(xbeta, ybeta) %>%
      mclust::dmvnorm(mean=c(0,0), sigma=sigma0)
    f0[f0<1e-300] = 1e-300
    
    # covariance matrix for the component (f1)
    sigma1 = matrix(c(sx1^2,sxy1,sxy1,sy1^2), 2, 2)
    
    # 2nd component
    f1 = gwas %>%
      dplyr::select(xbeta, ybeta) %>%
      mclust::dmvnorm(mean=c(0,0), sigma=sigma1)
    f1[f1<1e-300] = 1e-300
    
    # loglik of the mixture model: pi0 * f0 + (1-p0) * f1
    loglkl = sum(log(pi0*f0+(1-pi0)*f1))
    
    ## proportional contribution of density of f0 (for every point) to the total mixture
    pt = pi0*f0/(pi0*f0+(1-pi0)*f1)
    pt[pt>0.9999999] = 0.9999999
    
    #### The M step:
    # update pi0
    pi0 = mean(pt)
    if (pi0<0.0001) pi0 = 0.0001
    if (pi0>0.9999) pi0 = 0.9999
    
    # update sx0 & sy0
    sx0 = gwas %>% summarise(sqrt(sum(pt*(xbeta^2))/sum(pt))) %>% pull
    sy0 = gwas %>% summarise(sqrt(sum(pt*(ybeta^2))/sum(pt))) %>% pull
    dir0 = gwas %>% summarise(sum(pt*xbeta*ybeta)/sum(pt)) %>% pull %>% sign()
    if (dir0==0) dir0=sample(c(1,-1), 1)   # avoid slope = 0 (horizontal line)
    
    # update sx1, sy1 & sxy1
    sx1 = gwas %>% summarise(sqrt(sum((1-pt)*(xbeta^2))/(length(xbeta)-sum(pt)))) %>% pull
    sy1 = gwas %>% summarise(sqrt(sum((1-pt)*(ybeta^2))/(length(ybeta)-sum(pt)))) %>% pull
    sxy1 = gwas %>% summarise(sum((1-pt)*xbeta*ybeta)/(length(xbeta)-sum(pt))) %>% pull
    if (abs(sxy1) > 0.75*sx1*sy1)  sxy1 = sign(sxy1)*0.75*sx1*sy1     #===========
    
    ## Check convergence
    if (iter%%10==0){
      if ((loglkl - loglkl_ck)/loglkl < 1e-10){
        break
      } else {
        loglkl_ck = loglkl
      }
    }
  }
  
  # Diagnosis
  if (iter == 50000) warning("The algorithm may not have converged.\n")
  
  ### Results
  Fit = gwas %>% mutate(pt = pt) %>%
    mutate(po = 1-pt, clusters = factor(ifelse(pt >= 0.5, "Hunted", "Pleiotropic")))
  
  # slope of the eigenvector
  b = dir0*sy0/sx0
  bse = 0
  b_CI = c(b - 1.96*bse, b + 1.96*bse)
  entropy = Fit %>% filter(clusters == "Hunted") %>% summarise(mean(pt)) %>% pull
  
  return(list(b=b, bse=bse, b_CI=b_CI, iter=iter, pi0=pi0, entropy=entropy, Fit=Fit))
}

hunt = function(dat, snp_col="SNP", xbeta_col="BETA.incidence", xse_col="SE.incidence", xp_col="Pval.incidence",
                ybeta_col="BETA.prognosis", yse_col="SE.prognosis", yp_col="Pval.prognosis",
                xp_thresh=0.001, init_pi = 0.6, init_sigmaIP = 1e-5, Bootstrapping = TRUE, M = 100,
                show_adjustments = FALSE){
  
  # binding variable locally to the function:
  ## To avoid Notes: e.g. "hunt: no visible binding for global variable ‘xp’"; "hunt: no visible binding for global variable ‘xbeta’"
  xp <- xbeta <- ybeta <- clusters <- yse <- xse <- ybeta_adj <- yse_adj <- NULL
  
  all_cols <- c(snp_col, xbeta_col, xse_col, xp_col, ybeta_col, yse_col, yp_col)
  i <-  names(dat) %in% all_cols
  if (sum(i) == 0)
  {
    stop("None of the specified columns present!")
  }
  dat <- dat[, i]
  
  # Check if columns required for SH are present
  cols_req <- c(xbeta_col, xse_col, ybeta_col, yse_col)
  if (!all(cols_req %in% names(dat)))
  {
    stop("The following columns are not present and are required for the Slope-Hunter analysis:\n", paste(cols_req[!cols_req %in% names(dat)]), collapse="\n")
  }
  
  # generate p-values and SNP IDs if are not given
  cols_desired <- c(xp_col, yp_col, snp_col)
  i <- cols_desired %in% names(dat)
  if (!all(i))
  {
    message("The following column(s) is not present and will be generated:\n", paste(cols_desired[!i]))
    if(!i[1]){dat[[xp_col]] = pchisq((dat[[xbeta_col]]/dat[[xse_col]])^2, 1, lower.tail = FALSE)}
    if(!i[2]){dat[[yp_col]] = pchisq((dat[[ybeta_col]]/dat[[yse_col]])^2, 1, lower.tail = FALSE)}
    if(!i[3]){dat[[snp_col]] = paste0("snp", 1:nrow(dat))}
  }
  
  names(dat)[names(dat) == snp_col] <- "SNP"
  names(dat)[names(dat) == xbeta_col] <- "xbeta"
  names(dat)[names(dat) == xse_col] <- "xse"
  names(dat)[names(dat) == xp_col] <- "xp"
  names(dat)[names(dat) == ybeta_col] <- "ybeta"
  names(dat)[names(dat) == yse_col] <- "yse"
  names(dat)[names(dat) == yp_col] <- "yp"
  
  # filter in the variants associated with x
  gwas <- dat[, c("SNP", "xbeta", "xse", "xp", "ybeta", "yse", "yp")] %>%
    filter(xp <= xp_thresh)
  
  # fit slope-hunter model
  Model = shclust(gwas, pi0=init_pi, sxy1=init_sigmaIP)
  b = Model$b
  
  # estimate bse
  if (Bootstrapping){
    b.bts = vector("numeric", M)
    for(i in 1:M){
      print(paste("Bootstrap sample", i, "of", M, "samples ..."))   ### replace it by a progress bar ...
      Bts = sample(1:nrow(gwas), size = nrow(gwas), replace = TRUE)
      DT = gwas[Bts,]
      Model.DT = shclust(DT, pi0=init_pi, sxy1=init_sigmaIP)
      b.bts[i] = Model.DT$b
    }
    
    # If any of the model fits for the bootstrap samples generated NA
    if(any(is.na(b.bts))){
      b.bts = na.omit(b.bts)
      warning(paste("Only", length(b.bts), "bootstrap samples - out of", M, "- produced converged models used for estimating the standard error." ))
    }
    
    # calculate bse
    bse = sd(abs(b.bts))
    b_CI = c(b - 1.96*bse, b + 1.96*bse)
  } else{
    bse = Model$bse
    b_CI = Model$b_CI
  }
  
  # model fit
  Fit = Model$Fit
  
  if(show_adjustments)
  {
    ##### Adjust
    est = dat %>% mutate(ybeta_adj = ybeta - (b * xbeta),
                         yse_adj = sqrt(yse^2 + (b^2 * xse^2) + (xbeta^2 * bse^2) + (xse^2 * bse^2)),
                         yp_adj = pchisq((ybeta_adj/yse_adj)^2, 1, lower.tail = FALSE))
  }
  
  # Print main results
  print(paste0("Estimated slope: ", b))
  print(paste0("SE of the slope: ", bse))
  print(paste0("95% CI: ", b_CI[1], ", ", b_CI[2]))
  
  # Return
  if(Bootstrapping & show_adjustments){
    SH = list(est=est, b=b, bse=bse, b_CI=b_CI, pi=Model$pi0, entropy=Model$entropy, Fit=Model$Fit, iter = Model$iter, Bts.est = b.bts)
  }
  if(Bootstrapping & !show_adjustments){
    SH = list(b=b, bse=bse, b_CI=b_CI, pi=Model$pi0, entropy=Model$entropy, Fit=Model$Fit, iter = Model$iter, Bts.est = b.bts)
  }
  if(!Bootstrapping & show_adjustments){
    SH = list(est=est, b=b, bse=bse, b_CI=b_CI, pi=Model$pi0, entropy=Model$entropy, Fit=Model$Fit, iter = Model$iter)
  }
  if(!Bootstrapping & !show_adjustments){
    SH = list(b=b, bse=bse, b_CI=b_CI, pi=Model$pi0, entropy=Model$entropy, Fit=Model$Fit, iter = Model$iter)
  }
  class(SH) <- "SH"
  return(SH)
}

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

loop_methods <- vector('list', reps)
loop_incidence_GWAS <- vector('list', reps)
loop_progression_GWAS <- vector('list', reps)
incidence_rsq <- matrix(0, reps)
progression_rsq <- matrix(0, reps)
collider_rsq <- matrix(0, reps)
correlation_coef <- matrix(0, reps)
incidence_mean <- matrix(0, reps)
progression_mean <- matrix(0, reps)

for (n in 1:reps) {
  
  G <- matrix(0, n.ind, nSNPs)
  
  ##Simulate common cause (collider).
  U <- rnorm(n.ind, 0, 1)
  
  ##Simulate one-sample genetic data from individual-level data - number of SNPs, generate MAF, ratios for how many affect I and P.
  maf <- runif(nSNPs, 0.05, 0.5)
  for (j in 1:nSNPs) G[, j] <- rbinom(n.ind, 2, maf[j])
  incidence.betas <- rep(0, nSNPs)
  incidence.betas[incidence.SNPs] <- rnorm(length(incidence.SNPs), 0, 0.2)
  incidence.probs <- expit(-1 + as.vector(G %*% incidence.betas) + 1 * U)
  I <- rbinom(n.ind, 1, incidence.probs)
  
  progression.betas <- rep(0, nSNPs)
  progression.betas[progression.SNPs] <- rnorm(length(progression.SNPs), 0, 0.2)
  progression.probs <- expit(-1 + as.vector(G %*% progression.betas) + 1 * U)
  P <- rbinom(n.ind, 1, progression.probs)
  
  ##Simulate GWAS of incidence.
  
  incidence_GWAS <- data.frame()
  progression_GWAS <- data.frame()
  collider_bias_results <- data.frame()
  
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
  
  ##Run the methods
  
  ##Dudbridge method (updated CWLS from Cai et al. paper), based on April Hartley's code.
  
  ivw <- mr_ivw((mr_input(bx = incidence_GWAS[, 1], by = progression_GWAS[, 1], 
                          bxse = incidence_GWAS[, 2], byse = progression_GWAS[, 2])))
  dudbridgeweights <- 1/progression_GWAS[, 2]^2
  weighting <- (sum(dudbridgeweights*incidence_GWAS[, 1]^2))/
    ((sum(dudbridgeweights*incidence_GWAS[, 1]^2))-(sum(dudbridgeweights*incidence_GWAS[, 2]^2)))
  cf.db <- ivw$Estimate*weighting
  cf.se.db <- ivw$StdError*weighting
  
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
  
  loop_methods[[n]] <- list(collider_bias_results = collider_bias_results)
  loop_incidence_GWAS[[n]] <- list(incidence_GWAS = incidence_GWAS)
  loop_progression_GWAS[[n]] <- list(progression_GWAS = progression_GWAS)
  
  ##Diagnostics
  
  incidence_rsq[n] <- summary(lm(I ~ G))$r.squared
  progression_rsq[n] <- summary(lm(P ~ G))$r.squared
  collider_rsq[n] <- summary(lm(U ~ G))$r.squared
  correlation_coef[n] <- cor(I,P)
  incidence_mean[n] <- mean(I)
  progression_mean[n] <- mean(P)

  save(loop_incidence_GWAS, loop_progression_GWAS, loop_methods, incidence_rsq, progression_rsq, collider_rsq, correlation_coef, incidence_mean, progression_mean, file = 'sim_results.RData')

}

