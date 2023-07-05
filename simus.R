# R packages
library(mvtnorm); library(mice); library(MatchIt); library(miceadds); library(tidyverse); library(sandwich)

# data generation function
gene_data <- function(n, confounding, beta0){
  
  # gaussian independent confounders
  X <- rmvnorm(n, mean = rep(0,3), sigma = diag(3))
  
  # treatment ; beta0 = -1, -1.4, -2.2 for 30%, 20% and 10% treated
  eta_ttmt <- 0
  if (confounding == "strong") eta_ttmt <- beta0 - 0.5*X[,1] - 0.4*X[,2] - 0.7*X[,3]
  if (confounding == "moderate") eta_ttmt <- beta0 - 0.3*X[,1] - 0.4*X[,2] - 0.3*X[,3]
  if (confounding == "weak") eta_ttmt <- beta0 - 0.01*X[,1] - 0.05*X[,2] + 0.01*X[,3]
  
  prob_ttmt <- exp(eta_ttmt)/(1+exp(eta_ttmt))
  Z <- rbinom(n, 1, prob = prob_ttmt)
  
  # outcome 
  Y <- rep(NA, n)
  eta_out <- 0; prob_out <- 0
  
    
  if (confounding == "strong") eta_out <- -1 + 0.4*X[,1] + 0.5*X[,2] + 0.9*X[,3] + 1.2*Z
  if (confounding == "moderate") eta_out <- -1 + 0.4*X[,1] + 0.5*X[,2] + 0.3*X[,3] + 1.2*Z
  if (confounding == "weak") eta_out <- -1 + 0.1*X[,1] + 0.1*X[,2] - 0.1*X[,3] + 3*Z
    
  prob_out <- exp(eta_out)/(1+exp(eta_out))
  Y <- rbinom(n, 1, prob = prob_out)

  
  # missing data
  eta_na2 <-  -2 + 0.1*X[,1] + 0.9*X[,3] + 1.1*Z
  prob_na2 <- exp(eta_na2)/(1+exp(eta_na2))
  R2 <- rbinom(n, 1, prob_na2)
  X2full <- X[,2]
  X[,2] <- ifelse(R2 == 1, NA, X[,2])
  
  dat <- as_tibble(cbind("id" = seq(n), "y" = Y, "x1" = X[,1], "x2" = X[,2], "x3" = X[,3], "z" = Z))
}

# multiple imputation, matching and treatment effect estimation using Rubin's rules
missmatch <- function(data, nb_imp){
  
  out <- NULL
  
  # multiple imputation ==========================================================================================
  imp_data <- mice(data, m = nb_imp, printFlag = FALSE)
  
  # ps and treatment effect estimation ===========================================================================
  matched_ids <- list()
    
  # method independent : matches for each imputed dataset independenty (20 sets of matches)
  for (i in seq(nb_imp)){
        
    # ps estimation
    dat <- mice::complete(imp_data,i)
    psmodel <- glm(z ~ x1 + x2 + x3, family = binomial(link = "logit"), data = dat)
    psi <- fitted.values(psmodel)
        
    # imputed dataset i with ps
    completedata <- cbind(mice::complete(imp_data, i), psi)
    colnames(completedata) <- c("id", "y", "x1", "x2", "x3", "z", "psi")
        
    # ps matching for imputed dataset i
    matched <- matchit(z ~ psi, distance = completedata$psi, data = completedata, method = "nearest", caliper = 0.2)
    matched_ids[[i]] <- as_tibble(rbind(cbind("id" = as.numeric(rownames(na.omit(matched$match.matrix))), "idpair" = seq(length(rownames(na.omit(matched$match.matrix))))),
                                        cbind("id"= as.numeric(na.omit(matched$match.matrix)), "idpair" = seq(length(rownames(na.omit(matched$match.matrix)))))))
  }
      
  # treatment effect estimation =================================================================================
  models_OR <- list()
  models_diff <- list()
      
  for (i in seq(length(matched_ids))){
          
    # extract matched data
    matched_data <- data %>% dplyr::select(c("id", "y", "z")) %>% inner_join(matched_ids[[i]], by = "id")
          
    # estimate models for OR and RD
    models_OR[[i]] <- glm.cluster(matched_data, y ~ z, matched_data$idpair, family="binomial")
    models_diff[[i]] <- glm.cluster(matched_data, y ~ z, matched_data$idpair, family="gaussian")
  }
      
  # agregation using Rubin's rule ===============================================================================
  nb_matchs <- length(matched_ids)
  tteff_diff=0;var_tteff_diff_within=0;var_tteff_diff_between=0;var_tteff_diff=0;
      
  # risk difference
  tteff_diff <- mean(as.numeric(unlist(lapply(models_diff, function(x) return(x$glm_res$coefficients[2])))))
  var_tteff_diff_within <- mean(as.numeric(unlist(lapply(models_diff, function(x) return(x$vcov[2,2])))))
  var_tteff_diff_between <- sum((as.numeric(unlist(lapply(models_diff, function(x) return(x$glm_res$coefficients[2])))) - tteff_diff)^2)/(nb_matchs-1)
  var_tteff_diff <- var_tteff_diff_within + ifelse(nb_matchs > 1, (1 + 1/nb_matchs), 0)*var_tteff_diff_between
      
  # odds ratio
  tteff_OR <- mean(as.numeric(unlist(lapply(models_OR, function(x) return(x$glm_res$coefficients[2])))))
  var_tteff_OR_within <- mean(as.numeric(unlist(lapply(models_OR, function(x) return(x$vcov[2,2])))))
  var_tteff_OR_between <- sum((as.numeric(unlist(lapply(models_OR, function(x) return(x$glm_res$coefficients[2])))) - tteff_OR)^2)/(nb_matchs-1)
  var_tteff_OR <- var_tteff_OR_within + ifelse(nb_matchs > 1, (1 + 1/nb_matchs), 0)*var_tteff_OR_between
  
  out <- list(tteff_OR, var_tteff_OR, tteff_diff, var_tteff_diff)
      
  # output =======================================================================================================
  return(out)
  
}

# multiple imputation, matching and treatment effect estimation using Reiter's rules
missmatch_reiter <- function(data, nb_rep, nb_imp){
  
  out <- NULL
  
  # repeating data nb_rep times ==================================================================================
  datarep <- data[rep(seq_len(nrow(data)), nb_rep),]
  
  # multiple imputation (nb_rep times for each parameter draw from 1 to nb_imp) ==================================
  imp_data <- mice(datarep, m = nb_imp, printFlag = FALSE, ignore = c(rep(FALSE, nrow(data)), rep(TRUE, nrow(datarep)-nrow(data))))
  
  # ps estimation ================================================================================================
  ps <- NULL;
  for (i in seq(nb_imp)){ # for each parameter draw
    psi <- NULL;
    dati <- mice::complete(imp_data, i)
    for (j in seq(nb_rep)){ # for each imputed sub dataset
      datij <- dati[seq_len(nrow(data))+nrow(data)*(j-1), ]
      psmodelij <- glm(z ~ x1 + x2 + x3, family = binomial(link = "logit"), data = datij)
      psij <- fitted.values(psmodelij)
      psi <- c(psi, psij)
    }
    ps <- cbind(ps, psi)
  }
  
  # matching =============================================================================================
  matched_ids <- list()
      
  for (i in seq(nb_imp)){
        
    matched_ids[[i]] <- list()
        
    # imputed dataset i with ps
    completedatai <- cbind(mice::complete(imp_data, i), ps[,i])
    colnames(completedatai) <- c("id", "y", "x1", "x2", "x3", "z", "psi")
        
    for (j in seq(nb_rep)){
          
      # imputed subdataset ij with ps
      completedataij <- completedatai[seq_len(nrow(data))+nrow(data)*(j-1), ]
      matched <- matchit(z ~ psi, distance = completedataij$psi, data = completedataij, method = "nearest", caliper = 0.2)
      matched_ids[[i]][[j]] <- as_tibble(rbind(cbind("id" = as.numeric(rownames(na.omit(matched$match.matrix)))-((j-1)*nrow(data)), "idpair" = seq(length(rownames(na.omit(matched$match.matrix))))),
                                               cbind("id"= as.numeric(na.omit(matched$match.matrix))-((j-1)*nrow(data)), "idpair" = seq(length(rownames(na.omit(matched$match.matrix)))))))
      
    }
  }
    
  # treatment effect estimation =================================================================================
  models_OR <- list()
  models_diff <- list()
    
  for (i in seq(length(matched_ids))){
        
    models_OR[[i]] <- list()
    models_diff[[i]] <- list()
        
    for (j in seq(nb_rep)){
      # extract matched data
      matched_data <- data %>% dplyr::select(c("id", "y", "z")) %>% inner_join(matched_ids[[i]][[j]], by = "id")
          
      # estimate models for OR and RD
      models_OR[[i]][[j]] <- glm.cluster(matched_data, y ~ z, matched_data$idpair, family="binomial")
      models_diff[[i]][[j]] <- glm.cluster(matched_data, y ~ z, matched_data$idpair, family="gaussian")
    }
  }
    
  # agregation using Reiter's rule ===============================================================================
  
  # risk difference
  tteff_diff <- mean(as.numeric(unlist(lapply(models_diff, function(x) return(as.numeric(unlist(lapply(x, function(x) return(x$glm_res$coefficients[2])))))))))
  var_tteff_diff_w <- mean(as.numeric(unlist(lapply(models_diff, function(x) return(var(as.numeric(unlist(lapply(x, function(x) return(x$glm_res$coefficients[2]))))))))))
  var_tteff_diff_b <- var(as.numeric(unlist(lapply(models_diff, function(x) return(mean(as.numeric(unlist(lapply(x, function(x) return(x$glm_res$coefficients[2]))))))))))
  var_tteff_diff_u <- mean(as.numeric(unlist(lapply(models_diff, function(x) return(as.numeric(unlist(lapply(x, function(x) return(x$vcov[2,2])))))))))
  var_tteff_diff <- var_tteff_diff_u - (1+1/nb_rep)*var_tteff_diff_w + (1+1/nb_imp)*var_tteff_diff_b
    
  # odds ratio
  tteff_OR <- mean(as.numeric(unlist(lapply(models_OR, function(x) return(as.numeric(unlist(lapply(x, function(x) return(x$glm_res$coefficients[2])))))))))
  var_tteff_OR_w <- mean(as.numeric(unlist(lapply(models_OR, function(x) return(var(as.numeric(unlist(lapply(x, function(x) return(x$glm_res$coefficients[2]))))))))))
  var_tteff_OR_b <- var(as.numeric(unlist(lapply(models_OR, function(x) return(mean(as.numeric(unlist(lapply(x, function(x) return(x$glm_res$coefficients[2]))))))))))
  var_tteff_OR_u <- mean(as.numeric(unlist(lapply(models_OR, function(x) return(as.numeric(unlist(lapply(x, function(x) return(x$vcov[2,2])))))))))
  var_tteff_OR <-  var_tteff_OR_u - (1+1/nb_rep)*var_tteff_OR_w + (1+1/nb_imp)*var_tteff_OR_b
  
  out <- list(tteff_OR, var_tteff_OR, tteff_diff, var_tteff_diff)
  
  # output =======================================================================================================
  return(out)
}

# simulations
simus <- function(N, n, confounding, beta0, nb_rep, nb_imp, reiter){
  
  set.seed(1993)
  out <- NULL; data <- NULL;
  
  # loop over the N replicates
  for (i in seq(N)){
    
    # generate data
    data <- gene_data(n, confounding, beta0)
    
    # imputation and matching with Reiter's rules
    if (reiter){
      tryi <- try(missmatch_reiter(data, nb_rep, nb_imp))
      while (class(tryi) == "try-error"){
        data <- gene_data(n, confounding, beta0)
        tryi <- try(missmatch_reiter(data, nb_rep, nb_imp))
      }
    }
    
    # imputation and matching with Rubin's rules
    else {
      tryi <- try(missmatch(data, nb_imp))
      while (class(tryi) == "try-error"){
        data <- gene_data(n, confounding, beta0)
        tryi <- try(missmatch(data, nb_imp))
      }
    }
    
    out <- rbind(out, tryi)
  
  }
  
  # output
  return(out)
}

# example of a simulation run
res <- simus(10, 1000, "strong", -1, 10, 10, TRUE)