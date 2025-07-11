#' Function for one bootstrap replicate
#'
#' @param data dataset containing vaccine indicator, time variable, event type variable, and covariates
#' @param asymp_formula formula for asymptomatic Cox model
#' @param symp_formula formula for symptomatic Cox model
#' @param symp_lr_formula formula for LR GLM
#' @param n_boot number of bootstrap replicates for bootstrap CI
#' @param t0 numeric value or numeric vector of t0s
#' @param event_name name of event variable
#'
#' @returns dataframe with VE and P estimates
#' @keywords internal
one_boot <- function(data, asymp_formula, symp_formula, symp_lr_formula, n_boot, t0, event_name){
  
  # Create bootstrap data
  boot_id <- sample(1:nrow(data), replace = TRUE)
  boot_data <- data[boot_id, ,drop = FALSE]

  # Fit cox models for symptomatic & asymptomatic infections, symptomatic GLM
  symp_cox_fit <- survival::coxph(symp_formula, data = boot_data, model = TRUE)
  
  asymp_cox_fit <- survival::coxph(asymp_formula, data = boot_data, model = TRUE)
  
  symp_lr_fit <- glm(symp_lr_formula, 
                     # subset to infected participants only (given uninfected == 0)
                     data = boot_data[boot_data[[event_name]] > 0, ],
                     family = stats::binomial()
  )
  
  # Get point estimates
  ve_fit <- ve_from_fits(
    symp_cox_fit = symp_cox_fit, 
    asymp_cox_fit = asymp_cox_fit, 
    symp_lr_fit = symp_lr_fit, 
    data = boot_data, 
    t0 = t0
  )
  
  # Reformat results to get SE for VE_X and Ps more easily
  results_df <- data.frame(t0 = ve_fit$t0,
                           ve_i = ve_fit$ve_i,
                           ve_s = ve_fit$ve_s,
                           ve_ai = ve_fit$ve_ai,
                           ve_nai = ve_fit$ve_nai,
                           p_immune = ve_fit$p_immune,
                           p_doomed = ve_fit$p_doomed,
                           p_alwaysinf = ve_fit$p_alwaysinf,
                           p_converted = ve_fit$p_converted,
                           p_helped = ve_fit$p_helped,
                           p_helpedplus = ve_fit$p_helpedplus)
  
  return(results_df)
}

#' Function for n_boot bootstrap replicate and standard error/CIs
#'
#' @param data dataset containing vaccine indicator, time variable, event type variable, and covariates
#' @param asymp_formula formula for asymptomatic Cox model
#' @param symp_formula formula for symptomatic Cox model
#' @param symp_lr_formula formula for LR GLM
#' @param n_boot number of bootstrap replicates for bootstrap CI
#' @param t0 numeric value or numeric vector of t0s
#' @param event_name name of event variable
#'
#' @returns list of standard error, lower 95% CI bound, and upper 95% CI bound for each VE estimate and stratum probability
#' @keywords internal
bootstrap_estimates <- function(data, asymp_formula, symp_formula, symp_lr_formula, n_boot, t0, event_name){
  
  # Do n_boot bootstrap replicates
  boot_estimates <- replicate(n_boot, one_boot(data = data, 
                                               asymp_formula = asymp_formula, 
                                               symp_formula = symp_formula, 
                                               symp_lr_formula = symp_lr_formula,
                                               n_boot = n_boot,
                                               t0 = t0,
                                               event_name = event_name), simplify = FALSE) 
  
  boot_res <- data.frame(do.call(rbind, boot_estimates))
  
  # Full results list
  out <- list()
  
  # For each threshold:
  for(t in unique(boot_res$t0)){
    boot_res_t <- boot_res[boot_res$t0 == t,]
    
    t0_out <- list()
    
    # Get se, lower bound, upper bound for each column in boot_res for the given threshold
    for(col in colnames(boot_res_t)[-1]){
      se <- sd(boot_res_t[,col])
      lower <- quantile(boot_res_t[,col], p = 0.025, names = FALSE)
      upper <- quantile(boot_res_t[,col], p = 0.975, names = FALSE)
      
      t0_out[[col]] <- list(se = se, lower = lower, upper = upper)
    }
    
    # Add results for given threshold to full results list
    out[[paste0("t0_", t)]] <- t0_out
    
  }
  
  return(out)
  
}

