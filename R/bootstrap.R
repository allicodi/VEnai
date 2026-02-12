#' Function for one bootstrap replicate
#'
#' @param data dataset containing vaccine indicator, time variable, event type variable, and covariates
#' @param asymp_formula formula for asymptomatic Cox model
#' @param symp_formula formula for symptomatic Cox model
#' @param symp_lr_formula formula for LR GLM
#' @param n_boot number of bootstrap replicates for bootstrap CI
#' @param t0 numeric value or numeric vector of t0s
#' @param event_name name of event variable
#' @param weight_name name of weight variable, if applicable
#'
#' @returns dataframe with VE and P estimates
#' @keywords internal
one_boot <- function(data, asymp_formula, symp_formula, symp_lr_formula, 
                     symp_level, asymp_level, n_boot, t0, 
                     event_name, weight_name, vax_name, time_name, bounds,
                     sensitivity, delta, delta_X_variable){
  
  # Create bootstrap data
  boot_id <- sample(1:nrow(data), replace = TRUE)
  boot_data <- data[boot_id, ,drop = FALSE]

  # Fit cox models for symptomatic & asymptomatic infections, symptomatic GLM
  if(is.null(weight_name)){
    symp_cox_fit <- survival::coxph(symp_formula, data = boot_data, model = TRUE)
    
    asymp_cox_fit <- survival::coxph(asymp_formula, data = boot_data, model = TRUE)
    
    symp_lr_fit <- glm(symp_lr_formula, 
                       # subset to infected participants only (given uninfected == 0)
                       data = boot_data[boot_data[[event_name]] > 0, ],
                       family = stats::binomial()
    )
  } else{
    # get weights externally bc coxph env issue
    boot_weights <- boot_data[[weight_name]]
    
    #symp_cox_fit <- survival::coxph(symp_formula, data = boot_data, weights = boot_weights, model = TRUE)
    
    #asymp_cox_fit <- survival::coxph(asymp_formula, data = boot_data, weights = boot_weights, model = TRUE)
    
    # symp_lr_fit <- glm(symp_lr_formula, 
    #                    # subset to infected participants only (given uninfected == 0)
    #                    data = boot_data[boot_data[[event_name]] > 0, ],
    #                    weights = boot_weights[boot_data[[event_name]] > 0],
    #                    family = stats::binomial())
    
    # chatgpt trick bc environment issue with weights??
    symp_cox_fit <- eval(
      substitute(
        survival::coxph(symp_formula, data = boot_data, weights = W, model = TRUE),
        list(W = boot_weights)
      )
    )
    
    asymp_cox_fit <- eval(
      substitute(
        survival::coxph(asymp_formula, data = boot_data, weights = W, model = TRUE),
        list(W = boot_weights)
      )
    )

    symp_lr_fit <- eval(
      substitute(
        glm(symp_lr_formula, 
           # subset to infected participants only (given uninfected == 0)
           data = boot_data[boot_data[[event_name]] > 0, ],
           weights = W,
           family = stats::binomial()),
        list(W = boot_weights[boot_data[[event_name]] > 0])),
    )
  }
  
  # Get point estimates
  ve_fit <- ve_from_fits(
    symp_cox_fit = symp_cox_fit, 
    asymp_cox_fit = asymp_cox_fit, 
    symp_lr_fit = symp_lr_fit, 
    data = boot_data, 
    t0 = t0,
    vax_name = vax_name
  )
  
  if(bounds){
    ve_nai_bounds_res <- ve_nai_bounds(data = boot_data,
                                        vax_name = vax_name,
                                        time_name = time_name,
                                        event_name = event_name,
                                        symp_level = symp_level,
                                        asymp_level = asymp_level,
                                        t0 = t0)
  } else{
    ve_nai_bounds_res <- NULL
  }
  
  if(sensitivity){
    ve_nai_sens_res <- ve_nai_sensitivity(data = boot_data,
                                         symp_cox_fit = symp_cox_fit,
                                         asymp_cox_fit = asymp_cox_fit,
                                         symp_lr_fit = symp_lr_fit,
                                         t0 = t0,vax_name = vax_name,
                                         symp_ind_name = symp_ind_name,
                                         delta = delta,
                                         delta_X_variable = delta_X_variable)
  } else{
    ve_nai_sens_res <- NULL
  }
  
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
  
  return(list(results_df = results_df,
              ve_nai_bounds = ve_nai_bounds_res,
              ve_nai_sens = ve_nai_sens_res))
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
#' @param weight_name name of weight variable, if applicable
#'
#' @returns list of standard error, lower 95% CI bound, and upper 95% CI bound for each VE estimate and stratum probability
#' @keywords internal
bootstrap_estimates <- function(data, asymp_formula, symp_formula, symp_lr_formula, 
                                symp_level, asymp_level,n_boot, t0, event_name, 
                                weight_name, vax_name, time_name, bounds, sensitivity, delta, delta_X_variable){
  
  # Do n_boot bootstrap replicates
  boot_estimates <- replicate(n_boot, one_boot(data = data, 
                                               asymp_formula = asymp_formula, 
                                               symp_formula = symp_formula, 
                                               symp_lr_formula = symp_lr_formula,
                                               symp_level = symp_level,
                                               asymp_level = asymp_level, 
                                               n_boot = n_boot,
                                               t0 = t0,
                                               event_name = event_name,
                                               weight_name = weight_name,
                                               vax_name = vax_name, 
                                               time_name = time_name,
                                               bounds = bounds,
                                               sensitivity = sensitivity,
                                               delta = delta,
                                               delta_X_variable = delta_X_variable), simplify = FALSE) 
  
  # Extract and rbind 
  results_df_list <- lapply(boot_estimates, `[[`, "results_df")
  results_df_all <- do.call(rbind, results_df_list)

  if(bounds){
    ve_nai_bounds_list <- lapply(boot_estimates, `[[`, "ve_nai_bounds")
    ve_nai_bounds_all <- do.call(rbind, ve_nai_bounds_list)
  }
  
  if(sensitivity){
    ve_nai_sens_list <- lapply(boot_estimates, `[[`, "ve_nai_sens")
    ve_nai_sens_all <- do.call(rbind, ve_nai_sens_list)
  }
  
  # Full results list
  out <- list()
  
  # For each threshold:
  for(t in unique(t0)){
    boot_res_t <- results_df_all[results_df_all$t0 == t,]
    
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
    
    # Repeat above for lower and upper bound on each bound 
    if(bounds){
      bounds_t <- ve_nai_bounds_all[ve_nai_bounds_all$t0 == t,]
      bounds_out <- list()
      
      for(col in c("lower_bound", "upper_bound")){
        se <- sd(bounds_t[,col])
        lower <- quantile(bounds_t[,col], p = 0.025, names = FALSE)
        upper <- quantile(bounds_t[,col], p = 0.975, names = FALSE)
        
        bounds_out[[col]] <- list(se = se, lower = lower, upper = upper)
      }
      
      out[[paste0("ve_nai_bound_t0_", t)]] <- bounds_out
    }
    
    # Repeat above for lower and upper bound on VE_nai for each delta
    if(sensitivity){
      sens_t <- ve_nai_sens_all[ve_nai_sens_all$t0 == t,]
      sens_out <- data.frame()
      
      for(d in 1:length(delta)){
        d0 <- delta[d]
        se <- sd(sens_t$ve_nai[sens_t$delta == d0])
        lower <- quantile(sens_t$ve_nai[sens_t$delta == d0], p = 0.025, names = FALSE)
        upper <- quantile(sens_t$ve_nai[sens_t$delta == d0], p = 0.975, names = FALSE)
        
        sens_out <- rbind(sens_out, data.frame(delta = d0, se = se, lower = lower, upper = upper))
      }
      
      out[[paste0("ve_nai_sens_t0_", t)]] <- sens_out
    }
    
  }
  
  return(out)
  
}

