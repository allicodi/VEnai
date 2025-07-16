
#' Function to estimate the bounds on VE_NAI
#' 
#' @param data The \code{data} used to fit the above models. Should have columns \code{ftime}, \code{ftype}, 
#' 						 \code{vax}, and other columns for covariates.
#' @param vax_name name of vaccine variable
#' @param time_name name of time of event variable
#' @param event_name name of event variable
ve_nai_bounds <- function(data,
                          vax_name = "vax",
                          time_name = "ftime",
                          event_name = "ftype",
                          symp_level = 1, 
                          asymp_level = 2,
                          t0 = quantile(data[[time_name]], c(0.2, 0.5, 0.9))){
  
  # Estimates of cumulative incidence
  est <- cmprsk::cuminc(ftime = data[[time_name]],
                        fstatus = data[[event_name]],
                        group = data[[vax_name]])
  
  bound_est <- data.frame(t0 = t0, 
                          lower_bound = NA,
                          upper_bound = NA)
  
  # Get relevant probabilities and build bounds for each t in t0
  for(time in t0){
    # get estimates at given t0 = time
    est_t <- timepoints(est, time)$est
    
    # convert row names to columns & fill in values
    split_names <- do.call(rbind, strsplit(rownames(est_t), " ", fixed = TRUE))
    colnames(split_names) <- c("vax", "ftype")
    est_t_df <- cbind(data.frame(split_names), 
                      data.frame(cuminc = as.numeric(est_t)))
    
    # P(Y_I = 1, Y_S = 0 | V = 1) = proportion of asymptomatic infections in vaccinated
    P_Inf_1_Symp_0__V_1 <- est_t_df$cuminc[est_t_df$vax == 1 & est_t_df$ftype == asymp_level]
    
    # P(Y_I = 1, Y_S = 0 | V = 0) = proportion of asymptomatic infections in unvaccinated
    P_Inf_1_Symp_0__V_0 <- est_t_df$cuminc[est_t_df$vax == 0 & est_t_df$ftype == asymp_level]
    
    # P(Y_I = 0, Y_S = 0 | V = 1) = proportion of population that is AP, BP, or NI
    P_Inf_0_Symp_0__V_1 <- 1 - (est_t_df$cuminc[est_t_df$vax == 1 & est_t_df$ftype == asymp_level] + 
                                  est_t_df$cuminc[est_t_df$vax == 1 & est_t_df$ftype == symp_level])
    
    # P(Y_I = 0, Y_S = 0 | V = 0) = proportion of population that is NI
    P_Inf_0_Symp_0__V_0 <- 1 - (est_t_df$cuminc[est_t_df$vax == 0 & est_t_df$ftype == asymp_level] + 
                                  est_t_df$cuminc[est_t_df$vax == 0 & est_t_df$ftype == symp_level])
    
    ## Lower bound: VE_AI = 1 - [P(Y_I = 1, Y_S = 0 | V = 1) / P(Y_I = 1, Y_S = 0 | V = 0)]
    lower_bound <- 1 - (P_Inf_1_Symp_0__V_1 / P_Inf_1_Symp_0__V_0)
    
    ## Upper bound: [P(Y_I = 0, Y_S = 0 | V = 1) - P(Y_I = 0, Y_S = 0 | V = 0)] / P(Y_I = 1, Y_S = 0 | V = 0)
    upper_bound <- (P_Inf_0_Symp_0__V_1 - P_Inf_0_Symp_0__V_0) / P_Inf_1_Symp_0__V_0
    
    bound_est$lower_bound[bound_est$t0 == time] <- lower_bound
    bound_est$upper_bound[bound_est$t0 == time] <- upper_bound
    
  }
  
  return(bound_est)
  
  # --------
  # sub_V_1 <- data[data[[vax_name]] == 1]
  # sub_V_0 <- data[data[[vax_name]] == 0]
  # 
  # ## Lower bound: VE_AI = 1 - [P(Y_I = 1, Y_S = 0 | V = 1) / P(Y_I = 1, Y_S = 0 | V = 0)]
  # 
  # # P(Y_I = 1, Y_S = 0 | V = 1) = proportion of asymptomatic infections in vaccinated
  # P_Inf_1_Symp_0__V_1 <- length(which(sub_V_1[[event_name]] == asymp_level)) / nrow(sub_V_1)
  #   
  # # P(Y_I = 1, Y_S = 0 | V = 0) = proportion of asymptomatic infections in unvaccinated
  # P_Inf_1_Symp_0__V_0 <- length(which(sub_V_0[[event_name]] == asymp_level)) / nrow(sub_V_0)
  #   
  # lower_bound <- 1 - (P_Inf_1_Symp_0__V_1 / P_Inf_1_Symp_0__V_0)
  # 
  # # AP = asymptomatic prevented (asymptomatic under placebo, no inf vax)
  # # BP = both protected (symptomatic under placebo, no inf vax)
  # 
  # ## Upper bound: [P(Y_I = 0, Y_S = 0 | V = 1) - P(Y_I = 0, Y_S = 0 | V = 0)] / P(Y_I = 1, Y_S = 0 | V = 0)
  # P_Inf_0_Symp_0__V_1 <- length(which(sub_V_1[[event_name]] == 0)) / nrow(sub_V_1)
  # P_Inf_0_Symp_0__V_0 <- length(which(sub_V_0[[event_name]] == 0)) / nrow(sub_V_0)
  # 
  # upper_bound <- (P_Inf_0_Symp_0__V_1 - P_Inf_0_Symp_0__V_0) / P_Inf_1_Symp_0__V_0
  # 
  # return(list(lower_bound = lower_bound,
  #             upper_bound = upper_bound))
}