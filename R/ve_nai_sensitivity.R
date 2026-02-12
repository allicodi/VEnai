#' Function for VE NAI sensitivity analysis
#' 
#' @param tbd
#' @param delta_X_variable If non \code{NULL} then a separate value of delta is
#' applied for each level of delta_X_varaible. Currently only supports variables
#' with two unique values. If non \code{NULL}, then code will calculate ve_nai for
#' a two-way grid of all values of \code{delta}, where values of delta_X_variable
#' are assigned all unique combinations of the user-supplied delta.
#' @export
ve_nai_sensitivity <- function(
    symp_cox_fit,
    asymp_cox_fit,
    symp_lr_fit,
    data,
    t0 = quantile(data$ftime, c(0.2, 0.5, 0.9)),
    vax_name = "vax",
    symp_ind_name = "symp_ind",
    delta = 1,
    delta_X_variable = NULL,
    ...
){
  # P[Y_s(0) = 0 | Y_i(0) = 1, Y_i(1) = 1, X = x] /
  # P[Y_s(0) = 0 | Y_i(0) = 1, Y_i(1) = 0, X = x] = delta
  # for all x and for delta in (0, inf)
  
  # Need to estimate:
  # P(Y_S = 0 | V = 0, Y_I = 1, X = x) --> new? but just 1 - p_symp__inf_vax1? see ??? below
  # P(Y_I = 1 | V = 1, X = x) --> ci_inf_vax1 from original
  # P(Y_I = 1 | V = 0, X = x) --> ci_inf_vax0 from original
  if(!is.null(delta_X_variable)){
    delta_X_variable_levels <- unique(data[[delta_X_variable]])
    stopifnot(length(delta_X_variable_levels) == 2)
  }

  n <- dim(data)[1]

  symp_chaz <- survival::basehaz(symp_cox_fit)
  asymp_chaz <- survival::basehaz(asymp_cox_fit)
  # if the unique times are not the same, then it may indicate a problem with
  # how the models were fit
  stopifnot(all(unique(symp_chaz$time) == unique(asymp_chaz$time)))
  
  # find times closest to requested t0
  closest_indices <- sapply(t0, function(t) {
    which.min(abs(symp_chaz$time - t))
  })
  subset_symp_chaz <- symp_chaz[closest_indices,]
  subset_asymp_chaz <- asymp_chaz[closest_indices,]
  
  surv_data <- merge(subset_symp_chaz, subset_asymp_chaz, by = "time", suffixes = c("_symp", "_asymp"))
  
  surv_data$total_chaz <- surv_data$hazard_symp + surv_data$hazard_asymp
  # total baseline cumulative incidence of any infection
  surv_data$total_ci <- 1 - exp(-surv_data$total_chaz)
  
  # predict from models under vaccine/no vaccine
  data0 <- data; data0[[vax_name]] <- 0
  lp_symp_vax0 <- predict(symp_cox_fit, newdata = data0, type = "lp")
  lp_asymp_vax0 <- predict(asymp_cox_fit, newdata = data0, type = "lp")
  # P(Y_S = 1 | V = 0,  X = x)
  p_symp__inf_vax0 <- predict(
    symp_lr_fit, newdata = data0, type = "response"
  )
  
  data1 <- data; data1[[vax_name]] <- 1
  lp_symp_vax1 <- predict(symp_cox_fit, newdata = data1, type = "lp")
  lp_asymp_vax1 <- predict(asymp_cox_fit, newdata = data1, type = "lp")
  p_symp__inf_vax1 <- predict(
    symp_lr_fit, newdata = data1, type = "response"
  )
  
  ## NEW ##
  ## P(Y_S = 0 | V = 0, Y_I = 1, X = x) 
  
  # ??? symp_lr_fit was already fit in infected participants only? 
  # so P(Y_S = 1 | V = 1, X = X_i), i = 1,...,n = p_symp__inf_vax1 == P(Y_S = 1 | V = 1, Y_I = 1, X = X_i)??
  # in which case:
  p_symp0__inf_vax0 <- 1 - p_symp__inf_vax0
  
  # conditional cumulative incidence any infection under no vaccine/vaccine
  # P(Y_I = 1 | V = v, X = x)
  ci_inf_vax0 <- mapply(
    hazard_symp = surv_data$hazard_symp,
    hazard_asymp = surv_data$hazard_asymp,
    FUN = function(hazard_symp, hazard_asymp){
      1 - exp(- ( hazard_symp * exp(lp_symp_vax0) + hazard_asymp * exp(lp_asymp_vax0) ))
    }
  )
  ci_inf_vax1 <- mapply(
    hazard_symp = surv_data$hazard_symp,
    hazard_asymp = surv_data$hazard_asymp,
    FUN = function(hazard_symp, hazard_asymp){
      1 - exp(- ( hazard_symp * exp(lp_symp_vax1) + hazard_asymp * exp(lp_asymp_vax1) ))
    }
  )
  
  # Alternative VE_NAI:
  # VE_I(X_i) = 1 - 1 - ci_inf_vax1[i,] / ci_inf_vax0[i,]
  # Numerator: {p_symp0__inf_vax0 / (1 + VE_I(X_i) * [(1- delta) / delta]) } * ci_inf_vax1
  # Denominator: p_symp0__inf_vax0 * ci_inf_vax0
  
  # VE_I(X_i)
  
  # issue in some bootstrap samples where ci_inf_vax0 == 0, truncate temp?
  ci_inf_vax1 <- ifelse(ci_inf_vax1 < 0.001, 0.001, ci_inf_vax1)
  ci_inf_vax0 <- ifelse(ci_inf_vax0 < 0.001, 0.001, ci_inf_vax0)
  
  ve_any_inf <- 1 - ci_inf_vax1 / ci_inf_vax0
  
  # Denominator:
  ve_nai_denom <- colMeans(p_symp0__inf_vax0 * ci_inf_vax0)
  #ve_nai_denom <- colMeans(1 - ve_any_inf + (1/delta)*ve_any_inf)
  
  if(is.null(delta_X_variable)){
    delta_res <- vector(mode = "list", length = length(delta))

    for(d in 1:length(delta)){

      # Numerator:
      ve_nai_num <- colMeans((p_symp0__inf_vax0 / 
                                #(1 + (ve_any_inf * ((1 - delta[d]) / delta[d]))))
                                (1 - ve_any_inf + (1/delta[d])*ve_any_inf)
                             * ci_inf_vax1))
      
      ve_nai <- 1 - ve_nai_num / ve_nai_denom
      
      delta_res[[d]] <- data.frame(delta = delta[d],
                                   t0 = t0,
                                   ve_nai = ve_nai)
      
    }
  }else{
    delta_res <- vector(mode = "list", length = length(delta)^2)
    idx_level_1 <- which(data[[delta_X_variable]] == delta_X_variable_levels[1])
    idx_level_2 <- which(data[[delta_X_variable]] == delta_X_variable_levels[2])
    ct <- 0
    for(d1 in 1:length(delta)){
      for(d2 in 1:length(delta)){
        ct <- ct + 1
        delta_vec <- rep(NA, n)
        delta_vec[idx_level_1] <- delta[d1]
        delta_vec[idx_level_2] <- delta[d2]

        # Numerator:
        ve_nai_num <- colMeans((p_symp0__inf_vax0 / 
                                  (1 - ve_any_inf + (1/delta_vec)*ve_any_inf)
                               * ci_inf_vax1))
        
        ve_nai <- 1 - ve_nai_num / ve_nai_denom
        
        delta_res[[ct]] <- data.frame(
          delta1 = delta[d1],
          delta2 = delta[d2],
          t0 = t0,
          ve_nai = ve_nai
        )
      }
    }
    attr(delta_res, "delta_X_variable_levels") <- delta_X_variable_levels
  }
  
  delta_res <- do.call(rbind, delta_res)
  
  return(delta_res)
  
}


