#' Print the output of a \code{"ve_nai"} object.
#'
#' @param x A \code{"ve_nai"} object.
#' @param what type to print, default VE estimates, other option probability of belonging to strata
#' @param ... Other arguments (not used)
#'
#' @method print ve_nai
#' @export
print.ve_nai <- function(object, what = "ve", ...){
  
  if(what == "ve"){
    for(i in 1:length(object$ve_fit$t0)){
      boot_t0 <- object$boot_est[[paste0("t0_", object$ve_fit$t0[i])]]
      cat(paste0("------- Results for t0 = ", object$ve_fit$t0[i], " -------\n"))
      for(ve in c("ve_i", "ve_s", "ve_ai", "ve_nai")){
        cat(paste0(ve, ": ", round(object$ve_fit[[ve]][i], 4), 
              " (", round(boot_t0[[ve]]$lower,4), ", ", 
              round(boot_t0[[ve]]$upper,4), ") \n"))
      }
    }
  } else if(what == "ps"){
    for(i in 1:length(object$ve_fit$t0)){
      boot_t0 <- object$boot_est[[paste0("t0_", object$ve_fit$t0[i])]]
      cat(paste0("------- Results for t0 = ", object$ve_fit$t0[i], " -------\n"))
      for(p in c("p_immune", "p_doomed", "p_alwaysinf", "p_converted", "p_helped", "p_helpedplus")){
        cat(paste0(p, ": ", round(object$ve_fit[[p]][i], 4), 
                   " (", round(boot_t0[[p]]$lower,4), ", ", 
                   round(boot_t0[[p]]$upper,4), ") \n"))
      }
    }
  } else if(what == "ve_nai_bounds"){
    if(!is.null(object$ve_nai_bounds)){
      for(i in 1:length(object$ve_nai_bounds$t0)){
        boot_t0 <- object$boot_est[[paste0("ve_nai_bound_t0_", object$ve_nai_bounds$t0[i])]]
        cat(paste0("------- VE_nai bounds for t0 = ", object$ve_nai_bounds$t0[i], " -------\n"))
        cat(paste0("Lower bound: ", round(object$ve_nai_bounds$lower_bound[i],4),
                   " (", round(boot_t0$lower_bound$lower,4), ", ",
                   round(boot_t0$lower_bound$upper,4), ") \n"))
        cat(paste0("Upper bound: ", round(object$ve_nai_bounds$upper_bound[i],4),
                   " (", round(boot_t0$upper_bound$lower,4), ", ",
                   round(boot_t0$upper_bound$upper,4), ") \n"))
        
      }
    }
  }

}

#' Print the output of a \code{"ps_ve"} object.
#'
#' @param x A \code{"ps_ve"} object.
#' @param what type to print, default VE estimates, other option probability of belonging to strata
#'
#' @method print ps_ve
#' @export
print.ps_ve <- function(object, what = "ve"){
  if(what == "ve"){
    print(data.frame(
      object[c("t0", "ve_i", "ve_s", "ve_ai", "ve_nai")]
    ))
  }else if(what == "ps"){
    print(data.frame(
      object[c("t0", paste0("p_", c("immune", "doomed", "alwaysinf", "converted", "helped", "helpedplus")))]
    ))
  }
}