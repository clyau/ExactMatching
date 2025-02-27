## ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ##
## misc functions used for the simulation
## (1) num2cat(): convert numerical variables to categorical variables 
## (2) summ_results(): summarizes the simulation results
## (3) get_max(): find the max means from the list of mean dataframes
## (4) loop_wt_sd(): loop thru different variables with different weights
## ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ##
num2cat <- function(df, vars_to_cat, thresholds){
  ## thresholds needs to be a list of thresholds
  if(!is.list(thresholds)) 
    thresholds <- list(thresholds)
  ##
  for (i in (1 : length(vars_to_cat)) ) {
    df[[vars_to_cat[i]]] <- 
      cut(df[[vars_to_cat[i]]], 
          breaks = c(-Inf, thresholds[[i]], Inf),
          labels = LETTERS[1 : (length(thresholds[[i]]) + 1)])
  }
  return(df)
}
## ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ##
## :::::::::::::::::: summarize all simulation results ::::::::::::::::::::: ##
## simulations are done in script simAllDataRun.R, and the below function    ##
## ... is needed to summarize the results                                    ##
## ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ##
# summ_results <- function(results, seed_ipd1, seed_ipd2, cat_fn) {
#   ##                   ##
#   ## :::::: ess :::::: ##
#   ##                   ##
#   indices <- seq(3, length(results), by = 3)
#   selected_elements <- results[indices]
#   # Combine the selected elements into a single dataframe
#   ess_df <- do.call(rbind, selected_elements)
#   ess_df$randSeed.ipd1 <- seed_ipd1
#   ess_df$randSeed.ipd2 <- seed_ipd2
#   ##                              ##
#   ## :::::: weighted means :::::: ##
#   ##                              ##
#   indices <- seq(2, length(results), by = 3)
#   mean_list <- results[indices]
#   ##
#   if(cat_fn == 'cat201_minus1'){
#     derive_new_coln <- function(df, 
#                                 X3.B = "X3.B", X3.C = "X3.C", X3.D = "X3.D", 
#                                 X8.B = "X8.B", X8.C = "X8.C", 
#                                 X14.B = "X14.B", X14.C = "X14.C", 
#                                 X14.D = "X14.D", X14.E = "X14.E",
#                                 ## new colns
#                                 X3.Ad = "X3.Ad", X8.Ad = "X8.Ad", X14.Ad = "X14.Ad"
#     ) {
#       tdf <- data.frame(t(df)) 
#       tdf[[X3.Ad]] <- 1 - tdf[[X3.B]] - tdf[[X3.C]] - tdf[[X3.D]]
#       tdf[[X8.Ad]] <- 1 - tdf[[X8.B]] - tdf[[X8.C]]
#       tdf[[X14.Ad]] <- 1 - tdf[[X14.B]] - tdf[[X14.C]] - tdf[[X14.D]] - tdf[[X14.E]]
#       new.df <- data.frame(t(tdf))
#       return(new.df)
#     }
#     mean_list <- lapply(mean_list, derive_new_coln)
#   }
#   check_within_bounds <- function(df) {
#     is_within_bounds <- function(value, obs1, obs2){
#       value <- trunc(value * 10^9) / 10^9
#       obs1 <- trunc(obs1 * 10^9) / 10^9
#       osb2 <- trunc(obs2 * 10^9) / 10^9
#       return(value >= min(obs1, obs2) &
#                value <= max(obs1, obs2))
#     }
#     df$unc.in <- mapply(is_within_bounds, 
#                         value = df$ipd12.unc,
#                         obs1 = df$ipd1.obs, obs2 = df$ipd2.obs)
#     df$con.in <- mapply(is_within_bounds, 
#                         value = df$ipd12.con,
#                         obs1 = df$ipd1.obs, obs2 = df$ipd2.obs)
#     return(df)
#   }
#   mean_list <- lapply(mean_list, check_within_bounds)
#   
#   ##                               ##
#   ## :::::: simulated ipd's :::::: ##
#   ##                               ##
#   indices <- seq(1, length(results), by = 3)
#   ipd_list <- results[indices] ## need for max weights
#   ##                           ##
#   ## :::::: max weights :::::: ##
#   ##                           ##
#   colns_with_weights <- c("wt.unc.sdz", "wt.con.sdz", "wt.ps.sdz")
#   ##
#   # Function to get the max of each column in a dataframe
#   get_max <- function(df, cols) {
#     sapply(df[cols], function(column) {
#       if (all(is.na(column))) {
#         return(NA)  # or you could return some other value like NaN or 0
#       } else {
#         return(max(column, na.rm = TRUE))
#       }
#     })
#   }
#   ##
#   # Apply the function to each dataframe in the list
#   max_result <- lapply(ipd_list, get_max, cols = colns_with_weights)
#   max_result_df <- do.call(rbind, max_result)
#   ##
#   return(list(ess_df = ess_df,
#               max_df = max_result_df,
#               ipd_list = ipd_list,
#               mean_list = mean_list))
# }
## ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ##
## also part of the summ function
# Function to get the max of each column in a dataframe
get_max <- function(df, cols) {
  sapply(df[cols], function(column) {
    if (all(is.na(column))) {
      return(NA)  # or you could return some other value like NaN or 0
    } else {
      return(max(column, na.rm = TRUE))
    }
  })
}
## ::::::::::::::::::::::::::::::::::::::::::::::::::::::::
loop_wt_sd <- function(df, variables, weights) {
  # Initialize an empty list to store results
  wtd_sd <- list()
  # Loop over variables and weights
  for (variable in variables) {
    for (weight_name in names(weights)) {
      weight <- weights[[weight_name]]
      # Compute weighted variance for the current combination of variable and weight
      wtd_var <- Hmisc::wtd.var(df[[variable]], 
                                weight)
      # Store the result in the list
      wtd_sd[[paste0(paste(variable, 
                                weight_name, 
                                sep = "_"),  
                          '.sd')]] <- sqrt(wtd_var)
    }
  }
  # Convert the result list to a data frame for better readability
  wtd_sd_df <- data.frame(variable_weight = names(wtd_sd), 
                            wtd_sd = unlist(wtd_sd))
  wtd_sd_df <- separate(wtd_sd_df, 
                          col = variable_weight, 
                          into = c("vars", "wts"), 
                          sep = "_") 
  wtd_sd_df_w <- wtd_sd_df %>%
    pivot_wider(names_from = wts,
                values_from = wtd_sd)
  return(wtd_sd_df_w)
}
