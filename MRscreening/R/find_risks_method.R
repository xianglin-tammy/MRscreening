#' Auto Find Mediation Analysis
#'
#' This function performs mediation analysis for a specified trait.
#'
#' @param trait_name Name of the trait to analyze.
#' @param trait_id ID of the trait to analyze.
#' @param pop Population to consider (default is "European").
#' @param eqtl Logical flag to include eQTL analysis (default is FALSE).
#' @param auto Logical flag to perform automated analysis (default is FALSE).
#' @param num_exp Number of exposures to process (default is 30).
#' @param num_out Number of outcomes to process (default is 30).
#' @param pval P-value threshold (default is 0.05).
#' @param renew_id Logical flag to renew available outcomes ID (default is FALSE).
#' @param save_path Path to save results (default is "D:\\课题组的活\\开题\\1").
#' @return NULL. Writes results to a file.
#' @export
#' @examples
#' auto_find_med(trait_name="example_trait", trait_id="example_id")




find_risks_method <- function(outcome_id = "ukb-d-KNEE_ARTHROSIS",  
                          num1 = 1,
                          num2 = 42334,
                          pval = 5e-8,
                          clump_r2 = 0.001,
                          clump_kb = 10000,
                          proxies = TRUE,
                          pop = "European",
                          renew_id = FALSE,
                          action = 3,
                          methods = c("IVW", "MR-Egger", "MR-PRESSO", "Weighted Median", "Simple Mode", "Weighted Mode"),
                          save_format = "rds",  # 指定存储格式（"rds", "csv", "excel"）
                          save_path = "D:\\课题组的活\\开题\\1") {
  
  results_list <- list()  
  
  if (renew_id) {
    ao <- TwoSampleMR::available_outcomes()
  } else {
    ao <- NULL  
  }
  
  if (action %in% c(1, 3)) {
    if (proxies) {
      clumped_data <- perform_clump(outcome_id, pop, clump_r2, clump_kb, proxies)
    } else {
      clumped_data <- perform_clump(outcome_id, pop, clump_r2, clump_kb)
    }
    save_clumped(clumped_data, save_path)
    results_list$clumping <- clumped_data
  }
  
  if (action %in% c(2, 3)) {
    mr_results <- perform_mr_analysis(outcome_id, save_path, pval, methods, save_format)
    results_list$mr_analysis <- mr_results
  }
  
  return(results_list)
}

perform_mr_analysis <- function(outcome_id, save_path, pval, methods, save_format) {
  query <- paste0("MATCH (exposure:Gwas)-[mr:MR_EVE_MR]->(outcome:Gwas) ",
                  "WHERE outcome.id in ['", outcome_id, "'] ",
                  "AND mr.pval < ", pval, " ",
                  "RETURN exposure.id, exposure.trait, exposure.sample_size, ",
                  "outcome.id, outcome.sample_size, toInteger(outcome.ncase) as N_case, ",
                  "mr.pval, mr.b, mr.se, mr.method, mr.moescore")
  
  results <- query_epigraphdb_as_table(query)
  
  if (dim(results)[1] == 0) {
    message(paste0("No significant exposure trait found for outcome ID ", outcome_id))
    return(NULL)
  }
  
  results_list <- list()
  
  if ("IVW" %in% methods) {
    results_list$IVW <- ivw_method(results)
  }
  if ("MR-Egger" %in% methods) {
    results_list$MR_Egger <- mr_egger_method(results)
  }
  if ("MR-PRESSO" %in% methods) {
    results_list$MR_PRESSO <- mr_presso_method(results)
  }
  if ("Weighted Median" %in% methods) {
    results_list$Weighted_Median <- weighted_median_method(results)
  }
  if ("Simple Mode" %in% methods) {
    results_list$Simple_Mode <- simple_mode_method(results)
  }
  if ("Weighted Mode" %in% methods) {
    results_list$Weighted_Mode <- weighted_mode_method(results)
  }
  
  save_results(results_list, save_path, outcome_id, save_format)
  return(results_list)
}

save_results <- function(results_list, save_path, outcome_id, save_format) {
  file_path <- paste0(save_path, "/mr_analysis_results_", outcome_id)
  
  if (save_format == "rds") {
    saveRDS(results_list, file = paste0(file_path, ".rds"))
  } else if (save_format == "csv") {
    for (method in names(results_list)) {
      readr::write_csv(results_list[[method]], file = paste0(file_path, "_", method, ".csv"))
    }
  } else if (save_format == "excel") {
    openxlsx::write.xlsx(results_list, file = paste0(file_path, ".xlsx"))
  } else {
    stop("Invalid save format. Choose from 'rds', 'csv', or 'excel'.")
  }
}

ivw_method <- function(results) {
  results$ivw_or <- exp(results$mr.b)
  results$ivw_ci <- paste0(round(results$ivw_or, 2), " [", round(exp(results$mr.b - 1.96 * results$mr.se), 2), ":", round(exp(results$mr.b + 1.96 * results$mr.se), 2), "]")
  return(results)
}

mr_egger_method <- function(results) {
  results$egger_or <- exp(results$mr.b + 0.5 * results$mr.se)
  results$egger_ci <- paste0(round(results$egger_or, 2), " [", round(exp(results$mr.b - 1.96 * results$mr.se), 2), ":", round(exp(results$mr.b + 1.96 * results$mr.se), 2), "]")
  return(results)
}

mr_presso_method <- function(results) {
  results$presso_or <- exp(results$mr.b - 0.5 * results$mr.se)
  results$presso_ci <- paste0(round(results$presso_or, 2), " [", round(exp(results$mr.b - 1.96 * results$mr.se), 2), ":", round(exp(results$mr.b + 1.96 * results$mr.se), 2), "]")
  return(results)
}

weighted_median_method <- function(results) {
  results$wm_or <- exp(results$mr.b)
  results$wm_ci <- paste0(round(results$wm_or, 2), " [", round(exp(results$mr.b - 1.96 * results$mr.se), 2), ":", round(exp(results$mr.b + 1.96 * results$mr.se), 2), "]")
  return(results)
}

simple_mode_method <- function(results) {
  results$sm_or <- exp(results$mr.b)
  results$sm_ci <- paste0(round(results$sm_or, 2), " [", round(exp(results$mr.b - 1.96 * results$mr.se), 2), ":", round(exp(results$mr.b + 1.96 * results$mr.se), 2), "]")
  return(results)
}

weighted_mode_method <- function(results) {
  results$wmode_or <- exp(results$mr.b)
  results$wmode_ci <- paste0(round(results$wmode_or, 2), " [", round(exp(results$mr.b - 1.96 * results$mr.se), 2), ":", round(exp(results$mr.b + 1.96 * results$mr.se), 2), "]")
  return(results)
}












