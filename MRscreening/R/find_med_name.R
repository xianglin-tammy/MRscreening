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

find_med_name <- function(exposure_name,
                          outcome_name,
                          pop= "European",
                          eqtl= FALSE,
                          auto= FALSE,
                          num_exp=NULL,
                          num_out=NULL,
                          pval = 0.05,
                          renew_id = FALSE,
                          save_path = "D:\\课题组的活\\开题\\1") {

  if (renew_id) {
    ao <- TwoSampleMR::available_outcomes()
  } else {
    # Assuming ao is pre-loaded
  }

  if (!exists("ao") || !is.data.frame(ao)) {
    stop("ao must be a data frame and must exist.")
  }

  exposure_id <- ao %>%
    dplyr::filter(population == pop)

  if (!is.data.frame(exposure_id)) {
    stop("exposure_id must be a data frame.")
  }

  eqtl_id <- exposure_id[grep("eqtl-*", exposure_id$id), ]

  if (!eqtl) {
    exposure_id <- exposure_id[!exposure_id$id %in% eqtl_id$id, ]
  }

  if (is.null(num_exp)) {
    num_exp <- nrow(exposure_id)
  }

  write_results <- function(out, file_path) {
    tryCatch({
      if (file.exists(file_path)) {
        utils::write.table(out, file = paste0(save_path,"/",outcome_name,"_mediation_all.csv"), quote = FALSE, sep = ",", col.names = FALSE, append = TRUE, row.names = FALSE)
      } else {
        utils::write.table(out, file = paste0(save_path,"/","_mediation_all.csv"), quote = FALSE, sep = ",", col.names = TRUE, append = FALSE, row.names = FALSE)
      }
    }, error = function(e) {
      message("Error writing to file: ", e)
    })
  }

  if (!is.null(exposure_name) & !is.null(outcome_name)) {
    exposure_info <- ao %>%
      dplyr::filter(trait == exposure_name, population == pop) %>%
      slice(1)  # Take the first row if there are multiple matches

    outcome_info <- ao %>%
      dplyr::filter(trait == outcome_name, population == pop) %>%
      slice(1)  # Take the first row if there are multiple matches

    if (nrow(exposure_info) == 0) {
      stop(paste("Exposure", exposure_name, "not found in population", pop))
    }

    if (nrow(outcome_info) == 0) {
      stop(paste("Outcome", outcome_name, "not found in population", pop))
    }

    out <- query_and_tidy_conf(exposure = exposure_info$id,
                               outcome = outcome_info$id,
                               outcome_trait = outcome_name,
                               pval_threshold = pval,
                               mediator = TRUE)

    if (nrow(out) > 0) {
      print(out)
      file_path <- paste0(save_path, "/", exposure_name, "_", outcome_name, "_mediation_all.csv")
      write_results(out, file_path)
    } else {
      message("No data returned for exposure: ", exposure_name, " and outcome: ", outcome_name)
    }

  }

  if (auto) {
    outcome_id <- ao %>%
      dplyr::filter(population == pop)

    outcome_id <- outcome_id[!outcome_id$id %in% eqtl_id$id, ]

    if (is.null(num_out)) {
      num_out <- nrow(outcome_id)
    }

    for (i in 1:num_out) {
      print(outcome_id$trait[i])
      print(paste0(i, " :", outcome_id$id[i]))
      for (j in 1:num_exp) {
        if (outcome_id$id[i] != exposure_id$id[j]) {
          print(paste0(j, " :", exposure_id$id[j]))
          out <- query_and_tidy_conf(exposure = exposure_id$id[j],
                                     outcome = outcome_id$id[i],
                                     outcome_trait = outcome_id$trait[i],
                                     pval_threshold = pval,
                                     mediator = TRUE)

          if (nrow(out) != 0) {
            print(out)
            file_path <- paste0(save_path, "/auto_", exposure_id$id[j], "_", outcome_id$id[i], "_mediation_all.csv")
            write_results(out, file_path)
          } else {
            message("No data returned for exposure: ", exposure_id$id[j], " and outcome: ", outcome_id$id[i])
          }
        }
      }
    }
  }
}
