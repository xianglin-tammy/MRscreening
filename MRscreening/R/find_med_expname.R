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

find_med_expname <- function(exposure_name,
                             pop = "European",
                             eqtl = FALSE,
                             num_exp = 1000,
                             num_out = 1000,
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

  exposure_data <- ao %>%
    dplyr::filter(population == pop)

  if (!is.data.frame(exposure_data)) {
    stop("exposure_data must be a data frame.")
  }

  eqtl_data <- exposure_data[grep("eqtl-*", exposure_data$id), ]

  if (!eqtl) {
    exposure_data <- exposure_data[!exposure_data$id %in% eqtl_data$id, ]
  }

  if (is.null(num_exp)) {
    num_exp <- nrow(exposure_data)
  }

  write_results <- function(out, file_path) {
    tryCatch({
      if (file.exists(file_path)) {
        utils::write.table(out, file = file_path, quote = FALSE, sep = ",", col.names = FALSE, append = TRUE, row.names = FALSE)
      } else {
        utils::write.table(out, file = file_path, quote = FALSE, sep = ",", col.names = TRUE, append = FALSE, row.names = FALSE)
      }
    }, error = function(e) {
      message("Error writing to file: ", e)
    })
  }

  exposure_info <- ao %>%
    dplyr::filter(trait == exposure_name & population == pop)

  if (nrow(exposure_info) == 0) {
    stop(paste("Exposure name", exposure_name, "not found in population", pop))
  }

  exposure_id <- exposure_info$id[1]  # Assume using the first match if multiple exist

  outcome_data <- ao %>%
    dplyr::filter(population == pop)

  if (!is.data.frame(outcome_data)) {
    stop("outcome_data must be a data frame.")
  }

  outcome_data <- outcome_data[!outcome_data$id %in% eqtl_data$id, ]

  if (is.null(num_out)) {
    num_out <- nrow(outcome_data)
  }

  for (i in 1:num_out) {
    outcome_id <- outcome_data$id[i]
    outcome_name <- outcome_data$trait[i]

    print(outcome_name)
    print(paste0(i, " :", outcome_id))

    out <- query_and_tidy_conf(exposure = exposure_id,
                               outcome = outcome_id,
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
}
