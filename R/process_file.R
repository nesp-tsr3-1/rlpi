#' Process an LPI dataset (dataset_name using the refernce year ref_year)
#'
#' @param dataset_name - The name of the dataset to process
#' @param ref_year - Reference year for the LPI - when the index == 1
#' @param model_selection_flag - Default=0
#' @param gam_global_flag - 1 = process by GAM method, 0 = process by chain method. Default=1
#' @param data_length_min - Minimum data length to include in calculations. Default=2
#' @param avg_time_between_pts_max - Maximum time between datapoint to include. Default=100
#' @param global_gam_flag_short_data_flag - Set this to 1 GAM model is also to be generated for the short time series else the log linear model will be used. Default=0
#' @param auto_diagnostic_flag - 1=Automatically determine whether GAM models are good enough, 0=Manually ask for each. Default=1
#' @param lambda_min - Minimum lambda to include in calculations. Default=1
#' @param lambda_max - Minimum lambda to include in calculations. Default=-1
#' @param zero_replace_flag - 0 = +minimum value; 1 = +1\% of mean value; 2 = +1. Default=1
#' @param offset_all - 1 = Add offset to all values, to avoid log(0). Default=0
#' @param offset_none=FALSE - Does nothing (leaves 0 unaffected **used for testing will break if there are 0 values in the source data **)
#' @param offset_diff=FALSE - Offset time-series with 0 values adding 1\% of mean if max value in time-series<1 and 1 if max>=1
#' @param linear_model_short_flag - If=TRUE models short time-series with linear model
#' @return - Return length of lamda array (number of lambda values?) - results are saved to file
#' @export
#'
process_file <- function(dataset_name,
                        ref_year,
                        model_selection_flag,
                        gam_global_flag,
                        data_length_min,
                        avg_time_between_pts_max,
                        global_gam_flag_short_data_flag,
                        auto_diagnostic_flag,
                        lambda_min,
                        lambda_max,
                        zero_replace_flag,
                        offset_all,
                        offset_none,
                        offset_diff,
                        linear_model_short_flag,
                        cap_lambdas,
                        show_progress,
                        basedir) {
  md5val <- tools::md5sum(dataset_name)
  # Read data file
  df <- read.table(dataset_name, header = TRUE, col.names = c("species", "id", "year", "popvalue"))

  final_year <- max(df$year)
  initial_year <- ref_year # Note that data could start later than ref_year

  cat("Calculating LPI for species\n")

  # Write header for population lambdas file
  # This file gets data appended to it as each dataset is processed, in calc_lpi.R
  pop_lambda_filename <- file.path(basedir, gsub(".txt", "_PopLambda.txt", dataset_name))
  Pop_Headers <- t(c("population_id", as.vector(initial_year:final_year)))
  write.table(Pop_Headers, file = pop_lambda_filename, sep = ",", eol = "\n", quote = FALSE, col.names = FALSE, row.names = FALSE)

  speciesLambda <- calc_lpi(
    species = df['species'],
    id = df['id'],
    year = df['year'],
    popvalue = df['popvalue'],
    initial_year = initial_year,
    final_year = final_year,
    dataset_name = dataset_name,
    model_selection_flag = model_selection_flag,
    gam_global_flag = gam_global_flag,
    data_length_min = data_length_min,
    avg_time_between_pts_max = avg_time_between_pts_max,
    global_gam_flag_short_data_flag = global_gam_flag_short_data_flag,
    auto_diagnostic_flag = auto_diagnostic_flag,
    lambda_min = lambda_min,
    lambda_max = lambda_max,
    zero_replace_flag = zero_replace_flag,
    offset_all = offset_all,
    offset_none = offset_none,
    offset_diff = offset_diff,
    linear_model_short_flag = linear_model_short_flag,
    cap_lambdas = cap_lambdas,
    show_progress = show_progress,
    basedir = basedir
  )

  # Save species Lambda matrix into a file

  data_file_name <- file.path(basedir, "lpi_temp", paste0(md5val, "_splambda.csv"))
  cat(sprintf("Saving species lambda to file: %s\n", data_file_name))

  write.table(speciesLambda, data_file_name, sep = ",", col.names = FALSE, row.names = FALSE)

  # Adding count of pops per species:
  # *****
  sp.count <- as.data.frame(table(df['species']))

  data_file_name <- file.path(basedir, gsub(".txt", "_lambda.csv", dataset_name))
  rownames(speciesLambda) <- t(unique(df['species']))
  colnames(speciesLambda) <- initial_year:final_year
  sorted_lambdas <- speciesLambda[order(rownames(speciesLambda)), ]
  sorted_lambdas_count <- cbind(sp.count, sorted_lambdas)
  colnames(sorted_lambdas_count)[1] = 'SpeciesSSet'

  cat(sprintf("Saving species lambda to file: %s\n", data_file_name))
  write.table(sorted_lambdas_count, data_file_name, sep = ",", col.names = NA)

  cat("Calculating DTemp\n")

  DTemp <- colMeans(speciesLambda, na.rm = TRUE)
  DTemp[is.nan(DTemp)] <- -99
  DTemp <- matrix(DTemp, nrow = 1, dimnames = list(1, names(DTemp)))

  # JW: we are just saving the same file twice in different places.
  # Is one for diagnostic purposes and the other caching?

  data_file_name <- file.path(basedir, "lpi_temp", paste0(md5val, "_dtemp.csv"))
  cat("Saving DTemp to file: ", data_file_name, "\n")
  write.table(DTemp, data_file_name, sep = ",", row.names = FALSE)

  data_file_name <- file.path(basedir, gsub(".txt", "_dtemp.csv", dataset_name))
  cat("Saving DTemp to file: ", data_file_name, "\n")
  write.table(DTemp, data_file_name, sep = ",", row.names = FALSE)

  return(list("species_lambda" = speciesLambda, "dtemp" = DTemp))
}
