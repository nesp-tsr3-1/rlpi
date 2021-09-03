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
ProcessFile <- function(dataset_name,
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
  Data <- read.table(dataset_name, header = TRUE)

  # Get data from file as column vectors
  SpeciesSSet <- Data[1]
  idSSet <- Data[2]
  yearSSet <- Data[3]
  popvalueSSet <- Data[4]

  # Forget 'data' variable
  rm(Data)

  final_year <- max(yearSSet)

  # Note that data could be later than ref_year
  if (min(yearSSet) < ref_year) {
    initial_year <- ref_year
  } else {
    initial_year <- min(yearSSet)
  }
  initial_year <- ref_year

  # Call the LPI function sNames = unique(speciesSSet); nospecies = max(dim(sNames))
  # MethodFlag = matrix(0,1,nospecies)

  cat("Calculating LPI for species\n")

  pop_lambda_filename <- file.path(basedir, gsub(".txt", "_PopLambda.txt", dataset_name))
  Pop_Headers <- t(c("population_id", as.vector(initial_year:final_year)))
  write.table(Pop_Headers, file = pop_lambda_filename, sep = ",", eol = "\n", quote = FALSE, col.names = FALSE, row.names = FALSE)

  speciesLambda <- CalcLPI(
    species = SpeciesSSet,
    id = idSSet,
    year = yearSSet,
    popvalue = popvalueSSet,
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

  # FileNo removed - could use name?! *** #DataFileName = paste("lpi_temp/speciesLambda", FileNo, sep = "")
  # cat(sprintf("Saving species lambda to file: %s\n", DataFileName))
  # write.table(speciesLambda, DataFileName, sep = ",", col.names = FALSE, row.names = FALSE)

  # cat(dim(speciesLambda), "\n")
  # cat(dim(unique(speciesSSet)), "\n")

  DataFileName <- file.path(basedir, "lpi_temp", paste0(md5val, "_splambda.csv"))
  cat(sprintf("Saving species lambda to file: %s\n", DataFileName))
  write.table(speciesLambda, DataFileName, sep = ",", col.names = FALSE, row.names = FALSE)

  # Adding count of pops per species:
  # *****
  sp.count <- as.data.frame(table(SpeciesSSet))

  DataFileName <- file.path(basedir, gsub(".txt", "_lambda.csv", dataset_name))
  rownames(speciesLambda) <- t(unique(SpeciesSSet))
  colnames(speciesLambda) <- initial_year:final_year
  sorted_lambdas <- speciesLambda[order(rownames(speciesLambda)), ]
  sorted_lambdas_count <- cbind(sp.count, sorted_lambdas)

  cat(sprintf("Saving species lambda to file: %s\n", DataFileName))
  write.table(sorted_lambdas_count, DataFileName, sep = ",", col.names = NA)

  # rm(speciesSSet, idSSet, yearSSet, popvalueSSet)

  cat("Calculating DTemp\n")
  DTemp <- matrix(0, 1, dim(speciesLambda)[2])
  # For each year
  for (I in 1:dim(speciesLambda)[2]) {
    # Get data for this year 'I'
    yearData <- speciesLambda[, I]

    # Find populations that have data
    if (!cap_lambdas) {
      Index <- which(yearData != -1)
    } else {
      Index <- which(!is.na(yearData))
    }

    # If there are some populations
    if (length(Index) > 0) {
      # DTemp is mean lambda for those populations
      DTemp[I] <- mean(yearData[Index])
      # Otherwise -99
    } else {
      DTemp[I] <- -99
    }
  }

  # RF: Each file returns dimensions now
  # if (dim(speciesLambda)[2] > dsize)
  #  dsize = dim(speciesLambda)[2]

  # Save DTemp into file
  # FileNo removed - could use name?! *** #DataFileName = paste("lpi_temp/DTemp", FileNo, sep = "")
  # write.table(DTemp, DataFileName, sep = ",", col.names = FALSE, row.names = FALSE)

  DataFileName <- file.path(basedir, "lpi_temp", paste0(md5val, "_dtemp.csv"))

  cat("Saving DTemp to file: ", DataFileName, "\n")
  # write.table(DTemp, DataFileName, sep = ",", col.names = FALSE, row.names = FALSE)

  colnames(DTemp) <- initial_year:final_year

  write.table(DTemp, DataFileName, sep = ",", row.names = FALSE)

  DataFileName <- file.path(basedir, gsub(".txt", "_dtemp.csv", dataset_name))
  cat("Saving DTemp to file: ", DataFileName, "\n")
  # write.table(DTemp, DataFileName, sep = ",", col.names = FALSE, row.names = FALSE)
  write.table(DTemp, DataFileName, sep = ",", row.names = FALSE)

  # Return length of lamda array (number of lamda values?)
  # cat("Returning length of lambda\n")
  return(dim(speciesLambda)[2])
}
