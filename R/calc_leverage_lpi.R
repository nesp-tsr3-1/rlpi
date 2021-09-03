
#' calc_lambdas - calculate interannual differences from an index
#'
#' @param lpi - The index to calculate lambdas from
#'
#' @return Returns lamdbas for the lpi given
#' @export
#'
calc_lambdas <- function(lpi) {
  N <- length(lpi)
  lambdas <- matrix(0, N)
  lambdas[1] <- 1
  for (k in 2:N) {
    lambdas[k] <- log10(lpi[k]) - log10(lpi[k - 1])
  }
  return(lambdas)
}



#' #' Calculate the index using the annual differences in species_lambda_array
#' @param species_lambda_array Array of DTemps (annual differences)
#' @param fileindex The index of the file that this index is for
#' @param dsize The size of the data in DTemp
#' @param group Which group this file belongs to
#' @param weightings What the weightings are for this group
#' @param use_weightings Whether or not to use weightings (level 1)
#' @param use_weightings_b Whether or not to use weightings (level 1)
#' @param weightings_b What the weightingsB are for this group
#'
#' @return Returns an index
#' @export
#'
calc_leverage_lpi <- function(species_lambda_array, fileindex, dsize, group, weightings, use_weightings, use_weightings_b, weightings_b) {
  NoFiles <- length(unique(fileindex))
  Nogroups <- length(unique(group[[1]]))

  I <- matrix(0, dsize)
  I[1] <- 1

  cat(".")

  # For each year
  for (J in 2:dsize) {
    # Make two matrices of 0s of size 1xNogroups
    D <- matrix(0, 1, Nogroups)
    DI <- matrix(0, 1, Nogroups)

    # For each file (population) in this file/set of species lambdas
    for (FileNo in 1:NoFiles) {
      groupNo <- group[FileNo, 1]

      # Read speciesLambda from saved file FileName = paste('lpi_temp/speciesLambda',FileNo,sep='')
      # speciesLambda = read.table(FileName, header = FALSE, sep=',')

      speciesLambda <- species_lambda_array[fileindex == FileNo, ]

      if (J <= dim(speciesLambda)[2]) {

        # If there's still some left to get
        # Get the lamdas for this pop (species?)
        speciesLambdaVal <- speciesLambda[, J]

        # We shouldn't be sampling missing values....
        speciesLambdaVal <- speciesLambdaVal[!is.na(speciesLambdaVal)]
        # If we've got some meaningful data
        Index <- which(speciesLambdaVal != -1)
        if (length(Index) > 0) {
          if (use_weightings) {
            D[groupNo] <- D[groupNo] + mean(speciesLambdaVal[Index]) * weightings[[1]][FileNo]
          } else {
            D[groupNo] <- D[groupNo] + mean(speciesLambdaVal[Index])
          }

          DI[groupNo] <- DI[groupNo] + 1
        }
      }
    }

    # For each D
    # Take average if there's values (otherwise 0) - so this gives group average
    for (DIndex in 1:length(D)) {
      if (use_weightings == 1) {
        if (DI[DIndex] > 1) DI[DIndex] <- 1
      }
      if (DI[DIndex] > 0) {
        D[DIndex] <- D[DIndex] / DI[DIndex]
      } else {
        D[DIndex] <- 0
      }
    }

    DT <- 0
    DI <- 0
    # Sum over groups
    for (groupNo in 1:Nogroups) {
      # CHANGED AS D CAN BE 0 I.E ZERO GROWTH (AVERAGED)
      # if (D[groupNo] != 0) {
      if (use_weightings_b == 1) {
        # Catch any groups which have no data.
        if (!is.na(D[groupNo])) {
          DT <- DT + D[groupNo] * weightings_b[groupNo]
          # DI = DI + 1
        }
        DI <- 1
      } else {
        # Catch any groups which have no data.
        if (!is.na(D[groupNo])) {
          DT <- DT + D[groupNo]
          DI <- DI + 1
        }
      }
      # }
    }

    # Return the bootstrapped index
    if (DI == 0) {
      # If there was no data in this run, set to -1
      I[J] <- -1
      I[J] <- NA
    } else {
      if (is.na(I[J - 1])) {
        I[J] <- 1 * 10^(DT / DI)
      } else {
        I[J] <- I[J - 1] * 10^(DT / DI)
      }
    }
  }
  return(I)
}
