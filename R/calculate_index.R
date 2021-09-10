#' Calculate the index using the annual differences in dtemp_array
#'
#' @param dtemp_array Array of DTemps (annual differences)
#' @param fileindex The index of the file that this index is for
#' @param dsize The size of the data in DTemp
#' @param group Which group this file belongs to
#' @param weightings What the weightings are for this group
#' @param use_weightings Whether or not to use weightings (level 1)
#' @param use_weightings_b Whether or not to use weightings (level 1)
#' @param weightings_b What the weightingsB are for this group
#'
#' @return index - the calculated index
#' @export
#'
calculate_index <- function(dtemp_array, fileindex, dsize, group, weightings, use_weightings, use_weightings_b, weightings_b) {
  d <- data.frame(dtemp_array, fileindex = unique(fileindex)) %>%
    mutate(group = group[fileindex, 1]) %>%
    tidyr::pivot_longer(cols = starts_with("X"), names_to="year", values_to="d") %>%
    mutate(year = as.integer(substring(year, 2))) %>%
    filter(year > 1 & !is.na(d) & d != -99) %>%
    group_by(year, group) %>%
    summarise(d = ifelse(use_weightings, sum(d * weightings[fileindex, 1]), mean(d))) %>%
    group_by(year) %>%
    summarise(d = ifelse(use_weightings_b, sum(d * weightings_b[group, 1]), mean(d)))

  d2 <- rep(NA, dsize - 1)
  d2[d$year - 1] <- d$d

  d2 <- d2 %>%
    purrr::accumulate(function(a, x) { ifelse(is.na(a), 1, a) * 10^x }, .init = 1) %>%
    replace(is.na(.), -99)

  return(d2)
}

calculate_index_old <- function(dtemp_array, fileindex, dsize, group, weightings, use_weightings, use_weightings_b, weightings_b) {
#   # calculate LPI

  NoFiles <- length(unique(fileindex))
  Nogroups <- length(unique(group[[1]]))

  I <- matrix(0, dsize)
  I[1] <- 1

  # For each year indexed as 'J'
  for (J in 2:dsize) {
    # Create two vectors, one for summing 'lambdas' one for counting number summed
    D <- matrix(0, 1, Nogroups)
    DI <- matrix(0, 1, Nogroups)

    # cat("DI: ", DI, "\n")
    # cat("Nogroups: ", Nogroups, "\n")
    # For each file
    for (FileNo in 1:NoFiles) {
      groupNo <- group[FileNo, 1]
      # cat("groupNo: ", groupNo, "\n")
      # cat("group: \n")
      # print(group)
      # Read speciesLambda and DTemp from saved file
      # speciesLambdas are the annual differences in population for each species (row for each sp)
      # speciesLambda = species_lambda_array[fileindex == FileNo, ]

      # DTemps are the mean annual differences in population for each group/file
      # DTemp for this group/file
      DTemp <- as.matrix(dtemp_array[FileNo, ])

      # cat("DTemp: ", DTemp, "\n")

      if (J <= dim(DTemp)[2]) {
        # If it's not a missing value
        if (!is.na(DTemp[J])) {
          # or an earlier flag value
          if (DTemp[J] != -99) {
            if (use_weightings == 1) {
              # cat(sprintf("Using weighting %f for file number %d\n", weightings[[1]][FileNo], FileNo))
              D[groupNo] <- D[groupNo] + DTemp[J] * weightings[[1]][FileNo]
            } else {
              D[groupNo] <- D[groupNo] + DTemp[J]
            }
            DI[groupNo] <- DI[groupNo] + 1
          }
        }
      }
    }

    # cat("\n-D: \n", D, "\n")
    # cat("\n-DI: \n", DI, "\n")

    for (DIndex in 1:length(D)) {
      if (use_weightings == 1) {
        if (DI[DIndex] > 1) DI[DIndex] <- 1
      }
      # If more than one file contributed, take the average
      if (DI[DIndex] > 0) {
        D[DIndex] <- D[DIndex] / DI[DIndex]
      } else {
        # Otherwise, if no files contributed, set to 0??
        # D[DIndex] = 0
        # Surely better to flag missing value
        D[DIndex] <- NA
      }
    }

    # Average over groups
    DT <- 0
    DI <- 0
    # cat("DT: \n", DT, "\n")
    # cat("DI: \n", DI, "\n")
    # cat("weightings_b: \n", weightings_b, "\n")
    for (groupNo in seq(1, Nogroups)) {
      # CHANGED AS D CAN BE 0 I.E ZERO GROWTH (AVERAGED)
      # if (D[groupNo] != 0) {
      if (use_weightings_b == 1) {
        # cat(weightings_b, "\n")
        if (!is.na(D[groupNo])) {
          DT <- DT + D[groupNo] * weightings_b[groupNo]
          # DI = DI + 1
          # cat("A: ", DT, "\n")
        } else {
          cat(sprintf("group %d is NA in year %d\n", groupNo, J))
        }
        DI <- 1
      } else {
        if (!is.na(D[groupNo])) {
          DT <- DT + D[groupNo]
          DI <- DI + 1
          # cat("B: ", DT, "\n")
        } else {
          cat(sprintf("group %d is NA in year %d\n", groupNo, J))
        }
      }
      # }
    }

    # if (use_weightings == 1) {
    #  if (DI > 1) DI = 1
    # }

    # cat("DT: \n", DT, "\n")
    # cat("DI: \n", DI, "\n")

    # cat("I: ", I, "\n")

    if (is.na(DT)) {
      I[J] <- NA
      cat(sprintf("year %d is NA\n", J))
    } else {
      if (DI > 0) {
        if (is.na(I[J - 1]) | (I[J - 1] == -99)) {
          I[J] <- 1 * 10^(DT / DI)
          cat(sprintf("**** [year %d] Previous year data missing, assuming '1' **** \n", J))
        } else {
          I[J] <- I[J - 1] * 10^(DT / DI)
        }
      } else {
        I[J] <- -99
      }
    }
  }
  return(I)
}
