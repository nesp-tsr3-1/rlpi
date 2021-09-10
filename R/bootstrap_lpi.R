
#' bootstrap_lpi
#'
#' @param species_lambda_array - Array of DTemps (annual differences)
#' @param fileindex The index of the file that this index is for
#' @param dsize The size of the data in DTemp
#' @param group Which group this file belongs to
#' @param weightings What the weightings are for this group
#' @param use_weightings Whether or not to use weightings (level 1)
#' @param use_weightings_b Whether or not to use weightings (level 1)
#' @param weightings_b What the weightingsB are for this group
#'
#' @return Returns a bootstrapped LPI
#' @export
#'
bootstrap_lpi <- function(
  species_lambda_array,
  fileindex,
  dsize,
  group,
  weightings,
  use_weightings,
  use_weightings_b,
  weightings_b,
  cap_lambdas) {
  num_groups <- length(unique(group[[1]]))

  d <- data.frame(species_lambda_array, fileindex) %>%
    mutate(
      group = group[fileindex, 1],
      i = row_number()
    ) %>%
    tidyr::pivot_longer(cols = starts_with("V"), names_to="year") %>%
    mutate(year = as.integer(substring(year, 2))) %>%
    filter(year > 1 & !is.na(value)) %>%
    group_by(year, fileindex, group) %>%
    summarise(d = mean(sample(value, replace = TRUE))) %>%
    group_by(year, group) %>%
    summarise(d = ifelse(use_weightings, sum(d * weightings[fileindex, 1]), mean(d))) %>%
    group_by(year) %>%
    # ISSUE: maybe sum(d)/num_groups should be mean(d)
    summarise(d = ifelse(use_weightings_b, sum(d * weightings_b[group, 1]), sum(d) / num_groups)) %>%
    .$d

  d <- purrr::accumulate(d, function(a, x) { ifelse(is.na(a), 1, a) * 10^x }, .init = 1)

  return(d)
}

bootstrap_lpi_old <- function(
  species_lambda_array,
  fileindex,
  dsize,
  group,
  weightings,
  use_weightings,
  use_weightings_b,
  weightings_b,
  cap_lambdas) {
  no_files <- length(unique(fileindex))
  no_groups <- length(unique(group[[1]]))

  # Initialise first boot_i for this loop to 1
  boot_i <- matrix(0, dsize)
  boot_i[1] <- 1

  cat(".")

  # For each year
  for (J in 2:dsize) {
    # Make two matrices of 0s of size 1xno_groups
    D <- matrix(0, 1, no_groups)
    DI <- matrix(0, 1, no_groups)

    # For each file (population) in this file/set of species lambdas
    for (file_no in 1:no_files) {
      group_no <- group[file_no, 1]

      # Read species_lambda from saved file FileName = paste('lpi_temp/species_lambda',file_no,sep='')
      # species_lambda = read.table(FileName, header = FALSE, sep=',')
      species_lambda <- species_lambda_array[fileindex == file_no, J]

      if (!is.null(species_lambda)) {
        # We shouldn't be sampling missing values....
        species_lambda_val <- na.omit(species_lambda)

        # Create sample with replacement (single bootstrap instance)
        boot_val <- sample(species_lambda_val, replace = T)

        # If we've got some meaningful data
        if (!cap_lambdas) {
          index <- which(boot_val != -1)
        } else {
          index <- which(!is.na(boot_val))
        }
        if (length(index) > 0) {

          # Store sum of mean lamdas in D (summing over species within group)
          # D[group[file_no, 1]] = D[group[file_no, 1]] + mean(boot_val[index])

          if (use_weightings) {
            D[group_no] <- D[group_no] + mean(boot_val[index]) * weightings[[1]][file_no]
          } else {
            D[group_no] <- D[group_no] + mean(boot_val[index])
          }

          DI[group_no] <- DI[group_no] + 1
        }
      }
    }

    # For each D
    # Take average if there's values (otherwise 0) - so this gives group average
    for (Dindex in 1:length(D)) {
      if (use_weightings == 1) {
        if (DI[Dindex] > 1) DI[Dindex] <- 1
      }
      if (DI[Dindex] > 0) {
        D[Dindex] <- D[Dindex] / DI[Dindex]
      } else {
        D[Dindex] <- 0
      }
    }

    DT <- 0
    DI <- 0
    # Sum over groups
    for (group_no in 1:no_groups) {
      # CHANGED AS D CAN BE 0 I.E ZERO GROWTH (AVERAGED)
      # if (D[group_no] != 0) {
      if (use_weightings_b == 1) {
        # Catch any groups which have no data.
        if (!is.na(D[group_no])) {
          DT <- DT + D[group_no] * weightings_b[group_no]
          # DI = DI + 1
        }
        DI <- 1
      } else {
        # Catch any groups which have no data.
        if (!is.na(D[group_no])) {
          DT <- DT + D[group_no]
          DI <- DI + 1
        }
      }
      # }
    }

    # if (use_weightings == 1) {
    #  if (DI > 1) DI = 1
    # }
    # Return the bootstrapped index
    if (DI == 0) {
      # If there was no data in this run, set to -1
      boot_i[J] <- -1
      boot_i[J] <- NA
    } else {
      if (is.na(boot_i[J - 1])) {
        boot_i[J] <- 1 * 10^(DT / DI)
      } else {
        boot_i[J] <- boot_i[J - 1] * 10^(DT / DI)
      }
      # boot_i[J] = boot_i[J - 1] * 10^(DT/DI)
    }
  }
  # cat("boot_i: ", boot_i, "\n")
  return(boot_i)
}
