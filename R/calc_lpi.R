#' CalcLPI - Main funciton for calculating species lamdbas (interannual changes). Input data is a series
#' of 4-value rows: (species, id, year, popvalue). This function will model the species populations over time
#' using either the chain method (log-linear interpolation) or Generalised Additive Modelling. See gam_global_flag
#'
#'
#' @param species - Vector with name of species for each population value
#' @param id - Vector of ids for each population value
#' @param year - Vector of years for each population value
#' @param popvalue - Vector of population values
#' @param initial_year - Initial year to calculate the index from
#' @param final_year - Final year to calculate the index to
#' @param dataset_name - Name of the dataset that these value are from (for generating output files)
#' @param model_selection_flag - Default=0
#' @param gam_global_flag  - 1 = process by GAM method, 0 = process by chain method. Default=1
#' @param data_length_min - Minimum data length to include in calculations. Default=2
#' @param avg_time_between_pts_max - Maximum time between datapoint to include. Default=100
#' @param global_gam_flag_short_data_flag - Set this to 1 GAM model is also to be generated for the short time series else the log linear model will be used. Default=0
#' @param auto_diagnostic_flag - 1=Automatically determine whether GAM models are good enough, 0=Manually ask for each. Default=1
#' @param lambda_min - Minimum lambda to include in calculations. Default=1
#' @param lambda_max - Minimum lambda to include in calculations. Default=-1
#' @param zero_replace_flag - 0 = +minimum value; 1 = +1\% of mean value; 2 = +1. Default=2 Only for time-series that contain 0 values
#' @param offset_all - 1 = Add offset to all values in all time-series, to avoid log(0). Default=0
#' @param offset_none=FALSE - Does nothing (leaves 0 unaffected **used for testing will break if there are 0 values in the source data **)
#' @param offset_diff=FALSE - Offset time-series with 0 values adding 1\% of mean if max value in time-series<1 and 1 if max>=1
#' @param linear_model_short_flag - If=TRUE models short time-series with linear model
#' @return - Returns the species lambda array for the input species
#' @export
#'
CalcLPI <- function(species,
                    id,
                    year,
                    popvalue,
                    initial_year,
                    final_year,
                    dataset_name,
                    model_selection_flag, # determines whether we approve the models or the code does it automatically
                    gam_global_flag, # 1 = process by GAM method, 0 = process by chain method
                    data_length_min,
                    avg_time_between_pts_max,
                    global_gam_flag_short_data_flag, # set this if GAM model is also to be generated for the short time series else the log linear model will be used.
                    auto_diagnostic_flag, # is about how the smoothing parameter is chosen
                    lambda_min,
                    lambda_max,
                    zero_replace_flag, # 0 = +minimum value; 1 = +1% of mean value; 2 = +1
                    offset_all, # Add offset to all values, to avoid log(0)
                    offset_none, # Does nothing (leaves 0 unaffected **used for testing will break if there are 0 values in the source data **)
                    offset_diff, # Offset time-series with 0 values adding 1% of mean if max value in time-series<1 and 1 if max>=1
                    linear_model_short_flag, # if=TRUE models short time-series with linear model
                    cap_lambdas = FALSE,
                    show_progress = FALSE,
                    basedir = ".") {
  noRecs <- max(dim(popvalue))
  sNames <- unique(species)
  sid <- unique(id)
  nospecies <- max(dim(sNames))
  noPop <- max(dim(unique(id)))

  PopNotProcessed <- matrix(0, 1, noPop)
  MethodFlag <- matrix(0, 1, noPop)

  PopNotProcessedCounter <- 0
  PopProcessedGAMCounter <- 0
  sNamesCounter <- 0
  sNamesArray <- sNames
  sidArray <- sid
  PopProcessedGAM <- matrix(0, 1, noPop)

  cat(sprintf("Number of species: %s (in %s populations)\n", nospecies, noPop))
  # Calculate index value for each species
  speciesLambda <- matrix(0, nospecies, final_year - initial_year + 1)

  if (show_progress) {
    prog <- txtProgressBar(min = 0, max = nospecies, char = "*", style = 3)
  }
  # Here a file is created that includes the list of populations processed by chain
  MethodFlagLoop <- 0
  for (I in 1:nospecies) {
    # cat(".")
    Index <- 1

    # Delete sIndex if it exists
    if (length(which(objects() == "sIndex")) != 0) {
      rm(sIndex)
    }

    # Get indices of this species in 'species'
    sIndex <- which(species == toString(sNames[I, 1]))

    # Extract population data using that index
    Popid <- unique(id[sIndex, 1])

    # Delete PopLambda if it exists
    if (length(which(objects() == "PopLambda")) != 0) {
      rm(PopLambda)
    }

    PopIdsize <- length(Popid)
    # Blank matrix of -1s - default value of poplambd
    PopLambda <- matrix(-1, PopIdsize, final_year - initial_year + 1)
    JIndex <- 1
    # For each population of this species
    for (J in 1:PopIdsize) {
      IndexPop <- which(id == Popid[J]) # each population
      yearPop <- year[IndexPop, 1]
      PopN <- popvalue[IndexPop, 1]

      # Perform the data filtering (exclude time-series with two data points and more than 100 years between two points)
      DataTimeLength <- max(yearPop) - min(yearPop)
      AvgTimeBetweenPts <- DataTimeLength / length(PopN)
      # Apply data filter based on a series of criteria
      if ((length(PopN) >= data_length_min) & (AvgTimeBetweenPts < avg_time_between_pts_max)) {
        if (offset_all) {
          cat(sprintf("Offsetting all time-series by 1 to avoid log(0)\n"))
          PopN <- PopN + 1
        } else if (offset_none) {
          # Adds small value to zero values within a time series, To non-zero values it does nothing.
          IndexZero <- which(PopN == 0)
          if (length(IndexZero) > 0) {
            OffsetVal <- 1e-17
            PopN <- PopN + OffsetVal
          } else {
            PopN <- PopN
          }
        } else if (offset_diff) {
          # Offset different pops differentially
          IndexZero <- which(PopN == 0)
          # Check if this pop has 0 values
          if (length(IndexZero) > 0) {
            if (mean(PopN) == 0) {
              OffsetVal <- 1e-17
            } else {
              # If the maximum value of the population is more than to 1
              if (max(PopN) >= 1) {
                # Add 1
                PopN <- PopN + 1
              } else {
                # Otherwise add 1% of the mean
                IndexNonZero <- which(PopN != 0)
                OffsetVal <- mean(PopN[IndexNonZero]) * 0.01
                PopN <- PopN + OffsetVal
              }
            }
          }
        } else {
          # Replace zero values with 1 percent of average values
          IndexZero <- which(PopN == 0)

          # cat(sprintf("Number of zero values: %d\n", length(IndexZero)))

          if (zero_replace_flag == 1) {
            # cat(sprintf("Replacing zero values with 1 percent of average\n"))
            if (mean(PopN) == 0) {
              OffsetVal <- 1e-17
            } else {
              IndexNonZero <- which(PopN != 0)
              OffsetVal <- mean(PopN[IndexNonZero]) * 0.01
            }
          } else {
            # cat(sprintf("Replacing zero values with min of non-zero\n"))

            IndexNonZero <- which(PopN != 0)
            OffsetVal <- min(PopN[IndexNonZero])
          }

          if (zero_replace_flag == 2) {
            # cat(sprintf("Replacing zero values with 1\n"))
            if (mean(PopN) == 0) {
              OffsetVal <- 1e-17
            } else {
              IndexNonZero <- which(PopN != 0)
              OffsetVal <- 1
            }
          }
          if (length(IndexZero) > 0) {
            PopN <- PopN + OffsetVal
          }
        }


        # pdf()

        SortResults <- sort(yearPop, index.return = TRUE)
        yearPop <- SortResults$x
        TempI <- SortResults$ix
        PopN <- PopN[TempI]

        # Perform smoothing and interpolation

        if (length(which(objects() == "yearPopInt")) != 0) {
          rm(yearPopInt)
        }

        if (length(which(objects() == "PopNInt")) != 0) {
          rm(PopNInt)
        }

        yearPopInt <- yearPop[1]:yearPop[length(yearPop)]

        PopNLog <- log(PopN)
        Flag <- 0

        # if population has the same value then no need to build GAM model and interpolate
        GAMFlag <- gam_global_flag
        if (mean(PopN) == PopN[1]) {
          GAMFlag <- 0
        }

        if (GAMFlag == 1) {
          if (model_selection_flag == 0) {

            # Estimate the smoothing parameter to be half the population length
            SmoothParm <- round(length(PopN) / 2)
            # If this is 3 or more (The population is therefore 6+ datapoints)
            if (SmoothParm >= 3) {
              # Added this as was having trouble refering to mgcv::s in formulae
              s <- mgcv::`s`
              model <- mgcv::gam(PopNLog ~ s(yearPop, k = SmoothParm), fx = TRUE)
              # check if the model is ok
              if (auto_diagnostic_flag == 1) {
                rsd <- residuals(model)
                s <- mgcv::`s`
                modelres <- mgcv::gam(rsd ~ s(yearPop, k = length(PopN), bs = "cs"),
                  gamma = 1.4
                )
                if ((abs(sum(modelres$edf) - 1)) < 0.01) {
                  # 0.01
                  PopNInt <- predict(model, data.frame(yearPop = yearPopInt))
                  PopNInt <- exp(PopNInt)
                  Flag <- 1
                  PopProcessedGAMCounter <- PopProcessedGAMCounter + 1
                  PopProcessedGAM[PopProcessedGAMCounter] <- Popid[J]
                }
              } else {
                summary(model)
                readline(prompt = "Press any key to continue")
                plot(model,
                  pages = 1, residuals = TRUE, all.terms = TRUE, shade = TRUE,
                  shade.col = 2
                )
                readline(prompt = "Press any key to continue")
                mgcv::gam.check(model)
                Char <- readline(prompt = "Press 'Y' to accept model, 'N' to reject GAM model and use default method")
                while ((Char != "Y") & (Char != "N")) {
                  Char <- readline(prompt = "Press 'Y' to accept model, 'N' to reject GAM model and use default method")
                }
                if (Char == "Y") {
                  PopNInt <- predict(model, data.frame(yearPop = yearPopInt))
                  PopNInt <- exp(PopNInt)
                  Flag <- 1
                  PopProcessedGAMCounter <- PopProcessedGAMCounter + 1
                  PopProcessedGAM[PopProcessedGAMCounter] <- Popid[J]
                }
              }
            }
          } else {
            if (length(PopN) >= 6) {
              SmoothParm <- 3 # length(PopN) if K is set to max
              if (auto_diagnostic_flag == 1) {
                while ((length(PopN) >= SmoothParm) & (Flag == 0)) {
                  s <- mgcv::`s`
                  model <- mgcv::gam(PopNLog ~ s(yearPop, k = SmoothParm), fx = TRUE)
                  rsd <- residuals(model)
                  modelres <- mgcv::gam(rsd ~ s(yearPop, k = length(PopN), bs = "cs"),
                    gamma = 1.4
                  )
                  if ((abs(sum(modelres$edf) - 1)) < 0.01) {
                    Flag <- 1
                    PopNInt <- predict(model, data.frame(yearPop = yearPopInt))
                    PopNInt <- exp(PopNInt)
                    PopProcessedGAMCounter <- PopProcessedGAMCounter + 1
                    PopProcessedGAM[PopProcessedGAMCounter] <- Popid[J]
                  } else {
                    SmoothParm <- SmoothParm + 1
                  }
                }
              } else {
                while ((length(PopN) >= SmoothParm) & (Flag == 0)) {
                  s <- mgcv::`s`
                  model <- mgcv::gam(PopNLog ~ s(yearPop, k = SmoothParm), fx = TRUE)
                  summary(model)
                  readline(prompt = "Press any key to continue")
                  plot(model,
                    pages = 1, residuals = TRUE, all.terms = TRUE,
                    shade = TRUE, shade.col = 2
                  )
                  readline(prompt = "Press any key to continue")
                  mgcv::gam.check(model)
                  Char <- readline(prompt = "Press 'Y' to accept model, 'N' to reject model")
                  while ((Char != "Y") & (Char != "N")) {
                    Char <- readline(prompt = "Press 'Y' to accept model, 'N' to reject model")
                  }
                  if (Char == "Y") {
                    PopNInt <- predict(model, data.frame(yearPop = yearPopInt))
                    PopNInt <- exp(PopNInt)
                    Flag <- 1
                    PopProcessedGAMCounter <- PopProcessedGAMCounter + 1
                    PopProcessedGAM[PopProcessedGAMCounter] <- Popid[J]
                  } else {
                    SmoothParm <- SmoothParm + 1
                  }
                }
              }
            }
          }
        }
        # IF we get to here and the Flag is still 0, then no GAM has been made - either
        # the population is to short (<6) or the gam has failed the quality check
        if (Flag == 0) {

          # If this is true then try to GAM the short pops too
          if (global_gam_flag_short_data_flag == 1) {
            SmoothParm <- length(PopN)
            s <- mgcv::`s`
            model <- mgcv::gam(PopNLog ~ s(yearPop, k = SmoothParm), fx = TRUE)
            PopNInt <- predict(model, data.frame(yearPop = yearPopInt))
            PopNInt <- exp(PopNInt)
            PopProcessedGAMCounter <- PopProcessedGAMCounter + 1
            PopProcessedGAM[PopProcessedGAMCounter] <- Popid[J]
          } else {
            if (linear_model_short_flag == TRUE) {
              MethodFlagLoop <- MethodFlagLoop + 1
              MethodFlag[MethodFlagLoop] <- Popid[J]
              model <- lm(PopNLog ~ yearPop)
              # r2 <- summary(model)$r.squared
              # LM_R2_THRESH = 0.0
              # if (r2 > LM_R2_THRESH) {
              PopNInt <- predict(model, data.frame(yearPop = yearPopInt))
              PopNInt <- exp(PopNInt)
              # } else {
              # PopNotProcessedCounter = PopNotProcessedCounter + 1
              # PopNotProcessed[PopNotProcessedCounter] = Popid[J]
              # cat("R squared less than", LM_R2_THRESH, "\n")
              # next
              # }
            } else {
              # Apply the default approach (Chain)
              MethodFlagLoop <- MethodFlagLoop + 1
              MethodFlag[MethodFlagLoop] <- Popid[J]
              PopNInt <- matrix(-1, 1, length(yearPopInt))

              for (K in 1:length(yearPopInt)) {
                k <- which(yearPop == yearPopInt[K])
                if (length(k) > 0) {
                  PopNInt[K] <- PopN[k]
                } else {
                  # find the previous value
                  yearStart <- yearPopInt[K]
                  yearStart <- yearStart - 1
                  k <- which(yearPop == yearStart)
                  while (length(k) == 0) {
                    yearStart <- yearStart - 1
                    k <- which(yearPop == yearStart)
                  }
                  PopNStart <- PopN[k]
                  # find the next value
                  yearEnd <- yearPopInt[K]
                  yearEnd <- yearEnd + 1
                  k <- which(yearPop == yearEnd)
                  while (length(k) == 0) {
                    yearEnd <- yearEnd + 1
                    k <- which(yearPop == yearEnd)
                  }
                  PopNEnd <- PopN[k]
                  # Calculate the interpolated value
                  PopNInt[K] <- PopNStart * ((PopNEnd / PopNStart)^((yearPopInt[K] -
                    yearStart) / (yearEnd - yearStart)))
                }
              }
            }
          }
        }

        # only consider from initial_year onwards
        yearPop <- initial_year:final_year
        PopN <- matrix(0, 1, length(yearPop))

        k <- which(PopNInt == 0)
        k1 <- which(PopNInt > 0)
        TempVal <- 0
        if (length(k) > 0) {
          if (length(k1) > 0) {
            if (zero_replace_flag == 1) {
              TempVal <- mean(PopNInt[k1]) * 0.01
            } else {
              TempVal <- min(PopNInt[k1])
            }
            PopNInt <- PopNInt + TempVal
          }
        }
        for (K in initial_year:final_year) {
          k <- which(yearPopInt == K)
          if (length(k) > 0) {
            if (PopNInt[k] == 0) {
              PopN[K - initial_year + 1] <- -1
            } else {
              PopN[K - initial_year + 1] <- log10(PopNInt[k])
            }
          } else {
            PopN[K - initial_year + 1] <- -1
          }
        }

        # Calculate the growth rate
        PopLambda[JIndex, 1] <- 1
        Startyear <- initial_year + 1
        for (K in Startyear:final_year) {
          if ((PopN[K - initial_year + 1] != -1) & (PopN[K - initial_year] !=
            -1)) {
            PopLambda[JIndex, K - initial_year + 1] <- PopN[K - initial_year +
              1] - PopN[K - initial_year]
          } else {
            PopLambda[JIndex, K - initial_year + 1] <- -1
          }
        }
        JIndex <- JIndex + 1
      } else {
        PopNotProcessedCounter <- PopNotProcessedCounter + 1
        PopNotProcessed[PopNotProcessedCounter] <- Popid[J]
      }
    }

    # Save the population lamdas to a file:

    PopData <- cbind(as.vector(Popid), PopLambda)
    pop_lambda_filename <- file.path(basedir, gsub(".txt", "_PopLambda.txt", dataset_name))
    # Pop_Headers<-t(c("population_id", as.vector(initial_year:final_year)))
    # write.table(Pop_Headers,file=pop_lambda_filename, sep=",", eol="\n", quote=FALSE, col.names=FALSE, row.names = FALSE)
    write.table(PopData, sep = ",", eol = "\n", file = pop_lambda_filename, quote = FALSE, col.names = FALSE, row.names = FALSE, append = TRUE)


    # Save the species average lambda values

    Endyear <- final_year - initial_year + 1
    for (K in 1:Endyear) {
      k <- which(PopLambda[, K] != -1)
      if (length(k) > 0) {
        # Only get those lambdas that are not '-1'
        PopLambdaTemp <- PopLambda[k, K]
        # Fine which of these are less than our max
        IndexTemp <- which(PopLambdaTemp < lambda_max)
        IndexTempBad_max <- which(PopLambdaTemp > lambda_max)
        # If we have some...
        if (length(IndexTemp) > 0) {
          # Extract them as PopLambdaTemp1
          PopLambdaTemp1 <- PopLambdaTemp[IndexTemp]
          # Get those values that are also more then our min
          IndexTemp <- which(PopLambdaTemp1 > lambda_min)
          IndexTempBad_min <- which(PopLambdaTemp1 < lambda_min)
          # If there are some...
          if (length(IndexTemp) > 0) {
            # Then set species lambda to be their average...
            if (cap_lambdas) {
              speciesLambda[I, K] <- mean(c(PopLambdaTemp1[IndexTemp], rep(lambda_max, length(IndexTempBad_max)), rep(lambda_min, length(IndexTempBad_min))))
            } else {
              speciesLambda[I, K] <- mean(PopLambdaTemp1[IndexTemp])
            }
          } else {
            # Otherwise, if we have lambdas less than our max, but not more then min, set sp. av. to be NA
            if (cap_lambdas) {
              speciesLambda[I, K] <- lambda_min
            } else {
              speciesLambda[I, K] <- NA
            }
          }
        } else {
          # Otherwise, if we have no lambdas less than our max, set sp. av. to be NA
          if (cap_lambdas) {
            speciesLambda[I, K] <- lambda_max
          } else {
            speciesLambda[I, K] <- NA
          }
        }
      } else {
        # If all values are '-1' then set species lambda to be NA
        speciesLambda[I, K] <- NA
      }
    }

    sNamesCounter <- sNamesCounter + 1
    sNamesArray[sNamesCounter, 1] <- sNames[I, 1]
    sidArray[sNamesCounter] <- id[I, 1]

    # cat('Results saved\n')
    if (show_progress) setTxtProgressBar(prog, I)
  }
  if (show_progress) close(prog)
  cat("\n")

  PopNotProcessed1 <- PopNotProcessed[1, 1:PopNotProcessedCounter]
  write.table(PopNotProcessed1, file = file.path(basedir, "lpi_temp", "PopNotProcessed.txt")) # insert your desired file name here
  MethodFlag1 <- MethodFlag[1, 1:MethodFlagLoop]
  if (linear_model_short_flag == 1) {
    write.table(MethodFlag1, file = file.path(basedir, "lpi_temp", "PopProcessed_LM.txt")) # insert your desired file name here
  } else {
    write.table(MethodFlag1, file = file.path(basedir, "lpi_temp", "PopProcessed_Chain.txt")) # insert your desired file name here
  }
  PopProcessedGAM1 <- PopProcessedGAM[1, 1:PopProcessedGAMCounter]
  write.table(PopProcessedGAM1, file = file.path(basedir, "lpi_temp", "PopProcessedGAM.txt")) # insert your desired file name here
  sNamesArray1 <- sNamesArray[1:sNamesCounter, 1]
  write.table(sNamesArray1, file = file.path(basedir, "lpi_temp", "SpeciesName.txt"), quote = FALSE)

  Headers <- t(c("Species", as.vector(initial_year:final_year)))

  # sNamesT <- sNamesArray[sidArray, 1]

  # ids (sidArray) is per-population and sNamesArray is per species... currently exporting species lambdas
  # speciesData<-cbind(sidArray, as.vector(sNamesT), speciesLambda)
  speciesData <- cbind(as.vector(sNamesArray), speciesLambda)
  lambda_filename <- file.path(basedir, gsub(".txt", "_Lambda.txt", dataset_name))
  write.table(Headers, file = lambda_filename, sep = ",", eol = "\n", quote = FALSE, col.names = FALSE, row.names = FALSE)
  write.table(speciesData, sep = ",", eol = "\n", file = lambda_filename, quote = FALSE, col.names = FALSE, row.names = FALSE, append = TRUE)

  # sp_df = melt(speciesData)


  return(speciesLambda)
  # MethodFlag
}
