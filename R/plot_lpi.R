#' plot_lpi
#'
#' @param index - The index to plot
#' @param ref_year - The reference year of the plot (index == 1)
#' @param plot_max - The max y-value of the plot
#' @param ci_flag - whether confidence intervals are to be plotted
#' @param lower_ci - lower confidence interval values
#' @param upper_ci - upper confidence interval values
#' @param col - The color of the plot. Default is "black"
#'
#' @export
#'
plot_lpi <- function(index, ref_year, plot_max, ci_flag = 0, lower_ci = 0, upper_ci = 0, col = "black") {
  # plot the data
  year <- seq(ref_year, (ref_year + length(index)) - 1)
  plot(year, index, xlim = c(ref_year, plot_max), ylim = c(0, 2), ylab = paste("Index (", ref_year, "= 1.0)", sep = ""), col = col)
  # plot(year, index, xlim = c(ref_year, plot_max), ylab = paste("Index (", ref_year, " = 1.0)", sep=""))
  zeroEffectLine <- rep(1, (length(index)))
  lines(year, zeroEffectLine, col = "black")
  lines(year, index, col = col)

  if (ci_flag == 1) {
    lines(year, lower_ci, col = col)
    lines(year, upper_ci, col = col)
  }
}
