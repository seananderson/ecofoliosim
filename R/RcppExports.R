# This file was generated by Rcpp::compileAttributes
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Discrete-time Lotka-Volterra competition model from Loreau et al. 2013
#'
#' Simulate from an expanded discrete-time Lotka-Volterra competition model as
#' described in Loreau et al. 2013. This version uses C++ for speed.
#'
#' @details
#' Note that synchrony in environmental dynamics and demographic processes can
#' be added through the matrices \code{u_d} and \code{u_e}. Temporal
#' autocorrelation can be added through these matrices as well.
#' @useDynLib ecofoliosim
#'
#' @param rm A numeric vector of intrinsic maximum rate of natural increase for
#'   each species.
#' @param k A numeric vector of carrying capacities for each species.
#' @param b A numeric matrix of competition coefficients between species i
#'   (columns) and species j (rows).
#' @param n0 A numeric vector of starting biomass for each species.
#' @param sigma_d A vector of demographic standard deviation magnitudes for each
#'   species.
#' @param sigma_e A vector of environmental standard deviation magnitudes for
#'   each species.
#' @param u_d A matrix of standard normal N(0, 1) demographic deviations. Rows
#'   are time steps and columns are species.
#' @param u_e A matrix of standard normal N(0, 1) environmental deviations. Rows
#'   are time steps and columns are species.
#' @export
#' @references
#' Loreau, M. and de Mazancourt, C. (2008). Species synchrony and its drivers:
#' Neutral and nonneutral community dynamics in fluctuating environments. Amer.
#' Nat., 172(2):E48-E66.
#'
#' Loreau, M. (2010). From Populations to Ecosystems: Theoretical Foundations
#' for a New Ecological Synthesis. Princeton University Press, Princeton, NJ.
#'
#' Loreau, M. and de Mazancourt, C. (2013). Biodiversity and ecosystem
#' stability: a synthesis of underlying mechanisms. Ecol. Lett., 16(S1):106-115.
#'
#' @examples
#' # Figure 1b in Loreau et al. 2013:
#' species <- 2
#' time_steps <- 550
#' burnin <- 50
#' u_d <- matrix(ncol = species, nrow = time_steps,
#'   data = rnorm(species * time_steps, 0, 1))
#' u_e <- matrix(ncol = species, nrow = time_steps,
#'   data = rnorm(species * time_steps, 0, 1))
#' m <- lotvolt2(rm = c(0.5, 0.8), k = c(1000, 1500),
#'   b = rbind(c(NA, 0.7), c(0.9, NA)),
#'   n0 = c(1000, 1500), sigma_d = c(1, 1), sigma_e = c(0.02, 0.02),
#'   u_d = u_d, u_e = u_e)
#' matplot(m[-seq_len(burnin), ], type = "l", ylab = "Biomass",
#'   lty = 1, ylim = c(0, 1000), xlab = "Time")
lotvolt2 <- function(rm, k, b, n0, sigma_d, sigma_e, u_d, u_e) {
    .Call('ecofoliosim_lotvolt2', PACKAGE = 'ecofoliosim', rm, k, b, n0, sigma_d, sigma_e, u_d, u_e)
}

