#' Discrete-time Lotka-Volterra competition model from Loreau et al. 2013
#'
#' Simulate from an expanded discrete-time Lotka-Volterra competition model as
#' described in Loreau et al. 2013.
#'
#' @details
#' Note that synchrony in environmental dynamics and demographic processes can
#' be added through the matrices \code{u_d} and \code{u_e}. Temporal
#' autocorrelation can be added through these matrices as well.
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
#' @seealso \code{\link{lotvolt2}}
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
#' m <- lotvolt(rm = c(0.5, 0.8), k = c(1000, 1500),
#'   b = rbind(c(NA, 0.7), c(0.9, NA)),
#'   n0 = c(1000, 1500), sigma_d = c(1, 1), sigma_e = c(0.02, 0.02),
#'   u_d = u_d, u_e = u_e)
#' matplot(m[-seq_len(burnin), ], type = "l", ylab = "Biomass",
#'   lty = 1, ylim = c(0, 1000), xlab = "Time")

lotvolt <- function(rm, k, b, n0, sigma_d, sigma_e, u_d, u_e) {
  n_species <- length(rm)
  n_steps <- nrow(u_d)
  if(nrow(u_d) != nrow(u_e)) stop("Number of rows in u_d and u_e don't match")
  r <- matrix(ncol = n_species, nrow = n_steps)
  n <- matrix(ncol = n_species, nrow = n_steps)
  n[1, ] <- n0
  for(st in seq_len(n_steps - 1)) {
    for(i in seq_len(n_species)) {
      det_a <- n[st, i] / k[i]
      det_b <- 0
        for(j in seq_len(n_species)[-i]) {
          det_b <- det_b + (b[i, j] * n[st, j]) / k[j]
        }
      stoc_e <- sigma_e[i] * u_e[st, i]
      stoc_d <- (sigma_d[i] * u_d[st, i]) / (sqrt(n[st, i]))
      r[st, i] <- rm[i] * (1 - det_a - det_b) + stoc_e + stoc_d
      n[st + 1, i] <- exp(r[st, i] + log(n[st, i]))
    }
  }
  n
}


#' Run Lotka-Voltera comparison simulation
#'
#' A wrapper function to facilitate running simulations for two species.
#'
#' @param w1 Weight for species 1
#' @param b Competition matrix
#' @param obs_errors Either 0 (no observation errors) or a matrix of observation
#' errors
#' @param u_d Mean 0, standard deviation 1 demographic deviations.
#' @param u_e Mean 0, standard deviation 1 environmental deviations.
#' @param burnin How many iterations to discard as burnin.
#' @param sample_scale How often to sample the 'truth'. If set to 1 then all
#' observations are kept. If set to 2, for example, then every second
#' observation is discarded.
#' @param sigma_e Environmental standard deviation values to pass on.
#'
#' @export
lotvolt_comp <- function(w1, b, obs_errors = 0, u_d, u_e, burnin = 100, sample_scale = 1L, sigma_e = 0.03) {
  w2 <- 1 - w1
  out <- data.frame(w1 = w1, m = NA, v = NA)
  for(i in seq_along(w1)) {
    m <- lotvolt2(rm = c(0.8, 0.8), k = c(1000*w1[i], 1000*w2[i]),
      b = b,
      n0 = c(1000*w1[i], 1000*w2[i]), sigma_d = c(1, 1), sigma_e = c(sigma_e, sigma_e),
      u_d = u_d, u_e = u_e)
    m <- m[-seq_len(burnin), ]
    if(obs_errors[1] != 0) {
      m <- m * obs_errors
    }
    ret <- diff(log(rowSums(m)))
    if(sample_scale < 1) stop("sample_scale must be > 1 and be an integer")
    if(sample_scale > 1) {
      ret <- ret[seq(1, length(ret), sample_scale)]
    }
    out$m[i] <- mean(ret)
    out$v[i] <- sd(ret)
  }
  return(na.omit(out))
}

#' Run Lotka-Voltera comparison simulation for many species
#'
#' A wrapper function to facilitate running simulations where the number of
#' species can be more than 2.
#'
#' @param n_sim Number of Monte Carlo simulations to run
#' @param b Competition value (same for all species)
#' @param obs_errors Either 0 (no observation errors) or a matrix of observation
#' errors
#' @param u_d Mean 0, standard deviation 1 demographic deviations.
#' @param u_e Mean 0, standard deviation 1 environmental deviations.
#' @param burnin How many iterations to discard as burn in.
#' @param ef_interval How many bins to divide the efficient frontier into when
#' calculating the frontier based on the Monte Carlo results.
#' @importFrom metafolio create_asset_weights
#'
#' @export
lotvolt_sp <- function(n_sim, b, obs_errors = 0, u_d, u_e,
  burnin = 100, ef_interval = 20) {
  weights <- create_asset_weights(n_pop = ncol(u_d),
    n_sims = n_sim, weight_lower_limit = 0.001)
  N <- ncol(weights)
  junk <- data.frame(m = rep(NA, nrow(weights)), v = NA)
  for(i in seq_len(nrow(weights))) {
    m <- lotvolt2(rm = rep(0.8, N), k = 1000 * weights[i,],
      b = matrix(nrow = N, ncol = N, data = b),
      n0 = 1000 * weights[i,], sigma_d = rep(1, N),
      sigma_e = rep(0.03, N), u_d = u_d, u_e = u_e)
    m <- m[-seq_len(burnin), ]
    if(obs_errors[1] != 0) {
      m <- m * obs_errors
    }
    ret <- diff(log(rowSums(m)))
    junk$m[i] <- mean(ret)
    junk$v[i] <- sd(ret)
  }
  # let's pull out the efficient frontier:
  # [code from metafolio]
  junk <- na.omit(junk)
  m.bins <- seq(min(junk$m), max(junk$m), length.out = ef_interval)
  junk$m.bin <- findInterval(junk$m, m.bins)
  junk$optim_set <- 0
  for(i in unique(junk$m.bin)) {
    junk_to_check <- which(junk$m.bin == i)
    junk[junk_to_check, "optim_set"][which.min(junk[junk_to_check, ]$v)] <- 1
  }
  junk
}


#' Plot efficient frontiers
#'
#' @param dat A list of data frames containing columsn for variability (\code{v})
#' and mean (\code{m}).
#' @param xlim Possible x limits
#' @param ylim Possible y limits
#' @param show Vector of portfolios to show
#' @param colour_lines Should the lines be coloured?
#' @param ... Any extra arguments to pass to \code{plot.default}.
#'
#' @export
plot_frontiers <- function(dat, xlim = NULL, ylim = NULL, show = seq_along(dat),
  colour_lines = FALSE, ...) {
  if(is.null(xlim)) {
    xlim <- c(min(unlist(lapply(dat, function(x) min(x$v)))),
    max(unlist(lapply(dat, function(x) max(x$v)))))
  }
  if(is.null(ylim)) {
    ylim <- c(min(unlist(lapply(dat, function(x) min(x$m)))),
    max(unlist(lapply(dat, function(x) max(x$m)))))
  }
  plot(1, 1, xlim = xlim, ylim = ylim, type = "n", xlab = "Standard devation",
    ylab = "Mean", ...)
  if(is.character(show)) {
    if(show == "all") {
      to_show <- seq_along(dat)
    }
  }
  for(i in show) {
    if(colour_lines) {
      lines(dat[[i]]$v, dat[[i]]$m, lwd = 2, col = i)
    } else {
      lines(dat[[i]]$v, dat[[i]]$m, col = "#00000040", lwd = 3)
    }
  }
  invisible(list(xlim = xlim, ylim = ylim))
}
