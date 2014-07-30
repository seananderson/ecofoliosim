#include <Rcpp.h>
using namespace Rcpp;

//' Discrete-time Lotka-Volterra competition model from Loreau et al. 2013
//'
//' Simulate from an expanded discrete-time Lotka-Volterra competition model as
//' described in Loreau et al. 2013. This version uses C++ for speed.
//'
//' @details
//' Note that synchrony in environmental dynamics and demographic processes can
//' be added through the matrices \code{u_d} and \code{u_e}. Temporal
//' autocorrelation can be added through these matrices as well.
//' @useDynLib ecofoliosim
//'
//' @param rm A numeric vector of intrinsic maximum rate of natural increase for
//'   each species.
//' @param k A numeric vector of carrying capacities for each species.
//' @param b A numeric matrix of competition coefficients between species i
//'   (columns) and species j (rows).
//' @param n0 A numeric vector of starting biomass for each species.
//' @param sigma_d A vector of demographic standard deviation magnitudes for each
//'   species.
//' @param sigma_e A vector of environmental standard deviation magnitudes for
//'   each species.
//' @param u_d A matrix of standard normal N(0, 1) demographic deviations. Rows
//'   are time steps and columns are species.
//' @param u_e A matrix of standard normal N(0, 1) environmental deviations. Rows
//'   are time steps and columns are species.
//' @export
//' @references
//' Loreau, M. and de Mazancourt, C. (2008). Species synchrony and its drivers:
//' Neutral and nonneutral community dynamics in fluctuating environments. Amer.
//' Nat., 172(2):E48-E66.
//'
//' Loreau, M. (2010). From Populations to Ecosystems: Theoretical Foundations
//' for a New Ecological Synthesis. Princeton University Press, Princeton, NJ.
//'
//' Loreau, M. and de Mazancourt, C. (2013). Biodiversity and ecosystem
//' stability: a synthesis of underlying mechanisms. Ecol. Lett., 16(S1):106-115.
//'
//' @examples
//' # Figure 1b in Loreau et al. 2013:
//' species <- 2
//' time_steps <- 550
//' burnin <- 50
//' u_d <- matrix(ncol = species, nrow = time_steps,
//'   data = rnorm(species * time_steps, 0, 1))
//' u_e <- matrix(ncol = species, nrow = time_steps,
//'   data = rnorm(species * time_steps, 0, 1))
//' m <- lotvolt2(rm = c(0.5, 0.8), k = c(1000, 1500),
//'   b = rbind(c(NA, 0.7), c(0.9, NA)),
//'   n0 = c(1000, 1500), sigma_d = c(1, 1), sigma_e = c(0.02, 0.02),
//'   u_d = u_d, u_e = u_e)
//' matplot(m[-seq_len(burnin), ], type = "l", ylab = "Biomass",
//'   lty = 1, ylim = c(0, 1000), xlab = "Time")
// [[Rcpp::export]]
NumericMatrix lotvolt2(NumericVector rm, NumericVector k, NumericMatrix b,
    NumericVector n0, NumericVector sigma_d, NumericVector sigma_e,
    NumericMatrix u_d, NumericMatrix u_e) {

  int n_species;
  int n_steps;
  int n_steps_minus_one;
  double det_a;
  double det_b;
  double stoc_e;
  double stoc_d;

  n_species = u_d.ncol();
  //Rcpp::Rcout << "n_species" << n_species << std::endl;
  n_steps = u_d.nrow();
  //Rcpp::Rcout << "n_steps" << n_steps << std::endl;
  n_steps_minus_one = n_steps - 1;
  NumericMatrix r(n_steps, n_species);
  NumericMatrix n(n_steps, n_species);

  n(0, _) = n0;

  for(int st = 0; st < n_steps_minus_one; st++) {
  //Rcpp::Rcout << "det_b" << det_b << std::endl;
    for(int i = 0; i < n_species; i++) {
      det_a = n(st, i) / k(i);
      det_b = 0.0; // reset to 0 before adding
      for(int j = 0; j < n_species; j++)  {
        if(i != j) {
          det_b = det_b + (b(i, j) * n(st, j)) / k(j);
        }
      }
      stoc_e = sigma_e(i) * u_e(st, i);
      stoc_d = (sigma_d(i) * u_d(st, i)) / (sqrt(n(st, i)));
      r(st, i) = rm(i) * (1 - det_a - det_b) + stoc_e + stoc_d;
      n(st + 1, i) = exp(r(st, i) + log(n(st, i)));
    }
  }
  return n;
}
