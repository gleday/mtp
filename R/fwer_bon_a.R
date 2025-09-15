#' FWER control using adaptive Bonferroni (single-step)
#'
#' @inheritParams fwer_bon
#' @param lambda \[ `numeric(1)` \]\cr
#' A number between 0 and 1 specifying the value for
#' the tuning parameter of Storey's estimator.
#'
#' @description
#' Control FWER(k) using the generalized
#' adaptive Bonferroni procedure of Wang (2017).
#'
#' @details
#' The generalized adaptive Bonferroni procedure
#' (Wang, 2017; Definition 2.1)
#' yields:
#'
#' * the adjustment factor:
#' \eqn{\qquad\quad
#'  \displaystyle{
#'   a = \frac{\widehat{m}_0}{k}
#'  },
#' }\cr
#' where \eqn{\widehat{m}_0} is Storey's estimator of the
#' number \eqn{m_0} of true null hypotheses (see [m0()]).
#'
#' * the adjusted P-values:
#' \eqn{\qquad\quad
#'  \displaystyle{
#'   \widetilde{p}_{j} = \min\left( a\times p_{j}, 1\right),
#'   \ \text{for}\ j=1, \ldots, m.
#'  }
#' }
#'
#' * the adjusted critical value:
#' \eqn{\qquad
#'  \displaystyle{
#'   \widetilde{\alpha} = \frac{\alpha}{a}.
#'  }
#' }
#'
#' @inherit fwer_bon return
#'
#' @family FWER adaptive
#'
#' @author Gwenaël G.R. Leday
#'
#' @references
#' Storey, J. D. (2002). A direct approach to false discovery rates.
#' Journal of the Royal Statistical Society:
#' Series B (Statistical Methodology), 64(3), 479-498.\cr
#' Wang, L. (2017). Adaptive procedure for generalized
#' familywise error rate control. Communications in
#' Statistics-Simulation and Computation, 46(10), 8140-8151.
#'
#' @examples
#' # simulate p-values
#' set.seed(123)
#' sim <- simulate_p_values(m = 1000, m0 = 900, s1 = 0.1, s2 = 10)
#'
#' # 1. control FWER(k = 1) at 0.05 using adjusted p-values
#' adj_pv <- fwer_bon_a(sim$p_value, k = 1)
#' sum(adj_pv <= 0.05)
#'
#' # 2. control FWER(k = 1) at 0.05 using adjusted critical values
#' adj_cv <- fwer_bon_a(
#'   sim$p_value,
#'   k = 1,
#'   c_value = 0.05,
#'   output = "c_values"
#' )
#' sum(sim$p_value <= adj_cv)
#'
#' # 3. control FWER(k = 1) at 0.05 using decision indicator
#' ind <- fwer_bon_a(
#'   sim$p_value,
#'   k = 1,
#'   c_value = 0.05,
#'   output = "decisions"
#' )
#' sum(ind)
#'
#' # 4. adjustment factors to control FWER(k = 1) at 0.05
#' adj_fct <- fwer_bon_a(sim$p_value, k = 1, output = "factors")
#'
#' @export
fwer_bon_a <- function(
    p_values,
    k = 0,
    lambda = 0.5,
    c_value = NULL,
    output = c("p_values", "c_values", "decisions", "factors")) {
  # check arguments
  m <- length(p_values)
  .check_p_values()
  .check_k()
  .check_lambda()
  output <- arg_match(
    output,
    c("p_values", "c_values", "decisions", "factors")
  )
  .check_c_value()

  # get adjustment factor
  a <- m0(p_values = p_values, lambda = lambda) / (k + 1)

  # output
  adj_p_values <- pmin(a * p_values, 1)
  switch(output,
    p_values = adj_p_values,
    c_values = rep(c_value / a, m),
    decisions = (adj_p_values <= c_value) * 1L,
    factors = rep(a, m)
  )
}
