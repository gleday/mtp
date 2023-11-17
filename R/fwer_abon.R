#' FWER control using adaptive Bonferroni (single-step)
#'
#' @inheritParams fwer_bon
#' @inheritParams d0
#'
#' @description
#' Control \eqn{\text{FWER(k)}} using generalization of
#' adaptive Bonferroni's procedure (Wang, 2017)
#'
#' @details
#' The generalized adaptive Bonferroni procedure
#' (Wang, 2017; Definition 2.1)
#' consists in using the decision procedure
#' described in [mtp-package] with:
#'
#' * the adjustment factor:
#' \eqn{\qquad\quad
#'  \displaystyle{
#'   a = \frac{\widehat{d}_0}{k}
#'  },
#' }
#' where \eqn{\widehat{d}_0} is Storey's estimator of the
#' number \eqn{d_0} of true null hypotheses (see [d0()]).
#'
#' * the adjusted P-values:
#' \eqn{\qquad\quad
#'  \displaystyle{
#'   \widetilde{p}_{j} = \min\left( a\times p_{j}, 1\right),
#'   \ \text{for}\ j=1, \ldots, d.
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
#' The generalized adaptive Bonferroni procedure guarantees
#' that \eqn{\text{FWER(k)} \leq \alpha} without
#' assumptions on the dependence of P-values.
#'
#' @importFrom "assertthat" "assert_that" "is.count" "is.number" "is.string"
#' @importFrom "dplyr" "between"
#' @importFrom "purrr" "map_dfc" "map_lgl"
#'
#' @return
#' A [numeric] vector of:
#' * adjusted P-values when `.return = "p"`,
#' * adjustment factors when `.return = "a"`,
#' * adjusted critical values when `alpha` is provided.
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
#' @export
fwer_abon <- function(p_value, k = 1, lambda = 0.5, .return = "p", alpha = NULL) {

  # check arguments
  .check_p_value()
  .check_k()
  .check_lambda()
  .check_return()

  # get adjustment factor
  a <- d0(p_value = p_value, lambda = lambda) / k

  # output
  p <- pmin(a * p_value, 1)
  if (!is.null(alpha)) {
    return(rep(alpha / a, length(p_value)))
  } else {
    if (.return == "a") {
      return(rep(a, length(p_value)))
    }
  }
  p
}
