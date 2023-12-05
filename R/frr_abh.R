#' FRR control using adaptive Benjamini-Hochberg (step-up)
#'
#' @inheritParams frr_bh
#' @inheritParams d0
#'
#' @description
#' Control \eqn{\text{FRR}} using adaptive
#' Benjamini-Hochberg's step-up procedure and
#' Storey's plugin estimator
#'
#' @details
#' Storey's adaptive version of
#' Benjamini-Hochberg (BH) procedure
#' consists in using the decision procedure
#' described in [mtp-package] using:
#'
#' * the adjustment factors:
#' \eqn{\qquad\quad
#'  \displaystyle{
#'   a_j = \frac{\hat{d}_0}{j},
#'   \ \text{for}\ j=1, \ldots, d.
#'  },
#' }\cr
#' where \eqn{\widehat{d}_0} is Storey's estimator of the
#' number \eqn{d_0} of true null hypotheses (see [d0()]).
#'
#' * the adjusted P-values:
#' \eqn{\qquad\quad
#' \displaystyle{
#'   \begin{cases}
#'   \widetilde{p}_{1} = \min\left( a_j p_{j}, 1\right),\\
#'   \widetilde{p}_{j} = \min\left( a_j p_{j},
#'   \widetilde{p}_{j + 1}\right), \ \text{for}\ j = 1, \ldots, d-1
#'   \end{cases}
#'  }
#'  }
#'
#' * the adjusted critical values:
#' \eqn{\qquad
#' \displaystyle{
#'   \widetilde{\alpha}_{j} = \frac{\alpha}{a_j},
#'   \ \text{for}\ j=1, \ldots, d.
#' } }
#'
#' The BH procedure guarantees
#' that \eqn{\text{FRR} \leq \alpha} under
#' independence or positive dependence of P-values.
#'
#' @importFrom "assertthat" "assert_that" "is.count" "is.number" "is.string"
#' @importFrom "dplyr" "between"
#' @importFrom "purrr" "map_lgl"
#'
#' @return
#' A [numeric] vector of:
#' * adjusted P-values when `.return = "p"`,
#' * adjustment factors when `.return = "a"`,
#' * adjusted critical values when `alpha` is provided.
#'
#' @family FRR
#'
#' @references
#' Storey, J. D. (2002). A direct approach to false
#' discovery rates. Journal of the Royal Statistical
#' Society Series B: Statistical Methodology, 64(3), 479-498.\cr
#' Storey, J. D., Taylor, J. E., & Siegmund, D. (2004).
#' Strong control, conservative point estimation and
#' simultaneous conservative consistency of false
#' discovery rates: a unified approach.
#' Journal of the Royal Statistical Society Series B:
#' Statistical Methodology, 66(1), 187-205.
#'
#' @export
frr_abh <- function(p_value, lambda = 0.5, .return = "p", alpha = NULL) {

  # check arguments
  .check_p_value()
  .check_lambda()
  .check_return()

  # get adjustment factors
  d <- length(p_value)
  j <- d:1L
  o <- order(p_value, decreasing = TRUE)
  ro <- order(o)
  a <- d0(p_value = p_value, lambda = lambda) / j

  # output
  p <- pmin(cummin(a * p_value[o]), 1)[ro]
  if (!is.null(alpha)) {
    return((alpha * p_value) / p)
  } else {
    if (.return == "a") {
      return(p / p_value)
    }
  }
  p
}
