#' FRR control using adaptive Benjamini-Hochberg (step-up)
#'
#' @inheritParams frr_bh
#' @inheritParams m0
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
#'   a_j = \frac{\hat{m}_0}{j},
#'   \ \text{for}\ j=1, \ldots, m.
#'  },
#' }\cr
#' where \eqn{\widehat{m}_0} is Storey's estimator of the
#' number \eqn{m_0} of true null hypotheses (see [m0()]).
#'
#' * the adjusted p-values:
#' \eqn{\qquad\quad
#' \displaystyle{
#'   \begin{cases}
#'   \widetilde{p}_{1} = \min\left( a_j p_{j}, 1\right),\\
#'   \widetilde{p}_{j} = \min\left( a_j p_{j},
#'   \widetilde{p}_{j + 1}\right), \ \text{for}\ j = 1, \ldots, m-1
#'   \end{cases}
#'  }
#'  }
#'
#' * the adjusted critical values:
#' \eqn{\qquad
#' \displaystyle{
#'   \widetilde{\alpha}_{j} = \frac{\alpha}{a_j},
#'   \ \text{for}\ j=1, \ldots, m.
#' } }
#'
#' The BH procedure guarantees
#' that \eqn{\text{FRR} \leq \alpha} under some
#' assumptions on the dependence of p-values
#' (Simes inequality).
#'
#' @importFrom "assertthat" "assert_that" "is.count" "is.number" "is.string"
#' @importFrom "dplyr" "between"
#' @importFrom "purrr" "map_lgl"
#'
#' @return
#' A [numeric] vector of:
#' * adjusted p-values when `output = "p"`,
#' * adjustment factors when `output = "a"`,
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
frr_abh <- function(p, lambda = 0.5, alpha = NULL, output = "p") {

  # check arguments
  .check_p()
  .check_lambda()
  .check_alpha()
  .check_output()

  # get adjustment factors
  m <- length(p)
  j <- m:1L
  o <- order(p, decreasing = TRUE)
  ro <- order(o)
  a <- m0(p = p, lambda = lambda) / j

  # output
  p_adj <- pmin(cummin(a * p[o]), 1)[ro]
  if (!is.null(alpha)) {
    return((alpha * p) / p_adj)
  } else {
    if (output == "a") {
      return(p_adj / p)
    }
  }
  p_adj
}
