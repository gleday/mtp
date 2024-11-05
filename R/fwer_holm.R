#' FWER control using Holm (step-down)
#'
#' @inheritParams fwer_bon
#'
#' @description
#' Control \eqn{\text{FWER(k)}} using generalization of
#' Holm's step-down procedure
#'
#' @details
#' The generalized Holm procedure
#' (Lehmann & Romano, 2005; Theorem 2.2)
#' yields:
#'
#' * the adjustment factors:
#' \eqn{\qquad\quad
#'  \displaystyle{
#'   a_j = \frac{m - \max\left(j - k - 1, 0\right)}{k + 1},
#'   \ \text{for}\ j=1, \ldots, m.
#'  }
#' }
#'
#' * the adjusted P-values:
#' \eqn{\qquad\quad
#' \displaystyle{
#'   \begin{cases}
#'   \widetilde{p}_{(1)} = \min\left( a_1 p_{(1)}, 1\right),\\
#'   \widetilde{p}_{(j)} = \max\left( a_j p_{(j)},
#'   \widetilde{p}_{(j - 1)}\right), \ \text{for}\ j=2, \ldots, m
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
#' @importFrom "assertthat" "assert_that" "is.count" "is.number" "is.string"
#' @importFrom "dplyr" "between"
#' @importFrom "purrr" "map_lgl"
#'
#' @return
#' A [numeric] vector of:
#' * adjusted P-values when `output = "p"`,
#' * adjustment factors when `output = "a"`,
#' * adjusted critical values when `alpha` is provided.
#'
#' @family FWER
#'
#' @references
#' Holm, S. (1979). A simple sequentially rejective
#' multiple test procedure. Scandinavian journal of
#' statistics, 65-70.\cr
#' Lehmann, E. L., & Romano, J. P. (2005).
#' Generalizations of the familywise error rate.
#' The Annals of Statistics, 33(3), 1138-1154.
#'
#' @export
fwer_holm <- function(p, k = 0, alpha = NULL, output = "p") {

  # check input
  m <- length(p)
  .check_p()
  .check_k()
  .check_alpha()
  .check_output()

  # get adjustment factors
  j <- seq_len(m)
  o <- order(p)
  ro <- order(o)
  a <- (m - pmax(0, j - k - 1)) / (k + 1)

  # output
  p_adj <- pmin(cummax(a * p[o]), 1)[ro]
  if (!is.null(alpha)) {
    return((alpha * p) / p_adj)
  } else {
    if (output == "a") {
      return(p_adj / p)
    }
  }
  p_adj
}
