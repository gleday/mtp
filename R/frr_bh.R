#' FRR control using Benjamini-Hochberg (step-up)
#'
#' @inheritParams fwer_bon
#' @param alpha [numeric] scalar. Target \eqn{\text{FRR}}
#' level (between 0 and 1). Overrides `output` and
#' return adjusted critical values.
#'
#' @description
#' Control \eqn{\text{FRR}} using Benjamini-Hochberg's
#' step-up procedure
#'
#' @details
#' The Benjamini-Hochberg (BH) procedure
#' (Benjamini and Hochberg, 1995)
#' consists in using the decision procedure
#' described in [mtp-package] using:
#'
#' * the adjustment factors:
#' \eqn{\qquad\quad
#'  \displaystyle{
#'   a_j = \frac{m}{j},
#'   \ \text{for}\ j=1, \ldots, m.
#'  }
#' }
#'
#' * the adjusted P-values:
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
#' assumptions on the dependence of P-values
#' (Simes inequality).
#'
#' @importFrom "assertthat" "assert_that" "is.count" "is.number" "is.string"
#' @importFrom "dplyr" "between"
#' @importFrom "purrr" "map_lgl"
#'
#' @return
#' A [numeric] vector of:
#' * adjusted P-values when `ouput = "p"`,
#' * adjustment factors when `ouput = "a"`,
#' * adjusted critical values when `alpha` is provided.
#'
#' @family FRR
#'
#' @references
#' Benjamini, Y., & Hochberg, Y. (1995). Controlling the
#' false discovery rate: a practical and powerful approach to
#' multiple testing. Journal of the Royal statistical society:
#' series B (Methodological), 57(1), 289-300.
#'
#' @export
frr_bh <- function(p, alpha = NULL, output = "p") {

  # check arguments
  .check_p()
  .check_alpha()
  .check_output()

  # get adjustment factors
  m <- length(p)
  j <- m:1L
  o <- order(p, decreasing = TRUE)
  ro <- order(o)
  a <- m / j

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
