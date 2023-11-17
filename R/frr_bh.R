#' FRR control using Benjamini-Hochberg (step-up)
#'
#' @inheritParams fwer_bon
#'
#' @description
#' Control \eqn{\text{FRR}} using Benjamini-Hochberg's
#' step-up procedure
#'
#' @details
#' The Benjamini-Hochberg (BH) procedure
#' (Sarkar, 2008; Remark 4.2)
#' consists in using the decision procedure
#' described in [mtp-package] using:
#'
#' * the adjustment factors:
#' \eqn{\qquad\quad
#'  \displaystyle{
#'   a_j = \frac{d}{j},
#'   \ \text{for}\ j=1, \ldots, d.
#'  }
#' }
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
#' Benjamini, Y., & Hochberg, Y. (1995). Controlling the
#' false discovery rate: a practical and powerful approach to
#' multiple testing. Journal of the Royal statistical society:
#' series B (Methodological), 57(1), 289-300.
#'
#' @export
frr_bh <- function(p_value, .return = "p", alpha = NULL) {
  
  # check arguments
  .check_p_value()
  .check_return()
  
  # get adjustment factors
  d <- length(p_value)
  j <- d:1L
  o <- order(p_value, decreasing = TRUE)
  ro <- order(o)
  a <- d / j
  
  # output
  p <- pmin(cummin(a * p_value[o]), 1)[ro]
  if (!is.null(alpha)) {
    return( (alpha * p_value) / p )
  } else {
    if (.return == "a") {
      return(p / p_value)
    }
  }
  p
}
