#' FWER control using Hochberg (step-up)
#'
#' @inheritParams fwer_bon
#'
#' @description
#' Control \eqn{\text{FWER(k)}} using generalization of
#' Hochberg's step-up procedure
#'
#' @details
#' The generalized Holm procedure
#' (Sarkar, 2008; Remark 4.2)
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
#' The generalized Hochberg procedure uses the same
#' adjustment factors as the generalized Holm procedure
#' (see [fwer_holm()]) but yields a different
#' set of rejected hypotheses, as it is a step-up rather
#' than a step-down procedure.
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
#' @family FWER
#'
#' @references
#' Hochberg, Y. (1988). A sharper Bonferroni procedure for
#' multiple tests of significance. Biometrika, 75(4), 800-802.\cr
#' Sarkar, S. K. (2008). Generalizing Simes’ test and
#' Hochberg’s stepup procedure
#' The Annals of Statistics, 36(1), 337-363.
#'
#' @export
fwer_hoch <- function(p, k = 0, alpha = NULL, output = "p") {

  # check input
  m <- length(p)
  .check_p()
  .check_k()
  .check_alpha()
  .check_output()

  # get adjustment factors
  j <- m:1L
  o <- order(p)[j]
  ro <- order(o)
  a <- (m - pmax(0, j - k - 1)) / (k + 1)

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
