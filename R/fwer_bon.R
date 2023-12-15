#' FWER control using Bonferroni (single-step)
#'
#' @param p_value [numeric] vector. Observed P-values.
#' @param k [numeric] scalar. Exceedance level (positive integer).
#' @param .return [character] scalar. Return adjusted P-values
#' (`"p"`) or adjustment factors (`"a"`).
#' @param alpha [numeric] scalar. Target \eqn{\text{FWER(k)}}
#' level (between 0 and 1). Overrides `.return` and
#' return adjusted critical values.
#'
#' @description
#' Control \eqn{\text{FWER(k)}} using generalization of
#' Bonferroni's single-step procedure
#'
#' @details
#' The generalized Bonferroni procedure
#' (Lehmann & Romano, 2005; Theorem 2.1)
#' consists in using the decision procedure
#' described in [mtp-package] with:
#'
#' * the adjustment factor:
#' \eqn{\qquad\quad
#'  \displaystyle{
#'   a = \frac{m}{k}
#'  }.
#' }
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
#' The generalized Bonferroni procedure guarantees
#' that \eqn{\text{FWER(k)} \leq \alpha} without
#' assumptions on the dependence of P-values.
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
#' @family FWER
#'
#' @author Gwenaël G.R. Leday
#'
#' @references
#' Bonferroni, C. (1936). Teoria statistica delle classi e
#' calcolo delle probabilita. Pubblicazioni del R Istituto
#' Superiore di Scienze Economiche e Commericiali di Firenze, 8, 3-62.\cr
#' Lehmann, E. L., & Romano, J. P. (2005). Generalizations of
#' the familywise error rate. The Annals of Statistics, 33(3), 1138-1154.
#'
#' @export
fwer_bon <- function(p_value, k = 1, .return = "p", alpha = NULL) {

  # check arguments
  .check_p_value()
  .check_k()
  .check_return()

  # get adjustment factor
  m <- length(p_value)
  a <- m / k

  # output
  p <- pmin(a * p_value, 1)
  if (!is.null(alpha)) {
    return(rep(alpha / a, m))
  } else {
    if (.return == "a") {
      return(rep(a, m))
    }
  }
  p
}
