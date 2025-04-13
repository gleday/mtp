#' FWER control using Bonferroni (single-step)
#'
#' @param p [numeric] vector. Observed p-values.
#' @param k [numeric] scalar. Exceedance level (integer).
#' @param alpha [numeric] scalar. Target \eqn{\text{FWER(k)}}
#' level (between 0 and 1). Overrides `output` and
#' return adjusted critical values.
#' @param output [character] scalar. Return adjusted p-values
#' (`"p"`) or adjustment factors (`"a"`).
#'
#' @description
#' Control \eqn{\text{FWER(k)}} using generalization of
#' Bonferroni's single-step procedure
#'
#' @details
#' The generalized Bonferroni procedure
#' (Lehmann & Romano, 2005; Theorem 2.1)
#' yields:
#'
#' * the adjustment factor:
#' \eqn{\qquad\quad
#'  \displaystyle{
#'   a = \frac{m}{k + 1}
#'  }.
#' }
#'
#' * the adjusted p-values:
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
fwer_bon <- function(p, k = 0, alpha = NULL, output = "p") {

  # check input
  m <- length(p)
  .check_p()
  .check_k()
  .check_alpha()
  .check_output()

  # get adjustment factor
  a <- m / (k + 1)

  # output
  p_adj <- pmin(a * p, 1)
  if (!is.null(alpha)) {
    return(rep(alpha / a, m))
  } else {
    if (output == "a") {
      return(p_adj / p)
    }
  }
  p_adj
}
