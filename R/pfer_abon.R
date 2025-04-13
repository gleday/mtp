#' PFER control using adaptive Bonferroni (single-step)
#'
#' @inheritParams pfer_bon
#' @inheritParams m0
#'
#' @description
#' Control \eqn{\text{PFER}} using adaptive Bonferroni's
#' single-step procedure
#'
#' @details
#' The adaptive Bonferroni procedure
#' consists in using the decision procedure
#' described in [mtp-package] with:
#'
#' * the adjustment factor:
#' \eqn{\qquad\quad
#'  \displaystyle{
#'   a = \widehat{m}_0
#'  },
#' }
#' where \eqn{\widehat{m}_0} is Storey's estimator of the
#' number \eqn{m_0} of true null hypotheses (see [m0()]).
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
#'   \widetilde{\gamma} = \frac{\gamma}{a}.
#'  }
#' }
#'
#' The adaptive Bonferroni procedure guarantees
#' that \eqn{\text{PFER} \leq \gamma} under
#' independence of p-values.
#'
#' @importFrom "assertthat" "assert_that" "is.count" "is.number" "is.string"
#' @importFrom "dplyr" "between"
#' @importFrom "purrr" "map_lgl"
#'
#' @return
#' A [numeric] vector of:
#' * adjusted p-values when `output = "p"`,
#' * adjustment factors when `output = "a"`,
#' * adjusted critical values when `gamma` is provided.
#'
#' @family PFER
#'
#' @author Gwenaël G.R. Leday
#'
#' @references
#' Guo, W., Lynch, G., & Romano, J. P. (2018).
#' A new approach for large scale multiple testing with
#' application to FDR control for graphically
#' structured hypotheses. arXiv preprint arXiv:1812.00258.
#'
#' @export
pfer_abon <- function(p, lambda = 0.5, gamma = NULL, output = "p") {

  # check arguments
  .check_p()
  .check_lambda()
  .check_gamma()
  .check_output()

  # get adjustment factor
  m <- length(p)
  a <- m0(p = p, lambda = lambda)

  # output
  p_adj <- pmin(a * p, m)
  if (!is.null(gamma)) {
    return(rep(gamma / a, m))
  } else {
    if (output == "a") {
      return(rep(a, m))
    }
  }
  p_adj
}
