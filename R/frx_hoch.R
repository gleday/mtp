#' FRX control using modified Hochberg (step-down)
#'
#' @inheritParams frx_holm
#'
#' @description
#' Control \eqn{\text{FRX(d)}} using modification of
#' Hochberg's step-up procedure (Sarkar, 2008;
#' Guo et al., 2014)
#'
#' @details
#' The modified Holm procedure
#' (Guo et al., 2014; Theorem 3.1)
#' consists in using the decision procedure
#' described in [mtp-package] using:
#'
#' * the adjustment factors:
#' \eqn{\qquad\quad
#'  \displaystyle{
#'   a_j = \frac{d - j + \left\lfloor d j \right\rfloor + 1}
#'   {\left\lfloor d j \right\rfloor + 1},
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
#' Here, \eqn{\left\lfloor x \right\rfloor} denotes
#' the floor function returning the
#' largest integer that is smaller or equal to \eqn{x}.
#'
#' The modified Hochberg procedure guarantees
#' that \eqn{\text{FRX(d)} \leq \alpha} under some
#' assumptions on the dependence of P-values
#' (Simes inequality).
#'
#' The modified Hochberg procedure uses the same
#' adjustment factors as the modified Holm procedure
#' (see [frx_holm()]) but yields a different
#' set of rejected hypotheses.
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
#' @family FRX
#'
#' @references
#' Sarkar, S. K. (2008). Generalizing Simes’ test and
#' Hochberg’s stepup procedure.
#' The Annals of Statistics, 36(1), 337-363.\cr
#' Guo, W., He, L., & Sarkar, S. K. (2014). Further results on
#' controlling the false discovery proportion.
#'
#' @export
frx_hoch <- function(p, d = 0, alpha = NULL, output = "p") {

  # check arguments
  .check_p()
  .check_d()
  .check_alpha()
  .check_output()

  # get adjustment factors
  m <- length(p)
  j <- m:1L
  dj <- floor(d * j)
  o <- order(p)[j]
  ro <- order(o)
  a <- (m - j + dj + 1) / (dj + 1)

  # output
  p_adj <- pmin(1, cummin(a * p[o]))[ro]
  if (!is.null(alpha)) {
    return((alpha * p) / p_adj)
  } else {
    if (output == "a") {
      return(p_adj / p)
    }
  }
  p_adj
}
