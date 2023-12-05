#' FRX control using modified Hochberg (step-down)
#'
#' @inheritParams frx_holm
#'
#' @description
#' Control \eqn{\text{FRX(c)}} using modification of
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
#'   a_j = \frac{d - j + \left\lfloor c j \right\rfloor + 1}
#'   {\left\lfloor c j \right\rfloor + 1},
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
#' Here, \eqn{\left\lfloor x \right\rfloor} denotes
#' the floor function returning the
#' largest integer that is smaller or equal to \eqn{x}.
#'
#' The modified Hochberg procedure guarantees
#' that \eqn{\text{FRX(c)} \leq \alpha} under some
#' assumptions on the dependence of P-values.
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
#' * adjusted P-values when `.return = "p"`,
#' * adjustment factors when `.return = "a"`,
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
frx_hoch <- function(p_value, c = 0, .return = "p", alpha = NULL) {

  # check arguments
  .check_p_value()
  .check_c()
  .check_return()

  # get adjustment factors
  d <- length(p_value)
  j <- d:1L
  cj <- floor(c * j)
  o <- order(p_value)[j]
  ro <- order(o)
  a <- (d - j + cj + 1) / (cj + 1)

  # output
  p <- pmin(1, cummin(a * p_value[o]))[ro]
  if (!is.null(alpha)) {
    return((alpha * p_value) / p)
  } else {
    if (.return == "a") {
      return(p / p_value)
    }
  }
  p
}
