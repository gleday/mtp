#' FRX control using modified Holm (step-down)
#'
#' @param d [numeric] scalar. Exceedance level (between 0 and 1).
#' @param alpha [numeric] scalar. Target \eqn{\text{FRX(d)}}
#' level (between 0 and 1). Overrides `.return` and
#' return adjusted critical values.
#' @inheritParams fwer_bon
#'
#' @description
#' Control \eqn{\text{FRX(d)}} using modification of
#' Holm's step-down procedure (Lehmann & Romano, 2005)
#'
#' @details
#' The modified Holm procedure
#' (Lehmann & Romano, 2005; Theorem 3.1)
#' consists in using the decision procedure
#' described in [mtp-package] using:
#'
#' * the adjustment factors:
#' \eqn{\qquad\quad
#'  \displaystyle{
#'   a_j = \frac{m - j + \left\lfloor d j \right\rfloor + 1}
#'   {\left\lfloor d j \right\rfloor + 1},
#'   \ \text{for}\ j = 1, \ldots, m.
#'  }
#' }
#'
#' * the adjusted P-values:
#' \eqn{\qquad\quad
#' \displaystyle{
#'   \begin{cases}
#'   \widetilde{p}_{1} = \min\left( a_j p_{j}, 1\right),\\
#'   \widetilde{p}_{j} = \max\left( a_j p_{j},
#'   \widetilde{p}_{j - 1}\right), \ \text{for}\ j=2, \ldots, m
#'   \end{cases}
#'  }
#'  }
#'
#' * the adjusted critical values:
#' \eqn{\qquad
#' \displaystyle{
#'   \widetilde{\alpha}_{j} = \frac{\alpha}{a_j},
#'   \ \text{for}\ j = 1, \ldots, m.
#' } }
#'
#' Here, \eqn{\left\lfloor x \right\rfloor} denotes
#' the floor function returning the
#' largest integer that is smaller or equal to \eqn{x}.
#'
#' The modified Holm procedure guarantees
#' that \eqn{\text{FRX(d)} \leq \alpha} with some
#' assumptions on the dependence of P-values
#' (Simes inequality).
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
#' Lehmann, E. L., & Romano, J. P. (2005). Generalizations of
#' the familywise error rate. The Annals of Statistics, 33(3), 1138-1154.\cr
#' Guo, W., He, L., & Sarkar, S. K. (2014). Further results on
#' controlling the false discovery proportion.
#'
#' @export
frx_holm <- function(p_value, d = 0, .return = "p", alpha = NULL) {

  # check arguments
  .check_p_value()
  .check_d()
  .check_return()

  # get adjustment factors
  m <- length(p_value)
  j <- seq_len(m)
  dj <- floor(d * j)
  o <- order(p_value)
  ro <- order(o)
  a <- (m - j + dj + 1) / (dj + 1)

  # output
  p <- pmin(1, cummax(a * p_value[o]))[ro]
  if (!is.null(alpha)) {
    return((alpha * p_value) / p)
  } else {
    if (.return == "a") {
      return(p / p_value)
    }
  }
  p
}
