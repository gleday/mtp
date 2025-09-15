#' FRX control using modified Holm (step-down)
#'
#' @param d \[ `numeric(1)` \]\cr
#' A number between 0 and 1 specifying the exceedance level.
#' @param c_value \[ `numeric(1)` \]\cr
#' A number between 0 and 1 specifying the critical value for FRX(d).\cr
#' Required if `output = "c_values"` or `output = "decisions"`.
#' @inheritParams fwer_bon
#'
#' @description
#' Control FRX(d) using modification of
#' Holm's step-down procedure (Lehmann & Romano, 2005)
#'
#' @details
#' Denote by \eqn{p_{(1)} \leq \ldots \leq p_{(m)}}
#' the ordered p-values.
#'
#' The modified Holm procedure
#' (Lehmann & Romano, 2005; Theorem 3.1)
#' yields:
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
#' * the adjusted p-values:
#' \eqn{\qquad\quad
#'  \displaystyle{
#'   \begin{cases}
#'   \widetilde{p}_{(1)} = \min\left( a_j p_{(1)}, 1\right),\\
#'   \widetilde{p}_{(j)} = \max\left( a_j p_{(j)},
#'   \widetilde{p}_{(j - 1)}\right), \ \text{for}\ j = 2, \ldots, m
#'   \end{cases}
#'  }
#' }
#'
#' * the adjusted critical values:
#' \eqn{\qquad
#'  \displaystyle{
#'   \widetilde{\alpha}_{(j)} = \frac{\alpha}{a_j},
#'   \ \text{for}\ j = 1, \ldots, m.
#'  }
#' }
#'
#' Here, \eqn{p_{(1)} \leq \ldots \leq p_{(m)}}
#' denote the ordered p-values with
#' \eqn{\widetilde{p}_{(1)}, \ldots,
#' \widetilde{p}_{(m)}} and
#' \eqn{\widetilde{\alpha}_{(1)}, \ldots,
#' \widetilde{\alpha}_{(m)}}
#' their corresponding adjusted p-values and
#' critical values. Also,
#' \eqn{\left\lfloor x \right\rfloor} denotes
#' the floor function returning the
#' largest integer that is smaller or equal to \eqn{x}.
#'
#' @inherit fwer_bon return
#'
#' @family FRX
#'
#' @references
#' Lehmann, E. L., & Romano, J. P. (2005). Generalizations of
#' the familywise error rate. The Annals of Statistics, 33(3), 1138-1154.\cr
#' Guo, W., He, L., & Sarkar, S. K. (2014). Further results on
#' controlling the false discovery proportion.
#'
#' @examples
#' # simulate p-values
#' set.seed(123)
#' sim <- simulate_p_values(m = 1000, m0 = 900, s1 = 0.1, s2 = 10)
#'
#' # 1. control FRX(d = 0.1) at 0.05 using adjusted p-values
#' adj_pv <- frx_holm(sim$p_value, d = 0.1)
#' sum(adj_pv <= 0.05)
#'
#' # 2. control FRX(d = 0.1) at 0.05 using adjusted critical values
#' adj_cv <- frx_holm(
#'   sim$p_value,
#'   d = 0.1,
#'   c_value = 0.05,
#'   output = "c_values"
#' )
#' sum(sim$p_value <= adj_cv)
#'
#' # 3. control FRX(d = 0.1) at 0.05 using decision indicator
#' ind <- frx_holm(
#'   sim$p_value,
#'   d = 0.1,
#'   c_value = 0.05,
#'   output = "decisions"
#' )
#' sum(ind)
#'
#' # 4. adjustment factors to control FRX(d = 0.1) at 0.05
#' adj_fct <- frx_holm(sim$p_value, d = 0.1, output = "factors")
#'
#' @export
frx_holm <- function(
    p_values,
    d = 0,
    c_value = NULL,
    output = c("p_values", "c_values", "decisions", "factors")) {
  # check input
  .check_p_values()
  .check_d()
  output <- arg_match(
    output,
    c("p_values", "c_values", "decisions", "factors")
  )
  .check_c_value()

  # get adjustment factors
  m <- length(p_values)
  j <- seq_len(m)
  dj <- floor(d * j)
  o <- order(p_values)
  ro <- order(o)
  a <- (m - j + dj + 1) / (dj + 1)

  # output
  adj_p_values <- pmin(1, cummax(a * p_values[o]))[ro]
  switch(output,
    p_values = adj_p_values,
    c_values = (c_value * p_values) / adj_p_values,
    decisions = (adj_p_values <= c_value) * 1L,
    factors = a[ro]
  )
}
