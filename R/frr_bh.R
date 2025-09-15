#' FRR control using Benjamini-Hochberg (step-up)
#'
#' @inheritParams fwer_bon
#' @param c_value \[ `numeric(1)` \]\cr
#' A number between 0 and 1 specifying the critical value for FRR.\cr
#' Required if `output = "c_values"` or `output = "decisions"`.
#'
#' @description
#' Control FRR using Benjamini-Hochberg's
#' step-up procedure
#'
#' @details
#' The Benjamini-Hochberg (BH) procedure
#' (Benjamini and Hochberg, 1995)
#' yields:
#'
#' * the adjustment factors:
#' \eqn{\qquad\quad
#'  \displaystyle{
#'   a_j = \frac{m}{j},
#'   \ \text{for}\ j=1, \ldots, m.
#'  }
#' }
#'
#' * the adjusted p-values:
#' \eqn{\qquad\quad
#'  \displaystyle{
#'   \begin{cases}
#'   \widetilde{p}_{(m)} = \min\left( a_m p_{(m)}, 1\right),\\
#'   \widetilde{p}_{(j)} = \min\left( a_j p_{(j)},
#'   \widetilde{p}_{(j + 1)}\right), \ \text{for}\ j = m - 1, \ldots, 1
#'   \end{cases}
#'  }
#' }
#'
#' * the adjusted critical values:
#' \eqn{\qquad
#'  \displaystyle{
#'   \widetilde{\alpha}_{(j)} = \frac{\alpha}{a_j},
#'   \ \text{for}\ j=1, \ldots, m.
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
#' critical values.
#'
#' @inherit fwer_bon return
#'
#' @family FRR
#'
#' @references
#' Benjamini, Y., & Hochberg, Y. (1995). Controlling the
#' false discovery rate: a practical and powerful approach to
#' multiple testing. Journal of the Royal statistical society:
#' series B (Methodological), 57(1), 289-300.
#'
#' @examples
#' # simulate p-values
#' set.seed(123)
#' sim <- simulate_p_values(m = 1000, m0 = 900, s1 = 0.1, s2 = 10)
#'
#' # 1. control FRR at 0.1 using adjusted p-values
#' adj_pv <- frr_bh(sim$p_value)
#' sum(adj_pv <= 0.01)
#'
#' # 2. control FRR at 0.1 using adjusted critical values
#' adj_cv <- frr_bh(
#'   sim$p_value,
#'   c_value = 0.01,
#'   output = "c_values"
#' )
#' sum(sim$p_value <= adj_cv)
#'
#' # 3. control FRR at 0.05 using decision indicator
#' ind <- frr_bh(
#'   sim$p_value,
#'   c_value = 0.01,
#'   output = "decisions"
#' )
#' sum(ind)
#'
#' # 4. adjustment factors to control FRX(d = 0.1) at 0.05
#' adj_fct <- frr_bh(sim$p_value, output = "factors")
#'
#' @export
frr_bh <- function(
    p_values,
    c_value = NULL,
    output = c("p_values", "c_values", "decisions", "factors")) {
  # check arguments
  .check_p_values()
  output <- arg_match(
    output,
    c("p_values", "c_values", "decisions", "factors")
  )
  .check_c_value()

  # get adjustment factors
  m <- length(p_values)
  j <- m:1L
  o <- order(p_values, decreasing = TRUE)
  ro <- order(o)
  a <- m / j

  # output
  adj_p_values <- pmin(cummin(a * p_values[o]), 1)[ro]
  switch(output,
    p_values = adj_p_values,
    c_values = (c_value * p_values) / adj_p_values,
    decisions = (adj_p_values <= c_value) * 1L,
    factors = a[ro]
  )
}
