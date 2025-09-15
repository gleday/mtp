#' FWER control using Holm (step-down)
#'
#' @inheritParams fwer_bon
#'
#' @description
#' Control FWER(k) using generalization of
#' Holm's step-down procedure (also called the
#' Bonferroni-Holm or the sequentially rejective
#' Bonferroni procedure).
#'
#' @details
#' The generalized Holm procedure
#' (Lehmann & Romano, 2005; Theorem 2.2)
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
#'  \displaystyle{
#'   \begin{cases}
#'   \widetilde{p}_{(1)} = \min\left( a_1 p_{(1)}, 1\right),\\
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
#' @family FWER
#'
#' @references
#' Holm, S. (1979). A simple sequentially rejective
#' multiple test procedure. Scandinavian journal of
#' statistics, 65-70.\cr
#' Lehmann, E. L., & Romano, J. P. (2005).
#' Generalizations of the familywise error rate.
#' The Annals of Statistics, 33(3), 1138-1154.
#'
#' @examples
#' # simulate p-values
#' set.seed(123)
#' sim <- simulate_p_values(m = 1000, m0 = 900, s1 = 0.1, s2 = 10)
#'
#' # 1. control FWER(k = 1) at 0.05 using adjusted p-values
#' adj_pv <- fwer_holm(sim$p_value, k = 1)
#' sum(adj_pv <= 0.05)
#'
#' # 2. control FWER(k = 1) at 0.05 using adjusted critical values
#' adj_cv <- fwer_holm(
#'   sim$p_value,
#'   k = 1,
#'   c_value = 0.05,
#'   output = "c_values"
#' )
#' sum(sim$p_value <= adj_cv)
#'
#' # 3. control FWER(k = 1) at 0.05 using decision indicator
#' ind <- fwer_holm(
#'   sim$p_value,
#'   k = 1,
#'   c_value = 0.05,
#'   output = "decisions"
#' )
#' sum(ind)
#'
#' # 4. adjustment factors to control FWER(k = 1) at 0.05
#' adj_fct <- fwer_holm(sim$p_value, k = 1, output = "factors")
#'
#' @export
fwer_holm <- function(
    p_values,
    k = 0,
    c_value = NULL,
    output = c("p_values", "c_values", "decisions", "factors")) {
  # check input
  m <- length(p_values)
  .check_p_values()
  .check_k()
  output <- arg_match(
    output,
    c("p_values", "c_values", "decisions", "factors")
  )
  .check_c_value()

  # get adjustment factors
  j <- seq_len(m)
  o <- order(p_values)
  ro <- order(o)
  a <- (m - pmax(0, j - k - 1)) / (k + 1)

  # output
  adj_p_values <- pmin(cummax(a * p_values[o]), 1)[ro]
  switch(output,
    p_values = adj_p_values,
    c_values = (c_value * p_values) / adj_p_values,
    decisions = (adj_p_values <= c_value) * 1L,
    factors = a[ro]
  )
}
