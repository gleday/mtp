#' FRR control using adaptive Benjamini-Hochberg (step-up)
#'
#' @inheritParams frr_bh
#' @inheritParams fwer_bon_a
#'
#' @description
#' Control FRR using adaptive
#' Benjamini-Hochberg's step-up procedure and
#' Storey's plugin estimator
#'
#' @details
#' Storey's adaptive version of
#' Benjamini-Hochberg (BH) procedure
#' yields:
#'
#' * the adjustment factors:
#' \eqn{\qquad\quad
#'  \displaystyle{
#'   a_j = \frac{\hat{m}_0}{j},
#'   \ \text{for}\ j=1, \ldots, m.
#'  },
#' }\cr
#' where \eqn{\widehat{m}_0} is Storey's estimator of the
#' number \eqn{m_0} of true null hypotheses (see [m0()]).
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
#' Storey, J. D. (2002). A direct approach to false
#' discovery rates. Journal of the Royal Statistical
#' Society Series B: Statistical Methodology, 64(3), 479-498.\cr
#' Storey, J. D., Taylor, J. E., & Siegmund, D. (2004).
#' Strong control, conservative point estimation and
#' simultaneous conservative consistency of false
#' discovery rates: a unified approach.
#' Journal of the Royal Statistical Society Series B:
#' Statistical Methodology, 66(1), 187-205.
#'
#' @examples
#' # simulate p-values
#' set.seed(123)
#' sim <- simulate_p_values(m = 1000, m0 = 900, s1 = 0.1, s2 = 10)
#'
#' # 1. control FRR at 0.1 using adjusted p-values
#' adj_pv <- frr_bh_a(sim$p_value)
#' sum(adj_pv <= 0.01)
#'
#' # 2. control FRR at 0.1 using adjusted critical values
#' adj_cv <- frr_bh_a(
#'   sim$p_value,
#'   c_value = 0.01,
#'   output = "c_values"
#' )
#' sum(sim$p_value <= adj_cv)
#'
#' # 3. control FRR at 0.05 using decision indicator
#' ind <- frr_bh_a(
#'   sim$p_value,
#'   c_value = 0.01,
#'   output = "decisions"
#' )
#' sum(ind)
#'
#' # 4. adjustment factors to control FRX(d = 0.1) at 0.05
#' adj_fct <- frr_bh_a(sim$p_value, output = "factors")
#'
#' @export
frr_bh_a <- function(
    p_values,
    lambda = 0.5,
    c_value = NULL,
    output = c("p_values", "c_values", "decisions", "factors")) {
  # check arguments
  .check_p_values()
  .check_lambda()
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
  a <- m0(p_values = p_values, lambda = lambda) / j

  # output
  adj_p_values <- pmin(cummin(a * p_values[o]), 1)[ro]
  switch(output,
    p_values = adj_p_values,
    c_values = (c_value * p_values) / adj_p_values,
    decisions = (adj_p_values <= c_value) * 1L,
    factors = a[ro]
  )
}
