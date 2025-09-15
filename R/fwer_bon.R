#' FWER control using Bonferroni (single-step)
#'
#' @param p_values \[ `numeric()` \]\cr
#' A numeric vector of observed p-values.
#' @param k \[ `integer(1)` \]\cr
#' A non-negative integer smaller than `length(p_values)`
#' specifying the exceedance level.
#' @param c_value \[ `numeric(1)` \]\cr
#' A number between 0 and 1 specifying the critical value for FWER(k).\cr
#' Required if `output = "c_values"` or `output = "decisions"`.
#' @param output \[ `character(1)` \]\cr
#' The type of value returned by the function:
#' * `"p_values"` returns the adjusted p-values
#' * `"c_values"` returns the adjusted critical values
#' * `"decisions"` returns the rejection decisions
#' (1 = rejected and 0 = non-rejected)
#' * `"factors"` returns adjustment factors
#'
#' @description
#' Control FWER(k) using generalization of
#' Bonferroni's single-step procedure.
#'
#' @details
#' The generalized Bonferroni procedure
#' (Lehmann & Romano, 2005; Theorem 2.1)
#' yields:
#'
#' * adjustment factor:
#' \eqn{\qquad\quad
#'  \displaystyle{
#'   a = \frac{m}{k + 1}
#'  }.
#' }
#'
#' * adjusted p-values:
#' \eqn{\qquad\quad
#'  \displaystyle{
#'   \widetilde{p}_{j} = \min\left( a\times p_{j}, 1\right),
#'   \ \text{for}\ j=1, \ldots, m.
#'  }
#' }
#'
#' * adjusted critical value:
#' \eqn{\qquad
#'  \displaystyle{
#'   \widetilde{\alpha} = \frac{\alpha}{a}.
#'  }
#' }
#'
#' @return
#' \[ `numeric()` \] \cr \cr
#' A numeric vector the same length as `p_values` with:
#' * the adjusted p-values if `output = "p_values"`
#' * the adjusted critical values if `output = "c_values"`
#' * the rejection decisions (1 = rejected and 0 = non-rejected)
#' if `output = "decisions"`
#' * the adjustment factors if `output = "factors"`
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
#' @examples
#' # simulate p-values
#' set.seed(123)
#' sim <- simulate_p_values(m = 1000, m0 = 900, s1 = 0.1, s2 = 10)
#'
#' # 1. control FWER(k = 1) at 0.05 using adjusted p-values
#' adj_pv <- fwer_bon(sim$p_value, k = 1)
#' sum(adj_pv <= 0.05)
#'
#' # 2. control FWER(k = 1) at 0.05 using adjusted critical values
#' adj_cv <- fwer_bon(
#'   sim$p_value,
#'   k = 1,
#'   c_value = 0.05,
#'   output = "c_values"
#' )
#' sum(sim$p_value <= adj_cv)
#'
#' # 3. control FWER(k = 1) at 0.05 using decision indicator
#' ind <- fwer_bon(
#'   sim$p_value,
#'   k = 1,
#'   c_value = 0.05,
#'   output = "decisions"
#' )
#' sum(ind)
#'
#' # 4. adjustment factors to control FWER(k = 1) at 0.05
#' adj_fct <- fwer_bon(sim$p_value, k = 1, output = "factors")
#'
#' @export
fwer_bon <- function(
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

  # get adjustment factor
  a <- m / (k + 1)

  # output
  adj_p_values <- pmin(a * p_values, 1)
  switch(output,
    p_values = adj_p_values,
    c_values = rep(c_value / a, m),
    decisions = (adj_p_values <= c_value) * 1L,
    factors = rep(a, m)
  )
}
