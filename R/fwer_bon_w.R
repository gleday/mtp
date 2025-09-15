#' FWER control using weighted Bonferroni (single-step)
#'
#' @inheritParams fwer_bon
#' @param weights \[ `numeric()` \]\cr
#' A numeric vector of weights between 0 and 1
#' to be allocated to the hypotheses.
#'
#' @description
#' Control FWER(k) using the weighted generalized
#' Bonferroni procedure of Romano and Wolf (2010).
#'
#' @details
#' The weighted generalized Bonferroni procedure
#' (Romano & Romano, 2010; Theorem 6.1)
#' yields:
#'
#' * adjustment factor:
#' \eqn{\qquad\quad
#'  \displaystyle{
#'   a_j = \frac{1}{w_j (k + 1)}
#'  }.
#' }
#'
#' * adjusted p-values:
#' \eqn{\qquad\quad
#'  \displaystyle{
#'   \widetilde{p}_{j} = \min\left( a_j\times p_{j}, 1\right),
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
#' @inherit fwer_bon return
#'
#' @family FWER weighted
#'
#' @author Gwenaël G.R. Leday
#'
#' @references
#' Romano, J. P., & Wolf, M. (2010). Balanced control of
#' generalized error rates. The Annals of Statistics, 38(1), 598-633.
#'
#' @export
fwer_bon_w <- function(
    p_values,
    k = 0,
    weights = NULL,
    c_value = NULL,
    output = c("p_values", "c_values", "decisions", "factors")) {
  # check input
  m <- length(p_values)
  .check_p_values()
  .check_k()
  .check_weights()
  output <- arg_match(
    output,
    c("p_values", "c_values", "decisions", "factors")
  )
  .check_c_value()

  # get adjustment factor
  a <- 1 / (weights * (k + 1))

  # output
  adj_p_values <- pmin(a * p_values, 1)
  switch(output,
    p_values = adj_p_values,
    c_values = rep(c_value / a, m),
    decisions = (adj_p_values <= c_value) * 1,
    factors = rep(a, m)
  )
}
