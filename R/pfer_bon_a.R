#' PFER control using adaptive Bonferroni (single-step)
#'
#' @inheritParams pfer_bon
#' @inheritParams fwer_bon_a
#'
#' @description
#' Control PFER using adaptive Bonferroni's
#' single-step procedure
#'
#' @details
#' The Bonferroni procedure
#' (Meng et al, 2014; Theorem 2.1)
#' yields:
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
#' @inherit fwer_bon return
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
pfer_bon_a <- function(
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

  # get adjustment factor
  m <- length(p_values)
  a <- m0(p_values = p_values, lambda = lambda)

  # output
  adj_p_values <- pmin(a * p_values, m)
  switch(output,
    p_values = adj_p_values,
    c_values = rep(c_value / a, m),
    decisions = (adj_p_values <= c_value) * 1L,
    factors = rep(a, m)
  )
}
