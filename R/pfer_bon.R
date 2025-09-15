#' PFER control using Bonferroni (single-step)
#'
#' @inheritParams fwer_bon
#' @param c_value \[ `numeric(1)` \]\cr
#'  A non-negative integer smaller than `length(p_values)`
#'  specifying the critical value for FRR.\cr
#'  Required if `output = "c_values"` or `output = "decisions"`.
#'
#' @description
#' Control PFER using Bonferroni's
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
#'   a = m
#'  }.
#' }
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
#' Tukey, J. W. (1953). The problem of multiple comparison.
#' Unpublished manuscript. In The Collected Works of John W.
#' Tukey VIII. Multiple Comparisons: 1948–1983 1–300. Chapman
#' and Hall, New York.\cr
#' Gordon, A., Glazko, G., Qiu, X., & Yakovlev, A. (2007).
#' Control of the mean number of false discoveries,
#' Bonferroni and stability of multiple testing.
#' The Annals of Applied Statistics, 1(1), 179-190.\cr
#' Meng, X., Wang, J. & Wu, X. (2014) Multiple Comparisons
#' Controlling Expected Number of False Discoveries,
#' Communications in Statistics - Theory and Methods,
#' 43(13), 2830–2843.
#'
#' @export
pfer_bon <- function(
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

  # get adjustment factor
  m <- length(p_values)
  a <- m

  # output
  adj_p_values <- pmin(a * p_values, m)
  switch(output,
    p_values = adj_p_values,
    c_values = rep(c_value / a, m),
    decisions = (adj_p_values <= c_value) * 1L,
    factors = rep(a, m)
  )
}
