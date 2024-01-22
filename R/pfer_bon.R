#' PFER control using Bonferroni (single-step)
#'
#' @inheritParams fwer_bon
#' @param gamma [numeric] scalar. Target \eqn{\text{PFER}}
#' level (between 1 and m). Overrides `output` and
#' return adjusted critical values.
#'
#' @description
#' Control \eqn{\text{PFER}} using Bonferroni's
#' single-step procedure
#'
#' @details
#' The Bonferroni procedure
#' consists in using the decision procedure
#' described in [mtp-package] with:
#'
#' * the adjustment factor:
#' \eqn{\qquad\quad
#'  \displaystyle{
#'   a = m
#'  }.
#' }
#'
#' * the adjusted P-values:
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
#' The Bonferroni procedure guarantees
#' that \eqn{\text{PFER} \leq \gamma} without
#' assumptions on the dependence of P-values.
#'
#' @importFrom "assertthat" "assert_that" "is.count" "is.number" "is.string"
#' @importFrom "dplyr" "between"
#' @importFrom "purrr" "map_lgl"
#'
#' @return
#' A [numeric] vector of:
#' * adjusted P-values when `output = "p"`,
#' * adjustment factors when `output = "a"`,
#' * adjusted critical values when `gamma` is provided.
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
#' The Annals of Applied Statistics, 1(1), 179-190.
#'
#' @export
pfer_bon <- function(p, gamma = NULL, output = "p") {

  # check arguments
  .check_p()
  .check_gamma()
  .check_output()

  # get adjustment factor
  m <- length(p)
  a <- m

  # output
  p_adj <- pmin(a * p, m)
  if (!is.null(gamma)) {
    return(rep(gamma / a, m))
  } else {
    if (output == "a") {
      return(rep(a, m))
    }
  }
  p_adj
}
