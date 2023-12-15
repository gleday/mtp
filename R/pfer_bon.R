#' PFER control using Bonferroni (single-step)
#'
#' @inheritParams fwer_bon
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
#'   \widetilde{\alpha} = \frac{\alpha}{a}.
#'  }
#' }
#'
#' The Bonferroni procedure guarantees
#' that \eqn{\text{PFER} \leq \alpha} without
#' assumptions on the dependence of P-values.
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
pfer_bon <- function(p_value, .return = "p", alpha = NULL) {

  fwer_bon(
    p_value = p_value,
    k = 1,
    .return = .return,
    alpha = alpha
    )

}
