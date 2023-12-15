#' mtp: multiple hypothesis testing with error control
#'
#' @description
#' `mtp` offers functions to apply various forms of
#' statistical control (using different criteria and
#' methods) on the occurrence of Type I errors
#' (incorrect rejections) when
#' conducting multiple statistical tests.
#'
#' @section Setting:
#'
#' Testing of \eqn{m} null hypotheses
#' \eqn{H_{1}, \ldots, H_{m}} with
#' observed P-values \eqn{p_{1}, \ldots, p_{m}}.
#'
#' @section Notation:
#'
#' \eqn{p_{(1)} \leq \ldots \leq p_{(m)}}: ordered P-values.\cr
#' \eqn{H_{(1)}, \ldots, H_{(m)}}: ordered null hypotheses.\cr
#' \eqn{\text{R}}: number of rejected null hypotheses.\cr
#' \eqn{\text{V}}: number of false rejections.
#'
#' @section Type I error criteria:
#'
#' `mtp` allows to use the following criteria
#' to assess the likelihood of committing Type I errors:
#'
#' * the \strong{familywise error rate} (FWER),
#' defined as the probability that the number of
#' false rejections equals or exceeds a predefined
#' number \eqn{1 \leq \text{k} \leq \text{m}}:
#' \deqn{\text{FWER(k)} = \text{Pr}(\text{V} \geq \text{k}).}
#' See functions [fwer_bon()], [fwer_holm()],
#' [fwer_hoch()] and [fwer_abon()].
#' \cr\cr
#'
#' * the \strong{false rejection exceedance} (FRX),
#' defined as the probability that the proportion
#' of false rejections among all rejected hypotheses
#' exceeds a predefined proportion \eqn{0 \leq \text{d} \leq 1}:
#' \deqn{\text{FRX(d)} = \text{Pr}(\text{V}/\text{R} > \text{d}).}
#' See functions [frx_holm()] and [frx_hoch()].
#'
#' * the \strong{per family error rate} (PFER),
#' defined as the expected number of false rejections:
#' \deqn{\text{PFER} = \text{E}(\text{V}).}
#' See functions [pfer_bon()] and [pfer_abon()].
#'
#' * the \strong{false rejection rate} (FRR),
#' defined as the expected proportion
#' of false rejections among all rejected hypotheses:
#' \deqn{\text{FRR} = \text{E}(\text{V}/\text{R}).}
#' See functions [frr_bh()] and [frr_abh()].
#'
#' Note: the term "discovery" commonly
#' adopted in literature is replaced by the
#' more neutral term "rejection".
#'
#' @section Adjustment:
#'
#' Multiple testing procedures cap
#' the probability of producing Type I errors
#' (in the sense of the chosen criterion)
#' by adjusting the significance level of
#' each test.
#' The adjustment consists in modifying
#' the observed P-value upward or, equivalently,
#' the critical value downward with
#' an \strong{adjustment factor}.
#'
#' The functions in `mtp` perform these adjustments
#' and allow to obtain:
#' * the \strong{adjusted P-values}
#' \eqn{\widetilde{p}_{1}, \ldots, \widetilde{p}_{m}}.
#' when the argument `.return = "p"`.
#' * the \strong{adjustment factors} \eqn{a_{1}, \ldots, a_{m}}
#' when the argument `.return = "a"`.
#' * the \strong{adjusted critical values}
#' \eqn{\widetilde{\alpha}_{1}, \ldots, \widetilde{\alpha}_{m}}
#' when the argument `alpha` is provided.
#'
#' Adjustments in `mtp` preserve the ordering of P-values
#' and ensure that the adjusted P-values and critical values
#' are comprised between 0 and 1.
#'
#' @section Decision procedure:
#'
#' After adjustment, the decision procedure
#' consists in applying one of the following
#' (equivalent) rules:
#'
#' * reject hypotheses with \strong{adjusted P-values} smaller than or
#' equal to the critical value \eqn{\alpha}, i.e.:
#' \deqn{
#'  \displaystyle{
#'   \text{reject}\ \ H_{(j)} \qquad
#'   \text{if} \qquad \widetilde{p}_{(j)} \leq \alpha,
#'   \qquad\text{for}\ j=1, \ldots, m.
#'  }
#' }
#'
#' * reject hypotheses with observed P-values smaller than or equal to
#' the \strong{adjusted critical values}, i.e.:
#' \eqn{\widetilde{\alpha}_{1}, \ldots, \widetilde{\alpha}_{m}}:
#' \deqn{
#'  \displaystyle{
#'  \text{reject}\ \ H_{(j)} \qquad
#'   \text{if} \qquad p_{(j)} \leq \widetilde{\alpha}_{j},
#'   \qquad\text{for}\ j=1, \ldots, m.
#'  }
#' }
#'
#' Some adjustments produce
#' a different adjusted critical value for
#' each null hypothesis, while others
#' yield a single adjusted critical value
#' for all hypotheses
#' (\eqn{\widetilde{\alpha}_{1} = \ldots =
#' \widetilde{\alpha}_{m}}).
#'
#' @section Types of multiple testing procedures:
#'
#' `mtp` implements the following types of
#' multiple testing procedures:
#'
#' * \strong{one-step} procedures that adjust each test
#' regardless of the outcome of other tests.\cr
#' See for example [fwer_bon()] and [fwer_abon()].
#'
#' * \strong{step-down} procedures that adjust each test
#' sequentially, starting with the most significant,
#' (with smallest P-value).\cr
#' See for example [fwer_holm()] and [frx_holm()].
#'
#' * \strong{step-up} procedures that adjust each test
#' sequentially, starting with the least significant,
#' (with largest P-value).\cr
#' See for example [fwer_hoch()] and [frx_hoch()].
#'
#' Note that the code used for the step-down and
#' step-up procedures in `mtp` is largely based
#' on the code of `p.adjust()` from package `stat`.
#'
#' @docType package
#' @name mtp-package
#' @aliases mtp-package
#' @aliases mtp
#'
#' @references
#' Bonferroni, C. (1936). Teoria statistica delle
#' classi e calcolo delle probabilita. Pubblicazioni del R.
#' Istituto Superiore di Scienze Economiche e
#' Commericiali di Firenze, 8, 3-62.\cr
#' Holm, S. (1979). A simple sequentially rejective
#' multiple test procedure.
#' Scandinavian journal of statistics, 65-70.\cr
#' Hochberg, Y. (1988). A sharper Bonferroni procedure for
#' multiple tests of significance. Biometrika, 75(4), 800-802.\cr
#' Lehmann, E. L., & Romano, J. P. (2005). Generalizations of
#' the familywise error rate. The Annals of Statistics, 33(3), 1138-1154.\cr
#' Gordon, A., Glazko, G., Qiu, X., & Yakovlev, A. (2007).
#' Control of the mean number of false discoveries,
#' Bonferroni and stability of multiple testing.
#' The Annals of Applied Statistics, 1(1), 179-190.\cr
#' Sarkar, S. K. (2008). Generalizing Simes’ test and
#' Hochberg’s stepup procedure. The Annals of Statistics, 36(1), 337-363.\cr
#' Guo, W., He, L., & Sarkar, S. K. (2014). Further results on
#' controlling the false discovery proportion.\cr
#' Wang, L. (2017). Adaptive procedure for generalized
#' familywise error rate control. Communications in
#' Statistics-Simulation and Computation, 46(10), 8140-8151.
#'
#' @keywords internal
"_PACKAGE"
