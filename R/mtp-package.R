#' mtp: multiple hypothesis testing procedures
#'
#' @description
#' Control Type I errors (false positives) when
#' conducting multiple hypothesis tests using
#' a variety of error criteria and adjustment methods.
#'
#' Retrieve adjusted p-values, adjusted critical values,
#' adjustment factors, or rejection decisions.
#'
#' @section Control the (generalized) familywise error rate:
#' * [fwer_bon()]
#' * [fwer_holm()]
#' * [fwer_hoch()]
#' * [fwer_bon_a()]
#' * [fwer_bon_w()]
#'
#' @section Control the false rejection exceedance:
#' * [frx_holm()]
#' * [frx_hoch()]
#'
#' @section Control the false rejection rate:
#' * [frr_bh()]
#' * [frr_bh_a()]
#'
#' @section Control the per family error rate:
#' * [pfer_bon()]
#' * [pfer_bon_a()]
#'
#' @section Estimate the number/proportion of true null hypotheses:
#' * [m0()]
#' * [pi0()]
#'
#' @section Simulate p-values:
#' * [simulate_p_values()]
#'
#' @docType package
#' @name mtp-package
#' @aliases mtp-package
#' @aliases mtp
#'
#' @importFrom "assertthat" "assert_that" "is.count" "is.number"
#' @importFrom "assertthat" "is.string" "are_equal"
#' @importFrom "purrr" "map_lgl" "map_dbl"
#' @importFrom "stats" "rbeta"
#' @importFrom "tibble" "tibble"
#' @importFrom "rlang" "arg_match"
#'
#' @keywords internal
"_PACKAGE"
