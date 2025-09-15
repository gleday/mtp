#' @rdname m0
#'
#' @export
pi0 <- function(p_values, lambda = 0.5) {
  m0(p_values = p_values, lambda = lambda) / length(p_values)
}
