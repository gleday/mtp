#' @rdname d0
#'
#' @export
pi0 <- function(p_value, lambda = 0.5) {
  d0(p_value = p_value, lambda = lambda) / length(p_value)
}
