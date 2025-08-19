#' Title
#'
#' @param x x
#' @param ... x
#'
#' @return x
#' @export
summary.sRDA <- function(x, ...) {
    if (length(x) > 1) {
        for (ith in 1:length(x)) {

            cat(rep("-", 20), sep="")
            cat(paste0("\nComponent ", x[[ith]]$component))
            cat(paste0("\n1. Latent association : ", round(rh_test$estimate, 3)))
            cat(paste0("\n2. Absolute beta value: ", round(x[[ith]]$sum_abs_betas, 3)))
            cat("\n")

        }
    } else {
        cat("not yet finished")
    }

    return(invisible(x))
}
