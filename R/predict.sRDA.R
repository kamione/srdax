#' S3 methods predict
#' predict correlations of latent dimensions
#'
#' @param rda x
#' @param exploratory_new y
#' @param response_new z
#' @param ... x
#'
#' @return x
#' @export

predict.sRDA <- function(rda, exploratory_new, response_new, ...) {

    X_centered <- scale(as.matrix(exploratory_new))
    X_centered <- X_centered[, !colSums(!is.finite(X_centered))]

    Y_centered <- scale(as.matrix(response_new))
    Y_centered <- Y_centered[, !colSums(!is.finite(Y_centered))]

    Res_X <- X_centered
    results <- list()
    for (ith in 1:length(rda)) {
        XI <- scale(Res_X %*% rda[[ith]]$ALPHA)
        # update residual X
        Res_X = apply(Res_X, 2, function(Xcol) {
                Xcol - solve(t(XI) %*% XI) %*% t(XI) %*% Xcol %*% t(XI)
            }
        )
        ETA <- scale(Y_centered %*% rda[[ith]]$BETA)

        result <- list()

        result$XI <- XI
        result$ETA <- ETA
        result$cor_test <- cor.test(XI, ETA)

        results[[ith]] <- result
    }
    return(results)
}
