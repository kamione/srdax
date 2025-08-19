#' Bootstrapping test for the significance of features in sRDA
#'
#' @param rda x
#' @param method x
#' @param n_bootstrap Number of bootstrapping
#' @param ci Confidence interval
#'
#' @return a result data frame
#' @export
#'
#' @import parallel
#' @import doSNOW
#' @import foreach
#' @importFrom matrixStats rowSds
#' @importFrom purrr reduce
#' @importFrom dplyr bind_cols
#'
#' @examples
#' x
#'

run_stability_test <- function(
    rda,
    method,
    n_bootstrap = 1000,
    ci = 0.95
) {
    cat("Now accessing the stability of features: \n")

    n_cores <- parallel::detectCores() - 1 # avoid exhaustion of CPU cores
    cl <- parallel::makeCluster(n_cores, type = "SOCK")
    doSNOW::registerDoSNOW(cl)
    foreach::getDoParWorkers()

    # progress bar
    pb <- txtProgressBar(min = 0, max = n_bootstrap, width = 70, style = 3)
    progress <- function(n) {
        setTxtProgressBar(pb, n)
    }
    opts <- list(progress = progress)

    n_subjects <- dim(rda$explanatory)[1]

    results <- foreach::foreach(
        ith = 1:n_bootstrap,
        .combine = dplyr::bind_rows,
        .packages = c("srdax"),
        .options.snow = opts
    ) %dopar% {
        set.seed(ith)
        resampled_index <- sample(1:n_subjects, n_subjects, replace = TRUE)
        srda_result <- srda(
            explanatory = rda$explanatory[resampled_index, ],
            response = rda$response[resampled_index, ],
            lambdas = rda$selected_lambda,
            nonzeros = rda$selected_nonzero,
            penalization = method,
            max_iteration = 100
        )

        result <- data.frame(
            alphas = I(list(srda_result$ALPHA)),
            betas = I(list(srda_result$BETA))
        )
        return(result)
    }
    parallel::stopCluster(cl)
    close(pb)

    alphas <- list()
    betas <- list()
    for (ith in 1:n_bootstrap) {
        alphas[[ith]] <- as.data.frame(results[ith, 1][[1]])
        colnames(alphas[[ith]]) <- paste0("iteration_", ith)
        betas[[ith]] <- as.data.frame(results[ith, 2][[1]])
        colnames(betas[[ith]]) <- paste0("iteration_", ith)
    }

    # reverse the sign based on original XI
    for (ith in 1:n_bootstrap) {
        if(cor(rda$ALPHA, alphas[[ith]]) < 0) {
            alphas[[ith]] <- -alphas[[ith]]
            betas[[ith]] <- -betas[[ith]]
        }
    }

    combined_results <- list()

    combined_results$boots_alphas <- reduce(alphas, bind_cols)
    rownames(combined_results$boots_alphas) <- colnames(rda$explanatory)

    alpha_stability_results <- list()
    for (ith in 1:length(ci)) {
        alpha_stability_results[[ith]] <- .calculate_stability(
            boots = combined_results$boots_alphas,
            ci = ci[ith]
        )$is_stable
    }
    combined_results$boots_alpha_stats <- .calculate_stability(
        boots = combined_results$boots_alphas,
        ci = ci[length(ci)]
    )



    combined_results$boots_betas <- reduce(betas, bind_cols)
    rownames(combined_results$boots_betas) <- colnames(rda$response)
    combined_results$boots_betas_stats <- .calculate_stability(
        boots = combined_results$boots_betas,
        ci = ci
    )
    return(combined_results)
}

.calculate_stability <- function(boots, ci) {
    n_subjects <- dim(boots)[2]
    boots <- as.matrix(boots)
    stats <- data.frame(
        mean = rowMeans(boots, na.rm = TRUE),
        se = rowSds(boots, na.rm = TRUE)
    )
    stats$me <- stats$se * qt(ci, df = n_subjects - 1)
    stats$high <- stats$mean + stats$me
    stats$low <- stats$mean - stats$me
    stats$is_stable <- if_else(stats$high * stats$low > 0, 1, 0)
    stats$label = rownames(boots)
    rownames(stats) <- NULL
    return(stats[, c(7, 6, 1:5)])
}
