#' Permutation test for accessing the significance of latent variables in sRDA
#'
#' @param rda a sRDA class object
#' @param n_permutation number of permutations; default is 1,000 times
#' @return a result data frame
#' @export
#' @import parallel
#' @import doSNOW
#' @import foreach
#' @importFrom dplyr bind_rows
#' @importFrom purrr reduce

run_permutation_test <- function(rda, n_permutation = 1000) {
    # if (!.check_rda_class(rda)) {
    #     stop("The input is not a 'sRDA' class.")
    # }

    # if (length(rda) == 1) {
    #     rda <- list(rda)
    # }

    results <- list()
    for (ith in 1:length(rda)) {
        cat("Now accessing significance of Latent Variable", ith, "\n")
        permuted_ssr <- .permute(
            rda[[ith]],
            n_permutation = n_permutation
        )
        empirical_ssr <- rda[[ith]]$sum_squared_betas
        result <- data.frame(
            lvs = ith,
            empirical_ssr = empirical_ssr,
            permuted_ssr = I(list(permuted_ssr)),
            p_value = sum(permuted_ssr > empirical_ssr) / n_permutation
        )
        results[[ith]] <- result
    }
    return(reduce(results, bind_rows))
}

.permute <- function(rda, n_permutation) {
    n_cores <- parallel::detectCores() - 1 # avoid exhaustion of CPU cores
    cl <- parallel::makeCluster(n_cores, type = "SOCK")
    doSNOW::registerDoSNOW(cl)
    foreach::getDoParWorkers()

    # progress bar
    pb <- txtProgressBar(min = 0, max = n_permutation, width = 70, style = 3)
    progress <- function(n) {
        setTxtProgressBar(pb, n)
    }
    opts <- list(progress = progress)

    n_subjects <- dim(rda$explanatory)[1]

    permuted_ssr <- foreach::foreach(
        ith_perm = 1:n_permutation,
        .combine = c,
        .packages = c("srdax"),
        .options.snow = opts
    ) %dopar% {
        set.seed(ith_perm)
        resampled_index <- sample(1:n_subjects, n_subjects, replace = FALSE)
        result <- srda(
            explanatory = rda$explanatory[resampled_index, ],
            response = rda$response,
            lambdas = rda$selected_lambda,
            nonzeros = rda$selected_nonzero,
            penalization = "enet",
            max_iteration = 100
        )
        return(result$sum_squared_betas)
    }
    parallel::stopCluster(cl)
    close(pb)
    return(permuted_ssr)
}

.check_rda_class <- function(object) {
    if (class(object) == c("sRDA", "list")) {
        return(TRUE)
    } else {
        return(FALSE)
    }
}
