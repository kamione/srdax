#' Sparse Redundancy Analysis (sRDA)
#'
#' This package is a rewrite of the sRDA package available on CRAN. Please visit
#' the original package on https://cran.r-project.org/package=sRDA and their
#' paper for more details. If you use this package, please consider cite the
#' original paper of sRDA.
#'
#' @param explanatory explanatory matrix or data frame, n x p
#' @param response response matrix or data frame, n x q
#' @param penalization The penalization method: "enet", "ust", "none"
#' @param lambdas The ridge penalty parameter of the explanatory latent variable (Xi) used for enet
#' @param nonzeros The number of non-zero weights of the explanatory latent variable (Xi) used for enet or ust
#' @param max_iteration Maximum number of iterations
#' @param tolerance Convergence criteria
#' @param cv_n_folds x
#' @param parallel logical; run cross-validation in parallel or not
#' @param n_lvs Number of latent variables
#' @param seed seed number, default is 1234
#'
#' @return A sRDA object
#'
#' @importFrom elasticnet enet
#'
#' @export
#'
#' @examples
#' x
srda <- function(
    explanatory,
    response,
    penalization,
    lambdas = 1,
    nonzeros = 1,
    n_lvs = 1,
    max_iteration = 100,
    tolerance = 1e-20,
    parallel = FALSE,
    cv_n_folds,
    seed
) {

    penalization <- .check_penalization_method(penalization)
    cv <- .check_parameters(lambdas, nonzeros)

    if (cv && missing(cv_n_folds)) {
        cv_n_folds <- 5
    }

    if (missing(seed)) {
        seed <- 1234
    }

    # Multiple Latent Variables ------------------------------------------------
    if (n_lvs > 1) {
        result <- .run_multiple_lvs(
            X = explanatory,
            Y = response,
            penalization = penalization,
            lambdas = lambdas,
            nonzeros = nonzeros,
            n_lvs = n_lvs,
            tolerance = tolerance,
            max_iteration = max_iteration,
            cv_n_folds = cv_n_folds,
            parallel = parallel,
            seed = seed
        )
        class(result) <- c("sRDA", "list")
        return(result)
    }

    # Cross Validation ---------------------------------------------------------
    cv_results <- NULL
    cv_results$cv_results <- "Cross validation is not used."

    if (cv) {
        # run cross-validation to search for best lambda and nonzero
        cv_results <- .get_cv_penalty_parameters(
            explanatory = explanatory,
            response = response,
            penalization = penalization,
            lambdas = lambdas,
            nonzeros = nonzeros,
            cv_n_folds = cv_n_folds,
            max_iteration = max_iteration,
            tolerance = tolerance,
            parallel = parallel,
            seed = seed
        )
        # update parameters based on cv results
        lambda <- cv_results$best_lambda
        nonzero <- cv_results$best_nonzero
    } else {
        # update lambda and nonzero for penalized models
        lambda <- lambdas
        nonzero <- nonzeros
    }

    # Core Results -------------------------------------------------------------
    result <- .srda_core(
        X = explanatory,
        Y = response,
        penalization = penalization,
        lambda = lambda,
        nonzero = nonzero,
        max_iteration = max_iteration,
        tolerance = tolerance
    )
    result$selected_lambda <- lambda
    result$selected_nonzero <- nonzero
    result$cv_results <- cv_results$cv_results
    return(result)
}
