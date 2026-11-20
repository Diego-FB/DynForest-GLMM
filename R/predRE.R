# Reduced basis for a singular random-effect covariance:
# b = V z, z ~ N(0, diag(lambda)); directions without variance stay at 0.
re_singular_basis <- function(Sigma) {
  eg <- tryCatch(eigen((Sigma + t(Sigma)) / 2, symmetric = TRUE), error = function(e) NULL)
  if (is.null(eg) || any(!is.finite(eg$values))) return(NULL)
  keep <- eg$values > max(eg$values, 0) * 1e-10
  list(V = eg$vectors[, keep, drop = FALSE], lambda = eg$values[keep])
}

# Posterior mode of the random effects given -loglik(b). Uses Sigma_inv when it
# exists, otherwise (singular covariance) optimises in the reduced basis.
# Returns NULL if the optimisation does not converge.
re_map <- function(negloglik, Sigma_inv, basis, q) {
  if (is.null(basis)) {
    f <- function(b) negloglik(b) + 0.5 * as.numeric(t(b) %*% Sigma_inv %*% b)
    opt <- tryCatch(optim(par = rep(0, q), fn = f, method = "BFGS"), error = function(e) NULL)
    if (!is.null(opt) && opt$convergence == 0) return(opt$par)
    return(NULL)
  }
  r <- length(basis$lambda)
  if (r == 0) return(rep(0, q))
  f <- function(z) negloglik(as.vector(basis$V %*% z)) + 0.5 * sum(z^2 / basis$lambda)
  opt <- tryCatch(optim(par = rep(0, r), fn = f, method = "BFGS"), error = function(e) NULL)
  if (!is.null(opt) && opt$convergence == 0) return(as.vector(basis$V %*% opt$par))
  NULL
}

#' Function to compute individual random effects using mixed model output parameters
#'
#' @param model output object from longitudinal mixed model
#' @param formula list of formula for fixed and random part
#' @param data data to compute random effect
#'
#' @import lcmm
#' @importFrom splines ns
#'
#' @return A table of random-effects in column by subjects in row
#'
#' @keywords internal
#' @noRd
predRE <- function(model, formula, data){
  
  family <- model$family
  if (is.null(family)) family <- "gaussian"
  
  if (family == "gaussian") {
    
    subject <- "id"
    dataNA <- na.omit(data)
    beta <- model$beta
    
    # Variance-covariance matrix of the random-effects
    B <- matrix(0, ncol = sum(model$idea0), nrow = sum(model$idea0))
    B[upper.tri(B, diag = TRUE)] <- model$varcov
    B <- t(B)
    B[upper.tri(B, diag = TRUE)] <- model$varcov
    
    se <- model$stderr^2
    Z <- model.matrix(formula$random, dataNA)
    X <- model.matrix(formula$fixed, dataNA)
    
    Y <- model.matrix(
      reformulate(as.character(formula$fixed)[2], intercept = FALSE),
      dataNA
    )
    
    bi <- matrix(
      NA,
      nrow = length(unique(data$id)),
      ncol = ncol(B),
      dimnames = list(unique(data$id), colnames(Z))
    )
    
    for (id in unique(data$id)) {
      
      row.id <- which(dataNA$id == id)
      Zi <- Z[row.id, , drop = FALSE]
      Xi <- X[row.id, ]
      Yi <- Y[row.id, ]
      Vi <- Zi %*% B %*% t(Zi) + se * diag(length(row.id))
      
      b <- tryCatch(
        B %*% t(Zi) %*% solve(Vi) %*% (Yi - Xi %*% beta),
        error = function(e) return(rep(NA, ncol(B)))
      )
      
      bi[rownames(bi) == id, ] <- b
    }
    
    return(list(bi = bi))
  }
  
  if (family %in% c("binomial", "poisson", "negbin")) {
    
    dataNA <- na.omit(data)
    
    beta <- model$fixed
    Sigma <- model$varcov
    
    fixed_formula <- model$fixed_formula
    random_formula <- model$random_formula
    
    ids <- unique(data$id)
    
    Z_all <- model.matrix(random_formula, dataNA)
    q <- ncol(Z_all)
    
    bi <- matrix(
      NA,
      nrow = length(ids),
      ncol = q,
      dimnames = list(ids, colnames(Z_all))
    )
    
    Sigma_inv <- tryCatch(
      solve(Sigma),
      error = function(e) NULL
    )
    
    basis <- NULL
    if (is.null(Sigma_inv)) {
      basis <- re_singular_basis(Sigma)
      if (is.null(basis)) return(list(bi = bi))
    }
    
    response_name <- all.vars(fixed_formula)[1]
    
    for (id in ids) {
      
      data_i <- dataNA[dataNA$id == id, , drop = FALSE]
      
      if (nrow(data_i) == 0) next
      
      X_i <- model.matrix(fixed_formula, data_i)
      Z_i <- model.matrix(random_formula, data_i)
      y_i <- data_i[[response_name]]
      
      negll <- function(b) {
        
        eta <- as.vector(X_i %*% beta + Z_i %*% b)
        
        if (family == "binomial") {
          loglik <- sum(stats::dbinom(
            x = y_i,
            size = 1,
            prob = stats::plogis(eta),
            log = TRUE
          ))
        }
        
        if (family == "poisson") {
          loglik <- sum(stats::dpois(
            x = y_i,
            lambda = exp(eta),
            log = TRUE
          ))
        }
        
        if (family == "negbin") {
          theta_nb <- model$theta_nb
          
          loglik <- sum(stats::dnbinom(
            x = y_i,
            size = theta_nb,
            mu = exp(eta),
            log = TRUE
          ))
        }
        
        return(-loglik)
      }
      
      b_hat <- re_map(negll, Sigma_inv, basis, q)
      
      if (!is.null(b_hat)) {
        bi[rownames(bi) == id, ] <- b_hat
      }
    }
    
    return(list(bi = bi))
  }
  
  if (family == "ordinal") {

    dataNA <- na.omit(data)

    beta <- model$fixed
    alpha <- model$alpha
    Sigma <- model$varcov
    lev <- model$levels
    link <- if (is.null(model$link)) "logit" else model$link

    Fcdf <- switch(
      link,
      logit = stats::plogis,
      probit = stats::pnorm,
      cloglog = function(x) 1 - exp(-exp(x)),
      stop("Unsupported ordinal link in predRE().")
    )

    fixed_formula <- model$fixed_formula
    random_formula <- model$random_formula

    ids <- unique(data$id)

    Z_all <- model.matrix(random_formula, dataNA)
    q <- ncol(Z_all)

    bi <- matrix(
      NA,
      nrow = length(ids),
      ncol = q,
      dimnames = list(ids, colnames(Z_all))
    )

    Sigma_inv <- tryCatch(solve(Sigma), error = function(e) NULL)
    basis <- NULL
    if (is.null(Sigma_inv)) {
      basis <- re_singular_basis(Sigma)
      if (is.null(basis)) return(list(bi = bi))
    }

    response_name <- all.vars(fixed_formula)[1]
    K <- length(lev)
    alpha_ext <- c(-Inf, alpha, Inf)

    for (id in ids) {

      data_i <- dataNA[dataNA$id == id, , drop = FALSE]
      if (nrow(data_i) == 0) next

      # observations with a level not seen in the node are dropped
      k_i <- match(as.character(data_i[[response_name]]), lev)
      keep <- !is.na(k_i)
      if (!any(keep)) next
      data_i <- data_i[keep, , drop = FALSE]
      k_i <- k_i[keep]

      X_full <- model.matrix(fixed_formula[-2], data_i)
      X_i <- X_full[, setdiff(colnames(X_full), "(Intercept)"), drop = FALSE]
      Z_i <- model.matrix(random_formula, data_i)

      negll <- function(b) {

        eta <- as.vector(Z_i %*% b)
        if (length(beta) > 0) eta <- eta + as.vector(X_i %*% beta)

        p <- Fcdf(alpha_ext[k_i + 1] - eta) - Fcdf(alpha_ext[k_i] - eta)
        return(-sum(log(pmax(p, 1e-12))))
      }

      b_hat <- re_map(negll, Sigma_inv, basis, q)

      if (!is.null(b_hat)) {
        bi[rownames(bi) == id, ] <- b_hat
      }
    }

    return(list(bi = bi))
  }

  stop("Unsupported family in predRE().")
}
