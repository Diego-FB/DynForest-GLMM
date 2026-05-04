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
  
  if (family == "binomial") {
    
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
    
    if (is.null(Sigma_inv)) {
      return(list(bi = bi))
    }
    
    response_name <- all.vars(fixed_formula)[1]
    
    for (id in ids) {
      
      data_i <- dataNA[dataNA$id == id, , drop = FALSE]
      
      if (nrow(data_i) == 0) next
      
      X_i <- model.matrix(fixed_formula, data_i)
      Z_i <- model.matrix(random_formula, data_i)
      y_i <- data_i[[response_name]]
      
      nll <- function(b) {
        
        eta <- as.vector(X_i %*% beta + Z_i %*% b)
        
        loglik <- sum(dbinom(
          x = y_i,
          size = 1,
          prob = plogis(eta),
          log = TRUE
        ))
        
        prior <- -0.5 * as.numeric(t(b) %*% Sigma_inv %*% b)
        
        return(-(loglik + prior))
      }
      
      opt <- tryCatch(
        optim(
          par = rep(0, q),
          fn = nll,
          method = "BFGS"
        ),
        error = function(e) NULL
      )
      
      if (!is.null(opt) && opt$convergence == 0) {
        bi[rownames(bi) == id, ] <- opt$par
      }
    }
    
    return(list(bi = bi))
  }
  
  stop("Unsupported family in predRE().")
}
