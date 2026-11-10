#' Cumulative link mixed model fitted by Laplace approximation with RTMB
#'
#' Internal fitting routine used by \code{var_split_long()} for markers with
#' \code{family = "ordinal"}. The model is
#' \eqn{P(y_{ij} \le k) = F(\theta_k - x_{ij}^T\beta - z_{ij}^T b_i)}, with
#' ordered thresholds \eqn{\theta}, \eqn{b_i \sim N(0, \Sigma)} and \eqn{F} the
#' logit, probit or complementary log-log link. The random effects are
#' integrated out with the Laplace approximation of RTMB, so the gradient used
#' by the optimizer is exact. \eqn{\Sigma = LL^T} is parametrised through its
#' (unconstrained) Cholesky factor, which allows singular covariance matrices.
#'
#' @param y Integer response coded 1..K.
#' @param X Fixed-effect design matrix without intercept (the thresholds play
#'   that role).
#' @param Z Random-effect design matrix.
#' @param id Integer subject index 1..N.
#' @param K Number of response categories.
#' @param link \code{"logit"}, \code{"probit"} or \code{"cloglog"}.
#' @param start Optional list (\code{a1}, \code{dlt}, \code{beta}, \code{Lvec})
#'   with starting values, usually the parameters of the parent node.
#' @param iter.max,eval.max Limits passed to \code{nlminb()}; the defaults are
#'   those used by \code{glmmTMB} for the other families.
#'
#' @return \code{NULL} if the optimization failed, otherwise a list with the
#'   log-likelihood, thresholds (\code{alpha}), fixed effects (\code{beta}),
#'   random-effect covariance (\code{varcov}), conditional modes
#'   (\code{ranef}), gradient and Hessian of the fixed parameters, RDM
#'   (\code{rdm}), the optimizer convergence code (\code{conv}) and the fitted
#'   parameters (\code{par}) for warm starts.
#'
#' @keywords internal
#' @noRd
fit_ordinal_rtmb <- function(y, X, Z, id, K, link = "logit", start = NULL,
                             iter.max = 300, eval.max = 400){

  n <- length(y)
  q <- ncol(Z)
  p <- ncol(X)
  N <- max(id)
  nL <- q * (q + 1) / 2

  idx_low <- which(lower.tri(diag(q), diag = TRUE), arr.ind = TRUE)
  idx_low <- idx_low[order(idx_low[, 1], idx_low[, 2]), , drop = FALSE]

  y_up <- pmin(y, K - 1L)
  y_low <- pmax(y - 1L, 1L)
  is_K <- as.numeric(y == K)
  is_1 <- as.numeric(y == 1L)

  Fcdf <- switch(
    link,
    logit = function(x) RTMB::plogis(x),
    probit = function(x) RTMB::pnorm(x),
    cloglog = function(x) 1 - exp(-exp(x)),
    stop("Unsupported ordinal link.")
  )

  nll_fun <- function(par){

    # Fuera del entorno global RTMB no sobrecarga estos operadores por si solo
    "[<-" <- RTMB::ADoverload("[<-")
    "c" <- RTMB::ADoverload("c")

    RTMB::getAll(par)

    theta <- cumsum(c(a1, exp(dlt)))

    L <- matrix(0, q, q)
    for (k in seq_len(nrow(idx_low))){
      L[idx_low[k, 1], idx_low[k, 2]] <- Lvec[k]
    }

    B <- U %*% t(L)

    eta <- if (p > 0) as.vector(X %*% beta) else rep(0, n)
    for (k in seq_len(q)){
      eta <- eta + Z[, k] * B[id, k]
    }

    p_up <- is_K + (1 - is_K) * Fcdf(theta[y_up] - eta)
    p_low <- (1 - is_1) * Fcdf(theta[y_low] - eta)

    -sum(log(p_up - p_low + 1e-12)) - sum(RTMB::dnorm(U, 0, 1, log = TRUE))
  }

  if (is.null(start)){
    cum <- cumsum(tabulate(y, K))[-K] / n
    cum <- pmin(pmax(cum, 0.01), 0.99)
    th0 <- stats::qlogis(cum)
    start <- list(
      a1 = th0[1],
      dlt = if (K > 2) log(pmax(diff(th0), 0.1)) else numeric(0),
      beta = rep(0, p),
      Lvec = if (q == 2) c(0.7, 0, 0.3) else c(0.7, rep(0.05, nL - 1))
    )
  }

  par <- c(start, list(U = matrix(0, N, q)))

  obj <- RTMB::MakeADFun(nll_fun, par, random = "U", silent = TRUE)

  opt <- tryCatch(
    suppressWarnings(
      stats::nlminb(obj$par, obj$fn, obj$gr,
                    control = list(iter.max = iter.max, eval.max = eval.max))
    ),
    error = function(e) NULL
  )

  if (is.null(opt) || !is.finite(opt$objective)) return(NULL)

  par_hat <- obj$env$parList(opt$par)

  L <- matrix(0, q, q)
  for (k in seq_len(nrow(idx_low))){
    L[idx_low[k, 1], idx_low[k, 2]] <- par_hat$Lvec[k]
  }

  g <- tryCatch(as.numeric(obj$gr(opt$par)), error = function(e) rep(NA_real_, length(opt$par)))
  H <- tryCatch(stats::optimHess(opt$par, obj$fn, obj$gr), error = function(e) NULL)

  # RDM = g' H^-1 g; solo es interpretable con Hessiana definida positiva
  rdm <- tryCatch({
    if (any(eigen(H, symmetric = TRUE, only.values = TRUE)$values <= 0)) NA_real_
    else as.numeric(t(g) %*% solve(H) %*% g)
  }, error = function(e) NA_real_)

  ranef <- par_hat$U %*% t(L)
  colnames(ranef) <- colnames(Z)

  list(
    logLik = -opt$objective,
    alpha = cumsum(c(par_hat$a1, exp(par_hat$dlt))),
    beta = stats::setNames(par_hat$beta, colnames(X)),
    varcov = `dimnames<-`(L %*% t(L), list(colnames(Z), colnames(Z))),
    ranef = ranef,
    gradient = g,
    hessian = H,
    rdm = rdm,
    conv = opt$convergence,
    par = par_hat[c("a1", "dlt", "beta", "Lvec")]
  )
}
