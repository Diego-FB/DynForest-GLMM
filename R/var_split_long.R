#' Split function to build the two daughter nodes from longitudinal predictors
#'
#' @param X Input data
#' @param Y Outcome data
#' @param timeVar A character indicating the name of time variable
#' @param nsplit_option A character indicates how the values are chosen to build the two groups for the splitting rule (only for continuous predictors). Values are chosen using deciles (\code{nsplit_option}="quantile") or randomly (\code{nsplit_option}="sample"). Default value is "quantile".
#' @param cause (Only with competing events) Number indicates the event of interest.
#' @param nodesize Minimal number of subjects required in both child nodes to split. Cannot be smaller than 1.
#' @param init (Optional) Initial values for linear mixed models
#' @param gradTol Numeric or \code{NULL}. Binomial/poisson longitudinal glmmTMB
#'   fits that fail the formal convergence check (\code{conv_ok = FALSE}) but whose maximum
#'   fixed-effect gradient (\code{gradient_max}) is below \code{gradTol}, whose Hessian is
#'   positive-definite (\code{pdHess = TRUE}), whose random-effect variance is not collapsed to
#'   the boundary (minimum random-effect standard deviation > 1e-4, same default as
#'   \code{lme4::isSingular()}), and whose other diagnostic checks (\code{beta_ok}, \code{theta_ok},
#'   \code{re_ok}) are all TRUE, are accepted as converged (\code{conv_ok} is set back to
#'   \code{TRUE}) and flagged with \code{used_grad_rescue = TRUE} in \code{model_param}. This
#'   targets a false-negative pattern confirmed by simulation: a formally non-converged fit whose
#'   fixed-effect estimates are already stable because the gradient is essentially flat. The
#'   \code{pdHess}/boundary-variance checks were added after a diagnostic showed that, without
#'   them, a small fraction of binomial rescues (never poisson) had catastrophic standard errors
#'   from near-complete-separation-like fits; with both checks the catastrophic rate drops from
#'   ~4.8% to ~0.5% of rescued binomial nodes. Never applied to \code{negbin} fits (a flat gradient
#'   there can reflect a poorly identified overdispersion parameter on a flat likelihood surface
#'   rather than genuine convergence). Default is \code{0.01}; \code{NULL} disables it. Does not change the
#'   optimizer's iteration budget or negbin/gaussian behavior.
#'
#' @import lcmm
#' @importFrom splines ns
#'
#' @keywords internal
#' @noRd
var_split_long <- function(X, Y, timeVar = NULL, nsplit_option = "quantile",
                           cause = 1, nodesize = 1, init = NULL,
                           gradTol = 0.01){

  X_ncol <- ncol(X$X)
  all_imp_var <- split_var <- vector("list", X_ncol)
  impur_var <- rep(Inf, X_ncol)
  Pure <- FALSE
  model_param <- list()
  threshold_var <- var_sum <- rep(NA, X_ncol)
  conv_issue <- NULL
  conv_issue_detail <- list()

  for (i in 1:X_ncol){

    colnames_X_i <- colnames(X$X)[i]

    fixed_var <- all.vars(X$model[[i]]$fixed)
    random_var <- all.vars(X$model[[i]]$random)
    model_var <- unique(c(fixed_var,random_var))

    # compute features from mixed model
    data_model <- data.frame(id = as.numeric(X$id), time = X$time, X$X[,, drop = FALSE])
    colnames(data_model)[which(colnames(data_model)=="time")] <- timeVar
    data_model <- data_model[,c("id", model_var)]

    if (is.null(init[[colnames_X_i]][[1]])){
      init[[colnames_X_i]][[1]] <- NA
    }

    family <- X$model[[i]]$family
    if (is.null(family)) family <- "gaussian"

    if (!family %in% c("gaussian", "binomial", "poisson", "negbin")) {
      stop("Only family = 'gaussian', 'binomial', 'poisson' and 'negbin' are currently implemented.")
    }

    used_grad_rescue <- FALSE

    if (family == "gaussian") {

      # Mixed model with initial values for parameters ?
      if (!is.na(init[[colnames_X_i]][[1]])){
        model_output <- tryCatch(
          hlme(
            fixed = X$model[[i]]$fixed,
            random = X$model[[i]]$random,
            subject = "id",
            data = data_model,
            B = init[[colnames_X_i]],
            maxiter = 100,
            verbose = FALSE
          ),
          error = function(e){
            return(NULL)
          }
        )

        if (is.null(model_output)){
          model_output <- tryCatch(
            hlme(
              fixed = X$model[[i]]$fixed,
              random = X$model[[i]]$random,
              subject = "id",
              data = data_model,
              maxiter = 100,
              verbose = FALSE
            ),
            error = function(e){
              return(NULL)
            }
          )
        }

      } else {

        model_output <- tryCatch(
          hlme(
            fixed = X$model[[i]]$fixed,
            random = X$model[[i]]$random,
            subject = "id",
            data = data_model,
            maxiter = 100,
            verbose = FALSE
          ),
          error = function(e){
            return(NULL)
          }
        )
      }
    }

    if (family %in% c("binomial", "poisson", "negbin")) {

      fixed_formula <- X$model[[i]]$fixed
      random_formula <- X$model[[i]]$random

      ziformula <- X$model[[i]]$ziformula
      if (is.null(ziformula)) ziformula <- ~0

      random_rhs <- paste(deparse(random_formula[[2]]), collapse = "")

      glmm_formula <- as.formula(
        paste(
          deparse(fixed_formula),
          "+ (",
          random_rhs,
          "| id)"
        )
      )

      start_glmmTMB <- NULL

      if (is.list(init[[colnames_X_i]]) &&
          !is.null(init[[colnames_X_i]]$beta) &&
          !is.null(init[[colnames_X_i]]$theta)) {

        start_glmmTMB <- list(
          beta = init[[colnames_X_i]]$beta,
          theta = init[[colnames_X_i]]$theta
        )

        if (!is.null(init[[colnames_X_i]]$betadisp)) {
          start_glmmTMB$betadisp <- init[[colnames_X_i]]$betadisp
        }
      }

      family_glmmTMB <- switch(
        family,
        binomial = stats::binomial(),
        poisson = stats::poisson(),
        negbin = glmmTMB::nbinom2()
      )

      fit_error_msg <- NA_character_
      used_cold_retry <- FALSE

      model_output <- tryCatch(
        suppressWarnings(
          glmmTMB::glmmTMB(
            formula = glmm_formula,
            data = data_model,
            family = family_glmmTMB,
            ziformula = ziformula,
            start = start_glmmTMB
          )
        ),
        error = function(e){
          fit_error_msg <<- conditionMessage(e)
          return(NULL)
        }
      )

      # Fallback: si había warm start y el ajuste falló (o no pasa los checks), reintenta en frío
      needs_cold_retry <- is.null(model_output)

      if (!needs_cold_retry) {
        conv_ok_tmp <- tryCatch(model_output$fit$convergence == 0, error = function(e) FALSE)
        beta_ok_tmp <- tryCatch(all(is.finite(glmmTMB::fixef(model_output)$cond)), error = function(e) FALSE)
        theta_ok_tmp <- tryCatch(all(is.finite(glmmTMB::getME(model_output, "theta"))), error = function(e) FALSE)
        re_ok_tmp <- tryCatch({
          re_tmp <- glmmTMB::ranef(model_output)$cond$id
          all(is.finite(as.matrix(re_tmp)))
        }, error = function(e) FALSE)

        if (!isTRUE(conv_ok_tmp) || !isTRUE(beta_ok_tmp) || !isTRUE(theta_ok_tmp) || !isTRUE(re_ok_tmp)) {
          needs_cold_retry <- TRUE
        }
      }

      if (needs_cold_retry && !is.null(start_glmmTMB)) {

        model_output_cold <- tryCatch(
          suppressWarnings(
            glmmTMB::glmmTMB(
              formula = glmm_formula,
              data = data_model,
              family = family_glmmTMB,
              ziformula = ziformula,
              start = NULL
            )
          ),
          error = function(e){
            fit_error_msg <<- conditionMessage(e)
            return(NULL)
          }
        )

        if (!is.null(model_output_cold)) {
          model_output <- model_output_cold
          used_cold_retry <- TRUE
        }
      }
    }

    if (is.null(model_output)) {
      conv_issue <- c(conv_issue, colnames_X_i)

      conv_issue_detail[[colnames_X_i]] <- list(
        variable = colnames_X_i,
        family = family,
        error_msg = fit_error_msg,
        conv_ok = FALSE,
        pdHess = NA,
        beta_ok = FALSE,
        theta_ok = FALSE,
        re_ok = FALSE,
        theta_nb = NA_real_,
        used_cold_retry = used_cold_retry,
        used_grad_rescue = used_grad_rescue,
        theta_se_log = NA_real_,
        gradient_max = NA_real_

      )

      next()
    }

    if (family == "gaussian") {
      if (model_output$gconv[1] > 1e-04 | model_output$gconv[2] > 1e-04) {
        conv_issue <- c(conv_issue, colnames_X_i)
        next()
      }
    }

    if (family %in% c("binomial", "poisson", "negbin")) {

      conv_ok <- tryCatch(
        model_output$fit$convergence == 0,
        error = function(e) FALSE
      )

      pdHess <- tryCatch(
        model_output$sdr$pdHess,
        error = function(e) FALSE
      )

      beta_ok <- tryCatch(
        all(is.finite(glmmTMB::fixef(model_output)$cond)),
        error = function(e) FALSE
      )

      theta_ok <- tryCatch(
        all(is.finite(glmmTMB::getME(model_output, "theta"))),
        error = function(e) FALSE
      )

      re_ok <- tryCatch({
        re_tmp <- glmmTMB::ranef(model_output)$cond$id
        all(is.finite(as.matrix(re_tmp)))
      }, error = function(e) FALSE)

      # ------------------------------------------------------------
      # Rescate por gradiente (gradTol): un fallo formal de
      # convergencia (conv_ok = FALSE) con beta_ok/theta_ok/re_ok
      # TRUE, gradiente maximo por debajo de gradTol, Hessiano
      # positivo-definido (pdHess) y varianza de efecto aleatorio NO
      # colapsada (re_sd_min > 1e-4, mismo umbral por defecto que
      # lme4::isSingular(), usado antes para el fallback de blme) se
      # acepta como convergido. Solo binomial/poisson (NUNCA negbin,
      # por el riesgo de falso positivo en superficies de
      # verosimilitud planas con theta no identificado). No modifica
      # el numero de iteraciones del optimizador.
      #
      # Los chequeos de pdHess y re_sd_min se anadieron tras un
      # diagnostico sobre una prueba de equidad a gran escala (6
      # marcadores, N=400, ntree=100, 30 replicas): sin ellos, ~4.8%
      # de los nodos binomiales rescatados (nunca poisson) terminaban
      # con errores estandar catastroficos (hasta miles de millones),
      # tipicamente con intercepto en torno a -33, varianza de efecto
      # aleatorio colapsada (~0.00008) y un gradiente sospechosamente
      # exacto (~1e-14) - el patron clasico de una superficie de
      # verosimilitud plana por separacion (cuasi-)completa, no de
      # una convergencia genuina. Con ambos chequeos anadidos, esa
      # tasa cae a ~0.5% (validado sobre 2316 nodos binomiales
      # rescatados de 5 replicas reales), conservando el 71.5% de los
      # rescates originales. El resto de fallos de convergencia
      # (beta_ok/theta_ok/re_ok, gradTol, familia) no cambian.
      # ------------------------------------------------------------

      if (!isTRUE(conv_ok) && !is.null(gradTol) && family %in% c("binomial", "poisson") &&
          isTRUE(beta_ok) && isTRUE(theta_ok) && isTRUE(re_ok)) {

        gradient_max_tmp <- tryCatch(
          max(abs(model_output$sdr$gradient.fixed), na.rm = TRUE),
          error = function(e) NA_real_
        )

        re_sd_min_tmp <- tryCatch({
          vc <- glmmTMB::VarCorr(model_output)$cond$id
          min(attr(vc, "stddev"))
        }, error = function(e) NA_real_)

        if (is.finite(gradient_max_tmp) && gradient_max_tmp < gradTol &&
            isTRUE(pdHess) &&
            is.finite(re_sd_min_tmp) && re_sd_min_tmp > 1e-4) {
          conv_ok <- TRUE
          used_grad_rescue <- TRUE
        }
      }

      if (!isTRUE(conv_ok) || !isTRUE(beta_ok) || !isTRUE(theta_ok) || !isTRUE(re_ok)) {

        conv_issue <- c(conv_issue, colnames_X_i)

        conv_issue_detail[[colnames_X_i]] <- list(
          variable = colnames_X_i,
          family = family,
          error_msg = NA_character_,
          conv_ok = conv_ok,
          pdHess = pdHess,
          beta_ok = beta_ok,
          theta_ok = theta_ok,
          re_ok = re_ok,
          theta_nb = if (family == "negbin") {
            tryCatch(stats::sigma(model_output), error = function(e) NA_real_)
          } else {
            NA_real_
          },
          used_cold_retry = used_cold_retry,
          used_grad_rescue = used_grad_rescue,
          theta_se_log = if (family == "negbin") {
            tryCatch({
              summ <- summary(model_output$sdr)
              idx <- which(rownames(summ) == "betadisp")
              if (length(idx) == 0) NA_real_ else summ[idx, "Std. Error"]
            }, error = function(e) NA_real_)
          } else {
            NA_real_
          },
          gradient_max = tryCatch(
            max(abs(model_output$sdr$gradient.fixed), na.rm = TRUE),
            error = function(e) NA_real_
          )
        )

        next()
      }
    }

    if (family == "gaussian") {

      init[[colnames_X_i]] <- model_output$best

      model_param[[i]] <- list(
        family = "gaussian",
        beta = model_output$best[(model_output$N[1]+1):model_output$N[2]],
        varcov = model_output$best[(model_output$N[1]+model_output$N[2]+1):
                                     (model_output$N[1]+model_output$N[2]+model_output$N[3])],
        stderr = tail(model_output$best, n = 1),
        idea0 = model_output$idea0
      )

      RE <- merge(
        unique(Y$id),
        model_output$predRE,
        all.x = TRUE,
        by.x = "x",
        by.y = "id"
      )[,-1, drop = FALSE]
    }

    if (family %in% c("binomial", "poisson", "negbin")) {

      pl <- tryCatch(
        model_output$obj$env$parList(),
        error = function(e) NULL
      )

      init[[colnames_X_i]] <- list(
        beta = tryCatch(glmmTMB::getME(model_output, "beta"), error = function(e) NULL),
        theta = tryCatch(glmmTMB::getME(model_output, "theta"), error = function(e) NULL),
        betadisp = if (!is.null(pl) && "betadisp" %in% names(pl)) pl$betadisp else NULL
      )

      theta_nb <- NA_real_

      if (family == "negbin") {
        theta_nb <- tryCatch(
          stats::sigma(model_output),
          error = function(e) NA_real_
        )
      }

      poisson_diag <- list(theta_se_log = NA_real_, p_value_poisson = NA_real_, prefers_poisson = NA)

      if (family == "negbin") {
        theta_se_log <- tryCatch({
          summ <- summary(model_output$sdr)
          idx <- which(rownames(summ) == "betadisp")
          if (length(idx) == 0) NA_real_ else summ[idx, "Std. Error"]
        }, error = function(e) {
          NA_real_
        })

        if (is.finite(theta_se_log) && theta_se_log > 0) {
          p_val <- stats::pnorm(-1 / theta_se_log)
          poisson_diag <- list(
            theta_se_log = theta_se_log,
            p_value_poisson = p_val,
            prefers_poisson = p_val > 0.05
          )
        }
      }

      model_param[[i]] <- list(
        family = family,
        fixed = glmmTMB::fixef(model_output)$cond,
        varcov = as.matrix(glmmTMB::VarCorr(model_output)$cond$id),
        fixed_formula = X$model[[i]]$fixed,
        random_formula = X$model[[i]]$random,
        formula = glmm_formula,
        ziformula = ziformula,
        theta_nb = theta_nb,
        pdHess = pdHess,
        conv_ok = conv_ok,
        used_cold_retry = used_cold_retry,
        used_grad_rescue = used_grad_rescue,
        gradient_max = tryCatch(
          max(abs(model_output$sdr$gradient.fixed), na.rm = TRUE),
          error = function(e) NA_real_
        ),
        poisson_diag = poisson_diag
      )

      RE_raw <- glmmTMB::ranef(model_output)$cond$id
      RE_raw$id <- as.numeric(rownames(RE_raw))

      RE <- merge(
        unique(Y$id),
        RE_raw,
        all.x = TRUE,
        by.x = "x",
        by.y = "id"
      )[,-1, drop = FALSE]
    }

    ###########################

    # Other features to compute?

    ###########################

    data_summaries <- RE

    mtry2 <- ncol(data_summaries) # How many features do we draw?
    #mtry_sum <- sample(1:ncol(data_summaries), mtry2)
    mtry_sum <- seq(ncol(data_summaries))

    impur_sum <- rep(Inf, mtry2)
    all_imp_sum <- lapply(seq(mtry2), FUN = function(x) list(Inf, Inf))
    split_sum <- vector("list", mtry2)
    split_sum_threholds <- rep(NA, mtry2)

    for (i_sum in mtry_sum){

      if (!all(is.na(data_summaries[,i_sum]))){

        nsplit <- ifelse(length(unique(na.omit(data_summaries[,i_sum])))>10,
                         10, length(unique(na.omit(data_summaries[,i_sum]))))

        if (nsplit>1) {

          if (nsplit>2){
            if (nsplit_option == "quantile"){
              split_threholds <- unique(quantile(data_summaries[,i_sum], probs = seq(0,1,1/nsplit),
                                                 na.rm = T)[-c(1,nsplit+1)])
            }
            if (nsplit_option == "sample"){
              split_threholds <- unique(sample(data_summaries[,i_sum], nsplit))
            }
          }else{
            split_threholds <- mean(unique(data_summaries[,i_sum]))
          }

          # remove partition according to nodesize criteria
          group_length <- lapply(split_threholds, FUN = function(x){
            table(data_summaries[,i_sum]<=x)
          })

          split_nodesize_ok <- unlist(lapply(group_length, FUN = function(x) !any(x<nodesize)))
          split_threholds <- split_threholds[split_nodesize_ok]
          split_threholds_length <- length(split_threholds)

          if (split_threholds_length>0){ # could happened with tie values

            # Find best feature partition
            split_sum_list <- lapply(seq(split_threholds_length), FUN = function(x){

              split <- ifelse(data_summaries[,i_sum]<=split_threholds[x],1,2)

              if ((length(unique(split))>1)&(all(table(split)>=nodesize))){
                # Evaluate the partition
                impur_res <- impurity_split(Y, split, cause = cause)

                impur <- impur_res$impur
                imp_list <- impur_res$imp_list
              }else{
                impur <- Inf
                imp_list <- list(Inf, Inf)
              }

              return(list(split = split, impur = impur, imp_list = imp_list))

            })

            partition_sum_impur <- unlist(lapply(split_sum_list, function(x) return(x$impur)))

            if (any(partition_sum_impur!=Inf)){
              best_part_sum_nsplit <- which.min(partition_sum_impur)
              split_sum[[i_sum]] <- split_sum_list[[best_part_sum_nsplit]]$split
              impur_sum[i_sum] <- split_sum_list[[best_part_sum_nsplit]]$impur
              all_imp_sum[[i_sum]] <- split_sum_list[[best_part_sum_nsplit]]$imp_list
              split_sum_threholds[i_sum] <- split_threholds[best_part_sum_nsplit]
            }
          }
        }
      }
    }

    if (any(impur_sum!=Inf)){
      best_part_sum <- which.min(impur_sum)
      var_sum[i] <- mtry_sum[best_part_sum]
      split_var[[i]] <- split_sum[[best_part_sum]]
      impur_var[i] <- impur_sum[best_part_sum]
      all_imp_var[[i]] <- all_imp_sum[[best_part_sum]]
      threshold_var[i] <- split_sum_threholds[best_part_sum]
    }
  }

  if (all(impur_var==Inf)){
    return(list(Pure=TRUE))
  }

  var_split <- which.min(impur_var)

  return(list(split = split_var[[var_split]], impur = min(impur_var), impur_list = all_imp_var[[var_split]],
              variable = var_split, variable_summary = var_sum[var_split], threshold = threshold_var[var_split],
              model_param = list(model_param[[var_split]]), conv_issue = conv_issue, conv_issue_detail = conv_issue_detail,
              init = init, Pure = Pure))
}
