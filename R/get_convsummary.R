#' Convergence issue summary
#'
#' @param object A dynforest object.
#'
#' @return A data.frame summarising convergence issues by variable.
#'
#' @export
get_convsummary <- function(object) {
  
  if (!inherits(object, "dynforest")) {
    stop("object must be a dynforest object.")
  }
  
  if (is.null(object$rf) || !"conv_issue" %in% rownames(object$rf)) {
    stop("No conv_issue information found in object$rf.")
  }
  
  vars_long <- object$Inputs$Longitudinal
  
  families <- sapply(vars_long, function(v) {
    fam <- object$Longitudinal.model[[v]]$family
    if (is.null(fam)) fam <- "gaussian"
    fam
  })
  
  conv_df <- lapply(seq_len(object$param$ntree), function(b) {
    
    conv_b <- object$rf["conv_issue", b][[1]]
    
    if (is.null(conv_b) || length(conv_b) == 0) {
      return(NULL)
    }
    
    data.frame(
      tree = b,
      variable = as.character(unlist(conv_b)),
      stringsAsFactors = FALSE
    )
  })
  
  conv_df <- dplyr::bind_rows(conv_df)
  
  if (nrow(conv_df) == 0) {
    conv_summary <- data.frame(
      variable = vars_long,
      family = families,
      n_conv_issue = 0,
      n_trees_with_issue = 0,
      stringsAsFactors = FALSE
    )
  } else {
    conv_summary <- conv_df %>%
      dplyr::group_by(variable) %>%
      dplyr::summarise(
        n_conv_issue = dplyr::n(),
        n_trees_with_issue = dplyr::n_distinct(tree),
        .groups = "drop"
      )
    
    conv_summary <- data.frame(
      variable = vars_long,
      family = families,
      stringsAsFactors = FALSE
    ) %>%
      dplyr::left_join(conv_summary, by = "variable") %>%
      dplyr::mutate(
        n_conv_issue = ifelse(is.na(n_conv_issue), 0, n_conv_issue),
        n_trees_with_issue = ifelse(is.na(n_trees_with_issue), 0, n_trees_with_issue)
      )
  }
  
  split_df <- lapply(seq_len(object$param$ntree), function(b) {
    
    tree <- object$rf["V_split", b][[1]]
    
    if (is.null(tree) || nrow(tree) == 0) {
      return(NULL)
    }
    
    tree[tree$type == "Longitudinal" & !is.na(tree$var_split), ]
  })
  
  split_df <- dplyr::bind_rows(split_df)
  
  if (nrow(split_df) == 0) {
    split_summary <- data.frame(
      variable = vars_long,
      n_splits = 0,
      stringsAsFactors = FALSE
    )
  } else {
    split_summary <- split_df %>%
      dplyr::mutate(variable = vars_long[var_split]) %>%
      dplyr::group_by(variable) %>%
      dplyr::summarise(
        n_splits = dplyr::n(),
        .groups = "drop"
      )
    
    split_summary <- data.frame(
      variable = vars_long,
      stringsAsFactors = FALSE
    ) %>%
      dplyr::left_join(split_summary, by = "variable") %>%
      dplyr::mutate(
        n_splits = ifelse(is.na(n_splits), 0, n_splits)
      )
  }
  
  out <- conv_summary %>%
    dplyr::left_join(split_summary, by = "variable") %>%
    dplyr::mutate(
      conv_per_split = n_conv_issue / pmax(n_splits, 1)
    ) %>%
    dplyr::arrange(dplyr::desc(n_conv_issue))
  
  as.data.frame(out)
}