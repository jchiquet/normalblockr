library(combinat)

#' @description computes the Adjusted Rand Index (ARI) between two clustering
#' @returns ARI, value between 0 and 1 that describes how close the two clustering
#' are
#' @param true_clusters true clustering, given as a list of labels
#' @param inferred_clusters inferred clustering, given as a list of labels
ARI <- function(true_clusters, inferred_clusters){
  pdfCluster::adj.rand.index(true_clusters, inferred_clusters)}

#' @description computes the root mean squared error (RMSE) between two vectors
#' or matrices
#' @returns RMSE between x and y
#' @param x vector or matrix 1
#' @param y vector or matrix 2
rmse <- function(x, y){
  return(Metrics::rmse(x, y))
}

#' @description computes the root mean squared error (RMSE) between two precision matrices
#' specific to this issue because it can be required to permute the rows / columns
#' @returns RMSE between two precision matrices, one of them having potentially its rows /
#' columns permuted
#' @param omega_true true precision matrix
#' @param omega_estimate precision matrix estimator
#' @param permutation permutation to apply to omega_estimate before computing the RMSE
rmse_omega <- function(omega_true, omega_estimate, permutation = NULL) {
  if (is.null(permutation)) {
    res <- rmse(omega_true, omega_estimate)
  } else if (anyNA(permutation)) {
    res <- NA
  } else {
    res <- rmse(omega_true, omega_estimate[permutation,permutation])
  }
  res
}

#' @description computes ROC metrics (true positive, true negative, recall...)
#' @returns list of metrics to assess a network inference
#' @param omega_true true precision matrix
#' @param omega_estimate precision matrix estimator
#' @param permutation permutation to apply to omega_estimate before computing the metrics
roc_metrics <- function(omega_true, omega_estimate, permutation = NULL){
  if (!is.null(permutation)) omega_estimate <- omega_estimate[permutation, permutation]

  diag(omega_true) <- 0 ; diag(omega_estimate) <- 0

  true.nzero <- which(omega_true != 0)
  true.zero  <- which(omega_true == 0)

  nzero <- which(omega_estimate != 0)
  zero  <- which(omega_estimate == 0)

  TP <- sum(nzero %in% true.nzero)
  TN <- sum(zero %in%  true.zero) - nrow(omega_true) # removing diagonal values that do not count
  FP <- sum(nzero %in% true.zero)
  FN <- sum(zero %in%  true.nzero)

  recall    <- TP/(TP + FN) ## also recall and sensitivity
  fallout   <- FP/(FP + TN) ## also 1 - specificit
  precision <- TP/(TP + FP) ## also PPR
  recall[TP + FN == 0] <- NA
  fallout[TN + FP == 0] <- NA
  precision[TP + FP == 0] <- NA

  res <-  round(c(fallout,recall,precision),3)
  res[is.nan(res)] <- 0
  names(res) <- c("fallout","recall", "precision")

  return(res)
}

#' @description computes an Area Under the Curve (AUC) given a list of recall
#' and fallout values
#' @returns AUC values
#' @param recall list of recall values corresponding to each penalty
#' @param fallout list of fallout values corresponding to each penalty
auc <- function(recall, fallout){
  return(sum(diff(fallout) * (recall[-1] + recall[-length(recall)]) / 2))
}

#' @description computes an Area Under the Curve (AUC) given a precision matrix,
#' a normal-block model with various penalties and a permutation to apply to the
#' inferred precision matrix
#' @returns AUC value
#' @param omega_true true precision matrix
#' @param NB_sparse normal-block model with different penalties
#' @param permutation permutation to apply to omega_estimate before computing the metrics
get_auc <- function(omega_true, NB_sparse, permutation = NULL){
  fallout <- c() ; recall <- c()
  for(pen in NB_sparse$sparsity){
    omega_estimate <- NB_sparse$get_model(pen)$model_par$Omega
    res <- roc_metrics(omega_true, omega_estimate, permutation)
   if(!is.na(res[["fallout"]]) && !is.na(res[["recall"]])){
     fallout <- c(fallout, res[["fallout"]]) ; recall <- c(recall, res[["recall"]])
   }
  }
  if(pen == max( NB_sparse$sparsity)){recall <- rev(recall) ; fallout <- rev(fallout)}
  # One value of fallout may correspond to different recall values depending on the penalty
  fallout_unique <- unique(fallout) ; recall_unique <- c()
  for(x in fallout_unique){
    recall_unique <- c(recall_unique, max(recall[which(fallout == x)]))
  }
  return(auc(recall_unique, fallout_unique))
}

#' @description computes all the measures used to assess a Normal-Block model
#' @returns named list of measures to assess a Normal-Block model
#' @param NB_object optimized Normal-Block model with various penalties
#' @param param list of true parameters to compare the inferred parameters to
#' @param model_selection criterion (BIC, StARS...) to use to assess the model
#' @param fixed_blocks boolean to say whether the model should be assessed with fixed blocks
#' @param stability used only if model_selection = StARS, stability level for StARS
get_measures <- function(NB_object, param, model_selection, fixed_blocks = FALSE,
                         stability = 0.9) {
  # Select best sparsity level according to the chosen criterion
  if(is.numeric(model_selection)){
    model <- NB_object$get_model(model_selection)
    model_selection <- round(model_selection, 3)
  }else{model <- NB_object$get_best_model(model_selection, stability)}

  # Get best permutation of Omega according to rmse when possible
  omega_estimate <- model$model_par$Omega
  if (fixed_blocks) {
    best_perm <- NULL
    AUC = get_auc(param$Omega, NB_object, best_perm)
  }else{
    AUC = NA
    if (model$q < 7) {
      perms <- permn(model$q)
      ibest <- perms %>%
        purrr::map_dbl(\(P) rmse(param$Omega, omega_estimate[P,P])) %>%
        which.min()
      best_perm <- perms[[ibest]]
      AUC = get_auc(param$Omega, NB_object, best_perm)
    }else best_perm <- NA }

  ## get metrics
  res <- c(
    criterion = ifelse(model_selection == "StARS", paste0("StARS_", stability),
                       model_selection),
    fixed_blocks = fixed_blocks,
    AUC = AUC,
    ARI = ARI(get_clusters(param$C), model$clustering),
    rmse_B = rmse(param$B, model$model_par$B),
    rmse_D = rmse(diag(param$D), 1/model$model_par$dm1),
    rmse_kappa = ifelse(is.null(model$model_par$kappa), NA,
                        rmse(param$kappa, model$model_par$kappa)),
    rmse_fit = rmse(model$fitted, NB_object$data$Y),
    rmse_omega = rmse_omega(param$Omega, omega_estimate, best_perm),
    roc_metrics(param$Omega, omega_estimate, best_perm)
  )
  res
}
