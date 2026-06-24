## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NormalBlockBase ############################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' R6 abstract class for a generic sparse Normal Block model
NormalBlockBase <- R6::R6Class(
  classname = "NormalBlockBase",
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @field data object of NormalBlockData class, with responses and design matrix
    data  = NULL,

    #' @description Create a new [`NormalBlockBase`] object.
    #' @param data object of NormalBlockData class, with responses and design matrix
    #' @param q number of block/cluster
    #' @param sparsity sparsity penalty on the network density
    #' @param control structured list of more specific parameters, to generate with NB_control
    #' @param zero_inflation whether the concrete subclass models zero-inflation;
    #' set by the ZI subclasses themselves, not meant to be set by the end user.
    #' When `FALSE`, the (costly) zero-inflation probability fit (`kappa`/`B0`)
    #' is skipped entirely, since it would otherwise never be used downstream.
    #' @return A new [`NormalBlockBase`] object
    initialize = function(data, q, sparsity = 0, control = NB_control(), zero_inflation = FALSE) {
      self$data <- data

      stopifnot("There cannot be more blocks than there are entities to cluster" = q <= ncol(self$data$Y))

      ## variant (either diagonal or spherical residuals covariance)
      private$res_covariance <- control$noise_covariance

      ## pointer to the chosen optimization function
      private$optimizer <- ifelse(control$heuristic,
                                  private$heuristic_optimize,
                                  private$EM_optimize)
      ## name of the chosen clustering heuristic, looked up in
      ## private$clustering_methods by heuristic_clustering()
      private$approx <- control$heuristic
      private$clustering_approx <- control$clustering_approx

      ## penalty mask
      private$sparsity_ <- sparsity
      weights <- matrix(1, q, q)
      diag(weights) <- 0
      if (!is.null(control$sparsity_weights)) {
        weights <- control$sparsity_weights
      }
      private$weights <- weights

      cl0 <- control$clustering_init
      if (!is.null(cl0)) {
        if (!is.vector(cl0) & !is.matrix(cl0)) stop("Labels must be encoded in vector of labels or indicator matrix")
        if (is.vector(cl0)) {
          if (any(cl0 < 1 | cl0 > q))
            stop("Cluster labels must be between 1 and q")
          if (length(cl0) != self$p)
            stop("Cluster labels must match the number of Y's columns")
          if (length(unique(cl0)) != q)
            stop("The number of clusters in the initial clustering must be equal to q.")
          cl0 <- as_indicator(cl0)
        } else {
          if (nrow(cl0) != self$p)
            stop("Cluster-indicating matrix must have as many rows as Y has columns")
          if (ncol(cl0) != q)
            stop("Cluster-indicating matrix must have q columns")
          if (min(colSums(cl0)) < 1)
            stop("There cannot be empty clusters in the initial clustering matrix.")
        }
        private$C <- cl0
      } else {
        private$C <- matrix(NA, self$data$n, q)
      }

      ## Zero-inflation probabilities (kappa/B0) and the resulting fixed
      ## log-likelihood contribution (ZI_cond_mean) are only ever read by the
      ## ZI subclasses (see zi_diag_normal_inference(), and the ZI EM_optimize()/
      ## fitted/model_par methods) -- skip the p logistic regressions entirely
      ## for plain (non zero-inflated) models, where they would just be wasted
      ## work (private$kappa/B0/ZI_cond_mean stay at their NA default).
      if (zero_inflation) {
        if(self$data$npY < self$n * self$p){
          B0_list <- lapply(1:self$data$p,
                            f <- function(j){
                              df <- data.frame("zeros" = data$zeros[,j], self$data$X0)
                              model <- glm(zeros ~ 0 + ., family=binomial(link = "logit"), data=df)
                              return(model$coefficients)})
          private$B0 <- t(sapply(B0_list, unlist))
        }else{
          private$B0 <- matrix(rep(-Inf, self$data$p * self$data$d0), nrow = self$data$d0)
        }

        private$kappa <- apply(self$data$X0 %*% private$B0, MARGIN = c(1,2), FUN = sigmoid)
        private$ZI_cond_mean <-
          sum(xlogy(data$zeros, private$kappa)) +
          sum(xlogy(data$zeros_bar, 1 - private$kappa))
      }
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Setters    ------------------------
    #' @description
    #' Update a [`NormalBlockBase`] object
    #'
    #' All possible parameters of the child classes
    #' @param B regression matrix
    #' @param dm1 diagonal vector of inverse variance matrix (variables level)
    #' @param C the matrix of groups memberships (posterior probabilities)
    #' @param Omegaq groups inverse variance matrix
    #' @param gamma  variance of posterior distribution of W
    #' @param mu mean for posterior distribution of W
    #' @param kappa vector of zero-inflation probabilities
    #' @param alpha vector of groups probabilities
    #' @param M variational mean for posterior distribution of W
    #' @param S variational diagonal of variances for posterior distribution of W
    #' @param ll_list  list of log-lik (elbo) values
    #' @return Update the current [`normal`] object
    update = function(B = NA,
                      dm1 = NA,
                      C = NA,
                      Omegaq = NA,
                      gamma = NA,
                      mu = NA,
                      kappa = NA,
                      alpha = NA,
                      M = NA,
                      S = NA,
                      ll_list = NA) {
      if (!anyNA(B))       private$B       <- B
      if (!anyNA(dm1))     private$dm1     <- dm1
      if (!anyNA(C))       private$C       <- C
      if (!anyNA(Omegaq))  private$Omegaq  <- Omegaq
      if (!anyNA(gamma))   private$gamma   <- gamma
      if (!anyNA(kappa))   private$kappa   <- kappa
      if (!anyNA(mu))      private$mu      <- mu
      if (!anyNA(alpha))   private$alpha   <- alpha
      if (!anyNA(M))       private$M       <- M
      if (!anyNA(S))       private$S       <- S
      if (!anyNA(ll_list)) private$ll_list <- ll_list
    },

    #' @description calls optimization (EM or heuristic) and updates relevant fields
    #' @param control a list for controlling the optimization proces
    #' @return optimizes the model and updates its parameters
    optimize = function(control = list(niter = 100, threshold = 1e-4)) {
      optim_out <- private$optimizer(control)
      do.call(self$update, optim_out)
    },

    #' @description Create a clone of the current [`NormalBlockBase`] object after splitting cluster `cl`
    #' We split the cluster according to the species variances
    #' @param index index (integer) of the cluster to split
    #' @param in_place should the split applied to the object itself, or should a copy be sent?
    #' default FALSE (send a copy)
    #' @return A new [`NormalBlockBase`] object
    split = function(index, in_place = FALSE) {
      ## update private fields related to group parameters
      ## C, Omegaq, M, S, sparsity_weights

      ## indices of individuals split within the cluster
      cl  <- self$clustering == index
      var <- 1/private$dm1; var_median <- median(var[cl])
      split1 <- (var > var_median) & cl ;  split2 <- (var <= var_median) & cl

      ## Cluster split
      new_C <- cbind(private$C, .Machine$double.eps)
      new_C[split1, index] <- new_C[split1, index] - .Machine$double.eps
      new_C[split2, self$q + 1] <- new_C[split2, index]
      new_C[split2, index] <- .Machine$double.eps
      new_C <- new_C / rowSums(new_C)

      ## Variational means
      new_M <- cbind(private$M, 0)
      new_M[split2, self$q + 1] <- new_M[split2, index]
      new_M[split2, index] <- 0

      ## Variational variances
      if (is.matrix(private$S)) {
        new_S <- cbind(private$S, 0.1)
        new_S[split2, self$q + 1] <- new_C[split2, index]
        new_S[split2, index] <- 0.1
      } else {
        new_S <- c(private$S, mean(private$S))
      }

      ## Precision matrix
      new_Omegaq <- cbind(rbind(private$Omegaq,  0), 0)
      new_Omegaq[index, index] <- private$Omegaq[index, index]/2
      new_Omegaq[self$q + 1, self$q + 1] <- private$Omegaq[index, index]/2

      ## Sparsity weights
      if (self$q == 1) {
        new_weights <- matrix(c(0,1,1,0), 2, 2)
      } else {
        weights_cl <-  private$weights[index, setdiff(1:self$q, index)]
        weights_cl <-  c(weights_cl, mean(weights_cl))
        new_weights <- cbind(rbind(private$weights, weights_cl, deparse.level = 0),
                             c(weights_cl, 0))
      }

      if (in_place) {
        self$update(C = new_C, Omegaq = new_Omegaq, M = new_M, S = new_S)
        self$sparsity_weights <- new_weights
        return(invisible(self))
      } else {
        new_NB <- self$clone()
        new_NB$update(C = new_C, Omegaq = new_Omegaq, M = new_M, S = new_S)
        new_NB$sparsity_weights <- new_weights
        return(invisible(new_NB))
      }
    },

    #' @description generate and select a set of candidate models
    #' by splitting the clusters of the current model
    candidates_split = function() {
      # do not split groups with less than 2 guys
      candidates <- map((1:self$q)[self$cluster_sizes > 1], self$split)
      # keep candidates with at least 2 guys per cluster and non empty split
      clustering_sizes <- map(candidates, "clustering") %>% map(table)
      min_sizes  <- clustering_sizes %>% map_dbl(min)
      n_clusters <- clustering_sizes %>% map_dbl(length)
      candidates <- candidates[min_sizes > 1 & n_clusters == self$q + 1]

      for (i in seq_along(candidates))
        candidates[[i]]$optimize(list(niter = 5, threshold = 1e-4))
      candidates
    },

    #' @description generate and select a set of candidate models
    #' by merging the clusters of the current model
    candidates_merge = function() {
      stopifnot("need at least two clusters to merge them" = self$q > 1)
      candidates <- map(combn(self$q, 2, simplify = FALSE), self$merge)
      for (i in seq_along(candidates))
        candidates[[i]]$optimize(list(niter = 5, threshold = 1e-4))
      candidates
    },

    #' @description Create a clone of the current [`NormalBlockBase`] object after merging clusters `cl1` and `cl2`
    #' @param indices indices (couple of integer) of the clusters to merge
    #' @param in_place should the split applied to the object itself, or should a copy be sent?
    #' default FALSE (send a copy)
    #' @return A new [`NormalBlockBase`] object
    merge = function(indices, in_place=FALSE) {

      ## sorting by increasing group label
      indices <- sort(indices)

      ## Cluster merge
      new_C <- private$C[, -indices[2]]
      new_C[, indices[1]] <- private$C[, indices[1]] + private$C[, indices[2]]

      ## Variational means
      new_M <- private$M[, -indices[2]]
      new_M[, indices[1]] <- .5 * (private$M[, indices[1]] + private$M[, indices[2]])

      ## Variational variances
      if (is.matrix(private$S)) {
        new_S <- private$S[, -indices[2]]
        new_S[, indices[1]] <- .5 * (private$S[, indices[1]] + private$S[, indices[2]])
      } else {
        new_S <- private$S[-indices[2]]
        new_S[indices[1]] <- .5 * (private$S[indices[1]] + private$S[indices[2]])
      }

      ## Precision matrix
      new_Omegaq <- private$Omegaq[-indices[2], -indices[2]]
      new_Omegaq[indices[1], indices[1]] <-
        .5 * (private$Omegaq[indices[1], indices[1]] + private$Omegaq[indices[2], indices[2]])

      ## Sparsity weights
      new_weights <-  private$weights[-indices[2], -indices[2]]

      if (in_place) {
        self$update(C = new_C, Omegaq = new_Omegaq, M = new_M, S = new_S)
        self$sparsity_weights <- new_weights
        return(self)
      } else {
        new_NB <- self$clone()
        new_NB$update(C = new_C, Omegaq = new_Omegaq, M = new_M, S = new_S)
        new_NB$sparsity_weights <- new_weights
        return(new_NB)
      }
    },

    #' @description Predicts observations Y for new covariates X.
    #' @param new_X new set of covariates.
    #' @return A n*p prediction matrix for new observations
    predict = function(new_X){
      return(new_X %*% private$B)
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Extractors ------------------------
    #' @description Extract interaction network in the latent space
    #' @param type edge value in the network. Can be "support" (binary edges), "precision" (coefficient of the precision matrix) or "partial_cor" (partial correlation between species)
    #' @return a square matrix of size `self$q`
    latent_network = function(type = c("partial_cor", "support", "precision")) {
      net <- switch(
        match.arg(type),
        "support"     = 1 * (private$Omegaq != 0 & !diag(TRUE, ncol(private$Omegaq))),
        "precision"   = private$Omegaq,
        "partial_cor" = {
          tmp <- -private$Omegaq / tcrossprod(sqrt(diag(private$Omegaq))); diag(tmp) <- 1
          tmp
        }
      )
      ## Enforce sparse Matrix encoding to avoid downstream problems with igraph::graph_from_adjacency_matrix
      ## as it fails when given dsyMatrix objects
      Matrix::Matrix(net, sparse = TRUE)
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Graphical methods------------------
    #' @param type char for line type (see plot.default)
    #' @param log char for logarithmic axes (see plot.default)
    #' @param neg boolean plot negative log-likelihood (useful when log="y")
    #' @description plots log-likelihood values during model optimization
    plot_loglik = function(type = "b", log = "xy", neg = TRUE) {
      neg <- ifelse(neg, -1, 1)
      plot(seq_along(self$objective), neg * self$objective, type = type, log = log)
    },

    #' @description plot the latent network.
    #' @param type edge value in the network. Either "precision" (coefficient of the precision matrix) or "partial_cor" (partial correlation between species).
    #' @param output Output type. Either `igraph` (for the network) or `corrplot` (for the adjacency matrix)
    #' @param edge.color Length 2 color vector. Color for positive/negative edges. Default is `c("#F8766D", "#00BFC4")`. Only relevant for igraph output.
    #' @param node.labels vector of character. The labels of the nodes. The default will use the column names ot the response matrix.
    #' @param remove.isolated if `TRUE`, isolated node are remove before plotting. Only relevant for igraph output.
    #' @param layout an optional igraph layout. Only relevant for igraph output.
    #' @param plot logical. Should the final network be displayed or only sent back to the user. Default is `TRUE`.
    plot_network = function(type            = c("partial_cor", "support"),
                            output          = c("igraph", "corrplot"),
                            edge.color      = c("#F8766D", "#00BFC4"),
                            remove.isolated = FALSE,
                            node.labels     = NULL,
                            layout          = igraph::layout_in_circle,
                            plot = TRUE) {
      if(anyNA(private$Omegaq)) stop("NA in the precision matrix")

      type   <- match.arg(type)
      output <- match.arg(output)

      net <- self$latent_network(type)

      if (output == "igraph") {
        G <-  igraph::graph_from_adjacency_matrix(net, mode = "undirected", weighted = TRUE, diag = FALSE)

        if (!is.null(node.labels)) {
          igraph::V(G)$label <- node.labels
        } else {
          igraph::V(G)$label <- unlist(lapply(1:ncol(net), f <- function(x) paste0("Cluster_", x)))
        }
        ## Nice nodes
        V.deg <- igraph::degree(G)/sum(igraph::degree(G))
        igraph::V(G)$label.cex <- V.deg / max(V.deg) + .5
        igraph::V(G)$size <- tabulate(self$clustering, nbins = self$q) * 100 / self$p
        igraph::V(G)$label.color <- rgb(0, 0, .2, .8)
        igraph::V(G)$frame.color <- NA
        ## Nice edges
        igraph::E(G)$color <- ifelse(igraph::E(G)$weight > 0, edge.color[1], edge.color[2])
        if (type == "support")
          igraph::E(G)$width <- abs(igraph::E(G)$weight)
        else
          igraph::E(G)$width <- 15*abs(igraph::E(G)$weight)

        if (remove.isolated) {
          G <- igraph::delete.vertices(G, which(igraph::degree(G) == 0))
        }
        if (plot) plot(G, layout = layout)
      }
      if (output == "corrplot") {
        if (plot) {
          if (ncol(net) > 100)
            colnames(net) <- rownames(net) <- rep(" ", ncol(net))
          G <- net
          diag(net) <- 0
          corrplot::corrplot(as.matrix(net), method = "color", is.corr = FALSE, tl.pos = "td", cl.pos = "n", tl.cex = 0.5, type = "upper")
        } else  {
          G <- net
        }
      }
      invisible(G)
    },

    #' @description plot together latent network and log-likelihood values during model optimization
    plot = function(){
      self$plot_loglik(type = "b", log = "xy", neg = TRUE)
    },


    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## S3 methods ----------------------------
    #' @description User friendly print method
    #' @param model First line of the print output
    print = function(model = paste("A", self$who_am_I, ".\n")) {
      cat(model)
      cat("===========================================================================\n")
      print(as.data.frame(round(self$criteria, digits = 3), row.names = ""))
      cat("===========================================================================\n")
      cat("* Useful fields\n")
      cat("    $model_par, $posterior_par / $var_par, $clustering \n")
      cat("    $loglik, $BIC, $ICL, $objective, $nb_param, $criteria\n")
      cat("* Useful S3 methods\n")
      cat("    print(), coef(), sigma(), fitted(), predict() \n")
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ## PRIVATE MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(
    B                 = NA, # regression matrix
    dm1               = NA, # diagonal vector of inverse variance matrix (variables level)
    C                 = NA, # the matrix of posterior probabilities (tau) or group affectation
    Omegaq            = NA, # precision matrix for clusters
    kappa             = NA, # vector of zero-inflation probabilities
    B0                = NA, # vector of zero-inflation regression matrix
    alpha             = NA, # vector of groups probabilities
    gamma             = NA, # variance of  posterior distribution of W
    mu                = NA, # mean for posterior distribution of W
    M                 = NA, # variational mean for posterior distribution of W
    S                 = NA, # variational diagonal of variances for posterior distribution of W
    optimizer         = NA, # a link to the function that perform the optimization
    ll_list           = NA, # list of log-likelihoods or ELBOs
    sparsity_         = NA, # scalar controlling the overall sparsity
    weights           = NA, # sparsity weights specific to each pairs of group
    res_covariance    = NA, # shape of the residuals covariance (diagonal or spherical)
    approx            = NA, # use approximation/heuristic approach or not
    clustering_approx = NA, # name of the clustering heuristic, key into clustering_methods
    ZI_cond_mean      = NA, # conditional mean of the ZI component (fixed)
    niter             = NA, # number of EM iterations required by the inference, if applicable

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Methods for integrated (V)EM inference --------------
    ## Each concrete subclass overrides EM_optimize() to call its
    ## Rcpp/Armadillo core (see src/exports.cpp and
    ## inst/normal_block_models.qmd); EM_initialize() (heuristic
    ## initialization) stays in R and supplies the starting values.
    EM_initialize = function() {},

    get_Omegaq = function(Sigma) {
      if (private$sparsity_ == 0) {
        Omega <- solve(Sigma)
      } else {
        glasso_out <- glassoFast::glassoFast(Sigma, rho = private$sparsity_ * self$sparsity_weights)
        if (anyNA(glasso_out$wi)) {
          warning(
            "GLasso fails, the penalty is probably too small and the system badly conditionned \n reciprocal condition number =",
            rcond(Sigma), "\n We send back the original matrix and its inverse (unpenalized)."
          )
          Omega <- solve(Sigma)
        } else {
          Omega <- Matrix::symmpart(glasso_out$wi)
        }
      }
      Omega
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Methods for heuristic inference----------------------
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## MLE of MV Normal distribution
    multivariate_normal_inference = function(){
      B     <- self$data$XtXm1 %*% self$data$XtY
      R     <- ols_residuals(self$data)
      Sigma <- cov(R)
      list(B = B, R = R, Sigma = Sigma)
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## MLE of ZI Diagonal Normal distribution
    zi_diag_normal_inference = function(){
      B     <- self$data$XtXm1 %*% self$data$XtY
      dm1   <- self$data$nY / colSums(self$data$zeros_bar * (self$data$Y - self$data$X %*% B)^2)
      for (i in 1:3) { # a couple of iterates is enough
        B     <- private$zi_diag_normal_optim_B(dm1)
        dm1   <- self$data$nY / colSums(self$data$zeros_bar * (self$data$Y - self$data$X %*% B)^2)
      }
      R <- self$data$zeros_bar * (self$data$Y - self$data$X %*% B)
      list(B = B, dm1 = dm1, kappa = private$kappa, R = R)
    },

    ## Closed-form weighted least squares solve for B: the zero-inflation
    ## mask makes the weight matrix vary by both row and column, so each
    ## column of B is solved independently (mirrors nb_optim::solve_wls in
    ## src/zi_closed_form_solvers.h -- no iterative optimizer is needed since
    ## the objective is exactly quadratic in B). Uses a pseudo-inverse rather
    ## than solve() because, e.g. with a one-hot design and enough zero
    ## inflation, a whole design level can be all-zero for some variable,
    ## making XtWX exactly singular; ginv() falls back to the minimum-norm
    ## solution instead of erroring (mirroring arma::solve()'s automatic
    ## pinv fallback on the C++ side).
    zi_diag_normal_optim_B = function(dm1) {
      DM1 <- matrix(dm1, self$data$n, self$data$p, byrow = TRUE) * self$data$zeros_bar
      B <- matrix(0, self$d, self$p)
      for (j in seq_len(self$data$p)) {
        w <- DM1[, j]
        XtWX <- crossprod(self$data$X, self$data$X * w)
        XtWy <- crossprod(self$data$X, self$data$Y[, j] * w)
        B[, j] <- MASS::ginv(XtWX) %*% XtWy
      }
      B
    },

    heuristic_optimize = function(control){
      parameters <- private$get_heuristic_parameters()
      c(parameters, list(ll_list = NA))
    },

    heuristic_Sigmaq_from_Sigma = function(Sigma){
      Sigma_q <- (t(private$C) %*% Sigma %*% private$C) / outer(colSums(private$C), colSums(private$C))
      ### TODO: why is there any NA?
      if (anyNA(Sigma_q)) {
        diag(Sigma_q)[is.na(diag(Sigma_q))] <- mean(diag(Sigma_q)[!is.na(diag(Sigma_q))])
        Sigma_q[is.na(Sigma_q)] <- 0
      }
      Sigma_q
    },

    ## Registry of clustering heuristics used to turn the OLS/ZI residuals R
    ## (n x p) into an initial clustering of the p variables into self$q
    ## groups (a vector of length p with values in 1:q). One single table
    ## instead of one ad hoc private method per algorithm -- selectable via
    ## NB_control(clustering_approx = ...) ("kmeans"/"ward2"/"sbm"/"kmeansvar";
    ## "hclustvar" is reserved for the fallback below, not user-selectable).
    ## Benchmarked on two real datasets (inst/normal_block_models.qmd): no
    ## single method dominates everywhere -- e.g. kmeansvar is unremarkable on
    ## breast cancer proteomics data but wins decisively (20/24 q values) on
    ## university webpages text data, hence it being offered as a genuine
    ## alternative rather than dropping ClustOfVar from the dependencies.
    clustering_methods = list(
      kmeans    = function(R, q) kmeans(t(R), q, nstart = 30, iter.max = 50)$cluster,
      ward2     = function(R, q) cutree(hclust(dist(1 - cor(R)), method = "ward.D2"), q),
      sbm       = function(R, q) {
        options <- list(verbosity = 0, exploreMin = q, exploreMax = q, plot = FALSE, nbCores = 1)
        mySBM <- sbm::estimateSimpleSBM(cov(R), "gaussian", estimOptions = options)
        mySBM$setModel(q)
        mySBM$memberships
      },
      kmeansvar = function(R, q) ClustOfVar::kmeansvar(X.quanti = R, init = q, nstart = 30)$cluster,
      hclustvar = function(R, q) cutree(ClustOfVar::hclustvar(R), q)
    ),

    heuristic_clustering = function(R) {
      clustering <- private$clustering_methods[[private$clustering_approx]](R, self$q)
      if (length(unique(clustering)) < self$q) {
        clustering <- private$clustering_methods$hclustvar(R, self$q)
      }
      C <- as_indicator(clustering)
      if (min(colSums(C)) < 1) warning("Initialization failed to place elements in each cluster")
      C
    }

  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field inference_method inference procedure used (heuristic or integrated with EM)
    inference_method = function() ifelse(private$approx, "heuristic", "integrated"),
    #' @field n number of samples
    n = function() self$data$n,
    #' @field p number of responses per sample
    p = function() self$data$p,
    #' @field d number of variables (dimensions in X)
    d = function() self$data$d,
    #' @field d0 number of zi variables (dimensions in X0)
    d0 = function() self$data$d0,
    #' @field q number of blocks
    q = function() as.integer(ncol(private$C)),
    #' @field n_edges number of edges of the network (non null coefficient of the sparse precision matrix Omegaq)
    n_edges  = function() sum(private$Omegaq[upper.tri(private$Omegaq, diag = FALSE)] != 0),
    #' @field model_par a list with the matrices of the model parameters: B (covariates), dm1 (species variance), Omegaq (groups precision matrix))
    model_par = function() list(B = private$B, B0 = private$B0,
                                     dm1 = private$dm1, Omegaq = private$Omegaq),
    #' @field nb_param number of parameters in the model
    nb_param = function() {
      nb_param_D <- ifelse(private$res_covariance == "diagonal", self$p, 1)
      as.integer(self$p * self$d + self$q + self$n_edges + nb_param_D)
    },
    #' @field objective evolution of the objective function during (V)EM algorithm
    objective = function() private$ll_list[-1],
    #' @field loglik (or its variational lower bound)
    loglik = function() if (private$approx) NA else private$ll_list[[length(private$ll_list)]] + self$sparsity_term,
    #' @field deviance (or its variational lower bound)
    deviance = function() -2 * self$loglik,
    #' @field BIC (or its variational lower bound)
    #' @field entropy Entropy of the conditional distribution when applicable
    entropy    = function() 0,
    BIC = function() self$deviance + log(self$n) * self$nb_param,
    #' @field ICL variational lower bound of the ICL
    ICL        = function() self$BIC + 2 * self$entropy,
    #' @field EBIC variational lower bound of the EBIC
    EBIC   = function() self$BIC + 2 * ifelse(self$n_edges > 0, self$n_edges * log(self$q), 0),
    #' @field criteria a vector with loglik, BIC and number of parameters
    criteria   = function() {
      data.frame(nb_param = self$nb_param, q = self$q, n_edges = self$n_edges, sparsity = self$sparsity,
                 loglik = self$loglik, deviance = self$deviance, BIC = self$BIC, ICL = self$ICL, EBIC = self$EBIC,
                 niter = private$niter )
    },
    #' @field sparsity (overall sparsity parameter)
    sparsity = function(value) {
      if (missing(value)) {
        private$sparsity_
      } else {
        stopifnot("must be a positive scale" = value >= 0)
        private$sparsity_ <- value
      }
    },
    #' @field sparsity_weights (weights associated to each pair of groups)
    sparsity_weights = function(value) {
      if (missing(value)) {
        private$weights
      } else {
        stopifnot("must be a q x q matrix" =
                    all(is.matrix(value), nrow(value) == ncol(value), ncol(value) == self$q))
        private$weights <- value
      }
    },
    #' @field sparsity_term (sparsity_term term in log-likelihood due to sparsity)
    sparsity_term = function() self$sparsity * sum(abs(self$sparsity_weights * private$Omegaq)),
    #' @field get_res_covariance whether the residual covariance is diagonal or spherical
    get_res_covariance = function() private$res_covariance,
    #' @field memberships cluster memberships
    memberships = function() private$C,
    #' @field clustering given as the list of elements contained in each cluster
    clustering = function() {
      cl <- get_clusters(private$C)
      names(cl) <- colnames(self$data$Y)
      cl
    },
    #' @field cluster_sizes given as a vector of cluster sizes
    cluster_sizes = function() tabulate(self$clustering, nbins = self$q),
    #' @field elements_per_cluster given as the list of elements contained in each cluster
    elements_per_cluster = function() {
      if (is.null(names(self$clustering)))
        base::split(1:self$p, self$clustering)
      else
        base::split(names(self$clustering), self$clustering)
    }
  )
)
