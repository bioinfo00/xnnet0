#' Transcriptomics data from active vs latent tuberculosis patients
#'
#'
#' @format Sample of a GEO dataset containing labeled transcriptional profiles from active and latent tuberculosis patients
#' \describe{
#'   \item{X}{expression matrix with 50 samples (rows) and 31426 probes (columns)}
#'   \item{y}{a vector of binary labels for each sample: 0 = latent tuberculosis, 1 = active tuberculosis}
#' }
#' @source <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE37250>
"GSE37250"


#' Annotation libraries
#'
#' @format Each element of the list is an annotation library in the form of
#' a list of gene sets. The annotation libraries are a selection of the ones
#' provided by Enrichr (see link below).

#' @source <https://maayanlab.cloud/Enrichr/#stats>
"annotation_libraries"


#' Training/test split
#'
#' This function splits a dataset and corresponding
#' labels into a training and test set
#' @param X transcriptomics data with samples as rows and genes as features
#' @param y binary label for each sample in X
#' @param training_fraction fraction of data for training
#' @param rnd_seed random seed for reproducible results
#' @return A list containing \code{(X_train, y_train, X_test, y_test)}.
#' @export
#' @import caret
#' @import dplyr
#' @importFrom testthat expect_true
#' @examples
#' data("GSE37250") #load Tubercolosis dataset
#' GSE37250_split = train_test_split(GSE37250$X, GSE37250$y)

train_test_split = function(X,
                            y,
                            training_fraction = 0.7,
                            rnd_seed = 666) {
  #splits a labeled dataset (X, y) in training and test set
  #split_rnd_seed is the random seed for reproducibility

  #check if X and y have same size
  testthat::expect_true(nrow(X) == length(y))

  #check that y has two levels (0, 1)
  class_frequency = table(y)/length(y)
  y_levels = class_frequency %>% names %>% sort %>% as.numeric
  testthat::expect_true(sum(y_levels == c(0, 1)) == 2, info = 'y should take on two values: 0, 1')

  print(paste('####### class 0 frequency: ',  round(class_frequency[1], 2), '#######'))
  print(paste('####### class 1 frequency: ',  round(class_frequency[2], 2), '#######'))

  if (min(class_frequency) < 0.1) warning("high class imbalance")

  set.seed(rnd_seed)
  train_index = caret::createDataPartition(y, p = training_fraction, list = F)
  X_train = X[train_index,]
  y_train = y[train_index]
  X_test = X[-train_index,]
  y_test = y[-train_index]

  print(paste('####### training set has', nrow(X_train), 'samples and', ncol(X_train), 'variables', '#######'))
  print(paste('####### test set has', nrow(X_test), 'samples and', ncol(X_test), 'variables', '#######'))

  return(list(
    X_train = X_train,
    y_train = y_train,
    X_test = X_test,
    y_test = y_test
  ))
}


#' Building xnnet
#'
#' This function builds xnnet using a labeled training set of transcriptional profiles,
#' and a list of annotation libraries
#' @param X_train data matrix
#' @param y_train binary labels
#' @param annotation_libraries a list of annotation libraries. Each library is a list of gene sets
#' gene sets
#' @param n_hidden_nodes number of hidden nodes in the network
#' @param n_input_nodes maximum number of input nodes per hidden nodes
#' @param n_unassigned number of input genes not assigned to a specific gene set
#' @param number number of bootstrap samples for network training
#' @param min_decay minimum value of the decay parameter for regularization
#' @param max_decay maximum value of the decay parameter for regularization
#' @return A list containing an xnnet object for each annotation library
#' @export
#' @examples
#' data("GSE37250") #load Tubercolosis dataset
#' data("annotation_libraries")
#' GSE37250_split = train_test_split(GSE37250$X, GSE37250$y)
#' xnnet = build_xnnet(X_train = GSE37250_split$X_train, y_train = GSE37250_split$y_train,
#' annotation_libraries = annotation_libraries)
build_xnnet = function(X_train, y_train,
                       annotation_libraries,
                       n_input_nodes = 3,
                       n_hidden_nodes = 4,
                       n_unassigned = 3,
                       number = 10,
                       min_decay = 0.1,
                       max_decay = 1){

  print('step 1 of 3: performing Limma')
  limma_results = get_limma_results(X_train, y_train)
  print('done')

  print('step 2 of 3: processing GSEA results')
  GSEA_results = get_and_process_GSEA_results(annotation_libraries,
                                              limma_results)
  print('done')

  print('step 3 of 3: cross-validating neural networks')
  xnnet = lapply(GSEA_results, function(x) build_single_xnnet(X_train = X_train,
                                                              y_train = y_train,
                                                              annotation_libraries, x,
                                                              limma_results = limma_results,
                                                              n_input_nodes = n_input_nodes,
                                                              n_hidden_nodes = n_hidden_nodes,
                                                              n_unassigned = n_unassigned,
                                                              number = number,
                                                              min_decay = min_decay,
                                                              max_decay = max_decay))
  print('done')

  return(xnnet)
}

#' Limma differential analysis
#'
#' This function performs differential gene expression with Limma on
#' the training set
#' @param X_train data matrix
#' @param y_train binary labels
#' @param adj_pval_cutoff fraction of data for training (default = 0.05)
#' @import limma
#' @importFrom stats lm mad median model.matrix predict reorder rnorm terms

get_limma_results = function(X_train, y_train, adj_pval_cutoff = 0.05) {
  #performs Limma differential expression on labeled training set

  exp_design = factor(ifelse(y_train == 0, 'class_neg', 'class_pos'))
  exp_design_matrix = model.matrix( ~ 0 + exp_design)
  contrast_matrix = limma::makeContrasts(contrasts="exp_designclass_pos - exp_designclass_neg",
                                         levels = exp_design_matrix)

  fitTrtMean = limma::lmFit(t(X_train), exp_design_matrix)
  fit_contrast = limma::contrasts.fit(fitTrtMean, contrast_matrix)
  efit_contrast = limma::eBayes(fit_contrast)

  limma_results = limma::topTable(efit_contrast, adjust.method = "BH", n = Inf)
  limma_results$pi_val = -log10(limma_results$P.Val)*limma_results$logFC

  #sort by pi-values
  limma_results = limma_results %>% dplyr::arrange(desc(abs(.data$pi_val)))
  n_degs = sum(limma_results$adj.P.Val < adj_pval_cutoff)
  if (n_degs == 0){
    #warning("no significant genes found at adj pval < 0.05")
  } else {
    #print(paste('#######', n_degs, 'DEGS detected at adj pval < 0.05', '#######'))

  }

  return(limma_results)

}


#' Gene Set Enrichment Analysis (GSEA)
#'
#' This function performs GSEA using the Limma output and an annotation library
#' @param limma_results out of Limma
#' @param annotation_library list of gene sets contained in a given annotation library
#' @param min_term_size min size of annotation term (default = 10)
#' @param max_term_size max size of annotation term (default = Inf)
#' @param rank_metric
#' @importFrom fgsea fgsea
#' @importFrom testthat expect_true

get_GSEA_results = function(limma_results,
                            annotation_library,
                            min_term_size = 10,
                            max_term_size = Inf,
                            rank_metric = 'logFC'){
  #top_n_gene_sets = length(annotation_library)) {
  #performs GSEA to identify significant annotation library
  #annotation_database is a list of gene sets (e.g. Reactome)
  #top_n_gene_sets is the number of top gene sets retained for downstream analysis

  top_n_gene_sets = length(annotation_library)
  #check overlap between limma output genes and annotation genes
  measured_genes = row.names(limma_results)
  annotated_genes = unique(unlist(annotation_library))
  overlap = length(intersect(measured_genes, annotated_genes))

  switch(rank_metric,
         logFC = {ranks = limma_results$logFC},
         pval = {ranks = limma_results$P.Value},
         pi_val = {ranks = limma_results$pi_val}
  )
  names(ranks) = row.names(limma_results)

  GSEA_results = fgsea::fgsea(
    annotation_library,
    ranks,
    nperm = 1000,
    minSize = min_term_size,
    maxSize = max_term_size
  )

  n_significant_gene_sets = sum(GSEA_results$padj < 0.1)

  #print(paste('#######', n_significant_gene_sets, 'significance gene sets detected at adj pval < 0.1', '#######'))
  #if (n_significant_gene_sets == 0)
  # warning(paste('#######', "no significant genes found at adj pval < 0.1"), '#######')

  #print(paste('####### retaining top', top_n_gene_sets, 'gene sets for downstream analysis #######'))
  GSEA_results = GSEA_results %>%  arrange(desc(abs(.data$NES))) %>% dplyr::slice(1:top_n_gene_sets)

  return(GSEA_results)

}


#' Renormalize NES produced by GSEA
#'
#' @importFrom plyr ldply
#' @param GSEA_results

renormalize_NES = function(GSEA_results){

  all_GSEA_results = plyr::ldply(GSEA_results, rbind)

  model = lm(abs(NES) ~ log10(size) + .id, data = all_GSEA_results)
  all_GSEA_results$renormalized_NES = abs(all_GSEA_results$NES) - model$fitted.values

  all_GSEA_results = all_GSEA_results %>% dplyr::filter(.data$renormalized_NES > 0) %>%
    dplyr::arrange(desc(.data$renormalized_NES))

  GSEA_results = split(all_GSEA_results, all_GSEA_results$.id)

  return(GSEA_results)

}

#' Renormalize NES produced by GSEA
#'
#' @importFrom plyr ldply
#' @param annotation_libraries
#' @param limma_results
get_and_process_GSEA_results = function(annotation_libraries,
                                        limma_results){

  GSEA_results = lapply(annotation_libraries, function(x)
    get_GSEA_results(limma_results, x))

  GSEA_results = renormalize_NES(GSEA_results)
  GSEA_results = lapply(GSEA_results, function(x) reduce_GSEA_results(x))

  return(GSEA_results)

}



build_single_xnnet = function(X_train,
                              y_train,
                              annotation_libraries,
                              GSEA_results,
                              limma_results,
                              n_input_nodes = 5,
                              n_hidden_nodes = 5,
                              n_unassigned = 3,
                              number = 10,
                              min_decay = 0.1,
                              max_decay = 1){

  selected_xnnet_nodes = select_xnnet_nodes(GSEA_results, limma_results,
                                            n_input_nodes = n_input_nodes,
                                            n_hidden_nodes = n_hidden_nodes)

  hidden_renormalized_NES = GSEA_results[which(GSEA_results$pathway %in% selected_xnnet_nodes$hidden_nodes), 'renormalized_NES']
  interpretability_score = mean(hidden_renormalized_NES)

  xnnet_binary_matrix = create_xnnet_binary_matrix(annotation_library = annotation_libraries[[unique(GSEA_results$.id)]],
                                                   selected_xnnet_nodes = selected_xnnet_nodes)
  xnnet_binary_matrix = extend_with_unassigned_inputs(xnnet_binary_matrix, limma_results = limma_results, n_unassigned = n_unassigned)


  nnetFit_output = cross_validate_xnnet(
    X_train,
    y_train,
    xnnet_binary_matrix,
    number = number,
    min_decay = min_decay,
    max_decay = max_decay
  )

  initial_wt = nnetFit_output$initial_weights
  X_train_center = nnetFit_output$X_train_center
  X_train_scale = nnetFit_output$X_train_scale
  nnetFit = nnetFit_output$nnetFit

  #subset_limma_results

  xnnet_results = list(
    interpretability_score = interpretability_score,
    GSEA_results = GSEA_results,
    limma_results = limma_results,
    nnetFit = nnetFit,
    initial_wt = initial_wt,
    xnnet_binary_matrix = xnnet_binary_matrix,
    X_train_center = X_train_center,
    X_train_scale = X_train_scale
  )
  return(xnnet_results)
}


#' Reducing GSEA results with Weighted Set Cover
#'
#' Size constrained weighted set cover problem to find top N sets while
#' maximizing the coverage of all elements.
#'
#' @param GSEA_results GSEA results (output of get_GSEA_results)
#' @param topN The number of sets (or less when it completes early) to return.
#'
#' @return A list of \code{topSets} and \code{coverage}.
#' \describe{
#'  \item{topSets}{A list of set IDs.}
#'  \item{coverage}{The percentage of IDs covered in the top sets.}
#' }
#' @importFrom parallel mclapply

reduce_GSEA_results = function(GSEA_results, topN = 10, nThreads = 1) {

  idsInSet = GSEA_results$leadingEdge
  names(idsInSet) = GSEA_results$pathway
  #costs = scale(GSEA_results$size/abs(GSEA_results$NES))
  #costs = 1/GSEA_results$renormalized_NES

  costs = rep(1, nrow(GSEA_results))

  #cat("Begin weighted set cover...\n")
  names(costs) <- names(idsInSet)
  if (.Platform$OS.type == "windows") {
    nThreads = 1
  }
  multiplier <- 10
  # we only start with the top (multiplier * topN) most
  # significant sets to reduce computational cost
  max_num_set <- multiplier * topN
  if (length(idsInSet) > max_num_set) {
    # sort by absolute of cost (1/signedLogP)
    index <- order(abs(costs), decreasing=FALSE)
    costs <- costs[index][1:max_num_set]
    idsInSet <- idsInSet[index][1:max_num_set]
  }


  s.hat <- 1.0
  # get all unique genes in all enriched sets
  all.genes <- unique(unlist(idsInSet))
  remain <- s.hat * length(all.genes)

  # final results, contains a list of gene set names
  cur.res <- c()
  # current candidates with marginal gain and size
  all.set.names <- names(idsInSet)
  mc_results <- parallel::mclapply(all.set.names, function(cur_name, cur_res, idsInSet, costs) {
    cur_gain <- marginalGain(cur_name, cur_res, idsInSet, costs)
    cur_size <- length(idsInSet[[cur_name]])
    return(data.frame(geneset.name=cur_name, gain=cur_gain, size=cur_size, stringsAsFactors=FALSE))
  }, cur_res=cur.res, idsInSet=idsInSet, costs=costs, mc.cores=nThreads)
  candidates <- mc_results %>% dplyr::bind_rows()
  topN <- min(topN, nrow(candidates))
  for (i in seq(topN)) {
    # if there is no candidates, return
    if (nrow(candidates) == 0) {
      covered.genes <- unique(unlist(idsInSet[cur.res]))
      s.hat <- length(covered.genes) / length(all.genes)
      #cat("No more candidates, ending weighted set cover\n")
      return(list(topSets=cur.res, coverage=s.hat))
    }
    # find the set with maximum marginal gain
    # tie breaker: for two sets with sname marginal gain, pick the one with
    # larger size
    candidates <- candidates[order(-candidates$gain, -candidates$size), ]
    # update remain
    remain <- remain - length(marginalBenefit(candidates[1, "geneset.name"], cur.res, idsInSet))
    cur.res <- c(cur.res, candidates[1,"geneset.name"])
    if (remain == 0) {
      covered.genes <- unique(unlist(idsInSet[cur.res]))
      s.hat <- length(covered.genes) / length(all.genes)
      #cat("Remain is 0, ending weighted set cover\n")
      # full coverage solution
      return(list(topSets=cur.res, coverage=s.hat))
    }
    # update candidates
    # first remove the one just been selected
    candidates <- candidates[-1, ]
    # recalculate gain, remove rows with gain == 0
    if (nrow(candidates) > 0) {
      mc_results <- parallel::mclapply(seq(nrow(candidates)), function(row, candidates, cur_res, idsInSet, costs){
        cur_name <- candidates[row, "geneset.name"]
        cur_gain <- marginalGain(cur_name, cur_res, idsInSet, costs)
        if(cur_gain != 0) {
          candidates[candidates$geneset.name == cur_name, "gain"] <- cur_gain
          tmp_candidate <- candidates[candidates$geneset.name == cur_name,]
          return(tmp_candidate)
        }
      }, candidates=candidates, cur_res=cur.res, idsInSet=idsInSet, costs=costs, mc.cores=nThreads)

      new_candidates <- mc_results %>% dplyr::bind_rows()
      candidates <- new_candidates
    }
  }
  # not fully covered, compute the current coverage and return
  covered.genes <- unique(unlist(idsInSet[cur.res]))
  s.hat <- length(covered.genes) / length(all.genes)
  #cat("End weighted set cover...\n")

  reduced_GSEA_results = GSEA_results %>% dplyr::filter(.data$pathway %in% cur.res)
  return(reduced_GSEA_results)
}


# return a list of genes from all.genes that has not been
# covered so far
# cur.set.name: name of the candidate gene set
# cur.res: vector of names of gene sets in current result
marginalBenefit <- function(cur.set.name, cur.res, idsInSet) {
  all.genes <- unique(unlist(idsInSet))
  cur.genes <- idsInSet[[cur.set.name]]
  if(length(cur.res) == 0){
    not.covered.genes <- cur.genes
  } else{
    covered.genes <- unique(unlist(idsInSet[cur.res]))
    not.covered.genes <- setdiff(cur.genes, covered.genes)
  }
  return(not.covered.genes)
}

marginalGain <- function(cur.set.name, cur.res, idsInSet, costs) {
  abs_cur_cost <- abs(costs[cur.set.name])
  cur.mben <- marginalBenefit(cur.set.name, cur.res, idsInSet)
  return(length(cur.mben) / abs_cur_cost)
}


#' Defining xnnet nodes
#'
#' This function selects the xnnet input and hidden nodes
#' @param reduced_GSEA_results reduced GSEA results
#' @param limma_results Limma results
#' @param n_input_nodes number of input nodes (default = 5)
#' @param n_hidden_nodes number of input nodes (default = 4)
select_xnnet_nodes = function(reduced_GSEA_results, limma_results,
                              n_input_nodes, n_hidden_nodes){

  hidden_nodes = reduced_GSEA_results$pathway[1:n_hidden_nodes]

  reduced_GSEA_results = reduced_GSEA_results %>% filter(.data$pathway %in% hidden_nodes)
  input_nodes = sapply(reduced_GSEA_results$leadingEdge,
                       function(x) limma_results[x, ] %>%
                         arrange(desc(.data$pi_val)) %>%
                         slice(1:min(length(x), n_input_nodes)) %>%
                         row.names)

  return(list(input_nodes = input_nodes,
              hidden_nodes = hidden_nodes))

}


create_xnnet_binary_matrix = function(annotation_library, selected_xnnet_nodes){

  annotation_library = annotation_library[selected_xnnet_nodes$hidden_nodes]
  xnnet_input_nodes = as.character(selected_xnnet_nodes$input_nodes)
  annotation_library = lapply(annotation_library, function(x) intersect(x, xnnet_input_nodes))

  xnnet_binary_matrix = do.call(rbind, lapply(annotation_library, function(x) table(factor(x, levels = unique(xnnet_input_nodes)))))

  return(data.frame(t(xnnet_binary_matrix)))

}

extend_with_unassigned_inputs = function(xnnet_binary_matrix,
                                         limma_results,
                                         n_unassigned = 3) {
  #extend the network annotation DB to include unassigned inputs
  #these are top DEGS as identified by Limma which have not been
  #selected as xnnet nodes because not annotated or not
  #members of significant gene sets
  if (n_unassigned == 0) return(xnnet_binary_matrix)

  input_nodes = row.names(xnnet_binary_matrix)
  all_genes = row.names(limma_results)
  top_unassigned = all_genes[-which(all_genes %in% input_nodes)][1:n_unassigned]

  xnnet_binary_matrix$unassigned = 0
  unassigned_block = matrix(0,
                            nrow = length(top_unassigned),
                            ncol = ncol(xnnet_binary_matrix)
  )

  colnames(unassigned_block) = colnames(xnnet_binary_matrix)
  row.names(unassigned_block) = top_unassigned
  unassigned_block[, 'unassigned'] = 1
  xnnet_binary_matrix = rbind(xnnet_binary_matrix, unassigned_block)

  return(xnnet_binary_matrix)

}


generate_mask = function(xnnet_binary_matrix, n_classes = 2) {
  #this function generates a vector to mask all weights in the
  #neural network that do not express an association between
  #genes and annotation terms. These weights will be set to 0
  #during the training of the neural network.

  input_size = nrow(xnnet_binary_matrix)
  hidden_size = ncol(xnnet_binary_matrix)
  xnnet_binary_matrix_augmented = rbind(1, xnnet_binary_matrix)

  tot_mask_size = (1 + input_size) * hidden_size + (1 + hidden_size) * (n_classes - 1)
  mask = rep(1, tot_mask_size)
  xnnet_binary_matrix_augmented_flattened = as.vector(as.matrix(xnnet_binary_matrix_augmented))
  mask[1:length(xnnet_binary_matrix_augmented_flattened)] = xnnet_binary_matrix_augmented_flattened

  return(mask)

}

normalize_X = function(X, normalization = 'standard', z_cutoff = 10) {

  switch(
    normalization,

    standard = {
      X = scale(X, center = T)
      X_center = attr(X, "scaled:center")
      X_scale = attr(X, "scaled:scale")
    },

    robust = {
      X_center = apply(X, 2, median)
      X_scale = apply(X, 2, mad)
      X = scale(X, center = X_center, scale = X_scale)
    }
  )

  #threshold z-scores to deal with outliers
  X[X > z_cutoff] = z_cutoff
  X[X < -z_cutoff] = -z_cutoff

  return(list(X = X, X_center = X_center, X_scale = X_scale))
}

initialize_weights = function(mask){

  set.seed(123)
  initial_weights = rep(0, length(mask))
  initial_weights[which(mask == 1)] = rnorm(length(which(mask == 1)))
  return(initial_weights)
}

cross_validate_xnnet = function(X_train,
                                y_train,
                                xnnet_binary_matrix,
                                number = 20,
                                min_decay = 0.1,
                                max_decay = 1) {
  #finds the best model through repeated k-fold cross-validation
  #with optimization of the decay (regularization) parameter
  X_train = X_train[, row.names(xnnet_binary_matrix)]

  #grid to optimize decay parameter
  nnetGrid =  expand.grid(
    size = ncol(xnnet_binary_matrix),
    decay = seq(from = min_decay, to = max_decay, by = 0.1)
  )

  #set caret parameter for cross-validation
  fitControl = caret::trainControl(
    method = "boot",
    number = number,
    classProbs = TRUE,
    summaryFunction = caret::twoClassSummary
  )

  #additional neural network parameters
  mask = generate_mask(xnnet_binary_matrix)
  MaxNWts = length(mask)
  initial_weights = initialize_weights(mask)

  normalized_X = normalize_X(X_train)
  X_train = data.frame(normalized_X$X)

  augmented_data = augment_data(X_train, y_train, multiply = 10)
  X_train_augmented = augmented_data$augmented_X
  y_train_augmented = augmented_data$augmented_y
  y_train_augmented = factor(ifelse(y_train_augmented == 1, 'positive', 'negative'))

  X_train_center = normalized_X$X_center
  X_train_scale = normalized_X$X_scale

  nnetFit = caret::train(
    x = X_train_augmented,
    y = y_train_augmented,
    method = "nnet",
    metric = "ROC",
    MaxNWts = MaxNWts,
    maxit = 50,
    mask = mask,
    Wts = initial_weights,
    trControl = fitControl,
    tuneGrid = nnetGrid,
    trace = F
  )

  #remove training set because already saved
  #nnetFit$trainingData = NULL
  return(list(nnetFit = nnetFit,
              initial_weights = initial_weights,
              X_train_center = X_train_center,
              X_train_scale = X_train_scale))

}

compute_activation = function(xnnet_single, X, y, zscore_threshold = 3) {

  colnames(X) = make.names(colnames(X), unique = T)
  xnnet_input_nodes = xnnet_single$nnetFit$finalModel$coefnames
  X = X[, xnnet_input_nodes]
  X = scale(X, center = T)

  #threshold to handle potential outliers (treshold is the same as in training set)
  X[X > zscore_threshold] = zscore_threshold
  X[X < -zscore_threshold] = -zscore_threshold

  weights = xnnet_single$nnetFit$finalModel$wts
  xnnet_binary_matrix = xnnet_single$xnnet_binary_matrix

  mask = xnnet_single$nnetFit$finalModel$param$mask
  decoded_weights = decode_weights(weights, mask, xnnet_binary_matrix)
  hidden_activation = data.frame(sigmoid(cbind(1, X) %*% decoded_weights$first_layer))
  hidden_activation$sample = row.names(X)
  hidden_activation$class = factor(y)

  return(hidden_activation)
}

#' @importFrom  fMultivar mvFit
augment_data = function(X, y, multiply = 2){

  if (multiply == 0) return(list(augmented_X = X, augmented_y = y))

  samples_0 = sum(y == 0)
  fit_0 = fMultivar::mvFit(X[which(y == 0), ], method = 'st')
  synthetic_data_0 = fMultivar::rmvst(samples_0*multiply,
                                      dim = ncol(X),
                                      mu = as.numeric(fit_0@fit$estimated$beta),
                                      Omega = fit_0@fit$estimated$Omega,
                                      alpha = as.numeric(fit_0@fit$estimated$alpha))

  samples_1 = sum(y == 1)
  fit_1 = fMultivar::mvFit(X[which(y == 1), ], method = 'st')
  synthetic_data_1 = fMultivar::rmvst(samples_1*multiply,
                                      dim = ncol(X),
                                      mu = as.numeric(fit_1@fit$estimated$beta),
                                      Omega = fit_1@fit$estimated$Omega,
                                      alpha = as.numeric(fit_1@fit$estimated$alpha))

  synthetic_data = rbind(synthetic_data_0, synthetic_data_1)
  colnames(synthetic_data) = colnames(X)

  augmented_X =rbind(X, synthetic_data)
  augmented_y = c(y, rep(0, samples_0*multiply), rep(1, samples_1*multiply))

  return(list(augmented_X = augmented_X, augmented_y = augmented_y))

}




#' Predict sample labels using xnnet
#'
#' @param xnnet xxnet object return by the function build_xnnet
#' @param X_test test data
#' @export

threshold_scaled_data = function(scaled_data, z_threshold){
  scaled_data[scaled_data > z_threshold] = z_threshold
  scaled_data[scaled_data < -z_threshold] = -z_threshold
  return(scaled_data)
}

preprocess_X_test = function(xnnet_single, X_test, z_threshold){

  xnnet_input_nodes = colnames(xnnet_single$nnetFit$trainingData)
  xnnet_input_nodes = xnnet_input_nodes[-length(xnnet_input_nodes)]

  if (is.null(dim(X_test))){

    X_test = data.frame(rbind(X_test, X_test), row.names = NULL)
    X_test = X_test[, xnnet_input_nodes]
    X_test = scale(X_test, center = xnnet_single$X_train_center,
                   scale = xnnet_single$X_train_scale)
    X_test = threshold_scaled_data(X_test, z_threshold) %>% data.frame
    X_test$is_single = T

    } else {

      X_test = X_test[, xnnet_input_nodes]
      X_test = scale(X_test, center = xnnet_single$X_train_center,
                     scale = xnnet_single$X_train_scale)
      X_test = threshold_scaled_data(X_test, z_threshold) %>% data.frame
      X_test$is_single = F
    }

  return(X_test)
  }

xnnet_predict_single = function(xnnet_single, X_test, z_threshold) {

  preprocessed_X_test = preprocess_X_test(xnnet_single, X_test, z_threshold)
  xnnet_predictions_probs = caret::predict.train(xnnet_single$nnetFit, preprocessed_X_test, type = 'prob')[, 'positive']
  if(preprocessed_X_test$is_single[1] == T) xnnet_predictions_probs = xnnet_predictions_probs[1]

  return(xnnet_predictions_probs)
}

xnnet_predict = function(xnnet, X_test, z_threshold = 3){
  xnnet_predictions = lapply(xnnet, function(x) xnnet_predict_single(x, X_test, z_threshold = z_threshold))
  return(xnnet_predictions)
}

#' Assessing performance of xnnet predictions given labels
#'
#' This function plots xnnet as a heatmap with genes as rows and samples as columns.
#' Genes are grouped by the corresponding annotation and rows are grouped by sample class.
#' labels into a training and test set
#' @param xnnet_predictions output of build_xnnet
#' @param true_labels hidden node index
#' @import ggplot2
#' @importFrom pROC auc
#' @export
assess_xnnet_performance = function(xnnet, xnnet_predictions, true_labels){

  #pull out AUC on training set from xnnet object
  xnnet_train_AUC = sapply(xnnet, function(x) pROC::auc(x$nnetFit$trainingData %>% dplyr::select(.data$.outcome) %>% pull(),
                                                        x$nnetFit$finalModel$fitted.values %>% as.numeric, direction="<",
                                                        levels= c('negative', 'positive')))

  #compute AUC on set
  xnnet_test_AUC = sapply(xnnet_predictions, function(x) pROC::auc(true_labels, x, direction="<",
                                                                   levels= c('0', '1')))

  xnnet_interpretability =sapply(xnnet, function(x) x$interpretability_score)

  xnnet_performance = data.frame(annotation_library = names(xnnet_test_AUC),
                                 train_AUC = xnnet_train_AUC,
                                 test_AUC = xnnet_test_AUC,
                                 interpretability = xnnet_interpretability)
  row.names(xnnet_performance) = NULL
  AUC_plot = ggplot(xnnet_performance, aes(x = reorder(.data$annotation_library, .data$test_AUC), y = .data$test_AUC, .data$label)) +
    geom_col(col = 'white') + coord_flip() + theme_Publication() + xlab('') + ylab('test AUC')

  AUC_interpretability_plot = ggplot(xnnet_performance, aes(x = .data$test_AUC, y = .data$interpretability, label = .data$annotation_library)) +
    geom_point(size = 4) + theme_Publication() + xlab('test AUC') + ylab('interpretability') + ggrepel::geom_text_repel(size=6)


  return(list(xnnet_performance = xnnet_performance,
              AUC_plot = AUC_plot,
              AUC_interpretability_plot = AUC_interpretability_plot))

}


sigmoid = function(x) {
  return(1 / (1 + exp(-x)))
}
#threshold to handle potential outliers (treshold is the same as in training set)
#X[X > zscore_threshold] = zscore_threshold
#X[X < -zscore_threshold] = -zscore_threshold
decode_weights = function(W, mask, xnnet_binary_matrix, n_classes = 2) {
  #extract coefficients from results of neural network

  first_layer_size = (1 + nrow(xnnet_binary_matrix)) * ncol(xnnet_binary_matrix)

  W_first_layer = W[1:first_layer_size]
  W_first_layer = matrix(W_first_layer, ncol = ncol(xnnet_binary_matrix))
  row.names(W_first_layer) = c('b', row.names(xnnet_binary_matrix))
  colnames(W_first_layer) = colnames(xnnet_binary_matrix)

  W_second_layer = W[-(1:first_layer_size)]
  W_second_layer = matrix(W_second_layer, ncol = n_classes - 1)
  row.names(W_second_layer) = c('b1', colnames(xnnet_binary_matrix))

  return(list(first_layer = W_first_layer,
              second_layer = W_second_layer))
}




