# #split data
# GSE37250_split = train_test_split(GSE37250$X, GSE37250$y)
#
# #network cross-validation
# xnnet = build_xnnet(X_train = GSE37250_split$X_train,
#                     y_train = GSE37250_split$y_train,
#                     annotation_libraries = annotation_libraries,
#                     n_input_nodes = 3, n_hidden_nodes = 3)
# #
# # #predictions on test set
# xnnet_predictions = xnnet_predict(xnnet, X_test = GSE37250_split$X_test)
#
# # #assess model performance
# xnnet_performance = assess_xnnet_performance(xnnet, xnnet_predictions, true_labels = GSE37250_split$y_test)
#
# xnnet_performance = assess_xnnet_performance(xnnet, xnnet_predictions, true_labels = GSE37250_split$y_test, metric = 'accuracy')
#
#
#
# assess_xnnet_performance = function(xnnet, xnnet_predictions, true_labels, metric = 'AUC'){
#
#   training_true_labels = xnnet[[1]]$nnetFit$trainingData %>% dplyr::select(.data$.outcome)
#   end_of_real_data = which(row.names(training_true_labels) == 1) - 1
#   training_true_labels = training_true_labels %>% pull()
#   training_true_labels = training_true_labels[1:end_of_real_data]
#
#   if (match.arg(metric, c('AUC', 'accuracy')) == 'AUC'){
#
#     #pull out AUC on training set from xnnet object
#     train_performance = sapply(xnnet, function(x) pROC::auc(training_true_labels,
#                                                                   x$nnetFit$finalModel$fitted.values[1:end_of_real_data] %>% as.numeric, direction="<",
#                                                                   levels= c('negative', 'positive')))
#
#     #compute AUC on set
#     test_performance = sapply(xnnet_predictions, function(x) pROC::auc(true_labels, x, direction="<",
#                                                                             levels= c('0', '1')))
#
#   } else {
#
#     train_predicted_label = lapply(xnnet, function(x) ifelse(x$nnetFit$finalModel$fitted.values > 0.5, 'positive', 'negative'))
#     train_performance = sapply(train_predicted_label, function(x) sum(x[1:end_of_real_data] == training_true_labels)/length(training_true_labels))
#
#     test_predicted_label = lapply(xnnet_predictions, function(x) ifelse(x > 0.5, 1, 0))
#     test_performance = sapply(test_predicted_label, function(x) sum(x == true_labels)/length(true_labels))
#
#   }
#
#   interpretability = sapply(xnnet, function(x) x$interpretability_score)
#
#   global_performance = data.frame(
#     annotation_library = names(xnnet),
#     train_performance = train_performance,
#     test_performance = test_performance,
#     interpretability = interpretability,
#     metric = metric
#   )
#   row.names(global_performance) = NULL
#
#   performance_plot = ggplot(global_performance, aes(x = reorder(.data$annotation_library, .data$test_performance), y = .data$test_performance, .data$label)) +
#     geom_col(col = 'white') + coord_flip() + theme_Publication() + xlab('') + ylab(paste(metric, 'on test'))
#
#   performance_interpretability_plot = ggplot(global_performance, aes(x = .data$test_performance, y = .data$interpretability, label = .data$annotation_library)) +
#     geom_point(size = 4) + theme_Publication() + xlab(paste(metric, 'on test')) + ylab('interpretability') + ggrepel::geom_text_repel(size=6)
#
#
#   return(list(global_performance = global_performance,
#               performance_plot = performance_plot,
#               performance_interpretability_plot = performance_interpretability_plot))
#
# }




# # # #plotting network
# # # plot_xnnet(xnnet$Reactome_2016)
# # #
# # # #plot mean activation state in the two groups
# sample_activation = compute_hidden_activation(xnnet$BioCarta_2016,
#                                        X = GSE37250_split$X_train,
#                                        y = GSE37250_split$y_train)
#
# sample_mean_activation = sample_activation %>% dplyr::select(-sample) %>%
#   dplyr::group_by(class) %>% dplyr::summarise_all(mean)
#
# library(ggradar)
# ggradar(sample_mean_activation)


# X = X[, row.names(xnnet$BioCarta_2016$xnnet_binary_matrix)]
#
#
# augmented_data = augment_data(X, y, 10)
# pr = prcomp(augmented_data$augmented_X, center = T, scale. = T)
# pca_df = data.frame(pc1 = pr$x[, 1], pc2 = pr$x[, 2], col = augmented_data$augmented_y)
# ggplot(pca_df, aes(x = pc1, y = pc2, col = factor(col))) + geom_point()
#
# augment_data = function(X, y, multiply = 2){
#
#   if (multiply == 0) return(list(augmented_X = X, augmented_y = y))
#
#   samples_0 = sum(y == 0)
#   fit_0 = fMultivar::mvFit(X[which(y == 0), ], method = 'st')
#   synthetic_data_0 = fMultivar::rmvst(samples_0*multiply,
#                                     dim = ncol(X),
#                                     mu = as.numeric(fit_0@fit$estimated$beta),
#                                     Omega = fit_0@fit$estimated$Omega,
#                                     alpha = as.numeric(fit_0@fit$estimated$alpha))
#
#   samples_1 = sum(y == 1)
#   fit_1 = fMultivar::mvFit(X[which(y == 1), ], method = 'st')
#   synthetic_data_1 = fMultivar::rmvst(samples_1*multiply,
#                                       dim = ncol(X),
#                                       mu = as.numeric(fit_1@fit$estimated$beta),
#                                       Omega = fit_1@fit$estimated$Omega,
#                                       alpha = as.numeric(fit_1@fit$estimated$alpha))
#
#   synthetic_data = rbind(synthetic_data_0, synthetic_data_1)
#   colnames(synthetic_data) = colnames(X)
#
#   augmented_X =rbind(X, synthetic_data)
#   augmented_y = c(y, rep(0, samples_0*multiply), rep(1, samples_1*multiply))
#
#   return(list(augmented_X = augmented_X, augmented_y = augmented_y))
#
#   }




# X = GSE37250_split$X_test
# y = GSE37250_split$y_test
#
# bionnet = xnnet$BioCarta_2016
#
# xnnet_input_nodes =bionnet$nnetFit$finalModel$coefnames
# X = X[, xnnet_input_nodes]
# X = scale(X, center = T)
#
# weights = bionnet$nnetFit$finalModel$wts
# xnnet_binary_matrix = bionnet$xnnet_binary_matrix
#
# mask = bionnet$nnetFit$finalModel$param$mask
# decoded_weights = decode_weights(weights, mask, xnnet_binary_matrix)
# activation_scores = data.frame(sigmoid(cbind(1, X) %*% decoded_weights$first_layer))
# activation_scores$sample = row.names(X)
# activation_scores$class = factor(y)
#
# hidden_activation = lapply(xnnet, function(x) compute_hidden_activation(x, X, y))
#
# compute_hidden_activation(xnnet$BioCarta_2016, X, y)
#
# compute_hidden_activation = function(xnnet_single, X, y, zscore_threshold = 3) {
#
#   colnames(X) = make.names(colnames(X), unique = T)
#   xnnet_input_nodes = xnnet_single$nnetFit$finalModel$coefnames
#   X = X[, xnnet_input_nodes]
#   X = scale(X, center = T)
#
#   #threshold to handle potential outliers (treshold is the same as in training set)
#   #X[X > zscore_threshold] = zscore_threshold
#   #X[X < -zscore_threshold] = -zscore_threshold
#
#   weights = xnnet_single$nnetFit$finalModel$wts
#   xnnet_binary_matrix = xnnet_single$xnnet_binary_matrix
#
#   mask = xnnet_single$nnetFit$finalModel$param$mask
#   decoded_weights = decode_weights(weights, mask, xnnet_binary_matrix)
#   hidden_activation = data.frame(sigmoid(cbind(1, X) %*% decoded_weights$first_layer))
#   hidden_activation$sample = row.names(X)
#   hidden_activation$class = factor(y)
#
#   return(hidden_activation)
# }
#
#
#
#
#
# sigmoid = function(x) {
#   return(1 / (1 + exp(-x)))
# }
# #threshold to handle potential outliers (treshold is the same as in training set)
# #X[X > zscore_threshold] = zscore_threshold
# #X[X < -zscore_threshold] = -zscore_threshold
# decode_weights = function(W, mask, xnnet_binary_matrix, n_classes = 2) {
#   #extract coefficients from results of neural network
#
#   first_layer_size = (1 + nrow(xnnet_binary_matrix)) * ncol(xnnet_binary_matrix)
#
#   W_first_layer = W[1:first_layer_size]
#   W_first_layer = matrix(W_first_layer, ncol = ncol(xnnet_binary_matrix))
#   row.names(W_first_layer) = c('b', row.names(xnnet_binary_matrix))
#   colnames(W_first_layer) = colnames(xnnet_binary_matrix)
#
#   W_second_layer = W[-(1:first_layer_size)]
#   W_second_layer = matrix(W_second_layer, ncol = n_classes - 1)
#   row.names(W_second_layer) = c('b1', colnames(xnnet_binary_matrix))
#
#   return(list(first_layer = W_first_layer,
#               second_layer = W_second_layer))
# }
#
# compute_sample_activation_scores = function(bionnet, X, y, zscore_threshold = 3) {
#   bionnet_input_nodes = colnames(bionnet$X_train)
#   X = X[, bionnet_input_nodes]
#
#   #scale test set with respect to training set
#   X = scale(X, bionnet$mu_train, bionnet$scale_train)
#
#   #threshold to handle potential outliers (treshold is the same as in training set)
#   X[X > zscore_threshold] = zscore_threshold
#   X[X < -zscore_threshold] = -zscore_threshold
#
#   weights = bionnet$nnetFit$finalModel$wts
#   xnnet_binary_matrix = bionnet$xnnet_binary_matrix
#   mask = bionnet$nnetFit$finalModel$param$mask
#   decoded_weights = decode_weights(weights, mask, xnnet_binary_matrix, 2)
#   activation_scores = data.frame(sigmoid(cbind(1, X) %*% decoded_weights$first_layer))
#   activation_scores$sample = row.names(X)
#   activation_scores$class = factor(y)
#
#   return(activation_scores)
# }







#evaluate performance on training and test set
# xnnet_predictions = lapply(xnnet, function(x) xnnet_predict(x, GSE37250_split$X_test))
#
# xnnet$BioCarta_2016$nnetFit$finalModel$fitted.values
#
# xnnet_interpretability =sapply(xnnet, function(x) x$interpretability_score)
# xnnet_train_AUC = sapply(xnnet, function(x) pROC::auc(x$nnetFit$trainingData %>% dplyr::select(.outcome) %>% pull,
#                                                       x$nnetFit$finalModel$fitted.values %>% as.numeric))
#
# xnnet_test_AUC = sapply(xnnet_predictions, function(x) pROC::auc(y_test, x, ))
#
# xnnet_performance = data.frame(annotation_library = names(xnnet_test_AUC),
#                                train_AUC = xnnet_train_AUC,
#                                test_AUC = xnnet_test_AUC,
#                                interpretability = xnnet_interpretability)
# row.names(xnnet_performance) = NULL
# ggplot(xnnet_performance, aes(x = reorder(annotation_library, test_AUC), y = test_AUC, label)) +
#   geom_col(col = 'white') + coord_flip() + theme_Publication() + xlab('') + ylab('test AUC')
#
# ggplot(xnnet_performance, aes(x = test_AUC, y = interpretability, label = annotation_library)) +
#   geom_point(size = 4) + theme_Publication() + xlab('test AUC') + ylab('interpretability') + ggrepel::geom_label_repel(size=6)
#
#
# xnnet_performance = assess_xnnet_performance(xnnet_predictions, y_test)
#
# assess_xnnet_performance = function(xnnet_predictions, y_test){
#
#   #pull out AUC on training set from xnnet object
#   xnnet_train_AUC = sapply(xnnet, function(x) pROC::auc(x$nnetFit$trainingData %>% dplyr::select(.outcome) %>% pull,
#                                                         x$nnetFit$finalModel$fitted.values %>% as.numeric, direction="<",
#                                                         levels= c('negative', 'positive')))
#
#   #compute AUC on set
#   xnnet_test_AUC = sapply(xnnet_predictions, function(x) pROC::auc(y_test, x, direction="<",
#                                                                    levels= c('0', '1')))
#
#   xnnet_interpretability =sapply(xnnet, function(x) x$interpretability_score)
#
#   xnnet_performance = data.frame(annotation_library = names(xnnet_test_AUC),
#                                  train_AUC = xnnet_train_AUC,
#                                  test_AUC = xnnet_test_AUC,
#                                  interpretability = xnnet_interpretability)
#   row.names(xnnet_performance) = NULL
#   AUC_plot = ggplot(xnnet_performance, aes(x = reorder(annotation_library, test_AUC), y = test_AUC, label)) +
#     geom_col(col = 'white') + coord_flip() + theme_Publication() + xlab('') + ylab('test AUC')
#
#   AUC_interpretability_plot = ggplot(xnnet_performance, aes(x = test_AUC, y = interpretability, label = annotation_library)) +
#     geom_point(size = 4) + theme_Publication() + xlab('test AUC') + ylab('interpretability') + ggrepel::geom_text_repel(size=6)
#
#
#   return(list(xnnet_performance = xnnet_performance,
#               AUC_plot = AUC_plot,
#               AUC_interpretability_plot = AUC_interpretability_plot))
#
# }





# # GSE37250_split = train_test_split(GSE37250$X, GSE37250$y)
# xnnet = build_xnnet(GSE37250_split$X_train, GSE37250_split$y_train, annotation_libraries,
#                     n_input_nodes = 3, n_hidden_nodes = 3)
#
# training_data = xnnet$BioCarta_2016$nnetFit$trainingData
#
# training_data_split = split(training_data, training_data$.outcome)
#
# augmented_data = lapply(training_data_split,
#                         function(x) augment_data(x))
#
# d = ldply(augmented_data, rbind)
#
# X_train = X_train[, -1]
# y_train = ifelse(d$.outcome == 'positive', 1, 0)
#
# xnnet = build_xnnet(X_train, y_train, annotation_libraries,
#                     n_input_nodes = 3, n_hidden_nodes = 3)
#
#
#
#
# augment_data = function(X, multiply = 1){
#
#   n_samples = nrow(X)
#   class = X[, '.outcome']
#   X_red = X[, 1:(ncol(X) -1)]
#   fit = fMultivar::mvFit(X_red, method = 'st')
#   synthetic_data = fMultivar::rmvst(n_samples*multiply,
#                                     dim = 12,
#                                     mu = as.numeric(fit@fit$dp.complete$beta),
#                                     Omega = fit@fit$dp.complete$Omega,
#                                     alpha = as.numeric(fit@fit$dp.complete$alpha))
#
#   synthetic_data = data.frame(synthetic_data, .outcome = rep(class, multiply))
#   colnames(synthetic_data) = colnames(X)
#   X = rbind(X, synthetic_data)
#   return(X)
# }
#
#
#
#
# # # get_GSEA_results = function(limma_results,
# # #                             annotation_library,
# # #                             min_term_size = 10,
# # #                             max_term_size = Inf,
# # #                             rank_metric = 'logFC'){
# # #   #top_n_gene_sets = length(annotation_library)) {
# # #   #performs GSEA to identify significant annotation library
# # #   #annotation_database is a list of gene sets (e.g. Reactome)
# # #   #top_n_gene_sets is the number of top gene sets retained for downstream analysis
# # #
# # #   top_n_gene_sets = length(annotation_library)
# # #   #check overlap between limma output genes and annotation genes
# # #   measured_genes = row.names(limma_results)
# # #   annotated_genes = unique(unlist(annotation_library))
# # #   overlap = length(intersect(measured_genes, annotated_genes))
# # #
# # #   testthat::expect_true(overlap > 100, info = 'less than 100 measured genes annotated in annotation library')
# # #
# # #   print(paste('#######', overlap, 'out of ', length(measured_genes), 'measured genes are annotated in annotation library', '#######'))
# # #   switch(rank_metric,
# # #          logFC = {ranks = limma_results$logFC},
# # #          pval = {ranks = limma_results$P.Value},
# # #          pi_val = {ranks = limma_results$pi_val}
# # #   )
# # #   names(ranks) = row.names(limma_results)
# # #
# # #   GSEA_results = fgsea::fgsea(
# # #     annotation_library,
# # #     ranks,
# # #     nperm = 1000,
# # #     minSize = min_term_size,
# # #     maxSize = max_term_size
# # #   )
# # #
# # #   n_significant_gene_sets = sum(GSEA_results$padj < 0.1)
# # #   print(paste('#######', n_significant_gene_sets, 'significance gene sets detected at adj pval < 0.1', '#######'))
# # #   if (n_significant_gene_sets == 0)
# # #     warning(paste('#######', "no significant genes found at adj pval < 0.1"), '#######')
# # #
# # #   print(paste('####### retaining top', top_n_gene_sets, 'gene sets for downstream analysis #######'))
# # #   GSEA_results = GSEA_results %>%  arrange(desc(abs(NES))) %>% dplyr::slice(1:top_n_gene_sets)
# # #
# # #   return(GSEA_results)
# # #
# # # }
# # #


# # # split_data = train_test_split(GSE37250$X, GSE37250$y)
# # # X_train = split_data$X_train
# # # y_train = split_data$y_train
# # #
# # # limma_results = get_limma_results(X_train, y_train)
# # # gsea_results =get_GSEA_results(limma_results, reactome)
# # # gsea_results = normalize_NES(gsea_results)
# # #
# # #
# # #
# # # normalize_NES = function(GSEA_results){
# # #
# # #   x = log10(GSEA_results$size)
# # #   y = abs(GSEA_results$NES)
# # #   model <- lm(y ~ poly(x,3))
# # #   predicted.intervals <- predict(model, data.frame(x = x),
# # #                                  interval='confidence',
# # #                                  level=0.95)
# # #   fit_df = data.frame(x = x, y = y, upr_y_fit = predicted.intervals[,3])
# # #   normalized_NES = (fit_df$y - fit_df$upr_y_fit)
# # #   GSEA_results$normalized_NES = normalized_NES
# # #   print(head(GSEA_results))
# # #   GSEA_results = GSEA_results %>% dplyr::filter(normalized_NES > 0)
# # #
# # #   return(GSEA_results)
# # #
# # # }
# # #
# # # reduce_GSEA_results(gsea_results)
# # #
# # #
# # #
# # # # xnnet = build_xnnet(X_train = X_train,
# # # #                     y_train = y_train,
# # # #                     annotation_library = reactome)
# # # rescale_GSEA_results = function(GSEA_results, top_n_gene_sets = 30){
# # #
# # #   #rescales GSEA results taking into account the size of the different gene sets
# # #   x = gsea_results$size
# # #   y = abs(gsea_results$NES)
# # #   model <- lm(y ~ poly(x,3))
# # #
# # #
# # #
# # # }
# # #
# # #
# #
# #
#db_gene_sets = readRDS("../../xnnet/neural_nets/db_gene_sets")
# # annotation_libraries = lapply(db_gene_sets, function(x) x$db_gene_sets_list)
# # annotation_libraries = annotation_libraries[-c(5, 13, 14)]
# # GSE37250_split = train_test_split(GSE37250$X, GSE37250$y)
# # X_train = GSE37250_split$X_train
# # y_train = GSE37250_split$y_train
# #
# #
# # limma_results = get_limma_results(GSE37250_split$X_train,
# #                                   GSE37250_split$y_train)
# # GSEA_results = get_and_process_GSEA_results(annotation_libraries,
# #                                             limma_results)
# #
# # single_xnnet = build_single_xnnet(annotation_libraries, GSEA_results, limma_results)
# #
# # xnnet = lapply(GSEA_results, function(x) build_single_xnnet(annotation_libraries, x, limma_results))
# #
# #
# # build_single_xnnet = function(annotation_libraries, GSEA_results, limma_results,
# #                               n_input_nodes = 5,
# #                               n_hidden_nodes = 5,
# #                               n_unassigned = 3,
# #                               number = 10,
# #                               min_decay = 0.1,
# #                               max_decay = 1){
# #
# #   selected_xnnet_nodes = select_xnnet_nodes(GSEA_results, limma_results,
# #                                             n_input_nodes = n_input_nodes,
# #                                             n_hidden_nodes = n_hidden_nodes)
# #
# #   hidden_renormalized_NES = GSEA_results[which(GSEA_results$pathway %in% selected_xnnet_nodes$hidden_nodes), 'renormalized_NES']
# #   interpretability_score = mean(hidden_renormalized_NES)
# #
# #   xnnet_binary_matrix = create_xnnet_binary_matrix(annotation_library = annotation_libraries[[unique(GSEA_results$.id)]],
# #                                                    selected_xnnet_nodes = selected_xnnet_nodes)
# #   xnnet_binary_matrix = extend_with_unassigned_inputs(xnnet_binary_matrix, limma_results = limma_results, n_unassigned = n_unassigned)
# #
# #   nnetFit_output = cross_validate_xnnet(
# #     X_train,
# #     y_train,
# #     xnnet_binary_matrix,
# #     number = number,
# #     min_decay = min_decay,
# #     max_decay = max_decay
# #   )
# #
# #   initial_wt = nnetFit_output$initial_weights
# #   X_train_center = nnetFit_output$X_train_center
# #   X_train_scale = nnetFit_output$X_train_scale
# #   nnetFit = nnetFit_output$nnetFit
# #
# #   xnnet_results = list(
# #     interpretability_score = interpretability_score,
# #     GSEA_results = GSEA_results,
# #     limma_results = limma_results,
# #     nnetFit = nnetFit,
# #     initial_wt = initial_wt,
# #     xnnet_binary_matrix = xnnet_binary_matrix,
# #     X_train_center = X_train_center,
# #     X_train_scale = X_train_scale
# #   )
# #   return(xnnet_results)
# # }
# #
# #
# #
# # get_and_process_GSEA_results = function(annotation_libraries,
# #                                         limma_results){
# #
# #   GSEA_results = lapply(annotation_libraries, function(x)
# #     get_GSEA_results(limma_results, x))
# #
# #   GSEA_results = renormalize_NES(GSEA_results)
# #   GSEA_results = lapply(GSEA_results, function(x) reduce_GSEA_results(x))
# #
# #   return(GSEA_results)
# #
# # }
# #
# #
# #
# #
# #
# #
# #
# #
# #
# # selected_xnnet_nodes = select_xnnet_nodes(reduced_GSEA_results, limma_results,
# #                                           n_input_nodes = n_input_nodes,
# #                                           n_hidden_nodes = n_hidden_nodes)
# #
# # hidden_normalized_NES = reduced_GSEA_results[which(reduced_GSEA_results$pathway %in% selected_xnnet_nodes$hidden_nodes), 'normalized_NES']
# #
# # interpretability_score = reduced_GSEA_results %>% dplyr::filter(pathway %in% selected_xnnet_nodes$hidden_nodes) %>%
# #   dplyr::select(normalized_NES) %>% dplyr::pull() %>% mean
# #
# # xnnet_binary_matrix = create_xnnet_binary_matrix(annotation_library = annotation_library,
# #                                                  selected_xnnet_nodes = selected_xnnet_nodes)
# # xnnet_binary_matrix = extend_with_unassigned_inputs(xnnet_binary_matrix, limma_results = limma_results, n_unassigned = n_unassigned)
# #
# # nnetFit_output = cross_validate_xnnet(
# #   X_train,
# #   y_train,
# #   xnnet_binary_matrix,
# #   number = number,
# #   min_decay = min_decay,
# #   max_decay = max_decay
# # )
# #
# # initial_wt = nnetFit_output$initial_weights
# # X_train_center = nnetFit_output$X_train_center
# # X_train_scale = nnetFit_output$X_train_scale
# # nnetFit = nnetFit_output$nnetFit
# #
# # xnnet_results = list(
# #   interpretability_score = interpretability_score,
# #   GSEA_results = GSEA_results,
# #   limma_results = limma_results,
# #   nnetFit = nnetFit,
# #   initial_wt = initial_wt,
# #   xnnet_binary_matrix = xnnet_binary_matrix,
# #   X_train_center = X_train_center,
# #   X_train_scale = X_train_scale
# # )
# # return(xnnet_results)
# # }
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #


# db_gene_sets = readRDS("../../xnnet/neural_nets/db_gene_sets")
# annotation_libraries = lapply(db_gene_sets, function(x) x$db_gene_sets_list)
# #annotation_libraries = annotation_libraries[-c(5, 13, 14)]
# GSE37250_split = train_test_split(GSE37250$X, GSE37250$y)
# X_train = GSE37250_split$X_train
# y_train = GSE37250_split$y_train
#
# limma_results = get_limma_results(GSE37250_split$X_train,
#                                   GSE37250_split$y_train)
# GSEA_results = lapply(annotation_libraries, function(x)
#                       get_GSEA_results(limma_results, x))
#
# all_GSEA_results = plyr::ldply(GSEA_results, rbind)
# ggplot(all_GSEA_results, aes(x = .id, y = abs(NES))) +
#   geom_boxplot() + coord_flip()
#
# ggplot(all_GSEA_results, aes(x = log10(size), y = abs(NES), text = pathway)) +
#   geom_point(aes(col = .id))
# plotly::ggplotly()
#
#
# limma_results = get_limma_results(GSE37250_split$X_train,
#                                   GSE37250_split$y_train)
#
# GSEA_results = lapply(annotation_libraries, function(x)
#   get_GSEA_results(limma_results, x))
#
# all_GSEA_results = plyr::ldply(GSEA_results, rbind)
# ggplot(all_GSEA_results, aes(x = .id, y = abs(NES))) +
#   geom_boxplot() + coord_flip()
#
# ggplot(all_GSEA_results, aes(x = log10(size), y = abs(NES), text = pathway)) +
#   geom_point(aes(col = .id))
#
#
# GSEA_results = renormalize_NES(GSEA_results)
# all_GSEA_results = plyr::ldply(GSEA_results, rbind)
# ggplot(all_GSEA_results, aes(x = .id, y = abs(renormalized_NES))) +
#   geom_boxplot() + coord_flip() + theme_Publication() + theme(legend.position = "none")
#
#
# ggplot(all_GSEA_results, aes(x = log10(size), y = abs(renormalized_NES), text = pathway)) +
#   geom_point(aes(col = .id)) + theme_Publication() + theme(legend.position = "none")
#
# ggplot(all_GSEA_results, aes(x = log10(size), y = abs(NES), text = pathway)) +
#   geom_point(aes(col = .id)) + theme_Publication() + theme(legend.position = "none")
#
#
# plotly::ggplotly()


# # ggplot(all_GSEA_results, aes(x = log10(size), y = abs(NES), text = pathway)) +
# #   geom_point(aes(col = .id))
# # plotly::ggplotly()
# #
# # model = lm(abs(NES) ~ log10(size) + .id, data = all_GSEA_results)
# #
# # model = lm(abs(NES) ~ log10(size), data = all_GSEA_results)
# #
# # all_GSEA_results$dev = abs(all_GSEA_results$NES) - model$fitted.values
# #
# #
# # ggplot(all_GSEA_results, aes(x = .id, y = dev)) +
# #   geom_boxplot() + coord_flip()
# #
# #
# # ggplot(all_GSEA_results, aes(x = log10(size), y = dev, text = pathway)) +
# #   geom_point(aes(col = .id))
# # plotly::ggplotly()
# #
# #
# renormalize_NES = function(GSEA_results){
#
#   all_GSEA_results = plyr::ldply(GSEA_results, rbind)
#
#   model = lm(abs(NES) ~ log10(size) + .id, data = all_GSEA_results)
#   all_GSEA_results$renormalized_NES = abs(all_GSEA_results$NES) - model$fitted.values
#
#   all_GSEA_results = all_GSEA_results %>% dplyr::filter(renormalized_NES > 0) %>%
#     dplyr::arrange(desc(renormalized_NES))
#
#   GSEA_results = split(all_GSEA_results, all_GSEA_results$.id)
#
#   return(GSEA_results)
#
# }
# #
# #
# #
