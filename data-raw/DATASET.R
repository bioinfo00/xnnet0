## code to prepare `DATASET` dataset goes here

#for the expression dataset
GSE37250_processed = readRDS('../../xnnet/neural_nets/gse37250_processed')
X = GSE37250_processed$X #samples x genes
y = GSE37250_processed$y$infected #sample labels
# sample_idx = sample(1:nrow(X), 50)
# X = X[sample_idx, ]
# y = y[sample_idx]
# GSE37250 = list(X = X, y = y)
# usethis::use_data(GSE37250, compress = 'xz', overwrite = T)
#
# #for the annotation library
# db_gene_sets = readRDS("../../xnnet/neural_nets/db_gene_sets")
# annotation_libraries = lapply(db_gene_sets, function(x) x$db_gene_sets_list)
# annotation_libraries = annotation_libraries[c(1, 9, 18, 22)]
#
# #remove term containing non-ASCII characters (will throw a warning during check)
# annotation_libraries$Reactome_2016$`Loss of proteins required for interphase microtubule organization from the centrosome_Homo sapiens_R-HSA-380284` = NULL
# usethis::use_data(annotation_libraries, compress = 'xz', overwrite = T)

