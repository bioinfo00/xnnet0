data("GSE37250")
GSE37250_split = train_test_split(GSE37250$X, GSE37250$y)

context("testing xnnet helper functions")

test_that("testing Limma output",
         {
           limma_results = get_limma_results(GSE37250_split$X_train,
                                             GSE37250_split$y_train)

           expect_equal(nrow(limma_results), ncol(GSE37250_split$X_train))

         })
