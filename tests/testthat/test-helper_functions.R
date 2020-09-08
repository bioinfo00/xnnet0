context("testing weight initialization")

test_that("zeros in mask and initial weight
         vector are in the same position",
         {
           mask = sample(0:1, 10, replace = T)
           zeros_in_mask = which(mask == 0)

           initial_weights = initialize_weights(mask = mask)
           zeros_in_initial_weight_vector = which(initial_weights == 0)

           expect_equal(zeros_in_mask, zeros_in_initial_weight_vector)
         })
