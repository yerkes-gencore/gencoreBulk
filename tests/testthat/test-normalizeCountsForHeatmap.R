## DATA INIT

metadata <- expand.grid(
  individual = c(1,2,3),
  timepoint = c('D0', 'D2', 'D7')
)
metadata$sample <- paste0('S', metadata$individual, '_', metadata$timepoint)

counts <- matrix(data = c(
  0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
  7.13, 6.86, 6.57, 7.18, 6.80, 6.83, 6.80, 6.86, 7.51,
  0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
  7.53, 7.49, 7.40, 7.54, 7.64, 7.26, 7.37, 7.22, 7.13,
  6.88, 6.81, 6.87, 6.99, 7.00, 6.71, 7.17, 6.36, 6.34,
  2.35, 1.33, 1.61, 2.61, 1.32, 2.13, 1.99, 1.55, 1.72),
  nrow = 6, ncol=9, byrow=TRUE,
  dimnames = list(
    c('PPP1R3A', 'ENSM6801', 'ENSM8088', 'BMT2', 'TMEM168', 'LSMEM1'),
    c(metadata$sample))
)

## MEDIAN BASELINE

test_that("median baseline", {
  out <- normalizeCountsForHeatmap(counts, metadata,
                                   group_var = 'timepoint', baseline='D0')
  
  answer <- matrix(data = c(
    0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,
    0.27,  0.00, -0.29,  0.32, -0.06, -0.03, -0.06,  0.00,  0.65,
    0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,
    0.04,  0.00, -0.09,  0.05,  0.15, -0.23, -0.12, -0.27, -0.36,
    0.01, -0.06,  0.00,  0.12,  0.13, -0.16,  0.30, -0.51, -0.53,
    0.74, -0.28,  0.00,  1.00, -0.29,  0.52,  0.38, -0.06,  0.11
  ),
  nrow = 6, ncol=9, byrow=TRUE,
  dimnames = list(
    c('PPP1R3A', 'ENSM6801', 'ENSM8088', 'BMT2', 'TMEM168', 'LSMEM1'),
    c(metadata$sample))
  )
  expect_equal(out[5,3], answer[5,3]) ## 0
  expect_equal(out[6,8], answer[6,8]) ## -0.06
})

## INDIVIDUAL GROUPING

out <- normalizeCountsForHeatmapByIndividual(counts, metadata,
                                             group_var = 'timepoint', baseline='D0', 
                                             individual_var = 'individual')

answer <- matrix(data = c(
  0, 0, 0,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,
  0, 0, 0,  0.05, -0.06,  0.26, -0.33,  0.00,  0.94,
  0, 0, 0,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,
  0, 0, 0,  0.01,  0.15, -0.14, -0.16, -0.27, -0.27,
  0, 0, 0,  0.11,  0.19, -0.16,  0.29, -0.45, -0.53,
  0, 0, 0,  0.26, -0.01,  0.52, -0.36,  0.22,  0.11
),
nrow = 6, ncol=9, byrow=TRUE,
dimnames = list(
  c('PPP1R3A', 'ENSM6801', 'ENSM8088', 'BMT2', 'TMEM168', 'LSMEM1'),
  c(metadata$sample))
)

test_that("Within individual", {
  expect_equal(sum(out[,1]), 0)
  expect_equal(sum(out[,3]), 0)
  expect_equal(out[3,4], answer[3,4]) ## 0
  expect_equal(out[6,8], answer[6,8]) ## 0.22
})
  