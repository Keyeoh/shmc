context('double_beta()')

test_that('double_beta() fails on empty inputs', {
  expect_error(double_beta())
})

test_that('double_beta() works as expected on a simple case', {
  mock_ox = matrix(
    c(
      0.18, 0.22, 0.19, 0.19, 0.49, 0.50, 0.49, 0.49,
      0.19, 0.20, 0.22, 0.20, 0.49, 0.50, 0.48, 0.51
    ),
    nrow = 2,
    ncol = 8,
    byrow = TRUE
  )
  mock_bs = matrix(
    c(
      0.42, 0.41, 0.40, 0.38, 0.54, 0.54, 0.56, 0.56,
      0.39, 0.40, 0.39, 0.41, 0.56, 0.56, 0.56, 0.55
    ),
    nrow = 2,
    ncol = 8,
    byrow = TRUE
  )
  mock_groups = c(1, 1, 1, 1, 2, 2, 2, 2)

  result = double_beta(mock_bs, mock_ox, mock_groups)
})
