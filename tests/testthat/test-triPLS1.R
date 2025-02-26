test_that("triPLS1 throws no errors for one-component models", {
  Y = as.numeric(as.factor(Cornejo2025$Tongue$mode1$GenderID))
  Ycnt = Y - mean(Y)
  expect_no_error(triPLS1(Cornejo2025$Tongue$data, Ycnt, numComponents=1))
})

test_that("triPLS1 throws no errors for more than one component", {
  Y = as.numeric(as.factor(Cornejo2025$Tongue$mode1$GenderID))
  Ycnt = Y - mean(Y)
  expect_no_error(triPLS1(Cornejo2025$Tongue$data, Ycnt, numComponents=2))
})
