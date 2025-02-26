test_that("ncrossreg works with normal input", {
  Y = as.numeric(as.factor(Cornejo2025$Tongue$mode1$GenderID))
  Ycnt = Y - mean(Y)
  expect_no_error(ncrossreg(Cornejo2025$Tongue$data, Ycnt, cvFolds=2))
})

test_that("jack-knifing throws no errors or warnings", {
  Y = as.numeric(as.factor(Cornejo2025$Tongue$mode1$GenderID))
  Ycnt = Y - mean(Y)
  expect_no_error(ncrossreg(Cornejo2025$Tongue$data, Ycnt, maxNumComponents = 2))
})
