test_that("npred works with one new sample", {
  Y = as.numeric(as.factor(Cornejo2025$Tongue$mode1$GenderID))
  Ycnt = Y - mean(Y)
  model = triPLS1(Cornejo2025$Tongue$data, Ycnt, numComponents=1)
  expect_no_error(npred(model, Cornejo2025$Tongue$data[1,,]))
})

test_that("npred works with several new samples", {
  Y = as.numeric(as.factor(Cornejo2025$Tongue$mode1$GenderID))
  Ycnt = Y - mean(Y)
  model = triPLS1(Cornejo2025$Tongue$data, Ycnt, numComponents=1)
  expect_no_error(npred(model, Cornejo2025$Tongue$data[1:5,,]))
})
