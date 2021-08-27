# Test EI library
test_that("read_lib() reads EI library correctly", {
  EI_file <- system.file("EI.msp", package = "mspcompiler")

  expect_identical(read_lib(EI_file, type = "EI"), EI)
})

# Test MS2 msp library
test_that("read_lib() reads MS2 msp library correctly", {
  MS2_msp_file <- system.file("MS2.msp", package = "mspcompiler")

  expect_identical(read_lib(MS2_msp_file), MS2_msp)
})

# Test MS2 mgf library
test_that("read_lib() reads MS2 mgf library correctly", {
  MS2_mgf_file <- system.file("MS2.mgf", package = "mspcompiler")

  expect_identical(read_lib(MS2_mgf_file, format = "mgf"), MS2_mgf)
})
