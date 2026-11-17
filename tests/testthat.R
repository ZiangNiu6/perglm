# Activate the standard testthat dispatcher once test files are added.
if (length(list.files("testthat", pattern = "^test.*\\.[rR]$")) > 0L) {
  library(testthat)
  library(perglm)

  test_check("perglm")
}
