test_that("Look for lactate", {
  query <- list(list(1.25, 1.35, 'd'), list(4.05, 4.15, 'q'))
  result <- multiplet_query(query, hmdb_nmr_1h, urine_metabolites)
  expected <- list("accession"=c("HMDB0000190","HMDB0000030","HMDB0000701","HMDB0002802","HMDB0001348","HMDB0014352"),
    "name"=c("L-Lactic acid", "Biotin", "Hexanoylglycine", "Cortisone", "SM(d18:1/18:0)", "Azithromycin"),
    "similarity"=c(2.0000000, 0.7000000, 0.4000000, 0.1428571, 0.1333333, 0.0625000))
  expect_equal(result$accession, expected$accession)
  expect_equal(result$name, expected$name)
  expect_equal(abs(result$similarity[[1]] - expected$similarity[[1]]) < 1e-3, TRUE)
  expect_equal(abs(result$similarity[[2]] - expected$similarity[[2]]) < 1e-3, TRUE)
  expect_equal(abs(result$similarity[[3]] - expected$similarity[[3]]) < 1e-3, TRUE)
  
  result <- multiplet_query(query)
  
  expect_equal(result$accession, expected$accession)
  expect_equal(result$name, expected$name)
  expect_equal(abs(result$similarity[[1]] - expected$similarity[[1]]) < 1e-3, TRUE)
  expect_equal(abs(result$similarity[[2]] - expected$similarity[[2]]) < 1e-3, TRUE)
  expect_equal(abs(result$similarity[[3]] - expected$similarity[[3]]) < 1e-3, TRUE)
})
