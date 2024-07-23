test_that("Look for lactate", {
  matches = approximate_lookup("Lactate", "synonyms_cat", urine_metabolites, "lv", "lcs", limit=3)
  expect_equal(dim(matches)[[1]], 3)
  
  expect_equal(matches$name[[1]], "L-Lactic acid")
  expect_equal(matches$best_name[[1]],"Lactate")
  expect_equal(matches$index[[1]], 97)
  expect_equal(matches$best_score[[1]], 1)
  
  expect_equal(matches$name[[2]], "D-Lactic acid")
  expect_equal(matches$best_name[[2]], "LACTate")
  expect_equal(matches$index[[2]], 392)
  expect_equal(matches$best_score[[2]], 1) 
  
  expect_equal(matches$name[[3]], "Hydroxypropionic acid")
  expect_equal(matches$best_name[[3]], "b-Lactate")
  
  expect_equal(matches$index[[3]], 263)
  expect_equal(abs(matches$best_score[[3]] - 0.7777778) < 1e-5, TRUE)
  
  matches = approximate_lookup("Lactate", "synonyms_cat", urine_metabolites, "lv", "lcs", limit=4)
  expect_equal(dim(matches)[[1]], 4)
})
