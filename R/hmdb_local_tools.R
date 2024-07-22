library(dplyr)
library(stringdist)
library(stringr)
library(jsonlite)

hmdb_nmr_1h <- read.csv("R/spectral1hnmr.csv")
urine_metabolites = read.csv("R/metabolites.csv")

urine_metabolites$synonyms <- lapply(urine_metabolites$synonyms, function(x) {fromJSON(unescape_json(x))})


# Function to remove single quotes within double-quoted strings
unescape_json <- function(string) {
  in_double_quotes <- FALSE
  in_single_quotes <- FALSE
  
  result <- ""
  
  for (i in seq_along(strsplit(string, NULL)[[1]])) {
    char <- substr(string, i, i)
    if (in_double_quotes) {
      if (char == "\"") {
        in_double_quotes <- FALSE
        result <- paste0(result, char)
      } else if (char == "'") {
        next  # Skip adding the single quote if we are inside double quotes
      } else {
        result <- paste0(result, char)
      }
    } else {
      if (in_single_quotes) {
        if (char == "'") {
          in_single_quotes <- FALSE
          result <- paste0(result,  "\"")
        } else if (char == "\"") {
          next  # Skip adding the single quote if we are inside double quotes
        } else {
          result <- paste0(result, char)
        }
      } else {
        if (char == "\"") {
          in_double_quotes <- TRUE
          result <- paste0(result, char)
        } else {
          if (char == "'") {
            in_single_quotes <- TRUE
            result <- paste0(result, "\"")
          } else {
            result <- paste0(result, char)
          }
        } 
      }
    }
  }
  
  return(result)
}

#' Range Query Function
#'
#' This function performs range queries on NMR data and returns the set of matching accession IDs.
#'
#' @param ranges A list of ranges, where each range is a vector of the form \code{c(from, to)}.
#' @param nmrdb A data frame containing the NMR data with columns \code{shift1(ppm)}, \code{accession}, and \code{index}.
#' @return A vector containing the accession IDs that match the specified ranges.
#' @importFrom dplyr filter
#' @examples
#' \dontrun{
#' ranges <- list(c(1.0, 1.5), c(2.0, 2.5))
#' nmrdb <- data.frame(shift1.ppm = c(1.1, 2.1), accession = c("A1", "A2"), index = c(1, 2))
#' range_query_and(ranges, nmrdb)
#' }
#' @export
range_query_and <- function(ranges, nmrdb) {
  result_ids <- NULL
  for (rng in ranges) {
    range_index <- which(nmrdb[["shift1.ppm."]] > rng[1] & nmrdb[["shift1.ppm."]] < rng[2])
    res <- nmrdb[range_index, 'accession']
    if (is.null(result_ids)) {
      result_ids <- res
    } else {
      result_ids <- intersect(result_ids, res)
    }
  }
  return(result_ids)
}

#' Multiplet Query Function
#'
#' This function performs range queries on NMR data and returns a data frame with results ranked by similarity.
#'
#' @param query A list of queries, where each query is a vector of the form \code{c(from, to, type)}.
#' @param nmrdb A data frame containing the NMR data with columns \code{shift1.ppm.}, \code{accession}, and \code{type}.
#' @param metabolites A data frame containing the metabolites data with columns \code{accession} and \code{name}.
#' @return A data frame containing the matching metabolites, with columns for \code{accession}, \code{name}, and \code{similarity}, sorted by similarity in descending order.
#' @importFrom dplyr filter arrange desc
#' @examples
#' \dontrun{
#' query <- list(c(1.0, 1.5, "type1"), c(2.0, 2.5, "type2"))
#' nmrdb <- data.frame(shift1.ppm. = c(1.1, 2.1), accession = c("A1", "A2"), type = c("type1", "type2"))
#' metabolites <- data.frame(accession = c("A1", "A2"), name = c("Met1", "Met2"))
#' multiple_query(query, nmrdb, metabolites)
#' }
#' @export
multiplet_query <- function(query, nmrdb, metabolites) {
  res <- range_query_and(query, nmrdb)
  filtered_df <- metabolites %>% filter(accession %in% res)
  result <- data.frame(accession = character(), name = character(), similarity = numeric())
  
  for (i in 1:nrow(filtered_df)) {
    row <- filtered_df[i, ]
    signals <- nmrdb %>% filter(accession == row$accession)
    matches <- 0
    for (query_i in query) {
      from_i <- query_i[1]
      to_i <- query_i[2]
      mul_i <- query_i[3]
      for (j in 1:nrow(signals)) {
        signal <- signals[j, ]
        if (signal[['shift1.ppm.']] >= from_i && signal[['shift1.ppm.']] <= to_i && signal[['type']] == mul_i) {
          matches <- matches + 1
          break
        }
      }
    }
    similarity <- length(query) / nrow(signals) + matches / length(query)
    result <- rbind(result, data.frame(accession = row$accession, name = row$name, similarity = similarity))
  }
  
  result <- result %>% arrange(desc(similarity))
  return(result)
}


#' Approximate Lookup Function
#'
#' This function performs an approximate string matching lookup on a specified column of a data frame.
#'
#' @param df A data frame containing the data to search.
#' @param column A string specifying the name of the column to search.
#' @param search_string A string to search for in the specified column.
#' @param scorer1 A string specifying the method to use for the initial scoring (default is "lv" for Levenshtein distance).
#' @param scorer2 A string specifying the method to use for the secondary scoring (default is "lcs" for Longest common substring distance).
#' @param limit An integer specifying the maximum number of results to return (default is 10).
#' @return A data frame containing the best matches, with columns for the original name, the best matching synonym, the index of the match, and the best score.
#' @importFrom stringdist stringdist
#' @importFrom dplyr filter arrange desc
#' @examples
#' \dontrun{
#' df <- data.frame(name = c("Name1", "Name2"), 
#'                  shift1.ppm. = c(1.0, 2.0), 
#'                  accession = c("A1", "A2"), 
#'                  synonyms = c("Syn1, Syn1Alt", "Syn2, Syn2Alt"))
#' approximate_lookup(df, "name", "name1")
#' }
#' @export
approximate_lookup <- function(df, column, search_string, scorer1 = "lv", scorer2 = "lcs", limit = 10) {
  # Extract the column values as a list
  choices <- tolower(df[[column]])
  names <- df[['name']]
  search_string <- tolower(search_string)
  
  # Get the best matches using stringdist package
  distances <- stringdist::stringdist(search_string, choices, method = scorer1)
  scores <- 1 - distances / max(distances)
  #print(scores)
  matches <- data.frame(index = seq_along(choices), score = scores)
  matches <- matches %>% filter(score > 0.33)
  
  if (nrow(matches) > 0) {
    result <- data.frame(name = character(), best_name = character(), index = numeric(), best_score = numeric())
    
    for (i in 1:nrow(matches)) {
      match_index <- matches$index[i]
      matchings_names <- unlist(df[match_index, 'synonyms'])
      
      best_name <- ""
      best_score <- 0
      
      for (name in c(names[match_index], matchings_names)) {
        #print(name)
        score <- 1 - stringdist::stringdist(search_string, tolower(name), method = scorer2)/ max(nchar(search_string), nchar(name))
        if (score > best_score) {
          best_score <- score
          best_name <- name
        }
      }
      
      result <- rbind(result, data.frame(name = names[match_index], best_name = best_name, index = match_index, best_score = best_score))
    }
    
    result <- result %>% arrange(desc(best_score))
    return(head(result, limit))
  }
  return(NULL)
}
