#####################################################################
##########                  README                         ##########
#####################################################################

# Iniatially this project was performed using relative abundances. This script 
# illustrates by simulation that clr transformation is required to analyze 
# whether scfa counts in the sample are associated with bmi by comparing 
# different approaches in their performance. Also, that counts must be added
# BEFORE clr transformation. Alternatively, other log-ratio 
# approaches could have been used such as summed-log-ratios or balances.

options(tidyverse.quiet = TRUE)
options(dplyr.summarise.inform = FALSE)
library(tidyverse)
library(MCMCpack)
library(compositions)


#####################################################################
##########                   clr vs count                  ##########
#####################################################################

iterations <- as.list(1:1e3)

result <- map(iterations, function(i) {
  # 2 microbiota samples that we compare (a and b)
  n <- 2
  # equal readcounts (change to unequal in potential next step)
  readcounts_a <- 1e4
  readcounts_b <- 1e4
  # simulate probability for each genus draw (20 genera) for two samples
  a <- rep(1, 20)
  probs <- rdirichlet(2, a)

  # suppose that the columns c1 - c2 are SCFA producers
  c1 <- 11
  c2 <- 20
  # i store the ground truth just to compare if maybe clr sometimes is closer 
  # to that or not 
  scfa_prob_a <- sum(probs[1, c1:c2])
  scfa_prob_b <- sum(probs[2, c1:c2])
  ground_truth <- ifelse(scfa_prob_a > scfa_prob_b, "a", "b")

  # generate 2 samples 
  sample_a <- as.numeric(rmultinom(n = 1, size = readcounts_a, prob = probs[1, ]))
  sample_b <- as.numeric(rmultinom(n = 1, size = readcounts_b, prob = probs[2, ]))

  # indicate in tibble which are the SCFA producers
  df <- tibble(
    sample = rep(c("a", "b"), each = 20),
    counts = c(sample_a, sample_b),
    scfa = rep(c(rep(0, 10), rep(1, 10)), 2)) 
    

  # who has more scfa producers according to absolute counts?
  count_sum <- df %>% 
    filter(scfa == 1) %>% 
    group_by(sample) %>% 
    summarise(scfa_sum = sum(counts))

  bigger <- ifelse(as.numeric(count_sum[1, 2]) > as.numeric(count_sum[2, 2]), "a", "b")


  # now perform clr transformation and see if we get the same conclusion:
  sample_a_clr <- as.numeric(clr(sample_a))
  sample_b_clr <- as.numeric(clr(sample_b))

  df_clr <- tibble(
    sample = rep(c("a", "b"), each = 20),
    clr = c(sample_a_clr, sample_b_clr),
    scfa = rep(c(rep(0, 10), rep(1, 10)), 2))


  # who has more scfa producers according to clr sum?
  clr_sum <- df_clr %>% filter(scfa == 1) %>% 
    group_by(sample) %>% 
    summarise(scfa_sum_clr = sum(clr)) 

  bigger_clr <- ifelse(as.numeric(clr_sum[1, 2]) > as.numeric(clr_sum[2, 2]), "a", "b")

  data <- bind_cols(df, dplyr::select(df_clr, "clr")) %>%
          dplyr::select(sample, scfa, counts, clr)

  summmed_data <- bind_cols(count_sum, clr_sum[, 2])
  
  return(list(
    "equal_conclusion" = bigger == bigger_clr,
    "clr_correct" = bigger_clr == ground_truth,
    "count_correct" = bigger == ground_truth,
    "data" = data,
    "summed" = summmed_data
  ))
})  



# % agreement between both methods should be 100% but is:
map_dbl(result, ~.x$equal_conclusion) %>% mean()
# double check which method is closer to ground truth:
map_dbl(result, ~.x$clr_correct) %>% mean()
map_dbl(result, ~.x$count_correct) %>% mean()

# examples of disagreement:
ind <- c()
for (i in 1:1e3) {
  if (!result[[i]]$equal_conclusion) {
    ind <- c(ind, i)
  }
}
ind

examp <- ind[7]
  result[[examp]]




#####################################################################
##########                   rel vs count                  ##########
#####################################################################

iterations <- as.list(1:1e3)

result2 <- map(iterations, function(i) {
  # 2 microbiota samples that we compare (a and b)
  n <- 2
  # equal readcounts (change to unequal in potential next step)
  readcounts_a <- 1e4
  readcounts_b <- 1e4
  # simulate probability for each genus draw (20 genera) for two samples
  a <- rep(1, 20)
  probs <- rdirichlet(2, a)
  
  # i store the ground truth just to compare if maybe clr sometimes is closer 
  # to that or not 
  scfa_prob_a <- sum(probs[1, 11:20])
  scfa_prob_b <- sum(probs[2, 11:20])
  ground_truth <- ifelse(scfa_prob_a > scfa_prob_b, "a", "b")
  
  sample_a <- as.numeric(rmultinom(n = 1, size = readcounts_a, prob = probs[1, ]))
  sample_b <- as.numeric(rmultinom(n = 1, size = readcounts_b, prob = probs[2, ]))

  # lets say that half of these randomly generated genera are scfa producers
  df <- tibble(
    sample = rep(c("a", "b"), each = 20),
    counts = c(sample_a, sample_b),
    scfa = rep(c(rep(0, 10), rep(1, 10)), 2)) 
    

  # who has more scfa producers according to absolute counts?
  count_sum <- df %>% filter(scfa == 1) %>% 
    group_by(sample) %>% 
    summarise(scfa_sum = sum(counts))

  bigger <- ifelse(as.numeric(count_sum[1, 2]) > as.numeric(count_sum[2, 2]), "a", "b")


  # now perform clr transformation and see if we get the same conclusion:
  sample_a_rel <- as.numeric(sample_a/sum(as.numeric(sample_a)))
  sample_b_rel <- as.numeric(sample_b/sum(as.numeric(sample_b)))

  df_rel <- tibble(
    sample = rep(c("a", "b"), each = 20),
    rel = c(sample_a_rel, sample_b_rel),
    scfa = rep(c(rep(0, 10), rep(1, 10)), 2))


  # who has more scfa producers according to clr sum?
  rel_sum <- df_rel %>% filter(scfa == 1) %>% 
    group_by(sample) %>% 
    summarise(scfa_sum_rel = sum(rel)) 

  bigger_rel <- ifelse(as.numeric(rel_sum[1, 2]) > as.numeric(rel_sum[2, 2]), "a", "b")
  
  data <- bind_cols(df, dplyr::select(df_rel, "rel")) %>%
          dplyr::select(sample, scfa, counts, rel)
  
  summmed_data <- bind_cols(count_sum, rel_sum[, 2])
  
  return(list(
    "equal_conclusion" = bigger == bigger_rel,
    "rel_correct" = bigger_rel == ground_truth,
    "count_correct" = bigger == ground_truth,
    "data" = data,
    "summed" = summmed_data
  ))
  
})  



# % agreement between both methods should be 100% but is:
map_dbl(result2, ~.x$equal_conclusion) %>% mean()
# double check which method is closer to ground truth:
map_dbl(result2, ~.x$rel_correct) %>% mean()
map_dbl(result2, ~.x$count_correct) %>% mean()

# examples of disagreement:
ind <- c()
for (i in 1:1e3) {
  if (!result2[[i]]$equal_conclusion) {
    ind <- c(ind, i)
  }
}
ind

examp <- ind[2]
  result2[[examp]]
  



#####################################################################
##########              now with rarefied counts           ##########
#####################################################################

iterations <- as.list(1:1e3)

result <- map(iterations, function(i) {
  # 2 microbiota samples that we compare (a and b)
  n <- 2
  # equal readcounts (change to unequal in potential next step)
  readcounts_a <- ceiling(runif(1, 5000, 10000))
  readcounts_b <- ceiling(runif(1, 5000, 10000)) 

  # 20 genus names 
  genera <- glue::glue("genus{11:30}")
  # simulate probability for each genus draw (20 genera) for two samples
  a <- rep(1, 20)
  probs <- rdirichlet(2, a)

  # i store the ground truth just to compare if maybe clr sometimes is closer 
  # to that or not 
  scfa_prob_a <- sum(probs[1, 11:20])
  scfa_prob_b <- sum(probs[2, 11:20])
  ground_truth <- ifelse(scfa_prob_a > scfa_prob_b, "a", "b")

  sample_a <- as.numeric(rmultinom(n = 1, size = readcounts_a, prob = probs[1, ]))
  sample_b <- as.numeric(rmultinom(n = 1, size = readcounts_b, prob = probs[2, ]))

  genera_a <- map2(sample_a, genera, function(count, genus) {
    rep(genus, count)}) %>% unlist() 
  genera_b <- map2(sample_b, genera, function(count, genus) {
    rep(genus, count)}) %>% unlist()



  # perform rarefaction:
  rar_count <- min(sum(sample_a), sum(sample_b))
  rar_a <- sample(genera_a, rar_count, replace = T)
  rar_b <- sample(genera_b, rar_count, replace = T)
  sample_a_rar <- rar_a %>% sort() %>% plyr::count()
  sample_b_rar <- rar_b %>% sort() %>% plyr::count()

  sample_a_rar <- tibble(
    genus = genera, 
    freq = ifelse(genera %in% sample_a_rar$x, sample_a_rar$freq, 0)) %>%
    .$freq 
  sample_b_rar <- tibble(
    genus = genera, 
    freq = ifelse(genera %in% sample_b_rar$x, sample_b_rar$freq, 0)) %>%
    .$freq 


  # lets say that half of these randomly generated genera are scfa producers
  df <- tibble(
    sample = rep(c("a", "b"), each = 20),
    counts = c(sample_a_rar, sample_b_rar),
    scfa = rep(c(rep(0, 10), rep(1, 10)), 2)) 


  # who has more scfa producers according to absolute counts?
  count_sum <- df %>% filter(scfa == 1) %>% 
    group_by(sample) %>% 
    summarise(scfa_sum = sum(counts))

  bigger <- ifelse(as.numeric(count_sum[1, 2]) > as.numeric(count_sum[2, 2]), "a", "b")


  # now perform clr transformation and see if we get the same conclusion:
  sample_a_clr <- as.numeric(clr(sample_a))
  sample_b_clr <- as.numeric(clr(sample_b))

  df_clr <- tibble(
    sample = rep(c("a", "b"), each = 20),
    clr = c(sample_a_clr, sample_b_clr),
    scfa = rep(c(rep(0, 10), rep(1, 10)), 2))



  # who has more scfa producers according to clr sum?
  clr_sum <- df_clr %>% filter(scfa == 1) %>% 
    group_by(sample) %>% 
    summarise(scfa_sum_clr = sum(clr)) 

  bigger_clr <- ifelse(as.numeric(clr_sum[1, 2]) > as.numeric(clr_sum[2, 2]), "a", "b")

  data <- bind_cols(df, dplyr::select(df_clr, "clr")) %>%
          dplyr::select(sample, scfa, counts, clr)

  summmed_data <- bind_cols(count_sum, clr_sum[, 2])

  return(list(
    "equal_conclusion" = bigger == bigger_clr,
    "clr_correct" = bigger_clr == ground_truth,
    "rarefied_correct" = bigger == ground_truth,
    "data" = data,
    "summed" = summmed_data
  ))
})  





# % agreement between both methods should be 100% but is:
map_dbl(result, ~.x$equal_conclusion) %>% mean()
# double check which method is closer to ground truth:
map_dbl(result, ~.x$clr_correct) %>% mean()
map_dbl(result, ~.x$rarefied_correct) %>% mean()



# examples of disagreement:
ind <- c()
for (i in 1:1e3) {
  if (!result[[i]]$equal_conclusion) {
    ind <- c(ind, i)
  }
}
ind

examp <- ind[2]
  result[[examp]]






#####################################################################
##########              rarefied counts vs rel ab          ##########
#####################################################################

iterations <- as.list(1:1e3)

result <- map(iterations, function(i) {
  # 2 microbiota samples that we compare (a and b)
  n <- 2
  # equal readcounts (change to unequal in potential next step)
  readcounts_a <- ceiling(runif(1, 5000, 10000))
  readcounts_b <- ceiling(runif(1, 5000, 10000)) 

  # 20 genus names 
  genera <- glue::glue("genus{11:30}")
  # simulate probability for each genus draw (20 genera) for two samples
  a <- rep(1, 20)
  probs <- rdirichlet(2, a)

  # i store the ground truth just to compare if maybe clr sometimes is closer 
  # to that or not 
  scfa_prob_a <- sum(probs[1, 11:20])
  scfa_prob_b <- sum(probs[2, 11:20])
  ground_truth <- ifelse(scfa_prob_a > scfa_prob_b, "a", "b")

  sample_a <- as.numeric(rmultinom(n = 1, size = readcounts_a, prob = probs[1, ]))
  sample_b <- as.numeric(rmultinom(n = 1, size = readcounts_b, prob = probs[2, ]))

  genera_a <- map2(sample_a, genera, function(count, genus) {
    rep(genus, count)}) %>% unlist() 
  genera_b <- map2(sample_b, genera, function(count, genus) {
    rep(genus, count)}) %>% unlist()



  # perform rarefaction:
  rar_count <- min(sum(sample_a), sum(sample_b))
  rar_a <- sample(genera_a, rar_count, replace = T)
  rar_b <- sample(genera_b, rar_count, replace = T)
  sample_a_rar <- rar_a %>% sort() %>% plyr::count()
  sample_b_rar <- rar_b %>% sort() %>% plyr::count()

  sample_a_rar <- tibble(
    genus = genera, 
    freq = ifelse(genera %in% sample_a_rar$x, sample_a_rar$freq, 0)) %>%
    .$freq 
  sample_b_rar <- tibble(
    genus = genera, 
    freq = ifelse(genera %in% sample_b_rar$x, sample_b_rar$freq, 0)) %>%
    .$freq 


  # lets say that half of these randomly generated genera are scfa producers
  df <- tibble(
    sample = rep(c("a", "b"), each = 20),
    counts = c(sample_a_rar, sample_b_rar),
    scfa = rep(c(rep(0, 10), rep(1, 10)), 2)) 


  # who has more scfa producers according to absolute counts?
  count_sum <- df %>% filter(scfa == 1) %>% 
    group_by(sample) %>% 
    summarise(scfa_sum = sum(counts))

  bigger <- ifelse(as.numeric(count_sum[1, 2]) > as.numeric(count_sum[2, 2]), "a", "b")


  # now perform clr transformation and see if we get the same conclusion:
  sample_a_rel <- as.numeric(sample_a/sum(as.numeric(sample_a)))
  sample_b_rel <- as.numeric(sample_b/sum(as.numeric(sample_b)))

  df_rel <- tibble(
    sample = rep(c("a", "b"), each = 20),
    rel = c(sample_a_rel, sample_b_rel),
    scfa = rep(c(rep(0, 10), rep(1, 10)), 2))


  # who has more scfa producers according to clr sum?
  rel_sum <- df_rel %>% filter(scfa == 1) %>% 
    group_by(sample) %>% 
    summarise(scfa_sum_rel = sum(rel)) 

  bigger_rel <- ifelse(as.numeric(rel_sum[1, 2]) > as.numeric(rel_sum[2, 2]), "a", "b")

  data <- bind_cols(df, dplyr::select(df_rel, "rel")) %>%
          dplyr::select(sample, scfa, counts, rel)

  summmed_data <- bind_cols(count_sum, rel_sum[, 2])

  return(list(
    "equal_conclusion" = bigger == bigger_rel,
    "rel_correct" = bigger_rel == ground_truth,
    "rar_correct" = bigger == ground_truth,
    "data" = data,
    "summed" = summmed_data
  ))
})  





# % agreement between both methods should be 100% but is:
map_dbl(result, ~.x$equal_conclusion) %>% mean()
# double check which method is closer to ground truth:
map_dbl(result, ~.x$rel_correct) %>% mean()
map_dbl(result, ~.x$count_correct) %>% mean()

# examples of disagreement:
ind <- c()
for (i in 1:1e3) {
  if (!result[[i]]$equal_conclusion) {
    ind <- c(ind, i)
  }
}
ind

examp <- ind[2]
  result[[examp]]











#####################################################################
##########              count vs FIRST adding and then clr ##########
#####################################################################

iterations <- as.list(1:1e3)

result <- map(iterations, function(i) {
  # 2 microbiota samples that we compare (a and b)
  n <- 2
  # equal readcounts (change to unequal in potential next step)
  readcounts_a <- 1e4
  readcounts_b <- 1e4
  # simulate probability for each genus draw (30 genera) for two samples
  a <- rep(1, 30)
  probs <- rdirichlet(2, a)

  # suppose that the columns c1 - c2 are SCFA producers
  c1 <- 21
  c2 <- 30
  # i store the ground truth just to compare if maybe clr sometimes is closer 
  # to that or not 
  scfa_prob_a <- sum(probs[1, c1:c2])
  scfa_prob_b <- sum(probs[2, c1:c2])
  ground_truth <- ifelse(scfa_prob_a > scfa_prob_b, "a", "b")



  # generate 2 samples 
  sample_a <- as.numeric(rmultinom(n = 1, size = readcounts_a, prob = probs[1, ]))
  sample_b <- as.numeric(rmultinom(n = 1, size = readcounts_b, prob = probs[2, ]))

  # indicate in tibble which are the SCFA producers
  df <- tibble(
    sample = rep(c("a", "b"), each = 30),
    counts = c(sample_a, sample_b),
    scfa = rep(c(rep(0, 20), rep(1, 10)), 2)) 

  # who has more scfa producers according to absolute counts?
  count_sum <- df %>% 
    filter(scfa == 1) %>% 
    group_by(sample) %>% 
    summarise(scfa_sum = sum(counts))

  bigger <- ifelse(as.numeric(count_sum[1, 2]) > as.numeric(count_sum[2, 2]), "a", "b")


  # now perform clr transformation and see if we get the same conclusion:
  # first replace the SCFA counts by one SCFA count that is the sum:
  sample_a_manip <- c(sample_a[1:20], sum(sample_a[21:30]))
  sample_b_manip <- c(sample_b[1:20], sum(sample_b[21:30]))
  sample_a_manip
  sample_b_manip
  df
  count_sum
  sample_a_clr <- as.numeric(clr(sample_a_manip))
  sample_b_clr <- as.numeric(clr(sample_b_manip))

  sample_a_clr
  sample_b_clr


  bigger_clr <- ifelse(sample_a_manip[21] > sample_b_manip[21], "a", "b")
  
  return(list(
    "equal_conclusion" = bigger == bigger_clr,
    "clr_correct" = bigger_clr == ground_truth,
    "count_correct" = bigger == ground_truth,
    "count_sum" = count_sum,
    "a" = sample_a_clr,
    "b" = sample_b_clr
  ))
    
})  



# % agreement between both methods should be 100% but is:
map_dbl(result, ~.x$equal_conclusion) %>% mean()
# double check which method is closer to ground truth:
map_dbl(result, ~.x$clr_correct) %>% mean()
map_dbl(result, ~.x$count_correct) %>% mean()

# examples of disagreement:
ind <- c()
for (i in 1:1e3) {
  if (!result[[i]]$equal_conclusion) {
    ind <- c(ind, i)
  }
}
ind

examp <- ind[7]
  result[[examp]]







# conclusion:
# we should use clr if we sum up scfa before









