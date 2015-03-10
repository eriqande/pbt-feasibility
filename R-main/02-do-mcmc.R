

#### Load libraries, source files, load the data  ####
library(ggplot2)
library(lubridate)
library(plyr)
library(dplyr)
library(reshape2)
library(stringr)
library(parallel)

source("R/rmis_cleaning.R")
source("R/mcmc_funcs.R")


dat <- readRDS("data_for_mcmc.rds")



#### Gather up the names of the recovery groups we want to do now ####

# here are recovery groups that have at least some (and possibly all) electronic detection
electro <- dat$catch_sample %>%
  group_by(recovery_group, detection_method) %>%
  summarise_each(funs(sum(., na.rm = TRUE)), vars = mr_1st_partition_size:mr_2nd_sample_obs_adclips) %>%
  filter(detection_method == "E") %>%
  ungroup %>%
  select(recovery_group)

# then get a vector of the visual-detection-only fisheries.  
tmp <- unique(dat$catch_sample$recovery_group) 
viz_fisheries <- tmp[ !(tmp %in% electro$recovery_group)]
names(viz_fisheries) <- viz_fisheries


#### Then, do the runs for each visual fishery ####


# do each of them for a short run
results <- mclapply(viz_fisheries, function(x) { 
  s <- squash_data_for_estimation(r = dat$recovery,
                                  cs = dat$catch_sample,
                                  fp = dat$mark_and_tag_rate,
                                  rg = x)
  cwt_ppn_estimation_function(s, 10000, 10, no_ad_tag_ceiling = -0.05, 
                              unclipped_awts_to_U_minus = TRUE)
  },
  mc.cores = 8)



#### Some functions to assess the mcmc.  I will spin these out to another file at some point ####
extract_mcmc_uas <- function(x, thin = 1) {
   y <- as.data.frame(do.call(rbind, x$theta_uas))
   cbind(iteration = seq(from = thin, by = thin, length.out = nrow(y)), y)
}



#### Plot some traces with ggplot
tmp <-ldply(lapply(results, function(x) {extract_mcmc_uas(x, thin = 1)}), data.frame)
viz_uas <- melt(tmp, id.vars = c(".id", "iteration"), variable.name = "parameter")
ggplot(viz_uas, aes(x = iteration, y = value, colour = parameter)) + 
  geom_line() +
  facet_wrap(~ .id, ncol = 2)
