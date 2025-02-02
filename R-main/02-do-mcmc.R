

#### Load libraries, source files, load the data  ####
library(ggplot2)
library(lubridate)
library(dplyr)
library(reshape2)
library(stringr)
library(parallel)
library(GGally)

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
viz_results <- mclapply(viz_fisheries, function(x) { 
  s <- squash_data_for_estimation(r = dat$recovery,
                                  cs = dat$catch_sample,
                                  fp = dat$mark_and_tag_rate,
                                  rg = x)
  ret <- cwt_ppn_estimation_function(s, 10000, 10, no_ad_tag_ceiling = 0.05, 
                              unclipped_awts_to_U_minus = TRUE)
  ret
  },
  mc.cores = 8)

if(FALSE) {
#### Then, do the runs for each ETD fishery ####
# just do it for one now.  
}

#### Some functions to assess the mcmc.  I will spin these out to another file at some point ####
extract_mcmc_uas <- function(x, thin = 1) {
   y <- as.data.frame(do.call(rbind, x$theta_uas))
   cbind(iteration = seq(from = thin, by = thin, length.out = nrow(y)), y)
}


mean_mcmc_theta_gs <- function(x, toss_for_burn_in = 100) {
  tmp <- rowMeans(do.call(cbind, x$theta_gs[-(1:toss_for_burn_in)]))
  tbl_df(data.frame(tag_code = names(tmp), post_mean_theta_gs = tmp))
}

# here is one to extract the credible intervals
cred_int_mcmc_theta_gs <- function(x, quants = c(.05, 0.95), toss_for_burn_in = 100) {
  tmp <- do.call(cbind, x$theta_gs[-(1:toss_for_burn_in)]) # get a matrix of results, each iteration in one column
  tmp2 <- apply(tmp, 1, quantile, probs = quants) # get the quantiles
  tbl_df(data.frame(tag_code = colnames(tmp2), CI_lo = tmp2[1,], CI_hi = tmp2[2,]))
}


# this is for visually sampled fisheries
mean_mcmc_expect_pred <- function(x) {
  list(counts = colMeans(do.call(rbind, lapply(x$viz_expect_predict[-1], function(y) sapply(y[1:5], mean)))),
       tag_groups = colMeans(do.call(rbind, lapply(x$viz_expect_predict[-1], function(y) y$tag_group_counts)))
  )
}

# this is for ETD
mean_mcmc_expect_pred_ETD <- function(x) {
  list(counts = colMeans(do.call(rbind, lapply(x$etd_expect_predict[-1], function(y) sapply(y[1:10], mean)))),
       tag_groups = colMeans(do.call(rbind, lapply(x$etd_expect_predict[-1], function(y) y$tag_group_counts)))
  )
}

plot_uas_traces <- function(results, facet_cols = 2) {
  tmp <- plyr::ldply(lapply(results, function(x) {extract_mcmc_uas(x, thin = 1)}), data.frame)
  viz_uas <- melt(tmp, id.vars = c(".id", "iteration"), variable.name = "parameter")
  ggplot(viz_uas, aes_string(x = "iteration", y = "value", colour = "parameter")) + 
    geom_line() +
    facet_wrap(~ .id, ncol = facet_cols)
}

#### Plot some traces with ggplot  ####
plot_uas_traces(viz_results)
ggsave("uas_traces.pdf", width = 14, height = 10)

if(FALSE) {
plot_uas_traces(etd_results)
ggsave("uas_traces_etd.pdf", width = 14, height = 10)
}



#### make scatterplots of observed cwt counts versus the mean predictive ones from the mcmc simulation --VISUAL ####
### This is for visual fisheries
# first get the predicted quantities
tmp <- lapply(lapply(viz_results, mean_mcmc_expect_pred), function(y) 
  data.frame(tag_code = names(y$tag_groups), pred_counts = y$tag_groups, stringsAsFactors = FALSE))


tmp2 <- tbl_df(plyr::ldply(tmp, data.frame))
names(tmp2)[1] <- "recovery_group"

# then get the actual counts
tmp3 <- dat$recovery %>% 
  group_by(recovery_group, tag_code) %>%
  tally()
names(tmp3)[3] <- "obs_counts"

# left join them, and then set things to zero if they are NA in obs_counts
compare_counts <- left_join(tmp2, tmp3)
compare_counts$obs_counts[is.na(compare_counts$obs_counts)] <- 0

# now, also join in the marking and tagging fractions
compare_counts2 <- compare_counts %>% 
  inner_join(dat$mark_and_tag_rate)

# here we color it by the proportion of fish marked
ggplot(compare_counts2, aes(x = obs_counts, y = pred_counts, colour = f_marked)) +
  geom_point() + 
  geom_abline(intercept = 0, slope = 1, colour = "grey50", size = 0.3) +
  scale_colour_gradientn(colours = rev(rainbow(7))) +
  facet_wrap(~ recovery_group, ncol = 2, scales = "free")
ggsave("actual_vs_pred_cwt_recoveries_f_marked.pdf", width = 14, height = 8)

# here is the money shot!  The ones that are off the line are ones that
# apparently have no tags in the ad-clipped segment of the population.  So,
# clearly that is fishy...I'll have to figure out what is going on there.
# Aha! Those are mass-marked stocks whose proportions aren't being estimated correctly,
# which explains why the number of non-tagged ad-clipped fish is underestimated in those
# fisheries, too.
ggplot(compare_counts2, aes(x = obs_counts, y = pred_counts, colour = p_marked)) +
  geom_point() + 
  geom_abline(intercept = 0, slope = 1, colour = "grey50", size = 0.3) +
  scale_colour_gradientn(colours = rev(rainbow(7))) +
  facet_wrap(~ recovery_group, ncol = 2, scales = "free")
ggsave("actual_vs_pred_cwt_recoveries_p_marked.pdf", width = 14, height = 8)



if(FALSE) {
#### Make SCATTERPLOTS FOR THE ETD FISHERIES  ####
tmp <- lapply(lapply(etd_results, mean_mcmc_expect_pred_ETD), function(y) 
  data.frame(tag_code = names(y$tag_groups), pred_counts = y$tag_groups, stringsAsFactors = FALSE))

tmp2 <- tbl_df(plyr::ldply(tmp, data.frame))
names(tmp2)[1] <- "recovery_group"

# then get the actual counts
tmp3 <- dat$recovery %>% 
  group_by(recovery_group, tag_code) %>%
  tally()
names(tmp3)[3] <- "obs_counts"

# left join them, and then set things to zero if they are NA in obs_counts
compare_counts <- left_join(tmp2, tmp3)
compare_counts$obs_counts[is.na(compare_counts$obs_counts)] <- 0

# now, also join in the marking and tagging fractions
compare_counts2 <- compare_counts %>% 
  inner_join(dat$mark_and_tag_rate)

# here we color it by the proportion of fish marked
ggplot(compare_counts2, aes(x = obs_counts, y = pred_counts, colour = f_marked)) +
  geom_point() + 
  geom_abline(intercept = 0, slope = 1, colour = "grey50", size = 0.3) +
  scale_colour_gradientn(colours = rev(rainbow(7))) +
  facet_wrap(~ recovery_group, ncol = 2, scales = "free")
ggsave("etd_actual_vs_pred_cwt_recoveries_f_marked.pdf", width = 14, height = 8)

# The ones that are off the line are ones that
# apparently have few if any tags in the ad-clipped segment of the population. 
ggplot(compare_counts2, aes(x = obs_counts, y = pred_counts, colour = p_marked)) +
  geom_point() + 
  geom_abline(intercept = 0, slope = 1, colour = "grey50", size = 0.3) +
  scale_colour_gradientn(colours = rev(rainbow(7))) +
  facet_wrap(~ recovery_group, ncol = 2, scales = "free")
ggsave("etd_actual_vs_pred_cwt_recoveries_p_marked.pdf", width = 14, height = 8)
}


#### now, have a look at actual and predicted number of ad-clipped fish and ad-clip-cwt fish For VISUAL  ####
tmp <- lapply(lapply(viz_results, mean_mcmc_expect_pred), function(y) {
  mat <- matrix(y$counts, nrow=1, byrow = T)
  colnames(mat) <- names(y$counts)
  as.data.frame(mat)
})

# here are the predicted quantities
pred_counts <- plyr::ldply(tmp, data.frame)
names(pred_counts)[c(1, 4, 6)] <- c("recovery_group", "pred_ad_count", "pred_cwt_count")

# here are the actual numbers of ad_clipped fish with cwts ones
tmp2 <- dat$recovery %>% 
  filter(ad_clipped == "yes", cwt_status == "cwt") %>%
  group_by(recovery_group, cwt_status) %>%
  summarise(actual_ad_clipped_cwt_fish = n())

# here are the total numbers of ad_clipped fish (whether or not they have cwts)
tmp3 <- dat$recovery %>% 
  filter(ad_clipped == "yes") %>%
  group_by(recovery_group) %>%
  summarise(actual_all_ad_clipped_fish = n())

# and here we put them all together into a table
actual_v_pred_table <- pred_counts %>%
  inner_join(tmp2) %>%
  inner_join(tmp3) %>%
  select(recovery_group, 
         viz_samp_size, 
         actual_all_ad_clipped_fish, 
         pred_ad_count, 
         actual_ad_clipped_cwt_fish, 
         pred_cwt_count) 


if(FALSE) {
#### now, have a look at actual and predicted number of ad-clipped fish and ad-clip-cwt fish For ETD. INCOMPLETE  ####
tmp <- lapply(lapply(etd_results, mean_mcmc_expect_pred_ETD), function(y) {
  mat <- matrix(y$counts, nrow=1, byrow = T)
  colnames(mat) <- names(y$counts)
  as.data.frame(mat)
})

# here are the predicted quantities
pred_counts <- plyr::ldply(tmp, data.frame)
names(pred_counts)[c(1, 4, 6)] <- c("recovery_group", "pred_ad_count", "pred_cwt_count")

# here are the actual numbers of ad_clipped fish with cwts ones
tmp2 <- dat$recovery %>% 
  filter(ad_clipped == "yes", cwt_status == "cwt") %>%
  group_by(recovery_group, cwt_status) %>%
  summarise(actual_ad_clipped_cwt_fish = n())

# here are the total numbers of ad_clipped fish (whether or not they have cwts)
tmp3 <- dat$recovery %>% 
  filter(ad_clipped == "yes") %>%
  group_by(recovery_group) %>%
  summarise(actual_all_ad_clipped_fish = n())

# and here we put them all together into a table
actual_v_pred_table <- pred_counts %>%
  inner_join(tmp2) %>%
  inner_join(tmp3) %>%
  select(recovery_group, 
         viz_samp_size, 
         actual_all_ad_clipped_fish, 
         pred_ad_count, 
         actual_ad_clipped_cwt_fish, 
         pred_cwt_count) 
}

#### Another interesting plot to look at.  Estimated thetas versus recoveries!!  VISUAL ####

# this should show us how lopsided the different fisheries are with one or two release groups,
# and also how distorted the tagging rates are relative to proportion in fishery, etc.

# first get the mean theta_gs and add them to compare_counts_2
compare_counts_to_theta <- plyr::ldply(lapply(viz_results, mean_mcmc_theta_gs), data.frame) %>% 
  tbl_df %>%
  rename(recovery_group = .id) %>%
  inner_join(compare_counts2, .)


ggplot(compare_counts_to_theta, aes(x = post_mean_theta_gs, y = obs_counts, colour = f_marked)) + 
  geom_point() + 
  scale_colour_gradientn(colours = rev(rainbow(7))) +
  facet_wrap(~ recovery_group, ncol = 2, scales = "free")
ggsave("post_mean_theta_v_counts_f_marked.pdf", width = 14, height = 8)

ggplot(compare_counts_to_theta, aes(x = post_mean_theta_gs, y = obs_counts, colour = p_marked)) + 
  geom_point() + 
  scale_colour_gradientn(colours = rev(rainbow(7))) +
  facet_wrap(~ recovery_group, ncol = 2, scales = "free")
ggsave("post_mean_theta_v_counts_p_marked.pdf", width = 14, height = 8)


ggplot(compare_counts_to_theta, aes(x = post_mean_theta_gs, y = obs_counts, colour = release_location_state)) + 
  geom_point(size = 1.6) + 
  facet_wrap(~ recovery_group, ncol = 2, scales = "free") +
  scale_colour_discrete(name = "Release state") 
ggsave("post_mean_theta_v_counts_release_state.pdf", width = 14, height = 8)

# then do one coloured by the fraction of all the fish that are ad_clipped and tagged   
ggplot(compare_counts_to_theta, aes(x = post_mean_theta_gs, y = obs_counts, colour = f_marked * p_marked)) + 
  geom_point(size = 1.6) + 
  scale_colour_gradientn(name = "Ppn tagged and marked", colours = rev(rainbow(7))) +
  facet_wrap(~ recovery_group, ncol = 2, scales = "free") +
  ylab("Observed number of CWTs recovered") + 
  xlab("Posterior mean theta_g")
ggsave("post_mean_theta_v_counts_f_marked_times_p_marked.pdf", width = 14, height = 8)


if(FALSE) {
#### Another interesting plot to look at.  Estimated thetas versus recoveries!!  ETD   ####

# this should show us how lopsided the different fisheries are with one or two release groups,
# and also how distorted the tagging rates are relative to proportion in fishery, etc.

# first get the mean theta_gs and add them to compare_counts_2
compare_counts_to_theta <- plyr::ldply(lapply(etd_results, mean_mcmc_theta_gs), data.frame) %>% 
  tbl_df %>%
  rename(recovery_group = .id) %>%
  inner_join(compare_counts2, .)


ggplot(compare_counts_to_theta, aes(x = post_mean_theta_gs, y = obs_counts, colour = f_marked)) + 
  geom_point() + 
  scale_colour_gradientn(colours = rev(rainbow(7))) +
  facet_wrap(~ recovery_group, ncol = 2, scales = "free")
ggsave("etd_post_mean_theta_v_counts_f_marked.pdf", width = 14, height = 8)

ggplot(compare_counts_to_theta, aes(x = post_mean_theta_gs, y = obs_counts, colour = p_marked)) + 
  geom_point() + 
  scale_colour_gradientn(colours = rev(rainbow(7))) +
  facet_wrap(~ recovery_group, ncol = 2, scales = "free")
ggsave("etd_post_mean_theta_v_counts_p_marked.pdf", width = 14, height = 8)


ggplot(compare_counts_to_theta, aes(x = post_mean_theta_gs, y = obs_counts, colour = release_location_state)) + 
  geom_point(size = 1.6) + 
  facet_wrap(~ recovery_group, ncol = 2, scales = "free") +
  scale_colour_discrete(name = "Release state") 
ggsave("etd_post_mean_theta_v_counts_release_state.pdf", width = 14, height = 8)

# then do one coloured by the fraction of all the fish that are ad_clipped and tagged   
ggplot(compare_counts_to_theta, aes(x = post_mean_theta_gs, y = obs_counts, colour = f_marked * p_marked)) + 
  geom_point(size = 1.6) + 
  scale_colour_gradientn(name = "Ppn tagged and marked", colours = rev(rainbow(7))) +
  facet_wrap(~ recovery_group, ncol = 2, scales = "free") +
  ylab("Observed number of CWTs recovered") + 
  xlab("Posterior mean theta_g")
ggsave("etd_post_mean_theta_v_counts_f_marked_times_p_marked.pdf", width = 14, height = 8)
}


#### It might be interesting to add Bayesian credible intervals to the proportion estimates ####
cred_ints <- lapply(viz_results, cred_int_mcmc_theta_gs)


#### I want to make a scatterplot matrix like thing with the estimated fishery proportions
#### compared across fisheries---this will give us some idea of how flexible one might be in fiddling
#### with tagging and marking rates to improve recovery rates on a lot of stocks in a lot of fisheries
#### We will use the GGally library

# first I need to cast the thetas into a multi-column data frame grouped by tag_code and recovery_group 
# AND (this is important) if the count is zero in a recovery group, I will set the 
# estimated theta_g for that release group to NA, so we don't compare stocks that are not seen
# in one of the recovery groups.
tmp <- compare_counts_to_theta
tmp$post_mean_theta_gs[tmp$obs_counts == 0] <- NA
#tmp$post_mean_theta_gs[tmp$post_mean_theta_gs > 0.01] <- NA  # zap the few overly-large ones...

A <- tmp %>% 
  select(tag_code, recovery_group, post_mean_theta_gs) %>%
  mutate(recov_group = str_extract(recovery_group, "^([0-9][0-9])-([A-Z][A-Z])") 
         %>% str_replace("-", "_") %>% paste("f", ., sep = "")
         ) %>%  # make a shorter name for it
  select(-recovery_group) %>%
  dcast(., tag_code ~ recov_group, value.var = "post_mean_theta_gs") %>%
  tbl_df

# now, let's join the release location and the p_marked and f_marked on that 
B <- compare_counts_to_theta %>%
  select(tag_code, release_location_state, f_marked, p_marked) %>%
  inner_join(., A)


# now we can run ggpairs.  Gotta fiddle with the size and colors of dots still.
pm <- ggpairs(B, columns = 5:ncol(B),
              color = "release_location_state",
              alpha = 0.6)
pdf("theta_g_matrix.pdf", width = 28, height = 20)
print(pm)
dev.off()


