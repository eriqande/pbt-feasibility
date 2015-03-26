

# you gotta run 01 and 03 before this.


#### 0.  Drop the BC and OR recovery areas from our results
lean_results <- viz_results[!(names(viz_results) %in% c("13-OR  (6524)", "07-BC  (3509)"))]

#### 1. get the posterior mean theta_gs and theta_uas  ####

# we get them and store them is a list of lists.  Each component is a list with two
# components: uas and theta_gs, which are both named vectors.
post_means <- lapply(lean_results, function(x) {
  tmp <- mean_mcmc_theta_gs(x)
  theta_gs <- tmp$post_mean_theta_gs
  names(theta_gs) <- tmp$tag_code
  list(
    uas = colMeans(extract_mcmc_uas(x)[-(1:100),])[-1],
    theta_gs = theta_gs
  )
})


#### Write a function to deliver the predicted number of recoveries ####
# we will key this to a sample of 50,000 total fish that are investigated for ad-clips
# x is a component of post_means
# m_and_t is the mark_and_tag info to use
clean_predict <- function(n = 50000, x, scenario = "Current_CWT", m_and_t = dat$mark_and_tag_rate) {
  y <- viz_expect_predict(n, x$theta_gs, x$uas, m_and_t)
  tmp <- sapply(y[1:5], "[", 1)
  gross_counts <- tbl_df(data.frame(what_count = names(tmp), value = tmp, stringsAsFactors = FALSE))
  theta_g_counts <- tbl_df(data.frame(tag_code = names(y$tag_group_counts), value = y$tag_group_counts, stringsAsFactors = FALSE))
  names(gross_counts)[2] <- scenario
  names(theta_g_counts)[2] <- scenario
  list(gross_counts = gross_counts, theta_g_counts = theta_g_counts)
}


#### 2. Now, collect the expected number of recoveries given that under the current CWT tagging rates  ####
current_cwt <- lapply(post_means, function(x) clean_predict(x = x))
current_cwt_recovs <- tbl_df(plyr::ldply(current_cwt, function(x) as.data.frame(x$theta_g_counts), .id = "recovery_area"))

#### 3. Now collect expected number of recoveries when mark rate is the same, but the tag-rate is 100% for marked fish 

hundy_tag <- dat$mark_and_tag_rate
hundy_tag$p_marked <- 1.000
hundy_scenario <- lapply(post_means, function(x) clean_predict(x = x, m_and_t = hundy_tag, scenario = "Tag_100_Same_Mark_Rate"))
hundy_scenario_recovs <- tbl_df(plyr::ldply(hundy_scenario, function(x) as.data.frame(x$theta_g_counts), .id = "recovery_area"))


#### 4. Now plot those two for fun

comp1 <- inner_join(current_cwt_recovs, hundy_scenario_recovs, by = c("recovery_area", "tag_code")) %>%
  inner_join(dat$mark_and_tag_rate, by = "tag_code")

ggplot(comp1, aes(x = Current_CWT, y = (Tag_100_Same_Mark_Rate - Current_CWT) / Current_CWT, colour  = p_marked)) +
  geom_point() +
  facet_wrap(~ recovery_area, ncol = 2, scales = "free") +
  scale_colour_gradientn(colours = rev(rainbow(7))) 


#### 5. Now, consider the 100% mark-and-tag scenario, just to see how bad it could be ####
## We should probably exempt the DIT groups

mark_tag_all <- dat$mark_and_tag_rate
mark_tag_all <- inner_join(mark_tag_all, releases %>% select(tag_code_or_release_id, related_group_type),
           by = c("tag_code" = "tag_code_or_release_id"))

notDIT <- mark_tag_all$related_group_type != "D"


mark_tag_all$p_marked[notDIT] <- 1.000
mark_tag_all$f_marked[notDIT] <- 1.000
mark_tag_all$f_unmarked[notDIT] <- 0.000

mark_tag_all_scenario <- lapply(post_means, function(x) clean_predict(x = x, m_and_t = mark_tag_all, scenario = "Mark_and_Tag_All"))
mark_tag_all_recovs <- tbl_df(plyr::ldply(mark_tag_all_scenario, function(x) as.data.frame(x$theta_g_counts), .id = "recovery_area"))


#### 6. Now, go ahead and plot that assuming 50K fish still inspected for ad-clips ####

comp2 <- inner_join(current_cwt_recovs, mark_tag_all_recovs, by = c("recovery_area", "tag_code")) %>%
  inner_join(dat$mark_and_tag_rate, by = "tag_code")

ggplot(comp2, aes(x = Current_CWT, y = (Mark_and_Tag_All - Current_CWT) / Current_CWT, colour  = p_marked)) +
  geom_point() +
  facet_wrap(~ recovery_area, ncol = 2, scales = "free") +
  scale_colour_gradientn(colours = rev(rainbow(7)))


#### 7. That is all fine and well, but what about if we hold the number of recovered tags constant? ####
# i.e. how about if the sampling rate is adjusted so that we "decode" the same number of PBT tags?
mark_and_tag_all_gross <- tbl_df(plyr::ldply(mark_tag_all_scenario, function(x) as.data.frame(x$gross_counts), .id = "recovery_area"))
cwt_current_gross <- tbl_df(plyr::ldply(current_cwt, function(x) as.data.frame(x$gross_counts), .id = "recovery_area"))


# here is the fraction smaller that our sample should be to get the same number of ad-clipped fish recovered 
downsample_fracts <- inner_join(cwt_current_gross, mark_and_tag_all_gross, by = c("recovery_area", "what_count")) %>%
  filter(what_count == "ad_count") %>%
  mutate(downsample_fract = Current_CWT / Mark_and_Tag_All) %>%
  select(recovery_area, downsample_fract)

# now inner_join that so we can use it:
comp3 <- inner_join(comp2, downsample_fracts, by = "recovery_area") %>%
  mutate(Mark_and_Tag_All_ad_clip_count_const = Mark_and_Tag_All * downsample_fract)

ggplot(comp3, aes(x = Current_CWT, y = (Mark_and_Tag_All_ad_clip_count_const - Current_CWT) / Current_CWT, colour  = p_marked)) +
  geom_point() +
  facet_wrap(~ recovery_area, ncol = 2, scales = "free") +
  scale_colour_gradientn(colours = rev(rainbow(7)))

# That is really interesting.  Let's drill down on these and see
# what we really have.

#### Let's take all the different fisheries and cut on the expected number of recoveries under CWT ####

comp4 <- comp3 %>%
  mutate(recov_categ = cut(Current_CWT, breaks = c(-0.0001, 1, 5, 10, 20, 50, 100, 500, 2000)))

ggplot(comp4, aes(x = Current_CWT, y = Mark_and_Tag_All_ad_clip_count_const, colour = p_marked * f_marked)) + 
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + 
  facet_wrap(~ recov_categ, ncol = 2, scales = "free") +
  scale_colour_gradientn(colours = rev(rainbow(7)))




#### Now, what happens if we can increase the mark rate up to the inverse of the max(theta)? ####

theta_g_means <- tbl_df(plyr::ldply(post_means, 
                                    function(x)
                                      data.frame(tag_code = names(x$theta_gs), 
                                                 theta = x$theta_gs), 
                                    .id = "recovery_area"))

theta_g_maxes <- theta_g_means %>%
  group_by(tag_code) %>%
  summarise(max_theta = max(theta)) %>%


# but once I have these I have to come up with a good value of the theta below which 
# we will make all the mark-tag fractions 1, and above which we will make them 
# inversely proportional to theta.  It isn't clear how that value should be set.  
