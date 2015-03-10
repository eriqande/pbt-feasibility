

#' simulate a dirichlet R.V.
#' 
#' @param a  The parameter vector
rdirichlet <- function(a) {
  y <- rgamma(n = length(a), shape = a)
  y / sum(y)
}


#' This function takes something like rec2012 and cs2 and distinct_codes and picks out the required recovery group it
#' 
#'  It picks out the recovery_group and then it summarizes the catch sample data as we need it.  The resulting list
#'  is what gets passed to the estimation function. 
#'
#' @param r the full data frame of recovery data
#' @param cs the full data frame of catch sample data 
#' @pamam fp the data frame that has the f and p variables (the mark and tag rates of the different tag codes) and it includes the distinct codes that correspond
#' to the theta_g's that we will estimate.  
#' @param rg the recovery group (a string) to pull out of r and cs
squash_data_for_estimation <- function(r, cs, fp, rg) {
  ret <- list()
  
  ## first prep the recovery data
  ret$recovery <- r %>% filter(recovery_group == rg)
  
  ## then summarise the catch-sample data
  tmp <- cs %>% 
    filter(recovery_group == rg) %>%
    group_by(detection_method) %>%
    summarise_each(funs(sum(., na.rm = TRUE)), vars = mr_1st_partition_size:mr_2nd_sample_obs_adclips)
  
  # this is just some silly stuff to get the names
  tmp_summary_names <- cs %>%
    ungroup %>%
    group_by(detection_method) %>%
    select(mr_1st_partition_size:mr_2nd_sample_obs_adclips) %>%
    names()
  
  names(tmp) <- tmp_summary_names
  
  ret$catch_sample <- tmp
  
  ## then get the mark and tag rates
  ret$mark_and_tag <- fp
  
  # and then return it all
  ret
}


#' compute an "expected visual-detection fishery sample" given the current parameter values
#' 
#' You can change the mark and tag rates for the stocks and see the effects of that, too.
#' @param n  The number of fish to visually sample for adclips
#' @param theta_gs  The mixing proportions of different tag groups
#' @param theta_uas the four other parameters
#' @param mark_and_tag dataframe with (at least) the columns f_marked, f_unmarked, p_marked, p_unmarked
viz_expect_predict <- function(n, theta_gs, theta_uas, mark_and_tag) {
  
  n <- unlist(n)  # make sure it is not part of a table any longer.
  
  # make short variables
  fm <- mark_and_tag$f_marked
  fu <- mark_and_tag$f_unmarked
  pm <- mark_and_tag$p_marked
  pu <- mark_and_tag$p_unmarked
  
  # This is the expected fraction of fish with adclips
  ad_fract <- sum(theta_gs * fm) + sum(theta_uas[c("U_plus", "A_plus")])
  ad_count <- ad_fract * n  # the expected number of fish with adclips
  
  cwt_fract <- sum(theta_gs * fm * pm) / ad_fract   # fraction of ad-clipped fish bearing cwts
  cwt_count <- ad_count * cwt_fract  # expected number of ad-clipped fish carrying cwts
  
  tag_group_counts <- cwt_count * (theta_gs * fm * pm) / sum(theta_gs * fm * pm)  # expecte number of cwts from each tag_group
  
  list(viz_samp_size = n,
       ad_fract = ad_fract,
       ad_count = ad_count,
       cwt_fract = cwt_fract,
       cwt_count = cwt_count,
       tag_group_counts = tag_group_counts)
                                        
}

# this is here to set values while developing the function below
if(FALSE) {
  r <- dat$recovery
  cs <- dat$catch_sample
  fp <- dat$mark_and_tag_rate
  #  rg <- "06-BC  (4480)"
  #  rg <- "16-CA  (4183)"
  rg <- "01-AK-NW  (2111), 02-AK-NW  (1194)"
  s <- squash_data_for_estimation(r, cs, fp, rg)
}


#' Do the Bayesian mixed fishery proportion estimation
#' 
#' Takes as input the output of squash_data_for_estimation
#' @param  s  The "squashed" list of input data
#' @param unclipped_awts_to_U_minus  For cases when there are few awt fish and you want to assume that
#' there are now unclipped awt fish.
#' @param no_ad_tag_ceiling.  This is a hack necessary in fisheries that don't have electronic detection (like California's) 
#' because the estimated fraction of some release groups that are not-clipped but are tagged at a high rate can go through
#' the roof because there are no data on whether unclipped fish have cwts.  When this happens, lots of unknown fish get
#' allocated to the no_adclip with a cwt pile.  With no_ad_tag_ceiling set, fish in that category are tossed into the
#' U_minus pile no more then no_ad_tag_ceiling * (# of observed adclipped cwt fish) is less than no_ad_tag_ceiling.  Set it
#' to a negative value to not use it at all. 
cwt_ppn_estimation_function <- function(s, reps = 10000, thin = 10,
                                        unclipped_awts_to_U_minus = TRUE,
                                        no_ad_tag_ceiling = 0.02) {
  # preliminary stuff, define variables, etc.
  
  ## here are the parameters of the models.  The proportions of fish from the different tag_groups
  theta_gs <- rep(NA, length.out = nrow(s$mark_and_tag))  # set them to NA for now.  We will initialize them later
  names(theta_gs) <- s$mark_and_tag$tag_code
  
  theta_uas <- rep(NA, 4)
  names(theta_uas) <- c("A_plus", "A_minus", "U_plus", "U_minus")
  
  ## here are the directly observed counts of fish corresponding to the tag_codes in theta_gs
  theta_g_obs_n <- rep(0, length.out = nrow(s$mark_and_tag))
  names(theta_g_obs_n) <- names(theta_gs)
  tmp <- s$recovery %>%
    filter(cwt_status == "cwt") %>%
    group_by(tag_code) %>%
    tally()
  theta_g_obs_n[tmp$tag_code] <- tmp$n
  
  # and here are the directly observed fish amongst the recoveries in the "A_plus", "A_minus", "U_plus", "U_minus" categories
  # note that we can't actually observe U_plus and U_minus directly, so we won't be setting them here
  theta_uas_n_obs <- rep(0, 4)
  names(theta_uas_n_obs) <- names(theta_uas)
  tmp <- s$recovery %>%
    filter(cwt_status == "awt") %>%
    group_by(ad_clipped) %>%
    tally()
  if(nrow(tmp) > 0) {
    if(length(tmp$n[tmp$ad_clipped == "no"]) > 0) { theta_uas_n_obs["A_minus"] <- tmp$n[tmp$ad_clipped == "no"]}
    if(length(tmp$n[tmp$ad_clipped == "yes"]) > 0) { theta_uas_n_obs["A_plus"] <- tmp$n[tmp$ad_clipped == "yes"]}
  }
  
  ## and now we sum up the numbers of fish that have cwt and ad-clips.  
  recovs_factored <- s$recovery
  recovs_factored$ad_clipped <- factor(recovs_factored$ad_clipped, levels = c("yes", "no", "unknown"))
  recovs_factored$cwt_status <- factor(recovs_factored$cwt_status, levels = c("cwt", "no_read", "awt", "no_tag", "unknown"))
  recov_sums <- as.data.frame(table(recovs_factored$ad_clipped, recovs_factored$cwt_status))
  names(recov_sums) <- c("ad_clipped", "cwt_status", "n")
  
  # for now, just deal with the visual recoveries.  We do that by just adding the number
  # of no-adclip fish from the catch-sample information into the no, unknown category
  recov_sums$n[recov_sums$ad_clipped == "no" & recov_sums$cwt_status == "unknown"] <- 
    s$catch_sample[s$catch_sample$detection_method == "V", "mr_1st_sample_known_ad_status"] - 
    s$catch_sample[s$catch_sample$detection_method == "V", "mr_1st_sample_obs_adclips"]
  
   
  ## and this is the general shape of the values that we can compute probabilities for
  adc_cwt <- recov_sums %>% filter(ad_clipped != "unknown", cwt_status != "unknown", cwt_status != "no_read")  # I have to deal with the no_read category somehow later...
  
  ## and for the ones that have unknown values, we will store those numbers in a list.
  tmp <- recov_sums %>% filter((ad_clipped == "unknown" | cwt_status == "unknown") & cwt_status != "no_read")
  vis_counts <- lapply(1:nrow(tmp), function(x) {  # I am naming those "vis_counts" because it will be appropriate for visually sampled fisheries
    ret <- list()
    ret$vars <- c(as.character(tmp[x,1]), as.character(tmp[x,2]))
    names(ret$vars) <- names(tmp)[1:2]
    ret$n <- unlist(tmp[x,3])
    ret
  })
  # those are the counts of individuals in these latent categories that include an "unknown"
  
  # and here, after that we can compute for later use the total number of fish sampled (that is all those
  # that were checked for an ad-clip in visual fisheries.) We will use this later for predictive checks
  tot_fish_with_known_ad_status <- unname(unlist(s$catch_sample[s$catch_sample$detection_method == "V", "mr_1st_sample_known_ad_status"]))
  
  ## Now, here we will initialize the theta_gs to something silly---just set it to the observed proportion, plus a little bit of prior...
  ## I am just doing this while developing...
  theta_gs <- theta_g_obs_n + 0.5  # don't bother normalizing yet, as we will have to normalize things altogether at the end
  
  # and also initialize the theta_uas.  I'm just gonna set A_plus and and A_minus according to the observation and then let U_minus be about 50% and see
  # if it can get away from that value in the course of the MCMC.  And U_plus to be about 5%
  theta_uas <- theta_uas_n_obs + 1
  tmp <- sum(theta_gs, theta_uas)
  theta_uas[c("U_minus", "U_plus" )] <- c(0.8 * tmp, 0.05 * tmp) 
  
  totsum <- sum(theta_gs, theta_uas)
  theta_gs <- theta_gs / totsum
  theta_uas <- theta_uas / totsum
  
  
  ### before doing mcmc, we need a function that will take elements of vis_counts, and simulate all the 
  ### individuals in them into the rows of an imputed version of adc_count.  vp is the vis_count_table.
  ### This thing essentially does the conditioning that needs to get done.
  impute <- function(x, vp) {
    p <- vp$probs
    if(x$vars["ad_clipped"] != "unknown") {
      p[ vp$ad_clipped != x$vars["ad_clipped"]] <- 0
    }
    if(x$vars["cwt_status"] != "unknown") {
      p[ vp$cwt_status != x$vars["cwt_status"]] <- 0
    }
    p <- p / sum(p) 
    
    # then return a multinomial allocation of all these inividuals into "fully-observed" categories
    rmultinom(1, size = x$n, prob = p)
    
    #list(x$vars, p)
  }
  
  
  ## Now,  for the iterations of the MCMC we need to be able to compute a lot of things
  ## based on tagging and marking rate.  It will be convenient to name some convenience variables here
  fm <- s$mark_and_tag$f_marked
  fu <- s$mark_and_tag$f_unmarked
  pm <- s$mark_and_tag$p_marked
  pu <- s$mark_and_tag$p_unmarked
  
  #### Now that everything is initialized and we have our convenience variables, and 
  #### It is time to start doing the iterations of the MCMC.  Before we get started.
  #### Let's create two lists to store the output.  
  ret <- list()
  ret$theta_gs <- vector("list", ceiling(reps / thin)) 
  ret$theta_uas <- vector("list", ceiling(reps / thin))
  ret$predictive_viscounts <- vector("list", ceiling(reps / thin))
  ret$viz_expect_predict <- vector("list", ceiling(reps / thin))
  
  # and put the starting values in there:
  ret$theta_gs[[1]] <- theta_gs
  ret$theta_uas[[1]] <- theta_uas
  
  # initialize an index for storing the results
  store_idx <- 1
  
  #### This block does the MCMC iterations ####
  for(iterations in 1:reps) {
    # compute our main probability table for visually sampled fisheries
    vis_samp_prob_table <- adc_cwt
    vis_samp_prob_table$probs <- c(
      sum(theta_gs * fm * pm),  #   yes, cwt  (<--- ad_clipped, cwt_status)
      sum(theta_gs * fu * pu),  #    no, cwt
      theta_uas["A_plus"],      #   yes, awt
      theta_uas["A_minus"],     #    no, awt
      sum(theta_gs * fm * (1 - pm)) + theta_uas["U_plus"], #   yes, no_tag
      sum(theta_gs * fu * (1 - pu)) + theta_uas["U_minus"] #    no, no_tag
    )
    
    predictive_table <- vis_samp_prob_table
    predictive_table$pred_n <- rmultinom(1, tot_fish_with_known_ad_status, predictive_table$probs)
    
    
    # then we can lapply that function, bind the results and rowSum them to get the allocations of those individuals.  Cool.
    alloc <- rowSums(do.call(cbind, lapply(vis_counts, function(x) impute(x, vis_samp_prob_table))))
    
    # look at what that looks like:
    cbind(vis_samp_prob_table, alloc)
    
    if(no_ad_tag_ceiling >= 0) {
      limit <- ceiling(unlist(vis_samp_prob_table$n[1]) * no_ad_tag_ceiling)
      if(alloc[2] > limit) {
        transfer <- alloc[2] - limit
        alloc[2] <- limit
        alloc[6] <- alloc[6] + transfer
      }
    }
    
    #### Now, we have to allocate some of those individuals
    #### further to individual tag codes (or to U_plus and U_minus), so that
    #### we can use that information to simulate a new value for theta_gs and theta_uas. I am just going to sort of hard-wire
    #### this at the moment.  
    theta_uas_add <- rep(0, 4)
    names(theta_uas_add) <- names(theta_uas)
    theta_gs_add_yes_cwt <- rmultinom(1, alloc[1], theta_gs * fm * pm)  # yes, cwt
    theta_gs_add_no_cwt <- rmultinom(1, alloc[2], theta_gs * fu * pu)   #  no, cwt
    theta_uas_add["A_plus"] <- alloc[3]   # yes, awt
    theta_uas_add["A_minus"] <- alloc[4]  #  no, awt
    
    # for the categories that include U_plus and U_minus with the theta_gs, we add the U category at the front of the probs
    # then normalize, then move to the multinomial pull that first element off where it needs to to
    # for yes, no_tag
    tmp <- rmultinom(1, alloc[5], c(theta_uas["U_plus"], theta_gs * fm * (1 - pm)))
    theta_uas_add["U_plus"] <- tmp[1]
    theta_gs_add_yes_no_tag <- tmp[-1] 
    
    # now for no, no_tag
    tmp <- rmultinom(1, alloc[6], c(theta_uas["U_minus"], theta_gs * fu * (1 - pu)))
    theta_uas_add["U_minus"] <- tmp[1]
    theta_gs_add_no_no_tag <- tmp[-1]
    
    # and now we also have to be VERY careful to add the observed ad_clipped fish with no_tag to
    # the pile of fish we are allocating to individual tag_codes.
    tmp <- rmultinom(1, unlist(adc_cwt[5,"n"]), c(theta_uas["U_plus"], theta_gs * fm * (1 - pm)))
    theta_uas_add["U_plus"] <- tmp[1]
    theta_gs_obs_yes_no_tag <- tmp[-1] 
    
    #### Now that those guys are all allocated, we just need to simulate new values of theta_gs and 
    #### theta_uas from their Dirichlet full conditionals.  We will assume simple unit-information priors
    #### but we will give more prior weight to the four uas categories, just cuz we can.  Shouldn't make
    #### much difference anyway
    theta_gs_dirichlet_pars <- c( 1/length(theta_gs) +
                                    theta_g_obs_n +
                                    theta_gs_add_yes_cwt +
                                    theta_gs_add_no_cwt +
                                    theta_gs_add_yes_no_tag + 
                                    theta_gs_add_no_no_tag +
                                    theta_gs_obs_yes_no_tag)
    
    theta_uas_dirichlet_pars <- c( 1/length(theta_uas) +
                                     theta_uas_n_obs +
                                     theta_uas_add)
    
    # now we simulate a dirichlet random vector and put the first four elements back into theta_uas and the
    # rest into theta_gs
    tmp <- rdirichlet(c(theta_uas_dirichlet_pars, theta_gs_dirichlet_pars))
    theta_uas[] <- tmp[1:4]  # preserve names with the []'s
    if(unclipped_awts_to_U_minus == TRUE) {
      theta_uas[4] <- theta_uas[2] + theta_uas[4]
      theta_uas[2] <- 0
    }
    theta_gs[] <- tmp[-(1:4)]
    
    # and store these new values in the output variable if appropriate
    if(iterations %% thin == 0) {
      store_idx <- store_idx + 1
      ret$theta_gs[[store_idx]] <- theta_gs
      ret$theta_uas[[store_idx]] <- theta_uas
      ret$pred_table[[store_idx]] <- predictive_table
      ret$viz_expect_predict[[store_idx]] <- viz_expect_predict(tot_fish_with_known_ad_status, 
                                                                theta_gs,
                                                                theta_uas,
                                                                s$mark_and_tag)
    }
    
  }
  
  ret
  
}
