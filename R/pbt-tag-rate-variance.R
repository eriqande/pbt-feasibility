library(ggplot2)
library(stringr)
library(lubridate)
library(dplyr)
library(reshape2)
library(parallel)


#### Define Functions ####
#' function to simulate one rep for pbt tag rate variance simulation study
#' 
#' @param  pars A named list of the parameter values:
#' @param   S   number of families
#' @param    G   number of families successfully genotyped
#' @param    mu  mean number of juveniles per family at release stage
#' @param    r   dispersion parameter of neg binom for number of juvies at release
#' @param    A_F adult fraction: Take floor(2 * S * A_F) to get the number of adults at the time of adult samlpling
#' @param    a   the parameter of the Dirichlet distribution controlling variance in family-specific survival rate
#' @param    m   the tag sample rate. (Probability that a tagged (or should-be-tagged fish will be sampled). 
#' In the right up this is c * v
#' @param    n_fish_present  The number of fish to consider present in a hypothetical fishery that we 
#' will sample marked fish from.
#' @param    N_C  total number of fish from the release group in the catch
pbt_tag_rate_rep <- function(pars) {
  with(pars,{
    
    # set the number of adults
    N_A <- floor(2 * S * A_F)
    if(N_C > 0) {  # some stuff here to set N_C for simulation number 2.
      N_A <- N_C
    }
    # deal with equivalent CWT stuff
    p_cwt = G/S   # the comparable CWT tag rate
    n_cwt_a <- rbinom(n = 1, size = N_A, prob = p_cwt)
    p_cwt_a <- n_cwt_a / N_A # effective tagged fraction of adults
    
    # this is number of cwts recovered in simulation 2
    n_cwt2 <- rbinom(n = 1, size = N_C, prob = m)
    N_hat_cwt <- n_cwt2 / m
    
    # and here we want to get the correct (but never knowable) expansion factor for marking at a rate m
    # (in a case where p_cwt of the marked fish are tagged)
    p_cwt_am <- rbinom(n = 1, n_cwt_a, prob = m) / N_A
    
    # and we use that mark-and-tag rate to simulate some CWT recoveries from n_fish_present fish
    # in the fishery.  cwt_mt_recovered =  CWT-marked-tagged-recovered
    cwt_mt_recovered <- rbinom(n = 1, n_fish_present, prob = p_cwt_am)
    
    # simulate the PBT parts
    J <- rnbinom(n = S, mu = mu, size = r)   # family sizes at release stage
    B <- sample(x = 1:S, size = G)  # indices of families successfully genotyped
    p_pbt <- sum(J[B]) / sum(J)  #  PBT tagging rate of juvies at release
    
    # Adult survival part
    tmp <- rgamma(n = S, shape = a, scale = 1.0)  # gammas for dirichlet
    Y <- tmp / sum(tmp)  # the Dirichlet r.v.'s
    q <- Y * J / sum(Y * J)
    W <- rmultinom(n = 1, size = N_A, prob = q)
    p_pbt_a = sum(W[B]) / sum(W)
    
    # now we can simulate p_pbt2 for simulation #2
    n_pbt2 <- rbinom(n = 1, size = N_C, prob = p_pbt_a * m)
    N_hat_pbt <- n_pbt2 * S / (m * G)
    
    # now, we should simulate the observed fraction that are marked and tagged out of all of
    # the fish that are there, because that is the number that one would expand by to get the
    # total number of fish after expansion.
    p_pbt_am = rbinom(n = 1, sum(W[B]), prob = m) / sum(W)
    
    # and, now, given that tag-and-mark rate, we could simulate n_fish_present fish being
    # "present" in a fishery, and sample and obtain tags from the marked ones
    # and record that, and then expand it and see how close we get to n_fish_present
    pbt_mt_recovered <- rbinom(n = 1, n_fish_present, prob = p_pbt_am)
    
    # now, we estimate the effective sizes
    # for inbreeding effective size. 
    tmp <- rmultinom(n = 1, size = 2 * S, prob = q) # simulate offspring back as if const pop
    pIBD <- 0.25 * sum( (tmp / sum(tmp)) * (tmp-1) / (sum(tmp)-1) )
    NeI <- 1/pIBD * 0.5
    
    # and now we consider what the actual expanded number of impacts would be assuming a
    # mark rate of m and catching 10 marked fish with tags
    
    
    as.data.frame(c(pars,
                    #non_zero_fam_sd = non_zero_fam_sd, non_zero_fam_mean = non_zero_fam_mean,
                    N_A = N_A, N_C = N_C, J_sd = sd(J), p_cwt = p_cwt, p_cwt_a = p_cwt_a, p_pbt = p_pbt, 
                    p_pbt_a = p_pbt_a, pIBD = pIBD, NeI = NeI, p_pbt_am = p_pbt_am,
                    p_cwt_am = p_cwt_am, n_fish_present = n_fish_present, cwt_mt_recovered = cwt_mt_recovered,
                    pbt_mt_recovered = pbt_mt_recovered,
                    n_pbt2 = n_pbt2, n_cwt2 = n_cwt2,
                    N_hat_pbt = N_hat_pbt, N_hat_cwt = N_hat_cwt))
  })
}


#### Run the set of simulations at 24 combinations of values ####

# here are the values to run things at:
NumReps <- 20
Svals <- c(10, 20, 30, 50, 100, 200, 400, 1000)
A_Fvals <- 5
N_Cvals <- 50
Gprops <- seq(0.8, 1.0, by = 0.02)
disp_vals <- lapply(c(0.3, 0.5, 1.0, 2, 10), function(x) c(mu = 3000, r = 30, a = x))
disp_vals[[length(disp_vals) + 1]] <- c(mu = 3000, r = 10^8, a = 10000) # add the (nearly) wright-fisher case in there
mVals <- c(0.125, 0.25, 0.5, 0.75, 1.0)
NFP_vals <- 50


# Expand grid and mclapply over all the combinations, so we can do it in parallel:
combos <- expand.grid(S = Svals, A_F = A_Fvals, Gp = Gprops, disp = disp_vals, m = mVals,
                      n_fish_present = NFP_vals, N_C = N_Cvals)
set.seed(15)
big_list <- mclapply(1:nrow(combos), function(x) {
  spars <- list(
    S = combos[x, "S"],
    G = floor(combos[x, "S"] * combos[x, "Gp"]),
    mu = unname(combos[x, "disp"][[1]]["mu"]),
    r = unname(combos[x, "disp"][[1]]["r"]),
    A_F = combos[x, "A_F"],
    a = unname(combos[x, "disp"][[1]]["a"]),
    m = combos[x, "m"],
    N_C = combos[x, "N_C"],
    n_fish_present = combos[x, "n_fish_present"]
  )
  do.call(rbind, lapply(1:NumReps, function(y) pbt_tag_rate_rep(spars)))
},
mc.cores = 8)

big_frame <- tbl_df(do.call(rbind, big_list)) 

# first, make an NeN_factor column that takes that mean Ne over N for each a and calls it a factor
big_frame2 <- big_frame %>%
  group_by(a,S) %>%
  summarise(tmp = (1/(2*mean(pIBD))) / (2 * mean(S))) %>%
  summarise(Ne_over_N = mean(tmp)) %>%
  inner_join(big_frame, .) %>%
  mutate(NeN_factor = factor(sprintf("%.2f", Ne_over_N)))


# aha! I need to just put all the cwt ones, regardless of NeN factor into the pbt_a and pbt_am columns
# with the NeN factor = "CWT".  However, I don't want to put all of them there, because
# then they end up having too many outliers.  Just filter it down to ones that have NeN factor = 1.00
big_frame3 <- big_frame2 %>%
  filter(NeN_factor == "1.00") %>%
  group_by(S, G, N_C, m) %>%
  select(n_cwt2, N_hat_cwt) %>%
  rename(n_pbt2 = n_cwt2, N_hat_pbt = N_hat_cwt) %>%
  mutate(NeN_factor = "CWT") %>%
  rbind_list(., big_frame2)  %>%
  rename(N_hat = N_hat_pbt, Nb_over_2S = NeN_factor)

# %>%
#   mutate(expando_pbt = 10/p_pbt_am) %>%
#   mutate(expando_estimate = pbt_mt_recovered / (m * (G/S)))


#### Now, compute the mean and standard deviation of the expando_estimates ####
sds_and_means <- big_frame3 %>%
  group_by(S, G, N_C, Nb_over_2S, m) %>%
  summarise(mean_N_hat = mean(N_hat), sd_of_N_hat = sd(N_hat))


# Let's make a similar plot, but put lines top and bottom of different colors on the
# ends of the standard deviation range.
lapply(unique(sds_and_means$m), function(x) {
  
  tmp <- sds_and_means %>% filter(m == x) 
  nc_static <- 50  # this would have to be changed if N_C were changed
  ggplot(tmp, aes(x = G/S, y = N_hat, colour = Nb_over_2S)) + 
    geom_abline(intercept = nc_static, slope = 0) + # put this on the bottom cuz CWT will overwrite it.
    geom_line(aes(y = N_C - 2 * sd_of_N_hat)) +
    geom_line(aes(y = N_C + 2 * sd_of_N_hat)) +
    geom_point(aes(y = N_C - 2 * sd_of_N_hat), size = 1.4) +
    geom_point(aes(y = N_C + 2 * sd_of_N_hat), size = 1.4) +
    geom_point(aes(y = mean_N_hat), size = 1.6) +
    ylim(0, 100) +
    facet_wrap(~ S, ncol = 4) +
    scale_color_manual(values = rainbow(8))
  
  ggsave(paste("sd_line_horns_m_", x, ".pdf", sep = ""), width = 14, height = 10)
})



########   BONEYARD #######


# just some testing
bf_small <- big_frame3 %>% filter(S == 10, NeN_factor %in% c("0.34", "0.50", "0.66", "0.81", "1.00", "CWT"))
ggplot(bf_small, aes(x = factor(G/S), y = p_pbt_a, fill = NeN_factor)) + geom_boxplot(outlier.size = 1) 

ggplot(bf_small, aes(x = factor(G/S), y = expando_pbt, fill = NeN_factor)) + 
  geom_violin() +
  #geom_boxplot(outlier.size = 1) +
  expand_limits(y = 0)



# let's see if we can get all the standard deviations for the realized tagging rates
# onto a single page.  Note that m is not relevant to this.
tag_rate_var <- big_frame3 %>%
  mutate(G_over_S = G/S) %>%
  group_by(G_over_S, S, NeN_factor) %>%
  summarise(tag_rate_sd = sd(p_pbt_a))

ggplot(tag_rate_var, aes(x = G_over_S, y = tag_rate_sd, colour = NeN_factor)) +
  geom_line() +
  geom_point(size = 1.5) +
  facet_wrap(~ S, nrow = 5)

ggsave(file = "tag_rate_standard_devs.pdf", width = 14, height = 10)


lapply(unique(big_frame3$S), function(x) {
  bf <- big_frame3 %>% filter(S == x)
  ggplot(bf, aes(x = factor(G/S), y = p_pbt_a, fill = NeN_factor)) + 
  geom_boxplot(outlier.size = 1)
  
  ggsave(file = paste("pbt_var_boxplots_", x, ".pdf", sep = ""), width = 14, height = 10)
})


lapply(unique(big_frame3$S), function(x) {
  bf <- big_frame3 %>% filter(S == x)
  ggplot(bf, aes(x = factor(G/S), y = expando_estimate, fill = NeN_factor)) + 
    geom_boxplot(outlier.size = 1) +
    expand_limits(y = 0) + 
    facet_wrap(~ m, ncol = 2)
  
  ggsave(file = paste("pbt_expando_boxplots_", x, ".pdf", sep = ""), width = 14, height = 10)
})



#### Now, compute the mean and standard deviation of the expando_estimates ####
sds_and_means <- big_frame3 %>%
  group_by(S, G, NeN_factor, m) %>%
  summarise(mean_N_hat = mean(N_hat_pbt), sd_of_N_hat = sd(N_hat_pbt))


# Let's make a similar plot, but put lines top and bottom of different colors on the
# ends of the standard deviation range.
lapply(unique(sds_and_means$m), function(x) {
  
  tmp <- sds_and_means %>% filter(m == x) 
  
  ggplot(tmp, aes(x = G/S, y = sd_expansion, colour = NeN_factor)) + 
    geom_line(aes(y = 50 - 2 * sd_expansion)) +
    geom_line(aes(y = 50 + 2 * sd_expansion)) +
    ylim(0, 100) +
    geom_abline(intercept = 50, slope = 0) +
    facet_wrap(~ S, ncol = 5) 
  
  ggsave(paste("sd_line_horns_m_", x, ".pdf", sep = ""), width = 14, height = 10)
})

# and make a nice line plot
lapply(unique(sds_and_means$m), function(x) {
  
  tmp <- sds_and_means %>% filter(m == x) 
  
  ggplot(tmp, aes(x = G/S, y = sd_expansion, colour = NeN_factor)) + 
    geom_line() +
    facet_wrap(~ S, ncol = 5)
  
  ggsave(paste("standard_devs_m_", x, ".pdf", sep = ""), width = 14, height = 10)
})


# then, also show what it looks like if we  represent the standard deviation as 
# the true value +- 2 sds
lapply(unique(sds_and_means$m), function(x) {
  
  tmp <- sds_and_means %>% filter(m == x) 
  
  ggplot(tmp, aes(x = G/S, y = sd_expansion, colour = NeN_factor)) + 
    geom_linerange(aes(ymin = 50 - 2 * sd_expansion, ymax = 50 + 2 * sd_expansion),
                   position = position_dodge(width = .01)) +
    facet_wrap(~ S, ncol = 5)
  
  ggsave(paste("sd_lineranges_m_", x, ".pdf", sep = ""), width = 14, height = 10)
})



# Let's make a similar plot, but put lines top and bottom of different colors on the
# ends of the standard deviation range.
lapply(unique(sds_and_means$m), function(x) {
  
  tmp <- sds_and_means %>% filter(m == x) 
  
  ggplot(tmp, aes(x = G/S, y = sd_expansion, colour = NeN_factor)) + 
    geom_line(aes(y = 50 - 2 * sd_expansion)) +
    geom_line(aes(y = 50 + 2 * sd_expansion)) +
    ylim(0, 100) +
    geom_abline(intercept = 50, slope = 0) +
    facet_wrap(~ S, ncol = 5) 
  
  ggsave(paste("sd_line_horns_m_", x, ".pdf", sep = ""), width = 14, height = 10)
})





#### Process those simulations and make some plots/summaries ####
# combo24 <- readRDS("outputs/combo24.rds")  # do this if you don't want to re-do the sims
# add column for dispersion scenario
combo24 <- tbl_df(combo24)
combo24 <- combo24 %>%
  mutate(Dispersion_Scenario = ifelse(r>10000, "WrightFisher", "Ne50"))

# melt it so that we have p_pbt_a and p_cwt_a in long format
c24 <- tbl_df(melt(combo24, 
                   measure.vars = c("p_pbt", "p_pbt_a", "p_cwt", "p_cwt_a"),
                   value.name = "realized_tagged_fraction"))

# verify the effective sizes
combo24 %>%
  group_by(S, Dispersion_Scenario) %>%
  summarize(mean_pIBD = mean(pIBD), Ne = 1 / (2 * mean_pIBD))
# yep! they look right




# now we make four pages of plots.  Each is six panels (one for each S)
# and each page is a combination of dispersion scenario and G_fract

# first add a decent column for G_fract
c24 <- c24 %>% 
  mutate(G_fract = ifelse(G/S < .85, "80%", "96%"))


for(gf in c(80, 96)) for(ds in c("WrightFisher", "Ne50")) {
  
  cc <- c24 %>% 
    filter(G_fract == paste(gf, "%", sep = ""), 
           Dispersion_Scenario == ds,
           variable %in% c("p_pbt_a", "p_cwt_a")
    ) %>%
    mutate(floored_G = G/S) %>%
    droplevels
  
  
  
  g <- ggplot(cc, aes(x = realized_tagged_fraction, fill = variable)) + 
    geom_density(alpha = 0.5) + 
    facet_wrap(~ S, scales = "free", ncol = 2) + 
    geom_vline(aes(xintercept = floored_G)) + 
    scale_fill_manual(values = c("orange", "blue"))
  
  ggsave(filename =  paste("tag_fracts_", gf, "_", ds, ".pdf", sep=""), 
         plot = g,
         height = 8, 
         width = 10.5)
  
}


#### Count up how many release groups there are of different sizes ####

# get some necessary functions
source("R/rmis_cleaning.R")

# load it the releases data base
releases <- readRDS("data/releases.rds")
dim(releases)
releases$species <- int_to_factor(releases$species, species_levels())

make_rel_group_sizes_table <- function(Species = "Chinook", smolts_per_female = 3800) {
  # filter it down to tagged chinook from brood years >=2000, and pick out just the
  # columns with number of fish released (and agency)
  relcnts <-releases %>% 
    filter(record_code == "T", brood_year >= 2000, species == Species) %>% 
    select(reporting_agency, ends_with("count"), brood_year, tag_code_or_release_id)
  
  # change NAs in the counts to 0
  relcnts[, 2:5][is.na(relcnts[,2:5])] <- 0
  
  # count up the equivalent number of families that would have produced each
  # tag group if they were PBT-tagged.  
  parcnts <- relcnts %>% 
    mutate(total_smolts = cwt_1st_mark_count + cwt_2nd_mark_count + non_cwt_1st_mark_count + non_cwt_2nd_mark_count) %>%
    select(tag_code_or_release_id, reporting_agency, brood_year, total_smolts) %>%
    mutate(total_parents = total_smolts / smolts_per_female) %>%
    mutate(num_family_categories = cut(total_parents, breaks = c(0, 10, 25, 50, 100, 250, 500, 1000, 10^4)))
  
  # how many total fish released from each "number of parents category"
  # and also how many release groups?
  tab <- parcnts %>%
    group_by(num_family_categories) %>%
    summarise(tot_smolts = sum(total_smolts), tot_releases = n()) %>%
    mutate(rel_smolts = tot_smolts / sum(tot_smolts), rel_releases = tot_releases / sum(tot_releases))
  
  
  tab[2:3] <- lapply(tab[2:3], function(x) format(x, big.mark = ",", scientific = FALSE))
  tab[4:5] <- lapply(tab[4:5], function(x) sprintf("%.2f", 100 * x))
  tab$num_family_categories <- str_replace(tab$num_family_categories, pattern = "]", "") %>%
    str_replace(., "\\(", "") %>%
    str_replace(., ",", "--") %>%
    str_replace(., "1e\\+03", "1,000") %>% 
    str_replace(., "1e\\+04", "10,000")
  
  tab
}


CohoTable <- make_rel_group_sizes_table("Coho", 1800)
ChinookTable <- make_rel_group_sizes_table("Chinook", 3800)

CohoTable <- CohoTable[!is.na(CohoTable$num_family_categories), ]
ChinookTable <- ChinookTable[!is.na(ChinookTable$num_family_categories), ]

write.csv(CohoTable, file = "CohoTable.csv")
write.csv(ChinookTable, file = "ChinookTable.csv")


write.table(ChinookTable, sep = "  &  ", eol = "\\\\\n", row.names = FALSE, quote = FALSE)
write.table(CohoTable, sep = "  &  ", eol = "\\\\\n", row.names = FALSE, quote = FALSE)

#### LEFTOVERS!here is an example of running it ####
# here is an example pars:
pars50 <- list(
  S = 2500,
  G = 2400,
  mu = 3000,
  r = 30,
  A_F = 5,
  a = 1
)


# here is an example pars:
parsWF <- list(
  S = 25,
  G = 20,
  mu = 3000,
  r = 3000^2,
  A_F = 5,
  a = 1000
)

pars <- pars50
tmp <- do.call(rbind, lapply(1:500, function(x) pbt_tag_rate_rep(pars)))
summary(tmp)
(0.5 * 1 / mean(tmp$pIBD)) #/ pars$S 

# and here we plot it
ggplot(tmp, aes(x = p_cwt_a)) + 
  geom_density(fill = "orange", alpha = 0.5) + 
  geom_density(aes(x = p_pbt_a), fill = "blue", alpha = 0.5)



#### Here we get a normalized standard deviation of zero-censored family size

obs <- tbl_df(read.table("data/FRH_trios.txt", stringsAsFactors = FALSE)) %>%
  mutate(brood_year = str_match(pa, "_(2[0-9]{3})$")[,2])

# here we can get the means and sd of the families with at least one offspring:
obs %>% group_by(brood_year, pair) %>% summarise(num_offs = n()) %>% summarise(mean = mean(num_offs), sd = sd(num_offs))



#### Crunch out Mike Ackerman's Data to get the effective number of spawners at the Idaho Chinook hatcheries ####
mike <- read.csv("data/SteelheadRRSforEric.csv", skip = 2) %>% tbl_df

# clean up all the MS XL garbage and only keep columns we want
mike2 <- mike %>% 
  select(Dworshak:Wallowa) %>%
  filter(complete.cases(.))


mike_tot_spawn <- mike2[1,]

mike3 <- mike2[-(1:4),]

## Wait a minute! Now that I have read that in it is clear we want it by family!  So, we use what Craig sent instead.


#### Compute Nb/N ratios from Craig's data ####

# get the data:
craig <- read.csv("data/craig_chinook_fam_sizes.csv")

simpleNb <- function(n, i) {  
  tmp <- rep(i, n)
  pIBD <- 0.25 * sum( (tmp / sum(tmp)) * (tmp-1) / (sum(tmp)-1) )
  NeI <- 1/pIBD * 0.5
  NeI
}

simpleNb(craig$Sawtooth, craig$NumOffs)
