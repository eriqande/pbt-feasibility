library(ggplot2)
library(stringr)
library(lubridate)
library(dplyr)
library(reshape2)


#### Define Functions ####
#' function to simulate one rep for pbt tag rate variance simulation study
#' 
#' @param  pars A named list of the parameter values:
#'    S   number of families
#'    G   number of families successfully genotyped
#'    mu  mean number of juveniles per family at release stage
#'    r   dispersion parameter of neg binom for number of juvies at release
#'    A_F adult fraction: Take floor(2 * S * A_F) to get the number of adults at the time of adult samlpling
#'    a   the parameter of the Dirichlet distribution controlling variance in family-specific survival rate
pbt_tag_rate_rep <- function(pars) {
  with(pars,{
    
    # set the number of adults
    N_A <- floor(2 * S * A_F)
    # deal with equivalent CWT stuff
    p_cwt = G/S   # the comparable CWT tag rate
    p_cwt_a = rbinom(n = 1, size = N_A, prob = p_cwt) / N_A # effective tagged fraction of adults
    
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
    
    # now, we estimate the effective sizes
    # for inbreeding effective size. 
    tmp <- rmultinom(n = 1, size = 2 * S, prob = q) # simulate offspring back as if const pop
    pIBD <- 0.25 * sum( (tmp / sum(tmp)) * (tmp-1) / (sum(tmp)-1) )
    NeI <- 1/pIBD * 0.5
    
    # and now get the standard deviation of family size amongst families
    # that have at least one offspring if each family is expected to produce
    # 2 offspring.  Do this by simulating a lot of families to reduce variance
    #tmp <- rmultinom(n = 1000, size = 2 * S, prob = q)
    #non_zero_fam_sd <- sd(tmp[tmp>0])
    #non_zero_fam_mean <- mean(tmp[tmp>0])
    
    
    as.data.frame(c(pars,
                    #non_zero_fam_sd = non_zero_fam_sd, non_zero_fam_mean = non_zero_fam_mean,
                    N_A = N_A, J_sd = sd(J), p_cwt = p_cwt, p_cwt_a = p_cwt_a, p_pbt = p_pbt, 
                    p_pbt_a = p_pbt_a, pIBD = pIBD, NeI = NeI))
  })
}


#### Run the set of simulations at 24 combinations of values ####

# here are the values to run things at:
NumReps <- 1000
Svals <- c(10, 25, 50, 100, 250, 1000)
A_Fvals <- 5
Gprops <- c(.8, .96)
disp_vals <- list(
  WF = c(mu = 3000, r = 10^8, a = 2000),
  Ne50 = c(mu = 3000, r = 30, a = 1)
  )

# now just for-loop over them:
out_list <- list()
rep_num <- 1
for(S in Svals) for(A_F in A_Fvals) for(Gp in Gprops) for(disp in disp_vals) {
  spars <- list(
      S = S,
      G = floor(S * Gp),
      mu = unname(disp["mu"]),
      r = unname(disp["r"]),
      A_F = A_F,
      a = unname(disp["a"])
    )
    
  out_list[[rep_num]] <- do.call(rbind, lapply(1:NumReps, function(x) pbt_tag_rate_rep(spars)))
  rep_num <- rep_num + 1
}

combo24 <- do.call(rbind, out_list)


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

# filter it down to tagged chinook from brood years >=2000, and pick out just the
# columns with number of fish released (and agency)
relcnts <-releases %>% 
  filter(record_code == "T", brood_year >= 2000) %>% 
  select(reporting_agency, ends_with("count"), brood_year, tag_code_or_release_id)

# change NAs in the counts to 0
relcnts[, 2:5][is.na(relcnts[,2:5])] <- 0

# count up the equivalent number of families that would have produced each
# tag group if they were PBT-tagged.  
parcnts <- relcnts %>% 
  mutate(total_smolts = cwt_1st_mark_count + cwt_2nd_mark_count + non_cwt_1st_mark_count + non_cwt_2nd_mark_count) %>%
  select(tag_code_or_release_id, reporting_agency, brood_year, total_smolts) %>%
  mutate(total_parents = total_smolts / 3000) %>%
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

write.table(tab, sep = "  &  ", eol = "\\\\\n", row.names = FALSE, quote = FALSE)



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
