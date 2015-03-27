


# get the number of fish released in the wild with CWTs
source("R-main/01-prepare-for-mcmc.R")  # source this to get the "releases" variable

# order states as desired
releases$release_location_state <- factor(releases$release_location_state, levels = c("AK", "BC", "WA", "ID", "OR", "CA"))


#### Get the number of wild and hatchery fish in release groups by state and year ####
# group it
tmp <- releases %>%
  filter(rearing_type %in% c("W", "M"), brood_year >= 2005, species %in% 1:2) %>%
  group_by(species, release_location_state, brood_year, rearing_type) %>%
  select(starts_with("n_"))

# query and inner join results
num_fish <-  tmp %>% summarise_each(funs(sum))
num_releases <-  tmp %>% summarise(num_release_groups = n())
inner_join(num_releases, num_fish)

by_state_year <- inner_join(num_releases, num_fish)


by_state_year %>%
  group_by(species, brood_year) %>%
  summarise_each(funs(sum), starts_with("n_")) %>%
  filter(brood_year != 2012) %>%
  summarise_each(funs(mean), starts_with("n_"))



#### Now, count up the number of different tag_codes that don't look functionally different ####

# Here I am going to couunt up the number of tag_codes that have the same agency, state, region, 
# hatchery, run, release location and release date, just to see how many of these release groups
# are probably not different (or are only different because they were treated separately within
# the hatchery for some time, and hence could be kept separate from the family stage anyway.)
R <- releases


# Filter the release data to only chinook (1) and coho (2), and only the brood years from 2005 to 2010.  
# And only release groups that have actual tag_codes (record_code == "T") and are only from hatcheries.
R_filt <- R %>% 
  filter(species %in% 1:2, brood_year %in% 2005:2010, record_code == "T", rearing_type == "H")


R_filt_group <- R_filt %>%
  group_by(record_code, release_agency, species, run, brood_year, hatchery_location_code, stock_location_code, release_stage, release_location_code, first_release_date, last_release_date, release_location_state, hatchery_location_name, stock_location_name, release_location_name, related_group_type)

# here is the total number of tag_codes
R_filt %>% group_by(species) %>% tally()

# here we count up the number of tag codes in each functionally equivalent release group
num_tag_codes_by_functional_release_group <- R_filt_group %>%
  summarise(num_tag_codes = n())

# check that the numbers are right
num_tag_codes_by_functional_release_group %>%
  group_by(species) %>%
  summarise(tot = sum(num_tag_codes))


# now, add a column called "num_funct_equiv_release_groups" that is 1 for all that
num_tag_codes_by_functional_release_group$num_funct_equiv_release_groups <- 1


# now, count those up (as well as the number of redundant tag codes by hatchery)
final_summary <- num_tag_codes_by_functional_release_group %>% 
  group_by(species, brood_year, release_location_state, release_agency, hatchery_location_code, hatchery_location_name) %>%
  summarise(num_funct_equiv_rel_groups = sum(num_funct_equiv_release_groups), num_tag_codes = sum(num_tag_codes))


# check that the total number of tag codes is correct
final_summary %>% group_by(species) %>% summarise(tot = sum(num_tag_codes))


# write out a csv
write.csv(final_summary, file = "num_functional_release_groups_per_hatchery_DIT_accounted.csv")

# make a histogram
final_summary$SpeciesName = c("Chinook", "Coho")[final_summary$species]
ggplot(final_summary, aes(x = num_funct_equiv_rel_groups, fill = release_location_state)) +
  geom_histogram(breaks = seq(0.5, max(final_summary$num_funct_equiv_rel_groups) + 0.5, by = 1)) + 
  facet_wrap(~ SpeciesName, ncol = 2) +
  xlab("Number of Release-Equivalent Release Groups Per Hatchery Per Brood Year")

ggsave(file = "num_fe_rel_groups_hist_DIT_accounted.pdf", width = 14, height = 10)



#### Explore Sampling Rates  #####
cs <- readRDS("./data/catch_sample.rds")

cs %>%
  filter(sample_type == 1, species %in% 1:2, catch_year == 2012) %>%
  mutate(is_escapement_or_weird = fishery >= 50) %>%
  group_by(sampling_agency, species, is_escapement_or_weird) %>%
  select(number_caught, number_sampled, number_estimated) %>%
  filter(!(is.na(number_caught)), !is.na(number_sampled)) %>%
  mutate(samp_fract = number_sampled / number_caught) %>%
  filter(samp_fract <= 1.00) %>%
  summarise(mean_samp_fract = sum(number_sampled) / sum(number_caught), num_sampled = sum(number_sampled)) %>%
  ungroup %>%
  arrange(desc(num_sampled))

# the upshot is that the average sampling rate of most (non-escapement) fisheries
# is between 0.2 and 0.34.


marine2012 <- all_recov %>% 
  filter(run_year == 2012,   # filter it down to 2012
         str_sub(recovery_location_code, 2, 2) == "M" #  only take marine recoveries here
  )

# look at this.  It is the histogram of expansion factors for tagged fish
marine2012 %>%
  select(estimated_number)  %>%
  filter(estimated_number > 1) %>%
  ggplot(., aes(x = 1/estimated_number)) + geom_histogram()

# so, if it is .25 it means that for every tagged fish from the reporting group there
# were 4 TAGGED fish total from the release group that were expected to be in the catch.

# So, how about if you wanted to make claims about the total number of fish from the release
# group that were expected to be in the catch, whether they were tagged or not?  It depends
# on how fish were sampled.  If they were sampled visually for ad-clips, then you would have to
# further expand this number by the fraction of all the fish in the release group that had
# ad-clips.  If the fishery was electronically detected, then we would want to expand by the fraction
# of the release group that had coded wire tags.  We can figure this all out quite simply with a few
# inner joins. 

# first, I need to get a mark_and_tag that has coho in it:
coho_releases_too <- group_releases_hierarchically(BroodYears = 2004:2012, Species = 2)$Release_ID
tmp <- rbind(coho_releases_too, releases2004_2012$Release_ID)
# but I don't use that yet...


mark_and_tag_rates <- distinct_codes %>%
  rename(tag_code = tag_code_or_release_id) %>%
  mutate(n_total_fish = n_tag_ad + n_tag_noad + n_notag_ad + n_notag_noad,
         n_marked = n_tag_ad + n_notag_ad,
         n_unmarked = n_tag_noad + n_notag_noad,
         f_marked = n_marked / n_total_fish,
         f_unmarked = n_unmarked / n_total_fish,
         p_marked = ifelse(n_marked > 0, n_tag_ad / n_marked, 0),
         p_unmarked = ifelse(n_unmarked > 0, n_tag_noad / n_unmarked, 0)
  )  %>% 
  ungroup() %>%
  mutate(release_location_state = factor(release_location_state, levels = c("AK", "BC", "WA", "ID", "OR", "CA")))


# so, let's do it:
tmp <- marine2012 %>%
  select(tag_code, detection_method, sampling_agency, estimated_number)  %>%
  filter(estimated_number > 1) %>%  # this seems reasonable.  It has to at least represent one recovery!
  inner_join(., mark_and_tag_rates) %>%
  inner_join(.,
    releases %>%
      select(tag_code_or_release_id, related_group_type),
    by = c("tag_code" = "tag_code_or_release_id")) %>% 
  filter(!(related_group_type %in% c("D", "O"))) 
  
# now we are ready to compute the necessary expansion factors.  If estimated_number is
# the number of tagged fish IN THE CATCH represented by the tagged fish that were recovered,
# then to get the total number of fish from the release group IN THE CATCH, we just have to
# multiply estimated_number by 1 over the fraction of fish in the release group that were
# tagged.  The fraction tagged is f_marked * p_marked + f_unmarked * p_unmarked 
exp_factos <- tmp %>%
  mutate(fraction_tagged = f_marked * p_marked + f_unmarked * p_unmarked,
         total_expansion_factor = fraction_tagged / estimated_number)

exp_factos$sampling_agency[!(exp_factos$sampling_agency %in% c("ADFG", "CDFO", "WDFW","ODFW", "CDFW"))] <- "OTHER"
exp_factos$sampling_agency <- factor(exp_factos$sampling_agency,
                                     levels = c("ADFG", "CDFO", "WDFW","ODFW", "CDFW", "OTHER"))

ggplot(exp_factos, aes(x = total_expansion_factor, fill = release_location_state)) + 
  geom_histogram(binwidth = 0.04) +
  facet_wrap(~ sampling_agency, ncol = 3) +
  xlab("m")

ggsave(file = "m_histo_chinook.pdf", width = 9, height = 6)

# NOTE TO SELF, FILTER OUT THE DIT GROUPS



##### Make a figure of family tagging rate  given genotyping success rate ####
vals <- seq(0, 1.0, by = .002)
succ_rate <- rbind(
  data.frame(parentage_method = "Parent Pairs", individual_genotyping_success_rate = vals, family_tagging_rate = vals^2),
  data.frame(parentage_method = "Single Parent", individual_genotyping_success_rate = vals, family_tagging_rate = 1 - (1 - vals)^2)
)


ggplot(succ_rate, aes(x = individual_genotyping_success_rate, 
                      y = family_tagging_rate,
                      colour = parentage_method)) +
  geom_vline(xintercept = c(.9843, .9762), colour = "blue") +
  geom_line()

ggsave(file = "succ_rate.pdf", width = 9, height = 6)
