


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
  group_by(record_code, release_agency, species, run, brood_year, hatchery_location_code, stock_location_code, release_stage, release_location_code, first_release_date, last_release_date, release_location_state, hatchery_location_name, stock_location_name, release_location_name)

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
write.csv(final_summary, file = "num_functional_release_groups_per_hatchery.csv")

# make a histogram
ggplot(final_summary, aes(x = num_funct_equiv_rel_groups - 1, fill = release_location_state)) +
  geom_histogram(binwidth = 1) + 
  facet_wrap(~ species, ncol = 2)

ggsave(file = "num_fe_rel_groups_hist.pdf", width = 14, height = 10)
