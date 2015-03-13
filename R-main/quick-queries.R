


# get the number of fish released in the wild with CWTs
source("R-main/01-prepare-for-mcmc.R")  # source this to get the "releases" variable

# order states as desired
releases$release_location_state <- factor(releases$release_location_state, levels = c("AK", "BC", "WA", "ID", "OR", "CA"))

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