#### This script takes the release data base, then adds necessary information 
#### from the locations data base, and summarizes all that information into
#### a table useful for futher summaries and aggregations

#### Script by Eric Anderson



#### load libraries, etc  ####
library(ggplot2)
library(lubridate)
library(dplyr)
library(reshape2)
library(stringr)
source("R/rmis_cleaning.R")



#### get the entire releases data base  ####
releases <- readRDS("data/releases.rds")



#### Some useful helper function #####

# return true if code implies ad-clip and false otherwise (including if NA)
has_adclip <- function(x) {
  ifelse(!is.na(x) & str_sub(as.character(x), 1, 1) == 5,  TRUE, FALSE)
}



#### For each release group count up fish of different categories ####

# turn NA's in the mark_count columns into 0's for summing
count_cols <- c("cwt_1st_mark_count", "cwt_2nd_mark_count",
                "non_cwt_1st_mark_count", "non_cwt_2nd_mark_count")
releases[count_cols] <- lapply(count_cols, function(x) ifelse(is.na(releases[[x]]), 0, releases[[x]]))


# turn blanks in the release_location_rmis_region and release_location_rmis_basin into "Unspecified"
releases$release_location_rmis_region[releases$release_location_rmis_region == ""] <- "Unspecified"
releases$release_location_rmis_basin[releases$release_location_rmis_basin == ""] <- "Unspecified"


# this operation does a few things
# 1. make a column "tag_variety" that says whether the record is for a release with cwts, pseudo tags, or neither cwt or pseudo tags
# 2. create the columns:
#     n_tag_ad   number tagged and adclipped
#     n_tag_noad     number tagged but not adclipped, 
#     n_notag_ad      etc.
#     n_notag_noad
# NOTE: The pseudo-tagged and AWT fish end up in the notag columns
releases <- releases %>%
  mutate(tag_variety = ifelse(is.na(tag_type), "none", ifelse(tag_type==16, "pseudo", "cwt"))) %>%
  mutate(n_tag_ad = has_adclip(cwt_1st_mark) * cwt_1st_mark_count + has_adclip(cwt_2nd_mark) * cwt_2nd_mark_count,
         n_tag_noad = (1 - has_adclip(cwt_1st_mark)) * cwt_1st_mark_count + (1 - has_adclip(cwt_2nd_mark)) * cwt_2nd_mark_count,
         n_notag_ad = has_adclip(non_cwt_1st_mark) * non_cwt_1st_mark_count + has_adclip(non_cwt_2nd_mark) * non_cwt_2nd_mark_count,
         n_notag_noad = (1 - has_adclip(non_cwt_1st_mark)) * non_cwt_1st_mark_count + (1 - has_adclip(non_cwt_2nd_mark)) * non_cwt_2nd_mark_count
  )




#### a function to compute and store summaries for different levels of aggregation ####

# We want list of data frames, each one grouped by different (hierarchical) variables.
# We disregard "run" because there is a lot of missing data in it.
# Two variables that we want to include at every level of the hierarchy are: brood_year and tag_variety
# the others nest naturally as state -> region -> basin -> release_id

# So, we do this as a function which takes species and years as arguments to filter stuff 
# out first, then returns a list of all the differently grouped and summarised data.
# Note that we toss out release codes that have NA for state. (there are 29 of those across all years...)
group_releases_hierarchically <- function(Species = 1, BroodYears = 2000:2004) {
  
  # here are the different hierarchical levels at which we will aggregate;
  aggs <- list(State = list(~ release_location_state, ~ brood_year, ~ tag_variety),
               Region = list(~ release_location_state, ~ release_location_rmis_region, ~ brood_year,  ~ tag_variety),
               Basin = list(~ release_location_state, ~ release_location_rmis_region, ~ release_location_rmis_basin, ~ brood_year,  ~ tag_variety),
               Release_ID = list(~ release_location_state, ~ release_location_rmis_region, ~ release_location_rmis_basin,  ~ brood_year, ~ tag_variety, ~ tag_code_or_release_id)
  )
  
  # now, summarize at those different aggregations:
  
  grouped_list <- lapply(aggs, function(agg) {
    # make a grouped tbl_df
    rel_grouped <- releases %>% 
      filter(species == Species, brood_year %in% BroodYears, !is.na(release_location_state)) %>%   
      group_by_(.dots = agg)
    
    # then summarise it.  I didn't know how to sum each column and put n_distinct on there too, so I throw in the left join
    left_join(
      rel_grouped %>%
        select(starts_with("n_")) %>% 
        summarise_each(funs(sum)),
      
      rel_grouped %>%
        summarise(num_tag_codes = n_distinct(tag_code_or_release_id))
    )
  })
}



#### Then, for any particular species and years we can do this:

  
try_it <- group_releases_hierarchically()


