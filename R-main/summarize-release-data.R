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
# 1. make a column "tag" that says whether the record is for a release with cwts, pseudo tags, or neither cwt or pseudo tags
# 2. create the columns:
#     n_tag_ad   number tagged and adclipped
#     n_tag_noad     number tagged but not adclipped, 
#     n_notag_ad      etc.
#     n_notag_noad
# NOTE: The pseudo-tagged and AWT fish end up in the notag columns
releases <- releases %>%
  mutate(tag = ifelse(is.na(tag_type), "none", ifelse(tag_type==16, "pseudo", "cwt"))) %>%
  mutate(n_tag_ad = has_adclip(cwt_1st_mark) * cwt_1st_mark_count + has_adclip(cwt_2nd_mark) * cwt_2nd_mark_count,
         n_tag_noad = (1 - has_adclip(cwt_1st_mark)) * cwt_1st_mark_count + (1 - has_adclip(cwt_2nd_mark)) * cwt_2nd_mark_count,
         n_notag_ad = has_adclip(non_cwt_1st_mark) * non_cwt_1st_mark_count + has_adclip(non_cwt_2nd_mark) * non_cwt_2nd_mark_count,
         n_notag_noad = (1 - has_adclip(non_cwt_1st_mark)) * non_cwt_1st_mark_count + (1 - has_adclip(non_cwt_2nd_mark)) * non_cwt_2nd_mark_count
  )



# here is a quick sanity check:
rel_grouped <- releases %>%
  filter(species == 1, brood_year == 2000) %>% 
  group_by(release_location_state, release_location_rmis_region, release_location_rmis_basin, run, tag)

left_join(
  rel_grouped %>%
    select(starts_with("n_")) %>% 
    summarise_each(funs(sum)),
  
  rel_grouped %>%
    summarise(num_tag_codes = n_distinct(tag_code_or_release_id))
) %>% View


# Looks to me like I will just want to disregard "run" because there is a lot of missing data in it.
# I also should take the basins that are missing and just lump them together within region as "missing"
# and same for missing regions.  I could try to model different years and things, but since i am mostly
# aggregating like this to get a hierarchical prior, I think it will be best to just average over all
# possible years that could be contributing.
