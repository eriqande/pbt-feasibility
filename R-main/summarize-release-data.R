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


#### Some useful helper function #####

# return true if the mark code implies ad-clip and false otherwise (including if NA)
# this is used for counting up the release groups.  x is a vector of mark codes.
# They get read in as numeric, so they have to be padded back out to 4 digit strings
has_adclip <- function(x) {
  x[!is.na(x)] <- str_pad(x[!is.na(x)], 4, "left", "0")
  ifelse(!is.na(x) & str_sub(x, 1, 1) == 5,  TRUE, FALSE)
}


# this function returns "unknown", "yes", or "no" for the adclip status of recovered fish.
# It takes as input a vector of mark codes.
rec_adclip_status <- function(x) {
  x[!is.na(x)] <- str_pad(x[!is.na(x)], 4, "left", "0")
  x <- str_sub(x, 1, 1)  # get the first digit on it (this preserves NAs)
  
  # then get ready to switch the marks to our variables
  mycodes <- vector()
  mycodes["9"] <- "unknown"
  mycodes["5"] <- "yes"
  mycodes["0"] <- "no"
  
  ret <- unname(mycodes[x])
  
  ret[is.na(ret)] <- "unknown"  # this should never happen, but we include this in case it does
  
  ret
}

# this function takes a tag_status vector and condenses it into our categories
rec_tag_status <- function(x) {
  tagstats <- vector()
  tagstats[1] <- "cwt"
  tagstats[2] <- "no_tag"
  tagstats[c(3,4,7)] <- "no_read" 
  tagstats[8] <- "unknown"  # if head is not processed then it gets unknown
  tagstats[9] <- "awt"
  
  ret <- tagstats[x]
  
  ret[is.na(ret)] <- "unknown"
  
  ret
}


# this takes a vector of (supposedly) 19-character location codes and it parses them into a 
# data frame of hierarchically arranged location specifications.
parse_location_codes <- function(x) {
  tmp <- x %>%
    str_pad(width = 19, side = "right", pad = "-") %>%   # pad it out the right amount
    str_replace_all(" ", "*")  %>%  # replace internal spaces with *'s
    str_match("^([1-7])([MF])(.)(..)(....)(.{7})(...)$") %>%
    as.data.frame(stringsAsFactors = FALSE)
  
  names(tmp) <- c("full_loc_code",
                "state_or_province",
                "water_type",
                "sector",
                "region",
                "area",
                "location",
                "sub_location"
                )
  
  tmp
}



########### GETTING AND SUMMARIZING RELEASE DATA


#### get the entire releases data base  ####
releases <- readRDS("data/releases.rds")


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
  
  # here are the different hierarchical levels at which we will aggregate
  # NOTE: all the release_location_xxx fields are in the RMIS but are not part of the PSC specification
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


#### Then, for any particular species and years we can do this: ####

  
try_it <- group_releases_hierarchically()


############### GETTING AND SUMMARIZING RECOVERY DATA

# The general flow for this is as follows:
# 1. get the recovery data and condense it down to the tag_status and mark/beep status of each fish
#    and keep that associated with the catch_sample_id.
# 2. Use the location data base to get a hierarchical specification of all locations
# 2. Using the catch_sample_id and the location data base, give a hierarchy of locations for each recovered fish
# 3. Aggregate those cwt recoveries by counting them up grouped in different ways.
# 4. Add the catch-sample data so that we know how many untagged and un-beep fish were out there.

#### 1. Condensing and munging the recovery data  ####

# get the data (for now I just get the chinook and do one year of them.)
all_recov <- tbl_df(readRDS("data/chinook_recoveries.rds"))
tmp <- all_recov %>% 
  filter(run_year == 2012,   # filter it down to 2012
         str_sub(recovery_location_code, 2, 2) == "M" #  only take marine recoveries here
  )

# make three new columns:
# 1. "ad_clipped" gives the fish's ad-clip status (yes/no/unknown)
# 2. "beep" gives its beep status (yes/irrelevant).  We assume that every fish sampled in a 
#     fishery with electronic detection gets a positive beep status
# 3. "cwt_status" gives the status of the cwt recovered ("cwt", "no_tag", "no_read", "unknown", "awt")
# and then just pull out those columns along with the tag_code and the catch_sample_id and recovery location.
# All the other columns are just cruft that gets generated by referencing other tables.
rec <- tmp %>% 
  mutate(ad_clipped = rec_adclip_status(recorded_mark),
         beep = ifelse(detection_method == "E", "yes", "irrelevant"),
         cwt_status = rec_tag_status(tag_status)
         ) %>%
  select(tag_code, catch_sample_id, recovery_location_code, ad_clipped, cwt_status, beep)



#### 2. Split up location codes to hierarchical components and join them with the recovery data

# get the location codes and filter them to only those that have recoveries that we are focusing on
locs <- tbl_df(read.csv("data/locations.csv", stringsAsFactors = F, na.strings = "")) %>%
  filter(location_code %in% rec$recovery_location_code)

# then parse them out and keep just the columns that we want:
our_locs <- tbl_df(cbind(locs, parse_location_codes(locs$location_code))) %>%
  select(location_code, rmis_region, rmis_basin, full_loc_code:sub_location)

# and then join those to rec
rec_with_locs <- left_join(rec, our_locs, by = c("recovery_location_code" = "location_code"))

# now, just for fun and looking around, let's count these by different groupings:
rec_with_locs %>%
  group_by(state_or_province, sector, region, cwt_status, ad_clipped, beep) %>%
  tally() %>% View


