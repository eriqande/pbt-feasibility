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


#### Grab all the releases from 2004 to 2012 and order the release groups by state and region ####


# grab all these and then we will filter down later only to the ones that we've
# seen recovered in 2012
releases2004_2012 <- group_releases_hierarchically(BroodYears = 2004:2012)


# While we are at it, we want to make a data frame we can use to order things
# in some appropriate manner.  I could be more refined about ordering the regions
# within state, but i will worry about that if I have more time later
tmp <- releases2004_2012$Release_ID %>%  # make a frame of all the regions
  group_by(release_location_state, release_location_rmis_region) %>% 
  tally()


tmp$release_location_state <- factor(tmp$release_location_state, 
                                     levels = c("AK", "BC", "WA", "ID", "OR", "CA"))

release_order_info <- tmp %>%
  arrange(release_location_state) %>%
  ungroup %>%
  mutate(order_em = 1:nrow(tmp)) %>%
  inner_join(releases2004_2012$Release_ID, .)





############### GETTING AND SUMMARIZING RECOVERY DATA  ####

# The general flow for this is as follows:
# 1. get the recovery data and condense it down to the tag_status and mark/beep status of each fish
#    and keep that associated with the catch_sample_id.
# 2. Use the location data base to get a hierarchical specification of all locations

#### 1. Condensing and munging the recovery data. Add colums "ad_clipped", "beep", and "cwt_status"  ####

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


#### 2. Split up location codes into hierarchical components and join them with the recovery data  ####

# get the location codes and filter them to only those that have recoveries that we are focusing on
locs <- tbl_df(read.csv("data/locations.csv", stringsAsFactors = F, na.strings = "")) %>%
  filter(location_code %in% rec$recovery_location_code)

# then parse them out and keep just the columns that we want:
our_locs <- tbl_df(cbind(locs, parse_location_codes(locs$location_code))) %>%
  select(location_code, rmis_region, rmis_basin, full_loc_code:sub_location)

our_locs <- our_locs %>% distinct(location_code)  # for some reason there are duplicated location_codes!  We toss them out so there is just one of each here.


# and then join those to rec
rec_with_locs <- left_join(rec, our_locs, by = c("recovery_location_code" = "location_code"))

# NOTE: rmis_basin in the recoveries doesn't seem to be the basins --- more like a sub-region
# Though for freshwater recovery it probably means more.  Probably it is going to be best for us to
# use state_or_province, sector, region, area, location for agglomerating these.

# now, just for fun and looking around, let's count these by different groupings:
rec_with_locs %>%
  group_by(state_or_province, sector, region, cwt_status, ad_clipped, beep) %>%
  tally()


#### 2.1 Make some crazy "tree" plots that will guide us in aggregating recoveries ####
justcwts <- rec_with_locs %>%
  filter(cwt_status == "cwt")

rec2 <- justcwts %>%
  group_by(state_or_province) %>%
  summarise(state_n = n()) %>% 
  inner_join(justcwts)

rec3 <- rec2 %>%
  group_by(state_or_province, sector) %>%
  summarise(sector_n = n()) %>%
  inner_join(rec2)

rec4 <- rec3 %>%
  group_by(state_or_province, sector, region) %>%
  summarise(region_n = n()) %>%
  inner_join(rec3)

rec5 <- rec4 %>%
  group_by(state_or_province, sector, region, area) %>%
  summarise(area_n = n()) %>%
  inner_join(rec4)

rec6 <- rec5 %>%
  group_by(state_or_province, sector, region, area, location) %>%
  summarise(location_n = n()) %>%
  inner_join(rec5)


rec7 <- rec6 %>%
  select(state_or_province, state_n, sector, sector_n, region, region_n, area, area_n, location, location_n) %>%
  unique %>%
  ungroup %>%
  arrange(state_or_province, desc(sector_n), sector, 
          desc(region_n), region, desc(area_n), area, desc(location_n), location)

rec7$location_x <- 1:nrow(rec7)

rec8 <- rec7 %>%
  group_by(state_or_province, sector, region, area) %>%
  summarise(area_x = mean(location_x)) %>%
  left_join(rec7)

rec9 <- rec8 %>%
  group_by(state_or_province, sector, region) %>%
  summarise(region_x = mean(area_x)) %>%
  left_join(rec8)

rec10 <- rec9 %>%
  group_by(state_or_province, sector) %>%
  summarise(sector_x = mean(region_x)) %>%
  left_join(rec9)

rec11 <- rec10 %>%
  group_by(state_or_province) %>%
  summarise(state_x = mean(sector_x)) %>%
  left_join(rec10)

rec12 <- rec11 %>% 
  arrange(state_or_province, desc(sector_n), sector, 
          desc(region_n), region, desc(area_n), area, desc(location_n), location)

# and now we can ggplot that
number_ticks <- function(n) {function(limits) pretty(limits, n)}

rec12_4_plot <- rec12  # just make a copy cause i am going to tweak it
states_tmp <- c("AK", "BC/Yukon", "WA", "ID", "OR", "CA", "High Seas")  # fiddle to get it to plot North to South
rec12_4_plot$state_or_province <- factor(states_tmp[as.numeric(rec12_4_plot$state_or_province)], levels = states_tmp)
rec12_4_plot
g <- ggplot(data = rec12_4_plot) +
  geom_point(aes(x = state_x, y = state_n), color = "red", size = 0.2) + 
  geom_segment(aes(x = state_x, xend = sector_x, y = state_n, yend = sector_n), color = "red", size = 0.2) +
  geom_point(aes(x = sector_x, y = sector_n), color = "orange", size = 0.2) + 
  geom_segment(aes(x = sector_x, xend = region_x, y = sector_n, yend = region_n), color = "orange", size = 0.2) +
  geom_point(aes(x = region_x, y = region_n), color = "blue", size = 0.2) + 
  geom_segment(aes(x = region_x, xend = area_x, y = region_n, yend = area_n), color = "blue", size = 0.2) +
  geom_point(aes(x = area_x, y = area_n), color = "violet", size = 0.2) + 
  geom_segment(aes(x = area_x, xend = location_x, y = area_n, yend = location_n), color = "violet", size = 0.2) +
  scale_x_continuous(breaks=number_ticks(20)) +
  facet_wrap(~state_or_province, ncol = 1, scales = "free")
  
  
  
ggsave(file = "recovery_trees.pdf", width = 18.5, height = 20)


#### 2.2 Based on the recovery_trees.pdf plot, aggregate recovery locations, and name them, etc. ####
# This isn't super reproducible, but I am just going to focus on run_year = 2012 here
# and do it by hand...
fishery_breaks <- c(0, 22.1, 35.1, 78.1, 116.1, 235.1,  # 5 groups in alaska
                    500.1, 604.1, # 2 in BC
                    606.1, 610.1, 612.1, 745.1, 756.1, # 4 in WA, plus one small one
                    774.1, # one in oregon.
                    775.1, 776.1, 777.1, 778.1, 791.1)
rec13 <- rec12 %>%
  mutate(recovery_group = as.integer(cut(rec12$location_x, breaks = fishery_breaks))) %>%
  filter(!is.na(recovery_group))   # this tosses the high-seas recoveries

# now plot it with separators for the different recovery groups
frame_for_vline <- data.frame(state_or_province = rep(states_tmp[-c(4,7)], times = c(6, 2, 5, 1, 5)), fishery_breaks)
g + geom_vline(data = frame_for_vline[-1, ], mapping = aes(xintercept = fishery_breaks + 0.5))
ggsave(file = "recovery_trees_divided.pdf", width = 14, height = 10)



# now, let's also summarise to see what these fisheries are:
rec13 %>%
  group_by(state_or_province, region, recovery_group) %>%
  tally() %>%
  ungroup %>%
  arrange(recovery_group)

# and from that we see that we could name our recovery groups like this:
recov_group_names <- c("01-AK-NW", "02-AK-NW", "03-AK-SE", "04-AK-NE", "05-AK-Various",
  "06-BC", "07-BC",
  "08-WA", "09-WA", "10-WA", "11-WA", "12-WA",
  "13-OR",
  "14-CA", "15-CA", "16-CA", "17-CA", "18-CA")

rec13$recovery_group <- recov_group_names[rec13$recovery_group]

# let's count these up
tmp <- rec13 %>%
  group_by(recovery_group)  %>%
  summarise(num_cwts = sum(location_n)) %>%
  as.data.frame

# in fact, we could add these to the names of the recovery groups
new_names <- setNames(paste(tmp[,1], "  (", tmp[,2], ")", sep = ""), tmp[,1])

rec13$recovery_group <- new_names[rec13$recovery_group]


# now, just get the values columns we want
rec14 <- rec13 %>%
  select(state_or_province, sector, region, area, location, recovery_group)


# and finally, join this to rec_with_locs so that we have the recovery group for all the recoveries
rec2012_heavy <- inner_join(rec_with_locs, rec14)


# and here we just pick out the columns we are going to want to keep
rec2012 <- rec2012_heavy %>%
  select(tag_code:beep, recovery_group)


#### 2.3 Get the "distinct_codes" and order them  ####
# Now, filter it to just cwt recoveries, and store just the distinct codes,
# and on each of those reattach the order_me and the state, and basin, info,
# and give each a unique number that counts as their x-position
distinct_codes <- rec2012 %>% 
  filter(cwt_status == "cwt") %>%
  select(tag_code) %>%
  distinct %>%
  inner_join(release_order_info,  by = c("tag_code" = "tag_code_or_release_id")) %>%
  arrange(order_em, release_location_rmis_basin, brood_year) %>%
  mutate(tag_code_x_value = 1:nrow(.))


 


#### 2.4 Make some useful and important plots of the recoveries by tag group  ####

# then count up how many recoveries of each tag code in each recovery group
# and join the info in distinct codes to that
rec2012a <- rec2012 %>%
  group_by(recovery_group, tag_code) %>%
  summarise(num_tag_recovs_in_group = n()) %>%
  inner_join(distinct_codes)


# to make a rug we will want to have the x-values where state release groups start
state_starts <- distinct_codes %>% group_by(release_location_state) %>% summarise(xval = min(tag_code_x_value))
# and the same for where regions start too:
region_starts <- distinct_codes %>% group_by(release_location_state, release_location_rmis_region) %>% summarise(xval = min(tag_code_x_value))

# make a function to ggplot things:
recovery_histo_plot <- function(df, colorby =  "release_location_state") {
  ggplot(data = df) +
    geom_rug(data = region_starts, aes_string(x = "xval"), size = 0.2, colour = "grey50") +
    geom_vline(xintercept = state_starts$xval, size = 0.04, colour = "black") +
    facet_wrap(as.formula( "~ recovery_group"), ncol = 1, scales = "free_y") + 
    theme_bw() + theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
    geom_segment(aes_string(x = "tag_code_x_value", xend = "tag_code_x_value", y = 0, yend = "num_tag_recovs_in_group", colour = colorby),
                 size = 0.1) +
    scale_colour_discrete(drop = FALSE) +
    xlab("Tag code") + 
    ylab("Number of CWTs recovered")
}  

# then make the plot
recovery_histo_plot(rec2012a)

ggsave(file = "recovery_histo_panels.pdf", width = 8.5, height = 18)


# then make it and color it by factor(brood_year)
tmp <- rec2012a
tmp$brood_year <- factor(tmp$brood_year)
recovery_histo_plot(tmp, "brood_year")

ggsave(file = "recovery_histo_panels_brood_year.pdf", width = 8.5, height = 18)


# finally, break it down into two plots and get the states ordered correctly
rectmp <- rec2012a
rectmp$release_location_state <- factor(rectmp$release_location_state, levels = c("AK", "BC", "WA", "ID", "OR", "CA"))
tmp <- rectmp %>% filter(recovery_group < "10")
recovery_histo_plot(tmp)
ggsave(file = "recovery_histo_panel_1.pdf", width = 14, height = 10)

tmp <- rectmp %>% filter(recovery_group >= "10")
recovery_histo_plot(tmp)
ggsave(file = "recovery_histo_panel_2.pdf", width = 14, height = 10)



#### 2.5 Modify recovery data to correspond to the notation in the .tex document  ####

# cwt_status must be "cwt" even if it was not read, but the tag_code has to be unknown.
rec2012$tag_code[rec2012$cwt_status == "no_read"] <- "unknown"
rec2012$tag_code[rec2012$cwt_status == "unknown"] <- "pending"
rec2012$tag_code[rec2012$cwt_status %in% c("awt", "no_tag")] <- "irrelevant"


###################### FURTHER SUMMARIZE RELEASE DATA TO GET THE f's and p's ################


#### 1. Summarize the release data into the f and p values that we need for each tag code ####
mark_and_tag_rates <- distinct_codes %>%
  mutate(n_total_fish = n_tag_ad + n_tag_noad + n_notag_ad + n_notag_noad,
         n_marked = n_tag_ad + n_notag_ad,
         n_unmarked = n_tag_noad + n_notag_noad,
         f_marked = n_marked / n_total_fish,
         f_unmarked = n_unmarked / n_total_fish,
         p_marked = ifelse(n_marked > 0, n_tag_ad / n_marked, 0),
         p_unmarked = ifelse(n_unmarked > 0, n_tag_noad / n_unmarked, 0)
         )  %>% 
  ungroup() %>%
  mutate(release_location_state = factor(release_location_state, levels = c("AK", "BC", "WA", "ID", "OR", "CA"))) %>%
  select(tag_code, release_location_state, tag_code_x_value:p_unmarked)


#### 2. Make some plots to investigate sizes of release groups ####
ggplot(mark_and_tag_rates, aes(x = n_total_fish)) + geom_histogram()

# so there are some very large 
# now look at the proportion marked fish carrying cwts broken down by marking rate
ggplot(mark_and_tag_rates, aes(x = n_total_fish, y = f_marked)) + geom_point() + scale_x_log10()

tmp <- mark_and_tag_rates
tmp$mark_fract <- cut(mark_and_tag_rates$f_marked, c(-0.001, seq(0,1,by=0.2)))
ggplot(tmp, aes(x = n_total_fish, 
                               y = p_marked, 
                               colour = mark_fract
                               )) + 
  geom_point(alpha = 0.3, size = 1.5) + 
  scale_x_log10() +
  facet_wrap( ~ mark_fract)

ggsave("ppn_marked_that_have_tags_panels.pdf", width = 14, height = 10)


#### SUMMARIZE THE CATCH/SAMPLE DATA TO PREP FOR MCMC ####

#### 1. Read in the catch-sample data ####

# first, get the catch-sample data just for the catch-sample IDs that were involved 
cs1 <- readRDS("data/catch_sample.rds")


#### 2. Deal with catch-samples that occur in more than one recovery group ####

# pick out only those catch-samples that are relevant, but there is a problem.  It appears  some catch_sample_ids
# occur in more than one recovery group.  Blast!
# here are the catch_sample_ids
dupie_cs <- rec2012 %>%
  group_by(catch_sample_id, recovery_group) %>%
  distinct %>% 
  group_by(catch_sample_id) %>% 
  tally %>%
  filter(n>1, !is.na(catch_sample_id)) %>%
  select(catch_sample_id)

# now we can see which recovery groups they are in:
lame <- rec2012 %>% ungroup %>%
  group_by(catch_sample_id, recovery_group) %>%
  tally() %>%
  ungroup %>%
  filter(catch_sample_id %in% dupie_cs$catch_sample_id)
  
# and if we want to see just the recovery groups involved we can do:
lame %>% ungroup() %>%
  group_by(recovery_group) %>%
  distinct

# So, it appears that the way forward on this is to lump 14-CA, 15-CA, 17-CA, and 18-CA
# and also to lump 01-AK-NW and 02-AK-NW.
# I'll just go ahead and do that:
rec2012$recovery_group[str_detect(rec2012$recovery_group, "^1[4578]-CA") ] <-
  paste(unique(rec2012$recovery_group[str_detect(rec2012$recovery_group, "^1[4578]-CA") ]), collapse = ", ")

rec2012$recovery_group[str_detect(rec2012$recovery_group, "^0[12]-AK") ] <-
  paste(unique(rec2012$recovery_group[str_detect(rec2012$recovery_group, "^0[12]-AK") ]), collapse = ", ")

# also toss out those missing a catch_sample_id:  we dump about 6000 here.
rec2012 <- rec2012 %>% 
  filter(!is.na(catch_sample_id))

# So, that was lame, but necessary to keep from slowing us down here....
# check the number of distinct catch_sample_ids
rec2012 %>% 
  group_by(catch_sample_id, recovery_group) %>%
  tally()


#### 3. Attach recovery_groups to every catch_sample_id  ####
# continuing, we attach recovery_groups to every catch_sample_id, 
cs2 <- rec2012 %>%  
  group_by(catch_sample_id, recovery_group) %>% 
  summarise(tot_fish_in_recov = n()) %>%
  inner_join(cs1 %>% filter(catch_year == 2012))  # catch_sample_ids are not all unique across years...


#### 4. Explore the results and worry about some inconsistencies, especially with the BC data ####

# now we just look at some numbers sampled in different scenarios
cs2 %>%
  group_by(recovery_group, detection_method) %>%
  summarise_each(funs(sum(., na.rm = TRUE)), vars = mr_1st_partition_size:mr_2nd_sample_obs_adclips)

# see if we find NA's where maybe they oughtn't be
cs2 %>%
  group_by(recovery_group, detection_method) %>%
  summarise_each(funs(sum(.)), vars = mr_1st_partition_size:mr_2nd_sample_obs_adclips)
# this data base is sloppy.

# now, find rows that have NA's in "mr_1st_sample_known_ad_status" "mr_1st_sample_obs_adclips" 
# 'cuz the damn Canadians have some in there. Is that a problem?
bc_bozos <- cs2 %>%
  filter(detection_method == "V" & (is.na(mr_1st_sample_known_ad_status) | is.na(mr_1st_sample_obs_adclips) ))
# Well, those are just catches that they didn't sample for Marks.  OK....Whatever. The num_sampled is NA in it
# so it shouldn't throw things off.

# now, what about electronic WA fisheries that don't have "mr_2nd_sample_known_ad_status" "mr_2nd_sample_obs_adclips"
cs2 %>%
  filter(detection_method == "E" & (is.na(mr_2nd_sample_known_ad_status) | is.na(mr_2nd_sample_obs_adclips) ))
# that turns out to be just 2 rows, with almost no fish sampled.  Not a problem, really.


#### 5. Classify mark-sampling as beep-dependent or beep-independent in electronically detected fisheries ####

tmp <- as.data.frame(cs2)
tmp2 <- tmp[, c("mr_1st_partition_size", "mr_1st_sample_size", "mr_1st_sample_known_ad_status", "mr_1st_sample_obs_adclips", "mr_2nd_partition_size",
       "mr_2nd_sample_size", "mr_2nd_sample_known_ad_status", "mr_2nd_sample_obs_adclips")]
tmp2[is.na(tmp2)] <- 0
tmp3 <- as.matrix(tmp2)
first <- rowSums(tmp3[, 1:4])
second <- rowSums(tmp3[, 5:8])
beep_clip_sampling <- rep(NA, nrow(cs2))
beep_clip_sampling[cs2$detection_method == "E" & first > 0 & second == 0] <- "independent"
beep_clip_sampling[cs2$detection_method == "E" & first > 0 & second > 0] <- "dependent"


# investigate those
cs2$electronic_mark_sampling_type <- beep_clip_sampling
cs2 %>%
  group_by(recovery_group, detection_method, electronic_mark_sampling_type) %>%
  summarise_each(funs(sum(., na.rm = TRUE)), vars = mr_1st_partition_size:mr_2nd_sample_obs_adclips)

# from this, it is clear that the mark sampling is essentially all dependent (except for
# three fish which were meant to be dependent, no doubt)  So, we will treat them as such
# and just mix them all together.

# so we are back here:
cs2 %>%
  group_by(recovery_group, detection_method) %>%
  summarise_each(funs(sum(., na.rm = TRUE)), vars = mr_1st_partition_size:mr_2nd_sample_obs_adclips)
# So, the BC-06 sample is the only one with both electronic and visual sampling.


#### 6. Check that the total numbers of fish in the recovery file and the catch/sample information is consistent ####

tmp1 <- rec2012 %>% 
  filter(ad_clipped == "yes") %>%
  group_by(recovery_group) %>%
  summarise(adclipped_fish_in_recovery_data_base = n())

tmp2 <- cs2 %>%
  group_by(recovery_group) %>%
  summarise(adclipped_fish_reported_in_catch_sample = sum(mr_1st_sample_obs_adclips, na.rm = TRUE))
 
tmp3 <- rec2012 %>% 
  filter(ad_clipped == "yes" | ad_clipped == "unknown") %>%
  group_by(recovery_group) %>%
  summarise(adclipped_or_unknown_fish_in_recovery_data_base = n())


compare_totals <- inner_join(tmp1, tmp2)
ggplot(compare_totals, aes(y = adclipped_fish_in_recovery_data_base, x = adclipped_fish_reported_in_catch_sample)) +
  geom_point()  + 
  geom_abline(intercept = 0, slope = 1)

# which shows that we are pretty much right on except BC and OR.  OK....
# 07-BC really doesn't work!  All the others are OK, but 12-WA and 07-BC are too small anyway.

# and we can see what the picture looks like when we include the "unknown" ad-clipped fish
compare_totals <- inner_join(tmp3, tmp2)
ggplot(compare_totals, aes(y = adclipped_or_unknown_fish_in_recovery_data_base, x = adclipped_fish_reported_in_catch_sample)) +
  geom_point()  + 
  geom_abline(intercept = 0, slope = 1)
# that helps for some but not for others.


# So, we can continue!  

# let's just enumerate fish with different qualities/features/properties
rec2012 %>%
  group_by(recovery_group, cwt_status, ad_clipped, beep) %>%
  tally() %>%
  as.data.frame


#### Write out data to send off to 02-do-mcmc.R ####

# write it to a serialized data object list with  named components 
saveRDS(list(recovery = rec2012, 
             catch_sample = cs2, 
             mark_and_tag_rate = mark_and_tag_rates),
        compress = "xz",
        file = "data_for_mcmc.rds")



#### Explore the occurrence of unassociated and unassociated releaseses in CA #####
# this is just some fiddling I did after lab meeting.

rel_explore <- releases  %>%
  mutate(n_total_fish = n_tag_ad + n_tag_noad + n_notag_ad + n_notag_noad,
         n_marked = n_tag_ad + n_notag_ad,
         n_unmarked = n_tag_noad + n_notag_noad,
         f_marked = n_marked / n_total_fish,
         f_unmarked = n_unmarked / n_total_fish,
         p_marked = ifelse(n_marked > 0, n_tag_ad / n_marked, 0),
         p_unmarked = ifelse(n_unmarked > 0, n_tag_noad / n_unmarked, 0),
         f_marked_category = cut(f_marked, breaks = c(-.001, seq(0,1, by = .1)))
  ) %>% select(tag_code_or_release_id, brood_year, release_location_state, tag_variety, f_marked_category, n_total_fish, f_marked, p_marked)


# If I look at all the releases from Calfornia that have at least one individual recovered coastwide in 2012
# and count up the total number of fish released
mark_and_tag_rates %>%
  filter(release_location_state == "CA") %>% 
#  group_by(cut(f_marked, breaks = seq(0,1, by = 0.1))) %>% 
  summarise(tot_marked_released = sum(n_marked))

rel_explore %>%
  filter(release_location_state == "CA", brood_year %in% 2008:2010, tag_variety == "cwt") %>%
#  group_by(release_location_state, brood_year, tag_variety, f_marked_category) %>%
  summarise(total_marked_released = sum(n_total_fish * f_marked))
  

releases %>%
  filter(release_location_state == "CA") %>%
  group_by(release_location_rmis_basin) %>% 
  tally()

