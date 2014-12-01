# a bunch of functions to make it easier to clean up, format, and
# properly make factors of the RMIS data bases.

# take a simple vector of integers that have consecutive values
# with known levels and turn it into the appropriate factor vector
int_to_factor <- function(x, L) {
  f <- as.factor(x)
  levels(f) <- L
  f
}


# return species corresponding to integer codes
species_levels <- function() {
  scan("data_defs/psc-4.1/species.txt", what = "character")
}

# return tag_status levels corresponding to integer codes which
# are preserved here as character names of the vector that gets returned
tag_status_levels <- function() {
  tmp <- read.table("data_defs/psc-4.1/recovery_tag_status.txt", row.names = 1, sep = "=", strip.white = TRUE,
                    quote = "'")
  ret <- gsub(" *\\(.*$", "", tmp$V2)
  names(ret) <- gsub("[’‘]", "", rownames(tmp))
  ret
}

# return a vector of mark types with their 4-digit identifiers
# as names
mark_coding_hash <- function() {
  tab <- read.csv("data_defs/psc-4.1/mark_coding.csv", colClasses = "character")
  # now make the 9nnn categories
  nines <- tab[grep("^0", tab$Code), ]
  nines$Mark <- gsub("^No Adclip", "Adipose Clip Unknown", nines$Mark)
  nines$Code <- as.character(9000 + as.numeric(nines$Code))
  tabn <- rbind(tab, nines)
  ret <- tabn$Mark
  names(ret) <- tabn$Code
  ret
}

# given a numeric vector of the tag codes, this makes it a factor
mark_to_string <- function(x) {
  hash <- mark_coding_hash()
  char_codes <- sprintf("%04d", x) # make 4 digit characters out of them
  char_codes[grepl("NA", char_codes)] <- NA # reconstitute the missing data
  ret <- rep(NA_character_, length = length(x))
  ret[!is.na(char_codes)] <- hash[char_codes[!is.na(char_codes)]]
  ret
}