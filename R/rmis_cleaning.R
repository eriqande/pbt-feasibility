

# a bunch of functions to make it easier to clean up, format, and
# properly make factors of the RMIS data bases.

species_levels <- function() {
  c("Chinook", "Coho", "Steelhead", "Sockeye", "Chum", "Pink", "Masu", "Cutthroat", "Atlantic")
}

int_to_factor <- function(x, L) {
  f <- as.factor(x)
  levels(f) <- L
  f
}