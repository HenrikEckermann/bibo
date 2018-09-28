# In this script I change variable types and prepare the meta_long.csv to work
# in R. To see how those variables have been preprocessed, check the python jup
# yter notebook extract_meta.

# import metadata_long
library(dplyr)
library(here)
# change filepath to the metadata_long.csv
file <- here("data/meta_data/metadata_long.csv")
meta_long <- read.csv(file)
# convert to formats to work with
meta_long <- 
  select(meta_long,  -X,) %>%
  select(subject, week, bf, expressed_bf, formula, solid_food, type_solid_food, pet) %>%
  arrange(subject) %>%
  mutate(
    bf = as.numeric(bf),
    expressed_bf = as.numeric(expressed_bf),
    formula = as.numeric(formula),
    solid_food = as.numeric(solid_food),
    type_solid_food = ifelse(
      as.character(type_solid_food) == "", NA, 
      as.character(type_solid_food)
    )
)



