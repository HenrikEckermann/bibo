## ---- echo = FALSE, warning=FALSE, message=FALSE-------------------------
knitr::opts_chunk$set(echo = FALSE, warning=FALSE, message=FALSE)

## ---- fig.show = FALSE, results = "hide"---------------------------------
library(here)
# --------------------------- Data import from Leo
source(here("R/read.R"))

## ------------------------------------------------------------------------
# --------------------------- Check how many/which samples per day
library(knitr)
library(glue)
library(tidyverse)
# What information is there?
metadata <- sample_data(genus) %>% as.data.frame() %>% distinct(My.SQL.ID, .keep_all = T)
colnames(metadata)

## ------------------------------------------------------------------------
# select and rename relevant columns for checking
meta_sub <- 
  metadata %>% 
    mutate(subject_id = str_sub(unique.subjectID, 19,21)) %>%
    rename(age_days = time) %>%
    select(subject_id, age_days, sample)
    
# Are there subjects with more than one sample per time_point
multiple_samples <- 
  group_by(meta_sub, subject_id, age_days) %>% summarise(n = n()) %>%
  filter(n > 1)
kable(multiple_samples)

# I will have to ask why that is.
# For counting, I need to disregard all but one sample per time point.
# The selection is not based on any information and thus the resulting 
# dataframes should not be used for analyses without checking why we 
# have more than one sample and selecting accordingly. 

# show cases where we need to select from
# filter(meta_sub, subject_id %in% multiple_samples$subject_id) %>%
#   mutate(no = 1: 55)
# I write up all but one rownumber for each time point that has more than 1 n
v_disregard <- c(3, 6)
samples_disregard <- 
  filter(meta_sub, subject_id %in% multiple_samples$subject_id) %>%
  mutate(no = 1: 6) %>% 
  filter((no %in% v_disregard)) %>%
  select(sample)

## ------------------------------------------------------------------------
# Now we can count how many samples per time point
meta_sub %>% 
  filter(!(sample %in% samples_disregard$sample)) %>%
  group_by(age_days) %>% 
  summarise(n = n()) %>% kable()

## ------------------------------------------------------------------------
# To check how many subjects provided samples for certain time points, I
# can now create a new df so that we get NA where there is no sample 
# for a time point
meta_sub2 <- data_frame(
  subject_id = rep(unique(meta_sub$subject_id), each = 11),
  age_days = rep(sort(unique(meta_sub$age_days)), 186)) %>%
  left_join(meta_sub, by = c("subject_id", "age_days"))

# and now we can create the above mentioned df... 
sample_per_d <- 
  filter(meta_sub2, !(sample %in% samples_disregard$sample)) %>%
  mutate(age_days = glue("d{age_days}")) %>%
    spread(age_days, sample)
kable(sample_per_d)

## ---- echo = T-----------------------------------------------------------
# ... and select the time points of interest in the select function and omit 
# NA. Then by counting the rows we get the value of interest. E.g. for the time
# points 28 d, 75d, 105, 2193 we have...
sample_per_d %>%
  select(subject_id, d28, d75, d105, d2193) %>%
  na.omit() %>% dim()
# ... 68 complete cases. It may not be required to work only with complete cases
# but this might still be informative to decide what we can do.

## ------------------------------------------------------------------------
knitr::purl(here("R/sample_description.R"))

