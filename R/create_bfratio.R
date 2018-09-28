## ------------------------------------------------------------------------
# to view longer dfs I set print option
options(dplyr.print_max = 1e9)
library(microbiome)
library(tidyverse)
library(here)
# load HITChip data using Leo's script
source(here("R/read.R"))

## ------------------------------------------------------------------------
sample_data(genus) %>% as.data.frame() %>% distinct(My.SQL.ID, .keep_all = T)

## ------------------------------------------------------------------------
# load sample_data. disregard duplicate samples: %>% distinct(My.SQL.ID, .keep_all = TRUE)
# seems like this metadata is not complete. Therefore I use the excel from Gerben metadata <- sample_data(genus) 
# colnames(metadata)
metadata <- readxl::read_excel(here("data/meta_variables/my.metadata.xlsx"))
# to select the correct infants I use the resulting dataframe that arised from Gerbens selection process
select_samples <- read_csv(here("data/csv_gerben/complete.csv"))
# to confirm we have 49 pre/post for both groups
select_samples %>% group_by(childcarecenter, groupcode) %>% summarise(n = n())
metadata <- metadata %>%
    rename(subject_id = subject, sample = "Var.2") %>%
    filter(sample %in% select_samples$sample)

## ------------------------------------------------------------------------
# we have additional information stored in another excel sheet:
# I need to control for the age at CC entrace and the total length of CC (unless CC will extend)
# beyond post, which is just the indicator of the stool sample closest to post
# same for breastfeeding: I need the number or %BF until entrance and then during entrance,
# where % means breasfeedings/(breasfeedings + formula feedings), because those two variables will
# be highly correlated I assume, but I can check that...
# so the questions I need to aks are: have their been children who remained in CC after they gave their second sample?
# I assume that age at sample is in fact not the same among infants? (then 75, 105 etc. should be renamed really)
metadata2 <- readxl::read_excel(here("data/meta_variables/my.metadata2.xlsx"), sheet = "Data")
metadata2 <-
    metadata2 %>% 
        filter(ID %in% select_samples$subject) %>%
        rename(
            subject_id = ID, 
            cc = Childcare_yesno,
            per_bf_during_cc = `%BEFbetweenKDV1andKDV2`, 
            weeks_cc = `Tot#KDV`, 
            bf_until_end_cc = WeeksBEF,
            age_start = Begin_age_weeks) %>%
        select(subject_id, age_start, cc, weeks_cc, bf_until_end_cc, per_bf_during_cc)

## ------------------------------------------------------------------------
# how much variation is there in the # of week subjects were in CC?
metadata2 %>% summarise(mean_weeks_cc = mean(weeks_cc), sd_weeks_cc = sd(weeks_cc))

## ------------------------------------------------------------------------
# merge metadata, rename time (because it does not correspond to real age) 
meta <- metadata %>%
    left_join(metadata2, by = "subject_id") %>%
    mutate(time = ifelse(time == "75", "pre", "post")) %>%
    select(subject_id, time, cc, age_start, weeks_cc, bf_until_end_cc, per_bf_during_cc, sample)

## ------------------------------------------------------------------------
meta

## ------------------------------------------------------------------------
# Next we need to load additional metadata about the number of week breastfeeding and formula feeding
# as can be looked up in R/read_metadata_long.R, we get a df called meta_long
# ok I see a problem: above we have bf_before_cc, which is just +1 for each week an infant had bf, thereby
# disregarding the amount of bf. E.g. infant 204 had 3 weeks but we can see below 
# that it had on average 7.5 feedings a day for those 3 weeks, which is more informative
# I think about counting the number of bf until cc, use that as a predictor and then also count the bf during cc
# actually, I will then use bf/bf+formula. Then I avoid high negative corrlation of 2 predictors and still keep
# the info for when there is no correlation (e.g. very hungry infants that got additional formula)
source(here("R/read_metadata_long.R"))
meta_long %>% filter(subject == 204, week <=16)

## ------------------------------------------------------------------------
# After I had a look at the data above, I now try to replicate the bf_weeks and then change it to a count bf_before_cc
# and a count bf_weeks_cc
meta_long <-
    meta_long %>% 
        rename(subject_id = subject) %>%
        filter(subject_id %in% meta$subject_id) %>%
        left_join(meta, by = "subject_id")

## ------------------------------------------------------------------------
# intialize columns
meta$bf_count_pre <- NA
meta$bf_count_post <- NA
meta$bf_count_cc <- NA
meta$formula_count_pre <- NA
meta$formula_count_post <- NA
meta$formula_count_cc <- NA

# select only variables needed. I checked all the infants with NA for bf, formula and exp.bf
# most of the time there is NA for bf if e.g subject was formula fed or vice versa
# if not then these are only a few rows and given the consistency, most reasonably is to insert the value
# that was given before. The infant cannot starve + it received the same number of feedings before.
# only subjects 252, 448 have critically many NA. We would need to impute or drop those 2. They are both cc
#meta_long %>% 
#    select(subject_id, bf, expressed_bf, formula, cc) %>%
#    filter(subject_id %in% c(252, 448)) %>% print()

# select subset. Then for those we know NA will be very likely the that was there before. 
#bf_bool <- 
#    ifelse(
#        is.na(bf)&is.na(expressed_bf)&is.na(formula), mean(bf, na.rm = T),
#        ifelse(
#            is.na(bf)&(!is.na(expressed_bf)|!is.na(formula)), 0, bf))
# if there is NA but there is info in the other feedings columns, then NA = 0, 
# else most likely NA = mean unless subjects are 252 and 448 (why? because mothers throughout the 
# whole data are very consistent in feedings behavior)
bf_imputed <- 
    meta_long %>%
        select(subject_id, bf, expressed_bf, formula, week, age_start, weeks_cc) %>%
        filter(!subject_id %in% c(252, 448)) %>% 
        mutate(
            bf = ifelse(is.na(bf)&(!is.na(expressed_bf)|!is.na(formula)), 0,
                ifelse(is.na(bf)&is.na(expressed_bf)&is.na(formula), mean(bf, na.rm = T), bf))) %>%
        select(subject_id, week, age_start, weeks_cc, bf)
#
expressed_bf_imputed <- 
    meta_long %>%
        select(subject_id, bf, expressed_bf, formula) %>%
        filter(!subject_id %in% c(252, 448)) %>% 
        mutate(
            expressed_bf = ifelse(is.na(expressed_bf)&(!is.na(bf)|!is.na(formula)), 0,
                ifelse(is.na(bf)&is.na(expressed_bf)&is.na(formula), mean(expressed_bf, na.rm = T), expressed_bf))) %>%
        select(expressed_bf)

#
formula_imputed <- 
    meta_long %>%
        select(subject_id, bf, expressed_bf, formula) %>%
        filter(!subject_id %in% c(252, 448)) %>% 
        mutate(
            formula = ifelse(is.na(formula)&(!is.na(bf)|!is.na(expressed_bf)), 0, 
                ifelse(is.na(bf)&is.na(expressed_bf)&is.na(formula), mean(formula, na.rm = T), formula))) %>%
        select(formula)                               
feeding_imputed <- bind_cols(bf_imputed, expressed_bf_imputed, formula_imputed)
# add expressed bf and bf together because for our research question it does not matter
feeding_imputed <- feeding_imputed %>% mutate(bf = bf + expressed_bf)

## ------------------------------------------------------------------------
feeding_imputed

## ------------------------------------------------------------------------
# now for each subject and feedintype I need to aggregate the feedings
for (id in unique(feeding_imputed$subject_id)) {
    # pre CC bf
    bf_count_pre <- 
        feeding_imputed %>% filter(subject_id == id, week <= age_start) %>%
        select(bf) %>%
            sum()/2
    meta[meta$subject_id == id, "bf_count_pre"] = bf_count_pre
    # post CC bf
    bf_count_post <- 
        feeding_imputed %>% filter(subject_id == id, week <= age_start + weeks_cc) %>%
        select(bf) %>%
            sum()/2
    meta[meta$subject_id == id, "bf_count_post"] = bf_count_post
    # during CC bf
    meta[meta$subject_id == id, "bf_count_cc"] = bf_count_post - bf_count_pre

    # pre CC formula
    formula_count_pre <- 
        feeding_imputed %>% filter(subject_id == id, week <= age_start) %>%
        select(formula) %>%
            sum()/2
    meta[meta$subject_id == id, "formula_count_pre"] = formula_count_pre
    # post CC formula
    formula_count_post <- 
        feeding_imputed %>% filter(subject_id == id, week <= age_start + weeks_cc) %>%
        select(formula) %>%
            sum()/2
    meta[meta$subject_id == id, "formula_count_post"] = formula_count_post
    # during CC formula
    meta[meta$subject_id == id, "formula_count_cc"] = formula_count_post - formula_count_pre
}

## ------------------------------------------------------------------------
# those remaining NA are as mentioned completely NA for feedings variables OR there was not even
# an excel sheet for feeding.
meta %>% filter(is.na(bf_count_pre)|is.na(bf_count_post)) %>% distinct(subject_id)

## ------------------------------------------------------------------------
test <- meta %>% mutate(check_post = bf_count_post/(age_start + weeks_cc))
# they cannot correlate perfectly since I imputed NA where they likely dropped colons. I think
# it must be that the imputed values are closer to truth since the infants did not starve...
cor.test(test$bf_until_end_cc, test$check_post)
# lets see if formula and bf correlate as I suspected
cor.test(test$bf_count_pre, test$formula_count_pre)
cor.test(test$bf_count_post, test$formula_count_post)
cor.test(test$bf_count_cc, test$formula_count_cc)
# yes...And therefore it might be reasonable to create a ratio out of them. I loose the count then but
# I trust the ratio to be more correct than the count here since the count relies on the fact that mother count
# disciplines also minor deviances and I am not sure how realistic this is...I will see. bf_ratio will
# be either the bf_ratio for the time before CC for the pre sample or the ratio in the weeks during CC (post)
# pre and post correlate highly (.8) and therefore I cannot use separate predictors in regression. I currently
# think this is a reasonable solution. 
meta <- meta %>%
    mutate(
        bf_ratio_pre = bf_count_pre/(bf_count_pre + formula_count_pre),
        bf_ratio_cc = bf_count_cc/(bf_count_cc + formula_count_cc),
        bf_ratio = ifelse(time == "pre", bf_ratio_pre, bf_ratio_cc)
    )
meta   

## ------------------------------------------------------------------------
# I will now merge the data. I later will only need to slect genus for which abundance was not high enough in any sample
# But I don't know how to do this myself yet. Then the data is ready for model fitting. 
# Another question is to impute or to use lw deletion.
# I first use lw deletion to see if everything works as expected but if possible I want to impute
# because the information about cc/nocc and time is complete and that is most important!
data <- 
    otu_table(genus) %>% 
        as.data.frame() %>%
        rownames_to_column("genus") %>%
        mutate_if(is.numeric, log10) %>%
        gather(sample, value, -genus) %>%
        spread(genus, value) %>%
        filter(sample %in% meta$sample) %>%
        left_join(meta, by = "sample") %>%
        select(colnames(meta), everything()) %>% 
        select(- bf_until_end_cc, -per_bf_during_cc) # deselecting the old vars

data <- data %>% mutate(cc = ifelse(cc == 0, "no", "yes")) %>%
                  mutate(
                      cc = factor(cc, levels = c("no", "yes")), 
                      time = factor(time, levels = c("pre", "post")),
                      bf_count_pre_s = scale(bf_count_pre)[, 1],
                      bf_count_cc_s = scale(bf_count_cc)[, 1],
                      formula_count_pre_s = scale(formula_count_pre)[, 1],
                      formula_count_cc_s = scale(bf_count_cc)[, 1],
)

# rename problematic names 
colnames(data) <- gsub(" ", "_", colnames(data))
colnames(data) <- gsub("\\.", "", colnames(data))
colnames(data)[which(colnames(data) == "Clostridium_\\(sensu_stricto\\)")] <- "Clostridium_sensu_stricto"
