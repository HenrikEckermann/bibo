## ----setup, include=FALSE------------------------------------------------
library(knitr)
opts_chunk$set(echo = F, warning=F, message=F, results = "hide", fig.show = "hide", fig.pos = "H")
#for kableExtra:
#options(knitr.table.format = "html")

## ----echo = FALSE--------------------------------------------------------
# import libraries
library(readxl)
library(foreign)
library(tidyverse)
library(glue)
library(kableExtra)
library(here)
library(papaja)

# load data to have similar subjects 
load(here("data/cc_analyses_workspace.RData"))
# subjects ids in data$subject_id
## ------------------------------------------------------------------------
# BIBO confounders dataset 
confounders <- read.spss(here("data/meta_data/bibo_confounders.sav"), to.data.frame = T)
# select relevant columns
confounders <- select(confounders, -group193)
# rename for convenience
colnames(confounders) <- c("ID", "birthweight", "siblings", "pre_smoking",
 "pre_alcohol", "apgar5min", "maternal_age", "delivery", "gestat_length",
  "childsex", "maternal_education", "firstborn")
# Carolina sent me also this file because the confounders file is incomplete for usage:
june2017_gerben <- read_excel(here("data/meta_data/my.metadata2.xlsx"), sheet = "Data")
# rename colnames for convenience
colnames(june2017_gerben) <- c("ID", "begin_group", "group", "group_strict",
 "group_BF_noBF", "CC_noCC", "weeks_BF", "percent_BF_kdv1_kdv2", "CC",
 "age_CC_min2d", "age_CC_plus28d",  "age_begin_weeks", "age_plus_4weeks",
 "tot_kdv", "kdv1", "kdv2")
# merge confounders with df so that I can compare groups:
df_complete <- june2017_gerben %>% 
  select(-tot_kdv, - kdv1, -kdv2) %>%
  left_join(confounders, by = "ID")

## ------------------------------------------------------------------------
# # select final sample 
# df <- filter(df_complete, 
#   delivery == "natuurlijke bevalling" | delivery == "pomp", apgar5min > 7,
#   begin_group == 1, group_strict != 99999
# )

# I use Gerbens selection:
df_complete <- df_complete %>% filter(ID %in% data$subject_id)
df_complete %>% group_by(CC_noCC) %>% summarise(n = n())
# change numericals to factors
df_complete <- df_complete %>%
  mutate_at(
    vars(c("group_strict", "group_BF_noBF", "CC_noCC", "CC", "childsex", "firstborn")),
    funs(factor)
  )
data_pre <- data %>% 
  filter(time == "pre") %>% 
  select(subject_id, age_d, bf_count) %>% 
  rename(age_pre = age_d, bf_count_pre = bf_count, ID = subject_id)
data_post <- data %>% 
  filter(time == "post") %>% 
  select(subject_id, age_d, bf_count, sibling, csection, cc) %>% 
  rename(age_post = age_d, bf_count_post = bf_count, ID = subject_id)

data_add <- left_join(data_pre, data_post, by = "ID") %>%
  mutate(time_between_samples = age_post - age_pre)
data_add$ID <- as.numeric(data_add$ID)
df_complete <- df_complete %>% left_join(data_add, by = "ID") %>% rename(subject_id = "ID")

colnames(df_complete)

df_complete %>% group_by(cc) %>% summarise(sum(sibling, na.rm = T))


df_complete %>% group_by(CC) %>% summarise(mean = mean(time_between_samples), sd =  sd(time_between_samples))

# make df for table
df_table <- df_complete
colnames(df_table) <- c('ID', 'begin_group', 'group', 'Group',
 'group_BF_noBF', 'CC_noCC', 'weeks_BF', 'percent_BF_kdv1_kdv2', 'CC',
 'age_CC_min2d', 'age_CC_plus28d', 'age_begin_weeks', 'age_plus_4weeks',
 'birthweight', 'siblings', 'pre_smoking', 'pre_alcohol', 'apgar5min',
 'maternal_age', 'delivery', 'gestat_length', 'childsex',
 'maternal_education', 'firstborn','age_pre', 'breastfeeding_pre', 'age_post', 'breastfeeding_post', "sibling", "csection", "cc", 'time_between_samples')
head(df_table)

colnames(df_table)

df_table <-df_table %>% mutate(Childcare = ifelse(CC == "0", "HOME (n = 49)", "CC (n = 49)"))


# create df with descriptives 
d <- group_by(df_table, Childcare) %>%
  summarise(
    male = sum(childsex == 1),
    female = sum(childsex == 0),
    mean_age_pre = mean(age_pre, na.rm = T),
    sd_age_pre = sd(age_pre, na.rm = T),
    `mean (sd)_age_pre` = NA,
    min_age_pre = min(age_pre),
    max_age_pre = max(age_pre),
    mean_age_post = mean(age_post, na.rm = T),
    sd_age_post = sd(age_post, na.rm = T),
    `mean (sd)_age_post` = NA,
    min_age_post = min(age_post),
    max_age_post = max(age_post),
    mean_maternal_age = mean(maternal_age, na.rm = T),
    sd_maternal_age = sd(maternal_age, na.rm = T),
    `mean (sd)_maternal_age` = NA,
    min_maternal_age = min(maternal_age, na.rm = T), 
    max_maternal_age = max(maternal_age, na.rm = T), 
    
    mean_birthweight = mean(birthweight, na.rm = T),
    sd_birthweight = sd(birthweight, na.rm = T),
    `mean (sd)_birthweight` = NA,
    min_birthweight = min(birthweight, na.rm = T), 
    max_birthweight  = max(birthweight, na.rm = T),
    mean_breastfeeding_pre = mean(breastfeeding_pre, na.rm = T),
    sd_breastfeeding_pre = sd(breastfeeding_pre, na.rm = T),
    `mean (sd)_breastfeeding_pre` = NA,
    min_breastfeeding_pre = min(breastfeeding_pre, na.rm = T),
    max_breastfeeding_pre = max(breastfeeding_pre, na.rm = T),
    mean_breastfeeding_post = mean(breastfeeding_post, na.rm = T),
    sd_breastfeeding_post = sd(breastfeeding_post, na.rm = T),
    `mean (sd)_breastfeeding_post` = NA,
    min_breastfeeding_post = min(breastfeeding_post, na.rm = T),
    max_breastfeeding_post = max(breastfeeding_post, na.rm = T),
    yes = sum(sibling, na.rm = T),
    no =  sum(sibling == 0, na.rm = T),
    yes1 = sum(csection, na.rm = T),
    no1 = sum(csection == 0, na.rm = T)
  ) %>%
    mutate_if(is.numeric, funs(format(round(., 1), nsmall = 0))) %>%
    mutate(
      `mean (sd)_maternal_age` = glue("{mean_maternal_age} $\\pm$ {sd_maternal_age}"),
      `mean (sd)_birthweight` = glue("{mean_birthweight} $\\pm$ {sd_birthweight}"),
      `mean (sd)_breastfeeding_pre` = glue("{mean_breastfeeding_pre} $\\pm$ {sd_breastfeeding_pre}"),
      `mean (sd)_breastfeeding_post` = glue("{mean_breastfeeding_post} $\\pm$ {sd_breastfeeding_post}"),
      `mean (sd)_age_pre` = glue("{mean_age_pre} $\\pm$ {sd_age_pre}"),
      `mean (sd)_age_post` = glue("{mean_age_post} $\\pm$ {sd_age_post}"),
    ) %>%
    select(-mean_maternal_age, -sd_maternal_age, -mean_birthweight, -sd_birthweight, -mean_breastfeeding_pre, -sd_breastfeeding_pre,-mean_breastfeeding_post, -sd_breastfeeding_post, -mean_age_pre, -sd_age_pre, -mean_age_post, -sd_age_post)
d
sum(df_complete$sibling == 1, na.rm = T)
data$sibling
# obtain order of colnames since spread does not preserve it 
# but I need to keep the order for the rows laters
ordered_rows <- colnames(d)
d
# back to wide format, reorder, names ready for table, merge sd + mean to str

d <- 
  mutate_if(d, is.numeric, funs(toString(.))) %>%  
    gather(variable, value, -Childcare) %>%
    spread(Childcare, value) %>%
    arrange(match(variable, ordered_rows[-1])) %>%
    mutate(variable = gsub("_.+", "", variable))

#t.test(age_pre ~ CC, data = df_complete, var.equal = F)
df_complete %>% select(age_begin_weeks, age_pre) %>% mutate(age_pre = age_pre/7)
# chisquare
gmodels::CrossTable(df_complete$CC, df_complete$csection, chisq=T, fisher=T, expected=T, sresid=T, prop.c = FALSE, prop.r = FALSE, prop.t = FALSE, prop.chisq = FALSE, format='SPSS')
d %>% class()
d <- d %>% mutate(variable = ifelse(variable == "yes1", "yes", ifelse(variable == "no1", "no", variable)))
colnames(d)[1] <- ""
df_table
df_table %>% group_by(cc, csection) %>% summarise(n())
  
table_caption_descr <- "Descriptive statistics for demographic variables of infants and mothers included in the present study."
table_note_descr <- "CC = childcare. Breastfeeding refers to the average number of breast-feedings per day. Age refers to the age in days."
d

list_of_demographics <- list("**Gender**" = d[1:2,], "**Age PRE**" = d[3:5, ], "**Age POST**" = d[6:8, ], "**Maternal Education**" = d[9:11, ], "**Birthweight**" = d[12:14,], "**Breastfeeding (Birth - PRE)**" = d[15:17,], "**Breastfeeding (PRE - POST)**" = d[18:20,], "**Sibling**" = d[21:22, ], "**C-section**" = d[23:24, ])

dim(d)


final_table <- 
  kable(d, booktabs = T, caption = table_caption_descr,  escape = F, longtable = F) %>%
    kable_styling() %>%
    group_rows("Sample size", 1,1) %>%
    group_rows("Childsex", 2,3) %>%
    group_rows("Infant age at CC entrance (weeks)", 4,6) %>%
    group_rows("Maternal age", 7,9) %>%
    group_rows("Birthweight (gram)", 10, 12) %>%
    footnote(general = table_note_descr, general_title = "Note.", footnote_as_chunk = T, threeparttable = T) 
    


## ------------------------------------------------------------------------
# group_by(df_table, group_BF_noBF) %>% summarise(n= n())
# group_by(df_table, Group) %>% summarise(n= n())
# group_by(df_table, group_strict) %>% summarise(n= n())
# group_by(df_table, group) %>% summarise(n= n())


## ------------------------------------------------------------------------
pm_table
table_caption_pm <- "Model parameters for two-way PERMANOVAs"
table_note_pm <- "CC = Childcare entrance."
pm_table_apa <- 
  kable(pm_table, , align = c("l", rep("r", 5)), booktabs = T, caption = table_caption_pm,  escape = T) %>%
    kable_styling(latex_options = c("scale_down")) %>%
    footnote(general = table_note_pm, general_title = "Note.", footnote_as_chunk = T, threeparttable = T)
pm_table_apa


save.image(here("data/tables.RData"))
