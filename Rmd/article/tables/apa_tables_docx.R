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

## ------------------------------------------------------------------------
# BIBO confounders dataset 
confounders <- read.spss("~/workspace/research_master/dpb/minor_research_project/article/analyses/data/bibo_confounders.sav", to.data.frame = T)
# select relevant columns
confounders <- select(confounders, -group193)
# rename for convenience
colnames(confounders) <- c("ID", "birthweight", "siblings", "pre_smoking",
 "pre_alcohol", "apgar5min", "maternal_age", "delivery", "gestat_length",
  "childsex", "maternal_education", "firstborn")
# Carolina sent me also this file because the confounders file is incomplete for usage:
june2017_gerben <- read_excel("~/workspace/research_master/dpb/minor_research_project/article/analyses/data/june2017_gerben.xlsx", sheet = "Data")
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
# select final sample 
df <- filter(df_complete, 
  delivery == "natuurlijke bevalling" | delivery == "pomp", apgar5min > 7,
  begin_group == 1, group_strict != 99999
)
# change numericals to factors
df <- df %>%
  mutate_at(
    vars(c("group_strict", "group_BF_noBF", "CC_noCC", "CC", "childsex", "firstborn")),
    funs(factor)
  )

colnames(df)
# make df for table
df_table <- df 
colnames(df_table) <- c('ID', 'begin_group', 'group', 'Group',
 'group_BF_noBF', 'CC_noCC', 'weeks_BF', 'percent_BF_kdv1_kdv2', 'CC',
 'age_CC_min2d', 'age_CC_plus28d', 'age_begin_weeks', 'age_plus_4weeks',
 'birthweight', 'siblings', 'pre_smoking', 'pre_alcohol', 'apgar5min',
 'maternal_age', 'delivery', 'gestat_length', 'childsex',
 'maternal_education', 'firstborn')
levels(df_table$Group) <- c("BF + CC", "no BF + CC", "BF + no CC", "no BF + no CC")
head(df_table)
colnames(df_table)
df_table$age_begin_weeks

# create df with descriptives 
d <- group_by(df_table, Group) %>%
  summarise(
    N = n(),
    Male = sum(childsex == 1),
    Female = sum(childsex == 0),
    mean_age_begin_weeks = mean(age_begin_weeks, na.rm = T),
    sd_age_begin_weeks = sd(age_begin_weeks, na.rm = T),
    `mean (sd)_age_begin_weeks` = NA,
    min_age_begin_weeks = min(age_begin_weeks),
    max_age_begin_weeks = max(age_begin_weeks),
    mean_maternal_age = mean(maternal_age, na.rm = T),
    sd_maternal_age = sd(maternal_age, na.rm = T),
    `mean (sd)_maternal_age` = NA,
    min_maternal_age = min(maternal_age), 
    max_maternal_age = max(maternal_age), 
    mean_birthweight = mean(birthweight),
    sd_birthweight = sd(birthweight, na.rm = T),
    `mean (sd)_birthweight` = NA,
    min_birthweight = min(birthweight), 
    max_birthweight  = max(birthweight)) %>%
    mutate_if(is.numeric, funs(format(round(., 1), nsmall = 0))) %>%
    mutate(
      `mean (sd)_age_begin_weeks` = glue("{mean_age_begin_weeks} $\\pm$ {sd_age_begin_weeks}"),
      `mean (sd)_maternal_age` = glue("{mean_maternal_age} $\\pm$ {sd_maternal_age}"),
      `mean (sd)_birthweight` = glue("{mean_birthweight} $\\pm$ {sd_birthweight}"),
    ) %>%
    select(-mean_age_begin_weeks, -sd_age_begin_weeks, -mean_maternal_age, -sd_maternal_age, -mean_birthweight, -sd_birthweight)
d
    
# obtain order of colnames since spread does not preserve it 
# but I need to keep the order for the rows laters
ordered_rows <- colnames(d)
d
# back to wide format, reorder, names ready for table, merge sd + mean to str

d <- 
  mutate_if(d, is.numeric, funs(toString(.))) %>%  
    gather(variable, value, -Group) %>%
    spread(Group, value) %>%
    arrange(match(variable, ordered_rows[-1])) %>%
    mutate(variable = gsub("_.+", "", variable))

d  
colnames(d)[1] <- ""
table_caption_descr <- "Descriptive statistics for demographic variables of infants and mothers included in the present study."
table_note_descr <- "CC = childcare. BF = breastfeeding. Two way analyses of variances did not show significant differences in the mean between the groups of childcare center and breastfeeding for the continuous measures."
final_table <- 
  kable(d, booktabs = T, caption = table_caption_descr,  escape = F, longtable = F) %>%
    kable_styling() %>%
    group_rows("Sample size", 1,1) %>%
    group_rows("Childsex", 2,3) %>%
    group_rows("Infant age at CC entrance (weeks)", 4,6) %>%
    group_rows("Maternal age", 7,9) %>%
    group_rows("Birthweight (gram)", 10, 12) %>%
    footnote(general = table_note_descr, general_title = "Note.", footnote_as_chunk = T, threeparttable = T) 

d
## ------------------------------------------------------------------------
# group_by(df_table, group_BF_noBF) %>% summarise(n= n())
# group_by(df_table, Group) %>% summarise(n= n())
# group_by(df_table, group_strict) %>% summarise(n= n())
# group_by(df_table, group) %>% summarise(n= n())

## ------------------------------------------------------------------------
source("~/workspace/research_master/dpb/minor_research_project/article/analyses/gerben/adj_linear mixed model for Henrik to measure between groups changes.R")

## ------------------------------------------------------------------------
library(knitr)
library(kableExtra)

## ------------------------------------------------------------------------
permanova_table
table_caption_pm <- "Model parameters for two-way PERMANOVAs"
table_note_pm <- "CC = Childcare entrance."
pm_table <- 
  kable(permanova_table, , align = c("l", rep("r", 5)), booktabs = T, caption = table_caption_pm,  escape = T) %>%
    kable_styling(latex_options = c("scale_down")) %>%
    group_rows("Before CC", 1,5) %>%
    group_rows("After CC", 6,10) %>%
    footnote(general = table_note_pm, general_title = "Note.", footnote_as_chunk = T, threeparttable = T)

## ------------------------------------------------------------------------
# f table for LME4
table_caption_lme <- "Log fold change within childcare groups for significant linear mixed effects models. "
table_note_lme <- "CC = Childcare center. FDR =  Pvalue adjustment based on Benjamini Hochberg correction."
f_table <- 
  kable(my_f_table, booktabs = T, caption = table_caption_lme,  escape = T, longtable = F) %>%
    kable_styling() %>%
    footnote(general = table_note_lme, general_title = "Note.", footnote_as_chunk = T, threeparttable = T)
#latex_options = c("scale_down") 




library(papaja)







# there are extra steps required for papaja here because when using rowselection
# R adds and extra colum to tibble. Therefore, I make two sep table first 

perm_1 <- 
  mutate_if(pm1, is.numeric, funs(round(., 3))) %>%
    mutate(
      MeanSqs = ifelse(is.na(MeanSqs), "-", MeanSqs),
      F.Model = ifelse(is.na(F.Model), "-", F.Model),
      `Pr(>F)` = ifelse(is.na(`Pr(>F)`), "-", `Pr(>F)`) 
    ) %>%
    select(`Model Parameters`,SumsOfSqs, MeanSqs, F.Model, Df,`Pr(>F)`, R2)

colnames(perm_1) <- c("Model Parameter", "Sum of Squares", "Mean Sum of Squares", "F", "Df", "p", "R Square" )
# same for part 2
perm_2 <- 
  mutate_if(pm2, is.numeric, funs(round(., 3))) %>%
    mutate(
      MeanSqs = ifelse(is.na(MeanSqs), "-", MeanSqs),
      F.Model = ifelse(is.na(F.Model), "-", F.Model),
      `Pr(>F)` = ifelse(is.na(`Pr(>F)`), "-", `Pr(>F)`) 
    ) %>%
    select(`Model Parameters`,SumsOfSqs, MeanSqs, F.Model, Df,`Pr(>F)`, R2)

colnames(perm_2) <- c("Model Parameter", "Sum of Squares", "Mean Sum of Squares", "F", "Df", "p", "R Square" )



