## ------------------------------------------------------------------------
#knitr::purl(here("Rmd/childcare_article_analyses.Rmd"), here("R/childcare_article_analyses.R"))

## ------------------------------------------------------------------------
library(vegan)
library(microbiome)
library(tidyverse)
library(here)
library(future)
library(furrr)
library(broom)
library(ggpubr)
library(brms)
library(mice)

## ------------------------------------------------------------------------
# load data and helper functions
source("https://raw.githubusercontent.com/HenrikEckermann/in_use/master/bayesian_helper.R")
source("https://raw.githubusercontent.com/HenrikEckermann/in_use/master/mb_helper.R")
source("https://raw.githubusercontent.com/HenrikEckermann/in_use/master/reporting.R")

## ------------------------------------------------------------------------
load(here("data/data_transfer.RData"))
source(here("R/read.R"))

## ------------------------------------------------------------------------
# take over the meta variables I created in other docs
meta_new <- data_transfer[, 1:9] 

## ------------------------------------------------------------------------
# create catories for bf and childcare and specifically for ccyes vs rest
meta_new <- meta_new %>%
  mutate(
      groups = ifelse(time == "pre" & cc == "no", "noCCpre", ifelse(
          time == "pre" & cc == "yes", "CCpre", ifelse(
              time == "post" & cc == "no", "noCCpost", "CCpost"))),
      bf = ifelse(bf_ratio <= 0.25, "lowBF", ifelse(
          bf_ratio <0.75, "mediumBF", "highBF")),
      ccpost = ifelse(groups == "CCpost", 1, 0)) %>% 
  mutate(
      groups = as.factor(groups), 
      bf = as.factor(bf), 
      ccpost = as.factor(ccpost))

## ------------------------------------------------------------------------
siblings <- readxl::read_excel(here("data/meta_data/siblings_addon_yvonne.xlsx")) %>% 
    filter(subject_id %in% meta_new$subject_id)

## ------------------------------------------------------------------------
# add confounding variables
confounders <- foreign::read.spss(here("data/meta_data/bibo_confounders.sav"), to.data.frame = T) %>%
    select(subject_id = ID, childsex, delivery = DELIVERYmode, birthweight = BIRTHWEIGTH) %>%
    filter(subject_id %in% meta_new$subject_id) %>% 
    left_join(siblings, by = "subject_id") %>% 
    mutate(
        sibling = ifelse(sibling == 0, 0, 1),
        csection = ifelse(
        delivery == "natuurlijke bevalling", 0,
            ifelse(delivery == "6", NA, 
               ifelse(delivery == "999", NA,
                      ifelse(delivery == "pomp", 0, 1)))))
meta_new <- meta_new %>% left_join(confounders, by = "subject_id")

## ------------------------------------------------------------------------
# create new pseq object (read.R results in the object "genus" Leo created)
otu <- otu_to_df(genus, transpose = FALSE)
otu <- otu %>% 
    select(species, meta_new$sample_id) %>% 
    df_to_otu()
pseq <- phyloseq(otu, df_to_sd(meta_new), tax_table(genus))
# add diversity indeces to sample data
diversities <- 
    global(pseq, index = "all") %>% 
    select(contains("diversities")) %>% 
    rownames_to_column("sample_id")
colnames(diversities) <- gsub("diversities_", "", colnames(diversities))

sample_data(pseq) <- 
    sd_to_df(pseq) %>% 
    left_join(diversities, by = "sample_id") %>%
    df_to_sd()
meta <- sd_to_df(pseq)
# clr and relative abundance transformation to deal with compositionality of mb data
pseq.clr <- microbiome::transform(pseq, transform = "clr")
pseq.rel <- microbiome::transform(pseq, "compositional")

## ------------------------------------------------------------------------
otus.clr <- otu_to_df(pseq.clr)
colnames(otus.clr)[which(colnames(otus.clr) == "Clostridium \\(sensu stricto\\)")] <- "Clostridium_sensu_stricto"
colnames(otus.clr) <- c("sample_id", gsub("_", "", colnames(otus.clr)[-1]))
colnames(otus.clr) <- gsub("\\.", "", colnames(otus.clr))
colnames(otus.clr) <- gsub(" ", "", colnames(otus.clr))
genus <- colnames(otus.clr)[-1]
data <- sd_to_df(pseq.clr) %>%
    left_join(otus.clr, by = "sample_id")
data$childsex <- as.factor(data$childsex)
data$delivery <- as.factor(data$delivery)
data <- data %>% mutate(birthweight_s = scale(birthweight)[, 1])


# PCA with CLR values (euclidean distance of clr transformed values = Aitchison distance) 
pcx <- prcomp(otus.clr %>% column_to_rownames("sample_id"))
# extract loadings
pcx_rot <- 
    pcx$rotation %>%
        as.tibble() %>%
        mutate_all(function(x) x*10) %>%
        add_column(genus = rownames(pcx$rotation))

# add PCs to data
princomps <- pcx$x %>% as.data.frame() %>%
    rownames_to_column("sample_id") %>%
    select(PC1, PC2, PC3, PC4, PC5, sample_id)
data <- data %>% left_join(princomps, by = "sample_id") 
                   
# how much variance do pcs explain?
pc1 <- round(pcx$sdev[1]^2/sum(pcx$sdev^2),2)
pc2 <- round(pcx$sdev[2]^2/sum(pcx$sdev^2),2)
pc3 <- round(pcx$sdev[3]^2/sum(pcx$sdev^2),2)
pc4 <- round(pcx$sdev[4]^2/sum(pcx$sdev^2),2)
pc5 <- round(pcx$sdev[5]^2/sum(pcx$sdev^2),2)                   



## ------------------------------------------------------------------------
# I recode the contrasts so that I have the comparisons I want:
# the intercept will reflect our group we want to compare to others (cc post)
# the cc coefficient then compares to no cc post, the time coefficent 
# to cc pre and the interaction to no cc pre
# contrasts(data$cc)[1, 1] <- 1
# contrasts(data$cc)[2, 1] <- 0
# contrasts(data$time)[1, 1] <- 1
# contrasts(data$time)[2, 1] <- 0


# Next I impute data using predictive mean matching (PMM). PMM is
# less difficult to specify. I use the PCs to impute instead of all 
# genus abundances since these are correlated
data_imp <- data %>% 
    select(
        -everything(), 
        subject_id, 
        age_d_s, 
        time, 
        cc, 
        bf_count_s, 
        sibling,
        birthweight_s,
        csection,
        PC1,
        PC2,
        PC3,
        PC4,
        PC5) %>%
    mice(m = 10, method = "pmm", print = F, seed = 412) %>%
    mice::complete("all")
# This I use to join genus abundances again
deselect_col <- colnames(data_imp[[1]])
data_lj <- data %>% select(-deselect_col, subject_id, time)
data_imp <- map(data_imp, ~.x %>% left_join(data_lj, by = c("subject_id", "time")))

## ------------------------------------------------------------------------
# PMM does not understand that for one subject the sibling/csection must be similar, so I manually always use the pre
# sibling which is random.
data_imp <- map(data_imp, function(x) {
    #x[x$subject_id == "453" & x$time == "post", "sibling"] <- x[x$subject_id == "453" & x$time == "pre", "sibling"]
    sibling_pre <- x[x$subject_id == 453 & x$time == "pre", "sibling"]
    csection_pre_449 <- x[x$subject_id == 449 & x$time == "pre", "csection"]
    csection_pre_369 <- x[x$subject_id == 369 & x$time == "pre", "csection"]
    x <- x %>% mutate(
        sibling = ifelse(subject_id == "453", sibling_pre, sibling),
        csection = ifelse(subject_id == "449", csection_pre_449, 
                          ifelse(subject_id == "369", csection_pre_369, csection)),
        sibling = as.factor(sibling),
        csection = as.factor(csection)
        
    )

})

## ------------------------------------------------------------------------
data_imp$`1` %>% group_by(sibling, cc) %>% summarise(n= n())


prior <- c(
            set_prior("normal(0, 2)", class = "b"),
            set_prior("exponential(25)", class = "sd"),
            set_prior("normal(0, 10)", class = "Intercept"),
            set_prior("lkj(2)", class = "cor"))
prior_refit <- c(
            set_prior("normal(0, 2)", class = "b"),
            set_prior("exponential(55)", class = "sd"),
            set_prior("normal(0, 10)", class = "Intercept"),
            set_prior("lkj(2)", class = "cor"))


control <-  list(adapt_delta = 0.9999, max_treedepth = 15)
folder <- here("models/differential_abundance/skew_normal/mice/")

# define fitting function for fixed sigma for mice object
brm_sn_multiple_fut_comp <- function(genus) future({
    tryCatch(
        {
            # specify formula
            f_d <- as.formula(glue("{genus} ~ 1 + time*cc + age_d_s + bf_count_s + csection * sibling + (1+ time + age_d_s + bf_count_s|subject_id)"))
            formula <- bf(f_d)
            # give individual model name for storage
            model_file <- glue("{folder}csection/{genus}_full_multiple_comp")
            #fit model
            brm_multiple(
                family = skew_normal(), data = data_imp, formula = formula,
                chains = 4, warmup = 1000,
                control = control, prior = prior, file = model_file
                )
        },
        error=function(cond) {
            # Choose a return value in case of error
            return(cond)
        }
    )
})

## ------------------------------------------------------------------------
#plan(multiprocess)
#map(files[excluded_ids], brm_sn_multiple_fut_refit)
#models <- future_map(models_fut, value)

## ------------------------------------------------------------------------

finished <- list.files(here("models/differential_abundance/skew_normal/mice/csection"))
finished <- gsub("_full_multiple_comp.rds", "", finished)
finished %>% length()
plan(multiprocess)
models_fut <- future_map(genus[!genus %in% finished], brm_sn_multiple_fut_comp)
