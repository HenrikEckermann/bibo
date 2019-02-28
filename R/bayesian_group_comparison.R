## README: this script for written for torque cluster. You 
##         need to submit the qsub command within the 
##         orl_major_project folder. Otherwise, if you want to
##         reproduce, make sure to change {ORL}

#------ ENVIRONMENT PREPARATION

# for the torque cluster I need to tell where my library is...
libloc <- Sys.getenv("R_LIBS")
# dont use R_LIBS on local machine
if (libloc != "") {
  pkgs <- list("tidyverse", "here", "furrr", "glue", "brms", "future")
  lapply(pkgs, function(pkg) {
    library(pkg, character.only = TRUE, lib.loc = libloc)
  })
  } else {
  pkgs <- list("tidyverse", "here", "furrr", "glue", "brms")
  lapply(pkgs, function(pkg) {
    library(pkg, character.only = TRUE)
  })
}

# learn time 
start = Sys.time()
print(start)

BIBO <- Sys.getenv("BIBO")

#------ ACTUAL CODING

# load prepared data
load(glue("{BIBO}/rdata/data_mi.rds"))
# load own functions
source("https://raw.githubusercontent.com/HenrikEckermann/in_use/master/bayesian_helper.R")
source("https://raw.githubusercontent.com/HenrikEckermann/in_use/master/mb_helper.R")
source("https://raw.githubusercontent.com/HenrikEckermann/in_use/master/reporting.R")



# define the formula we use at the example of Bifidobacterium
formula_mu <- glue("Bifidobacterium ~ time*cc + age_d_s + bf_count_s + sibling*csection + (1|subject_id)") %>% as.formula()
formula_sigma <- sigma ~ 1 + time*cc
formula <- bf(formula_mu, formula_sigma)
# see default priors brms
get_prior(formula, data = data_imp[[1]], family = student)
# we use default prior except for using a normal(0, 1) prior for the b
# which makes the model more skeptical of large effect sizes

genera <- colnames(data_imp[[1]])[29:158]

# fit models 
models <- map(genera, function(genus) {
    formula_mu <- glue("{genus} ~ time*cc + age_d_s + bf_count_s + sibling*csection + (1|subject_id)") %>% as.formula()
    formula_sigma <- sigma ~ 1 + time*cc
    formula <- bf(formula_mu, formula_sigma)
    
    # fit example model
    fit <- brm_multiple(
        data = data_imp,
        family = student(),
        formula = bf(formula_mu, formula_sigma),
        prior = c(prior(normal(0, 1), class = b),
                  prior(normal(0, 4), class = Intercept, dpar = sigma),
                  prior(normal(0, 0.2), class = b, dpar = sigma), 
                  prior(gamma(2, 0.2), class = nu)),
        cores = 4,
        file = glue("{BIBO}/models/bayesian_student/student_{genus}")
    ) 
})

# store vector of succeful fits 
diagnosed <- map2(models, genera, return_diag)

save(models, diagnosed, file = glue{"{BIBO}/rdata/bayesian_student.rds"})