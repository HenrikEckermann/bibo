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
formula_sigma <- sigma ~ 1 + time*cc + sibling*csection
formula <- bf(formula_mu, formula_sigma)
# see default priors brms
get_prior(formula, data = data_imp[[1]], family = student)
# we use default prior except for using a normal(0, 1) prior for the b
# which makes the model more skeptical of large effect sizes

genera <- colnames(data_imp[[1]])[29:158]
exclude_vector <- c()
# fit models 
comparisons <- map(genera, function(genus) {
  
    # specify model
    formula_mu <- glue("{genus} ~ time*cc + age_d_s + bf_count_s + sibling*csection + (1|subject_id)") %>% as.formula()
    formula_sigma <- sigma ~ 1 + time*cc
    formula <- bf(formula_mu, formula_sigma)
    
    # fit model
    fit <- brm_multiple(
        data = data_imp,
        family = student(),
        formula = bf(formula_mu, formula_sigma),
        prior = c(prior(normal(0, 1), class = b),
                  prior(normal(0, 4), class = Intercept, dpar = sigma),
                  prior(normal(0, 0.2), class = b, dpar = sigma), 
                  prior(gamma(2, 0.2), class = nu)),
        cores = 4,
        file = glue("{BIBO}/models/bayesian_student/student_sigma_{genus}")
    ) 
    
    # create a grid so we can get the mu for each subgroup
    nd <- with(data_imp[[1]], 
               expand.grid(cc = levels(cc), time = levels(time),
                           sibling = levels(sibling), csection = levels(csection)))
    # we use median value for bf
    nd$bf_count_s <- median(model.frame(fit)$bf_count_s)
    # it makes no sense that age remains similar when changing pre to post
    # so we use the median for pre and post
    nd <- nd %>% mutate(age_d_s = ifelse(time == "pre", -0.835, 0.626))
    
    # only add mu for genera if diagnostics are OK
    if (return_diag(fit, genus)) {
      
        # calulate mu for each subgroup
        plinp_mu <- posterior_linpred(fit, nd, transform = T,  re.form = ~ 0) %>% as_tibble()
        colnames(plinp_mu) <- c(
            "home_pre_nosib_nocsec",
            "cc_pre_nosib_nocsec",
            "home_post_nosib_nocsec",
            "cc_post_nosib_nocsec",
            "home_pre_sib_nocsec",
            "cc_pre_sib_nocsec",
            "home_post_sib_nocsec", 
            "cc_post_sib_nocsec",
            "home_pre_nosib_csec",
            "cc_pre_nosib_csec", 
            "home_post_nosib_csec",
            "cc_post_nosib_csec",
            "home_pre_sib_csec",
            "cc_pre_sib_csec",
            "home_post_sib_csec",
            "cc_post_sib_csec")
            
        # mu comparisons (csec group separately)
        mu_group <- plinp_mu %>%
            mutate(
                genus = genus, 
                
                # ccpre vs homepre 
                homepre = (home_pre_nosib_nocsec + home_pre_sib_nocsec)/2,
                ccpre = (cc_pre_nosib_nocsec + cc_pre_sib_nocsec)/2,
                ccpre_homepre = ccpre - homepre,
                
                # ccpost vs ccpre
                ccpost = (cc_post_nosib_nocsec + cc_post_sib_nocsec)/2,
                ccpost_ccpre = ccpost -ccpre,
                
                # ccpost vs homepost
                homepost = (home_post_nosib_nocsec + home_post_sib_nocsec)/2,
                ccpost_homepost = ccpost - homepost,
                
                # homepost vs homepre
                homepost_homepre = homepost - homepre,
                
                # sibling (if csec: home_pre_nosib_csec + cc_pre_nosib_csec + home_post_nosib_csec + cc_post_nosib_csec)
                nosib = (home_pre_nosib_nocsec + cc_pre_nosib_nocsec + home_post_nosib_nocsec + cc_post_nosib_nocsec)/4,
                sib = (home_pre_sib_nocsec + cc_pre_sib_nocsec + home_post_sib_nocsec + cc_post_sib_nocsec)/4,
                sib_nosib = sib - nosib,
                
                # csec
                nocsec = (home_pre_nosib_nocsec + cc_pre_nosib_nocsec + home_post_nosib_nocsec + cc_post_nosib_nocsec +
                          home_pre_sib_nocsec + cc_pre_sib_nocsec + home_post_sib_nocsec + cc_post_sib_nocsec)/8,
                csec = (home_pre_nosib_csec + cc_pre_nosib_csec + home_post_nosib_csec + cc_post_nosib_csec +
                          home_pre_sib_csec + cc_pre_sib_csec + home_post_sib_csec + cc_post_sib_csec)/8,
                csec_nocsec = csec - nocsec,
                
                # sib * csec
                csec_sib = (home_pre_sib_csec + cc_pre_sib_csec + home_post_sib_csec + cc_post_sib_csec)/4,
                csec_nosib = (home_pre_nosib_csec + cc_pre_nosib_csec + home_post_nosib_csec + cc_post_nosib_csec)/4,
                csecsib_csecnosib = csec_sib - csec_nosib
            ) %>%
        select(genus, ccpre_homepre, ccpost_ccpre, ccpost_homepost, homepost_homepre, sib_nosib, csec_nocsec, csecsib_csecnosib)
            
        
        
        # calulate sigma for each subgroup
        plinp_sigma <- posterior_linpred(fit, nd, transform = T,  re.form = ~ 0, dpar = "sigma") %>% as_tibble()
        colnames(plinp_sigma) <- c(
            "home_pre_nosib_nocsec",
            "cc_pre_nosib_nocsec",
            "home_post_nosib_nocsec",
            "cc_post_nosib_nocsec",
            "home_pre_sib_nocsec",
            "cc_pre_sib_nocsec",
            "home_post_sib_nocsec", 
            "cc_post_sib_nocsec",
            "home_pre_nosib_csec",
            "cc_pre_nosib_csec", 
            "home_post_nosib_csec",
            "cc_post_nosib_csec",
            "home_pre_sib_csec",
            "cc_pre_sib_csec",
            "home_post_sib_csec",
            "cc_post_sib_csec")
            
        # sigma comparisons (csec group separately)
        sigma_group <- plinp_sigma %>%
            mutate(
                genus = genus, 
                
                # ccpre vs homepre 
                homepre = (home_pre_nosib_nocsec + home_pre_sib_nocsec)/2,
                ccpre = (cc_pre_nosib_nocsec + cc_pre_sib_nocsec)/2,
                ccpre_homepre = ccpre - homepre,
                
                # ccpost vs ccpre
                ccpost = (cc_post_nosib_nocsec + cc_post_sib_nocsec)/2,
                ccpost_ccpre = ccpost -ccpre,
                
                # ccpost vs homepost
                homepost = (home_post_nosib_nocsec + home_post_sib_nocsec)/2,
                ccpost_homepost = ccpost - homepost,
                
                # homepost vs homepre
                homepost_homepre = homepost - homepre,
                
                # sibling (if csec: home_pre_nosib_csec + cc_pre_nosib_csec + home_post_nosib_csec + cc_post_nosib_csec)
                nosib = (home_pre_nosib_nocsec + cc_pre_nosib_nocsec + home_post_nosib_nocsec + cc_post_nosib_nocsec)/4,
                sib = (home_pre_sib_nocsec + cc_pre_sib_nocsec + home_post_sib_nocsec + cc_post_sib_nocsec)/4,
                sib_nosib = sib - nosib,
                
                # csec
                nocsec = (home_pre_nosib_nocsec + cc_pre_nosib_nocsec + home_post_nosib_nocsec + cc_post_nosib_nocsec +
                          home_pre_sib_nocsec + cc_pre_sib_nocsec + home_post_sib_nocsec + cc_post_sib_nocsec)/8,
                csec = (home_pre_nosib_csec + cc_pre_nosib_csec + home_post_nosib_csec + cc_post_nosib_csec +
                          home_pre_sib_csec + cc_pre_sib_csec + home_post_sib_csec + cc_post_sib_csec)/8,
                csec_nocsec = csec - nocsec,
                
                # sib * csec
                csec_sib = (home_pre_sib_csec + cc_pre_sib_csec + home_post_sib_csec + cc_post_sib_csec)/4,
                csec_nosib = (home_pre_nosib_csec + cc_pre_nosib_csec + home_post_nosib_csec + cc_post_nosib_csec)/4,
                csecsib_csecnosib = csec_sib - csec_nosib
            ) %>%
        select(genus, ccpre_homepre, ccpost_ccpre, ccpost_homepost, homepost_homepre, sib_nosib, csec_nocsec, csecsib_csecnosib)
        
        return(list(mu_group, sigma_group))
        
        
          
    } else {
        print(glue("{genus} excluded (diagnostics)"))
        exclude_vector <- c(exclude_vector, genus)
        return(NA)
    }
})


save(comparisons, exclude_vector, file = glue("{BIBO}/rdata/bayesian_student.rds"))

