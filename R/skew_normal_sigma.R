## README: this script for written for torque cluster. You 
##         need to submit the qsub command within the 
##         orl_major_project folder. Otherwise, if you want to
##         reproduce, make sure to change {ORL}

#------ ENVIRONMENT PREPARATION

# for the torque cluster I need to tell where my library is...
libloc <- Sys.getenv("R_LIBS")
# dont use R_LIBS on local machine
if (libloc != "") {
  pkgs <- list("tidyverse", "here", "future", "furrr", "glue", "brms")
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


# load prepared data
load(glue("{BIBO}/rdata/test_data2.rds"))
# load own functions
source("https://raw.githubusercontent.com/HenrikEckermann/in_use/master/mb_helper.R")
source("https://raw.githubusercontent.com/HenrikEckermann/in_use/master/reporting.R")


