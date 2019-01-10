library(purrr)
library(mice)



est5 <- nhanes %>%
  mice(seed = 123, print = FALSE) %>%
  mice::complete("all") %>%
  map(lm, formula = chl ~ age + bmi + hyp) 

est5 %>% pool()
imp <- mice(nhanes, seed = 123, print = FALSE)
fit <- with(imp, lm(chl ~ age + bmi + hyp))
est1 <- pool(fit)
est1

nhanes
library(mice)
imp <- mice(nhanes, seed = 123, print = FALSE)
fit <- with(imp, lme4::lmer(chl ~ + bmi + hyp + (1|age)))
est1 <- summary(pool(fit))
est1

sessionInfo()

imp <- mice(nhanes, m=5, print = FALSE, seed = 55152)
fit <- with(imp, lm(bmi~hyp+chl))
summary(pool(fit))
library(purrr)
test <- nhanes %>%
  mice(seed = 123, print = FALSE) %>%
  mice::complete("all") %>%
  map(lm, formula = chl ~ age + bmi + hyp) %>%
  pool() %>%
  summary()
test 
nhanes %>%
  mice(seed = 123, print = FALSE) %>%
  mice::complete("all") %>%
  map(lm, formula = chl ~ age + bmi + hyp) %>%
  map(summary)


summary(lm(formula = chl ~ age + bmi + hyp, data = nhanes))
