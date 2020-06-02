-   [Aim](#aim)
-   [Packages](#packages)
    -   [Install packages](#install-packages)
    -   [Load Packages](#load-packages)
    -   [A function to calculate the pseudo
        R-squared](#a-function-to-calculate-the-pseudo-r-squared)
-   [Results obtained using the Stockholm CT-score in the core
    model](#results-obtained-using-the-stockholm-ct-score-in-the-core-model)
    -   [results](#results)
    -   [Results for each individual
        variable.](#results-for-each-individual-variable.)
-   [ordinal](#ordinal)

Aim
===

The aim of this script is to use SVM on the MRI data in order to predict
outcome.

Packages
========

Install packages
----------------

This step only has to be done once!

``` r
install.packages("tidymodels")
install.packages("tidyverse")
install.packages("knitr")
install.packages("AppliedPredictiveModeling")
install.packages("caret")
install.packages("MLmetrics")
install.packages("ggfortify")
install.packages("kableExtra")
install.packages("summarytools")
install.packages("readxl")
install.packages("corrr", dependencies = T)
install.packages("RColorBrewer")
install.packages("corrplot")
#This may be required for the Support Vector Machines analysis
install.packages("kernlab")
install.packages('modelgrid')
install.packages('e1071')
install.packages('MASS')
install.packages('descr')
install.packages("DescTools")
install.packages("CaTools")
install.packages("table1", dependencies = T)
install.packages("mice")
install.packages("pscl")
install.packages("naniar")
install.packages("glmnet")
install.packages("pROC")
install.packages("stringi")
install.packages("randomForest")
install.packages("mltools")
install.packages("data.table")
install.packages("cowplot")
install.packages("GA")
install.packages("funModeling")
install.packages("desirability")
install.packages("doParallel")
```

Load Packages
-------------

``` r
library(MASS)
library(data.table)
library(tidymodels)
library(tidyverse)
library(knitr)
library(AppliedPredictiveModeling)
library(caret)
library(MLmetrics)
library(ggfortify)
library(kableExtra)
library(summarytools)
library(corrplot)
library(readxl)
library(modelgrid)
library(e1071)
library(descr)
library(DescTools)
library(caTools)
library(table1)
library(corrr)
library(lubridate)
library(mice)
library(pscl)
library(naniar)
library(stringi)
library(glmnet)
library(pROC)
library(mltools)
library(cowplot)
library(GA)
library(funModeling)
library(doParallel)


options(knitr.kable.NA = '')
```

A function to calculate the pseudo R-squared
--------------------------------------------

This function will used to calculate Nagelkerke’s Pseudo-*R*<sup>2</sup>
in the multivariate logistic regression models, which will be fitted
later on in this report.

``` r
# function used to calculate R squared

r2_calc <- function(data, model, nullmodel){
  CoxSnell <- 1 - exp(-(2/nrow(data) * (logLik(model) -  logLik(nullmodel))))
  r2_full <- CoxSnell/(1 - exp((2 * nrow(data)^(-1)) * logLik(nullmodel)))
  return(r2_full)
}
```

Results obtained using the Stockholm CT-score in the core model
===============================================================

``` r
core_model <- c("Alder", "GCSKS", "hypoxia", "hypotension", "Pupiller")

imputed_log <- readRDS('imputed_data.rds') %>% 
  mice::complete("long") %>%
    group_by(.imp) %>%
  select("Imputation" =.imp, id, dich_gos,
         adams, stockholm, rotterdam, firsching, hamdeh, 
          stockholm_mri,
         core_model, mortality) %>%
   group_by(Imputation) %>%  nest()

core_vars_formula <- paste("dich_gos ~", paste(core_model, 
                                               collapse = " + "),
                            sep = " " ) %>% 
   as.formula()

core_stockholm <-  paste("dich_gos ~", paste(core_model, 
                                             collapse = " + "),
                            sep = " " , "+", "stockholm") %>% 
   as.formula()

core_adams <-  paste("dich_gos ~", paste(core_model, 
                                             collapse = " + "),
                            sep = " " , "+", "adams") %>% 
   as.formula()

core_karolinska <-  paste("dich_gos ~", paste(core_model, 
                                             collapse = " + "),
                            sep = " " , "+", "stockholm_mri") %>% 
   as.formula()

core_firsching <-  paste("dich_gos ~", paste(core_model, 
                                             collapse = " + "),
                            sep = " " , "+", "firsching") %>% 
   as.formula()

core_hamdeh <-  paste("dich_gos ~", paste(core_model, 
                                             collapse = " + "),
                            sep = " " , "+", "hamdeh") %>% 
   as.formula()

var_names <- readRDS('imputed_data.rds') %>% 
  mice::complete("long") %>% 
   select(adams, stockholm_mri, hamdeh, firsching) %>%
   names() 

formulas_stockholm <- paste('dich_gos ~ ', paste(core_model, 
                                                 collapse = " + "),
                            sep = " " , "+", "stockholm", "+",var_names) #formulas to be used in the the glm function

formulas_univariate <- paste('dich_gos ~ ', var_names)
```

fitting the models

``` r
estimates <- imputed_log %>%
mutate(adams = map(data, ~ glm(formulas_stockholm[1], 
                               family = binomial(link="logit"), 
                               data = .))) %>%
  mutate(hamdeh = map(data, ~ glm(formulas_stockholm[3], 
                               family = binomial(link="logit"), 
                               data = .))) %>% 
mutate(firsching = map(data, ~ glm(formulas_stockholm[4], 
                               family = binomial(link="logit"), 
                               data = .))) %>% 
  mutate(karolinska = map(data, ~ glm(formulas_stockholm[2], 
                                      family = binomial(link="logit"), data = .))) %>%
  mutate(null = map(data, ~ glm(dich_gos ~ 1, 
                              family = binomial(link="logit"), 
                              data = .))) %>%
    mutate(stockholm = map(data, ~ glm(dich_gos ~ stockholm, 
                              family = binomial(link="logit"), 
                              data = .))) %>% 
   mutate(core = map(data, ~ glm(core_vars_formula, family = binomial(link="logit"), data = .))) %>%
   mutate(core_stockholm = map(data, ~ glm(core_stockholm, family = binomial(link="logit"), data = .))) %>%
   mutate(core_adams = map(data, ~ glm(core_adams, family = binomial(link="logit"), data = .))) %>%
   mutate(core_karolinska = map(data, ~ glm(core_karolinska, family = binomial(link="logit"), data = .))) %>%
   mutate(core_firsching = map(data, ~ glm(core_firsching, family = binomial(link="logit"), data = .))) %>%
   mutate(core_hamdeh = map(data, ~ glm(core_hamdeh, family = binomial(link="logit"), data = .))) %>%
   mutate(r_sq_karolinska = pmap_dbl(list(data, karolinska, null), r2_calc)) %>%
  mutate(r_sq_hamdeh = pmap_dbl(list(data, hamdeh, null), 
                                r2_calc)) %>%   
  mutate(r_sq_firsching = pmap_dbl(list(data, firsching, null), 
                                 r2_calc)) %>% 
   mutate(r_sq_core = pmap_dbl(list(data, core, null), r2_calc)) %>%
  mutate(r_sq_stockholm = pmap_dbl(list(data, stockholm, null), r2_calc)) %>%
  mutate(r_sq_core_stockholm = pmap_dbl(list(data, core_stockholm, null), r2_calc)) %>%
  mutate(r_sq_core_adams = pmap_dbl(list(data, core_adams, null), 
                                    r2_calc)) %>%
   mutate(r_sq_core_firsching = pmap_dbl(list(data, core_firsching, null), 
                                    r2_calc)) %>%
   mutate(r_sq_core_hamdeh = pmap_dbl(list(data, core_hamdeh, null), 
                                    r2_calc)) %>%
  mutate(r_sq_core_karolinska = pmap_dbl(list(data, core_karolinska, null), r2_calc)) %>%
  mutate(r_sq_adams = pmap_dbl(list(data, adams, null), r2_calc)) %>%
mutate(glance = map(karolinska, ~ glance(.))) %>%
mutate(aug = map2(karolinska, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>% 
mutate(tid = map(karolinska, ~ tidy(.))) %>%
mutate(glance_core = map(core, ~ glance(.))) %>% 
mutate(glance_core_stockholm = map(core_stockholm,
                                   ~ glance(.))) %>% 
  mutate(glance_stockholm = map(stockholm,
                                   ~ glance(.))) %>% 
mutate(glance_core_adams = map(core_adams, ~ glance(.))) %>% 
mutate(glance_core_hamdeh = map(core_hamdeh, ~ glance(.))) %>% 
mutate(glance_core_firsching = map(core_firsching, ~ glance(.))) %>% 
mutate(glance_core_karolinska = map(core_karolinska, ~ glance(.))) %>%
mutate(glance_adams = map(adams, ~ glance(.))) %>%
  mutate(glance_hamdeh = map(hamdeh, ~ glance(.))) %>%
   mutate(glance_firsching = map(firsching, ~ glance(.))) %>%
   mutate(aug_core = map2(core, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>% 
mutate(aug_core_stockholm = map2(core_stockholm, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>% 
  mutate(aug_stockholm = map2(stockholm, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>% 
mutate(aug_core_adams = map2(core_adams, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>% 
mutate(aug_core_hamdeh = map2(core_hamdeh, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>% 
mutate(aug_core_firsching = map2(core_firsching, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>% 
mutate(aug_core_karolinska = map2(core_karolinska, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>%
mutate(aug_adams = map2(adams, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>%
  mutate(aug_hamdeh = map2(hamdeh, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>%
  mutate(aug_firsching = map2(firsching, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>%
  select(- c(adams, core_adams, core_stockholm, 
             core_karolinska,
       core, hamdeh, firsching, karolinska, core_firsching, core_hamdeh, stockholm)) %>% 
    mutate(karolinska_only = map(data, ~ glm(formulas_univariate[2], 
                                      family = binomial(link="logit"), data = .))) %>% 
  mutate(adams_only = map(data, ~ glm(formulas_univariate[1], 
                                      family = binomial(link="logit"), data = .))) %>% 
  mutate(hamdeh_only = map(data, ~ glm(formulas_univariate[3], 
                                      family = binomial(link="logit"), data = .))) %>% 
  mutate(firsching_only = map(data, ~ glm(formulas_univariate[4], 
                                      family = binomial(link="logit"), data = .))) %>% 
  mutate(r_sq_karolinska_only = pmap_dbl(list(data, karolinska_only, null), r2_calc)) %>%
  mutate(r_sq_adams_only = pmap_dbl(list(data, adams_only, null), r2_calc)) %>%
  mutate(r_sq_hamdeh_only = pmap_dbl(list(data, hamdeh_only, null), r2_calc)) %>%
  mutate(r_sq_firsching_only = pmap_dbl(list(data, firsching_only, null), r2_calc)) %>%
  mutate(glance_hamdeh_only = map(hamdeh_only, ~ glance(.))) %>% 
  mutate(glance_firsching_only = map(firsching_only, ~ glance(.))) %>% 
   mutate(glance_karolinska_only = map(karolinska_only,
                                       ~ glance(.))) %>%
  mutate(glance_adams_only = map(adams_only, ~ glance(.))) %>% 
  mutate(aug_firsching_only = map2(firsching_only, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>%
  mutate(aug_karolinska_only = map2(karolinska_only, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>%
  mutate(aug_hamdeh_only = map2(hamdeh_only, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>%
  mutate(aug_adams_only = map2(adams_only, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>%
  mutate(aug_karolinska_only = map2(karolinska_only, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>%
  select(-c(null, karolinska_only, hamdeh_only, 
        adams_only, firsching_only)) %>% 
     mutate_at(vars(contains("aug")), funs(map(., ~ roc(.x,dich_gos,
                                              .fitted)))) %>% 
  mutate_at(vars(contains("aug")), funs(map_dbl(., ~ auc(.x))))

names(estimates) <- gsub("aug", "auc", names(estimates))

model_vars <- estimates %>%  ungroup() %>% 
  sample_n(size =  1, replace = F) %>% 
   unnest(tid, .drop = T) %>% 
      pull(term) %>% 
     as_vector()
```

Likelihood ratio test to compare the addition of MRI to cora variables
and the Stockholm CT-score

``` r
data_table <- imputed_log %>% 
  filter(Imputation == 3) %>% 
  select(data) %>% 
  unnest()

adams <- glm(formulas_stockholm[1], 
                               family = binomial(link="logit"), 
                               data = data_table)


hamdeh <- glm(formulas_stockholm[3], 
                               family = binomial(link="logit"), 
                               data = data_table)

firsching <- glm(formulas_stockholm[4], 
                               family = binomial(link="logit"), 
                               data = data_table) 

karolinska <- glm(formulas_stockholm[2], 
                                      family = binomial(link="logit"), data = data_table)

core_stockholm <- glm(core_stockholm, family = binomial(link="logit"), data = data_table) 


karolinska_only <- glm(formulas_univariate[2], 
                                      family = binomial(link="logit"), data = data_table) 

adams_only <- glm(formulas_univariate[1], 
                                      family = binomial(link="logit"),
                    data = data_table)

hamdeh_only <- glm(formulas_univariate[3], 
                                      family = binomial(link="logit"),
                    data = data_table)

firsching_only <- glm(formulas_univariate[4], 
                                      family = binomial(link="logit"),
                    data = data_table)


core_karolinska <- anova(core_stockholm, karolinska , test="LRT")

tidy(core_karolinska) %>% 
  kable(digits = 3, caption = "comparison between CT score and Stockholm mri")
```

<table>
<caption>
comparison between CT score and Stockholm mri
</caption>
<thead>
<tr>
<th style="text-align:right;">
Resid..Df
</th>
<th style="text-align:right;">
Resid..Dev
</th>
<th style="text-align:right;">
df
</th>
<th style="text-align:right;">
Deviance
</th>
<th style="text-align:right;">
p.value
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:right;">
343
</td>
<td style="text-align:right;">
362.632
</td>
<td style="text-align:right;">
</td>
<td style="text-align:right;">
</td>
<td style="text-align:right;">
</td>
</tr>
<tr>
<td style="text-align:right;">
340
</td>
<td style="text-align:right;">
331.372
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
31.261
</td>
<td style="text-align:right;">
0
</td>
</tr>
</tbody>
</table>

``` r
core_hamdeh <- anova(core_stockholm, hamdeh , test="LRT")

tidy(core_hamdeh) %>% 
 kable(digits = 3, caption = "comparison between CT score and hamdeh")
```

<table>
<caption>
comparison between CT score and hamdeh
</caption>
<thead>
<tr>
<th style="text-align:right;">
Resid..Df
</th>
<th style="text-align:right;">
Resid..Dev
</th>
<th style="text-align:right;">
df
</th>
<th style="text-align:right;">
Deviance
</th>
<th style="text-align:right;">
p.value
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:right;">
343
</td>
<td style="text-align:right;">
362.632
</td>
<td style="text-align:right;">
</td>
<td style="text-align:right;">
</td>
<td style="text-align:right;">
</td>
</tr>
<tr>
<td style="text-align:right;">
339
</td>
<td style="text-align:right;">
346.220
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
16.412
</td>
<td style="text-align:right;">
0.003
</td>
</tr>
</tbody>
</table>

``` r
core_firsching <- anova(core_stockholm, firsching , test="LRT")

tidy(core_firsching) %>% 
  kable(digits = 3, caption = "comparison between CT score and firsching")
```

<table>
<caption>
comparison between CT score and firsching
</caption>
<thead>
<tr>
<th style="text-align:right;">
Resid..Df
</th>
<th style="text-align:right;">
Resid..Dev
</th>
<th style="text-align:right;">
df
</th>
<th style="text-align:right;">
Deviance
</th>
<th style="text-align:right;">
p.value
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:right;">
343
</td>
<td style="text-align:right;">
362.632
</td>
<td style="text-align:right;">
</td>
<td style="text-align:right;">
</td>
<td style="text-align:right;">
</td>
</tr>
<tr>
<td style="text-align:right;">
339
</td>
<td style="text-align:right;">
342.423
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
20.209
</td>
<td style="text-align:right;">
0
</td>
</tr>
</tbody>
</table>

``` r
core_adams <- anova(core_stockholm, adams , test="LRT")

tidy(core_adams) %>% 
 kable(digits = 3, caption = "comparison between CT score and adams")
```

<table>
<caption>
comparison between CT score and adams
</caption>
<thead>
<tr>
<th style="text-align:right;">
Resid..Df
</th>
<th style="text-align:right;">
Resid..Dev
</th>
<th style="text-align:right;">
df
</th>
<th style="text-align:right;">
Deviance
</th>
<th style="text-align:right;">
p.value
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:right;">
343
</td>
<td style="text-align:right;">
362.632
</td>
<td style="text-align:right;">
</td>
<td style="text-align:right;">
</td>
<td style="text-align:right;">
</td>
</tr>
<tr>
<td style="text-align:right;">
340
</td>
<td style="text-align:right;">
345.493
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
17.14
</td>
<td style="text-align:right;">
0.001
</td>
</tr>
</tbody>
</table>

DeLong test to compare AUC values of the different MRI-based grading
systems

``` r
roc_karo <- augment(karolinska_only) %>% 
  roc(., dich_gos,.fitted)


roc_adams <- augment(adams_only) %>% 
  roc(., dich_gos,.fitted)

roc_hamdeh <- augment(hamdeh_only) %>% 
  roc(., dich_gos,.fitted)

roc_firsching <- augment(firsching_only) %>% 
  roc(., dich_gos,.fitted)


roc.test(roc_adams, roc_karo, method=c("delong"), boot.n=2000, boot.stratified=TRUE )
```

    ## 
    ##  DeLong's test for two correlated ROC curves
    ## 
    ## data:  roc_adams and roc_karo
    ## Z = -1.6268, p-value = 0.1038
    ## alternative hypothesis: true difference in AUC is not equal to 0
    ## sample estimates:
    ## AUC of roc1 AUC of roc2 
    ##   0.6812224   0.7130009

``` r
roc.test(roc_adams, roc_firsching, method=c("delong"), boot.n=2000, boot.stratified=TRUE )
```

    ## 
    ##  DeLong's test for two correlated ROC curves
    ## 
    ## data:  roc_adams and roc_firsching
    ## Z = 0, p-value = 1
    ## alternative hypothesis: true difference in AUC is not equal to 0
    ## sample estimates:
    ## AUC of roc1 AUC of roc2 
    ##   0.6812224   0.6812224

``` r
roc.test(roc_adams, roc_hamdeh, method=c("delong"), boot.n=2000, boot.stratified=TRUE )
```

    ## 
    ##  DeLong's test for two correlated ROC curves
    ## 
    ## data:  roc_adams and roc_hamdeh
    ## Z = -0.10397, p-value = 0.9172
    ## alternative hypothesis: true difference in AUC is not equal to 0
    ## sample estimates:
    ## AUC of roc1 AUC of roc2 
    ##   0.6812224   0.6823120

``` r
roc.test(roc_firsching, roc_hamdeh, method=c("delong"), boot.n=2000, boot.stratified=TRUE )
```

    ## 
    ##  DeLong's test for two correlated ROC curves
    ## 
    ## data:  roc_firsching and roc_hamdeh
    ## Z = -0.071385, p-value = 0.9431
    ## alternative hypothesis: true difference in AUC is not equal to 0
    ## sample estimates:
    ## AUC of roc1 AUC of roc2 
    ##   0.6812224   0.6823120

``` r
roc.test(roc_firsching, roc_karo, method=c("delong"), boot.n=2000, boot.stratified=TRUE )
```

    ## 
    ##  DeLong's test for two correlated ROC curves
    ## 
    ## data:  roc_firsching and roc_karo
    ## Z = -1.5206, p-value = 0.1284
    ## alternative hypothesis: true difference in AUC is not equal to 0
    ## sample estimates:
    ## AUC of roc1 AUC of roc2 
    ##   0.6812224   0.7130009

``` r
roc.test(roc_hamdeh, roc_karo, method=c("delong"), boot.n=2000, boot.stratified=TRUE )
```

    ## 
    ##  DeLong's test for two correlated ROC curves
    ## 
    ## data:  roc_hamdeh and roc_karo
    ## Z = -1.504, p-value = 0.1326
    ## alternative hypothesis: true difference in AUC is not equal to 0
    ## sample estimates:
    ## AUC of roc1 AUC of roc2 
    ##   0.6823120   0.7130009

Bootstrap auc confidence intervals

``` r
order_vars <- c("core", "stockholm", "adams_only", 
               "firsching_only", "hamdeh_only", "karolinska_only",
               "core_stockholm", "core_adams", "core_firsching",
               "core_hamdeh", "core_karolinska","adams", 
               "firsching", "hamdeh", "karolinska")

set.seed(123)

auc_ci <- estimates %>% 
  mutate(adams = map(data, ~ glm(formulas_stockholm[1], 
                               family = binomial(link="logit"), 
                               data = .))) %>%
  mutate(hamdeh = map(data, ~ glm(formulas_stockholm[3], 
                               family = binomial(link="logit"), 
                               data = .))) %>% 
mutate(firsching = map(data, ~ glm(formulas_stockholm[4], 
                               family = binomial(link="logit"), 
                               data = .))) %>% 
  mutate(karolinska = map(data, ~ glm(formulas_stockholm[2], 
                                      family = binomial(link="logit"), data = .))) %>%
    mutate(stockholm = map(data, ~ glm(dich_gos ~ stockholm, 
                              family = binomial(link="logit"), 
                              data = .))) %>% 
   mutate(core = map(data, ~ glm(core_vars_formula, family = binomial(link="logit"), data = .))) %>%
   mutate(core_stockholm = map(data, ~ glm(core_stockholm, family = binomial(link="logit"), data = .))) %>%
   mutate(core_adams = map(data, ~ glm(core_adams, family = binomial(link="logit"), data = .))) %>%
   mutate(core_karolinska = map(data, ~ glm(core_karolinska, family = binomial(link="logit"), data = .))) %>%
   mutate(core_firsching = map(data, ~ glm(core_firsching, family = binomial(link="logit"), data = .))) %>%
   mutate(core_hamdeh = map(data, ~ glm(core_hamdeh, family = binomial(link="logit"), data = .))) %>%
  mutate(aug = map2(karolinska, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>% 
   mutate(aug_core = map2(core, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>% 
mutate(aug_core_stockholm = map2(core_stockholm, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>% 
  mutate(aug_stockholm = map2(stockholm, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>% 
mutate(aug_core_adams = map2(core_adams, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>% 
mutate(aug_core_hamdeh = map2(core_hamdeh, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>% 
mutate(aug_core_firsching = map2(core_firsching, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>% 
mutate(aug_core_karolinska = map2(core_karolinska, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>%
mutate(aug_adams = map2(adams, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>%
  mutate(aug_hamdeh = map2(hamdeh, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>%
  mutate(aug_firsching = map2(firsching, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>%
  select(- c(adams, core_adams, core_stockholm, 
             core_karolinska,
       core, hamdeh, firsching, karolinska, core_firsching, core_hamdeh, stockholm)) %>% 
    mutate(karolinska_only = map(data, ~ glm(formulas_univariate[2], 
                                      family = binomial(link="logit"), data = .))) %>% 
  mutate(adams_only = map(data, ~ glm(formulas_univariate[1], 
                                      family = binomial(link="logit"), data = .))) %>% 
  mutate(hamdeh_only = map(data, ~ glm(formulas_univariate[3], 
                                      family = binomial(link="logit"), data = .))) %>% 
  mutate(firsching_only = map(data, ~ glm(formulas_univariate[4], 
                                      family = binomial(link="logit"), data = .))) %>% 
  mutate(aug_firsching_only = map2(firsching_only, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>%
  mutate(aug_karolinska_only = map2(karolinska_only, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>%
  mutate(aug_hamdeh_only = map2(hamdeh_only, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>%
  mutate(aug_adams_only = map2(adams_only, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>%
  mutate(aug_karolinska_only = map2(karolinska_only, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>%
  select(-c(karolinska_only, hamdeh_only, 
        adams_only, firsching_only)) %>% 
     mutate_at(vars(contains("aug")), funs(map(., ~ roc(.x,dich_gos,
                                              .fitted)))) %>% 
  mutate_at(vars(contains("aug")), funs(map(., ~ ci.auc(.x, method = "bootstrap")))) %>% 
  ungroup() %>% 
  select(Imputation, contains("aug")) %>% 
  gather(key = "model", value = "ci", -Imputation) %>% 
  mutate(model = str_remove(model, "aug_")) %>% 
  mutate(model = str_replace(model, "aug", "karolinska")) %>% 
  separate_rows(ci) %>%
    filter(!ci %in% c('c', '', 2.5,50,97.5)) %>% 
  mutate(ci = as.numeric(ci)) %>% 
  group_by(Imputation, model) %>% 
  summarise(low_ci = min(ci), high_ci = max(ci)) %>% 
  ungroup() %>% 
  mutate(low_ci = as.numeric(low_ci), 
         high_ci = as.numeric(high_ci)) %>%
  group_by(model) %>% 
  summarise(low = mean(low_ci), high = mean(high_ci)) %>% 
  ungroup() %>% 
   mutate(low = round(low, 3), high = round(high, 3)) %>% 
  mutate(auc_ci = paste("(", low, "-", high, ")", sep = "")) %>%
  select(model, auc_ci) %>% 
  mutate(model = as.factor(model)) %>% 
  arrange(factor(model, levels = order_vars))
```

Bootstrapping pseudo R squared

``` r
set.seed(123)

bootstrapped <- imputed_log %>% 
  filter(Imputation == 3) %>% 
  ungroup() %>% 
  select(data) %>% 
  unnest(data) %>% 
  bootstraps(.,times = 2000) %>% 
  mutate(data = map(splits, ~ analysis(.x))) %>% 
  select(-splits)
```

running all the glm models and calculating nagelkerke

``` r
bootstrapped_results <- bootstrapped %>%
mutate(adams = map(data, ~ glm(formulas_stockholm[1], 
                               family = binomial(link="logit"), 
                               data = .))) %>%
  mutate(hamdeh = map(data, ~ glm(formulas_stockholm[3], 
                               family = binomial(link="logit"), 
                               data = .))) %>% 
mutate(firsching = map(data, ~ glm(formulas_stockholm[4], 
                               family = binomial(link="logit"), 
                               data = .))) %>% 
  mutate(karolinska = map(data, ~ glm(formulas_stockholm[2], 
                                      family = binomial(link="logit"), data = .))) %>%
  mutate(null = map(data, ~ glm(dich_gos ~ 1, 
                              family = binomial(link="logit"), 
                              data = .))) %>%
    mutate(stockholm = map(data, ~ glm(dich_gos ~ stockholm, 
                              family = binomial(link="logit"), 
                              data = .))) %>% 
   mutate(core = map(data, ~ glm(core_vars_formula, family = binomial(link="logit"), data = .))) %>%
   mutate(core_stockholm = map(data, ~ glm(core_stockholm, family = binomial(link="logit"), data = .))) %>%
   mutate(core_adams = map(data, ~ glm(core_adams, family = binomial(link="logit"), data = .))) %>%
   mutate(core_karolinska = map(data, ~ glm(core_karolinska, family = binomial(link="logit"), data = .))) %>%
   mutate(core_firsching = map(data, ~ glm(core_firsching, family = binomial(link="logit"), data = .))) %>%
   mutate(core_hamdeh = map(data, ~ glm(core_hamdeh, family = binomial(link="logit"), data = .))) %>%
   mutate(r_sq_karolinska = pmap_dbl(list(data, karolinska, null), r2_calc)) %>%
  mutate(r_sq_hamdeh = pmap_dbl(list(data, hamdeh, null), 
                                r2_calc)) %>%   
  mutate(r_sq_firsching = pmap_dbl(list(data, firsching, null), 
                                 r2_calc)) %>% 
   mutate(r_sq_core = pmap_dbl(list(data, core, null), r2_calc)) %>%
  mutate(r_sq_stockholm = pmap_dbl(list(data, stockholm, null), r2_calc)) %>%
  mutate(r_sq_core_stockholm = pmap_dbl(list(data, core_stockholm, null), r2_calc)) %>%
  mutate(r_sq_core_adams = pmap_dbl(list(data, core_adams, null), 
                                    r2_calc)) %>%
   mutate(r_sq_core_firsching = pmap_dbl(list(data, core_firsching, null), 
                                    r2_calc)) %>%
   mutate(r_sq_core_hamdeh = pmap_dbl(list(data, core_hamdeh, null), 
                                    r2_calc)) %>%
  mutate(r_sq_core_karolinska = pmap_dbl(list(data, core_karolinska, null), r2_calc)) %>% 
  mutate(r_sq_adams = pmap_dbl(list(data, adams, null), r2_calc)) %>%
  select(- c(adams, core_adams, core_stockholm, 
             core_karolinska,
       core, hamdeh, firsching, karolinska, core_firsching, core_hamdeh, stockholm)) %>% 
    mutate(karolinska_only = map(data, ~ glm(formulas_univariate[2], 
                                      family = binomial(link="logit"), data = .))) %>% 
  mutate(adams_only = map(data, ~ glm(formulas_univariate[1], 
                                      family = binomial(link="logit"), data = .))) %>% 
  mutate(hamdeh_only = map(data, ~ glm(formulas_univariate[3], 
                                      family = binomial(link="logit"), data = .))) %>% 
  mutate(firsching_only = map(data, ~ glm(formulas_univariate[4], 
                                      family = binomial(link="logit"), data = .))) %>% 
  mutate(r_sq_karolinska_only = pmap_dbl(list(data, karolinska_only, null), r2_calc)) %>%
  mutate(r_sq_adams_only = pmap_dbl(list(data, adams_only, null), r2_calc)) %>%
  mutate(r_sq_hamdeh_only = pmap_dbl(list(data, hamdeh_only, null), r2_calc)) %>%
  mutate(r_sq_firsching_only = pmap_dbl(list(data, firsching_only, null), r2_calc)) %>% 
  select(id, contains("r_sq")) %>% 
  gather(key = "model", value = "r_sq", -id) %>% 
  group_by(model) %>% 
  summarise(low = round(quantile(r_sq, 0.025),3),
              high = round(quantile(r_sq, 0.975),3), sd = sd(r_sq), mean = mean(r_sq)) %>% 
  ungroup() %>% 
  mutate(model = str_remove(model, "r_sq_")) %>% 
  mutate(me = 1.96 * mean) %>% 
  mutate(rsq_ci = paste("(", low, "-", high, ")", sep = "")) %>%
  select(model, rsq_ci) %>% 
  mutate(model = as.factor(model)) %>% 
  arrange(factor(model, levels = order_vars))
```

results
-------

Results for all logistic regression models that were fitted to all
imputed datasets, respectively.

``` r
estimates %>%
  unnest(glance) %>%
  dplyr::select(Imputation, r_sq_karolinska, 
                null.deviance:df.residual) %>% 
  kable(format = "latex", escape = FALSE, booktabs = T, linesep = "",
        caption = "Comparing the results for all imputed datasets",
        digits = 2, align = 'c', col.names = c("Imputation", "$R^2$",
                                               "Null deviance",
                                                   "df.null", "logLik",
                                                   "AIC", "BIC",
                    "deviance", "df.residuals")) %>%
   kable_styling(full_width = F)
```

average (CI)

``` r
core_comparison <- estimates %>%
  unnest(glance) %>%  
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_karolinska, AIC, 
                "AUC" = auc) %>% 
summarize_all(funs(mean)) %>% 
  mutate_at(vars(c("r_sq", "AUC", "AIC")), funs(round(., 3))) %>% 
mutate(Model = "Core + Stockholm CT-score + Stockholm MRI grading system") 

core_comparison <- estimates %>%
  unnest(glance_hamdeh) %>% 
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_hamdeh,AIC, 
                "AUC" = auc_hamdeh) %>% 
summarize_all(funs(mean)) %>% 
    mutate_at(vars(c("r_sq", "AUC", "AIC")), funs(round(., 3))) %>% 
   mutate(Model = "Core + Stockholm CT-score + MRI grading system of Abu Hamdeh et al.") %>% 
   bind_rows(core_comparison) %>% 
   select(Model, everything())

core_comparison <- estimates %>%
  unnest(glance_firsching) %>% 
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_firsching,AIC, 
                "AUC" = auc_firsching) %>% 
summarize_all(funs(mean)) %>% 
    mutate_at(vars(c("r_sq", "AUC", "AIC")), funs(round(., 3))) %>% 
   mutate(Model = "Core + Stockholm CT-score + MRI grading system of Firsching et al.") %>% 
   bind_rows(core_comparison) %>% 
   select(Model, everything())

core_comparison <- estimates %>%
  unnest(glance_adams) %>% 
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_adams, AIC, "AUC" = auc_adams) %>% 
summarize_all(funs(mean)) %>% 
    mutate_at(vars(c("r_sq", "AUC", "AIC")), funs(round(., 3))) %>% 
mutate(Model = "Core + Stockholm CT-score + MRI grading system of Adams et al.") %>% 
   bind_rows(core_comparison) %>% 
   select(Model, everything())

core_comparison <- estimates %>%
  unnest(glance_core_karolinska) %>% 
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_core_karolinska,AIC, 
                "AUC" = auc_core_karolinska) %>% 
summarize_all(funs(mean)) %>% 
    mutate_at(vars(c("r_sq", "AUC", "AIC")), funs(round(., 3))) %>% 
   mutate(Model = "Core + Stockholm MRI grading system") %>% 
   bind_rows(core_comparison) %>% 
   select(Model, everything())

core_comparison <- estimates %>%
  unnest(glance_core_hamdeh) %>% 
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_core_hamdeh,AIC, 
                "AUC" = auc_core_hamdeh) %>% 
summarize_all(funs(mean)) %>% 
    mutate_at(vars(c("r_sq", "AUC", "AIC")), funs(round(., 3))) %>% 
mutate(Model = "Core + MRI grading system of Abu Hamdeh et al.") %>% 
   bind_rows(core_comparison) %>% 
   select(Model, everything())

core_comparison <- estimates %>%
  unnest(glance_core_firsching) %>% 
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_core_firsching,AIC, 
                "AUC" = auc_core_firsching) %>% 
summarize_all(funs(mean)) %>% 
    mutate_at(vars(c("r_sq", "AUC", "AIC")), funs(round(., 3))) %>% 
mutate(Model = "Core + MRI grading system of Firsching et al.") %>% 
   bind_rows(core_comparison) %>% 
   select(Model, everything())

core_comparison <- estimates %>%
  unnest(glance_core_adams) %>% 
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_core_adams, AIC, 
                "AUC" = auc_core_adams) %>% 
summarize_all(funs(mean)) %>% 
    mutate_at(vars(c("r_sq", "AUC", "AIC")), funs(round(., 3))) %>% 
   mutate(Model = "Core + MRI grading system of Adams et al.") %>% 
   bind_rows(core_comparison) %>% 
   select(Model, everything())

core_comparison <- estimates %>%
  unnest(glance_core_stockholm) %>%
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_core_stockholm, AIC, 
                "AUC" = auc_core_stockholm) %>% 
summarize_all(funs(mean)) %>% 
    mutate_at(vars(c("r_sq", "AUC", "AIC")), funs(round(., 3))) %>% 
   mutate(Model = "Core + Stockholm CT-score") %>% 
   bind_rows(core_comparison) %>% 
   select(Model, everything())

core_comparison <- estimates %>%
  unnest(glance_karolinska_only) %>%  
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_karolinska_only, AIC, 
                "AUC" = auc_karolinska_only) %>% 
summarize_all(funs(mean)) %>% 
    mutate_at(vars(c("r_sq", "AUC", "AIC")), funs(round(., 3))) %>% 
   mutate(Model = "Stockholm MRI grading system") %>% 
   bind_rows(core_comparison) %>% 
   select(Model, everything())

core_comparison <- estimates %>%
  unnest(glance_hamdeh_only) %>%  
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_hamdeh_only, AIC, 
                "AUC" = auc_hamdeh_only) %>% 
summarize_all(funs(mean)) %>% 
    mutate_at(vars(c("r_sq", "AUC", "AIC")), funs(round(., 3))) %>% 
   mutate(Model = "MRI grading system of Abu Hamdeh et al.") %>% 
   bind_rows(core_comparison) %>% 
   select(Model, everything())

core_comparison <- estimates %>%
  unnest(glance_firsching_only) %>%  
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_firsching_only, AIC, 
                "AUC" = auc_firsching_only) %>% 
summarize_all(funs(mean)) %>% 
    mutate_at(vars(c("r_sq", "AUC", "AIC")), funs(round(., 3))) %>% 
   mutate(Model = "MRI grading system of Firsching et al.") %>% 
   bind_rows(core_comparison) %>% 
   select(Model, everything())

core_comparison <- estimates %>%
  unnest(glance_adams_only) %>%  
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_adams_only, AIC, 
                "AUC" = auc_adams_only) %>% 
summarize_all(funs(mean)) %>% 
    mutate_at(vars(c("r_sq", "AUC", "AIC")), funs(round(., 3))) %>% 
   mutate(Model = "MRI grading system of Adams et al.") %>% 
   bind_rows(core_comparison) %>% 
   select(Model, everything())

core_comparison <- estimates %>%
  unnest(glance_stockholm) %>%  
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_stockholm, AIC, 
                "AUC" = auc_stockholm) %>% 
summarize_all(funs(mean)) %>% 
    mutate_at(vars(c("r_sq", "AUC", "AIC")), funs(round(., 3))) %>% 
   mutate(Model = "Stockholm CT-score") %>% 
   bind_rows(core_comparison) %>% 
   select(Model, everything())

core_comparison <- estimates %>%
  unnest(glance_core) %>%  
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_core, AIC, 
                "AUC" = auc_core) %>% 
summarize_all(funs(mean)) %>% 
    mutate_at(vars(c("r_sq", "AUC", "AIC")), funs(round(., 3))) %>% 
   mutate(Model = "Core") %>% 
   bind_rows(core_comparison) %>% 
   select(Model, everything())

core_comparison <- core_comparison %>% 
  mutate_at(vars(c("r_sq", "AUC", "AIC")), funs(round(., 3))) %>%
  select(-AIC) %>% 
  mutate(rsq_ci = bootstrapped_results$rsq_ci, 
         auc_ci = auc_ci$auc_ci) %>% 
  unite("r_sq", c("r_sq", "rsq_ci"), sep = " ") %>% 
   unite("AUC", c("AUC", "auc_ci"), sep = " ")
  

core_comparison %>% 
  kable(format = "latex", escape = FALSE, booktabs = T, linesep = "",
        caption = "Average result for all imputed datasets",
        digits = 2, align = 'c', col.names = c("Model","$R^2$ (95% CI)",
                                                "AUC (95% CI)")) %>%
   kable_styling(full_width = F)

results_table_stockholm <- core_comparison %>%
  kable(format = "html", booktabs = T, linesep = "",
        caption = "Average result for all imputed datasets",
        digits = 3, align = 'c', col.names = c("Model","R^2 (95% CI)",
                                                "AUC (95% CI)")) %>%
   kable_styling(full_width = F)



# save_kable(results_table_stockholm, "stockholm_results.html")
```

average

``` r
core_comparison <- estimates %>%
  unnest(glance) %>%  
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_karolinska, AIC, 
                "AUC" = auc) %>% 
summarize_all(funs(mean)) %>% 
mutate(Model = "Core + Stockholm CT-score + Stockholm MRI grading system") 

core_comparison <- estimates %>%
  unnest(glance_hamdeh) %>% 
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_hamdeh,AIC, 
                "AUC" = auc_hamdeh) %>% 
summarize_all(funs(mean)) %>% 
   mutate(Model = "Core + Stockholm CT-score + MRI grading system of Abu Hamdeh et al.") %>% 
   bind_rows(core_comparison) %>% 
   select(Model, everything())

core_comparison <- estimates %>%
  unnest(glance_firsching) %>% 
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_firsching,AIC, 
                "AUC" = auc_firsching) %>% 
summarize_all(funs(mean)) %>% 
   mutate(Model = "Core + Stockholm CT-score + MRI grading system of Firsching et al.") %>% 
   bind_rows(core_comparison) %>% 
   select(Model, everything())

core_comparison <- estimates %>%
  unnest(glance_adams) %>% 
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_adams, AIC, "AUC" = auc_adams) %>% 
summarize_all(funs(mean)) %>% 
mutate(Model = "Core + Stockholm CT-score + MRI grading system of Adams et al.") %>% 
   bind_rows(core_comparison) %>% 
   select(Model, everything())

core_comparison <- estimates %>%
  unnest(glance_core_karolinska) %>% 
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_core_karolinska,AIC, 
                "AUC" = auc_core_karolinska) %>% 
summarize_all(funs(mean)) %>% 
   mutate(Model = "Core + Stockholm MRI grading system") %>% 
   bind_rows(core_comparison) %>% 
   select(Model, everything())

core_comparison <- estimates %>%
  unnest(glance_core_hamdeh) %>% 
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_core_hamdeh,AIC, 
                "AUC" = auc_core_hamdeh) %>% 
summarize_all(funs(mean)) %>% 
mutate(Model = "Core + MRI grading system of Abu Hamdeh et al.") %>% 
   bind_rows(core_comparison) %>% 
   select(Model, everything())

core_comparison <- estimates %>%
  unnest(glance_core_firsching) %>% 
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_core_firsching,AIC, 
                "AUC" = auc_core_firsching) %>% 
summarize_all(funs(mean)) %>% 
mutate(Model = "Core + MRI grading system of Firsching et al.") %>% 
   bind_rows(core_comparison) %>% 
   select(Model, everything())

core_comparison <- estimates %>%
  unnest(glance_core_adams) %>% 
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_core_adams, AIC, 
                "AUC" = auc_core_adams) %>% 
summarize_all(funs(mean)) %>% 
   mutate(Model = "Core + MRI grading system of Adams et al.") %>% 
   bind_rows(core_comparison) %>% 
   select(Model, everything())

core_comparison <- estimates %>%
  unnest(glance_core_stockholm) %>%
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_core_stockholm, AIC, 
                "AUC" = auc_core_stockholm) %>% 
summarize_all(funs(mean)) %>% 
   mutate(Model = "Core + Stockholm CT-score") %>% 
   bind_rows(core_comparison) %>% 
   select(Model, everything())

core_comparison <- estimates %>%
  unnest(glance_karolinska_only) %>%  
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_karolinska_only, AIC, 
                "AUC" = auc_karolinska_only) %>% 
summarize_all(funs(mean)) %>% 
   mutate(Model = "Stockholm MRI grading system") %>% 
   bind_rows(core_comparison) %>% 
   select(Model, everything())

core_comparison <- estimates %>%
  unnest(glance_hamdeh_only) %>%  
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_hamdeh_only, AIC, 
                "AUC" = auc_hamdeh_only) %>% 
summarize_all(funs(mean)) %>% 
   mutate(Model = "MRI grading system of Abu Hamdeh et al.") %>% 
   bind_rows(core_comparison) %>% 
   select(Model, everything())

core_comparison <- estimates %>%
  unnest(glance_firsching_only) %>%  
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_firsching_only, AIC, 
                "AUC" = auc_firsching_only) %>% 
summarize_all(funs(mean)) %>%
   mutate(Model = "MRI grading system of Firsching et al.") %>% 
   bind_rows(core_comparison) %>% 
   select(Model, everything())

core_comparison <- estimates %>%
  unnest(glance_adams_only) %>%  
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_adams_only, AIC, 
                "AUC" = auc_adams_only) %>% 
summarize_all(funs(mean)) %>% 
   mutate(Model = "MRI grading system of Adams et al.") %>% 
   bind_rows(core_comparison) %>% 
   select(Model, everything())

core_comparison <- estimates %>%
  unnest(glance_stockholm) %>%  
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_stockholm, AIC, 
                "AUC" = auc_stockholm) %>% 
summarize_all(funs(mean)) %>% 
   mutate(Model = "Stockholm CT-score") %>% 
   bind_rows(core_comparison) %>% 
   select(Model, everything())

core_comparison <- estimates %>%
  unnest(glance_core) %>%  
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_core, AIC, 
                "AUC" = auc_core) %>% 
summarize_all(funs(mean)) %>% 
   mutate(Model = "Core") %>% 
   bind_rows(core_comparison) %>% 
   select(Model, everything())
  

core_comparison %>% 
  kable(format = "latex", escape = FALSE, booktabs = T, linesep = "",
        caption = "Average result for all imputed datasets",
        digits = 2, align = 'c', col.names = c("Model","$R^2$", 
                                               "AIC",
                                                "AUC")) %>%
   kable_styling(full_width = F)
```

``` r
# results_table_stockholm_orig <- core_comparison %>%
#   kable(format = "html", booktabs = T, linesep = "",
#         caption = "Average result for all imputed datasets",
#         digits = 3, align = 'c', col.names = c("Model","R^2", 
#         "AIC",
#                                                 "AUC")) %>%
#    kable_styling(full_width = F)


# save_kable(results_table_stockholm_orig, "stockholm_results_orig.html")
```

Results for each individual variable.
-------------------------------------

Results for all the imputed datasets combined

``` r
by_variable <- imputed_log %>%
  mutate(karolinska = map(data, ~ glm(formulas_stockholm[2], 
                                      family = binomial(link="logit"), data = .)))

vars_table <- summary(mice::pool(by_variable$karolinska)) %>%
as_tibble() %>%
mutate(variable = model_vars) %>%
mutate(OR = exp(estimate)) %>%
filter(variable != "(Intercept)") %>%
dplyr::select(variable, everything(), -c(statistic, estimate, df)) 

# Two different coefficients are presented for the Pupils variable,
# as this variable has three levels. 
#Thus, a new vector containing two different 
# names for the pupils is required

vars_table %>% mutate(variable = str_replace_all(string = variable,
                                                 pattern = "\\d",
                                  replacement = "")) %>%
   select(-term) %>% 
   arrange(desc(OR)) %>% 
  kable(caption = "Combined results of all imputated datasets", align = 'c',
        format = "latex", booktabs = T, linesep = "", longtable = TRUE,  digits = 2) %>%
   kable_styling(full_width = F, latex_options = "repeat_header")
```

``` r
# vars_latex <- vars_table %>% mutate(variable = str_replace_all(string = variable,
#                                                  pattern = "\\d",
#                                   replacement = "")) %>%
#    select(-term) %>%
#    arrange(desc(OR)) %>%
#   kable(caption = "Combined results of all imputated datasets", align = 'c',
#         format = "latex", booktabs = T, linesep = "", longtable = TRUE,  digits = 2) %>%
#    kable_styling(full_width = F)
# 
# save_kable(vars_latex, "stockholm_vars.png")
```

ordinal
=======

``` r
dw_combine_imputations<- function( imputeddata,formula,reg_type,outcome_levels,imputations){
  
  
  model<-formula
  
  library(rms)
  
  
  Stat_sumc_SE<-NULL
  Stat_sumc_coef<-NULL
  Stat_sumc_var<-NULL
  Stat_psr2<-NULL
  for (i in 1:imputations){       
    datado<-complete(imputeddata,action=i)
    
    if(  reg_type==1){
      mod<-lm(formula(model),data=datado)
      preddys<-2  
    } else if ( reg_type==2){
      datado<-droplevels(datado)
      mod<-lrm(formula(model),data=datado)
      Stat_psr2<-c(Stat_psr2,mod$stats[10])
      preddys<-outcome_levels
    } else {
      print( "no such type yet supported")
    } 
    Stat_sumc_SE<-rbind(Stat_sumc_SE,sqrt(diag(vcov(mod))))
    Stat_sumc_var<-rbind(Stat_sumc_var,(diag(vcov(mod))))
    Stat_sumc_coef<-rbind(Stat_sumc_coef,mod$coef)
    
  }
  # ender<-length(attr(mod$terms,'predvars'))
  ender<-length(Stat_sumc_SE[1,])
  
  Stat_sumc_SE<-Stat_sumc_SE[,preddys:ender]
  Stat_sumc_coef<- Stat_sumc_coef[,preddys:ender]
  Stat_sumc_var<- Stat_sumc_var[,preddys:ender]
  
  
  
  m<-length(datado[,1])
  qmeans<-colMeans(Stat_sumc_coef)
  within_var<-colMeans(Stat_sumc_var)
  
  between_var<-(t(Stat_sumc_coef)-qmeans)^2
  bewteen_var_full<-apply(between_var,1,sum)*(1/(m-1))
  full_var<-within_var+bewteen_var_full*(1+1/m)
  full_SE<-sqrt(full_var)
  pvalues<-format(2*pt(-abs(qmeans/full_SE),df=m-1),digits=3,scientific=FALSE)
  
  if (  reg_type==1){
    print('last imputation')
    print(summary(mod))
    print('IMPORTANT ! The above is last imputation only ! ')
  } else if (  reg_type==2){
    print('Last imputaion')
    print(mod)
    print('IMPORTANT ! The above is last imputation only !')
    print('Pseudo-R2 by imputation ')
    print( Stat_psr2)
    print('mean Pseudo-R2')
    print( mean(Stat_psr2)) }
  
  print('Combined p-values')
  print(pvalues)
  
}
```

Stockholm MRI score

``` r
print(dw_combine_imputations(imputeddata = readRDS('imputed_data.rds'),
                             formula = formulas_stockholm[2],
                             reg_type = 2, outcome_levels = 2, 
                             imputations = 7))
```

    ## [1] "Last imputaion"
    ## Logistic Regression Model
    ##  
    ##  lrm(formula = formula(model), data = datado)
    ##  
    ##                        Model Likelihood     Discrimination    Rank Discrim.    
    ##                           Ratio Test           Indexes           Indexes       
    ##  Obs           351    LR chi2     135.66    R2       0.428    C       0.836    
    ##   unfavourable 183    d.f.            10    g        1.848    Dxy     0.672    
    ##   favourable   168    Pr(> chi2) <0.0001    gr       6.348    gamma   0.672    
    ##  max |deriv| 3e-09                          gp       0.337    tau-a   0.336    
    ##                                             Brier    0.165                     
    ##  
    ##                  Coef    S.E.   Wald Z Pr(>|Z|)
    ##  Intercept        1.4972 0.9043  1.66  0.0978  
    ##  Alder           -0.0451 0.0079 -5.73  <0.0001 
    ##  GCSKS            0.1273 0.0419  3.04  0.0024  
    ##  hypoxia=1       -0.5625 0.3289 -1.71  0.0872  
    ##  hypotension=0    1.1628 0.4036  2.88  0.0040  
    ##  Pupiller         0.3916 0.5326  0.74  0.4622  
    ##  Pupiller=2      -0.0700 0.9179 -0.08  0.9392  
    ##  stockholm       -0.2808 0.1273 -2.21  0.0274  
    ##  stockholm_mri   -0.1915 0.3707 -0.52  0.6055  
    ##  stockholm_mri=3 -1.4779 0.6292 -2.35  0.0188  
    ##  stockholm_mri=4 -0.9689 0.9911 -0.98  0.3282  
    ##  
    ## [1] "IMPORTANT ! The above is last imputation only !"
    ## [1] "Pseudo-R2 by imputation "
    ##        R2        R2        R2        R2        R2        R2        R2 
    ## 0.5030049 0.4798932 0.4752411 0.4388076 0.5089004 0.4905219 0.4276832 
    ## [1] "mean Pseudo-R2"
    ## [1] 0.4748646
    ## [1] "Combined p-values"
    ##           Alder           GCSKS       hypoxia=1   hypotension=0        Pupiller 
    ##     "0.0000001"     "0.0107441"     "0.2739787"     "0.0009432"     "0.5371685" 
    ##      Pupiller=2       stockholm   stockholm_mri stockholm_mri=3 stockholm_mri=4 
    ##     "0.9900172"     "0.0000314"     "0.6223121"     "0.0282624"     "0.3669397" 
    ##           Alder           GCSKS       hypoxia=1   hypotension=0        Pupiller 
    ##     "0.0000001"     "0.0107441"     "0.2739787"     "0.0009432"     "0.5371685" 
    ##      Pupiller=2       stockholm   stockholm_mri stockholm_mri=3 stockholm_mri=4 
    ##     "0.9900172"     "0.0000314"     "0.6223121"     "0.0282624"     "0.3669397"

Adams

``` r
print(dw_combine_imputations(imputeddata = readRDS('imputed_data.rds'),
                             formula = formulas_stockholm[1],
                             reg_type = 2, outcome_levels = 2, 
                             imputations = 7))
```

    ## [1] "Last imputaion"
    ## Logistic Regression Model
    ##  
    ##  lrm(formula = formula(model), data = datado)
    ##  
    ##                        Model Likelihood     Discrimination    Rank Discrim.    
    ##                           Ratio Test           Indexes           Indexes       
    ##  Obs           351    LR chi2     117.63    R2       0.380    C       0.815    
    ##   unfavourable 183    d.f.            10    g        1.690    Dxy     0.629    
    ##   favourable   168    Pr(> chi2) <0.0001    gr       5.419    gamma   0.629    
    ##  max |deriv| 3e-10                          gp       0.315    tau-a   0.315    
    ##                                             Brier    0.174                     
    ##  
    ##                Coef    S.E.   Wald Z Pr(>|Z|)
    ##  Intercept      0.8695 0.6981  1.25  0.2129  
    ##  Alder         -0.0424 0.0076 -5.55  <0.0001 
    ##  GCSKS          0.1347 0.0411  3.28  0.0010  
    ##  hypoxia=1     -0.4945 0.3168 -1.56  0.1185  
    ##  hypotension=0  1.2654 0.4118  3.07  0.0021  
    ##  Pupiller       0.3389 0.5050  0.67  0.5022  
    ##  Pupiller=2    -0.1401 0.8631 -0.16  0.8711  
    ##  stockholm     -0.3398 0.1259 -2.70  0.0070  
    ##  adams          1.1749 0.5596  2.10  0.0358  
    ##  adams=2       -2.7656 1.0363 -2.67  0.0076  
    ##  adams=3       -4.4454 1.5574 -2.85  0.0043  
    ##  
    ## [1] "IMPORTANT ! The above is last imputation only !"
    ## [1] "Pseudo-R2 by imputation "
    ##        R2        R2        R2        R2        R2        R2        R2 
    ## 0.4713654 0.4355970 0.4399829 0.4021710 0.4731195 0.4424800 0.3799156 
    ## [1] "mean Pseudo-R2"
    ## [1] 0.4349473
    ## [1] "Combined p-values"
    ##         Alder         GCSKS     hypoxia=1 hypotension=0      Pupiller 
    ## "0.000000215" "0.006373661" "0.268865218" "0.000616376" "0.602910180" 
    ##    Pupiller=2     stockholm         adams       adams=2       adams=3 
    ## "0.981456373" "0.000010429" "0.076079529" "0.017124694" "0.012380363" 
    ##         Alder         GCSKS     hypoxia=1 hypotension=0      Pupiller 
    ## "0.000000215" "0.006373661" "0.268865218" "0.000616376" "0.602910180" 
    ##    Pupiller=2     stockholm         adams       adams=2       adams=3 
    ## "0.981456373" "0.000010429" "0.076079529" "0.017124694" "0.012380363"

Abu Hamdeh

``` r
print(dw_combine_imputations(imputeddata = readRDS('imputed_data.rds'),
                             formula = formulas_stockholm[3],
                             reg_type = 2, outcome_levels = 2, 
                             imputations = 7))
```

    ## [1] "Last imputaion"
    ## Logistic Regression Model
    ##  
    ##  lrm(formula = formula(model), data = datado)
    ##  
    ##                        Model Likelihood     Discrimination    Rank Discrim.    
    ##                           Ratio Test           Indexes           Indexes       
    ##  Obs           351    LR chi2     121.90    R2       0.391    C       0.821    
    ##   unfavourable 183    d.f.            11    g        1.724    Dxy     0.641    
    ##   favourable   168    Pr(> chi2) <0.0001    gr       5.609    gamma   0.641    
    ##  max |deriv| 3e-10                          gp       0.321    tau-a   0.321    
    ##                                             Brier    0.172                     
    ##  
    ##                Coef    S.E.   Wald Z Pr(>|Z|)
    ##  Intercept      0.8408 0.7074  1.19  0.2346  
    ##  Alder         -0.0439 0.0078 -5.63  <0.0001 
    ##  GCSKS          0.1292 0.0414  3.12  0.0018  
    ##  hypoxia=1     -0.4304 0.3206 -1.34  0.1795  
    ##  hypotension=0  1.2897 0.4148  3.11  0.0019  
    ##  Pupiller       0.5069 0.5165  0.98  0.3264  
    ##  Pupiller=2    -0.3683 0.8829 -0.42  0.6765  
    ##  stockholm     -0.3321 0.1268 -2.62  0.0088  
    ##  hamdeh         1.1228 0.5592  2.01  0.0447  
    ##  hamdeh=2      -2.6570 1.0364 -2.56  0.0104  
    ##  hamdeh=3      -3.9710 1.5665 -2.53  0.0112  
    ##  hamdeh=4      -5.8815 2.1174 -2.78  0.0055  
    ##  
    ## [1] "IMPORTANT ! The above is last imputation only !"
    ## [1] "Pseudo-R2 by imputation "
    ##        R2        R2        R2        R2        R2        R2        R2 
    ## 0.4740625 0.4332985 0.4381281 0.4072318 0.4663109 0.4428310 0.3914405 
    ## [1] "mean Pseudo-R2"
    ## [1] 0.4361862
    ## [1] "Combined p-values"
    ##         Alder         GCSKS     hypoxia=1 hypotension=0      Pupiller 
    ## "0.000000145" "0.009514523" "0.339132331" "0.000652263" "0.494829937" 
    ##    Pupiller=2     stockholm        hamdeh      hamdeh=2      hamdeh=3 
    ## "0.958910291" "0.000011625" "0.118674527" "0.030358226" "0.039896792" 
    ##      hamdeh=4 
    ## "0.023513396" 
    ##         Alder         GCSKS     hypoxia=1 hypotension=0      Pupiller 
    ## "0.000000145" "0.009514523" "0.339132331" "0.000652263" "0.494829937" 
    ##    Pupiller=2     stockholm        hamdeh      hamdeh=2      hamdeh=3 
    ## "0.958910291" "0.000011625" "0.118674527" "0.030358226" "0.039896792" 
    ##      hamdeh=4 
    ## "0.023513396"

Firsching

``` r
print(dw_combine_imputations(imputeddata = readRDS('imputed_data.rds'),
                             formula = formulas_stockholm[4],
                             reg_type = 2, outcome_levels = 2, 
                             imputations = 7))
```

    ## [1] "Last imputaion"
    ## Logistic Regression Model
    ##  
    ##  lrm(formula = formula(model), data = datado)
    ##  
    ##                        Model Likelihood     Discrimination    Rank Discrim.    
    ##                           Ratio Test           Indexes           Indexes       
    ##  Obs           351    LR chi2     115.41    R2       0.374    C       0.814    
    ##   unfavourable 183    d.f.            11    g        1.635    Dxy     0.628    
    ##   favourable   168    Pr(> chi2) <0.0001    gr       5.129    gamma   0.628    
    ##  max |deriv| 8e-11                          gp       0.314    tau-a   0.314    
    ##                                             Brier    0.175                     
    ##  
    ##                Coef    S.E.   Wald Z Pr(>|Z|)
    ##  Intercept      1.1843 0.7067  1.68  0.0938  
    ##  Alder         -0.0417 0.0075 -5.56  <0.0001 
    ##  GCSKS          0.1356 0.0413  3.28  0.0010  
    ##  hypoxia=1     -0.5546 0.3166 -1.75  0.0798  
    ##  hypotension=0  1.0318 0.4104  2.51  0.0119  
    ##  Pupiller       0.2434 0.5100  0.48  0.6332  
    ##  Pupiller=2     0.1342 0.8753  0.15  0.8781  
    ##  stockholm     -0.3297 0.1250 -2.64  0.0084  
    ##  firsching     -0.2773 0.3855 -0.72  0.4720  
    ##  firsching=2   -0.1144 0.5916 -0.19  0.8466  
    ##  firsching=3   -1.1078 1.0148 -1.09  0.2750  
    ##  firsching=4   -0.2583 1.3191 -0.20  0.8448  
    ##  
    ## [1] "IMPORTANT ! The above is last imputation only !"
    ## [1] "Pseudo-R2 by imputation "
    ##        R2        R2        R2        R2        R2        R2        R2 
    ## 0.4673385 0.4378791 0.4477671 0.4053320 0.4707025 0.4542602 0.3738448 
    ## [1] "mean Pseudo-R2"
    ## [1] 0.436732
    ## [1] "Combined p-values"
    ##         Alder         GCSKS     hypoxia=1 hypotension=0      Pupiller 
    ## "0.000000159" "0.007524190" "0.289707005" "0.002145554" "0.781637460" 
    ##    Pupiller=2     stockholm     firsching   firsching=2   firsching=3 
    ## "0.726191849" "0.000005274" "0.331612756" "0.847045757" "0.336100193" 
    ##   firsching=4 
    ## "0.850126023" 
    ##         Alder         GCSKS     hypoxia=1 hypotension=0      Pupiller 
    ## "0.000000159" "0.007524190" "0.289707005" "0.002145554" "0.781637460" 
    ##    Pupiller=2     stockholm     firsching   firsching=2   firsching=3 
    ## "0.726191849" "0.000005274" "0.331612756" "0.847045757" "0.336100193" 
    ##   firsching=4 
    ## "0.850126023"

``` r
sessionInfo()
```

    ## R version 3.6.2 (2019-12-12)
    ## Platform: x86_64-w64-mingw32/x64 (64-bit)
    ## Running under: Windows 10 x64 (build 18363)
    ## 
    ## Matrix products: default
    ## 
    ## locale:
    ## [1] LC_COLLATE=English_Sweden.1252  LC_CTYPE=English_Sweden.1252   
    ## [3] LC_MONETARY=English_Sweden.1252 LC_NUMERIC=C                   
    ## [5] LC_TIME=English_Sweden.1252    
    ## 
    ## attached base packages:
    ## [1] parallel  stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] rms_5.1-4                       SparseM_1.78                   
    ##  [3] doParallel_1.0.15               funModeling_1.9.3              
    ##  [5] Hmisc_4.3-1                     Formula_1.2-3                  
    ##  [7] survival_3.1-8                  GA_3.2                         
    ##  [9] iterators_1.0.12                foreach_1.4.8                  
    ## [11] cowplot_1.0.0                   mltools_0.3.5                  
    ## [13] pROC_1.16.1                     glmnet_3.0-2                   
    ## [15] Matrix_1.2-18                   stringi_1.4.6                  
    ## [17] naniar_0.5.0                    pscl_1.5.2                     
    ## [19] mice_3.8.0                      lubridate_1.7.8                
    ## [21] corrr_0.4.1                     table1_1.1                     
    ## [23] caTools_1.18.0                  DescTools_0.99.32              
    ## [25] descr_1.1.4                     e1071_1.7-3                    
    ## [27] modelgrid_1.1.1.0               readxl_1.3.1                   
    ## [29] corrplot_0.84                   summarytools_0.9.5             
    ## [31] kableExtra_1.1.0                ggfortify_0.4.8                
    ## [33] MLmetrics_1.1.1                 caret_6.0-85                   
    ## [35] lattice_0.20-38                 AppliedPredictiveModeling_1.1-7
    ## [37] knitr_1.28                      forcats_0.5.0                  
    ## [39] stringr_1.4.0                   readr_1.3.1                    
    ## [41] tidyverse_1.3.0                 yardstick_0.0.5                
    ## [43] workflows_0.1.0                 tune_0.0.1                     
    ## [45] tibble_3.0.0                    rsample_0.0.5                  
    ## [47] tidyr_1.0.2                     recipes_0.1.9                  
    ## [49] purrr_0.3.3                     parsnip_0.0.5                  
    ## [51] infer_0.5.1                     ggplot2_3.3.0                  
    ## [53] dplyr_0.8.5                     dials_0.0.4                    
    ## [55] scales_1.1.0                    broom_0.5.5                    
    ## [57] tidymodels_0.1.0                data.table_1.12.8              
    ## [59] MASS_7.3-51.5                  
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] SnowballC_0.6.0      ModelMetrics_1.2.2.1 visdat_0.5.3        
    ##   [4] acepack_1.4.1        multcomp_1.4-12      dygraphs_1.1.1.6    
    ##   [7] rpart_4.1-15         inline_0.3.15        generics_0.0.2      
    ##  [10] GPfit_1.0-8          TH.data_1.0-10       callr_3.4.3         
    ##  [13] polspline_1.1.17     future_1.16.0        tokenizers_0.2.1    
    ##  [16] webshot_0.5.2        xml2_1.2.2           httpuv_1.5.2        
    ##  [19] StanHeaders_2.21.0-1 assertthat_0.2.1     gower_0.2.1         
    ##  [22] xfun_0.12            hms_0.5.3            bayesplot_1.7.1     
    ##  [25] evaluate_0.14        promises_1.1.0       fansi_0.4.1         
    ##  [28] dbplyr_1.4.2         igraph_1.2.4.2       DBI_1.1.0           
    ##  [31] htmlwidgets_1.5.1    stats4_3.6.2         ellipsis_0.3.0      
    ##  [34] crosstalk_1.0.0      backports_1.1.5      markdown_1.1        
    ##  [37] rapportools_1.0      moments_0.14         vctrs_0.2.4         
    ##  [40] quantreg_5.55        ROCR_1.0-7           entropy_1.2.1       
    ##  [43] withr_2.1.2          pryr_0.1.4           checkmate_2.0.0     
    ##  [46] xts_0.12-0           prettyunits_1.1.1    cluster_2.1.0       
    ##  [49] lazyeval_0.2.2       crayon_1.3.4         ellipse_0.4.1       
    ##  [52] pkgconfig_2.0.3      nlme_3.1-142         nnet_7.3-12         
    ##  [55] rlang_0.4.5          globals_0.12.5       lifecycle_0.2.0     
    ##  [58] miniUI_0.1.1.1       sandwich_2.5-1       MatrixModels_0.4-1  
    ##  [61] colourpicker_1.0     modelr_0.1.5         tidytext_0.2.2      
    ##  [64] cellranger_1.1.0     tcltk_3.6.2          matrixStats_0.55.0  
    ##  [67] loo_2.2.0            boot_1.3-23          zoo_1.8-7           
    ##  [70] reprex_0.3.0         base64enc_0.1-3      ggridges_0.5.2      
    ##  [73] processx_3.4.2       png_0.1-7            viridisLite_0.3.0   
    ##  [76] bitops_1.0-6         KernSmooth_2.23-16   pander_0.6.3        
    ##  [79] shape_1.4.4          jpeg_0.1-8.1         shinystan_2.5.0     
    ##  [82] magrittr_1.5         plyr_1.8.6           gplots_3.0.1.2      
    ##  [85] gdata_2.18.0         threejs_0.3.3        compiler_3.6.2      
    ##  [88] rstantools_2.0.0     RColorBrewer_1.1-2   plotrix_3.7-7       
    ##  [91] lme4_1.1-23          cli_2.0.2            DiceDesign_1.8-1    
    ##  [94] listenv_0.8.0        janeaustenr_0.1.5    ps_1.3.2            
    ##  [97] htmlTable_1.13.3     tidyselect_1.0.0     highr_0.8           
    ## [100] yaml_2.2.1           latticeExtra_0.6-29  grid_3.6.2          
    ## [103] tidypredict_0.4.5    tools_3.6.2          rstudioapi_0.11     
    ## [106] foreign_0.8-72       gridExtra_2.3        prodlim_2019.11.13  
    ## [109] rpart.plot_3.0.8     digest_0.6.25        shiny_1.4.0         
    ## [112] lava_1.6.6           Rcpp_1.0.3           later_1.0.0         
    ## [115] httr_1.4.1           rsconnect_0.8.16     colorspace_1.4-1    
    ## [118] rvest_0.3.5          fs_1.3.1             splines_3.6.2       
    ## [121] statmod_1.4.34       expm_0.999-4         shinythemes_1.1.2   
    ## [124] xtable_1.8-4         rstanarm_2.19.3      jsonlite_1.6.1      
    ## [127] nloptr_1.2.2.1       timeDate_3043.102    rstan_2.19.3        
    ## [130] CORElearn_1.54.2     ipred_0.9-9          R6_2.4.1            
    ## [133] lhs_1.0.1            pillar_1.4.3         htmltools_0.4.0     
    ## [136] mime_0.9             glue_1.3.1           fastmap_1.0.1       
    ## [139] minqa_1.2.4          DT_0.13              class_7.3-15        
    ## [142] codetools_0.2-16     pkgbuild_1.0.6       mvtnorm_1.1-0       
    ## [145] furrr_0.1.0          gtools_3.8.1         tidyposterior_0.0.2 
    ## [148] magick_2.3           shinyjs_1.1          rmarkdown_2.1       
    ## [151] munsell_0.5.0        haven_2.2.0          reshape2_1.4.3      
    ## [154] gtable_0.3.0
