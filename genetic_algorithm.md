-   [Aim](#aim)
-   [Packages](#packages)
    -   [Install packages](#install-packages)
    -   [Load Packages](#load-packages)
-   [Data preparation and splitting](#data-preparation-and-splitting)
    -   [Standardizing numeric
        variables](#standardizing-numeric-variables)
    -   [split](#split)
-   [Genetic algorithm](#genetic-algorithm)
    -   [Preparing functions that will be used later
        on](#preparing-functions-that-will-be-used-later-on)
    -   [Running the GA](#running-the-ga)
    -   [results GA](#results-ga)
-   [SVM](#svm)
    -   [tune control](#tune-control)
-   [Model performance](#model-performance)
    -   [GA SVM](#ga-svm)
    -   [Core SVM](#core-svm)
-   [ggplot graph](#ggplot-graph)

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

Data preparation and splitting
==============================

Load the data, which hass already undergone pre-processing and multiple
imputation (please examine the pdf document which contains the stepwise
selection procedure).

``` r
# Selecting a single imputed dataset

set.seed(123)

Data <- readRDS('imputed_data_ml.rds') %>% 
  mice::complete("long") %>%
  group_by(.imp) %>%
  nest() %>% 
  ungroup() %>% 
  sample_n(size =  1, replace = F) %>% 
  unnest() %>% 
 select(- c(.imp, .id)) 
```

check how many there are of each category

``` r
Data %>% 
  group_by(dich_gos) %>% 
  summarise(Frequency = n()) %>% 
  mutate(Percentage = round(100 * (Frequency/sum(Frequency))))  %>% 
  kable()
```

Examine the structure of the data

``` r
str(Data, list.len=ncol(Data))
```

Standardizing numeric variables
-------------------------------

``` r
pre_processing_recipe <- recipe(dich_gos ~ ., data = Data) %>% 
  step_normalize(all_numeric()) %>% 
  step_dummy(all_nominal(), - all_outcomes()) 

Data <- pre_processing_recipe %>% 
  prep(training = Data) %>% 
  juice(all_predictors()) %>% 
  mutate(dich_gos = Data$dich_gos)
```

split
-----

``` r
set.seed(123)

train_test_split <- initial_split(Data, strata = "dich_gos", prop = .66)

train_data <- training(train_test_split)
test_data <- testing(train_test_split)
```

Genetic algorithm
=================

Preparing functions that will be used later on
----------------------------------------------

The function “custom fitness” makes sure that the model tries to
optimize the area under the ROC curve, divided by the number of
variables used. thus, it’s a bit like the AIC in the sence that it
penalizes “overly-complicated” models. Although, I have since changed it
to just the AUC, as the penalty for adding variables proved to be far
too great.

The function “get\_roc\_metric” just helps evaluate the results in the
TRAINING dataset

``` r
cvIndex <- createMultiFolds(train_data$dich_gos, k = 10, times = 5)

ctrl <- trainControl(method = "repeatedcv", 
                     number = 10,
                     repeats = 5, 
                     classProbs = TRUE, 
                     summaryFunction = twoClassSummary, 
                     allowParallel = FALSE, 
                     index = cvIndex)

Dcv <- function(ind, x, y, cntrl) 
{
    library(caret)
    library(MASS)
    library(desirability)
    ind <- which(ind == 1)
    if (length(ind) == 0) return(0)
    
     mtry = sqrt(ncol(x[, ind]))
tunegrid = expand.grid(.mtry=round(mtry))
    
    out <- train(x[, ind], y, 
                 method = "rf", 
                 metric = "ROC", 
                 trControl = cntrl,
                 tuneGrid = tunegrid)
    
    
    rocVal <- caret:::getTrainPerf(out)[, "TrainROC"]
    dROC <- dMax(0.5, 1)
    dPreds <- dMin(11, ncol(x))
    ## Comnined the two with a geometirc mean
    allD <- dOverall(dROC, dPreds)
    ## Calculate the overall desirability value
    predict(allD, data.frame(ROC = rocVal, NumPred = length(ind)))
}

initialSmall <- function(object, ...)
{
    population <- sample(0:1,
                         replace = TRUE,
                         size = object@nBits * object@popSize,
                         prob = c(0.7, 0.3))
    population <- matrix(population,
                         nrow = object@popSize,
                         ncol = object@nBits)
    return(population)
}
```

preparing for the GA

``` r
data_x <- train_data %>% 
  select(-dich_gos)

param_nBits=ncol(data_x)
col_names=colnames(data_x)

y <- train_data %>% 
  pull(dich_gos) 
```

Running the GA
--------------

``` r
GA <- ga(type = "binary", 
           fitness = Dcv, 
           lower = 0, upper = 1, 
           run = 200,
           maxiter = 2000, 
           nBits = param_nBits, 
           names = col_names, 
           x = data_x, 
           y = y, 
           cntrl = ctrl, 
           population = initialSmall,
           keepBest = TRUE, 
           parallel = FALSE,
           monitor = FALSE,
           seed=123)

saveRDS(GA, "GA_11.rds")
```

results GA
----------

``` r
GA <- readRDS("GA_12.rds")

plot(GA)
```

![](genetic_algorithm_files/figure-markdown_github/unnamed-chunk-9-1.png)

``` r
only_letters <- function(x){      
 
  
  y <- strsplit(unlist(x), "[^a-zA-Z]+") 
  z <- y %>% map(~paste(., collapse = "_")) %>% simplify()
  return(z)}

chosen_vars_ga=col_names[GA@solution[1,]==1]

chosen_vars_ga %>% 
  tibble() %>% 
  rename("Variable" = ".") %>% 
  mutate(Variable = only_letters(Variable)) %>% 
  mutate(Variable = str_replace(Variable, "_X", "")) %>% 
  kable(caption = "Variables selected by the GA")
```

<table>
<caption>
Variables selected by the GA
</caption>
<thead>
<tr>
<th style="text-align:left;">
Variable
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Alder
</td>
</tr>
<tr>
<td style="text-align:left;">
GCSKS
</td>
</tr>
<tr>
<td style="text-align:left;">
hem\_cc\_splenium
</td>
</tr>
<tr>
<td style="text-align:left;">
hem\_bg\_bilateral
</td>
</tr>
<tr>
<td style="text-align:left;">
hem\_tha\_unilat
</td>
</tr>
<tr>
<td style="text-align:left;">
hem\_mes\_crus
</td>
</tr>
<tr>
<td style="text-align:left;">
hem\_pons\_ventral
</td>
</tr>
<tr>
<td style="text-align:left;">
diff\_mes\_tegmentum
</td>
</tr>
<tr>
<td style="text-align:left;">
diff\_pons\_unilat
</td>
</tr>
<tr>
<td style="text-align:left;">
flair\_bg\_unilat
</td>
</tr>
<tr>
<td style="text-align:left;">
flair\_mes\_bilateral
</td>
</tr>
<tr>
<td style="text-align:left;">
hypotension
</td>
</tr>
</tbody>
</table>

``` r
# core_model <- c("Alder", "GCSKS", "hypotension_X0", 
#                 "Pupiller_2", "Pupiller_1")

core_model <- c("Alder", "GCSKS", "hypotension_X0")

chosen_vars_ga <- setdiff(chosen_vars_ga, core_model)

chosen_vars_formula_ga <- paste("dich_gos","~", paste(chosen_vars_ga, 
                                                       collapse = " + "), 
                                 sep = " " ) %>% 
   paste("+", paste(core_model, collapse = " + "), sep = " " ) %>% 
   as.formula()
```

SVM
===

``` r
rec_ga <-recipe(chosen_vars_formula_ga, data = train_data)

core_formula <- paste("dich_gos","~", paste(core_model, 
                                                       collapse = " + "), 
                                 sep = " " ) %>% 
  as.formula()

rec_core <- recipe(core_formula, data = train_data) 
```

tune control
------------

In preparation for the repeated 10-fold cross validation

``` r
control_random <- trainControl(
 method = "repeatedcv",
 number = 10,
 repeats = 10,
 classProbs = TRUE, 
  summaryFunction = twoClassSummary, search = "random")

control_grid <- trainControl(
 method = "repeatedcv",
 number = 10,
 repeats = 10,
 classProbs = TRUE, 
  summaryFunction = twoClassSummary, search = "grid")
```

``` r
set.seed(123)

svm_core <- train(rec_core, data = train_data,
method = "svmRadial",
tuneLength = 12,
trControl = control_random,
metric = "ROC")

set.seed(123)

svm_ga <- train(rec_ga, data = train_data,
method = "svmRadial",
tuneLength = 12,
trControl = control_random,
  metric = "ROC")
```

``` r
ga_c <-  as_vector(svm_ga$finalModel@param[[1]])

ga_sigma <- as_vector(svm_ga$finalModel@kernelf@kpar[["sigma"]])

sigmalist_ga <- as_vector(seq(ga_sigma-0.1, ga_sigma+0.1,
                                by = .025))

costlist_ga <- as_vector(seq(ga_c-1, (ga_c+1), by = .25))

core_c <-  as_vector(svm_core$finalModel@param[[1]])

core_sigma <- as_vector(svm_core$finalModel@kernelf@kpar[["sigma"]])

sigmalist_core <- as_vector(seq(core_sigma-0.1, core_sigma+0.1,
                                by = .025))

costlist_core <- as_vector(seq(core_c-1, (core_c+1), by = .25))
```

``` r
grid_ga <- expand.grid(sigma = sigmalist_ga,
  C = costlist_ga)

grid_core <- expand.grid(sigma = sigmalist_core,
  C = costlist_core)

set.seed(123)

svm_core <- train(rec_core, data = train_data,
                    method = "svmRadial",
                    tuneGrid = grid_core,
                    trControl = control_grid,
                      metric = "ROC")

saveRDS(svm_core, file = "svm_core_12.rds")

set.seed(123)

svm_ga <- train(rec_ga, data = train_data,
                    method = "svmRadial",
                    tuneGrid = grid_core,
                    trControl = control_grid,
                      metric = "ROC")

saveRDS(svm_ga, file = "svm_ga_12.rds")
```

Model performance
=================

read previously saved rds files

``` r
svm_ga <- readRDS("svm_ga_12.rds")

svm_core <- readRDS("svm_core_12.rds")
```

GA SVM
------

``` r
testing_results <- test_data %>%
        mutate("ga" = predict(svm_ga, test_data),
               "svm_core" = predict(svm_core, test_data))
  
ga <- conf_mat(testing_results, truth = dich_gos, estimate = ga)

autoplot(ga, type = "heatmap") +
  ggtitle("confusion matrix, SVM with GA") +
  theme(plot.title = element_text(hjust = 0.5))
```

![](genetic_algorithm_files/figure-markdown_github/unnamed-chunk-18-1.png)

``` r
ga <- (summary(ga)) %>% 
  select("Metric" = .metric, "Estimate" = .estimate) %>% 
  mutate(Metric = str_replace(string=Metric, pattern='f_meas', replacement='F-score')) %>% 
 mutate(Metric = str_replace(string=Metric, pattern='j_index', replacement="Youden's index")) %>% 
  mutate(Metric = str_replace(string=Metric, pattern='kap', 
                              replacement="kappa")) %>% 
  filter(Metric %in% c("specificity", "sensitivity", "F-score", 
                       "Youden's index", "accuracy", "kappa", "precision",                              "recall"))

kable(ga,digits = 3, caption = "SVM with GA")
```

<table>
<caption>
SVM with GA
</caption>
<thead>
<tr>
<th style="text-align:left;">
Metric
</th>
<th style="text-align:right;">
Estimate
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
accuracy
</td>
<td style="text-align:right;">
0.748
</td>
</tr>
<tr>
<td style="text-align:left;">
kappa
</td>
<td style="text-align:right;">
0.496
</td>
</tr>
<tr>
<td style="text-align:left;">
Youden’s index
</td>
<td style="text-align:right;">
0.496
</td>
</tr>
<tr>
<td style="text-align:left;">
precision
</td>
<td style="text-align:right;">
0.767
</td>
</tr>
<tr>
<td style="text-align:left;">
recall
</td>
<td style="text-align:right;">
0.742
</td>
</tr>
<tr>
<td style="text-align:left;">
F-score
</td>
<td style="text-align:right;">
0.754
</td>
</tr>
</tbody>
</table>

ROC

``` r
prob_ga <- predict(svm_ga, newdata = test_data %>% select(-dich_gos), 
                      type = "prob")
roc_ga <- roc(response = test_data$dich_gos, predictor = prob_ga[, "unfavourable"])

roc_ga_graph <- ggroc(roc_ga) + theme_minimal() + ggtitle("ROC: SVM with GA") + 
    geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="grey", linetype="dashed") +  annotate("text", x = .9, y = .2, 
       label = paste("AUC =", round(roc_ga$auc, 2))) + xlab("False Positive Rate (1-Specificity)") + ylab("True Positive Rate (Sensitivity)") + ggtitle("ROC: SVM with GA")
```

Core SVM
--------

``` r
core_svm <- conf_mat(testing_results, truth = dich_gos, estimate = svm_core)

autoplot(core_svm, type = "heatmap") +
  ggtitle("confusion matrix, SVM with the core variables") +
  theme(plot.title = element_text(hjust = 0.5))
```

![](genetic_algorithm_files/figure-markdown_github/unnamed-chunk-20-1.png)

``` r
core_svm <- (summary(core_svm)) %>% 
  select("Metric" = .metric, "Estimate" = .estimate) %>% 
  mutate(Metric = str_replace(string=Metric, pattern='f_meas', replacement='F-score')) %>% 
 mutate(Metric = str_replace(string=Metric, pattern='j_index', replacement="Youden's index")) %>% 
  mutate(Metric = str_replace(string=Metric, pattern='kap', 
                              replacement="kappa")) %>% 
  filter(Metric %in% c("specificity", "sensitivity", "F-score", 
                       "Youden's index", "accuracy", "kappa", "precision",                              "recall"))

kable(core_svm,digits = 3, caption = "SVM with the core variables")
```

<table>
<caption>
SVM with the core variables
</caption>
<thead>
<tr>
<th style="text-align:left;">
Metric
</th>
<th style="text-align:right;">
Estimate
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
accuracy
</td>
<td style="text-align:right;">
0.647
</td>
</tr>
<tr>
<td style="text-align:left;">
kappa
</td>
<td style="text-align:right;">
0.305
</td>
</tr>
<tr>
<td style="text-align:left;">
Youden’s index
</td>
<td style="text-align:right;">
0.310
</td>
</tr>
<tr>
<td style="text-align:left;">
precision
</td>
<td style="text-align:right;">
0.763
</td>
</tr>
<tr>
<td style="text-align:left;">
recall
</td>
<td style="text-align:right;">
0.468
</td>
</tr>
<tr>
<td style="text-align:left;">
F-score
</td>
<td style="text-align:right;">
0.580
</td>
</tr>
</tbody>
</table>

ROC

``` r
prob_core_svm <- predict(svm_core, newdata = test_data %>% select(-dich_gos), 
                      type = "prob")
roc_core_svm <- roc(response = test_data$dich_gos, predictor = prob_core_svm[, "unfavourable"])

roc_core_svm_graph <- ggroc(roc_core_svm) + theme_minimal() + ggtitle("ROC: SVM with the core variables") + 
    geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="grey", linetype="dashed") +  annotate("text", x = .9, y = .2, 
       label = paste("AUC =", round(roc_core_svm$auc, 2))) + xlab("False Positive Rate (1-Specificity)") + ylab("True Positive Rate (Sensitivity)")
```

``` r
print(roc_plots <- plot_grid(roc_core_svm_graph,
                       roc_ga_graph, 
                       rows = 1, cols = 2))
```

![](genetic_algorithm_files/figure-markdown_github/unnamed-chunk-22-1.png)

ggplot graph
============

``` r
ggroc(roc_core_svm, alpha = 0.7, colour = "Blue",
linetype = 1, size = 2) +
theme_minimal() +
geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="grey", linetype="dashed") + annotate("text", x= 0.68, y = 0.5, label = paste("AUC =", round(auc(roc_core_svm), digits = 2))) + xlab("False Positive Rate (1-Specificity)") + ylab("True Positive Rate (Sensitivity)")
```

![](genetic_algorithm_files/figure-markdown_github/unnamed-chunk-23-1.png)

``` r
# roc_list <- roc(test_data$dich_gos ~ prob_ga[, "unfavourable"] + prob_core_svm[, "unfavourable"])

results <- test_data %>%
        mutate("ga" = predict(svm_ga, test_data, type = "prob")[,"unfavourable"],
    "svm_core" = predict(svm_core, test_data, type = "prob")[,"unfavourable"]) %>% 
  select(dich_gos, "Clinical variables + MRI characteristics" = ga, "Clinical variables" = svm_core) %>% 
  gather(key = "Model", value = "fitted", - dich_gos) 


results %>% 
  group_by(Model) %>% 
  roc_curve(truth = dich_gos, fitted) %>% 
  ggplot(
    aes(
      x = 1 - specificity, 
      y = sensitivity, 
      color = Model
    )
  ) + # plot with 2 ROC curves for each model
  geom_line(size = 1.1) +
  geom_abline(slope = 1, intercept = 0, size = 0.4) +
  scale_color_manual(values = c("#48466D", "#3D84A8")) +
  coord_fixed() +
  theme_cowplot() 
```

![](genetic_algorithm_files/figure-markdown_github/unnamed-chunk-24-1.png)

``` r
# ggsave("auc_plot_initial.png")
```

``` r
results %>% 
  group_by(Model) %>%
roc_auc(truth = dich_gos, fitted) %>% 
  select(Model, "AUC" = .estimate) %>%
  kable(digits = 2, format = "latex", booktabs = T, linesep = "", align = 'c') %>% 
  row_spec(0,  bold = T, color = "Black") %>%
  row_spec(2, bold = T, color = "#3D84A8") %>% 
  row_spec(1, bold = T, color = "#48466D") %>% 
  save_kable(., "auc_table.png")
```

``` r
p <- results %>% 
  group_by(Model) %>% 
  roc_curve(truth = dich_gos, fitted) %>% 
  ggplot(
    aes(
      x = 1 - specificity, 
      y = sensitivity, 
      color = Model
    )
  ) + # plot with 2 ROC curves for each model
  geom_line(size = 1.1) +
  geom_abline(slope = 1, intercept = 0, size = 0.4) +
  scale_color_manual(values = c("#48466D", "#3D84A8")) +
  coord_fixed() +
  theme_cowplot() 
  
library(png)

auc_table <- readPNG("auc_table.png")

roc_plot <- ggdraw(p) + 
  draw_image(auc_table, scale = 0.45, x = 0.1, y = - 0.1)  

ggsave2("roc_plot.png", plot = roc_plot,
        dpi = 300)
```
