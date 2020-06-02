-   [Aim](#aim)
-   [Packages](#packages)
    -   [Install packages](#install-packages)
    -   [Load Packages](#load-packages)
-   [Data preparation](#data-preparation)
    -   [Load the data](#load-the-data)
    -   [Calculate the time difference between the MRI-scan and the
        trauma](#calculate-the-time-difference-between-the-mri-scan-and-the-trauma)
    -   [Structure](#structure)
    -   [Studying the data more
        thoroughly](#studying-the-data-more-thoroughly)
    -   [continue pre-processing the
        data](#continue-pre-processing-the-data)
    -   [David Nelson’s pupil fix](#david-nelsons-pupil-fix)
    -   [mutate MRI scores as ordered
        “factors”](#mutate-mri-scores-as-ordered-factors)
    -   [Data summary](#data-summary)
-   [Missing](#missing)
-   [Univariate logistic regression](#univariate-logistic-regression)
    -   [A function to calculate the pseudo
        R-squared](#a-function-to-calculate-the-pseudo-r-squared)
    -   [Fitting univariate logistic
        regression](#fitting-univariate-logistic-regression)
-   [Multiple imputation](#multiple-imputation)
-   [Fitting glm models to imputed
    datasets](#fitting-glm-models-to-imputed-datasets)
    -   [Fitting the GLM models and extracting their
        coefficients](#fitting-the-glm-models-and-extracting-their-coefficients)
    -   [results](#results)
    -   [Results for each individual
        variable..](#results-for-each-individual-variable..)
    -   [ROC](#roc)
-   [Session information](#session-information)

Aim
===

The aim of this script is to run the data analysis for the MRI
variables. In this case, a combination of univariate logistic regression
and stepwise variable selection were used to select variables for a
multivariate logistic regression. The same procedure was also used for
ordinal logistic regression.

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
install.packages("CaTools")
install.packages("table1", dependencies = T)
install.packages("mice")
install.packages("pscl")
install.packages("naniar")
install.packages("cowplot")
install.packages("glmnet")
install.packages("pROC")
install.packages("rcompanion")
install.packages("DescTools")
```

Load Packages
-------------

``` r
library(MASS)
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
library(caTools)
library(table1)
library(corrr)
library(lubridate)
library(mice)
library(pscl)
library(naniar)
library(cowplot)
library(glmnet)
library(pROC)
library(DescTools)
library(rcompanion)

options(knitr.kable.NA = '')
```

Data preparation
================

Load the data
-------------

``` r
Data <- read_excel('dai_results_mri_hem.xls', na = "NA") 
```

Calculate the time difference between the MRI-scan and the trauma
-----------------------------------------------------------------

calculate the difference, in days, between the date of the trauma and
the date of the MRI-scan

``` r
Data <- Data %>% 
   mutate(mri_date = as_date(mri_date)) %>% 
   mutate(trauma_date = as_date(trauma_date)) %>% 
   select(-GCSOpl)

Data <- Data %>% 
   mutate(mri_time = difftime(time1 = mri_date, time2 = trauma_date, 
                                     units = "days"))
```

Structure
---------

keep only the first 406 rows, corresponding to all patients admitted
prior to December 31 2018, and only patients who have had their
MRI-scans performed within the first 28 days after their trauma and were
above the age 15

``` r
Data <- Data %>%
  filter(id < 1593) %>% 
  drop_na(mri_date) %>%
   select(-contains("_date")) %>% 
   filter(mri_time <= 28) %>% 
  filter(Alder >= 15)
```

Examine the structure of the data

``` r
str(Data, list.len=ncol(Data))
```

change the class of certain variables, based on that which was seen
using the str() function above.

``` r
Data <- Data %>% 
   mutate(FinalGOS = factor(FinalGOS, ordered = TRUE)) %>% 
   mutate(Pupreak = as.factor(Pupreak)) %>% 
   mutate(rotterdam = factor(rotterdam, ordered = TRUE)) %>%
   mutate(marshall = as.factor(marshall)) %>% 
   mutate(adams = factor(adams, ordered = TRUE)) %>% 
   select(-c(GOSB,GOSC, infarct_bg))
```

Blood pressure and oxygen saturation

``` r
# replace "0" with "NA"

Data <- Data %>% 
  replace_with_na(replace = list(SaO2Opl = 0,
                             BltrOpl = 0, FinalGOS = 0))

# Dichotomize the saturation into hypoxia (1) or no hypoxia (0) 
# Also dichotomize the bloodpressure into hypotension (1) 
# or no hypotension (0)

Data <- Data%>% 
   mutate(hypoxia = as_factor(ifelse(SaO2Opl < 90, "1", "0"))) %>% 
   mutate(hypotension = as_factor(ifelse(BltrOpl < 90, "1", "0")))

# removing the variables for prehospital oxygen saturation and prehospital blood pressure

Data <- Data %>% 
   select(-c(SaO2Opl, BltrOpl))
```

Create a dichotomized GOS variable and a mortality variable

``` r
Data <- Data %>% 
   mutate(dich_gos = as_factor(ifelse(FinalGOS == 4 | FinalGOS == 5, 
                                      "favourable", 
                                      "unfavourable"))) %>% 
   mutate(mortality = as_factor(ifelse(FinalGOS == 1, "Dead", 
                                      "Alive"))) 
```

Studying the data more thoroughly
---------------------------------

dich\_gos

``` r
Data %>% 
   select(dich_gos,hem_subcort_unilat:infarct_cerebellum) %>% 
   select(-contains("pons")) %>% 
   select(-contains("infarct")) %>%
   select(-contains("mes")) %>% 
   group_by(dich_gos) %>% 
   summarize_if(is.numeric, funs(sum(., na.rm = TRUE))) %>% 
   gather(key = "Variable", value = "Frequency", -dich_gos) %>% 
   na.omit() %>%
   ungroup() %>% 
   group_by(Variable) %>% 
   spread(dich_gos, Frequency) %>% 
   mutate(Total = unfavourable+favourable) %>% 
   mutate("Percent (unfavourable)" = round(unfavourable/Total
          *100, digits = 0)) %>% 
   arrange(desc(unfavourable)) %>% 
   kable(caption = "Percent of unfavourable outcomes following lesions outside of the brainstem")
```

<table>
<caption>
Percent of unfavourable outcomes following lesions outside of the
brainstem
</caption>
<thead>
<tr>
<th style="text-align:left;">
Variable
</th>
<th style="text-align:right;">
unfavourable
</th>
<th style="text-align:right;">
favourable
</th>
<th style="text-align:right;">
Total
</th>
<th style="text-align:right;">
Percent (unfavourable)
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
flair\_cc\_splenium
</td>
<td style="text-align:right;">
75
</td>
<td style="text-align:right;">
36
</td>
<td style="text-align:right;">
111
</td>
<td style="text-align:right;">
68
</td>
</tr>
<tr>
<td style="text-align:left;">
hem\_subcort\_bilateral
</td>
<td style="text-align:right;">
62
</td>
<td style="text-align:right;">
42
</td>
<td style="text-align:right;">
104
</td>
<td style="text-align:right;">
60
</td>
</tr>
<tr>
<td style="text-align:left;">
hem\_cc\_splenium
</td>
<td style="text-align:right;">
60
</td>
<td style="text-align:right;">
27
</td>
<td style="text-align:right;">
87
</td>
<td style="text-align:right;">
69
</td>
</tr>
<tr>
<td style="text-align:left;">
diff\_cc\_splenium
</td>
<td style="text-align:right;">
52
</td>
<td style="text-align:right;">
30
</td>
<td style="text-align:right;">
82
</td>
<td style="text-align:right;">
63
</td>
</tr>
<tr>
<td style="text-align:left;">
hem\_cc\_corpus
</td>
<td style="text-align:right;">
46
</td>
<td style="text-align:right;">
20
</td>
<td style="text-align:right;">
66
</td>
<td style="text-align:right;">
70
</td>
</tr>
<tr>
<td style="text-align:left;">
flair\_subcort\_bilateral
</td>
<td style="text-align:right;">
41
</td>
<td style="text-align:right;">
16
</td>
<td style="text-align:right;">
57
</td>
<td style="text-align:right;">
72
</td>
</tr>
<tr>
<td style="text-align:left;">
flair\_bg\_unilat
</td>
<td style="text-align:right;">
40
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
47
</td>
<td style="text-align:right;">
85
</td>
</tr>
<tr>
<td style="text-align:left;">
flair\_cc\_corpus
</td>
<td style="text-align:right;">
38
</td>
<td style="text-align:right;">
19
</td>
<td style="text-align:right;">
57
</td>
<td style="text-align:right;">
67
</td>
</tr>
<tr>
<td style="text-align:left;">
hem\_bg\_unilat
</td>
<td style="text-align:right;">
32
</td>
<td style="text-align:right;">
26
</td>
<td style="text-align:right;">
58
</td>
<td style="text-align:right;">
55
</td>
</tr>
<tr>
<td style="text-align:left;">
hem\_subcort\_unilat
</td>
<td style="text-align:right;">
31
</td>
<td style="text-align:right;">
36
</td>
<td style="text-align:right;">
67
</td>
<td style="text-align:right;">
46
</td>
</tr>
<tr>
<td style="text-align:left;">
flair\_tha\_unilat
</td>
<td style="text-align:right;">
28
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
39
</td>
<td style="text-align:right;">
72
</td>
</tr>
<tr>
<td style="text-align:left;">
hem\_tha\_unilat
</td>
<td style="text-align:right;">
26
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
38
</td>
<td style="text-align:right;">
68
</td>
</tr>
<tr>
<td style="text-align:left;">
flair\_subcort\_unilat
</td>
<td style="text-align:right;">
25
</td>
<td style="text-align:right;">
26
</td>
<td style="text-align:right;">
51
</td>
<td style="text-align:right;">
49
</td>
</tr>
<tr>
<td style="text-align:left;">
hem\_tha\_bilateral
</td>
<td style="text-align:right;">
25
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
28
</td>
<td style="text-align:right;">
89
</td>
</tr>
<tr>
<td style="text-align:left;">
hem\_cc\_genu
</td>
<td style="text-align:right;">
24
</td>
<td style="text-align:right;">
18
</td>
<td style="text-align:right;">
42
</td>
<td style="text-align:right;">
57
</td>
</tr>
<tr>
<td style="text-align:left;">
diff\_bg\_unilat
</td>
<td style="text-align:right;">
22
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
24
</td>
<td style="text-align:right;">
92
</td>
</tr>
<tr>
<td style="text-align:left;">
flair\_cc\_genu
</td>
<td style="text-align:right;">
21
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
27
</td>
<td style="text-align:right;">
78
</td>
</tr>
<tr>
<td style="text-align:left;">
diff\_cc\_corpus
</td>
<td style="text-align:right;">
19
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
23
</td>
<td style="text-align:right;">
83
</td>
</tr>
<tr>
<td style="text-align:left;">
diff\_cc\_genu
</td>
<td style="text-align:right;">
18
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
23
</td>
<td style="text-align:right;">
78
</td>
</tr>
<tr>
<td style="text-align:left;">
diff\_tha\_unilat
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
18
</td>
<td style="text-align:right;">
83
</td>
</tr>
<tr>
<td style="text-align:left;">
hem\_bg\_bilateral
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
19
</td>
<td style="text-align:right;">
79
</td>
</tr>
<tr>
<td style="text-align:left;">
flair\_tha\_bilateral
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
16
</td>
<td style="text-align:right;">
88
</td>
</tr>
<tr>
<td style="text-align:left;">
flair\_bg\_bilateral
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
92
</td>
</tr>
<tr>
<td style="text-align:left;">
diff\_subcort\_unilat
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
16
</td>
<td style="text-align:right;">
62
</td>
</tr>
<tr>
<td style="text-align:left;">
diff\_subcort\_bilateral
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
62
</td>
</tr>
<tr>
<td style="text-align:left;">
diff\_tha\_bilateral
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
100
</td>
</tr>
<tr>
<td style="text-align:left;">
diff\_bg\_bilateral
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
67
</td>
</tr>
</tbody>
</table>

``` r
Data %>% 
   select(dich_gos,hem_subcort_unilat:infarct_cerebellum) %>% 
   select(contains("pons"), dich_gos) %>% 
   group_by(dich_gos) %>% 
   summarize_if(is.numeric, funs(sum(., na.rm = TRUE))) %>% 
   gather(key = "Variable", value = "Frequency", -dich_gos) %>% 
   na.omit() %>%
   ungroup() %>% 
   group_by(Variable) %>% 
   spread(dich_gos, Frequency) %>% 
   mutate(Total = unfavourable+favourable) %>% 
   mutate("Percent (unfavourable)" = round(unfavourable/Total
          *100, digits = 0)) %>% 
   arrange(desc(unfavourable)) %>% 
   kable(caption = "Percent of unfavourable outcomes following lesions in the pons")
```

<table>
<caption>
Percent of unfavourable outcomes following lesions in the pons
</caption>
<thead>
<tr>
<th style="text-align:left;">
Variable
</th>
<th style="text-align:right;">
unfavourable
</th>
<th style="text-align:right;">
favourable
</th>
<th style="text-align:right;">
Total
</th>
<th style="text-align:right;">
Percent (unfavourable)
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
hem\_pons\_dorsal
</td>
<td style="text-align:right;">
27
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
35
</td>
<td style="text-align:right;">
77
</td>
</tr>
<tr>
<td style="text-align:left;">
flair\_pons\_dorsal
</td>
<td style="text-align:right;">
25
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
33
</td>
<td style="text-align:right;">
76
</td>
</tr>
<tr>
<td style="text-align:left;">
hem\_pons\_ventral
</td>
<td style="text-align:right;">
25
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
26
</td>
<td style="text-align:right;">
96
</td>
</tr>
<tr>
<td style="text-align:left;">
flair\_pons\_unilat
</td>
<td style="text-align:right;">
22
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
33
</td>
<td style="text-align:right;">
67
</td>
</tr>
<tr>
<td style="text-align:left;">
flair\_pons\_ventral
</td>
<td style="text-align:right;">
21
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
25
</td>
<td style="text-align:right;">
84
</td>
</tr>
<tr>
<td style="text-align:left;">
hem\_pons\_bilateral
</td>
<td style="text-align:right;">
21
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
22
</td>
<td style="text-align:right;">
95
</td>
</tr>
<tr>
<td style="text-align:left;">
hem\_pons\_unilat
</td>
<td style="text-align:right;">
17
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
25
</td>
<td style="text-align:right;">
68
</td>
</tr>
<tr>
<td style="text-align:left;">
diff\_pons\_unilat
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
17
</td>
<td style="text-align:right;">
88
</td>
</tr>
<tr>
<td style="text-align:left;">
flair\_pons\_bilateral
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:right;">
100
</td>
</tr>
<tr>
<td style="text-align:left;">
diff\_pons\_dorsal
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
85
</td>
</tr>
<tr>
<td style="text-align:left;">
diff\_pons\_ventral
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
100
</td>
</tr>
<tr>
<td style="text-align:left;">
diff\_pons\_bilateral
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
100
</td>
</tr>
</tbody>
</table>

``` r
Data %>% 
   select(dich_gos,hem_subcort_unilat:infarct_cerebellum) %>% 
   select(contains("mes"), dich_gos) %>% 
   group_by(dich_gos) %>% 
   summarize_if(is.numeric, funs(sum(., na.rm = TRUE))) %>% 
   gather(key = "Variable", value = "Frequency", -dich_gos) %>% 
   na.omit() %>%
   ungroup() %>% 
   group_by(Variable) %>% 
   spread(dich_gos, Frequency) %>% 
   mutate(Total = unfavourable+favourable) %>% 
   mutate("Percent (unfavourable)" = round(unfavourable/Total
          *100, digits = 0)) %>%
   arrange(desc(unfavourable)) %>% 
   kable(caption = "Percent of unfavourable outcomes following lesions in the midbrain")
```

<table>
<caption>
Percent of unfavourable outcomes following lesions in the midbrain
</caption>
<thead>
<tr>
<th style="text-align:left;">
Variable
</th>
<th style="text-align:right;">
unfavourable
</th>
<th style="text-align:right;">
favourable
</th>
<th style="text-align:right;">
Total
</th>
<th style="text-align:right;">
Percent (unfavourable)
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
hem\_mes\_tegmentum
</td>
<td style="text-align:right;">
47
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
57
</td>
<td style="text-align:right;">
82
</td>
</tr>
<tr>
<td style="text-align:left;">
flair\_mes\_tegmentum
</td>
<td style="text-align:right;">
46
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
57
</td>
<td style="text-align:right;">
81
</td>
</tr>
<tr>
<td style="text-align:left;">
flair\_mes\_unilat
</td>
<td style="text-align:right;">
38
</td>
<td style="text-align:right;">
17
</td>
<td style="text-align:right;">
55
</td>
<td style="text-align:right;">
69
</td>
</tr>
<tr>
<td style="text-align:left;">
hem\_mes\_bilateral
</td>
<td style="text-align:right;">
33
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
38
</td>
<td style="text-align:right;">
87
</td>
</tr>
<tr>
<td style="text-align:left;">
diff\_mes\_unilat
</td>
<td style="text-align:right;">
29
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
34
</td>
<td style="text-align:right;">
85
</td>
</tr>
<tr>
<td style="text-align:left;">
hem\_mes\_unilat
</td>
<td style="text-align:right;">
29
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
42
</td>
<td style="text-align:right;">
69
</td>
</tr>
<tr>
<td style="text-align:left;">
flair\_mes\_crus
</td>
<td style="text-align:right;">
27
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
32
</td>
<td style="text-align:right;">
84
</td>
</tr>
<tr>
<td style="text-align:left;">
diff\_mes\_tegmentum
</td>
<td style="text-align:right;">
23
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
23
</td>
<td style="text-align:right;">
100
</td>
</tr>
<tr>
<td style="text-align:left;">
flair\_mes\_bilateral
</td>
<td style="text-align:right;">
23
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
24
</td>
<td style="text-align:right;">
96
</td>
</tr>
<tr>
<td style="text-align:left;">
hem\_mes\_crus
</td>
<td style="text-align:right;">
21
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
29
</td>
<td style="text-align:right;">
72
</td>
</tr>
<tr>
<td style="text-align:left;">
flair\_mes\_tectum
</td>
<td style="text-align:right;">
18
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
21
</td>
<td style="text-align:right;">
86
</td>
</tr>
<tr>
<td style="text-align:left;">
diff\_mes\_crus
</td>
<td style="text-align:right;">
16
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
18
</td>
<td style="text-align:right;">
89
</td>
</tr>
<tr>
<td style="text-align:left;">
diff\_mes\_bilateral
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
100
</td>
</tr>
<tr>
<td style="text-align:left;">
hem\_mes\_tectum
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:right;">
67
</td>
</tr>
<tr>
<td style="text-align:left;">
diff\_mes\_tectum
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
89
</td>
</tr>
</tbody>
</table>

``` r
Data %>% 
   select(dich_gos,hem_subcort_unilat:infarct_cerebellum) %>% 
   select(contains("infarct"), dich_gos) %>% 
   group_by(dich_gos) %>% 
   summarize_if(is.numeric, funs(sum(., na.rm = TRUE))) %>% 
   gather(key = "Variable", value = "Frequency", -dich_gos) %>% 
   na.omit() %>%
   ungroup() %>% 
   group_by(Variable) %>% 
   spread(dich_gos, Frequency) %>% 
   mutate(Total = unfavourable+favourable) %>% 
   mutate("Percent (unfavourable)" = round(unfavourable/Total
          *100, digits = 0)) %>% 
   arrange(desc(unfavourable)) %>% 
   kable(caption = "Percent of unfavourable outcomes following infarcts")
```

<table>
<caption>
Percent of unfavourable outcomes following infarcts
</caption>
<thead>
<tr>
<th style="text-align:left;">
Variable
</th>
<th style="text-align:right;">
unfavourable
</th>
<th style="text-align:right;">
favourable
</th>
<th style="text-align:right;">
Total
</th>
<th style="text-align:right;">
Percent (unfavourable)
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
infarct\_posterior
</td>
<td style="text-align:right;">
31
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
37
</td>
<td style="text-align:right;">
84
</td>
</tr>
<tr>
<td style="text-align:left;">
infarct\_media
</td>
<td style="text-align:right;">
26
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
31
</td>
<td style="text-align:right;">
84
</td>
</tr>
<tr>
<td style="text-align:left;">
infarct\_anterior
</td>
<td style="text-align:right;">
22
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
23
</td>
<td style="text-align:right;">
96
</td>
</tr>
<tr>
<td style="text-align:left;">
infarct\_cerebellum
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
100
</td>
</tr>
</tbody>
</table>

mortality

``` r
Data %>% 
   select(mortality,hem_subcort_unilat:infarct_cerebellum) %>% 
   select(-contains("pons")) %>% 
   select(-contains("infarct")) %>%
   select(-contains("mes")) %>% 
   group_by(mortality) %>% 
   summarize_if(is.numeric, funs(sum(., na.rm = TRUE))) %>% 
   gather(key = "Variable", value = "Frequency", -mortality) %>% 
   na.omit() %>%
   ungroup() %>% 
   group_by(Variable) %>% 
   spread(mortality, Frequency) %>% 
   mutate(Total = Dead+Alive) %>% 
   mutate("Percent (Dead)" = round(Dead/Total
          *100, digits = 0)) %>% 
   arrange(desc(Dead)) %>% 
   kable(caption = "Percent mortalities following lesions outside of the brainstem")
```

<table>
<caption>
Percent mortalities following lesions outside of the brainstem
</caption>
<thead>
<tr>
<th style="text-align:left;">
Variable
</th>
<th style="text-align:right;">
Alive
</th>
<th style="text-align:right;">
Dead
</th>
<th style="text-align:right;">
Total
</th>
<th style="text-align:right;">
Percent (Dead)
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
flair\_cc\_splenium
</td>
<td style="text-align:right;">
95
</td>
<td style="text-align:right;">
16
</td>
<td style="text-align:right;">
111
</td>
<td style="text-align:right;">
14
</td>
</tr>
<tr>
<td style="text-align:left;">
diff\_cc\_splenium
</td>
<td style="text-align:right;">
68
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:right;">
82
</td>
<td style="text-align:right;">
17
</td>
</tr>
<tr>
<td style="text-align:left;">
hem\_cc\_splenium
</td>
<td style="text-align:right;">
73
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:right;">
87
</td>
<td style="text-align:right;">
16
</td>
</tr>
<tr>
<td style="text-align:left;">
hem\_subcort\_bilateral
</td>
<td style="text-align:right;">
91
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
104
</td>
<td style="text-align:right;">
12
</td>
</tr>
<tr>
<td style="text-align:left;">
flair\_bg\_unilat
</td>
<td style="text-align:right;">
38
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
47
</td>
<td style="text-align:right;">
19
</td>
</tr>
<tr>
<td style="text-align:left;">
hem\_cc\_corpus
</td>
<td style="text-align:right;">
57
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
66
</td>
<td style="text-align:right;">
14
</td>
</tr>
<tr>
<td style="text-align:left;">
flair\_cc\_corpus
</td>
<td style="text-align:right;">
49
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
57
</td>
<td style="text-align:right;">
14
</td>
</tr>
<tr>
<td style="text-align:left;">
flair\_subcort\_bilateral
</td>
<td style="text-align:right;">
49
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
57
</td>
<td style="text-align:right;">
14
</td>
</tr>
<tr>
<td style="text-align:left;">
hem\_tha\_unilat
</td>
<td style="text-align:right;">
30
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
38
</td>
<td style="text-align:right;">
21
</td>
</tr>
<tr>
<td style="text-align:left;">
diff\_cc\_corpus
</td>
<td style="text-align:right;">
16
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
23
</td>
<td style="text-align:right;">
30
</td>
</tr>
<tr>
<td style="text-align:left;">
flair\_tha\_bilateral
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
16
</td>
<td style="text-align:right;">
44
</td>
</tr>
<tr>
<td style="text-align:left;">
hem\_bg\_unilat
</td>
<td style="text-align:right;">
51
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
58
</td>
<td style="text-align:right;">
12
</td>
</tr>
<tr>
<td style="text-align:left;">
hem\_subcort\_unilat
</td>
<td style="text-align:right;">
60
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
67
</td>
<td style="text-align:right;">
10
</td>
</tr>
<tr>
<td style="text-align:left;">
diff\_cc\_genu
</td>
<td style="text-align:right;">
17
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
23
</td>
<td style="text-align:right;">
26
</td>
</tr>
<tr>
<td style="text-align:left;">
flair\_tha\_unilat
</td>
<td style="text-align:right;">
33
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
39
</td>
<td style="text-align:right;">
15
</td>
</tr>
<tr>
<td style="text-align:left;">
hem\_cc\_genu
</td>
<td style="text-align:right;">
36
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
42
</td>
<td style="text-align:right;">
14
</td>
</tr>
<tr>
<td style="text-align:left;">
hem\_tha\_bilateral
</td>
<td style="text-align:right;">
22
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
28
</td>
<td style="text-align:right;">
21
</td>
</tr>
<tr>
<td style="text-align:left;">
diff\_bg\_unilat
</td>
<td style="text-align:right;">
19
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
24
</td>
<td style="text-align:right;">
21
</td>
</tr>
<tr>
<td style="text-align:left;">
flair\_cc\_genu
</td>
<td style="text-align:right;">
22
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
27
</td>
<td style="text-align:right;">
19
</td>
</tr>
<tr>
<td style="text-align:left;">
diff\_tha\_unilat
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
18
</td>
<td style="text-align:right;">
22
</td>
</tr>
<tr>
<td style="text-align:left;">
flair\_bg\_bilateral
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
31
</td>
</tr>
<tr>
<td style="text-align:left;">
diff\_tha\_bilateral
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
75
</td>
</tr>
<tr>
<td style="text-align:left;">
hem\_bg\_bilateral
</td>
<td style="text-align:right;">
16
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
19
</td>
<td style="text-align:right;">
16
</td>
</tr>
<tr>
<td style="text-align:left;">
flair\_subcort\_unilat
</td>
<td style="text-align:right;">
49
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
51
</td>
<td style="text-align:right;">
4
</td>
</tr>
<tr>
<td style="text-align:left;">
diff\_bg\_bilateral
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
33
</td>
</tr>
<tr>
<td style="text-align:left;">
diff\_subcort\_bilateral
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
12
</td>
</tr>
<tr>
<td style="text-align:left;">
diff\_subcort\_unilat
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
16
</td>
<td style="text-align:right;">
6
</td>
</tr>
</tbody>
</table>

``` r
Data %>% 
   select(mortality,hem_subcort_unilat:infarct_cerebellum) %>% 
   select(contains("pons"), mortality) %>% 
   group_by(mortality) %>% 
   summarize_if(is.numeric, funs(sum(., na.rm = TRUE))) %>% 
   gather(key = "Variable", value = "Frequency", -mortality) %>% 
   na.omit() %>%
   ungroup() %>% 
   group_by(Variable) %>% 
   spread(mortality, Frequency) %>% 
   mutate(Total = Dead+Alive) %>% 
   mutate("Percent (Dead)" = round(Dead/Total
          *100, digits = 0)) %>% 
   arrange(desc(Dead)) %>% 
   kable(caption = "Percent mortalities following lesions in the pons")
```

<table>
<caption>
Percent mortalities following lesions in the pons
</caption>
<thead>
<tr>
<th style="text-align:left;">
Variable
</th>
<th style="text-align:right;">
Alive
</th>
<th style="text-align:right;">
Dead
</th>
<th style="text-align:right;">
Total
</th>
<th style="text-align:right;">
Percent (Dead)
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
hem\_pons\_ventral
</td>
<td style="text-align:right;">
17
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
26
</td>
<td style="text-align:right;">
35
</td>
</tr>
<tr>
<td style="text-align:left;">
flair\_pons\_bilateral
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:right;">
53
</td>
</tr>
<tr>
<td style="text-align:left;">
flair\_pons\_dorsal
</td>
<td style="text-align:right;">
25
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
33
</td>
<td style="text-align:right;">
24
</td>
</tr>
<tr>
<td style="text-align:left;">
flair\_pons\_ventral
</td>
<td style="text-align:right;">
18
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
25
</td>
<td style="text-align:right;">
28
</td>
</tr>
<tr>
<td style="text-align:left;">
hem\_pons\_dorsal
</td>
<td style="text-align:right;">
28
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
35
</td>
<td style="text-align:right;">
20
</td>
</tr>
<tr>
<td style="text-align:left;">
hem\_pons\_unilat
</td>
<td style="text-align:right;">
19
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
25
</td>
<td style="text-align:right;">
24
</td>
</tr>
<tr>
<td style="text-align:left;">
hem\_pons\_bilateral
</td>
<td style="text-align:right;">
17
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
22
</td>
<td style="text-align:right;">
23
</td>
</tr>
<tr>
<td style="text-align:left;">
diff\_pons\_unilat
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
17
</td>
<td style="text-align:right;">
24
</td>
</tr>
<tr>
<td style="text-align:left;">
diff\_pons\_ventral
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
50
</td>
</tr>
<tr>
<td style="text-align:left;">
diff\_pons\_bilateral
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
50
</td>
</tr>
<tr>
<td style="text-align:left;">
diff\_pons\_dorsal
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
15
</td>
</tr>
<tr>
<td style="text-align:left;">
flair\_pons\_unilat
</td>
<td style="text-align:right;">
31
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
33
</td>
<td style="text-align:right;">
6
</td>
</tr>
</tbody>
</table>

``` r
Data %>% 
   select(mortality,hem_subcort_unilat:infarct_cerebellum) %>% 
   select(contains("mes"), mortality) %>% 
   group_by(mortality) %>% 
   summarize_if(is.numeric, funs(sum(., na.rm = TRUE))) %>% 
   gather(key = "Variable", value = "Frequency", -mortality) %>% 
   na.omit() %>%
   ungroup() %>% 
   group_by(Variable) %>% 
   spread(mortality, Frequency) %>% 
   mutate(Total = Dead+Alive) %>% 
   mutate("Percent (Dead)" = round(Dead/Total
          *100, digits = 0)) %>%
   arrange(desc(Dead)) %>% 
   kable(caption = "Percent of unfavourable outcomes following lesions in the midbrain")
```

<table>
<caption>
Percent of unfavourable outcomes following lesions in the midbrain
</caption>
<thead>
<tr>
<th style="text-align:left;">
Variable
</th>
<th style="text-align:right;">
Alive
</th>
<th style="text-align:right;">
Dead
</th>
<th style="text-align:right;">
Total
</th>
<th style="text-align:right;">
Percent (Dead)
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
hem\_mes\_tegmentum
</td>
<td style="text-align:right;">
46
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
57
</td>
<td style="text-align:right;">
19
</td>
</tr>
<tr>
<td style="text-align:left;">
flair\_mes\_tegmentum
</td>
<td style="text-align:right;">
47
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
57
</td>
<td style="text-align:right;">
18
</td>
</tr>
<tr>
<td style="text-align:left;">
hem\_mes\_bilateral
</td>
<td style="text-align:right;">
30
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
38
</td>
<td style="text-align:right;">
21
</td>
</tr>
<tr>
<td style="text-align:left;">
hem\_mes\_unilat
</td>
<td style="text-align:right;">
34
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
42
</td>
<td style="text-align:right;">
19
</td>
</tr>
<tr>
<td style="text-align:left;">
flair\_mes\_bilateral
</td>
<td style="text-align:right;">
17
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
24
</td>
<td style="text-align:right;">
29
</td>
</tr>
<tr>
<td style="text-align:left;">
flair\_mes\_crus
</td>
<td style="text-align:right;">
25
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
32
</td>
<td style="text-align:right;">
22
</td>
</tr>
<tr>
<td style="text-align:left;">
diff\_mes\_unilat
</td>
<td style="text-align:right;">
28
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
34
</td>
<td style="text-align:right;">
18
</td>
</tr>
<tr>
<td style="text-align:left;">
flair\_mes\_tectum
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
21
</td>
<td style="text-align:right;">
29
</td>
</tr>
<tr>
<td style="text-align:left;">
diff\_mes\_tegmentum
</td>
<td style="text-align:right;">
18
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
23
</td>
<td style="text-align:right;">
22
</td>
</tr>
<tr>
<td style="text-align:left;">
flair\_mes\_unilat
</td>
<td style="text-align:right;">
50
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
55
</td>
<td style="text-align:right;">
9
</td>
</tr>
<tr>
<td style="text-align:left;">
hem\_mes\_tectum
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:right;">
33
</td>
</tr>
<tr>
<td style="text-align:left;">
diff\_mes\_crus
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
18
</td>
<td style="text-align:right;">
22
</td>
</tr>
<tr>
<td style="text-align:left;">
diff\_mes\_bilateral
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
27
</td>
</tr>
<tr>
<td style="text-align:left;">
diff\_mes\_tectum
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
33
</td>
</tr>
<tr>
<td style="text-align:left;">
hem\_mes\_crus
</td>
<td style="text-align:right;">
27
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
29
</td>
<td style="text-align:right;">
7
</td>
</tr>
</tbody>
</table>

``` r
# Data %>% 
#    select(dich_gos,hem_subcort_unilat:infarct_cerebellum) %>% 
#    select(contains("infarct"), dich_gos) %>% 
#    group_by(dich_gos) %>% 
#    summarize_if(is.numeric, funs(sum(., na.rm = TRUE))) %>% 
#    gather(key = "Variable", value = "Frequency", -dich_gos) %>% 
#    na.omit() %>%
#    ungroup() %>% 
#    group_by(Variable) %>% 
#    spread(dich_gos, Frequency) %>% 
#    mutate(Total = unfavourable+favourable) %>% 
#    mutate("Percent (unfavourable)" = round(unfavourable/Total
#           *100, digits = 0)) %>% 
#    arrange(desc(unfavourable)) %>% 
#    kable(caption = "Percent of unfavourable outcomes following infarcts")
```

continue pre-processing the data
--------------------------------

Identify all the MRI variables

``` r
imaging_vars <- Data %>% 
   select(hem_subcort_unilat:infarct_cerebellum) %>% 
   names() %>% 
   as_vector()
```

remove variables with fewer than 10 obvervations, as these have
essentially too few observations to be of any value in the analysis
(currently not run, but the code itself is kept for possible future
uses)

``` r
cols_to_remove <- imaging_vars[colSums(Data[imaging_vars], na.rm = TRUE) <= 12]

cols_to_keep <- imaging_vars[colSums(Data[imaging_vars], na.rm = TRUE) >= 12]
Data[!names(Data) %in% imaging_vars || names(Data) %in% cols_to_keep]
```

Coerce all binary imaging variables into factors instead of numeric
variables

``` r
Data <- Data %>% 
   mutate_at(imaging_vars, funs(factor(.)))
```

Drop variables with no identified lesions

``` r
Data <- Data[sapply(Data, function(x) !is.factor(x) | nlevels(x) > 1)]

Data <- Data %>% select_if(~sum(!is.na(.)) > 0)
```

David Nelson’s pupil fix
------------------------

``` r
# Fix pupils

# levels(Pupreak)
 

# Princip -> om en sida saknas imupteras det frön den andra
# om en sida saknas pga extern skada räknas som normal

attach(Data)

Pupiller<-rep(4,length(Data$Pupreak))

levels(Data$Pupreak)
```

    ##  [1] "NANorm"   "NAStel"   "NATrög"   "NK"       "Norm"     "NormStel"
    ##  [7] "NormTrög" "Stel"     "StelNorm" "StelTrög" "Trög"     "TrögNA"  
    ## [13] "TrögNorm" "TrögStel"

``` r
Pupreak<-Data$Pupreak

Pupiller[which(Pupreak=='Norm')]<-2
Pupiller[which(Pupreak=='TrögNorm')]<-2

Pupiller[which(Pupreak=='Trög')]<-2
Pupiller[which(Pupreak=='NormTrög')]<-2
Pupiller[which(Pupreak=='TrögNA')]<-2
Pupiller[which(Pupreak=='NormNorm')]<-2
Pupiller[which(Pupreak=='NATrög')]<-2

Pupiller[which(Pupreak=='TRÖG')]<-2
Pupiller[which(Pupreak=='TRÖ')]<-2
Pupiller[which(Pupreak=='TRÖ')]<-2
Pupiller[which(Pupreak=='TÖG')]<-2
Pupiller[which(Pupreak=='TRÖD')]<-2
Pupiller[which(Pupreak=='norm/stel (skada vä)')]<-2
Pupiller[which(Pupreak=='norm/stel (defekt vä)')]<-2


Pupiller[which(Pupreak=='TrögStel')]<-1
Pupiller[which(Pupreak=='NAStel')]<-1
Pupiller[which(Pupreak=='NormStel')]<-1
Pupiller[which(Pupreak=='StelTrög')]<-1
Pupiller[which(Pupreak=='StelNA')]<-1
Pupiller[which(Pupreak=='StelNorm')]<-1

Pupiller[which(Pupreak=='Stel')]<-0
Pupiller[which(Pupreak=='StelStel')]<-0
Pupiller[which(Pupreak=='Norm')]<-2

Pupiller[which(is.na(Pupreak))]<-NA
   
Pupiller[which(Pupreak=="NA")]<-NA
Pupiller[which(Pupreak=="u/norm(emaljöga hö)")]<-2
Pupiller[which(Pupreak== "" )]<-NA
Pupiller[which(Pupreak== "u" )]<-NA
Pupiller[which(Pupreak== "TöG"  )]<-2
Pupiller[which(Pupreak==  "Trög/Stel" )]<-1
Pupiller[which(Pupreak== "Trög/Norm"  )]<-2
Pupiller[which(Pupreak==  "Trög"  )]<-2
Pupiller[which(Pupreak=="TRöG")]<-2
Pupiller[which(Pupreak=="TRöD")]<-2
Pupiller[which(Pupreak== "StelNK" )]<-0
Pupiller[which(Pupreak=="Stel/Trög")]<-1
Pupiller[which(Pupreak=="stel/trög")]<-1
Pupiller[which(Pupreak=="stel/normal")]<-1
Pupiller[which(Pupreak=="stel/norm.\n" )]<-1
Pupiller[which(Pupreak=="Stel/Norm (skada hö)" )]<-2
Pupiller[which(Pupreak=="Stel/Norm")]<-1
Pupiller[which(Pupreak=="Stel/norm")]<-1
Pupiller[which(Pupreak=="stel/ej bedömbar")]<-0
Pupiller[which(Pupreak== "Stel" )]<-0
Pupiller[which(Pupreak== "stel")]<-0
Pupiller[which(Pupreak== "stel")]<-0
Pupiller[which(Pupreak== "saknar pupiller" )]<-NA
Pupiller[which(Pupreak=="op. bil.")]<-NA
Pupiller[which(Pupreak=="ONOR")]<-NA
Pupiller[which(Pupreak=="NormTrög")]<-2
Pupiller[which(Pupreak=="NormTrög")]<-2

Pupiller[which(Pupreak=="NATr\xf6g")]<-2

Pupiller[which(Pupreak=="StelTr\xf6g")]<-1
Pupiller[which(Pupreak=="Tr\xf6g")]<-2
Pupiller[which(Pupreak=="Tr\xf6gNA")]<-2
Pupiller[which(Pupreak=="Tr\xf6gNorm")]<-2
Pupiller[which(Pupreak=="Tr\xf6gStel")]<-1

Pupiller[which(Pupreak=="NormNK")]<-2
Pupiller[which(Pupreak=="NormNA")]<-2
Pupiller[which(Pupreak=="Norm/Trög")]<-2
Pupiller[which(Pupreak== "norm/trög")]<-2
Pupiller[which(Pupreak=="norm/stel(skada)")]<-2
Pupiller[which(Pupreak== "Norm/Stel(starr)")]<-2
Pupiller[which(Pupreak== "norm/stel (skada vö)" )]<-2
Pupiller[which(Pupreak=="norm/stel (defekt vö)")]<-2
Pupiller[which(Pupreak=="Norm/Stel")]<-1
Pupiller[which(Pupreak=="norm/stel" )]<-1
Pupiller[which(Pupreak=="Norm/NA" )]<-2
Pupiller[which(Pupreak=="NANorm" )]<-2
Pupiller[which(Pupreak=="Norm" )]<-2
Pupiller[which(Pupreak=="NORM" )]<-2
Pupiller[which(Pupreak=="Nej" )]<-0
Pupiller[which(Pupreak=="NK")]<-NA
Pupiller[which(Pupreak=="NKNorm")]<-2
Pupiller[which(Pupreak=="NATrög")]<-2
Pupiller[which(Pupreak=="NAStel")]<-0
Pupiller[which(Pupreak=="NANorm")]<-2
Pupiller[which(Pupreak==  "linsop. bil.  stela"  )]<-NA
Pupiller[which(Pupreak== "NANorm")]<-2
Pupiller[which(Pupreak==  "glaucomop. bil." )]<-NA
Pupiller[which(Pupreak==  "" )]<-NA
Pupiller[which(Pupreak=="0" )]<-0
Pupiller[which(Pupreak=="0,TR" )]<-1
Pupiller[which(Pupreak=="-/NORM" )]<-1
Pupiller[which(Pupreak=="stel/norm" )]<-1
Pupiller[which(Pupreak=="TRö" )]<-2
Pupiller[which(Pupreak=="Norm/Stel (starr)" )]<-2
Pupiller[which(Pupreak=="Normtrög)" )]<-2
Pupiller[which(Pupreak=="Normtrög" )]<-2
Pupiller[which(Pupreak=="NormTr\xf6g" )]<-2
Pupiller[which(Pupreak=="Normtrög" )]<-2
Pupiller[which(Pupreak=="Normtrög" )]<-2
Pupiller[which(Pupreak=="Steltrög" )]<-1
Pupiller[which(Pupreak=="stel/norm.\r\n")]<-1



Pupiller<-as.factor(Pupiller)
Data$Pupiller_old<-Data$Pupiller
Data$Pupiller<-Pupiller


ifelse(any(Pupiller==4),print('Some Pupil names not caught in program'), return(Data$Pupiller<-Pupiller))
```

    ## [1] NA

``` r
remove (Pupiller)

# To check
# levels(Pupiller)
# which(Pupiller==4)
# Pupreak[which(Pupiller==4)]
# levels(Pupreak[which(Pupiller==4)])
# levels(Data$Pupiller[which(Pupiller==4)])
# 
detach(Data)

Data <- Data %>% 
   select(-Pupreak) %>% 
   mutate(Pupiller = factor(Pupiller, ordered = TRUE))
```

mutate MRI scores as ordered “factors”
--------------------------------------

``` r
Data <- Data %>% 
   mutate(adams = factor(adams, ordered = TRUE)) %>% 
   mutate(firsching = factor(firsching, ordered = TRUE)) %>% 
   mutate(hamdeh = factor(hamdeh, ordered = TRUE)) %>%  
   mutate(stockholm_mri = factor(stockholm_mri, ordered = TRUE)) %>% 
   select(-c(helsinki, marshall, mri_time, stockholm_no_dai, karolinska, contains("infarct")))
```

Data summary
------------

``` r
 Data %>% map_df(~(data.frame(class = class(.x),
                              unique_values = n_distinct(.x, na.rm = TRUE))),
                     .id = "variable") %>% 
   kable(caption = 
            "The variables, their class and the number of unique values",
         format = "latex", booktabs = T,linesep = "", 
         longtable = T) %>% 
   kable_styling(full_width = F, latex_options = "repeat_header")
```

Missing
=======

``` r
Data %>%
    gather(key = "key", value = "val") %>%
    mutate(is_missing = is.na(val)) %>%
    group_by(key, is_missing) %>%
    summarise(Missing = n()) %>%
    filter(is_missing==T) %>%
    select(-is_missing) %>%
    arrange(desc(Missing)) %>% 
   kable(caption = "Missing data",
         format = "latex", booktabs = T,linesep = "", 
         longtable = T) %>% 
   kable_styling(full_width = F, latex_options = "repeat_header")
```

``` r
Data %>%
    gather(key = "key", value = "val") %>% 
   summarize(average_missing_percent = round(mean(is.na(val))*100)) %>% kable()
```

<table>
<thead>
<tr>
<th style="text-align:right;">
average\_missing\_percent
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:right;">
5
</td>
</tr>
</tbody>
</table>

Row plot

``` r
missing_values <- Data %>%
  gather(key = "key", value = "val") %>%
  mutate(is_na = is.na(val)) %>%
  group_by(key) %>%
  mutate(total = n()) %>%
  group_by(key, total, is_na) %>%
  summarise(num.isna = n()) %>%
  mutate(pct = num.isna / total * 100)


levels <-
    (missing_values  %>% filter(is_na == T) %>% arrange(desc(pct)))$key

print(Data %>%
  mutate(id = row_number()) %>%
  gather(-id, key = "key", value = "val") %>%
  mutate(isna = is.na(val)) %>%
  ggplot(aes(key, id, fill = isna)) +
    geom_raster(alpha=0.8) +
    scale_fill_manual(name = "",
        values = c('steelblue', 'tomato3'),
        labels = c("Present", "Missing")) +
    scale_x_discrete(limits = levels) +
    labs(x = "Variable",
           y = "Row Number", title = "Missing values in rows") +
    coord_flip())
```

![](mri_analysis_files/figure-markdown_github/unnamed-chunk-18-1.png)

Univariate logistic regression
==============================

list of variable names

``` r
var_names <- Data %>% 
   select(-c(dich_gos, FinalGOS, id)) %>%
   names() 
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

Fitting univariate logistic regression
--------------------------------------

Running univariate logistic regression for all numeric variables
included in this study

``` r
formulas <- paste('dich_gos ~ ', var_names) #formulas to be used in the the glm function

model_null <- glm(dich_gos ~ 1, family = binomial(link="logit"), 
                     data = Data)

univariate_analysis <- var_names %>% 
   tibble() %>% 
   rename(Variables = ".") %>% 
   mutate(fitted_models = map(formulas, ~ glm(formula = .x, 
                                  family = binomial(link="logit"),
                                                    data= Data))) %>%
   group_by(Variables) %>%
   mutate(tid = map(fitted_models, ~ tidy(.))) %>% #extracting the coefficients
   mutate(tid = map(tid, ~select(.x, term, estimate, p.value))) %>% 
   unnest(tid, .drop = F) %>% 
   rename(p_value = "p.value", Coefficient = "estimate") %>% 
   filter(!str_detect(term, "Intercept")) %>% 
   select(-term) %>% 
   mutate(OR = exp(Coefficient)) %>% 
   arrange((p_value)) %>% 
   mutate_at(vars(p_value), funs(round(.,4))) %>% 
   mutate(glance = map(fitted_models, ~ glance(.))) %>% 
   mutate(aug = map(fitted_models, ~ augment(.))) %>% 
   mutate(r_sq = map_dbl(fitted_models, 
                         ~ r2_calc(data = Data, model = .,
                                   nullmodel = model_null))) %>%
   select(-c(fitted_models,glance,aug, Coefficient)) 

univariate_analysis %>% 
   ungroup() %>% 
   mutate(Variables = str_replace_all(string = Variables, pattern = "_",
                                      replacement = " ")) %>% 
  kable(format = "latex", escape = FALSE, booktabs = TRUE, linesep = "",
        caption = "Results form the univariate logistic regression",
        digits = 2, align = 'c',longtable = TRUE, col.names = c("Variables", "p-value","OR","$R^2$")) %>%
   kable_styling(full_width = F, latex_options = "repeat_header")
```

Select variables with a p value that is below 0.05 in the univariate
logistic regression

``` r
best_vars <- univariate_analysis %>% 
   filter(p_value < 0.05) %>% select(Variables) %>% .[["Variables"]]

core_model <- c("Alder", "GCSKS", "hypoxia", "hypotension", "Pupiller")

# Remove variables that are part of the core model from the "best_vars" vector
# As these will be added anyway and in case e.g. GCS was significant
# in the univariate analysis, I would like for there not to be two such
# variables in the final model.

best_vars <- setdiff(best_vars, core_model)

best_vars %>% 
   tidy() %>% 
   rename("Variables" = x) %>% 
   kable(format = "latex", booktabs = T, linesep = "",
        caption = "MRI variables that were related to outcome in the univariate logistic regression analysis",
        digits = 2, align = 'c',longtable = TRUE) %>%
   kable_styling(full_width = F, latex_options = "repeat_header")
```

check how many variables have no lesions at all

``` r
setdiff(var_names, univariate_analysis$Variables) %>% 
   tibble() %>% 
   rename("Variables with no lesions identified" = "." ) %>% 
   kable()
```

Multiple imputation
===================

``` r
imputed <- mice(data = Data, seed = 123, m = 10, maxit = 30) 

# in order to save the data so that it may be used for e.g. SVM.
saveRDS(imputed,'imputed_data.rds')
```

The code chunk above is only run the first time. For all subsequent
coding sessions, the code below will be run.

``` r
#this is just so that I don't have to run mice() all the time

imputed <- readRDS('imputed_data.rds')
```

We now have 7 imputed datasets. If we wish to fit a logistic regression
model using this data, we would need to do so 7 times.

However, by nesting all 7 datasets and fitting the models in an
iterative manner, this can be achieved rather simply.

Fitting glm models to imputed datasets
======================================

Fitting the GLM models and extracting their coefficients
--------------------------------------------------------

Iteratively fitting the logistic regression model for all seven imputed
datasets, simultaneously.

``` r
imputed_log <- imputed %>%
  mice::complete("long") %>%
  group_by(.imp) %>%
  select("Imputation" =.imp, dich_gos,
         all_of(best_vars), all_of(core_model)) %>%
   group_by(Imputation) %>%  nest() 
 
best_vars_formula <- paste("~", paste(best_vars, collapse = " + "),
                            sep = " " ) %>%
   paste("+", paste(core_model, collapse = " + "), sep = " " ) %>%
   as.formula()

core_vars_formula <- paste("dich_gos ~", paste(core_model, collapse = " + "),
                            sep = " " ) %>% 
   as.formula()

# scope <- list(upper = best_vars_formula,
#              lower = ~1)

# expression <- expression(f1 <- glm(dich_gos ~ 1, family =  binomial(link="logit")),
#                       f2 <- stats::step(f1, scope = scope))
# 
# 
# fit <- with(imputed, expression)
# 
#  formulas <- lapply(fit$analyses, formula)
#  terms <- lapply(formulas, terms)
#  votes <- unlist(lapply(terms, labels))
#  table(votes)
# 
# fit.without <- with(imputed, glm(dich_gos ~hem_mes_tegmentum + flair_bg_unilat +
#      infarct_media +
#     Alder +
#     GCSKS + hypoxia + hypotension + Pupiller, family =  binomial(link="logit")))
# 
# fit.with <- with(imputed, exp = expression(glm(dich_gos ~hem_mes_tegmentum +  flair_bg_unilat +
#    infarct_media +
#      Alder +
#     GCSKS + hypoxia + hypotension + Pupiller, family =  binomial(link="logit"))))
# D1(fit.with, fit.without)
# 
#  scope <- list(upper = best_vars_formula_two,
#               lower = ~1)

best_vars_formula_two <- as.formula(dich_gos ~ hem_mes_tegmentum +  flair_bg_unilat + 
    Alder + 
    GCSKS + hypoxia + hypotension + Pupiller)  


estimates <- imputed_log %>%
mutate(glm = map(data, ~ glm(best_vars_formula_two, family = binomial(link="logit"), data = .))) %>%
mutate(null = map(data, ~ glm(dich_gos ~ 1, family = binomial(link="logit"), data = .))) %>%
   mutate(core = map(data, ~ glm(core_vars_formula, family = binomial(link="logit"), data = .))) %>%
   mutate(r_sq = pmap_dbl(list(data, glm, null), r2_calc)) %>%
   mutate(r_sq_core = pmap_dbl(list(data, core, null), r2_calc)) %>%
mutate(glance = map(glm, ~ glance(.))) %>%
mutate(aug = map(glm, ~ augment(.))) %>%
mutate(tid = map(glm, ~ tidy(.))) %>%
mutate(glance_core = map(core, ~ glance(.))) %>% 
mutate(aug_core = map(core, ~ augment(.))) %>%
select(-null)

best_vars_two <- estimates %>%  ungroup() %>% 
   sample_n(size =  1, replace = F) %>% 
   unnest(tid, .drop = T) %>% 
   ungroup() %>% 
      pull(term) %>% 
     as_vector()
```

results
-------

Results for all logistic regression models that were fitted to all
imputed datasets, respectively.

``` r
estimates %>%
  unnest(glance) %>%
  dplyr::select(Imputation, r_sq, null.deviance:df.residual) %>% 
  kable(format = "latex", escape = FALSE, booktabs = T, linesep = "",
        caption = "Comparing the results for all imputed datasets",
        digits = 2, align = 'c', col.names = c("Imputation", "$R^2$",
                                               "Null deviance",
                                                   "df.null", "logLik",
                                                   "AIC", "BIC",
                    "deviance", "df.residuals")) %>%
   kable_styling(full_width = F)
```

average

``` r
core_comparison <- estimates %>%
  unnest(glance) %>%  
  dplyr::select(r_sq,logLik:deviance, null.deviance) %>% 
   ungroup() %>% 
summarize_all(funs(mean)) %>% 
   mutate(Model = "Core + MRI")

core_comparison <- estimates %>%
  unnest(glance_core) %>%  
   ungroup() %>% 
  dplyr::select("r_sq" = r_sq_core,logLik:deviance, null.deviance) %>% 
summarize_all(funs(mean)) %>% 
   mutate(Model = "Core") %>% 
   bind_rows(core_comparison) %>% 
   select(Model, everything()) %>% 
   select(-Imputation)

core_comparison %>% 
  kable(format = "latex", escape = FALSE, booktabs = T, linesep = "",
        caption = "Average result for all imputed datasets",
        digits = 2, align = 'c', col.names = c("Model","$R^2$", "logLik",
                                                   "AIC", "BIC",
                    "deviance", "Null deviance")) %>%
   kable_styling(full_width = F)
```

Results for each individual variable..
--------------------------------------

Results for all the imputed datasets combined

``` r
vars_table <- summary(mice::pool(estimates$glm)) %>%
as_tibble() %>%
mutate(variable = best_vars_two) %>%
mutate(OR = exp(estimate)) %>%
filter(variable != "(Intercept)") %>%
dplyr::select(variable, everything(), -c(statistic, estimate, df)) 

# Two different coefficients are presented for the Pupils variable,
# as this variable has three levels. 
#Thus, a new vector containing two different 
# names for the pupils is required

pupils <- vars_table %>% 
   filter(variable == "Pupiller1" | variable == "Pupiller2") 

vars_table %>% mutate(variable = str_replace_all(string = variable,
                                                 pattern = "\\d",
                                  replacement = "")) %>%
   filter(variable != "Pupiller") %>%
   bind_rows(pupils) %>%
   arrange(desc(OR)) %>% 
  kable(caption = "Combined results of all imputated datasets", align = 'c',
        format = "latex", booktabs = T, linesep = "", longtable = TRUE,  digits = 2) %>%
   kable_styling(full_width = F, latex_options = "repeat_header")
```

ROC
---

``` r
auroc <- estimates %>% 
   select(-c(glm,glance, tid, r_sq)) %>% 
    ungroup()

auroc_1 <- auroc %>% 
   filter(Imputation == "1") %>% 
   unnest(aug) %>% 
   roc(dich_gos,.fitted)
   
 
p1 <- ggroc(auroc_1, alpha = 0.7, colour = "black", 
                                       linetype = 1, size = 2) + 
   theme_minimal() +  
    geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="grey", linetype="dashed") +  annotate("text", x = .9, y = .2, 
       label = paste("AUC =", round(auc(auroc_1), digits = 3))) + xlab("False Positive Rate (1-Specificity)") + ylab("True Positive Rate (Sensitivity)")

auroc_2 <- auroc %>% 
   filter(Imputation == "2") %>% 
   unnest(aug) %>% 
   roc(dich_gos,.fitted)
   
 
p2 <- ggroc(auroc_2, alpha = 0.7, colour = "black", 
                                       linetype = 1, size = 2) + 
   theme_minimal() + 
    geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="grey", linetype="dashed") +  annotate("text", x = .9, y = .2, 
       label = paste("AUC =", round(auc(auroc_2), digits = 3))) + xlab("False Positive Rate (1-Specificity)") + ylab("True Positive Rate (Sensitivity)")

auroc_3 <- auroc %>% 
   filter(Imputation == "3") %>% 
   unnest(aug) %>% 
   roc(dich_gos,.fitted)
   
 
p3 <- ggroc(auroc_3, alpha = 0.7, colour = "black", 
                                       linetype = 1, size = 2) + 
   theme_minimal() + 
    geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="grey", linetype="dashed") +  annotate("text", x = .9, y = .2, 
       label = paste("AUC =", round(auc(auroc_3), digits = 3))) + xlab("False Positive Rate (1-Specificity)") + ylab("True Positive Rate (Sensitivity)")

auroc_4 <- auroc %>% 
   filter(Imputation == "4") %>% 
   unnest(aug) %>% 
   roc(dich_gos,.fitted)
   
 
p4 <- ggroc(auroc_4, alpha = 0.7, colour = "black", 
                                       linetype = 1, size = 2) + 
   theme_minimal() + 
    geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="grey", linetype="dashed") +  annotate("text", x = .9, y = .2, 
       label = paste("AUC =", round(auc(auroc_4), digits = 3))) + xlab("False Positive Rate (1-Specificity)") + ylab("True Positive Rate (Sensitivity)")

auroc_5 <- auroc %>% 
   filter(Imputation == "5") %>% 
   unnest(aug) %>% 
   roc(dich_gos,.fitted)
   
 
p5 <- ggroc(auroc_5, alpha = 0.7, colour = "black", 
                                       linetype = 1, size = 2) + 
   theme_minimal() + 
    geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="grey", linetype="dashed") +  annotate("text", x = .9, y = .2, 
       label = paste("AUC =", round(auc(auroc_5), digits = 3))) + xlab("False Positive Rate (1-Specificity)") + ylab("True Positive Rate (Sensitivity)")

auroc_6 <- auroc %>% 
   filter(Imputation == "6") %>% 
   unnest(aug) %>% 
   roc(dich_gos,.fitted)
   
 
p6 <- ggroc(auroc_6, alpha = 0.7, colour = "black", 
                                       linetype = 1, size = 2) + 
   theme_minimal() +  
    geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="grey", linetype="dashed") +  annotate("text", x = .9, y = .2, 
       label = paste("AUC =", round(auc(auroc_6), digits = 3))) + xlab("False Positive Rate (1-Specificity)") + ylab("True Positive Rate (Sensitivity)")

auroc_7 <- auroc %>% 
   filter(Imputation == "7") %>% 
   unnest(aug) %>% 
   roc(dich_gos,.fitted)
   
 
p7 <- ggroc(auroc_7, alpha = 0.7, colour = "black", 
                                       linetype = 1, size = 2) + 
   theme_minimal() + 
    geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="grey", linetype="dashed") +  annotate("text", x = .9, y = .2, 
       label = paste("AUC =", round(auc(auroc_7), digits = 3))) + xlab("False Positive Rate (1-Specificity)") + ylab("True Positive Rate (Sensitivity)")

auroc_core <- auroc %>% 
   filter(Imputation == "7") %>% 
   unnest(aug_core) %>% 
   roc(dich_gos,.fitted)
   
 
p8 <- ggroc(auroc_core, alpha = 0.7, colour = "black", 
                                       linetype = 1, size = 2) + 
   theme_minimal() + 
    geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="grey", linetype="dashed") +  annotate("text", x = .9, y = .2, 
       label = paste("AUC =", round(auc(auroc_core), digits = 3))) + xlab("False Positive Rate (1-Specificity)") + ylab("True Positive Rate (Sensitivity)")

title <- ggdraw() + 
  draw_label(
    "ROC curves for each imputed dataset",
    fontface = 'bold',
    x = 0.5, size = 24
  )

labs <- c(auroc$Imputation, "Core")

roc_plots <- plot_grid(p1,p2,p3,p4,p5,p6,p7,p8, rows = 4, cols = 2, labels = labs,
  label_size = 18,
  align = "v",label_x = 0.15, label_y = 1.05
)
```

plot

``` r
plot_grid(title, roc_plots, ncol = 1,
  rel_heights = c(0.1, 1.2))
```

![](mri_analysis_files/figure-markdown_github/unnamed-chunk-31-1.png)

Session information
===================

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
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] rcompanion_2.3.25               DescTools_0.99.32              
    ##  [3] pROC_1.16.1                     glmnet_3.0-2                   
    ##  [5] Matrix_1.2-18                   cowplot_1.0.0                  
    ##  [7] naniar_0.5.0                    pscl_1.5.2                     
    ##  [9] mice_3.8.0                      lubridate_1.7.8                
    ## [11] corrr_0.4.1                     table1_1.1                     
    ## [13] caTools_1.18.0                  descr_1.1.4                    
    ## [15] e1071_1.7-3                     modelgrid_1.1.1.0              
    ## [17] readxl_1.3.1                    corrplot_0.84                  
    ## [19] summarytools_0.9.5              kableExtra_1.1.0               
    ## [21] ggfortify_0.4.8                 MLmetrics_1.1.1                
    ## [23] caret_6.0-85                    lattice_0.20-38                
    ## [25] AppliedPredictiveModeling_1.1-7 knitr_1.28                     
    ## [27] forcats_0.5.0                   stringr_1.4.0                  
    ## [29] readr_1.3.1                     tidyverse_1.3.0                
    ## [31] yardstick_0.0.5                 workflows_0.1.0                
    ## [33] tune_0.0.1                      tibble_3.0.0                   
    ## [35] rsample_0.0.5                   tidyr_1.0.2                    
    ## [37] recipes_0.1.9                   purrr_0.3.3                    
    ## [39] parsnip_0.0.5                   infer_0.5.1                    
    ## [41] ggplot2_3.3.0                   dplyr_0.8.5                    
    ## [43] dials_0.0.4                     scales_1.1.0                   
    ## [45] broom_0.5.5                     tidymodels_0.1.0               
    ## [47] MASS_7.3-51.5                  
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] tidyselect_1.0.0     lme4_1.1-23          htmlwidgets_1.5.1   
    ##   [4] grid_3.6.2           CORElearn_1.54.2     munsell_0.5.0       
    ##   [7] codetools_0.2-16     statmod_1.4.34       DT_0.13             
    ##  [10] future_1.16.0        miniUI_0.1.1.1       withr_2.1.2         
    ##  [13] colorspace_1.4-1     highr_0.8            rstudioapi_0.11     
    ##  [16] stats4_3.6.2         bayesplot_1.7.1      listenv_0.8.0       
    ##  [19] labeling_0.3         rstan_2.19.3         farver_2.0.3        
    ##  [22] DiceDesign_1.8-1     TH.data_1.0-10       vctrs_0.2.4         
    ##  [25] generics_0.0.2       ipred_0.9-9          xfun_0.12           
    ##  [28] R6_2.4.1             markdown_1.1         rstanarm_2.19.3     
    ##  [31] bitops_1.0-6         lhs_1.0.1            assertthat_0.2.1    
    ##  [34] promises_1.1.0       multcomp_1.4-12      nnet_7.3-12         
    ##  [37] gtable_0.3.0         multcompView_0.1-8   globals_0.12.5      
    ##  [40] processx_3.4.2       sandwich_2.5-1       timeDate_3043.102   
    ##  [43] rlang_0.4.5          EMT_1.1              splines_3.6.2       
    ##  [46] ModelMetrics_1.2.2.1 rapportools_1.0      checkmate_2.0.0     
    ##  [49] inline_0.3.15        yaml_2.2.1           reshape2_1.4.3      
    ##  [52] modelr_0.1.5         tidytext_0.2.2       threejs_0.3.3       
    ##  [55] crosstalk_1.0.0      backports_1.1.5      httpuv_1.5.2        
    ##  [58] rsconnect_0.8.16     tokenizers_0.2.1     tools_3.6.2         
    ##  [61] lava_1.6.6           tcltk_3.6.2          ellipsis_0.3.0      
    ##  [64] ggridges_0.5.2       Rcpp_1.0.3           plyr_1.8.6          
    ##  [67] base64enc_0.1-3      ps_1.3.2             prettyunits_1.1.1   
    ##  [70] rpart_4.1-15         zoo_1.8-7            haven_2.2.0         
    ##  [73] cluster_2.1.0        fs_1.3.1             furrr_0.1.0         
    ##  [76] magrittr_1.5         data.table_1.12.8    magick_2.3          
    ##  [79] lmtest_0.9-37        colourpicker_1.0     reprex_0.3.0        
    ##  [82] GPfit_1.0-8          mvtnorm_1.1-0        SnowballC_0.6.0     
    ##  [85] matrixStats_0.55.0   tidyposterior_0.0.2  hms_0.5.3           
    ##  [88] shinyjs_1.1          mime_0.9             evaluate_0.14       
    ##  [91] xtable_1.8-4         tidypredict_0.4.5    shinystan_2.5.0     
    ##  [94] shape_1.4.4          gridExtra_2.3        rstantools_2.0.0    
    ##  [97] compiler_3.6.2       ellipse_0.4.1        rpart.plot_3.0.8    
    ## [100] crayon_1.3.4         minqa_1.2.4          StanHeaders_2.21.0-1
    ## [103] htmltools_0.4.0      later_1.0.0          Formula_1.2-3       
    ## [106] libcoin_1.0-5        visdat_0.5.3         expm_0.999-4        
    ## [109] DBI_1.1.0            dbplyr_1.4.2         boot_1.3-23         
    ## [112] cli_2.0.2            pryr_0.1.4           parallel_3.6.2      
    ## [115] gower_0.2.1          igraph_1.2.4.2       pkgconfig_2.0.3     
    ## [118] coin_1.3-1           xml2_1.2.2           foreach_1.4.8       
    ## [121] dygraphs_1.1.1.6     webshot_0.5.2        prodlim_2019.11.13  
    ## [124] rvest_0.3.5          janeaustenr_0.1.5    callr_3.4.3         
    ## [127] digest_0.6.25        rmarkdown_2.1        cellranger_1.1.0    
    ## [130] nortest_1.0-4        modeltools_0.2-22    shiny_1.4.0         
    ## [133] gtools_3.8.1         nloptr_1.2.2.1       lifecycle_0.2.0     
    ## [136] nlme_3.1-142         jsonlite_1.6.1       viridisLite_0.3.0   
    ## [139] fansi_0.4.1          pillar_1.4.3         loo_2.2.0           
    ## [142] fastmap_1.0.1        httr_1.4.1           plotrix_3.7-7       
    ## [145] pkgbuild_1.0.6       survival_3.1-8       glue_1.3.1          
    ## [148] xts_0.12-0           shinythemes_1.1.2    iterators_1.0.12    
    ## [151] pander_0.6.3         class_7.3-15         stringi_1.4.6
