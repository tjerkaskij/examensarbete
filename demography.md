-   [Aim](#aim)
-   [Packages](#packages)
    -   [Install packages](#install-packages)
    -   [Load Packages](#load-packages)
-   [The data](#the-data)
-   [Pre-process the data table](#pre-process-the-data-table)
    -   [select variables](#select-variables)
    -   [structure of the data](#structure-of-the-data)
        -   [Changing characters to
            factors](#changing-characters-to-factors)
-   [Shapiro-wilk test](#shapiro-wilk-test)
-   [Mann–Whitney U-test](#mannwhitney-u-test)
-   [*χ*<sup>2</sup> - test](#chi2---test)

Aim
===

The aim of this assignment is to investigate of cortical thickness can
be used  
to succesfully predict the correct diagnosis of patients with movement
disorders.

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
```

Load Packages
-------------

``` r
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

options(knitr.kable.NA = '')
```

The data
========

``` r
data <- read_excel('AV-id dataset till Jonathan 191020.xlsx', na = "NA") 
```

Pre-process the data table
==========================

select variables
----------------

``` r
data <- data %>% 
  select(Patientnummer, MRT, Alder, GCSKS, ISS, Pupreak, Pupstrl, GCSOpl, Sed, Ethyl, BltrOpl, ICP, Kon) 
```

structure of the data
---------------------

``` r
str(data)
```

    ## tibble [1,652 x 13] (S3: tbl_df/tbl/data.frame)
    ##  $ Patientnummer: num [1:1652] 1 2 3 4 5 6 7 8 9 10 ...
    ##  $ MRT          : chr [1:1652] "Nej" "Nej" "Ja" "Nej" ...
    ##  $ Alder        : num [1:1652] 66 34 25 72 16 62 57 63 32 80 ...
    ##  $ GCSKS        : num [1:1652] 3 15 3 12 3 7 14 4 3 5 ...
    ##  $ ISS          : num [1:1652] 0 0 0 0 0 0 0 0 0 0 ...
    ##  $ Pupreak      : chr [1:1652] "TrögStel" "Norm" "Stel" "Norm" ...
    ##  $ Pupstrl      : num [1:1652] 12 22 33 22 11 11 22 13 11 22 ...
    ##  $ GCSOpl       : num [1:1652] 10 15 3 12 6 3 14 5 3 14 ...
    ##  $ Sed          : chr [1:1652] "Ja" "Nej" "Ja" "Nej" ...
    ##  $ Ethyl        : chr [1:1652] "Ja" "Ja" "Ja" "Nej" ...
    ##  $ BltrOpl      : num [1:1652] 200 120 0 130 80 80 0 100 110 130 ...
    ##  $ ICP          : chr [1:1652] "Ja" "Nej" "Ja" "Nej" ...
    ##  $ Kon          : chr [1:1652] "Man" "Man" "Man" "Man" ...

### Changing characters to factors

``` r
data <- data %>% 
  mutate(Pupreak = str_replace(Pupreak, "NA", "")) %>% 
  mutate_if(is.character, as.factor) 
```

Shapiro-wilk test
=================

A graphical analysis of the data is perhaps a bette roption than this.
However, I currently do not have the time for that.

``` r
numerics <- names(data)[sapply(data, is.numeric)] 

data %>% 
  select(MRT, numerics) %>% 
  select(-Patientnummer) %>% 
    gather(key = "variable_name", value = "value", Alder:BltrOpl) %>% 
    group_by(variable_name, MRT)  %>% 
    do(tidy(shapiro.test(.$value))) %>% 
    ungroup() %>% 
    select(-method) %>% 
    filter(!is.na(MRT)) %>% 
    kable(digits = 3, caption = "Shapiro-wilt test for normality", booktabs = TRUE)
```

<table>
<caption>
Shapiro-wilt test for normality
</caption>
<thead>
<tr>
<th style="text-align:left;">
variable\_name
</th>
<th style="text-align:left;">
MRT
</th>
<th style="text-align:right;">
statistic
</th>
<th style="text-align:right;">
p.value
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Alder
</td>
<td style="text-align:left;">
Ja
</td>
<td style="text-align:right;">
0.934
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Alder
</td>
<td style="text-align:left;">
Nej
</td>
<td style="text-align:right;">
0.972
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
BltrOpl
</td>
<td style="text-align:left;">
Ja
</td>
<td style="text-align:right;">
0.852
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
BltrOpl
</td>
<td style="text-align:left;">
Nej
</td>
<td style="text-align:right;">
0.830
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
GCSKS
</td>
<td style="text-align:left;">
Ja
</td>
<td style="text-align:right;">
0.786
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
GCSKS
</td>
<td style="text-align:left;">
Nej
</td>
<td style="text-align:right;">
0.829
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
GCSOpl
</td>
<td style="text-align:left;">
Ja
</td>
<td style="text-align:right;">
0.908
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
GCSOpl
</td>
<td style="text-align:left;">
Nej
</td>
<td style="text-align:right;">
0.849
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
ISS
</td>
<td style="text-align:left;">
Ja
</td>
<td style="text-align:right;">
0.874
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
ISS
</td>
<td style="text-align:left;">
Nej
</td>
<td style="text-align:right;">
0.804
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Pupstrl
</td>
<td style="text-align:left;">
Ja
</td>
<td style="text-align:right;">
0.866
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Pupstrl
</td>
<td style="text-align:left;">
Nej
</td>
<td style="text-align:right;">
0.787
</td>
<td style="text-align:right;">
0
</td>
</tr>
</tbody>
</table>

``` r
# Initial attempt, which did not look as nice as 
# the current version:
# other <- data %>% 
#   select(MRT, numerics) %>% 
#   select(-Patientnummer) %>% 
# group_by(MRT) %>%
#   summarise_all(.funs = funs(statistic = shapiro.test(.)$statistic, 
#                              p.value = shapiro.test(.)$p.value))
```

Mann–Whitney U-test
===================

For *numeric* data

``` r
data %>% 
  select_if(is.numeric) %>%
  select(- Patientnummer) %>% 
  map_df(~wilcox.test(. ~ data$MRT)$p.value) %>% 
  kable(digits = 3, caption = "Mann-Whitney U test", booktabs = TRUE)
```

<table>
<caption>
Mann-Whitney U test
</caption>
<thead>
<tr>
<th style="text-align:right;">
Alder
</th>
<th style="text-align:right;">
GCSKS
</th>
<th style="text-align:right;">
ISS
</th>
<th style="text-align:right;">
Pupstrl
</th>
<th style="text-align:right;">
GCSOpl
</th>
<th style="text-align:right;">
BltrOpl
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.981
</td>
</tr>
</tbody>
</table>

*χ*<sup>2</sup> - test
======================

$$\\chi^2 = \\sum \\frac {(O - E)^2}{E}$$

Used for *categorical* data

``` r
chi_sq_p_value <- data %>%
  select_if(is.factor) %>%
  drop_na() %>% 
    summarise_at(2:6,funs(chisq.test(., 
           MRT, simulate.p.value = TRUE)$p.value)) %>% 
  mutate(Estimate = "p_value")
 

chi_sq_statistic <- data %>%
  select_if(is.factor) %>% 
  drop_na() %>% 
    summarise_at(2:6,funs(chisq.test(., 
           MRT)$statistic)) %>% 
  mutate(Estimate = "Statistic")

chi_sq_resultss <- bind_rows(chi_sq_statistic, chi_sq_p_value) %>% 
  select(Estimate, everything())

chi_sq_resultss %>% 
  kable(digits = 3, booktabs = TRUE, caption = "$\\chi^2$ - test")
```

<table>
<caption>
*χ*<sup>2</sup> - test
</caption>
<thead>
<tr>
<th style="text-align:left;">
Estimate
</th>
<th style="text-align:right;">
Pupreak
</th>
<th style="text-align:right;">
Sed
</th>
<th style="text-align:right;">
Ethyl
</th>
<th style="text-align:right;">
ICP
</th>
<th style="text-align:right;">
Kon
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Statistic
</td>
<td style="text-align:right;">
175.535
</td>
<td style="text-align:right;">
109.987
</td>
<td style="text-align:right;">
0.625
</td>
<td style="text-align:right;">
300.926
</td>
<td style="text-align:right;">
0.449
</td>
</tr>
<tr>
<td style="text-align:left;">
p\_value
</td>
<td style="text-align:right;">
0.000
</td>
<td style="text-align:right;">
0.000
</td>
<td style="text-align:right;">
0.724
</td>
<td style="text-align:right;">
0.000
</td>
<td style="text-align:right;">
0.520
</td>
</tr>
</tbody>
</table>
