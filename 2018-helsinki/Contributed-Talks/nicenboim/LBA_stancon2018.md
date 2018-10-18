---
title: "The implementation of a model of choice: the (truncated) linear ballistic accumulator."
author: 
  - name: B. Nicenboim ^[Thanks to Henrik Singmann and Andrew Heathcote for their assistance in reverse engineering the R package rtdists and thanks to Gabriel Tillman for his helpful advice on the parameters of the LBA. Thanks also to Shravan Vasishth for his comments. ]
date: 
output:
  html_document:
    toc: true
    number_sections: true
    fig_caption: true
    css: styles.css
bibliography: bibliography.bib
---

<script type="text/x-mathjax-config">
MathJax.Hub.Config({
  TeX: { equationNumbers: { autoNumber: "AMS" } }
});
</script>










## Abstract {-}

It is very common in cognitive science and psychology to use experimental tasks that involve making a fast choice among a restricted number of alternatives. While a class of choice models, the so-called sequential-sampling models, can give an accurate account of the relationship between the accuracy of the choice and the time it took to respond, it is fairly common to ignore the tradeoff between accuracy and response times and analyze them as if they were independent. The reason for this is that sequential-sampling models are mathematically and computationally difficult to apply. In this notebook, I focus on one influential and relatively simple model that belongs to the class of sequential-sampling models: the linear ballistic accumulator with a drift rate drawn from a normal distribution (restricted to positive values) [@BrownHeathcote2008;@HeathcoteLove2012]. Even though this model has been proved to be  well-suited for tasks that elicit a speeded response (out of any number of possible choices), its hierarchical version is  difficult to implement and fit. First, I discuss the motivation for fitting this model using the Stroop task [@Stroop1935] as a case study. Then, I discuss the challenges in the implementation of the model in (R)Stan [@RStan2018], which might also apply to other hierarchical models with complex likelihood functions. Finally, I show some results that exemplify how the linear ballistic accumulator can be used for examining individual differences.

# Introduction

Examining fast decisions that do not require long deliberation can uncover automatic cognitive processes, since these decisions are taken before more complex strategies have time to kick in. This is a reason for the popularity of tasks that involve making decisions under time pressure in cognitive science and psychology; some examples are the Stroop task [@Stroop1935], Eriksen flanker task [@Eriksen1974], lexical decision tasks, speeded acceptability judgments, and so forth. This type of task is characterized by yielding two measures: the accuracy of the decision made and the time it took to respond. These two measures, however, are not independent since participants can answer faster sacrificing accuracy, or be more cautious decreasing response times but increasing accuracy.

## The problems of ignoring the relationship between response times and accuracy: Individual differences in the Stroop task

I illustrate the problems of ignoring the relationship between response times and accuracy using the Stroop task as a case study. I will focus on a subset of the data of 3337 participants that undertook one variant of the Stroop task, as part of the battery of tasks run in @ManyLabs3. There are many variations of the Stroop task, in this particular case, participants were presented with one word at the center of the screen, which was either "red", "blue", and "green" (`word`) written in either red, blue, or green font (`color`). In one third of the trials the word matched the color of the text ("congruent" condition) and in the rest of the trials it did not match ("incongruent" condition). Participants were instructed to only pay attention to the color, and press `1` if the color of the word was red, `2` if it was blue, and `3` if it was green. The Stroop effect, that is, the difficulty in identifying the color when it mismatches the word in the incongruent condition (e.g., green in color blue) in comparison to a baseline condition, here, the congruent condition (e.g., green in color green) is extremely robust across variations of the task. The most robust finding is that correct responses in the  incongruent condition take longer than incongruent conditions (or any baseline condition). 

To speed up computation, I subset 25 participants of the original dataset (containing 3337 participants). 


```r
set.seed(42)
library(tidyverse)
library(magrittr)
library(stringr)
library(ggplot2)

dstroop <- read_csv("data/StroopCleanSet.csv")  %>%
            # rename unintuitive column names:
            rename(correct= trial_error, 
                    RT = trial_latency,
                    condition = congruent) %>% 
            # Create word and color column names
              mutate(subject = session_id %>% as.factor %>% as.numeric, 
                     word = str_match(trial_name,"(.*?)\\|")[,2] %>%
                            str_to_lower,
                     color = str_match(block_name,"(.*?)\\|")[,2] %>%
                             str_to_lower,
                     cond = if_else(condition == "Congruent", -1, 1))

# To speed thing up, I select the first 25 subjects out
# the 3337 of the original study.
dstroop_25 <- dstroop %>% filter(subject <= 25)

dstroop_25 %>% filter(correct==1) %>% 
            group_by(condition) %>% 
            summarize(mean(RT))
```
I show the estimate of the delay between incongruent and congruent conditions calculated with a log-linear mixed model using the `brms` package [@brms]:


```r
library(brms)
options(mc.cores = parallel::detectCores())

priors <- c(set_prior("normal(0,10)", class = "Intercept"),
           set_prior("normal(0,2)", class = "b"),
           set_prior("normal(0,1)", class = "sigma"),
           set_prior("lkj(2)", class = "cor"),
           set_prior("normal(0,1)", class = "sd"))


m_rt <- brm(RT ~ cond + (cond | subject), filter(dstroop_25, correct==1), 
  family=lognormal(), prior = priors, control = list (adapt_delta = .9))
m_rt
```


```
##  Family: lognormal 
##   Links: mu = identity; sigma = identity 
## Formula: RT ~ cond + (cond | subject) 
##    Data: filter(dstroop_25, correct == 1) (Number of observations: 1515) 
## Samples: 4 chains, each with iter = 2000; warmup = 1000; thin = 1; 
##          total post-warmup samples = 4000
##     ICs: LOO = NA; WAIC = NA; R2 = NA
##  
## Group-Level Effects: 
## ~subject (Number of levels: 25) 
##                     Estimate Est.Error l-95% CI u-95% CI Eff.Sample Rhat
## sd(Intercept)           0.12      0.02     0.09     0.16       1056 1.00
## sd(cond)                0.02      0.01     0.00     0.04       1594 1.00
## cor(Intercept,cond)     0.38      0.36    -0.48     0.90       4000 1.00
## 
## Population-Level Effects: 
##           Estimate Est.Error l-95% CI u-95% CI Eff.Sample Rhat
## Intercept     6.34      0.02     6.29     6.39        765 1.00
## cond          0.04      0.01     0.02     0.06       4000 1.00
## 
## Family Specific Parameters: 
##       Estimate Est.Error l-95% CI u-95% CI Eff.Sample Rhat
## sigma     0.29      0.01     0.28     0.30       4000 1.00
## 
## Samples were drawn using sampling(NUTS). For each parameter, Eff.Sample 
## is a crude measure of effective sample size, and Rhat is the potential 
## scale reduction factor on split chains (at convergence, Rhat = 1).
```

Below I show the effect on milliseconds:


```r
estimates_rt <- posterior_samples(m_rt, pars= c("b") ) %>%
# Stroop effect in milliseconds:
  mutate(stroop_effect_ms = exp(b_Intercept + b_cond) - 
          exp(b_Intercept - b_cond))

mean(estimates_rt$stroop_effect_ms)
```

```
## [1] 44.86322
```

```r
quantile(estimates_rt$stroop_effect_ms,c(0.025,.5,.975))
```

```
##     2.5%      50%    97.5% 
## 25.75824 44.69775 64.39404
```




















































