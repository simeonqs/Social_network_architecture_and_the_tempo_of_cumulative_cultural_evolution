# packages
if(!require(tidyverse)){install.packages('tidyverse'); library(tidyverse)}
if(!require(lme4)){install.packages('lme4'); library(lme4)}
if(!require(stargazer)){install.packages('stargazer'); library(stargazer)}

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


#### Table S1 ####
load("../../2_Python_agent_based_models/data/df_TTC_m1.Rda")
df = df_TTC_m1
df$pop_size = as.factor(df$pop_size)
M1_TTC_glm = glm(log(epoch+1) ~ graph*pop_size, family=gaussian(link="log"),data=df)

load("../../2_Python_agent_based_models/data/df_TTC_m2.Rda")
df = df_TTC_m2
df$pop_size = as.factor(df$pop_size)
M2_TTC_glm = glm(log(epoch+1) ~ graph*pop_size, family=gaussian(link="log"),data=df)
stargazer(M1_TTC_glm,M2_TTC_glm,apply.coef= function(x) exp(x), summary = FALSE, digits=3, type="html", intercept.bottom = FALSE,single.row = TRUE, font.size = "small",report="vcstp*",out="../results/tableS1.html")

#### Table S2 ####
df = df_TTC_m1 %>% filter(graph!="full") %>% droplevels()
df$pop_size = as.factor(df$pop_size)
df$degree = as.factor(df$degree)
M1_TTC_glm = glm(log(df$epoch+1) ~ graph*pop_size*degree, family=gaussian(link="log"),data=df)

#remove fully connected
df = df_TTC_m2 %>% filter(graph!="full") %>% droplevels()
df$pop_size = as.factor(df$pop_size)
df$degree = as.factor(df$degree)
M2_TTC_glm = glm(log(df$epoch+1) ~ graph*pop_size*degree, family=gaussian(link="log"),data=df)
stargazer(M1_TTC_glm,M2_TTC_glm,apply.coef= function(x) exp(x),summary = FALSE, digits=3, type="html", intercept.bottom = FALSE,single.row = TRUE, font.size = "small",report="vcstp*",out="../results/tableS2.html")
