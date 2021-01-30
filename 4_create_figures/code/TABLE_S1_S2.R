library(tidyverse)
library(lme4)
library(stargazer)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

load("../data/df_TTC_m1.Rda")
summary(df_TTC_m1)
df_TTC_m1 %>% filter(graph!="full") %>% group_by(pop_size) %>% summarize(sd(epoch)/mean(epoch))
df_TTC_m1 %>% filter(graph!="full") %>% group_by(pop_size) %>% summarize(IQR(epoch)/(quantile(epoch,.25)+quantile(epoch,.75)))

#### Table S1 ####
#M1 Full network to others comparison
df = df_TTC_m1
df$pop_size = as.factor(df$pop_size)
M1_TTC_glm = glm(log(epoch+1) ~ graph*pop_size, family=gaussian(link="log"),data=df)

#M2 Full network to others comparison
load("../data/df_TTC_m2.Rda")
df = df_TTC_m2
df$pop_size = as.factor(df$pop_size)
M2_TTC_glm = glm(log(epoch+1) ~ graph*pop_size, family=gaussian(link="log"),data=df)
stargazer(M1_TTC_glm,M2_TTC_glm,apply.coef= function(x) exp(x), summary = FALSE, digits=3, type="html", intercept.bottom = FALSE,single.row = TRUE, font.size = "small",report="vcstp*",out="../output/tableS1.html")

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
stargazer(M1_TTC_glm,M2_TTC_glm,apply.coef= function(x) exp(x),summary = FALSE, digits=3, type="html", intercept.bottom = FALSE,single.row = TRUE, font.size = "small",report="vcstp*",out="../output/tableS2.html")