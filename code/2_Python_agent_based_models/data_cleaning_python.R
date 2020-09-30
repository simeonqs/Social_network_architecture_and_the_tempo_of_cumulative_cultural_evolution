library(tidyverse)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#create df_TTC_m1
{
  df_1 = read.csv("m1/TTC_clustered.csv")
  df_2 = read.csv("m1/TTC_degree.csv")
  df_3 = read.csv("m1/TTC_full.csv")
  df_4 = read.csv("m1/TTC_modularclustered.csv")
  df_5 = read.csv("m1/TTC_modular.csv")
  df_6 = read.csv("m1/TTC_multilevel.csv")
  df_7 = read.csv("m1/TTC_smallworld.csv")
  
  df = rbind(df_1,df_2,df_3,df_4,df_5,df_6,df_7)
  df = df[df$discovery %in% c(7,14),]
  df %>% group_by(graph_type) %>% summarize(n()) %>% print(n=Inf)
  df$graph_type = as.character(df$graph_type)
  df = df %>% mutate(graph = strsplit(graph_type, '_'))
  df$graph = sapply(df$graph, `[`, 2)
  df$graph = as.factor(df$graph)
  df = df %>% mutate(degree = strsplit(graph_type, '_'))
  df$degree = sapply(df$degree, function(x) x[length(x)]) %>% str_remove("d")
  df = df[df$degree != "16",]
  df$degree = factor(df$degree, levels=c("8","12", "18","24", "30", "full"))
  df$epoch = df$epoch+1
  df_TTC_m1 = df
  save(df_TTC_m1, file="df_TTC_m1.Rda")
}

#create df_TTC_m2
{
  df_1 = read.csv("m2/TTC_clustered.csv")
  df_2 = read.csv("m2/TTC_degree.csv")
  df_3 = read.csv("m2/TTC_full.csv")
  df_4 = read.csv("m2/TTC_modularclustered.csv")
  df_5 = read.csv("m2/TTC_modular.csv")
  df_6 = read.csv("m2/TTC_multilevel.csv")
  df_7 = read.csv("m2/TTC_smallworld.csv")
  
  df = rbind(df_1,df_2,df_3,df_4,df_5,df_6,df_7)
  df = df[df$discovery %in% c(7,14),]
  df = df %>% group_by(graph_type,sim) %>% arrange(timestep) %>% slice_head(n = 1)
  df %>% group_by(graph_type) %>% summarize(n()) %>% print(n=Inf)
  df$graph_type = as.character(df$graph_type)
  df = df %>% mutate(graph = strsplit(graph_type, '_'))
  df$graph = sapply(df$graph, `[`, 2)
  df$graph = as.factor(df$graph)
  df = df %>% mutate(degree = strsplit(graph_type, '_'))
  df$degree = sapply(df$degree, function(x) x[length(x)]) %>% str_remove("d")
  df = df[df$degree != "16",]
  df$degree = factor(df$degree, levels=c("8","12", "18","24", "30", "full"))
  df$epoch = df$epoch+1
  df_TTC_m2 = df
  save(df_TTC_m2, file="df_TTC_m2.Rda")
}

#create_df_TTD_m2
{
  df_1 = read.csv("m2/TTD_clustered.csv")
  df_2 = read.csv("m2/TTD_degree.csv")
  df_3 = read.csv("m2/TTD_full.csv")
  df_4 = read.csv("m2/TTD_modularclustered.csv")
  df_5 = read.csv("m2/TTD_modular.csv")
  df_6 = read.csv("m2/TTD_multilevel.csv")
  df_7 = read.csv("m2/TTD_smallworld.csv")
  df = rbind(df_1,df_2,df_3,df_4,df_5,df_6,df_7)
  df$graph_type = as.character(df$graph_type)
  df = df %>% mutate(graph = strsplit(graph_type, '_'))
  df$graph = sapply(df$graph, `[`, 2)
  df$graph = as.factor(df$graph)
  df = df %>% mutate(degree = strsplit(graph_type, '_'))
  df$degree = sapply(df$degree, function(x) x[length(x)]) %>% str_remove("d")
  df = df[df$degree != "16",]
  df$degree = factor(df$degree, levels=c("8","12","18","24","30","full"))
  summary(df)
  df$epoch = df$epoch+1
  df_TTD_m2 = df
  save(df_TTD_m2, file="df_TTD_m2.Rda")
}

#create TTD-TTC
{
  df1 = df_TTC_m2 %>% select(sim,graph_type,pop_size,graph,degree,epoch) %>% rename(epoch_TTC=epoch)
  df2 = df_TTD_m2 %>% select(sim,graph_type,pop_size,graph,degree,epoch) %>% rename(epoch_TTD=epoch)
  df_TTC_TTD_m2 = left_join(df1,df2)
  df_TTC_TTD_m2 = df_TTC_TTD_m2 %>% mutate(delta_TTC_TTD=epoch_TTD-epoch_TTC)
  save(df_TTC_TTD_m2,file="df_TTC_TTD_m2.Rda")
}


#prep df with network information
{
  df1 = read.csv("network_properties/DETAILS_d8_N_64.csv")
  df1$pop_size=64
  df1$degree=8
  df2 = read.csv("network_properties/DETAILS_d8_N_144.csv")
  df2$pop_size=144
  df2$degree=8
  df3 = read.csv("network_properties/DETAILS_d8_N_324.csv")
  df3$pop_size=324
  df3$degree=8
  df4 = read.csv("network_properties/DETAILS_d12_N_64.csv")
  df4$pop_size=64
  df4$degree=12
  df5 = read.csv("network_properties/DETAILS_d12_N_144.csv")
  df5$pop_size=144
  df5$degree=12
  df6 = read.csv("network_properties/DETAILS_d12_N_324.csv")
  df6$pop_size=324
  df6$degree=12
  df7 = read.csv("network_properties/DETAILS_d18_N_144.csv")
  df7$pop_size=144
  df7$degree=18
  df8 = read.csv("network_properties/DETAILS_d18_N_324.csv")
  df8$pop_size=324
  df8$degree=18
  df9 = read.csv("network_properties/DETAILS_d24_N_144.csv")
  df9$pop_size=144
  df9$degree=24
  df10 = read.csv("network_properties/DETAILS_d24_N_324.csv")
  df10$pop_size=324
  df10$degree=24
  df11 = read.csv("network_properties/DETAILS_d30_N_324.csv")
  df11$pop_size=324
  df11$degree=30
  
  df = rbind(df1,df2,df3,df4,df5,df6,df7,df8,df9,df10,df11) %>% 
    select(-X) %>% 
    rename(graph = type)
  levels(df$graph) <- c("lattice","random","full","modular_lattice","modular","multilevel","small_world")
  summary(df)
  df$degree = as.character(df$degree)
  df = df %>% mutate(degree = if_else(graph=="full","full",as.character(degree)))
  head(df)
  df$degree = factor(df$degree, levels=c("8","12", "18","24", "30", "full"))
  df = df %>% unique()
}


#load data
load("df_TTC_m1.Rda")
df_TTC_m1 = df_TTC_m1 %>% droplevels()
load("df_TTC_m2.Rda")
df_TTC_m2 = df_TTC_m2 %>% droplevels()
load("df_TTC_TTD_m2.Rda")
df_TTC_TTD_m2 = df_TTC_TTD_m2 %>% droplevels()

#rename to those used in figures
levels(df_TTC_m1$graph) <- c("lattice","random","full","modular_lattice","modular","multilevel","small_world") 
levels(df_TTC_m2$graph) <- c("lattice","random","full","modular_lattice","modular","multilevel","small_world") 
levels(df_TTC_TTD_m2$graph) <- c("lattice","random","full","modular_lattice","modular","multilevel","small_world") 

#join datasets
df_TTC_m1 = left_join(df_TTC_m1,df)
df_TTC_m2 = left_join(df_TTC_m2,df)
df_TTC_TTD_m2 = left_join(df_TTC_TTD_m2,df)

summary(df_TTC_m1)
#reorder to decreasing complexity
df_TTC_m1$graph = factor(df_TTC_m1$graph, 
                         levels = c("random","small_world","lattice","modular","modular_lattice","multilevel","full"))
df_TTC_m2$graph = factor(df_TTC_m2$graph, 
                         levels = c("random","small_world","lattice","modular","modular_lattice","multilevel","full"))
df_TTC_TTD_m2$graph = factor(df_TTC_TTD_m2$graph, 
                             levels = c("random","small_world","lattice","modular","modular_lattice","multilevel","full"))

save(df_TTC_m1,file="../data/df_TTC_m1.Rda")
save(df_TTC_m2,file="../data/df_TTC_m2.Rda")
save(df_TTC_TTD_m2,file="../data/df_TTC_TTD_m2.Rda")

