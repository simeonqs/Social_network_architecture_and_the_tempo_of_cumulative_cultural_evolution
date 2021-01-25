######## FIGURE 2 ###########



# 1. LOADING --------------------------------------------------------------


# packages
if(!require(igraph)){install.packages('igraph'); library(igraph)}
if(!require(GGally)){install.packages('GGally'); library(GGally)}
if(!require(sna)){install.packages('sna'); library(sna)}
if(!require(viridis)){install.packages('viridis'); library(viridis)}
if(!require(gridExtra)){install.packages('gridExtra'); library(gridExtra)}
if(!require(tidyverse)){install.packages('tidyverse'); library(tidyverse)}
if(!require(survival)){install.packages('survival'); library(survival)}
if(!require(survminer)){install.packages('survminer'); library(survminer)}
if(!require(ggridges)){install.packages('ggridges'); library(ggridges)}
if(!require(grid)){install.packages('grid'); library(grid)}


# set colour palette
colviridis = viridis(7, begin = 0, end = 1, direction = 1, option = 'viridis')
alphaviridis = seq(from=0.1, to=1, by=0.1)
matviridis = matrix(NA, nrow = length(colviridis), ncol=length(alphaviridis))
rownames(matviridis) = c("FULL","SMALLWORLD","DEGREE","CLUSTERED","MODULAR","MOD_CLUST","MULTILEVEL")
for(i in 1:length(colviridis)){
  for(j in 1:length(alphaviridis)){
    matviridis[i,j] = alpha(colviridis[i], alphaviridis[j])
  }
}  

# reset full to greyscale
for(j in 1:length(alphaviridis)){
  matviridis[1,j]  = alpha('grey75', alphaviridis[j])
}




# 2. INPUT DATA -----------------------------------------------------------


# LOAD AND PICK DATA FOR MODEL1 OR MODEL 2 

load("/Users/maucantor/Documents/Google Drive/Postdoctoral/Articles/MLS_CulturalEvolution/Code/ABM/ABM_python/ABM_Final/df_TTC_m1.Rda")
load("/Users/maucantor/Documents/Google Drive/Postdoctoral/Articles/MLS_CulturalEvolution/Code/ABM/ABM_python/ABM_Final/df_TTC_m2.Rda")
load("/Users/maucantor/Documents/Google Drive/Postdoctoral/Articles/MLS_CulturalEvolution/Code/ABM/ABM_python/ABM_Final/df_TTD_m2.Rda")
load("/Users/maucantor/Documents/Google Drive/Postdoctoral/Articles/MLS_CulturalEvolution/Code/ABM/ABM_python/ABM_Final/df_TTC_TTD_m2.Rda")







# 3. RIDGES ---------------------------------------------------------------


### N= 64 ####
## TTC_TTD_m2 ##

cu = df_TTC_TTD_m2

# Set small population size N=64, low degree K=12, include fully connected networks
cu = cu[which(cu$degree==12 | cu$degree=='full'),]
cu = cu[which(cu$pop_size==64),]
table(cu$degree, cu$graph, cu$pop_size)

# reorder nets by inverse structural complexity
cu$graph = factor(cu$graph, levels = c("multilevel", 'modularclustered',"modular", "clustered", "small", "degree",  'full'))

# PLOT
rids_ttcttd =cu %>%
ggplot(aes(x = log(delta_TTC_TTD+1), y = graph)) +
     geom_density_ridges(aes(fill = graph), scale = 3, size = 0.3, jittered_points = F,
                         bandwidth = 0.18#,
                        #quantile_lines = T, quantiles = 2
                        ) +
  scale_fill_manual(values = as.character(matviridis[7:1,6])) +  
  scale_x_continuous(limits = c(0, 3.5), expand = c(0, 0), labels = c('1', '2.5', '5',  '10'), 
                     breaks = c(0, 0.92, 1.61, 2.31))+
  labs(x = '', y = '') +
  theme_ridges() + theme(legend.position = "none",
                         axis.title.x = element_text(vjust = 0.5, hjust = 0.5),
                         axis.text.x = element_blank(),
                         axis.text.y = element_blank(),
                         plot.margin=unit(c(0.05,0.01,0.05,0.01),"cm")) +
  ggtitle(label='C')





## TTC_m2 ##

cu = df_TTC_m2

# Set small population size N=64, low degree K=12, include fully connected networks
cu = cu[which(cu$degree==12 | cu$degree=='full'),]
cu = cu[which(cu$pop_size==64),]
table(cu$degree, cu$graph, cu$pop_size)

# reorder nets by inverse structural complexity
cu$graph = factor(cu$graph, levels = c("multilevel", 'modularclustered',"modular", "clustered", "small", "degree",  'full'))

# PLOT
rids_ttcm2 =
cu %>%
ggplot(aes(x = log(epoch+1), y = graph)) +
     geom_density_ridges(aes(fill = graph), scale = 3, size = 0.3, jittered_points = F,
                         bandwidth = 0.18
                        ) +
  scale_fill_manual(values = as.character(matviridis[7:1,6])) +  
  scale_x_continuous(limits = c(-1, 8), expand = c(0, 0), labels = c('1', '10',  '100', '1000'), breaks = c(0, 2.31, 4.608, 6.908)) +
  labs(x = '', y = '') +
  theme_ridges() + theme(legend.position = "none",
                         axis.title.x = element_text(vjust = 0.5, hjust = 0.5),
                         axis.text.x = element_blank(),
                         axis.text.y = element_blank(),
                         plot.margin=unit(c(0.05,0.01,0.05,0.01),"cm")) +
  ggtitle(label='B')




## TTC_m1 ##

cu = df_TTC_m1

# Set small population size N=64, low degree K=12, include fully connected networks
cu = cu[which(cu$degree==12 | cu$degree=='full'),]
cu = cu[which(cu$pop_size==64),]
table(cu$degree, cu$graph, cu$pop_size)

# reorder nets by inverse structural complexity
cu$graph = factor(cu$graph, levels = c("multilevel", 'modularclustered',"modular", "clustered", "small", "degree",  'full'), labels = c('Multilevel', 'Modular lattice', 'Modular','Lattice','Small-world','Random', "Fully connected"))


# PLOT
rids_ttcm1 = cu %>%
ggplot(aes(x = log(epoch+1), y = graph)) +
     geom_density_ridges(aes(fill = graph), scale = 3, size = 0.3, jittered_points = F,
                         bandwidth = 0.18
                        ) +
    scale_fill_manual(values = as.character(matviridis[7:1,6])) +  
  scale_x_continuous(limits = c(-1, 8), expand = c(0, 0), labels = c('1', '10',  '100', '1000'), breaks = c(0, 2.31, 4.608, 6.908)) +
  labs(x = '', y = 'Network size N=64') +
  theme_ridges() + theme(legend.position = "none",
                         axis.title.x = element_text(vjust = 0.5, hjust = 0.5),
                         axis.title.y = element_text(vjust = 0.5, hjust = 0.5),
                                                  axis.text.x = element_blank(),
                         plot.margin=unit(c(0.05,0.01,0.05,0.01),"cm")) +
  ggtitle(label='A')










#### N= 144 ####

## TTC_TTD_m2 ##

cu = df_TTC_TTD_m2

cu = cu[which(cu$degree==12 | cu$degree=='full'),]
cu = cu[which(cu$pop_size==144),]
table(cu$degree, cu$graph, cu$pop_size)

# reorder nets by inverse structural complexity
cu$graph = factor(cu$graph, levels = c("multilevel", 'modularclustered',"modular", "clustered", "small", "degree",  'full'))

# PLOT
rids_ttcttd_144 =cu %>%
ggplot(aes(x = log(delta_TTC_TTD+1), y = graph)) +
     geom_density_ridges(aes(fill = graph), scale = 3, size = 0.3, jittered_points = F,
                         bandwidth = 0.18
                        ) +
  scale_fill_manual(values = as.character(matviridis[7:1,6])) +  
  scale_x_continuous(limits = c(0, 3.5), expand = c(0, 0), labels = c('1', '2.5', '5',  '10'), 
                     breaks = c(0, 0.92, 1.61, 2.31))+
  labs(x = '', y = '') +
  theme_ridges() + theme(legend.position = "none",
                         axis.title.x = element_text(vjust = 0.5, hjust = 0.5),
                         axis.text.x = element_blank(),
                         axis.text.y = element_blank(),
                         plot.margin=unit(c(0.05,0.01,0.05,0.01),"cm")) +
  ggtitle(label='')





## TTC_m2 ##

cu = df_TTC_m2

cu = cu[which(cu$degree==12 | cu$degree=='full'),]
cu = cu[which(cu$pop_size==144),]
table(cu$degree, cu$graph, cu$pop_size)

# reorder nets by inverse structural complexity
cu$graph = factor(cu$graph, levels = c("multilevel", 'modularclustered',"modular", "clustered", "small", "degree",  'full'))


# PLOT
rids_ttcm2_144 =
cu %>%
ggplot(aes(x = log(epoch+1), y = graph)) +
     geom_density_ridges(aes(fill = graph), scale = 3, size = 0.3, jittered_points = F,
                         bandwidth = 0.18
                        ) +
  scale_fill_manual(values = as.character(matviridis[7:1,6])) +  
  scale_x_continuous(limits = c(-1, 8), expand = c(0, 0), labels = c('1', '10',  '100', '1000'), breaks = c(0, 2.31, 4.608, 6.908)) +
  labs(x = '', y = '') +
  scale_y_discrete(labels = c(' ', ' ', ' ', ' ', ' ', ' ', ' ')) +
  theme_ridges() + theme(legend.position = "none",
                         axis.title.x = element_text(vjust = 0.5, hjust = 0.5),
                         axis.text.x = element_blank(),
                         axis.text.y = element_blank(),
                         plot.margin=unit(c(0.05,0.01,0.05,0.01),"cm")) +
  ggtitle(label='')




## TTC_m1 ##

cu = df_TTC_m1

cu = cu[which(cu$degree==12 | cu$degree=='full'),]
cu = cu[which(cu$pop_size==144),]
table(cu$degree, cu$graph, cu$pop_size)

# reorder nets by inverse structural complexity
cu$graph = factor(cu$graph, levels = c("multilevel", 'modularclustered',"modular", "clustered", "small", "degree",  'full'), labels = c('Multilevel', 'Modular lattice', 'Modular','Lattice','Small-world','Random', "Fully connected"))


# PLOT
rids_ttcm1_144 = cu %>%
ggplot(aes(x = log(epoch+1), y = graph)) +
     geom_density_ridges(aes(fill = graph), scale = 3, size = 0.3, jittered_points = F,
                         bandwidth = 0.18
                        ) +
  scale_fill_manual(values = as.character(matviridis[7:1,6])) +  
  scale_x_continuous(limits = c(-1, 8), expand = c(0, 0), labels = c('1', '10',  '100', '1000'), breaks = c(0, 2.31, 4.608, 6.908)) +
  labs(x = '', y = 'Network size N=144') +
  theme_ridges() + theme(legend.position = "none",
                         axis.title.x = element_text(vjust = 0.5, hjust = 0.5),
                         axis.title.y = element_text(vjust = 0.5, hjust = 0.5),
                         axis.text.x = element_blank(),
                         plot.margin=unit(c(0.05,0.01,0.05,0.01),"cm")) +
  ggtitle(label='')






#### N=324 ####

## TTC_TTD_m2 ##

cu = df_TTC_TTD_m2

# Set small population size N=64, low degree K=12, include fully connected networks
cu = cu[which(cu$degree==12 | cu$degree=='full'),]
cu = cu[which(cu$pop_size==324),]
table(cu$degree, cu$graph, cu$pop_size)

# reorder nets by inverse structural complexity
cu$graph = factor(cu$graph, levels = c("multilevel", 'modularclustered',"modular", "clustered", "small", "degree",  'full'))


# PLOT
rids_ttcttd_324 =cu %>%
ggplot(aes(x = log(delta_TTC_TTD+1), y = graph)) +
     geom_density_ridges(aes(fill = graph), scale = 3, size = 0.3, jittered_points = F,
                         bandwidth = 0.18
                        ) +
  scale_fill_manual(values = as.character(matviridis[7:1,6])) +  
  scale_x_continuous(limits = c(0, 3.5), expand = c(0, 0), labels = c('1', '2.5', '5',  '10'), 
                     breaks = c(0, 0.92, 1.61, 2.31))+
  labs(x = 'Time from recombination \n to diffusion (epoch)', y = '') +
  theme_ridges() + theme(legend.position = "none",
                         axis.title.x = element_text(vjust = 0.5, hjust = 0.5),
                         axis.text.y = element_blank(),
                         plot.margin=unit(c(0.05,0.01,0.05,0.01),"cm")) +
  ggtitle(label='')





## TTC_m2 ##

cu = df_TTC_m2

# Set small population size N=64, low degree K=12, include fully connected networks
cu = cu[which(cu$degree==12 | cu$degree=='full'),]
cu = cu[which(cu$pop_size==324),]
table(cu$degree, cu$graph, cu$pop_size)

# reorder nets by inverse structural complexity
cu$graph = factor(cu$graph, levels = c("multilevel", 'modularclustered',"modular", "clustered", "small", "degree",  'full'))

# PLOT
rids_ttcm2_324 =
cu %>%
ggplot(aes(x = log(epoch+1), y = graph)) +
     geom_density_ridges(aes(fill = graph), scale = 3, size = 0.3, jittered_points = F,
                         bandwidth = 0.18
                        ) +
  scale_fill_manual(values = as.character(matviridis[7:1,6])) +  
  scale_x_continuous(limits = c(-1, 8), expand = c(0, 0), labels = c('1', '10',  '100', '1000'), breaks = c(0, 2.31, 4.608, 6.908)) +
  labs(x = 'Time to recombination with \n one-to-one diffusion (epoch)', y = '') +
  theme_ridges() + theme(legend.position = "none",
                         axis.title.x = element_text(vjust = 0.5, hjust = 0.5),
                         axis.text.y = element_blank(),
                         plot.margin=unit(c(0.05,0.01,0.05,0.01),"cm")) +
  ggtitle(label='')




## TTC_m1 ##

cu = df_TTC_m1

# Set small population size N=64, low degree K=12, include fully connected networks
cu = cu[which(cu$degree==12 | cu$degree=='full'),]
cu = cu[which(cu$pop_size==324),]
table(cu$degree, cu$graph, cu$pop_size)

# reorder nets by inverse structural complexity
cu$graph = factor(cu$graph, levels = c("multilevel", 'modularclustered',"modular", "clustered", "small", "degree",  'full'), labels = c('Multilevel', 'Modular lattice', 'Modular','Lattice','Small-world','Random', "Fully connected"))


# PLOT
rids_ttcm1_324 = cu %>%
ggplot(aes(x = log(epoch+1), y = graph)) +
     geom_density_ridges(aes(fill = graph), scale = 3, size = 0.3, jittered_points = F,
                         bandwidth = 0.18
                        ) +
  scale_fill_manual(values = as.character(matviridis[7:1,6])) +  
  scale_x_continuous(limits = c(-1, 8), expand = c(0, 0), labels = c('1', '10',  '100', '1000'), breaks = c(0, 2.31, 4.608, 6.908)) +
  labs(x = 'Time to recombination with \n one-to-many diffusion (epoch)', y = 'Network size N=324') +
  theme_ridges() + theme(legend.position = "none",
                         axis.title.x = element_text(vjust = 0.5, hjust = 0.5),
                         axis.title.y = element_text(vjust = 0.5, hjust = 0.5),
                         plot.margin=unit(c(0.05,0.01,0.05,0.01),"cm")) +
  ggtitle(label='')







# 4. PLOTTING -------------------------------------------------------------

lay <- rbind(c(1,1,1,1,2,2,2,3,3,3),
             c(4,4,4,4,5,5,5,6,6,6),
             c(7,7,7,7,8,8,8,9,9,9))
grid.arrange(rids_ttcm1, rids_ttcm2, rids_ttcttd,
             rids_ttcm1_144, rids_ttcm2_144, rids_ttcttd_144,
            rids_ttcm1_324, rids_ttcm2_324, rids_ttcttd_324,
             layout_matrix = lay)