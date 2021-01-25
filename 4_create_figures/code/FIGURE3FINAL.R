######## FIGURE 2 ###########



# 1. LOADING --------------------------------------------------------------


# packages
if(!require(viridis)){install.packages('viridis'); library(viridis)}
if(!require(tidyverse)){install.packages('tidyverse'); library(tidyverse)}
if(!require(survival)){install.packages('survival'); library(survival)}
if(!require(survminer)){install.packages('survminer'); library(survminer)}
if(!require(ggpubr)){install.packages('ggpubr'); library(ggpubr)}



# 2. INPUT DATA -----------------------------------------------------------


# LOAD 

load("/Users/maucantor/Documents/Google Drive/Postdoctoral/Articles/MLS_CulturalEvolution/Code/ABM/ABM_python/ABM_Final/df_TTC_m1.Rda")

df_TTC_m1 %>% filter(graph!="full") %>% group_by(pop_size) %>% summarize(sd(epoch)/mean(epoch))
df_TTC_m1 %>% filter(graph!="full") %>% group_by(pop_size) %>% summarize(IQR(epoch)/(quantile(epoch,.25)+quantile(epoch,.75)))
summary(df_TTC_m1$epoch)

#M1 by degree and pop size

plot(df$graph,log(df$epoch))
plot(density(df$epoch[df$graph=="random"] %>% log()))
{
  df = df_TTC_m1 %>% filter(graph!="full",pop_size==64,degree=="8") %>% droplevels()
  df$pop_size = as.factor(df$pop_size)
  fit.coxph = glm(log(epoch+1)~graph, family=gaussian(link="log"), data=df)
  summary(fit.coxph)
  r1= exp(fit.coxph$coefficients)
  r1
  
  df = df_TTC_m1 %>% filter(graph!="full",pop_size==64,degree=="12") %>% droplevels()
  df$pop_size = as.factor(df$pop_size)
  fit.coxph = glm(log(epoch+1)~graph, data=df, family=gaussian(link="log"))
  r2= exp(fit.coxph$coefficients)
  r2
  
  df = df_TTC_m1 %>% filter(graph!="full",pop_size==144,degree=="8") %>% droplevels()
  df$pop_size = as.factor(df$pop_size)
  fit.coxph = glm(log(epoch+1)~graph, family=gaussian(link="log"), data=df)
  r3= exp(fit.coxph$coefficients)
  
  df = df_TTC_m1 %>% filter(graph!="full",pop_size==144,degree=="12") %>% droplevels()
  df$pop_size = as.factor(df$pop_size)
  fit.coxph = glm(log(epoch+1)~graph, family=gaussian(link="log"), data=df)
  r4= exp(fit.coxph$coefficients)
  
  df = df_TTC_m1 %>% filter(graph!="full",pop_size==144,degree=="18") %>% droplevels()
  df$pop_size = as.factor(df$pop_size)
  fit.coxph = glm(log(epoch+1)~graph, family=gaussian(link="log"), data=df)
  r5 = exp(fit.coxph$coefficients)
  
  df = df_TTC_m1 %>% filter(graph!="full",pop_size==144,degree=="24") %>% droplevels()
  df$pop_size = as.factor(df$pop_size)
  fit.coxph = glm(log(epoch+1)~graph, family=gaussian(link="log"), data=df)
  r6 = exp(fit.coxph$coefficients)
  
  df = df_TTC_m1 %>% filter(graph!="full",pop_size==324,degree=="8") %>% droplevels()
  df$pop_size = as.factor(df$pop_size)
  fit.coxph = glm(log(epoch+1)~graph, family=gaussian(link="log"), data=df)
  r7= exp(fit.coxph$coefficients)
  
  df = df_TTC_m1 %>% filter(graph!="full",pop_size==324,degree=="12") %>% droplevels()
  df$pop_size = as.factor(df$pop_size)
  fit.coxph = glm(log(epoch+1)~graph, family=gaussian(link="log"), data=df)
  r8= exp(fit.coxph$coefficients)
  
  df = df_TTC_m1 %>% filter(graph!="full",pop_size==324,degree=="18") %>% droplevels()
  df$pop_size = as.factor(df$pop_size)
  fit.coxph = glm(log(epoch+1)~graph, family=gaussian(link="log"), data=df)
  r9 = exp(fit.coxph$coefficients)
  
  df = df_TTC_m1 %>% filter(graph!="full",pop_size==324,degree=="24") %>% droplevels()
  df$pop_size = as.factor(df$pop_size)
  fit.coxph = glm(log(epoch+1)~graph, family=gaussian(link="log"), data=df)
  r10 = exp(fit.coxph$coefficients)
  
  df = df_TTC_m1 %>% filter(graph!="full",pop_size==324,degree=="30") %>% droplevels()
  df$pop_size = as.factor(df$pop_size)
  fit.coxph = glm(log(epoch+1)~graph, family=gaussian(link="log"), data=df)
  r11 = exp(fit.coxph$coefficients)
  
  df_degree_popsize = as.data.frame(rbind(r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11),row.names = F)
  pop_size = c(64,64,144,144,144,144,324,324,324,324,324)
  degree = c(8,12,8,12,18,24,8,12,18,24,30)
  df_degree_popsize$N = pop_size
  df_degree_popsize$K = degree
  
  df_degree_popsize$graphrandom = 1
  df_degree_popsize = df_degree_popsize %>% mutate(row=paste(paste("N=",N,sep=""),paste("K=",K,sep="")))
  df_output = df_degree_popsize %>% select(row,K,N, graphrandom,graphsmall_world,graphlattice,graphmodular,graphmodular_lattice,graphmultilevel)
  df_output = df_output %>% pivot_longer(cols=c(graphrandom,graphsmall_world,graphlattice,graphmodular,graphmodular_lattice,graphmultilevel))
  df_output$name = factor(df_output$name,levels = c("graphrandom","graphsmall_world","graphlattice","graphmodular","graphmodular_lattice","graphmultilevel"))
  levels(df_output$name) <- c("Random","Small world","Lattice","Modular","Modular lattice","Multilevel") 
  
  N.labs <- c("Network size N=64","Network size N=144","Network size N=324")
  names(N.labs) <- c(64,144,324)
  m1_heatmap = ggplot(data=df_output,aes(x=name,y=as.factor(K)))+
    facet_wrap(~N,scales="free",labeller = labeller(N=N.labs))+
    geom_tile(aes(fill=value),size=.2,color="black")+
    geom_text(aes(label = sprintf("%0.3f", round(value, digits = 3)),color=value<.89),size=3,show.legend = F) +
    scale_color_manual(values=c("black","white"))+
    scale_fill_gradient2(high='red',low='darkblue',mid='white',midpoint=1)+
    labs(x="",y="",fill=expression(paste("exp(", beta, ")")))+
    theme_minimal()+
    theme(axis.text.x = element_blank())
  m1_heatmap
  ggsave("cox_ph_output/arch_performance_heatmap_GLM.pdf",height=7,width=17.8,units="cm",scale=1.5)
}


#M2 by degree and pop size
{
  load("df_TTC_m2.Rda")
  df = df_TTC_m2 %>% filter(graph!="full",pop_size==64,degree=="8") %>% droplevels()
  df$pop_size = as.factor(df$pop_size)
  fit.coxph = glm(log(epoch+1)~graph, family=gaussian(link="log"), data=df)
  r1= exp(fit.coxph$coefficients)
  
  df = df_TTC_m2 %>% filter(graph!="full",pop_size==64,degree=="12") %>% droplevels()
  df$pop_size = as.factor(df$pop_size)
  fit.coxph = glm(log(epoch+1)~graph, family=gaussian(link="log"), data=df)
  r2= exp(fit.coxph$coefficients)
  
  df = df_TTC_m2 %>% filter(graph!="full",pop_size==144,degree=="8") %>% droplevels()
  df$pop_size = as.factor(df$pop_size)
  fit.coxph = glm(log(epoch+1)~graph, family=gaussian(link="log"), data=df)
  r3= exp(fit.coxph$coefficients)
  
  df = df_TTC_m2 %>% filter(graph!="full",pop_size==144,degree=="12") %>% droplevels()
  df$pop_size = as.factor(df$pop_size)
  fit.coxph = glm(log(epoch+1)~graph, family=gaussian(link="log"), data=df)
  r4= exp(fit.coxph$coefficients)
  
  df = df_TTC_m2 %>% filter(graph!="full",pop_size==144,degree=="18") %>% droplevels()
  df$pop_size = as.factor(df$pop_size)
  fit.coxph = glm(log(epoch+1)~graph, family=gaussian(link="log"), data=df)
  r5 = exp(fit.coxph$coefficients)
  
  df = df_TTC_m2 %>% filter(graph!="full",pop_size==144,degree=="24") %>% droplevels()
  df$pop_size = as.factor(df$pop_size)
  fit.coxph = glm(log(epoch+1)~graph, family=gaussian(link="log"), data=df)
  r6 = exp(fit.coxph$coefficients)
  
  df = df_TTC_m2 %>% filter(graph!="full",pop_size==324,degree=="8") %>% droplevels()
  df$pop_size = as.factor(df$pop_size)
  fit.coxph = glm(log(epoch+1)~graph, family=gaussian(link="log"), data=df)
  r7= exp(fit.coxph$coefficients)
  
  df = df_TTC_m2 %>% filter(graph!="full",pop_size==324,degree=="12") %>% droplevels()
  df$pop_size = as.factor(df$pop_size)
  fit.coxph = glm(log(epoch+1)~graph, family=gaussian(link="log"), data=df)
  r8= exp(fit.coxph$coefficients)
  
  df = df_TTC_m2 %>% filter(graph!="full",pop_size==324,degree=="18") %>% droplevels()
  fit.coxph = glm(log(epoch+1)~graph, family=gaussian(link="log"), data=df)
  r9 = exp(fit.coxph$coefficients)
  
  df = df_TTC_m2 %>% filter(graph!="full",pop_size==324,degree=="24") %>% droplevels()
  df$pop_size = as.factor(df$pop_size)
  fit.coxph = glm(log(epoch+1)~graph, family=gaussian(link="log"), data=df)
  r10 = exp(fit.coxph$coefficients)
  
  df = df_TTC_m2 %>% filter(graph!="full",pop_size==324,degree=="30") %>% droplevels()
  df$pop_size = as.factor(df$pop_size)
  fit.coxph = glm(log(epoch+1)~graph, family=gaussian(link="log"), data=df)
  r11 = exp(fit.coxph$coefficients)
  
  df_degree_popsize = as.data.frame(rbind(r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11),row.names = F)
  pop_size = c(64,64,144,144,144,144,324,324,324,324,324)
  degree = c(8,12,8,12,18,24,8,12,18,24,30)
  df_degree_popsize$N = pop_size
  df_degree_popsize$K = degree
  df_degree_popsize$graphrandom = 1.000
  df_degree_popsize = df_degree_popsize %>% mutate(row=paste(paste("N=",N,sep=""),paste("K=",K,sep="")))
  df_output = df_degree_popsize %>% select(row,K,N, graphrandom,graphsmall_world,graphlattice,graphmodular,graphmodular_lattice,graphmultilevel)
  df_output = df_output %>% pivot_longer(cols=c(graphrandom,graphsmall_world,graphlattice,graphmodular,graphmodular_lattice,graphmultilevel))
  df_output$name = factor(df_output$name,levels = c("graphrandom","graphsmall_world","graphlattice","graphmodular","graphmodular_lattice","graphmultilevel"))
  levels(df_output$name) <- c("Random","Small world","Lattice","Modular","Modular lattice","Multilevel") 
  N.labs <- c("Network size N=64","Network size N=144","Network size N=324")
  names(N.labs) <- c(64,144,324)
  m2_heatmap = ggplot(data=df_output,aes(x=name,y=as.factor(K)))+
    facet_wrap(~N,scales="free",labeller = labeller(N=N.labs))+
    geom_tile(aes(fill=value),size=.2,color="black")+
    geom_text(aes(label = sprintf("%0.3f", round(value, digits = 3)),color=value<.98),size=3,show.legend = F) +
    scale_color_manual(values=c("black","white"))+
    scale_fill_gradient2(high='red',low='darkblue',mid='white',midpoint=1)+
    labs(x="",y="Network connectivity (degree)",fill=expression(paste("exp(", beta, ")")))+
    theme_minimal()+
    theme(axis.text.x = element_blank(),strip.background = element_blank(),
          strip.text.x = element_blank())
  m2_heatmap
  ggsave("cox_ph_output/arch_performance_heatmap_m2.pdf",height=7,width=17.8,units="cm",scale=1.5)
  
}

#M2 raw TTD by degree and pop size
{
  load("df_TTC_TTD_m2.Rda")
  df = df_TTC_TTD_m2 %>% filter(graph!="full",pop_size==64,degree=="8") %>% droplevels()
  fit.coxph = glm(log(epoch_TTD+1)~graph, family=gaussian(link="log"), data=df)
  r1= exp(fit.coxph$coefficients)
  
  df = df_TTC_TTD_m2 %>% filter(graph!="full",pop_size==64,degree=="12") %>% droplevels()
  fit.coxph = glm(log(epoch_TTD+1)~graph, family=gaussian(link="log"), data=df)
  r2= exp(fit.coxph$coefficients)
  
  df = df_TTC_TTD_m2 %>% filter(graph!="full",pop_size==144,degree=="8") %>% droplevels()
  fit.coxph = glm(log(epoch_TTD+1)~graph, family=gaussian(link="log"), data=df)
  r3= exp(fit.coxph$coefficients)
  
  df = df_TTC_TTD_m2 %>% filter(graph!="full",pop_size==144,degree=="12") %>% droplevels()
  fit.coxph = glm(log(epoch_TTD+1)~graph, family=gaussian(link="log"), data=df)
  r4= exp(fit.coxph$coefficients)

  
  df = df_TTC_TTD_m2 %>% filter(graph!="full",pop_size==144,degree=="18") %>% droplevels()
  fit.coxph = glm(log(epoch_TTD+1)~graph, family=gaussian(link="log"), data=df)
  r5 = exp(fit.coxph$coefficients)
  
  df = df_TTC_TTD_m2 %>% filter(graph!="full",pop_size==144,degree=="24") %>% droplevels()
  fit.coxph = glm(log(epoch_TTD+1)~graph, family=gaussian(link="log"), data=df)
  r6 = exp(fit.coxph$coefficients)
  
  df = df_TTC_TTD_m2 %>% filter(graph!="full",pop_size==324,degree=="8") %>% droplevels()
  fit.coxph = glm(log(epoch_TTD+1)~graph, family=gaussian(link="log"), data=df)
  r7= exp(fit.coxph$coefficients)
  
  df = df_TTC_TTD_m2 %>% filter(graph!="full",pop_size==324,degree=="12") %>% droplevels()
  fit.coxph = glm(log(epoch_TTD+1)~graph, family=gaussian(link="log"), data=df)
  r8= exp(fit.coxph$coefficients)
  
  df = df_TTC_TTD_m2 %>% filter(graph!="full",pop_size==324,degree=="18") %>% droplevels()
  fit.coxph = glm(log(epoch_TTD+1)~graph, family=gaussian(link="log"), data=df)
  r9 = exp(fit.coxph$coefficients)
  
  df = df_TTC_TTD_m2 %>% filter(graph!="full",pop_size==324,degree=="24") %>% droplevels()
  fit.coxph = glm(log(epoch_TTD+1)~graph, family=gaussian(link="log"), data=df)
  r10 = exp(fit.coxph$coefficients)
  
  df = df_TTC_TTD_m2 %>% filter(graph!="full",pop_size==324,degree=="30") %>% droplevels()
  fit.coxph = glm(log(epoch_TTD+1)~graph, family=gaussian(link="log"), data=df)
  r11 = exp(fit.coxph$coefficients)
  
  df_degree_popsize = as.data.frame(rbind(r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11),row.names = F)
  pop_size = c(64,64,144,144,144,144,324,324,324,324,324)
  degree = c(8,12,8,12,18,24,8,12,18,24,30)
  df_degree_popsize$N = pop_size
  df_degree_popsize$K = degree
  df_degree_popsize$graphrandom = 1
  df_degree_popsize = df_degree_popsize %>% mutate(row=paste(paste("N=",N,sep=""),paste("K=",K,sep="")))
  df_output = df_degree_popsize %>% select(row,K,N, graphrandom,graphsmall_world,graphlattice,graphmodular,graphmodular_lattice,graphmultilevel)
  df_output = df_output %>% pivot_longer(cols=c(graphrandom,graphsmall_world,graphlattice,graphmodular,graphmodular_lattice,graphmultilevel))
  df_output$name = factor(df_output$name,levels = c("graphrandom","graphsmall_world","graphlattice","graphmodular","graphmodular_lattice","graphmultilevel"))
  levels(df_output$name) <- c("Random","Small world","Lattice","Modular","Modular lattice","Multilevel") 
  
  N.labs <- c("Network size N=64","Network size N=144","Network size N=324")
  names(N.labs) <- c(64,144,324)
  m2_heatmap_TTD = ggplot(data=df_output,aes(x=name,y=as.factor(K)))+
    facet_wrap(~N,scales="free",labeller = labeller(N=N.labs))+
    geom_tile(aes(fill=value),size=.2,color="black")+
    geom_text(aes(label = sprintf("%0.3f", round(value, digits = 3)),color=value>1.12),size=3,show.legend = F) +
    scale_color_manual(values=c("black","white"))+
    scale_fill_gradient2(high='red',low='darkblue',mid='white',midpoint=1)+
    labs(x="Network architecture",y="",fill=expression(paste("exp(", beta, ")")))+
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 20,hjust=1),strip.background = element_blank(),
          strip.text.x = element_blank())
  m2_heatmap_TTD
  ggsave("cox_ph_output/arch_performance_heatmap_m2_TTD.pdf",height=7,width=17.8,units="cm",scale=1.5)
}

ggarrange(m1_heatmap,m2_heatmap,m2_heatmap_TTD,labels=c("A","B","C"), ncol=1)
ggsave("cox_ph_output/arch_performance_heatmap_rawTTD_GLM.png",height=15,width=17.8,units="cm",scale=1.5)
ggsave("cox_ph_output/arch_performance_heatmap_rawTTD_GLM.pdf",height=15,width=17.8,units="cm",scale=1.5)

#M2 TTD-TTC by degree and pop size
{
  load("df_TTC_TTD_m2.Rda")
  summary(df_TTC_TTD_m2)
  df = df_TTC_TTD_m2 %>% filter(graph!="full",pop_size==64,degree=="8") %>% droplevels()
  df$pop_size = as.factor(df$pop_size)
  surv_object <- Surv(time = df$delta_TTC_TTD)
  summary(df)
  fit.coxph = coxph(surv_object ~ graph, data=df)
  summary(fit.coxph)
  r1= exp(fit.coxph$coefficients)
  p = ggforest(fit.coxph, data = df)
  p
  #ggsave("cox_ph_output/subset_models/M2_N64_K8_TTD_TTC.png")
  
  df = df_TTC_TTD_m2 %>% filter(graph!="full",pop_size==64,degree=="12") %>% droplevels()
  df$pop_size = as.factor(df$pop_size)
  surv_object <- Surv(time = df$delta_TTC_TTD)
  summary(df)
  fit.coxph = coxph(surv_object ~ graph, data=df)
  summary(fit.coxph)
  r2= exp(fit.coxph$coefficients)
  p = ggforest(fit.coxph, data = df)
  p
  #ggsave("cox_ph_output/subset_models/M2_N64_K12_TTD_TTC.png")
  
  df = df_TTC_TTD_m2 %>% filter(graph!="full",pop_size==144,degree=="8") %>% droplevels()
  df$pop_size = as.factor(df$pop_size)
  surv_object <- Surv(time = df$delta_TTC_TTD)
  summary(df)
  fit.coxph = coxph(surv_object ~ graph, data=df)
  summary(fit.coxph)
  r3= exp(fit.coxph$coefficients)
  p = ggforest(fit.coxph, data = df)
  p
  #ggsave("cox_ph_output/subset_models/M2_N144_K8_TTD_TTC.png")
  
  df = df_TTC_TTD_m2 %>% filter(graph!="full",pop_size==144,degree=="12") %>% droplevels()
  df$pop_size = as.factor(df$pop_size)
  surv_object <- Surv(time = df$delta_TTC_TTD)
  summary(df)
  fit.coxph = coxph(surv_object ~ graph, data=df)
  summary(fit.coxph)
  r4= exp(fit.coxph$coefficients)
  p = ggforest(fit.coxph, data = df)
  p
  #ggsave("cox_ph_output/subset_models/M2_N144_K12_TTD_TTC.png")
  
  df = df_TTC_TTD_m2 %>% filter(graph!="full",pop_size==144,degree=="18") %>% droplevels()
  df$pop_size = as.factor(df$pop_size)
  surv_object <- Surv(time = df$delta_TTC_TTD)
  summary(df)
  fit.coxph = coxph(surv_object ~ graph, data=df)
  summary(fit.coxph)
  r5 = exp(fit.coxph$coefficients)
  p = ggforest(fit.coxph, data = df)
  p
  #ggsave("cox_ph_output/subset_models/M2_N144_K18_TTD_TTC.png")
  
  df = df_TTC_TTD_m2 %>% filter(graph!="full",pop_size==144,degree=="24") %>% droplevels()
  df$pop_size = as.factor(df$pop_size)
  surv_object <- Surv(time = df$delta_TTC_TTD)
  summary(df)
  fit.coxph = coxph(surv_object ~ graph, data=df)
  summary(fit.coxph)
  r6 = exp(fit.coxph$coefficients)
  p = ggforest(fit.coxph, data = df)
  p
  #ggsave("cox_ph_output/subset_models/M2_N144_K24_TTD_TTC.png")
  
  df = df_TTC_TTD_m2 %>% filter(graph!="full",pop_size==324,degree=="8") %>% droplevels()
  df$pop_size = as.factor(df$pop_size)
  surv_object <- Surv(time = df$delta_TTC_TTD)
  summary(df)
  fit.coxph = coxph(surv_object ~ graph, data=df)
  summary(fit.coxph)
  r7= exp(fit.coxph$coefficients)
  p = ggforest(fit.coxph, data = df)
  p
  #ggsave("cox_ph_output/subset_models/M2_N324_K8_TTD_TTC.png")
  
  df = df_TTC_TTD_m2 %>% filter(graph!="full",pop_size==324,degree=="12") %>% droplevels()
  df$pop_size = as.factor(df$pop_size)
  surv_object <- Surv(time = df$delta_TTC_TTD)
  summary(df)
  fit.coxph = coxph(surv_object ~ graph, data=df)
  summary(fit.coxph)
  r8= exp(fit.coxph$coefficients)
  p = ggforest(fit.coxph, data = df)
  p
  #ggsave("cox_ph_output/subset_models/M2_N324_K12_TTD_TTC.png")
  
  df = df_TTC_TTD_m2 %>% filter(graph!="full",pop_size==324,degree=="18") %>% droplevels()
  df$pop_size = as.factor(df$pop_size)
  surv_object <- Surv(time = df$delta_TTC_TTD)
  summary(df)
  fit.coxph = coxph(surv_object ~ graph, data=df)
  summary(fit.coxph)
  r9 = exp(fit.coxph$coefficients)
  p = ggforest(fit.coxph, data = df)
  p
  #ggsave("cox_ph_output/subset_models/M2_N324_K18_TTD_TTC.png")
  
  df = df_TTC_TTD_m2 %>% filter(graph!="full",pop_size==324,degree=="24") %>% droplevels()
  df$pop_size = as.factor(df$pop_size)
  surv_object <- Surv(time = df$delta_TTC_TTD)
  summary(df)
  fit.coxph = coxph(surv_object ~ graph, data=df)
  summary(fit.coxph)
  r10 = exp(fit.coxph$coefficients)
  p = ggforest(fit.coxph, data = df)
  p
  #ggsave("cox_ph_output/subset_models/M2_N324_K24_TTD_TTC.png")
  
  df = df_TTC_TTD_m2 %>% filter(graph!="full",pop_size==324,degree=="30") %>% droplevels()
  df$pop_size = as.factor(df$pop_size)
  surv_object <- Surv(time = df$delta_TTC_TTD)
  summary(df)
  fit.coxph = coxph(surv_object ~ graph, data=df)
  summary(fit.coxph)
  r11 = exp(fit.coxph$coefficients)
  p = ggforest(fit.coxph, data = df)
  p
  #ggsave("cox_ph_output/subset_models/M2_N324_K30_TTD_TTC.png")
  paste("a","b")
  
  df_degree_popsize = as.data.frame(rbind(r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11),row.names = F)
  pop_size = c(64,64,144,144,144,144,324,324,324,324,324)
  degree = c(8,12,8,12,18,24,8,12,18,24,30)
  df_degree_popsize$N = pop_size
  df_degree_popsize$K = degree
  df_degree_popsize$graphrandom = 1
  df_degree_popsize = df_degree_popsize %>% mutate(row=paste(paste("N=",N,sep=""),paste("K=",K,sep="")))
  df_output = df_degree_popsize %>% select(row,K,N, graphrandom,graphsmall_world,graphlattice,graphmodular,graphmodular_lattice,graphmultilevel)
  df_output = df_output %>% pivot_longer(cols=c(graphrandom,graphsmall_world,graphlattice,graphmodular,graphmodular_lattice,graphmultilevel))
  df_output$name = factor(df_output$name,levels = c("graphrandom","graphsmall_world","graphlattice","graphmodular","graphmodular_lattice","graphmultilevel"))
  levels(df_output$name) <- c("Random","Small world","Lattice","Modular","Modular lattice","Multilevel") 
  
  N.labs <- c("Network size N=64","Network size N=144","Network size N=324")
  names(N.labs) <- c(64,144,324)
  m2_heatmap_TTD = ggplot(data=df_output,aes(x=name,y=as.factor(K)))+
    facet_wrap(~N,scales="free",labeller = labeller(N=N.labs))+
    geom_tile(aes(fill=value),size=.2,color="black")+
    geom_text(aes(label = round(value, 3),color=value>2),size=3,show.legend = F) +
    scale_color_manual(values=c("black","white"))+
    scale_fill_gradient2(high='darkblue',low='red',mid='white',midpoint=1)+
    labs(x="Network architecture",y="Network connectivity (degree)",fill="Hazard\nratio")+
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 20,hjust=1))
  m2_heatmap_TTD
  ggsave("cox_ph_output/arch_performance_heatmap_m2_TTD.pdf",height=7,width=17.8,units="cm",scale=1.5)
}
ggarrange(m1_heatmap,m2_heatmap,m2_heatmap_TTD,labels=c("A","B","C"), ncol=1)
ggsave("cox_ph_output/arch_performance_heatmap_TTD-TTC.png",height=15,width=17.8,units="cm",scale=1.5)
ggsave("cox_ph_output/arch_performance_heatmap_TTD-TTC.pdf",height=15,width=17.8,units="cm",scale=1.5)

#FOR SUPPLEMENTALS
#M1 by degree and architecture
{
  df = df_TTC_m1 %>% filter(graph=="random",pop_size==64) %>% droplevels()
  df$pop_size = as.factor(df$pop_size)
  surv_object <- Surv(time = df$epoch)
  fit.coxph = coxph(surv_object ~ degree, data=df)
  r1= exp(fit.coxph$coefficients)
  p = ggforest(fit.coxph, data = df)
  p
  ggsave("cox_ph_output/subset_models/M1_N64_random.png")
  
  df = df_TTC_m1 %>% filter(graph=="small_world",pop_size==64) %>% droplevels()
  df$pop_size = as.factor(df$pop_size)
  surv_object <- Surv(time = df$epoch)
  fit.coxph = coxph(surv_object ~ degree, data=df)
  r2= exp(fit.coxph$coefficients)
  p = ggforest(fit.coxph, data = df)
  p
  ggsave("cox_ph_output/subset_models/M1_N64_small_world.png")
  
  df = df_TTC_m1 %>% filter(graph=="lattice",pop_size==64) %>% droplevels()
  df$pop_size = as.factor(df$pop_size)
  surv_object <- Surv(time = df$epoch)
  fit.coxph = coxph(surv_object ~ degree, data=df)
  r3= exp(fit.coxph$coefficients)
  p = ggforest(fit.coxph, data = df)
  p
  ggsave("cox_ph_output/subset_models/M1_N64_lattice.png")
  
  df = df_TTC_m1 %>% filter(graph=="modular",pop_size==64) %>% droplevels()
  df$pop_size = as.factor(df$pop_size)
  surv_object <- Surv(time = df$epoch)
  fit.coxph = coxph(surv_object ~ degree, data=df)
  r4= exp(fit.coxph$coefficients)
  p = ggforest(fit.coxph, data = df)
  p
  ggsave("cox_ph_output/subset_models/M1_N64_modular.png")
  
  df = df_TTC_m1 %>% filter(graph=="modular_lattice",pop_size==64) %>% droplevels()
  df$pop_size = as.factor(df$pop_size)
  surv_object <- Surv(time = df$epoch)
  fit.coxph = coxph(surv_object ~ degree, data=df)
  r5= exp(fit.coxph$coefficients)
  p = ggforest(fit.coxph, data = df)
  p
  ggsave("cox_ph_output/subset_models/M1_N64_modular_lattice.png")
  
  df = df_TTC_m1 %>% filter(graph=="multilevel",pop_size==64) %>% droplevels()
  df$pop_size = as.factor(df$pop_size)
  surv_object <- Surv(time = df$epoch)
  fit.coxph = coxph(surv_object ~ degree, data=df)
  r6= exp(fit.coxph$coefficients)
  p = ggforest(fit.coxph, data = df)
  p
  ggsave("cox_ph_output/subset_models/M1_N64_multilevel.png")
  
  #N 144
  df = df_TTC_m1 %>% filter(graph=="random",pop_size==144) %>% droplevels()
  df$pop_size = as.factor(df$pop_size)
  surv_object <- Surv(time = df$epoch)
  fit.coxph = coxph(surv_object ~ degree, data=df)
  r7= exp(fit.coxph$coefficients)
  p = ggforest(fit.coxph, data = df)
  p
  ggsave("cox_ph_output/subset_models/M1_N144_random.png")
  
  df = df_TTC_m1 %>% filter(graph=="small_world",pop_size==144) %>% droplevels()
  df$pop_size = as.factor(df$pop_size)
  surv_object <- Surv(time = df$epoch)
  fit.coxph = coxph(surv_object ~ degree, data=df)
  r8= exp(fit.coxph$coefficients)
  p = ggforest(fit.coxph, data = df)
  p
  ggsave("cox_ph_output/subset_models/M1_N144_small_world.png")
  
  df = df_TTC_m1 %>% filter(graph=="lattice",pop_size==144) %>% droplevels()
  df$pop_size = as.factor(df$pop_size)
  surv_object <- Surv(time = df$epoch)
  fit.coxph = coxph(surv_object ~ degree, data=df)
  r9= exp(fit.coxph$coefficients)
  p = ggforest(fit.coxph, data = df)
  p
  ggsave("cox_ph_output/subset_models/M1_N144_lattice.png")
  
  df = df_TTC_m1 %>% filter(graph=="modular",pop_size==144) %>% droplevels()
  df$pop_size = as.factor(df$pop_size)
  surv_object <- Surv(time = df$epoch)
  fit.coxph = coxph(surv_object ~ degree, data=df)
  r10= exp(fit.coxph$coefficients)
  p = ggforest(fit.coxph, data = df)
  p
  ggsave("cox_ph_output/subset_models/M1_N144_modular.png")
  
  df = df_TTC_m1 %>% filter(graph=="modular_lattice",pop_size==144) %>% droplevels()
  df$pop_size = as.factor(df$pop_size)
  surv_object <- Surv(time = df$epoch)
  fit.coxph = coxph(surv_object ~ degree, data=df)
  r11= exp(fit.coxph$coefficients)
  p = ggforest(fit.coxph, data = df)
  p
  ggsave("cox_ph_output/subset_models/M1_N144_modular_lattice.png")
  
  df = df_TTC_m1 %>% filter(graph=="multilevel",pop_size==144) %>% droplevels()
  df$pop_size = as.factor(df$pop_size)
  surv_object <- Surv(time = df$epoch)
  fit.coxph = coxph(surv_object ~ degree, data=df)
  r12= exp(fit.coxph$coefficients)
  p = ggforest(fit.coxph, data = df)
  p
  ggsave("cox_ph_output/subset_models/M1_N144_multilevel.png")
  
  #N 324
  df = df_TTC_m1 %>% filter(graph=="random",pop_size==324) %>% droplevels()
  df$pop_size = as.factor(df$pop_size)
  surv_object <- Surv(time = df$epoch)
  fit.coxph = coxph(surv_object ~ degree, data=df)
  r13= exp(fit.coxph$coefficients)
  p = ggforest(fit.coxph, data = df)
  p
  ggsave("cox_ph_output/subset_models/M1_N324_random.png")
  
  df = df_TTC_m1 %>% filter(graph=="small_world",pop_size==324) %>% droplevels()
  df$pop_size = as.factor(df$pop_size)
  surv_object <- Surv(time = df$epoch)
  fit.coxph = coxph(surv_object ~ degree, data=df)
  r14= exp(fit.coxph$coefficients)
  p = ggforest(fit.coxph, data = df)
  p
  ggsave("cox_ph_output/subset_models/M1_N324_small_world.png")
  
  df = df_TTC_m1 %>% filter(graph=="lattice",pop_size==324) %>% droplevels()
  df$pop_size = as.factor(df$pop_size)
  surv_object <- Surv(time = df$epoch)
  fit.coxph = coxph(surv_object ~ degree, data=df)
  r15= exp(fit.coxph$coefficients)
  p = ggforest(fit.coxph, data = df)
  p
  ggsave("cox_ph_output/subset_models/M1_N324_lattice.png")
  
  df = df_TTC_m1 %>% filter(graph=="modular",pop_size==324) %>% droplevels()
  df$pop_size = as.factor(df$pop_size)
  surv_object <- Surv(time = df$epoch)
  fit.coxph = coxph(surv_object ~ degree, data=df)
  r16= exp(fit.coxph$coefficients)
  p = ggforest(fit.coxph, data = df)
  p
  ggsave("cox_ph_output/subset_models/M1_N324_modular.png")
  
  df = df_TTC_m1 %>% filter(graph=="modular_lattice",pop_size==324) %>% droplevels()
  df$pop_size = as.factor(df$pop_size)
  surv_object <- Surv(time = df$epoch)
  fit.coxph = coxph(surv_object ~ degree, data=df)
  r17= exp(fit.coxph$coefficients)
  p = ggforest(fit.coxph, data = df)
  p
  ggsave("cox_ph_output/subset_models/M1_N324_modular_lattice.png")
  
  df = df_TTC_m1 %>% filter(graph=="multilevel",pop_size==324) %>% droplevels()
  df$pop_size = as.factor(df$pop_size)
  surv_object <- Surv(time = df$epoch)
  fit.coxph = coxph(surv_object ~ degree, data=df)
  r18= exp(fit.coxph$coefficients)
  p = ggforest(fit.coxph, data = df)
  p
  ggsave("cox_ph_output/subset_models/M1_N324_multilevel.png")
  
  df_degree_arch = as.data.frame(rbind(r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,r14,r15,r16,r17,r18),row.names = F)
  
}


