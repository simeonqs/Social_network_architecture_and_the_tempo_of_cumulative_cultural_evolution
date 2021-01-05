library(tidyverse)
library(ggpubr)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
load("../../2_Python_agent_based_models/data/df_TTC_m1.Rda")

#IQR
df_TTC_m1 %>% filter(graph!="full") %>% group_by(pop_size) %>% summarize(sd(epoch)/mean(epoch))
df_TTC_m1 %>% filter(graph!="full") %>% group_by(pop_size) %>% summarize(IQR(epoch)/(quantile(epoch,.25)+quantile(epoch,.75)))
summary(df_TTC_m1$epoch)

#Model 1 by degree and pop size
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
  ggsave("results/cox_ph_output/arch_performance_heatmap_GLM.pdf",height=7,width=17.8,units="cm",scale=1.5)
}


#Model 2 by degree and pop size
{
  load("../data/df_TTC_m2.Rda")
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

#Model 2 raw TTD by degree and pop size
{
  load("../data/df_TTC_TTD_m2.Rda")
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
ggsave("../output/arch_performance_heatmap_rawTTD_GLM.png",height=15,width=17.8,units="cm",scale=1.5)
ggsave("../output/arch_performance_heatmap_rawTTD_GLM.pdf",height=15,width=17.8,units="cm",scale=1.5)


