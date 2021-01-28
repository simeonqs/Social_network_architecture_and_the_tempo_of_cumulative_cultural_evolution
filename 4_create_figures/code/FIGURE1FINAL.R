####### FIGURE 1 ###########



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




# 2. INPUT DATA -----------------------------------------------------------


# LOAD AND PICK DATA FOR MODEL1 OR MODEL 2 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
load("../../2_Python_agent_based_models/data/df_TTC_m1.Rda")

df = df_TTC_m1

# 3. CREATE NETWORKS ------------------------------------------------------


# CREATING NETWORKS N=324, K=12


set.seed(123)

N <- 324
degree <- 12

# Details
details <- data.frame(type=c("FULL","SMALLWORLD","DEGREE","CLUSTERED","MODULAR","MOD_CLUST","MULTILEVEL"), 
                      DENSITY=NA, DEGREE=NA, CLUSTERING=NA, MODULARITY=NA, MEAN_DISTANCE=NA)


# Function to make a lattice (with no boundaries)
lattice <- function(N1, N2, degree) {
	N <- N1*N2
	xy <- expand.grid(x=1:N1, y=1:N2)
	network <- matrix(0, N, N)
	for (i in 1:nrow(xy)) {
		xy.tmp <- xy
		xy.tmp$x <- xy.tmp$x[i] - xy.tmp$x
		xy.tmp$y <- xy.tmp$y[i]- xy.tmp$y
		xy.tmp$x[xy.tmp$x <= -(max(xy$x)-2)] <- max(xy$x) + xy.tmp$x[xy.tmp$x <= -(max(xy$x)-2)]
		xy.tmp$y[xy.tmp$y <= -(max(xy$x)-2)] <- max(xy$y) + xy.tmp$y[xy.tmp$y <= -(max(xy$x)-2)]
		xy.tmp$x[xy.tmp$x >= (max(xy$x)-3)] <- max(xy$x) - xy.tmp$x[xy.tmp$x >= (max(xy$x)-3)]
		xy.tmp$y[xy.tmp$y >= (max(xy$x)-3)] <- max(xy$y) - xy.tmp$y[xy.tmp$y >= (max(xy$x)-3)]
		dist <- sqrt((xy.tmp$x-xy.tmp$x[i])^2 + (xy.tmp$y-xy.tmp$y[i])^2)
		network[i,order(dist)[1:((degree)+1)]] <- 1
		network[order(dist)[1:((degree)+1)],i] <- 1
	}
	diag(network) <- 0
	return(network)
}



## FULLY CONNECTED ##

network <- matrix(1, nrow=N, ncol=N)
diag(network) <- 0
d <- graph.adjacency(network,mode="undirected",diag=FALSE)
zz <- 1
details$CLUSTERING[zz] <- transitivity(d)
details$DEGREE[zz] <- mean(igraph::degree(d))
details$MODULARITY[zz] <- modularity(fastgreedy.community(d))
details$DENSITY[zz] <- edge_density(d)
details$MEAN_DISTANCE[zz] <- mean_distance(d)

## new layout
net <- as.network(network, 
                 matrix.type='adjacency',
                 directed = F,
                 ignore.eval=FALSE,
                 names.eval='weight')
full = ggnet2(net,
       color = matviridis[1,][fastgreedy.community(d)$membership],
       # color = sort(viridis(10, begin = 0.1, end = 1, direction = 1, option = 'plasma')[fastgreedy.community(d)$membership]),
       # alpha = 0.9, 
       size = .7, 
       edge.size = 0.5,
       edge.color = 'grey50',
       edge.alpha = 0.2,
       mode = 'kamadakawai',
       label= NA,
       label.size=0,
       legend.position = 0,
       legend.size = 0
       ) + ggtitle('Fully connected')

rm(d, network)



## RANDOM ##

# Degree controlled only (low transitivity, low modularity)
d <- sample_k_regular(N,degree)
zz <- 3
details$CLUSTERING[zz] <- transitivity(d)
details$DEGREE[zz] <- mean(igraph::degree(d))
details$MODULARITY[zz] <- modularity(fastgreedy.community(d))
details$DENSITY[zz] <- edge_density(d)
details$MEAN_DISTANCE[zz] <- mean_distance(d)

#d <- graph.adjacency(network, mode="undirected")
#network <- graph.adjacency(d, mode="undirected")
network <- as.matrix(as_adjacency_matrix(d))

## new layout
net <- as.network(network, 
                 matrix.type='adjacency',
                 directed = F,
                 ignore.eval=FALSE,
                 names.eval='weight')
degr = ggnet2(net,
        color = matviridis[2,][fastgreedy.community(d)$membership],
        #color = sort(viridis(10, begin = 0.1, end = 1, direction = 1, option = 'plasma')[fastgreedy.community(d)$membership]), # with 'modules'
       #color = sort(viridis(10, begin = 0.1, end = 1, direction = 1, option = 'plasma')[1]),
       #alpha = 0.9, 
       size = .7, 
       edge.size = 0.5,
       edge.color = 'grey50',
       edge.alpha = 0.2,
       mode = 'kamadakawai',
       label= NA,
       label.size=0,
       legend.position = 0,
       legend.size = 0
       )+ ggtitle('Random')

rm(d, network)



## SMALL WORLD ##


## [Small world]
d <- sample_smallworld(1, N, degree/2, 0.05)

network <- as.matrix(as_adjacency_matrix(d))
zz <- 2
details$CLUSTERING[zz] <- transitivity(d)
details$DEGREE[zz] <- mean(igraph::degree(d))
details$MODULARITY[zz] <- modularity(fastgreedy.community(d))
details$DENSITY[zz] <- edge_density(d)
details$MEAN_DISTANCE[zz] <- mean_distance(d)

## new layout
net <- as.network(network, 
                 matrix.type='adjacency',
                 directed = F,
                 ignore.eval=FALSE,
                 names.eval='weight')
small = ggnet2(net,
       color = matviridis[3,][fastgreedy.community(d)$membership],
       # color = sort(viridis(10, begin = 0.1, end = 1, direction = 1, option = 'plasma')[fastgreedy.community(d)$membership]),
       # alpha = 0.9, 
       size = .7, 
       edge.size = 0.5,
       edge.color = 'grey50',
       edge.alpha = 0.2,
       mode = 'kamadakawai',
       label= NA,
       label.size=0,
       legend.position = 0,
       legend.size = 0
       )+ ggtitle('Small-world')

rm(d, network)



## LATTICE ##

network <- lattice(sqrt(N),sqrt(N), degree)

# fix having too many edges from the lattice
	if (mean(rowSums(network)) > degree) {
		edges <- which(network>0, arr.ind=TRUE)
		edges <- edges[which(edges[,1] < edges[,2]),]
		edges.to.remove <- nrow(edges)-degree*N/2
		remove <- round(seq(1,nrow(edges),length.out=edges.to.remove))
		edges <- edges[remove,]
		network[cbind(edges[,1],edges[,2])]<-0
		network[cbind(edges[,2],edges[,1])]<-0
	}

d <- graph.adjacency(network, mode="undirected")
zz <- 4
details$CLUSTERING[zz] <- transitivity(d)
details$DEGREE[zz] <- mean(igraph::degree(d))
details$MODULARITY[zz] <- modularity(fastgreedy.community(d))
details$DENSITY[zz] <- edge_density(d)
details$MEAN_DISTANCE[zz] <- mean_distance(d)

## new layout
net <- as.network(network, 
                 matrix.type='adjacency',
                 directed = F,
                 ignore.eval=FALSE,
                 names.eval='weight')
clus = ggnet2(net,
      color = matviridis[4,][fastgreedy.community(d)$membership],
       # color = sort(viridis(10, begin = 0.1, end = 1, direction = 1, option = 'plasma')[fastgreedy.community(d)$membership]), # with 'modules'
       # #color = sort(viridis(10, begin = 0.1, end = 1, direction = 1, option = 'plasma')[1]),
       # alpha = 0.9, 
       size = .7, 
       edge.size = 0.5,
       edge.color = 'grey50',
       edge.alpha = 0.2,
       mode = layout_on_grid(d),
       label= NA,
       label.size=0,
       legend.position = 0,
       legend.size = 0
       )+ ggtitle('Lattice')

rm(d, network)



## MODULAR ##

# modular but random within (10 modules, all connected to one-another)

mods <- 9
N.per.mod <- N/mods
network <- matrix(0, nrow=N, ncol=N)

for (i in 1:mods) {
	
	d <- sample_k_regular(N.per.mod,degree-1)
	indices <- (i-1)*N.per.mod + (1:N.per.mod)
	network[indices, indices] <- as.matrix(as_adjacency_matrix(d))

	# connect to neighbouring community
	mod.indices <- rep(1:mods, each= N.per.mod/mods)
	for (j in 1:mods) {
		if (i != j) {
			indices.mod.other <- (j-1)*N.per.mod + (1:N.per.mod)
			network[cbind(indices[which(mod.indices==j)], indices.mod.other[which(mod.indices==i)])] <- 1
			network[cbind(indices.mod.other[which(mod.indices==i)], indices[which(mod.indices==j)])] <- 1
		} else {
			network[indices[which(mod.indices == i)][1:2],indices[which(mod.indices == i)][1:2]] <- 1
			network[indices[which(mod.indices == i)][3:4],indices[which(mod.indices == i)][3:4]] <- 1
		}			
	}
}

diag(network) <- FALSE

d <- graph.adjacency(network, mode="undirected", diag=FALSE)
zz <- 5
details$CLUSTERING[zz] <- transitivity(d)
details$DEGREE[zz] <- mean(igraph::degree(d))
details$MODULARITY[zz] <- modularity(fastgreedy.community(d))
details$DENSITY[zz] <- edge_density(d)
details$MEAN_DISTANCE[zz] <- mean_distance(d)

## new layout
net <- as.network(network, 
                 matrix.type='adjacency',
                 directed = F,
                 ignore.eval=FALSE,
                 names.eval='weight')
modu = ggnet2(net,
       color = matviridis[5,][fastgreedy.community(d)$membership],
       # color = sort(viridis(10, begin = 0.1, end = 1, direction = 1, option = 'plasma')[fastgreedy.community(d)$membership]), 
       # alpha = 0.9, 
       size = .7, 
       edge.size = 0.5,
       edge.color = 'grey50',
       edge.alpha = 0.2,
       mode = 'fruchtermanreingold',
       label= NA,
       label.size=0,
       legend.position = 0,
       legend.size = 0
       )+ ggtitle('Modular')

rm(d, network)



## MODULAR LATTICE ##

# modular and clustered (10 modules, all connected to one-another)
mods <- 9
N.per.mod <- N/mods
network <- matrix(0, nrow=N, ncol=N)

for (i in 1:mods) {
	
	net <- lattice(sqrt(N.per.mod),sqrt(N.per.mod), degree-2)
	indices <- (i-1)*N.per.mod + (1:N.per.mod)
	network[indices, indices] <- net

	# connect to neighbouring community
	mod.indices <- rep(1:mods, each= N.per.mod/mods)
	for (j in 1:mods) {
		if (i != j) {
			indices.mod.other <- (j-1)*N.per.mod + (1:N.per.mod)
			network[cbind(indices[which(mod.indices==j)], indices.mod.other[which(mod.indices==i)])] <- 1
			network[cbind(indices.mod.other[which(mod.indices==i)], indices[which(mod.indices==j)])] <- 1
		} else {
			network[indices[which(mod.indices == i)][1:2],indices[which(mod.indices == i)][1:2]] <- 1
			network[indices[which(mod.indices == i)][3:4],indices[which(mod.indices == i)][3:4]] <- 1
		}			
	}
}

diag(network) <- FALSE

d <- graph.adjacency(network, mode="undirected", diag=FALSE)
zz <- 6
details$CLUSTERING[zz] <- transitivity(d)
details$DEGREE[zz] <- mean(igraph::degree(d))
details$MODULARITY[zz] <- modularity(fastgreedy.community(d))
details$DENSITY[zz] <- edge_density(d)
details$MEAN_DISTANCE[zz] <- mean_distance(d)

## new layout
net <- as.network(network, 
                 matrix.type='adjacency',
                 directed = F,
                 ignore.eval=FALSE,
                 names.eval='weight')
mocl = ggnet2(net,
               color = matviridis[6,][fastgreedy.community(d)$membership],
       # color = sort(viridis(10, begin = 0.1, end = 1, direction = 1, option = 'plasma')[fastgreedy.community(d)$membership]), 
       # alpha = 0.9, 
       size = .7, 
       edge.size = 0.5,
       edge.color = 'grey50',
       edge.alpha = 0.2,
       mode = 'fruchtermanreingold',
       label= NA,
       label.size=0,
       legend.position = 0,
       legend.size = 0
       )+ ggtitle('Modular lattice')




## MULTILEVEL ##

# Multilevel (3 modules of 3)
mods.inner <- 3
mods.outer <- 3

N2 <- N/mods.outer

N.per.mod <- N2/mods.inner
network.inner <- matrix(0, nrow=N2, ncol=N2)
network <- matrix(0, nrow=N, ncol=N)

for (z in 1:mods.outer) {
	indices.outer <- (z-1)*N.per.mod*mods.outer + (1:(N.per.mod*mods.outer))

	for (i in 1:mods.inner) {	
		net <- lattice(sqrt(N.per.mod),sqrt(N.per.mod), degree-2)
		indices <- (i-1)*N.per.mod + (1:N.per.mod)
		network.inner[indices, indices] <- net

		# connect to neighbouring community
		mod.indices <- rep(1:(mods.inner+1), each= N.per.mod/(mods.inner+1))
		for (j in 1:mods.inner) {
			if (i != j) {
				indices.mod.other <- (j-1)*N.per.mod + (1:N.per.mod)
				network.inner[cbind(indices[which(mod.indices==j)], indices.mod.other[which(mod.indices==i)])] <- 1
				network.inner[cbind(indices.mod.other[which(mod.indices==i)], indices[which(mod.indices==j)])] <- 1
			} else {
				network.inner[indices[which(mod.indices == i)][1:2],indices[which(mod.indices == i)][1:2]] <- 1
				network.inner[indices[which(mod.indices == i)][3:4],indices[which(mod.indices == i)][3:4]] <- 1
			}			
		}
	}

	# build full network
	network[indices.outer,indices.outer] <- network.inner
	
}

# now connect modules
mods.inner.indices <- 1:(mods.inner*mods.outer)
mod.outer.to.inner <- rep(1:mods.inner,each=mods.outer)
mod.indices.to.connect <- vector()
for (i in 1:mods.outer) {	
	mod.indices <- rep(1:(mods.inner+1), each= N.per.mod/(mods.inner+1))
	mod.indices[mod.indices <= mods.outer] <- 0
	mod.indices[mod.indices > 0] <- rep(mods.inner.indices, each=(sum(mod.indices > 0)/length(mods.inner.indices)))
	mod.indices <- rep(mod.indices, mods.inner)
	mod.indices.to.connect <- c(mod.indices.to.connect,mod.indices)
}

indices <- 1:N
for (i in 1:(mods.outer*mods.inner)) {
for (j in 1:(mods.outer*mods.inner)) {
	indices.mod.other.i <- (i-1)*N.per.mod + (1:N.per.mod)
	indices.mod.other.j <- (j-1)*N.per.mod + (1:N.per.mod)
	network[which(mod.indices.to.connect == j)[which(which(mod.indices.to.connect == j) %in% indices.mod.other.i)], which(mod.indices.to.connect == i)[which(which(mod.indices.to.connect == i) %in% indices.mod.other.j)]] <- 1
}
}

diag(network) <- FALSE

d <- graph.adjacency(network, mode="undirected", diag=FALSE)
zz <- 7
details$CLUSTERING[zz] <- transitivity(d)
details$DEGREE[zz] <- mean(igraph::degree(d))
details$MODULARITY[zz] <- modularity(fastgreedy.community(d))
details$DENSITY[zz] <- edge_density(d)
details$MEAN_DISTANCE[zz] <- mean_distance(d)

## new layout
net <- as.network(network, 
                 matrix.type='adjacency',
                 directed = F,
                 ignore.eval=FALSE,
                 names.eval='weight')
mls = ggnet2(net,
      color = matviridis[7,][fastgreedy.community(d)$membership],
       # color = sort(viridis(10, begin = 0.1, end = 1, direction = 1, option = 'plasma')[fastgreedy.community(d)$membership]),
       # alpha = 0.9, 
       size = .7, 
       edge.size = 0.5,
       edge.color = 'grey50',
       edge.alpha = 0.2,
       mode = 'fruchtermanreingold',
       label= NA,
       label.size=0,
       legend.position = 0,
       legend.size = 0
       )+ ggtitle('Multilevel')









# 4. SURVIVAL ANALYSES ----------------------------------------------------


### CALCULATING THE CUMULATIVE INCIDENCE ###

### facet by degree and type, colmuns by network size


## Network size N=324, varying degree ##

cu = df[which(df$pop_size==324 & df$graph!='full'),]
cu = transform(cu, graph=factor(graph,levels = (c("full", "degree", "small", "clustered", "modular", 'modularclustered', 'multilevel'))))

fit9 <- survfit(Surv(log(epoch+1)) ~ graph + degree , data = cu)

strata <- rep(1:length(fit9$strata),times=fit9$strata)
medians <- rep(NA, length(fit9$strata))
survmedian <- rep(NA, length(fit9$strata))
survmedian05 <- rep(0.5, length(fit9$strata))
for (i in 1:length(medians)) {
  fit.time <- fit9$time[which(strata == i)]
  fit.surv <- fit9$surv[which(strata == i)]
  medians[i] <- median(fit.time[which.min(abs(0.5-fit.surv))])
  survmedian[i] <- fit.surv[which.min(abs(0.5-fit.surv))]
}


## Visualize survival curves
surv_typedeg324 = ggsurvplot(
   fit9,                     # survfit object with calculated statistics.
  fun = "event", ### “event”: plots cumulative events (f(y) = 1-y). It’s also known as the cumulative incidence##
   #pval = TRUE,             # show p-value of log-rank test.
   conf.int = TRUE,         # show confidence intervals for point estimates of survival curves.
   #conf.int.style = "step", # customize style of confidence intervals
   xlab = "Time (epoch)",  
   ylab = "", 
    xlim = c(0, 4.67), # truncate by 1000 epochs
   ggtheme = theme_minimal() , # customize plot and risk table with a theme.
  surv.median.line = "v",  # add the median survival pointer.
  legend='none',
   palette = as.character(rev(c(#'grey',
              matviridis[7,c(10,8,6,4,2)],
              matviridis[6,c(10,8,6,4,2)],
              matviridis[5,c(10,8,6,4,2)],
              matviridis[4,c(10,8,6,4,2)],
              matviridis[3,c(10,8,6,4,2)],
              matviridis[2,c(10,8,6,4,2)])))
)

# adding median lines
dfmed <- data.frame(medians, survmedian, survmedian05, fit9$strata,
                    graph = rep(c("degree", "small", "clustered", "modular", 'modularclustered', 'multilevel'), each=5))

surv_typedeg324$plot <- surv_typedeg324$plot + 
  scale_x_continuous(limits = c(0, 8), expand = c(0, 0), labels = c('1', '10',  '100', '1000'), breaks = c(0, 2.31, 4.608, 6.908)) +
  facet_wrap(~graph, ncol=1) +
  geom_segment( data = dfmed[which(dfmed$graph=='degree'),], aes(y = 0, yend = survmedian05, x=medians, xend=medians), colour = "black", linetype='dashed')+
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.text.y = element_blank()) +
  ggtitle("Network size N=324")







## Network size N=144, varying degree ##

cu = df[which(df$pop_size==144 & df$graph!='full'),]
cu = transform(cu, graph=factor(graph,levels = (c("full", "degree", "small", "clustered", "modular", 'modularclustered', 'multilevel'))))

fit11 <- survfit(Surv(log(epoch+1)) ~ graph + degree + degree, data = cu)

strata <- rep(1:length(fit11$strata),times=fit11$strata)
medians <- rep(NA, length(fit11$strata))
survmedian <- rep(NA, length(fit11$strata))
survmedian05 <- rep(0.5, length(fit11$strata))
for (i in 1:length(medians)) {
  fit.time <- fit11$time[which(strata == i)]
  fit.surv <- fit11$surv[which(strata == i)]
  medians[i] <- median(fit.time[which.min(abs(0.5-fit.surv))])
  survmedian[i] <- fit.surv[which.min(abs(0.5-fit.surv))]
}


## Visualize survival curves
surv_typedeg144 = ggsurvplot(
   fit11,                     # survfit object with calculated statistics.
  fun = "event", ### “event”: plots cumulative events (f(y) = 1-y). It’s also known as the cumulative incidence##
   #pval = TRUE,             # show p-value of log-rank test.
   conf.int = TRUE,         # show confidence intervals for point estimates of survival curves.
   #conf.int.style = "step", # customize style of confidence intervals
   xlab = "Time (epoch)",  
   ylab = "", 
    xlim = c(0, 4.67), # truncate by 1000 epochs
   ggtheme = theme_minimal() , # customize plot and risk table with a theme.
  surv.median.line = "v",  # add the median survival pointer.
  legend='none',
   palette = as.character(rev(c(#'grey',
              matviridis[7,c(10,8,6,4,2)],
              matviridis[6,c(10,8,6,4,2)],
              matviridis[5,c(10,8,6,4,2)],
              matviridis[4,c(10,8,6,4,2)],
              matviridis[3,c(10,8,6,4,2)],
              matviridis[2,c(10,8,6,4,2)])))
)

# adding median lines
dfmed <- data.frame(medians, survmedian, survmedian05, fit11$strata,
                    graph = rep(c("degree", "small", "clustered", "modular", 'modularclustered', 'multilevel'), each=4))

surv_typedeg144$plot <- surv_typedeg144$plot + 
  scale_x_continuous(limits = c(0, 8), expand = c(0, 0), labels = c('1', '10',  '100', '1000'), breaks = c(0, 2.31, 4.608, 6.908)) +
  facet_wrap(~graph, ncol=1) +
   geom_segment( data = dfmed[which(dfmed$graph=='degree'),], aes(y = 0, yend = survmedian05, x=medians, xend=medians), colour = "black", linetype='dashed')+
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.text.y = element_blank()) +
  ggtitle("Network size N=144")








## Network size N=64, varying degree ### 

cu = df[which(df$pop_size==64 & df$graph!='full'),]
cu = transform(cu, graph=factor(graph,levels = (c("full", "degree", "small", "clustered", "modular", 'modularclustered', 'multilevel'))))

fit8 <- survfit(Surv(log(epoch+1)) ~ graph + degree, data = cu)

strata <- rep(1:length(fit8$strata),times=fit8$strata)
medians <- rep(NA, length(fit8$strata))
survmedian <- rep(NA, length(fit8$strata))
survmedian05 <- rep(0.5, length(fit8$strata))
for (i in 1:length(medians)) {
  fit.time <- fit8$time[which(strata == i)]
  fit.surv <- fit8$surv[which(strata == i)]
  medians[i] <- median(fit.time[which.min(abs(0.5-fit.surv))])
  survmedian[i] <- fit.surv[which.min(abs(0.5-fit.surv))]
}


## Visualize survival curves
surv_typedeg64 = ggsurvplot(
   fit8,                     # survfit object with calculated statistics.
  fun = "event", ### “event”: plots cumulative events (f(y) = 1-y). It’s also known as the cumulative incidence##
  #pval = TRUE,             # show p-value of log-rank test.
   conf.int = TRUE,         # show confidence intervals for point estimates of survival curves.
   #conf.int.style = "step", # customize style of confidence intervals
   xlab = "Time (epoch)",  
   ylab = "Probability of reaching recombination", 
    xlim = c(0, 4.67), # truncate by 1000 epochs
   ggtheme = theme_minimal() , # customize plot and risk table with a theme.
  surv.median.line = "v",  # add the median survival pointer.
  legend='none',
   palette = as.character(rev(c(#'grey',
              matviridis[7,c(4,2)],
              matviridis[6,c(4,2)],
              matviridis[5,c(4,2)],
              matviridis[4,c(4,2)],
              matviridis[3,c(4,2)],
              matviridis[2,c(4,2)])))
)

# adding median lines
dfmed <- data.frame(medians, survmedian, survmedian05, fit8$strata,
                    graph = rep(c("degree", "small", "clustered", "modular", 'modularclustered', 'multilevel'), each=2))

surv_typedeg64$plot <- surv_typedeg64$plot + 
  scale_x_continuous(limits = c(0, 8), expand = c(0, 0), labels = c('1', '10',  '100', '1000'), breaks = c(0, 2.31, 4.608, 6.908)) +
  geom_segment( data = dfmed[which(dfmed$graph=='degree'),], aes(y = 0, yend = survmedian05, x=medians, xend=medians), colour = "black", linetype='dashed')+
  facet_wrap(~graph, ncol=1) +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank()) +
  ggtitle("Network size N=64") 





# 5. PLOTTING -------------------------------------------------------------


degr <- arrangeGrob(degr, top = textGrob(expression(bold("A")), x = unit(0, "npc")
        , y   = unit(1, "npc"), just=c("left","top"),
        gp=gpar(col="black", fontsize=14, fontfamily="Arial", face='bold')))

surv_typedeg64$plot <- arrangeGrob(surv_typedeg64$plot,
                                   top = textGrob(expression(bold("B")), x = unit(0, "npc"),
                                                  y   = unit(1, "npc"), just=c("left","top"),
                                                  gp=gpar(col="black", fontsize=14, fontfamily="Arial", face='bold')))
surv_typedeg144$plot <- arrangeGrob(surv_typedeg144$plot, top = textGrob(" ", x = unit(0, "npc")
        , y   = unit(1, "npc"), just=c("left","top"),
        gp=gpar(col="black", fontsize=14, fontfamily="Arial", face='bold')))
surv_typedeg324$plot <- arrangeGrob(surv_typedeg324$plot, top = textGrob(" ", x = unit(0, "npc")
        , y   = unit(1, "npc"), just=c("left","top"),
        gp=gpar(col="black", fontsize=14, fontfamily="Arial", face='bold')))




        
## adding nets
lay <- rbind(c(1,7,7,8,8,9,9),
             c(2,7,7,8,8,9,9),
             c(3,7,7,8,8,9,9),
             c(4,7,7,8,8,9,9),
             c(5,7,7,8,8,9,9),
             c(6,7,7,8,8,9,9))
grid.arrange(degr, small, clus, modu, mocl, mls,
              surv_typedeg64$plot, surv_typedeg144$plot, surv_typedeg324$plot,
             layout_matrix = lay)
