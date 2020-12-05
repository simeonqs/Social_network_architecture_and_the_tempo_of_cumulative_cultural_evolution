if(!require(igraph)){install.packages('igraph'); library(igraph)}
if(!require(GGally)){install.packages('GGally'); library(GGally)}
if(!require(sna)){install.packages('sna'); library(sna)}
if(!require(viridis)){install.packages('viridis'); library(viridis)}
if(!require(gridExtra)){install.packages('gridExtra'); library(gridExtra)}
set.seed(123)

# Change this path
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
output_folder <- "/Users/ssmeele/Desktop/Social_network_architecture_and_the_tempo_of_cumulative_cultural_evolution-master"
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

for(N in c(64, 144, 324)){
  
  degrees = c(8, 12, 18, 24, 30)
  if(N == 64) degrees = c(8, 12)
  if(N == 144) degrees = c(8, 12, 18, 24)

  mods.inner <- 3
  mods.outer <- 3
  degree.fix <- 2
  mods <- 9
  if(N < 324) { mods.inner <- 2; mods.outer <- 2; degree.fix <- 1; mods <- 4}


  
  for(degree in degrees){
    
    # Details
    details <- data.frame(type=c("FULL","SMALLWORLD","DEGREE","CLUSTERED","MODULAR","MOD_CLUST","MULTILEVEL"), DENSITY=NA, DEGREE=NA, CLUSTERING=NA, MODULARITY=NA, MEAN_DISTANCE=NA)
    
    # Function to make a lattice (with no boundaries)
    lattice <- function(N1, N2, degree, fix.degree=FALSE, plot=FALSE) {
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
	
	if (mean(rowSums(network)) > degree & fix.degree) {
		edges <- which(network>0, arr.ind=TRUE)
		edges <- edges[which(edges[,1] < edges[,2]),]
		edges.to.remove <- nrow(edges)-degree*N/2
		remove <- round(seq(1,nrow(edges),length.out=edges.to.remove))
		edges <- edges[remove,]
		network[cbind(edges[,1],edges[,2])]<-0
		network[cbind(edges[,2],edges[,1])]<-0
	}

	if (plot) plot(graph.adjacency(network,mode="undirected"), layout=as.matrix(xy.tmp),vertex.size=1)
	
	return(network)
    }
    
    
    
    ## Full (high degree)
    
    network <- matrix(1, nrow=N, ncol=N)
    diag(network) <- 0
    d <- graph.adjacency(network,mode="undirected",diag=FALSE)
    zz <- 1
    details$CLUSTERING[zz] <- transitivity(d)
    details$DEGREE[zz] <- mean(igraph::degree(d))
    details$MODULARITY[zz] <- modularity(fastgreedy.community(d))
    details$DENSITY[zz] <- edge_density(d)
    details$MEAN_DISTANCE[zz] <- mean_distance(d)
    
    # Don't plot
    
    write_graph(d, file=paste(output_folder, "/3_R_agent_based_models/networks/Networks_N", N, "/full_d",degree,".pajek",sep=""), format="pajek")
    write_graph(d, file=paste(output_folder, "/2_Python_agent_based_models/networks/full/", N, "_full.txt",sep=""), format="pajek")
    d = read.table(paste(output_folder, "/2_Python_agent_based_models/networks/full/", N, "_full.txt",sep=""), 
                   skip = 2, header = F, sep = ' ')
    write.table(d, paste(output_folder, "/2_Python_agent_based_models/networks/full/", N, "_full.txt",sep=""), row.names = F, col.names = F)
    
    ## Small world
    
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
                   color = sort(viridis(10, begin = 0.1, end = 1, direction = 1, option = 'plasma')[fastgreedy.community(d)$membership]),
                   alpha = 0.9, 
                   size = 2.5, 
                   edge.size = 0.5,
                   edge.color = 'grey32',
                   edge.alpha = 0.4,
                   mode = 'kamadakawai',
                   label= NA,
                   label.size=0,
                   legend.position = 0,
                   legend.size = 0
    )
    #small <- arrangeGrob(small, top = textGrob("A", x = unit(0, "npc")
    #         , y   = unit(1, "npc"), just=c("left","top"),
    #         gp=gpar(col="black", fontsize=18, fontfamily="Arial")))
    
    
    write_graph(d, file=paste(output_folder, "/3_R_agent_based_models/networks/Networks_N", N, "/small_world_d",degree,".pajek",sep=""), format="pajek")
    write_graph(d, file=paste(output_folder, "/2_Python_agent_based_models/networks/smallworld/", N, "_small_world_d", degree, ".txt",sep=""), format="pajek")
    d = read.table(paste(output_folder, "/2_Python_agent_based_models/networks/smallworld/", N, "_small_world_d", degree, ".txt",sep=""), 
                   skip = 2, header = F, sep = ' ')
    write.table(d, paste(output_folder, "/2_Python_agent_based_models/networks/smallworld/", N, "_small_world_d", degree, ".txt",sep=""), row.names = F, col.names = F)
    
    ## Degree controlled only (low transitivity, low modularity)
    
    d <- sample_k_regular(N,degree)
    zz <- 3
    details$CLUSTERING[zz] <- transitivity(d)
    details$DEGREE[zz] <- mean(igraph::degree(d))
    details$MODULARITY[zz] <- modularity(fastgreedy.community(d))
    details$DENSITY[zz] <- edge_density(d)
    details$MEAN_DISTANCE[zz] <- mean_distance(d)
    
    d <- graph.adjacency(network, mode="undirected")
    
    
    ## new layout
    net <- as.network(network, 
                      matrix.type='adjacency',
                      directed = F,
                      ignore.eval=FALSE,
                      names.eval='weight')
    degr = ggnet2(net,
                  color = sort(viridis(10, begin = 0.1, end = 1, direction = 1, option = 'plasma')[fastgreedy.community(d)$membership]), # with 'modules'
                  #color = sort(viridis(10, begin = 0.1, end = 1, direction = 1, option = 'plasma')[1]),
                  alpha = 0.9, 
                  size = 2.5, 
                  edge.size = 0.5,
                  edge.color = 'grey32',
                  edge.alpha = 0.4,
                  mode = 'kamadakawai',
                  label= NA,
                  label.size=0,
                  legend.position = 0,
                  legend.size = 0
    )
    #degr <- arrangeGrob(degr, top = textGrob("B", x = unit(0, "npc")
    #         , y   = unit(1, "npc"), just=c("left","top"),
    #         gp=gpar(col="black", fontsize=18, fontfamily="Arial")))
    
    write_graph(d, file=paste(output_folder, "/3_R_agent_based_models/networks/Networks_N", N, "/degree_d",degree,".pajek",sep=""), format="pajek")
    write_graph(d, file=paste(output_folder, "/2_Python_agent_based_models/networks/degree/", N, "_degree_d", degree, ".txt",sep=""), format="pajek")
    d = read.table(paste(output_folder, "/2_Python_agent_based_models/networks/degree/", N, "_degree_d", degree, ".txt",sep=""), 
                   skip = 2, header = F, sep = ' ')
    write.table(d, paste(output_folder, "/2_Python_agent_based_models/networks/degree/", N, "_degree_d", degree, ".txt",sep=""), row.names = F, col.names = F)
    
    ## clustered
    
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
                  color = sort(viridis(10, begin = 0.1, end = 1, direction = 1, option = 'plasma')[fastgreedy.community(d)$membership]), # with 'modules'
                  #color = sort(viridis(10, begin = 0.1, end = 1, direction = 1, option = 'plasma')[1]),
                  alpha = 0.9, 
                  size = 2.5, 
                  edge.size = 0.5,
                  edge.color = 'grey32',
                  edge.alpha = 0.4,
                  mode = layout_on_grid(d),
                  label= NA,
                  label.size=0,
                  legend.position = 0,
                  legend.size = 0
    )
    #clus <- arrangeGrob(clus, top = textGrob("C", x = unit(0, "npc")
    #         , y   = unit(1, "npc"), just=c("left","top"),
    #         gp=gpar(col="black", fontsize=18, fontfamily="Arial")))
    
    
    write_graph(d, file=paste(output_folder, "/3_R_agent_based_models/networks/Networks_N", N, "/clustered_d",degree,".pajek",sep=""), format="pajek")
    write_graph(d, file=paste(output_folder, "/2_Python_agent_based_models/networks/clustered/", N, "_clustered_d", degree, ".txt",sep=""), format="pajek")
    d = read.table(paste(output_folder, "/2_Python_agent_based_models/networks/clustered/", N, "_clustered_d", degree, ".txt",sep=""), 
                   skip = 2, header = F, sep = ' ')
    write.table(d, paste(output_folder, "/2_Python_agent_based_models/networks/clustered/", N, "_clustered_d", degree, ".txt",sep=""), row.names = F, col.names = F)
    
    
    ## modular but random within (10 modules, all connected to one-another)
    
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
                  color = sort(viridis(10, begin = 0.1, end = 1, direction = 1, option = 'plasma')[fastgreedy.community(d)$membership]), 
                  alpha = 0.9, 
                  size = 2.5, 
                  edge.size = 0.5,
                  edge.color = 'grey32',
                  edge.alpha = 0.4,
                  mode = 'fruchtermanreingold',
                  label= NA,
                  label.size=0,
                  legend.position = 0,
                  legend.size = 0
    )
    #modu <- arrangeGrob(modu, top = textGrob("D", x = unit(0, "npc")
    #         , y   = unit(1, "npc"), just=c("left","top"),
    #         gp=gpar(col="black", fontsize=18, fontfamily="Arial")))
    
    
    write_graph(d, file=paste(output_folder, "/3_R_agent_based_models/networks/Networks_N", N, "/modular_d",degree,".pajek",sep=""), format="pajek")
    write_graph(d, file=paste(output_folder, "/2_Python_agent_based_models/networks/modular/", N, "_modular_d", degree, ".txt",sep=""), format="pajek")
    d = read.table(paste(output_folder, "/2_Python_agent_based_models/networks/modular/", N, "_modular_d", degree, ".txt",sep=""), 
                   skip = 2, header = F, sep = ' ')
    write.table(d, paste(output_folder, "/2_Python_agent_based_models/networks/modular/", N, "_modular_d", degree, ".txt",sep=""), row.names = F, col.names = F)
    
    
    ## modular and clustered (10 modules, all connected to one-another)
    
    N.per.mod <- N/mods
    network <- matrix(0, nrow=N, ncol=N)
    
    for (i in 1:mods) {
      
      net <- lattice(sqrt(N.per.mod),sqrt(N.per.mod), degree-degree.fix, N<144)
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
                  color = sort(viridis(10, begin = 0.1, end = 1, direction = 1, option = 'plasma')[fastgreedy.community(d)$membership]), 
                  alpha = 0.9, 
                  size = 2.5, 
                  edge.size = 0.5,
                  edge.color = 'grey32',
                  edge.alpha = 0.4,
                  mode = 'fruchtermanreingold',
                  label= NA,
                  label.size=0,
                  legend.position = 0,
                  legend.size = 0
    )
    #mocl <- arrangeGrob(mocl, top = textGrob("E", x = unit(0, "npc")
    #         , y   = unit(1, "npc"), just=c("left","top"),
    #         gp=gpar(col="black", fontsize=18, fontfamily="Arial")))
    
    
    write_graph(d, file=paste(output_folder, "/3_R_agent_based_models/networks/Networks_N", N, "/modular_clustered_d",degree,".pajek",sep=""), format="pajek")
    write_graph(d, file=paste(output_folder, "/2_Python_agent_based_models/networks/modularclustered/", N, "_modular_clustered_d", degree, ".txt",sep=""), format="pajek")
    d = read.table(paste(output_folder, "/2_Python_agent_based_models/networks/modularclustered/", N, "_modular_clustered_d", degree, ".txt",sep=""), 
                   skip = 2, header = F, sep = ' ')
    write.table(d, paste(output_folder, "/2_Python_agent_based_models/networks/modularclustered/", N, "_modular_clustered_d", degree, ".txt",sep=""), row.names = F, col.names = F)
    
    
    # Multilevel
    
    N2 <- N/mods.outer
    
    N.per.mod <- N2/mods.inner
    network.inner <- matrix(0, nrow=N2, ncol=N2)
    network <- matrix(0, nrow=N, ncol=N)
    
    for (z in 1:mods.outer) {
      indices.outer <- (z-1)*N.per.mod*mods.inner + (1:(N.per.mod*mods.inner))
      
      for (i in 1:mods.inner) {	
        net <- lattice(sqrt(N.per.mod),sqrt(N.per.mod), degree-degree.fix, N<144)
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
                 color = sort(viridis(10, begin = 0.1, end = 1, direction = 1, option = 'plasma')[fastgreedy.community(d)$membership]),
                 alpha = 0.9, 
                 size = 2.5, 
                 edge.size = 0.5,
                 edge.color = 'grey32',
                 edge.alpha = 0.4,
                 mode = 'fruchtermanreingold',
                 label= NA,
                 label.size=0,
                 legend.position = 0,
                 legend.size = 0
    )
    #mls <- arrangeGrob(mls, top = textGrob("F", x = unit(0, "npc")
    #         , y   = unit(1, "npc"), just=c("left","top"),
    #         gp=gpar(col="black", fontsize=18, fontfamily="Arial")))
    
    
    write_graph(d, file=paste(output_folder, "/3_R_agent_based_models/networks/Networks_N", N, "/multilevel_d",degree,".pajek",sep=""), format="pajek")
    write_graph(d, file=paste(output_folder, "/2_Python_agent_based_models/networks/multilevel/", N, "_multilevel_d", degree, ".txt",sep=""), format="pajek")
    d = read.table(paste(output_folder, "/2_Python_agent_based_models/networks/multilevel/", N, "_multilevel_d", degree, ".txt",sep=""), 
                   skip = 2, header = F, sep = ' ')
    write.table(d, paste(output_folder, "/2_Python_agent_based_models/networks/multilevel/", N, "_multilevel_d", degree, ".txt",sep=""), row.names = F, col.names = F)
    
    # Plotting all networks
    pdf(paste(output_folder, "/1_create_networks/results/Network_degree_",degree,"_N_",N,".pdf",sep=""), height=8, width=14)
    grid.arrange(small, degr, clus, modu, mocl, mls, ncol=3)
    dev.off()
    
    # Write out details
    write.csv(details, file=paste(output_folder, "/1_create_networks/results/DETAILS_d",degree,"_N_",N,".csv",sep=""))
  
  } # End degree loop
  
} # End N loop
