hub_detection <- function(num_iterations, adjacency, moduleColors, colors, betweenness = T, hubscore = T, pagerank = T, closeness = T, MM = T, edge = T, limit_edge = T, sig_value = 0.05, clusters = 5, datExpr = NULL) {
	
	require(WGCNA)
	require(igraph)
	require(data.table)
	require(WGCNA)
	require(foreach)
	require(doParallel)
	require(R.utils)
	
	#initialise blank lists to store results
	last = list()
	edge_values = list()

	#if the module colors being investigated aren't given, takes them from the moduleColors variable
	if (missing(colors)) {
		colors <- unique(moduleColors)
	}

	############
	#cycles through each module given to find hubs
	############
	for (i in 1:length(colors)){
		
		#print module name and number of genes in module
		print(paste0("Starting the ", colors[i], " module, module ", i, " of ", length(colors)))
		print(paste0("Genes in module: ",  length(moduleColors[moduleColors == colors[i]])))
		
		#initalise the variable (last) list with names of genes in module
		genes <- names(moduleColors[moduleColors == colors[i]])
		temp <- data.table(genes) #convert to data.table object - this has one coloumn which has gene names in it
		last[[i]] <- temp
		names(last)[i] <- colors[i] #changes name of the member in list to the module name being investigated  

		#create a graph object from the adjacency matrix using igraph
		network <- adjacency[genes, genes]
		network <- network[genes, genes] ##makes sure the order of cols and rows are correct
		network_graph <- graph_from_adjacency_matrix(network, mode = "undirected", weighted = TRUE, diag = FALSE,  add.colnames = NULL, add.rownames = NA)


		################################
		#calculate the hub scores in the actual data
		################################
		##if  betweenness and stuff ot her ones
		if(betweenness){
			network_betweenness <- betweenness(network_graph, v = V(network_graph), directed = FALSE, nobigint = TRUE, normalized = FALSE, weights = E(network_graph)$weight)
		}
		
		if(hubscore){
			network_hubscore <- hub_score(network_graph, scale = FALSE,weights = E(network_graph)$weight)$vector
		}

		if(pagerank){
			network_pagerank <- page_rank(network_graph, algo = "prpack", vids = V(network_graph), directed = FALSE, damping = 0.85, personalized = NULL, weights = E(network_graph)$weight,  options = NULL)$vector
		}

		if(closeness){
			network_closeness <- closeness(network_graph, vids = V(network_graph),  weights = E(network_graph)$weight, normalized = FALSE)
		}

		if(MM){
			geneModuleMembership1 = signedKME(datExpr, MEs)
			colnames(geneModuleMembership1) <- colnames(MEs)
			network_MM <- geneModuleMembership1[genes,colors[i], drop = FALSE]
			network_MM <- unlist(t(network_MM)[1,])
		}
		if(edge){
			network_edgebetweenness <- edge_betweenness(network_graph, e = E(network_graph), directed = FALSE,  weights = E(network_graph)$weight)
			names(network_edgebetweenness) <- as_ids(E(network_graph))
		}

		

		#########################
		#do the iteration test
		###########################

		set.seed(10)
		genes <- colnames(network)
		
		#set up clusters for parallel computing
		cl<-makeCluster(clusters)
		registerDoParallel(cl)

		if(betweenness){
			print(paste0("Finding betweenness values"))
			results_between <- foreach(icount(num_iterations), .packages = c("igraph")) %dopar% {
			    	#randomise the list of genes that are colnames of genes
				rand_genes <- sample(genes)
				#add the randomised genes as row and colnames to randomise the adjacency matrix
				colnames(network) <- rand_genes
				rownames(network) <- rand_genes
				#create the igraph object from new randomised adjacency matrix
				network_graph <- graph_from_adjacency_matrix(network, mode = "undirected", weighted = TRUE, diag = FALSE,  add.colnames = NULL, add.rownames = NA)
				#calculate the betweenness values for this new random labelled graph
				temp_network_betweenness <- betweenness(network_graph, v = V(network_graph), directed = FALSE, nobigint = TRUE, normalized = FALSE, weights = E(network_graph)$weight)
			}
		}

		if(hubscore){
			print(paste0("Finding hubscore values"))
			results_hubscore <- foreach(icount(num_iterations), .packages = c("igraph")) %dopar% {
			    	#randomise the list of genes that are colnames of genes
				rand_genes <- sample(genes)
				#create a dupicalte of the original adjacency matrix
				#add the randomised genes as row and colnames to randomise the adjacency matrix
				colnames(network) <- rand_genes
				rownames(network) <- rand_genes
				network_graph <- graph_from_adjacency_matrix(network, mode = "undirected", weighted = TRUE, diag = FALSE,  add.colnames = NULL, add.rownames = NA)
				#create the igraph object
				temp_network_hubscore <- hub_score(network_graph, scale = FALSE,weights = E(network_graph)$weight)$vector
				
			}
		}

		if(pagerank){
			print(paste0("Finding pagerank values"))
			results_pagerank <- foreach(icount(num_iterations), .packages = c("igraph")) %dopar% {
			    	#randomise the list of genes that are colnames of genes
				rand_genes <- sample(genes)
				#create a dupicalte of the original adjacency matrix
				#add the randomised genes as row and colnames to randomise the adjacency matrix
				colnames(network) <- rand_genes
				rownames(network) <- rand_genes
				network_graph <- graph_from_adjacency_matrix(network, mode = "undirected", weighted = TRUE, diag = FALSE,  add.colnames = NULL, add.rownames = NA)
				#create the igraph object
				temp_network_pagerank <- page_rank(network_graph, algo = "prpack", vids = V(network_graph), directed = FALSE, damping = 0.85, personalized = NULL, weights = E(network_graph)$weight,  options = NULL)$vector
			}
		}
		
		if(closeness){
			print(paste0("Finding closeness values"))
			results_closeness <- foreach(icount(num_iterations), .packages = c("igraph")) %dopar% {
			  	  #randomise the list of genes that are colnames of genes
				rand_genes <- sample(genes)
				#create a dupicalte of the original adjacency matrix
				#add the randomised genes as row and colnames to randomise the adjacency matrix
				colnames(network) <- rand_genes
				rownames(network) <- rand_genes
				network_graph <- graph_from_adjacency_matrix(network, mode = "undirected", weighted = TRUE, diag = FALSE,  add.colnames = NULL, add.rownames = NA)
				#create the igraph object
				temp_network_closeness <- closeness(network_graph, vids = V(network_graph),  weights = E(network_graph)$weight, normalized = FALSE)
			}
		}
		
		if(MM){
			print(paste0("Finding module membership values"))
			results_MM <- foreach(icount(num_iterations)) %dopar% {
			   	#randomise the list of genes that are colnames of genes
				rand_genes <- sample(genes)
				temp_network_MM <- network_MM
				names(temp_network_MM) <- rand_genes
				temp_network_MM
			}
		}

		
		if(edge){
			results_edgebetweenness <- foreach(icount(num_iterations), .packages = c("igraph")) %dopar% {
			    	#randomise the list of genes that are colnames of genes
				rand_genes <- sample(genes)
				#create a duplicate of the original adjacency matrix
				#add the randomised genes as row and colnames to randomise the adjacency matrix
				colnames(network) <- rand_genes
				rownames(network) <- rand_genes
				###as edgebetweenness uses the first gene it pulls from the network to report results to remove any repetitions (ie. remove SCNA|DCLK1 and DCLK1|SNCA, only using which gene comes up first in the column names)
				network <- network[genes, genes]
				network_graph <- graph_from_adjacency_matrix(network, mode = "undirected", weighted = TRUE, diag = FALSE,  add.colnames = NULL, add.rownames = NA)
				#create the igraph object
				temp_network_edgebetweenness <- edge_betweenness(network_graph, e = E(network_graph), directed = TRUE,  weights = E(network_graph)$weight)
				names(temp_network_edgebetweenness) <- as_ids(E(network_graph))
				temp_network_edgebetweenness
			}
		}

		stopCluster(cl)
		#################################
		#finding pvalues                #
		#################################
		#create a blank numeric vector for results
		pvalue=vector('numeric')
		print(paste0("Calculating pvalues"))
		##########################
		#if betweenness
		##########################
		#for each gene do following:
		if(betweenness){
			pb <- txtProgressBar(min = 0, max = length(genes), width = 60, style = 3)
			for (p in 1:length(genes)){
				setTxtProgressBar(pb, p)
				count = 0
				#check if the hub score is higher than the original score and if is then add a count
				for (l in 1:num_iterations){
					if (results_between[[l]][genes[p]] >= network_betweenness[genes[p]]){
						count = count + 1
					}
				}
				#divide by number of iterations to get pvalue
				pvalue[p] <- (count + 1)/(num_iterations + 1)
			}
			close(pb)
			names(pvalue) <- genes
			data <- t(rbind(network_betweenness, pvalue))
			data[,1] <- as.numeric(format(round(data[,1], 5), nsmall = 5))
			colnames(data)[2] <- "betweenness_pvalue"
			last[[i]] <- cbind(last[[i]], data)
			
			if(any(data[,2] <= sig_value)){
				sig_genes <- rownames(data[data[,2] <= sig_value,])
				print(paste("The", colors[i], "module has", length(sig_genes) ,"significant betweenness hubs"))
				plotting <- list()
				if(length(sig_genes > 0)){
					for(a in 1:length(sig_genes)){
						plotting[[a]] <- network_betweenness[sig_genes[a]]
						for(z in 1:num_iterations){
							point <- results_between[[z]][sig_genes[a]]
							plotting[[a]] <- c(plotting[[a]], point)
							#table(plotting)
						}
					}
					assign(paste(colors[i], "_sig_genes_betweenness", sep = ""), plotting)
					nam <- paste(colors[i], "_sig_genes_betweenness", sep = "")
					save(list = nam, file = paste(nam, "RData", sep = "."))
				}
			}
		}

		##########################
		#if hub score
		##########################
		#for each gene do following:
		if(hubscore){
			pb <- txtProgressBar(min = 0, max = length(genes), width = 60, style = 3)
			for (p in 1:length(genes)){
				setTxtProgressBar(pb, p)
				count = 0
				#check if the hub score is higher than the original score and if is then add a count
				for (l in 1:num_iterations){
					if (results_hubscore[[l]][genes[p]] >= network_hubscore[genes[p]]){
						count = count + 1
					}
				}
				#divide by number of iterations to get pvalue
				pvalue[p] <- (count + 1)/(num_iterations + 1)
			}
			close(pb)
			names(pvalue) <- genes
			data <- t(rbind(network_hubscore, pvalue))
			data[,1] <- as.numeric(format(round(data[,1], 5), nsmall = 5))
			colnames(data)[2] <- "hubscore_pvalue"
			last[[i]] <- cbind(last[[i]], data)

			if(any(data[,2] < sig_value)){
				sig_genes <- rownames(data[data[,2] < sig_value,])
				print(paste("The", colors[i], "module has", length(sig_genes) ,"significant hub score hubs"))
				plotting <- list()
				if(length(sig_genes > 0)){
					for(a in 1:length(sig_genes)){
						plotting[[a]] <- network_hubscore[sig_genes[a]]
						for(z in 1:num_iterations){
							point <- results_hubscore[[z]][sig_genes[a]]
							plotting[[a]] <- c(plotting[[a]], point)
							#table(plotting)
						}
					}
					assign(paste(colors[i], "_sig_genes_hubscore", sep = ""), plotting)
					nam <- paste(colors[i], "_sig_genes_hubscore", sep = "")
					save(list = nam, file = paste(nam, "RData", sep = "."))
				}
			}			
		}
		######################
		#pagerank
		######################
		#for each gene do following:
		if(pagerank){
			pb <- txtProgressBar(min = 0, max = length(genes), width = 60, style = 3)
			for (p in 1:length(genes)){
				setTxtProgressBar(pb, p)
				count = 0
				#check if the hub score is higher than the original score and if is then add a count
				for (l in 1:num_iterations){
					if (results_pagerank[[l]][genes[p]] >= network_pagerank[genes[p]]){
						count = count + 1
					}
				}
				#divide by number of iterations to get pvalue
				pvalue[p] <- (count + 1)/(num_iterations + 1)
			}
			close(pb)
			names(pvalue) <- genes
			data <- t(rbind(network_pagerank, pvalue))
			data[,1] <- as.numeric(format(round(data[,1], 5), nsmall = 5))
			colnames(data)[2] <- "pagerank_pvalue"
			last[[i]] <- cbind(last[[i]], data)
			if(any(data[,2] <= sig_value)){
				sig_genes <- rownames(data[data[,2] <= sig_value,])
				print(paste("The", colors[i], "module has", length(sig_genes) ,"significant pagerank hubs"))
				plotting <- list()
				if(length(sig_genes > 0)){
					for(a in 1:length(sig_genes)){
						plotting[[a]] <- network_pagerank[sig_genes[a]]
						for(z in 1:num_iterations){
							point <- results_pagerank[[z]][sig_genes[a]]
							plotting[[a]] <- c(plotting[[a]], point)
							#table(plotting)
						}
					}
					assign(paste(colors[i], "_sig_genes_pagerank", sep = ""), plotting)
					nam <- paste(colors[i], "_sig_genes_pagerank", sep = "")
					save(list = nam, file = paste(nam, "RData", sep = "."))
				}
			}
		}

		#######################
		#closeness
		#######################
		#for each gene do following:
		if(closeness){
			pb <- txtProgressBar(min = 0, max = length(genes), width = 60, style = 3)
			for (p in 1:length(genes)){
				setTxtProgressBar(pb, p)
				count = 0
				#check if the hub score is higher than the original score and if is then add a count
				for (l in 1:num_iterations){
					if (results_closeness[[l]][genes[p]] >= network_closeness[genes[p]]){
						count = count + 1
					}
				}
				#divide by number of iterations to get pvalue
				pvalue[p] <- (count + 1)/(num_iterations + 1)
			}
			close(pb)
			names(pvalue) <- genes
			data <- t(rbind(network_closeness, pvalue))
			data[,1] <- as.numeric(format(round(data[,1], 5), nsmall = 5))
			colnames(data)[2] <- "closeness_pvalue"
			last[[i]] <- cbind(last[[i]], data)

			if(any(data[,2] <= sig_value)){
				sig_genes <- rownames(data[data[,2] <= sig_value,])
				print(paste("The", colors[i], "module has", length(sig_genes) ,"significant closeness hubs"))
				plotting <- list()
				if(length(sig_genes > 0)){
					for(a in 1:length(sig_genes)){
						plotting[[a]] <- network_closeness[sig_genes[a]]
						for(z in 1:num_iterations){
							point <- results_closeness[[z]][sig_genes[a]]
							plotting[[a]] <- c(plotting[[a]], point)
							#table(plotting)
						}
					}
				}
				assign(paste(colors[i], "_sig_genes_closeness", sep = ""), plotting)
				nam <- paste(colors[i], "_sig_genes_closeness", sep = "")
				save(list = nam, file = paste(nam, "RData", sep = "."))
			}
		}

		#######################
		#module membership
		########################
		#for each gene do following:
		if(MM){
			pb <- txtProgressBar(min = 0, max = length(genes), width = 60, style = 3)
			for (p in 1:length(genes)){
				setTxtProgressBar(pb, p)
				count = 0
				#check if the hub score is higher than the original score and if is then add a count
				for (l in 1:num_iterations){
					if (results_MM[[l]][genes[p]] >= network_MM[genes[p]]){
						count = count + 1
					}
				}
				#divide by number of iterations to get pvalue
				pvalue[p] <- (count + 1)/(num_iterations + 1)
			}
			close(pb)
			names(pvalue) <- genes
			data <- t(rbind(network_MM, pvalue))
			data[,1] <- as.numeric(format(round(data[,1], 5), nsmall = 5))
			colnames(data)[2] <- "MM_pvalue"
			last[[i]] <- cbind(last[[i]], data)

			if(any(data[,2] <= sig_value)){
				sig_genes <- rownames(data[data[,2] <= sig_value,])
				print(paste("The", colors[i], "module has", length(sig_genes) ,"significant MM hubs"))
				plotting <- list()
				if(length(sig_genes > 0)){
					for(a in 1:length(sig_genes)){
						plotting[[a]] <- network_MM[sig_genes[a]]
						for(z in 1:num_iterations){
							point <- results_MM[[z]][sig_genes[a]]
							plotting[[a]] <- c(plotting[[a]], point)
							#table(plotting)
						}
					}
				}
				assign(paste(colors[i], "_sig_genes_MM", sep = ""), plotting)
				nam <- paste(colors[i], "_sig_genes_MM", sep = "")
				save(list = nam, file = paste(nam, "RData", sep = "."))
			}
		}
		#########################
		#edge betweenness
		#########################
		if(edge){
			if(limit_edge){
				pairs <- names(network_edgebetweenness[order(network_edgebetweenness, decreasing = TRUE)][1:length(genes)])
			}else{
				pairs <- names(network_edgebetweenness)
			}
			pb <- txtProgressBar(min = 0, max = length(pairs), width = 60, style = 3)
			for (p in 1:length(pairs)){
				setTxtProgressBar(pb, p)
				count = 0
				#check if the hub score is higher than the original score and if is then add a count
				for (l in 1:num_iterations){
					if (results_edgebetweenness[[l]][pairs[p]] >= network_edgebetweenness[pairs[p]]){
						count = count + 1
					}
				}
				#divide by number of iterations to get pvalue
				pvalue[p] <- (count + 1)/(num_iterations +  1)
			}
			close(pb)
			names(pvalue) <- pairs
			if(limit_edge){
				data <- t(rbind(network_edgebetweenness, pvalue))
				data[,2] <- as.numeric(format(round(data[,2], 5), nsmall = 5))
				data <- data[pairs,]
				colnames(data)[2] <- "edge_pvalue"
				last[[i]] <- cbind(last[[i]], pairs, data)

				if(any(data[,2] <= sig_value)){
					sig_genes <- rownames(data[data[,2] <= sig_value,])
					print(paste("The", colors[i], "module has", length(sig_genes) ,"significant edges"))
					plotting <- list()
					if(length(sig_genes > 0)){
						for(a in 1:length(sig_genes)){
							plotting[[a]] <- network_edgebetweenness[sig_genes[a]]
							for(z in 1:num_iterations){
								point <- results_edgebetweenness[[z]][sig_genes[a]]
								plotting[[a]] <- c(plotting[[a]], point)
								#table(plotting)
							}
						}
					}
					assign(paste(colors[i], "_sig_genes_edgebetweenness", sep = ""), plotting)
					nam <- paste(colors[i], "_sig_genes_edgebetweenness", sep = "")
					save(list = nam, file = paste(nam, "RData", sep = "."))
				}

			}else{
				edge_results <- t(rbind(network_edgebetweenness, pvalue))
				edge_results[,2] <- as.numeric(format(round(edge_results[,2], 5), nsmall = 5))
				colnames(edge_results)[2] <- "edge_pvalue"
				edge_values[[i]] <- edge_results
			}
		}
	}
	if (edge == T & limit_edge == F){
		last <- list(hub_detection = last, edge_values = edge_values)
	}
	return(last)
}
