igraph2nea.net <- function (NET) {
	require("igraph");
	edges <- E(NET); 
	elist <- incident_edges(NET, v = V(NET));
	Net <- as.list(NULL);
	Net$links <- as.list(NULL);
	for (vv in 1:length(V(NET))) {
		# print(V(NET)[.nei(c(vv))]);
		# Net$links[[vv]] <- as.character(unique(V(NET)[.nei(c(vv))]));
		Net$links[[vv]] <- unique(V(NET)[.nei(c(vv))]);
	}
	Net$Ntotal <- sum(sapply(Net$links, length))/2;
	return(Net);
}