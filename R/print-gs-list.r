print.gs.list <- function(gs.list, file = "gs.list.groups") {
	t1 <- NULL;
	for (gs in names(gs.list)) {
		if (length(gs.list[[gs]]) > 0) {
			t1 <- rbind(t1, cbind(gs.list[[gs]], gs));
		}
	}
	write.table(t1, file=file, append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE,  col.names = FALSE)
}