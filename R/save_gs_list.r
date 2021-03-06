#' Create a TAB-delimited text file from AGS or FGS
#'
#' Each line in this file represents one gene/protein from an AGS/FGS and is accompanied with respective AGS/FGS ID. This format can be used e.g. as input at web site \href{https://www.evinet.org/}{EviNet}

#' @param gs.list a list created with \code{\link{samples2ags}}, \code{\link{mutations2ags}}, \code{\link{as_genes_fgs}}, or \code{\link{import.gs}}.
#' @param File output file name.

#' @seealso \code{\link{samples2ags}}, \code{\link{mutations2ags}}, \code{\link{as_genes_fgs}}, \code{\link{import.gs}}
#' @references \url{http://www.biomedcentral.com/1471-2105/13/226}
#' @references \url{https://www.evinet.org/}

#' @examples
#' netpath <- system.file("extdata", "merged6_and_wir1_HC2", package="NEArender")
#' net <- import.net(netpath);
#' fgs.genes <- as_genes_fgs(net);
#' save_gs_list(fgs.genes, File = "single_gene_ags.groups.tsv");
#' @importFrom utils write.table
#' @export


save.gs.list <- function(gs.list, ID = "gs", attributes = FALSE, shape = NA) { 
	t1 <- NULL;
	for (gs in names(gs.list)) {
		if (attributes) {
			Values <- gs.list[[gs]];
			# shape = "diamond";
			t1 <- rbind(t1, cbind(names(Values), 
				sign(Values) * Values / ifelse(Values >= 0, max(Values, na.rm = TRUE), min(Values, na.rm = TRUE)),
				shape,
				gs)); 
		} else {
			t1 <- rbind(t1, cbind(gs.list[[gs]], gs));
		}
	}
	write.table(t1, file = paste("NEA", ID, ifelse(attributes, 'attributes', ''), "groups", sep = "."), append = FALSE, quote = FALSE, 
		sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE,  col.names = FALSE);
}