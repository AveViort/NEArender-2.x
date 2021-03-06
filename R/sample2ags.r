#' Create AGS from a raw data matrix.
#'
#' This function creates a number of new AGSs from a given dataset, such as gene copy number, gene/protein expression, gene methylation etc. Such matrices M typically have size Ngenes x Nsamples, so that the current function returns a list of \code{length=ncol(M)}. The AGSs for each of the Nsamples are created with one of the five available methods (see parameter \code{method}).

#' @param method Method to select sample-specific genes. One of
#' \itemize{
#' \item "significant" : Using one-sided z-test, select gene values of which in the given sample \emph{i} deviate from the mean over all the samples, requiring q-value (Benjamini-Hochberg FDR) be below \emph{cutoff.q}.(default)
#' \item "topnorm"     : similar to "significant", i.e. calculates \eqn{(x[i] - mean(x)) / sd(x)} but does not evaluate significance. Instead, top N ranked genes \code{Ntop} are taken into AGS.
#' \item "top"         : similar to "topnorm", but \eqn{x[i] - mean(x)} is not divided with \emph{sd(x)}. This might help to prioritize genes with higher \emph{mean(x)} and ignore ones with low signal. Consider also that AGSs from "top" overlap much more with each other than those from "topnorm", i.e. would be less sample-specific.
#' \item "toppos"     : similar to "top", but retrives only genes with positive values of \eqn{x[i] - mean(x)}. This might be useful when the gene expression values are small counts (such as in single-cell RNA sequencing), so that considering the left part of the distribution would not bring high-quality AGS.
#' \item "toprandom"     : generates lists of Ntop random genes for each AGS.
#' }
#' @param m0 input matrix.
#' @param col.mask To include only columns with IDs that contain the specified mask. Follows the regular expression syntax.
#' @param namesFromColumn Number of the column (if any) that contains the gene/protein names. Note that it is only necessary of the latter  are NOT the unique rownames of the matrix. This could be sometimes needed to be able to process redundant expression etc. profiles.
#' @param cutoff.q cutoff value. Mutually exclusive with "Ntop" (default: 0.05)
#' @param Ntop Number of top ranking genes to include into each sample-specific AGS. Mutually exclusive with "cutoff.q". A practically recommended value of Ntop could be in the range 30...300. Ntop>1000 might decrease the analysis specificity.
#'
#' @examples
#' data("fantom5.43samples", package = "NEArender")
#' ags.list <- samples2ags(fantom5.43samples, cutoff.q = 0.01, method = "significant")
#' @export



samples2ags <- function(m0, Ntop = NA, col.mask = NA, namesFromColumn = NA, method = c("significant", "top", "toppos", "topnorm", "bestp", "toprandom"), cutoff.q = 0.05, 
	return.value = FALSE, gs.list = NA #these two options enable a special mode: extract expression values for pre-specified gene set(s) and use them as gene attributes on evinet.org for node coloring
) {
	if (!method %in% c("topnorm") & return.value == TRUE) {
		stop(paste("If 'return.value' is TRUE, then parameter 'method' should only be: [", paste(c("topnorm"), collapse = ", "), "]. Terminated...", sep = ""));
	}
	
	if (return.value == TRUE & !is.list(gs.list)) {
		stop(paste("If 'return.value' is TRUE, then parameter 'gs.list' should refer to a pre-compiled list of gene sets. Terminated..."));
	}
	
	if (return.value == FALSE & is.list(gs.list)) {
		stop(paste("If 'return.value' is FALSE, then parameter 'gs.list' is irrelevant. Terminated..."));
	}
	
	if (!method %in% c("bestp", "significant", "top", "toppos", "topnorm", "toprandom")) {
		print(paste("Parameter 'method' should be one of: ", paste(c("significant", "top", "topnorm", "toprandom"), collapse=", "), ". The default is 'method = significant'...", sep=""));
	}
	
	if (length(method) < 1 ) {
		print("Parameter 'method' is undefined. The default 'method = significant' will be used.");
		method = "significant";
	}

	if (method == "significant" & !is.na(Ntop)) {
		print("Parameter 'Ntop' is irrelevant when 'method = significant'. Terminated...");
	}
	
	if (return.value == TRUE & !is.na(Ntop)) {
		print("Parameter 'Ntop' is irrelevant when 'return.value' is TRUE.");
	}
	
	if (grepl("top", method, ignore.case = TRUE) & is.na(Ntop) & !return.value) {
		stop("Parameter 'Ntop' is missing...");
	}
	
	ags.list <- NULL; 
	if (is.na(namesFromColumn)) {m1 <- m0;} else {m1 <- m0[,(namesFromColumn+1):ncol(m0)];}
	if (!is.na(col.mask)) {m1 <- m1[,colnames(m1)[grep(col.mask,colnames(m1))]];}
	
	if (method != "bestp") {
		uc <- sweep(m1, 1, rowMeans(m1, na.rm = TRUE), FUN = "-");

		if (method == "toprandom") {
			uc <- m1;
			for (label in colnames(uc)) {
				x = uc[,label];
				ags.list[[label]] <- sample(tolower(names(x)), Ntop);
			}
		}
	} else {
		for (label in colnames(m1)) {
			x = m1[,label];
			ags.list[[label]] <- tolower(names(x))[order(x, decreasing = FALSE)][1:Ntop];
		}
	}

	if (method %in% c("significant", "topnorm")) {
		SD <- apply(m1, 1, sd, na.rm = TRUE);
		uc <- sweep(uc, 1, SD, FUN = "/");
		if (method == "significant") {
			if (!return.value) {
				p1 <- 2 * pnorm(abs(uc), lower.tail = FALSE);
				q1 <- apply(p1, 2, function (x) p.adjust(x, method = "BH"));
				ags.list <- apply(q1, 2, function (x) tolower(names(x))[which(x < cutoff.q)]);
			}
		}
	}
	
	if (method %in% c("top", "toppos", "topnorm")) {
		if (!return.value) {
			for (label in colnames(uc)) {
				if (method == "toppos") {
					x = uc[,label];
				} else {
					x = abs(uc[,label]);
				}
				ags.list[[label]] <- tolower(names(x))[order(x, decreasing = TRUE)][1:Ntop];
			}
		}
	}

	if (!return.value) {
		if (length(ags.list) < 1) {
			print("Could not produce any gene sets under the current parameters: method=", method, "; Ntop=", Ntop, "; cutoff.q=", cutoff.q, sep = "");
			return(NULL);
		}
		return(ags.list); 
	}  else {
		gs.values <- as.list(NULL);
		rownames(uc) <- tolower(rownames(uc));
		if ("FALSE" %in% names(table(names(gs.list) %in% colnames(uc)))) {
			for (label in colnames(uc)) {
				gs.values[[label]] <- as.list(NULL);
				for (gs in names(gs.list)) {
					gs.values[[label]][[gs]] <- rep(NA, times = length(gs.list[[gs]]));
					names(gs.values[[label]][[gs]]) <- gs.list[[gs]];
					sharedGenes <- intersect(rownames(uc), names(gs.values[[label]][[gs]]));
					gs.values[[label]][[gs]][sharedGenes] <- uc[sharedGenes,label];
				}
			} 
		} else {
			for (gs in names(gs.list)) {
				gs.values[[gs]] <- rep(NA, times = length(gs.list[[gs]]));
				names(gs.values[[gs]]) <- gs.list[[gs]];
				sharedGenes <- intersect(rownames(uc), names(gs.values[[gs]]));
				gs.values[[gs]][sharedGenes] <- uc[sharedGenes,gs];
			}
		}
		return(gs.values); 
	}
}
