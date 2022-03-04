#' ROC for NEA benchmarks
#' 
#' Plot ROC curve(s) for benchmarked network(s)

#' @param tpvsfp a list contaning data for the ROC curve(s) from \code{\link{benchmark}}: \code{cross.z}, \code{cutoffs}, \code{ne}, \code{nv},\code{tp}, \code{fp}, \code{tn}, \code{fn}.
#' @param mode either ROC curve (roc) or precision recall curve (prc). default: roc
#' @param tpr.stop upper limit for the Y-axis (true positives). default: 1
#' @param fpr.stop upper limit for the X-axis (false prositives), default: 1
#' @param coff.z the point where to stop the ROC curve (in order to disable, set coff.z to below the minimal possible level)
#' @param coff.fdr the point where to put the circle ('cross') at the ROC curve 
#' @param main title name for the plot, default: none
#' @param cex.main font size for the main title
#' @param cex.leg font size for the legend
#' @param use.lty if different line types should be used for the curves, useful when the number of curves is very big (>10). default:"FALSE"
#' @param print.stats add N(edges) and N(vertices) to the legend lines; only works if non-empty $ne and $nv elements are submitted in tpvsfp. default:"TRUE"
#' @param sort_by_letter_and_remove enable a specific order for the curves in the legend, set by prefix letters \code{gsub("[A-Z]\\:", "", p1, fixed=F, ignore.case=T)}, such as "a:network1", "b:network-2" etc., whereby a: and b: are removed from the text strings. default:"FALSE"
#' 
#' @details The function \code{\link{roc}} can be called either from inside \code{\link{benchmark}}, or separately. In the latter case, it should receive the first argument as an object of the same type as \code{\link{benchmark}} can return. Generally, a ROC curve encompasses the whole interval from 0\% to 100\% at both axes. However in case of predictions made via network analysis, only the interval with (at least) formally significant scores is of interest (\href{https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-15-308}{Merid et.al 2014}).
#' Therefore the benchmark results are presented as ROC curves that end at minimal acceptable confidence, set via either a z-score cut-off. Additionally, an extra point is denoted that might correspond to \code{cross.z} in \code{\link{roc}}. The scale is given in no. of tested member genes rather then as percentage of the total no. of tests. The scale ranges \code{xlim} and \code{ylim} of the plot can be reduced by submitting parameters \code{tpr.stop}  and \code{fpr.stop} with values <1.
#' 
#' @seealso \link{benchmark}
#' @references \url{http://www.biomedcentral.com/1471-2105/15/308}
#' @keywords benchmark

#' @examples
#' # Benchmark and plot one networks on the whole set of test GSs, using no mask:
#' fpath <- system.file("extdata", "CAN_SIG_GO.34.txt", package="NEArender")
#' gs.list <- import.gs(fpath, Lowercase = 1, col.gene = 2, col.set = 3);
#' netpath <- system.file("extdata", "merged6_and_wir1_HC2", package="NEArender")
#' net <- import.net(netpath)
#' \donttest{
#' b0 <- benchmark (NET = net,
#'  GS = gs.list,
#'  echo=1, na.replace = 0, mask = ".", minN = 0,
#'  coff.z = 1.965, coff.fdr = 0.1, Parallelize=2);
#'  roc(b0, coff.z = 1.64);
#' }
#' \dontrun{
#' ## Benchmark and plot a number of networks on GO terms and KEGG pathways separately, using masks
#' b1 <- NULL;
#' for (mask in c("kegg_", "go_")) {
#' b1[[mask]] <- NULL;
#' for (file.net in c("netpath")) {
#' # a series of networks can be put here: c("netpath1", "netpath2", "netpath3")
#' net <- import.net(netpath, col.1 = 1, col.2 = 2, Lowercase = 1, echo = 1)
#' b1[[mask]][[file.net]] <- benchmark (NET = net, GS = gs.list,
#' gs.gene.col = 2, gs.group.col = 3, net.gene1.col = 1, net.gene2.col = 2,
#' echo=1, na.replace = 0, mask = mask, minN = 0, Parallelize=2);
#' }}
#' par(mfrow=c(2,1));
#' roc(b1[["kegg_"]], coff.z = 2.57,main="kegg_");
#' roc(b1[["go_"]], coff.z = 2.57,main="go_");
#' }
#'
#' @importFrom grDevices rainbow
#' @importFrom graphics abline legend points
#' @export


roc <- function(tpvsfp, # object of a special type; must contain "cross.z", "cutoffs", "ne", "nv", "tp", "fp", "tn", "fn"
	mode = "roc", # either ROC curve or precision recall curve (mode='prc'), default: roc
	tpr.stop = 1, # upper limit for the Y-axis (true positives), default: 1
	fpr.stop = 1, # upper limit for the X-axis (false positives), default: 1
	coff.z = 1.965, # the point where to stop the ROC curve (in order to disable, set coff.z to below the minimal possible level)
	# coff.fdr = 0.1, # the point where to put the circle ('cross') at the ROC curve 
	cex.main = 1, # font size for the main title
	cex.leg = 0.75, # font size for the legend
	main = NA, # title text, default: none
	use.lty = FALSE, # if different line types should be used for the curves, useful when the number of curves is very big (>10)
	print.stats = TRUE, # add N(edges) and N(vertices) to the legend lines; only works if non-empty $ne and $nv elements are submitted in tpvsfp
	sort_by_letter_and_remove = FALSE # enable a specific sorting order
) {
	if (mode(tpvsfp) != "list") {
		print('Submitted first agrument is wrong. Should be either a single list with entries c("cross.z", "cutoffs", "fp", "tp", "ne", "nv") or a list of such lists...');
		return();
	}

	# if (suppressWarnings(all(sort(names(tpvsfp)) == c("cross.z", "cutoffs", "ne", "nv", "tp", "fp", "tn", "fn")))) {
	if (!FALSE %in% (names(tpvsfp) %in% c("cross.z", "cutoffs", "ne", "nv", "tp", "fp", "tn", "fn"))) {
		objs <- NULL; objs[["1"]]=tpvsfp;
	} else {
		objs=tpvsfp;
	} 
	Nets <- sort(names(objs));

	if (length(Nets) > 1 | Nets[1] != "1") {
		Col <- rainbow(length(Nets)); names(Col) <- Nets;
	} else {
		Col <- NULL; Col["1"] <- c("black");
	} 
	Max.fp <- 0; Max.tp <- 0;
	for (p1 in Nets) {
		pred <- objs[[p1]];
		Max.signif.n = max(pred$fp[[1]][which(pred$cutoffs[[1]] > coff.z)], na.rm = TRUE);
		if (Max.fp < Max.signif.n ) {Max.fp <- Max.signif.n;}
		Max.signif.p = max(pred$tp[[1]][which(pred$cutoffs[[1]] > coff.z)], na.rm = TRUE);
		if (Max.tp < Max.signif.p ) {Max.tp <- Max.signif.p;}
	}
	
	if (mode == "roc") {
		Xlim = c(0, Max.fp * fpr.stop); 
		Ylim = c(0, Max.tp * tpr.stop);
		Xlab = "False predictions";
		Ylab = "True predictions";
	} else {
		Xlim = c(0, 1); 
		Ylim = c(0, 1);
		Xlab = 'Recall, TP / (TP + FN)';
		Ylab = 'Precision, TP / (TP + FP)';
	}

	plot(1, 1, xlim=Xlim, ylim=Ylim, cex.main = cex.main, type = "n", xlab = Xlab, ylab = Ylab, main = main);
	Lty = rep(1:4, times = length(Nets));
	i = 1; Le <- NULL;
	for (p1 in Nets) {
		if (print.stats) {
			# Le = c(Le, paste(substr(p1, 1, 15), "; n(V)=", sapply(objs, function (x) x$nv), "; n(E)=", sapply(objs, function (x) x$ne), sep=""));
			Le = c(Le, paste(substr(p1, 1, 15), "; n(V)=", objs[[p1]]$nv, "; n(E)=", objs[[p1]]$ne, sep = ""));
		} else {
			if (sort_by_letter_and_remove) {
				Le = c(Le, paste(gsub("[A-Z]\\:", "", p1, fixed = FALSE,  ignore.case = TRUE), sep = ""));
			} else {
				Le = c(Le, paste(p1, sep = ""));
			}
		}
		pred <- objs[[p1]];
	
		if (mode == "roc") {
			pick <- which(pred$cutoffs[[1]] > coff.z);
			X = pred$fp[[1]][pick];
			Y = pred$tp[[1]][pick];
		} else {
			if (mode == "prc") {
				X = pred$tp[[1]] / (pred$tp[[1]] + pred$fn[[1]]);
				Y = pred$tp[[1]] / (pred$fp[[1]] + pred$tp[[1]]);
			} else {
				stop("Curve mode undefined...");
			}
		}

		points(X,Y, pch=".", type="l", col=Col[p1], lty=ifelse(use.lty, Lty[i], 1));
		abline(0, 1,col = "gray60",lwd=0.5, lty=2);
		Pos <- pred$cutoffs[[1]][which(pred$cutoffs[[1]] > pred$cross.z)];
		points(
			pred$fp[[1]][which((Pos - pred$cross.z) == min(Pos - pred$cross.z, na.rm = TRUE))],
			pred$tp[[1]][which((Pos - pred$cross.z) == min(Pos - pred$cross.z, na.rm = TRUE))], 
			pch = "O", type = "p", col = Col[p1]);
		i = i + 1;
	}
	llist <- 1;
	if (use.lty) {llist <- Lty[1:i];}
	legend(x = "bottomright", legend = Le, title = paste( "", sep = ""), col = Col, bty = "n", pch = 15, lty = llist, cex = cex.leg);
}