require(BiocGenerics)
require(AnnotationDbi)
require(GSEABase)
require(parallel)
require(Category)
require(GOstats)

################################################################################
# INSTRUCTIONS
# source("module.enrichment.analysis.withDependencies.R")
# see examples above module.enrichment.analysis
# Note this assumes that the GeneSetCollection object, uses the primary 'gene ID' type.
# for human's, this is the Entrez Gene ID 
# To convert from gene symbols to entrez gene ID's:
# library(GSEABase)
# gsc.eg <- mapIdentifiers(gsc, SymbolIdentifier("org.Mm.eg.db"), EntrezIdentifier("org.Mm.eg.db"))
# 
# To convert from ENSEMBL ID's to entrez gene ID's:
 require(BiocGenerics)
 require(AnnotationDbi)
 require(GSEABase)
 require(parallel)
 require(Category)
 require(GOstats)
# ens.ids <- readLines("SP0ENSG_cluster1.txt")
# ens.ids.clean <- sub("\\..+", "", ens.ids)
# # cbind(ens.ids, ens.ids.clean)
# gs <- GeneSet(ens.ids.clean, geneIdType=ENSEMBLIdentifier("org.Mm.eg.db"), setName="SP0ENSG_cluster1")
# # gs
# # didn't work.... why not?
# # mapIdentifiers(gs, ENSEMBLIdentifier("org.Mm.eg.db"), EntrezIdentifier("org.Mm.eg.db"))
# library(org.Mm.eg.db)
# get("ENSG00000213558", org.Mm.egENSEMBL2EG)
# get("ENSG00000213558", org.Mm.egENSEMBL2EG)
# get("ENSG00000102053", org.Mm.egENSEMBL2EG)
# mget(ens.ids.clean, org.Mm.egENSEMBL2EG, ifnotfound=NA)
# gs
# geneIdType(gs)
# geneIdType(gs) <- EntrezIdentifier("org.Mm.eg.db")
# gs
# history(100)
# 
################################################################################

###
### from /Users/marcow/src/R/enrichR/dev/bin/module.enrichment.analysis.R
###

#' From a list of modules, run an enrichment analysis.
#' 
#' The input is a GeneSetCollection object, where each GeneSet has a unique name, and
#' uses the primary annotation type, which is mostly the Entrez Gene ID, or the sgd ID for 
#' yeast, tair ID for arabidopsis.
#' Given an outdir, this function writes tsv files, 1 per GeneSet, within files called:\cr
#' <outdir>/<test.type>/<GeneSet.name>-<test.type>-<test.direction>-enrichment.tsv. It's thus in your interest
#' to name you GeneSet's usefully (hint: \code{GSEABase::\link{setName}}).
#' 
#' @section Development/Speedup:
#' To speed things up, you can specify a snow, or parallel cluster, which will run 1 GeneSet per
#' node. At the time of writing, on OSX, I get intermittent reliability using snow SOCK cluster
#' (ie I get unserialize errors); rpvm is not available for OSX, so snow PVM cluster is out of the
#' question. Rmpi needs systemwide MPI installation, ... parallel's SOCK cluster just uses snow's
#' SOCK cluster, so that doesn't work\dots Only parallel's FORK cluster is working for me.
#'
#' @param modules a GeneSetCollection where each GeneSet uses the primary ID type for that species.
#'   This is usually a numerical Entrez Gene ID, but for yeast, it's an SGD ID.
#' @param outdir the output directory
#' @param snowCluster a snow cluster. To use multiple cores, we implemented snow cluster. 
#'  This can be \code{NULL}, in which case it will run single threaded.
#' @param test.type one of \dQuote{BP}, \dQuote{CC}, \dQuote{MF}, \dQuote{KEGG}, \dQuote{PFAM},
#'   or one of the Broad collections, eg \dQuote{C2CGP}, \dots
#' @param test.direction one of \dQuote{over}, \dQuote{under}
#' @param minGenesPerModule the minimum number of genes in each GeneSet to do enrichment upon
#' @param minGenesPerTerm the minimum number of genes within a GO/KEGG/PFAM/C2CGP/etc term to be considered.
#' @param p.adjust.method FDR/FWER correction method. default = \dQuote{BH}. see \code{\link{p.adjust}}
#' @param fdr.thresh fdr threshold for term enrichment, based upon \code{p.adjust.method}
#' @param minResultTerms if not enough terms are significant, then output the most 
#'    \code{minResultTerms} significant terms
#' @param msigdb.gsc a GeneSetCollection object from loading Broad genesets. 
#'     Only needed if \code{test.type} is a Broad type.
#' @return 
#'  Creates N files, containing enrichment results for each of the N qualifying modules within the AllegroFile.
#' a module qualifies if it has enough genes in it (see minGenesPerModule); 
#' Within each file, there's the results from M terms, where M is at least the top minResultTerms, or all 
#'  significant results which have FDR < \code{fdr.thresh}. 
#' Note, only the terms with sufficient sizes are
#' tested for enrichment (see minGenesPerTerm)
#' @author Mark Cowley, 2011-07-01
#' @examples
#' g1 <- c("1", "10", "100", "1000", "10000", "100008586", "10001", "10002", "10003", "100037264")
#' gs1 <- GeneSet(g1, geneIdType=EntrezIdentifier("org.Mm.eg.db"), setName="gs1")
#' g2 <- c("100037265", "10004", "100049542", "100049716", "10005", "10006", "10007", "10008", "10009", "100093630")
#' gs2 <- GeneSet(g2, geneIdType=EntrezIdentifier("org.Mm.eg.db"), setName="gs2")
#' gsc <- GeneSetCollection(list(gs1, gs2))
#' 
#' # the single-threaded way
#' module.enrichment.analysis(gsc, ".", NULL, test.type="BP")
#' 
#' # the multi-threaded way, via snow (NB: I get unserialize errors)
#' library(snow)
#' snowCluster <- makeSOCKcluster(rep("localhost", 2))
#' # snowCluster <- makeSOCKcluster(rep("localhost", getOption("Ncpus", 4)-2))
#' clusterExport(snowCluster, ls())
#' module.enrichment.analysis(gsc, ".", snowCluster, test.type="BP")
#' 
#' # the parallel way, via parallel (NB: seems to work better than snow's SOCK cluster.)
#' require(parallel)
#' # clu <- makeCluster(rep("localhost", 2))
#' clu <- makeForkCluster(4)
#' clusterExport(clu, ls())
#' module.enrichment.analysis(gsc, ".", clu, test.type="BP")
#' stopCluster(clu)
#' 
module.enrichment.analysis <- function(
	modules, 
	outdir=".", 
	snowCluster=NULL, 
	test.type=c(
		"BP", "MF", "CC", 
		"KEGG", 
		"PFAM", 
		"C2CGP",
		# "C1ALL", 
		# "C2ALL", "C2CGP", "C2CP", "C2BIOCARTA", "C2KEGG", "C2REACTOME", 
		# "C3ALL", "C3MIR", "C3TFT", 
		# "C4ALL", "C4CGN", "C4CM"
	),
	test.direction="over",
	minGenesPerModule=10, 
	minGenesPerTerm=10,
	p.adjust.method="BH",
	fdr.thresh=0.05,
	minResultTerms=10,
	DEBUG=FALSE,
	verbose=TRUE
	) {
		
	#
	# ERROR CHECKING
	#
	!missing(modules) || stop("You need to supply modules - a list of modules")
	is(modules, "GeneSetCollection") || stop("modules must be a GeneSetCollection")
	length(modules) > 0 || stop("the modules GeneSetCollection object was empty")
	all(annotation(modules) == annotation(modules)[1]) || stop("All GeneSet modules must use the same annotation package")
	!is.null(outdir) || stop("outdir must be specified")
	# !is.null(snowCluster) || stop('You need a snowCluster object. here\'s the simplest possible snowCluster: require(snow); snowCluster <- makeCluster(rep("localhost", 1), "SOCK")')
	test.type <- test.type[1]
	!missing(outdir)      || stop("You need to supply outdir - the output directory")
	test.direction %in% c("over", "under") || stop("test.direction must be one of over or under")
	is.numeric(minGenesPerModule) && minGenesPerModule >= 0 || stop("minGenesPerModule must be a positive integer or 0")
	is.numeric(minGenesPerTerm) && minGenesPerTerm >= 0 || stop("minGenesPerTerm must be a positive integer or 0")
	p.adjust.method %in% p.adjust.methods || stop("p.adjust.method must be one of the values within p.adjust.methods")
	is.numeric(fdr.thresh) && fdr.thresh >= 0.0 && fdr.thresh < 1.0 || stop("fdr.thresh must be in [0,1]")
	is.numeric(minResultTerms) && minResultTerms >= 0 || stop("minResultTerms must be a positive integer or 0")
	
	hgt <- list()
	EnrichmentResults <- list()
	files <- sprintf("%s/%s/%s-%s-%s-enrichment.tsv", outdir, test.type, setName(modules), test.type, test.direction)
	dir.create(dirname(files[1]), recursive=TRUE, mode="0770", showWarnings=FALSE)
	
	# for( i in seq(along=modules)) {
	a <- clusterApplyLB2(snowCluster, doit=!is.null(snowCluster), seq(along=modules), function(i) {
		module.enrichment.analysis.GeneSet(
			modules[[i]],
			outdir=outdir, 
			test.type=test.type,
			test.direction=test.direction,
			minGenesPerModule=minGenesPerModule, 
			minGenesPerTerm=minGenesPerTerm,
			p.adjust.method=p.adjust.method,
			fdr.thresh=fdr.thresh,
			minResultTerms=minResultTerms,
			DEBUG=DEBUG,
			verbose=verbose
		)
	})
	hgt <- sapply(a, "[", "hgt")
	EnrichmentResults <- sapply(a, "[", "EnrichmentResults")
	names(EnrichmentResults) <- names(hgt) <- names(files) <- setName(modules)

	idx <- !is.na(EnrichmentResults)
	EnrichmentResults <- subset(EnrichmentResults, idx)
	hgt               <- subset(hgt  , idx)
	files             <- subset(files, idx)
	
	res <- data.frame(geneSetName=character(), numTerms=numeric(), numSignificantTerms=numeric(), files=character(), stringsAsFactors=FALSE)
	if( length(EnrichmentResults) > 0 ) {
		res <- data.frame(
			geneSetName=names(files), 
			files=files, 
			test.type=test.type,
			numTerms=sapply(EnrichmentResults, nrow),
			numSignificantTerms=sapply(EnrichmentResults, function(x) sum(x$FDR < 0.05)),
			stringsAsFactors=FALSE
		)
	}
	invisible( res )
	###############################################################################
}

###
### from /Users/marcow/src/R/enrichR/dev/bin/module.enrichment.analysis.R
###

#' 
#' @examples
#' g1 <- c("1", "10", "100", "1000", "10000", "100008586", "10001", "10002", "10003", "100037264")
#' gs1 <- GeneSet(g1, geneIdType=EntrezIdentifier("org.Mm.eg.db"), setName="gs1")
#' g2 <- c("100037265", "10004", "100049542", "100049716", "10005", "10006", "10007", "10008", "10009", "100093630")
#' gs2 <- GeneSet(g2, geneIdType=EntrezIdentifier("org.Mm.eg.db"), setName="gs2")
#' module.enrichment.analysis.GeneSet(gs1, ".", test.type="BP")
module.enrichment.analysis.GeneSet <- function(
	gs,
	outdir=NULL, 
	test.type=c(
		"BP", "MF", "CC", 
		"KEGG", 
		"PFAM", 
		"C2CGP",
		# "C1ALL", 
		# "C2ALL", "C2CGP", "C2CP", "C2BIOCARTA", "C2KEGG", "C2REACTOME", 
		# "C3ALL", "C3MIR", "C3TFT", 
		# "C4ALL", "C4CGN", "C4CM"
	),
	test.direction="over",
	minGenesPerModule=10, 
	minGenesPerTerm=10,
	p.adjust.method="BH",
	fdr.thresh=0.05,
	minResultTerms=10,
	DEBUG=FALSE,
	verbose=TRUE) {
		
	require(Category) || stop("required package Category is not installed")
	hgt <- EnrichmentResults <- NA
	
	#
	# EXTRACT SOME DEFAULT ARGUMENTS
	#
	org.db = annotation(geneIdType(gs))
	require(org.db, character.only=TRUE) || stop("Required organism.db package not installed")

	debug2(GeneSet.hyperGTest, doit=DEBUG)
	
	outfile <- sprintf("%s/%s/%s-%s-enrichment.tsv", outdir, test.type, setName(gs), test.type)
	dir.create(dirname(outfile), recursive=TRUE, mode="0770", showWarnings=FALSE)
	
	if( numGenes(gs) < minGenesPerModule ) {
		if(verbose) {
			cat(sprintf("Skipping %s as it has too few genes (%d/%d)\n", setName(gs), numGenes(gs), minGenesPerModule))
		}
	} else {
		# 
		# calculate the enrichment vs all the terms in the term collection
		# 
		
		# errors can be thrown if there are no genes in any annotation terms (eg 5 genes with no GO term.)
		hgt <- tryCatch(
			GeneSet.hyperGTest(gs, 
				org.db, test.type=test.type, p.thresh=1.0, 
				test.direction=test.direction, 
				DEBUG=DEBUG
			),
			error=function(e) e # do nothing with the error.
		)
		# cat("Class of result: ", class(hgt), "\n")
		if( inherits(hgt, "error") ) {
			cat(sprintf("%s could not be annotated: %s\n", setName(gs), hgt$message))
			hgt <- NA
		} else if( identical(hgt, NA) ) {
			if(verbose) {
				cat(sprintf("%s could not be functionally annotated (ie had an NA hyperGTest result), even though it had these %d genes: %s\n", 
						setName(gs), numGenes(gs), paste(geneIds(gs), collapse=", ")))
			}
		} else {
			# calc FDR and filter results: include the top 10 enriched terms, and any others that have fdr<0.05
			EnrichmentResults <- HyperGResult2pina(
				hgt, 
				minGenesPerTerm=minGenesPerTerm, 
				minResults=minResultTerms, 
				fdr.thresh=fdr.thresh, 
				p.adjust.method=p.adjust.method
			)
			write.table(EnrichmentResults, outfile, sep="\t", col.names=T, quote=F, row.names=F)
			if(verbose) cat(sprintf("Wrote enrichment results for %s vs %s to %s\n", setName(gs), test.type, outfile))
		}
	}
	
	res <- list(hgt=hgt, EnrichmentResults=EnrichmentResults)
	return(res)
}



###
### from /Users/marcow/src/R/enrichR/dev/bin/clusterApply2.R
###

#' clusterApplyLB which can be turned on/off at will.
#' 
#' It's nigh impossible to debug code when its inside a clusterApplyLB
#' block. So, you end up editing code, flipping back & forth between
#' lapply and clusterApplyLB. This function flips the mode from using
#' clusterApplyLB or lapply with just a logical variable: \code{doit}
#' @param cl a SNOW cluster object
#' @param x an array
#' @param doit logical: if \code{TRUE}, use the snow cluster \code{cl}. if \code{FALSE},
#'   use good old \code{\link{lapply}}
#' @param fun function or character string naming a function
#' @param \dots further arguments passed to \code{fun}
#' @return the result of applying fun to each element of \code{x}
#' @author Mark Cowley, 2011-08-12
#' @export
clusterApplyLB2 <- function(cl, doit, x, fun, ...) {
	if( doit ) {
		clusterApplyLB(cl, x, fun, ...)
	}
	else {
		lapply(x, fun, ...)
	}
}

###
### from /Users/marcow/src/R/enrichR/dev/bin/debug2.R
###

#' Debug a function, controlled by a boolean
#' The standard \code{\link[base]{debug}} function works great, but this function lets you turn
#' debuging on and off via a logical global parameter, meaning no more code
#' commenting
#' 
#' @param fun: any interpreted R function. see \code{\link[base]{debug}}
#' @param text a text string that can be retrieved when the browser is entered.
#'   see \code{\link[base]{debug}}
#' @param condition a condition that can be retrieved when the browser is entered.
#'   see \code{\link[base]{debug}}
#' @param doit logical: if TRUE, then debugging mode will be entered, if \code{FALSE},
#'   no debugging will occur. Default=\code{TRUE} if running an interactive session.
#' @return
#' @seealso \code{\link[base]{debug}}
#' @author Mark Cowley, 2011-08-12
#' @export
debug2 <- function(fun, text="", condition=NULL, doit=interactive()) {
	if( doit )
		debug(fun=fun, text=text, condition=condition)
}

###
### from /Users/marcow/src/R/enrichR/dev/bin/GeneSet.hyperGtest.R
###

#' Run an enrichment analysis upon a GeneSet for GOBP, GOMF, GOCC, KEGG, PFAM
#'
#' This supports testing GeneSet's for enrichment vs Gene Ontology Biological Processes,
#' Molecular Function, Cellular Component, or KEGG pathways, or PFAM domains.
#' HyperG test on a GeneSet to either KEGG, GO, PFAM, or Broad GeneSetCollections
#' 
#' Broad's MSigDB\cr
#' i've created environments org.Mm.egC2CGP, org.Mm.egC2CGP, org.Rn.egC2CGP, and their
#' reverses: org.Mm.egC2CGP2EG, org.Mm.egC2CGP2EG, org.Rn.egC2CGP2EG.
#' I've also extended Category to provide a C2CGPHyperGParams-class.
#' So, once you've loaded the mapping environments into the current workspace, and
#' laoded my version of Category_2.18.0-2, you can run tests against C2CGP like you
#' would any of the other tests.\cr
#' I've also created the generic org.Mm.egMSIGDB mappings, but still working on
#' the Category code to allow for a MSIGDBHyperGParams-class, which can be subsetted
#' by geneset catgegory (eg C2CGP, like you can for the different clades of GO.)
#'
#' @param gs a GeneSet object from GSEABase
#' @param annot.db eg org.Mm.eg.db
#' @param test.clade one of "BP", "MF", "CC", "KEGG", "PFAM"
#' @param p.thresh a p-value threshold for the hyperGTest
#' @param test.direction are you looking for "over" or "under" enriched terms
#' @param GOtest.conditional A logical indicating whether the calculation should condition on the GO structure.
#' @param msigdb.gsc The MSigDB XML File as a GeneSetCollection, for use if test.type is a Broad test.
#' @param \dots: other arguments passed to hyperGTest.
#' @author Mark Cowley, 2011-07-01
#' @examples
#' 	library(GSEABase)
#' 	e1 <- GeneSet(c("9440","637","5713","9441","355","5430","9412","5433","6874","4478","8795","356","8772","841","5432","5518","5704","6881","5434","843","6837","11235","112950","29079","29888","7430","54797","6877","8428","5440","6882","6872","6873","4792","5431","4791","396","5435","85369","5599","6879","5436","84246","5441","5437","9861","80143","5515","81857","55588","9969","6908","7132","6878","6883","387","5710","219541","129685","5439","1147"), 
#' 	    geneIdType=EntrezIdentifier("org.Mm.eg.db"))
#' 	e1
#' 	hgt.gobp <- GeneSet.hyperGTest(e1, "org.Mm.eg.db", test.type="BP", p.thresh=1.0, test.direction="over")
#' 	hgt.kegg <- GeneSet.hyperGTest(e1, "org.Mm.eg.db", test.type="KEGG", p.thresh=1.0, test.direction="over")
#' 	hgt.pfam <- GeneSet.hyperGTest(e1, "org.Mm.eg.db", test.type="PFAM", p.thresh=1.0, test.direction="over")
#' 	hgt.c2cgp <- GeneSet.hyperGTest(e1, "org.Mm.eg.db", test.type="C2CGP", p.thresh=1.0, test.direction="over")
#' 	HyperGResult2pina(hgt.gobp)
#' 	HyperGResult2pina(hgt.kegg)
#' 	HyperGResult2pina(hgt.pfam)
#'  HyperGResult2pina(hgt.c2cgp)
#'  # f <- "msigdb_v3.xml"
#'  # msigdb3 <- getBroadSets("msigdb_v3.xml")
#'  # hgt.c2cgp <- GeneSet.c2cgp(e1, "", "c2cgp", msigdb.gsc=msigdbv3)
#' @export
GeneSet.hyperGTest <- function(gs, annot.db=NULL,
	test.type=c("BP", "MF", "CC", "KEGG", "PFAM",
		"C2CGP",
		# "C1ALL", "C2ALL", "C2CGP", "C2CP", "C2BIOCARTA", "C2KEGG", "C2REACTOME",
		# "C3ALL", "C3MIR", "C3TFT", "C4ALL", "C4CGN", "C4CM"
	),
	p.thresh=1.0, test.direction=c("over", "under"), GOtest.conditional=FALSE,
	msigdb.gsc, # @TODO: drop this parameter, when Category is updated.
	DEBUG=FALSE,
	...) {
	!is.null(annot.db) || stop("Must specify an annot.db. eg org.Mm.eg.db")
	require(GSEABase) || stop("required package 'GSEABase' is not installed")
	inherits(gs, "GeneSet") || stop("gs must be a GeneSet object (from GSEABase package)")
	require(GOstats) || stop("required package 'GOstats' is not installed")
	# require(annot.db, character=T)

	test.type <- test.type[1]
	test.direction <- test.direction[1]

	.get.genes.from.GeneSet <- function(x, is.entrezGeneID) {
		res <- x@geneIds[!is.na(x@geneIds)]
		res <- unique(res)
		# if these are Entrez Gene ID's, then as.numeric will work & this is the prefered class for these symbols;
		# if not Entrez Gene ID's (eg Yeast SGD ID's) then there will be NA's, and we should use the character IDs
		int <- as.numeric(res)
		if(all(!is.na(int))) res <- int
		res
	}
		
	res <- NA
	
	if( test.type %in% c("BP", "CC", "MF") ) { # @seealso help("GOHyperGParams-class")
		params <- new("GOHyperGParams",
			geneIds=.get.genes.from.GeneSet(gs),
			annotation=annot.db,
			ontology=test.type,
			pvalueCutoff=p.thresh,
			conditional=GOtest.conditional,
			testDirection=test.direction
		)
		res <- hyperGTest(params, ...)
		# @seealso help("GOHyperGResult-class"), help("HyperGResult-class")
	}
	else if( test.type == "KEGG" ) { # @seealso help("KEGGHyperGParams-class")
		params <- new("KEGGHyperGParams",
			geneIds=.get.genes.from.GeneSet(gs),
			annotation=annot.db,
			pvalueCutoff=p.thresh,
			testDirection=test.direction
		)
		res <- hyperGTest(params, ...)
		# @seealso help("GOHyperGResult-class"), help("HyperGResult-class")
	}
	else if( test.type == "PFAM" ) { # @seealso help("PFAMHyperGParams-class")
		params <- new("PFAMHyperGParams",
			geneIds=.get.genes.from.GeneSet(gs),
			annotation=annot.db,
			pvalueCutoff=p.thresh,
			testDirection=test.direction
		)
		res <- hyperGTest(params, ...)
		# @seealso help("GOHyperGResult-class"), help("HyperGResult-class")
	}
	else if( test.type == "C2CGP" ) { # @seealso help("C2CGPHyperGParams-class")
	params <- new("C2CGPHyperGParams",
			geneIds=.get.genes.from.GeneSet(gs),
			annotation=annot.db,
			pvalueCutoff=p.thresh,
			testDirection=test.direction
		)
		res <- hyperGTest(params, ...)
		# @seealso help("GOHyperGResult-class"), help("HyperGResult-class")
	}
	else if( substring(test.type,1,1) == "C" ) {
		# "C1ALL", "C2ALL", "C2CGP", "C2CP", "C2BIOCARTA", "C2KEGG", "C2REACTOME",
		# "C3ALL", "C3MIR", "C3TFT", "C4ALL", "C4CGN", "C4CM"
		types <- sapply(msigdb.gsc, function(elt) bcCategory(collectionType(elt)))
		stypes <- sapply(msigdb.gsc, function(elt) bcSubCategory(collectionType(elt)))
		idx <- switch(test.type,
			c1all=     types=="C1",
			c2all=     types=="C2",
			c2cgp=     types=="C2" & stypes=="CGP",
			c2cp=      types=="C2" & stypes=="CP",
			c2biocarta=types=="C2" & stypes=="CP:BIOCARTA",
			c2kegg=    types=="C2" & stypes=="CP:KEGG",
			c2reactome=types=="C2" & stypes=="CP:REACTOME",
			c3all=     types=="C3",          
			c3mir=     types=="C3" & stypes=="MIR",
			c3tft=     types=="C3" & stypes=="TFT",
			c4all=     types=="C4",          
			c4cgn=     types=="C4" & stypes=="CGN",
			c4cm=      types=="C4" & stypes=="CM",
			stop("error, unsupported test.type"))
		gsc <- msigdb.gsc[idx]
		debug2(hyperGtest.GeneSet.vs.GeneSetCollection, doit=DEBUG)
		try(res <- hyperGtest.GeneSet.vs.GeneSetCollection(gs, gsc, universeSize=NULL, testName="MSigDB"))
	}
	
	return( res )
}

###
### from /Users/marcow/src/R/enrichR/dev/bin/HyperGResult2pina.R
###

#' Convert a GOHyperGResult or HyperGResult object to a table
#'
#' Convert a GOHyperGResult or HyperGResult object to a table,
#'  with an FDR column, and limiting the number of significant terms.
#'
#' FDR correction: Since you generally test 1 geneset for significant enrichment vs
#' an entire collection of genesets, you should perform adjustments for multiple
#' testing. 
#' This is a bit trickier for GO-testing, since GO terms are hierarchical
#' and often parents and child terms can be extremely similar - adjusting for 2 
#' non-independent tests here may be a bit too strict. 
#' Within a HyperGResult object,
#' only those terms that had some overlap with \textit{your} genesets are included, 
#' ie, not all of the genesets that were tested. This is particularly important when 
#' doing a user query of say
#' 5 proteins vs a collection of modules, where we typically find only very few
#' categories with any overlap. In these cases, you should set \code{numTermsTested}
#' equal to the number of terms that were tested (after filtering for GeneSet size).
#'
#' @param x a GOHyperGResult, or HyperGResult object
#' @param minGenesPerTerm The minimum number of genes in a term. esp important for PFAM
#'    where many terms have few genes. Set to 0 to bypass this.
#' @param p.adjust.method The method to adjust for multiple testing. Default=\dQuote{BH}
#' @param numTermsTested The number of Terms that were tested for enrichment. This is
#'   occasionally larger than the number of results in \code{x}. 
#' If \code{NULL}, then this is 
#'  set to the number of terms in the \code{x} object. Otherwise, set to a large positive
#' integer. See details.
#' @param fdr.thresh the FDR threshold to use. set to \dQuote{NULL} to ignore.
#' @param minResults Non-NULL value of the minimum number of significant GO terms to
#'    return. useful if none of the results are < fdr.thresh
#' @param digits the number of significant figure to use. default=4.
#' @return a \code{data.frame} with these columns:
#' \item{Term}{the GO Term}
#' \item{Count}{the number of genes in your geneset and the GO term}
#' \item{Size}{The number of genes annotated to the GO term}
#' \item{Enrichment}{The ratio of enrichment vs random expectation}
#' \item{P}{the enrichment Pvalue}
#' \item{FDR}{the FDR of each GO term}
#' @examples
#' library(GSEABase)
#' e1 <- GeneSet(c("9440","637","5713","9441","355","5430","9412","5433","6874","4478","8795","356","8772","841","5432","5518","5704","6881","5434","843","6837","11235","112950","29079","29888","7430","54797","6877","8428","5440","6882","6872","6873","4792","5431","4791","396","5435","85369","5599","6879","5436","84246","5441","5437","9861","80143","5515","81857","55588","9969","6908","7132","6878","6883","387","5710","219541","129685","5439","1147"), 
#'     geneIdType=EntrezIdentifier("org.Mm.eg.db"))
#' e1
#' hgt.gobp <- GeneSet.hyperGTest(e1, "org.Mm.eg.db", test.type="BP",   p.thresh=1.0, test.direction="over")
#' hgt.kegg <- GeneSet.hyperGTest(e1, "org.Mm.eg.db", test.type="KEGG", p.thresh=1.0, test.direction="over")
#' hgt.pfam <- GeneSet.hyperGTest(e1, "org.Mm.eg.db", test.type="PFAM", p.thresh=1.0, test.direction="over")
#' HyperGResult2pina(hgt.gobp)
#' HyperGResult2pina(hgt.kegg)
#' HyperGResult2pina(hgt.pfam)
#' @export
HyperGResult2pina <- function(x, minGenesPerTerm=10, 
		p.adjust.method="BH", numTermsTested=NULL, fdr.thresh=0.05, minResults=10, digits=4) {
	inherits(x, "HyperGResultBase") || stop(sprintf("x must be an object that inherits from HyperGResultBase, and not: %s", class(x)))
	
	res <- summary(x)
	if( nrow(res) == 0 ) return(NA)
	res <- res[res$Size >= minGenesPerTerm, ]
	if( nrow(res) == 0 ) return(NA)
	
	if( is.null(numTermsTested) ) numTermsTested <- length(res$Pvalue)
	res$FDR <- p.adjust(res$Pvalue, p.adjust.method, n=numTermsTested)
	if( is.null(fdr.thresh) ) fdr.thresh <- 1.01

	# keep at least minResults rows, but include enough rows until the FDR threshold has been breached
	idx <- which(res$FDR<fdr.thresh | 1:nrow(res) <= minResults)
	res <- res[idx, ]

	
	terms <- ""
	if (colnames(res)[1] == "KEGGID") {
		# This was probably a KEGG analysis; prepend KEGG: to this column
		terms <- sprintf("%s (%s)", res$Term, paste("KEGG:", res[,1], sep=""))
	}
	else if (colnames(res)[1] == "PFAMID") {
		require(PFAM.db) || stop("required package PFAM.db is not installed")
		pfam.terms <- mget(res$PFAMID, PFAMID, ifnotfound=NA)
		terms <- sprintf("%s (PFAM:%s)", pfam.terms, res[,1])
		terms <- sub("^NA ", "unknown domain ", terms)
	}
	else if (colnames(res)[1] %in% c("GOBPID", "GOMFID", "GOCCID")) {
		# GO term.
		terms <- sprintf("%s (%s)", res$Term, res[,1])
	}
	else {
		# CUSTOM hyperGtest result
		terms <- res[,1]
	}
	# help("HyperGResult-accessors")
	res <- data.frame(
		Term=terms, Count=res$Count, Size=res$Size, 
		Enrichment=signif(res$OddsRatio, digits), 
		P=signif(res$Pvalue, digits), 
		FDR=signif(res$FDR, digits), 
		stringsAsFactors=FALSE 
	)
	
	res
}

###
### from /Users/marcow/src/R/enrichR/dev/bin/hyperGtest.GeneSet.vs.GeneSetCollection.R
###

#' Do a Hypergeometic enrichment test on the genes in a GeneSet, vs all GeneSets in a GeneSetCollection
#'
#' This function takes a GeneSet and GeneSetCollection, both of which use the same Gene ID's eg Entrez Gene,
#' or UniProt, and performs enrichment analysis. It returns a HyperGResult object, from the Category package
#' I'm still not 100% sure this calculates the stat's right - it uses phyper.
#'
#' @param gs a GeneSet
#' @param gsc a GeneSetCollection with >= 1 GeneSet's
#' @param universeSize The number of proteins in the universe. if NULL, this is set to the number of unique
#'	 genes in gsc.
#' @param testName the name of the test. This becomes the colname[1] in the result summary table.
#' @return An instance of HyperGResult. the annotation is forced to org.Mm.eg
#' @author Mark Cowley, 2011-07-01
#' @usage
#' require(GSEABase) || stop("required package GSEABase is not installed")
#' fl <- system.file("extdata", "Broad.xml", package="GSEABase")
#' gsc <- getBroadSets(fl) # GeneSetCollection of 2 sets
#' gs <- gsc[[1]]
#' hgt <- hyperGtest.GeneSet.vs.GeneSetCollection(gs, gsc, NULL)
#' summary(hgt)
#'
hyperGtest.GeneSet.vs.GeneSetCollection <- function(gs, gsc, universeSize=NULL, testName="CUSTOM") {
	!missing(gs) || stop("query GeneSet must be provided")
	!missing(gsc) || stop("GeneSetCollection must be provided")
	!missing(universeSize) || stop("universeSize must be provided")
	require(GSEABase) || stop("required package GSEABase is not installed")
	require(Category) || stop("required package Category is not installed")
	
	if( is.null(universeSize) ) {
		universe <- lapply(gsc, geneIds)
		universe <- sort(unique(unlist(universe)))
		universeSize <- length(universe)
	}
	
	query <- geneIds(gs)
	
	# any(query %in% unlist(lapply(gsc, geneIds))) || stop("There are no ID's in common between query GeneSet, and the GeneSetCollection")
	
	#
	# manually do a HyperGTest
	#
	pvals <- list()
	catToGeneId <- list()
	oddsRatios <- list()
	expectedCounts <- list()
	for(geneset in gsc) {
		idx <- length(pvals) + 1
		
		# q: vector of quantiles representing the number of white balls drawn without replacement from an urn which contains both
		#	 black and white balls.
		# m: the number of white balls in the urn.
		# n: the number of black balls in the urn.
		# k: the number of balls drawn from the urn.
		# 10% bg freq, and we have 2/20
		# eg phyper(q=2,m=100,n=900,k=20, lower.tail=F) # = 0.323
		# eg phyper(q=0,m=100,n=900,k=20, lower.tail=F) # = 0.88
		# eg phyper(q=18,m=100,n=900,k=20, lower.tail=F) # = 3.52e-19
	
		q <- length(intersect(query, geneIds(geneset)))
		m <- length(geneIds(geneset))
		n <- universeSize-m
		k <- length(query)
	
		pvals[[idx]] <- phyper(q-1, m, n, k, lower.tail=FALSE)
		catToGeneId[[idx]] <- intersect(query, geneIds(geneset))
		expectedCounts[[idx]] <- m/universeSize*k
		oddsRatios[[idx]] <- q/expectedCounts[[idx]]
	}
	names(pvals) <- names(catToGeneId) <- names(oddsRatios) <- names(expectedCounts) <- sapply(gsc, setName)
	o <- order(unlist(pvals), decreasing=FALSE)
	pvals <- unlist(pvals)[o]; catToGeneId <- catToGeneId[o]; oddsRatios <- unlist(oddsRatios)[o]; expectedCounts <- unlist(expectedCounts)[o]
	gsc <- gsc[o]
	
	# help("HyperGResult-class")
	hgt <- new("HyperGResult", pvalues=pvals, catToGeneId=lapply(gsc, geneIds), annotation="org.Mm.eg", 
		geneIds=query, testName=testName, pvalueCutoff=1.0, testDirection="over", oddsRatios=oddsRatios,
		expectedCounts=expectedCounts
	)
	
	return(hgt)
}

