import_supercluster <- function(file, geneIdType=ENSEMBLIdentifier("org.Mm.eg.db")) {
	txt <- readLines(file)
	
	# if( is.null(geneIdType) ) {
	# 	geneIdType <- ENSEMBLIdentifier("org.Mm.eg.db")
	# }
	
	symbols <- list()
	geneset.names <- list()
	
	i <- 1
	geneset.idx <- 0
	done <- FALSE
	while (!done) {
		if( i == 1 || txt[i] == "<EndOfFile>::::::::::::::" ) {
			# then new record has started
			geneset.idx <- geneset.idx + 1
			geneset.names[[geneset.idx]] <- sub(".csv", "", txt[i+1])
			numGenes <- as.numeric(sub(",.*", "", txt[i+4]))
			symbols[[geneset.idx]] <- txt[(i+5):(i+4+numGenes)]
			i <- i+4+numGenes+1
		}
		else if ( txt[i] == "<EndOfFile>" ) {
			done = TRUE
		}
		else {
			i <- i + 1
		}
	}
	names(symbols) <- geneset.names
	
	# trim the trailing version numbers from the ID's
	for(i in seq(along=symbols)) {
		symbols[[i]] <- sub("\\.[0-9]+", "", symbols[[i]])
	}
	
	#return(symbols)
	
	genesetlist <- list()
	# convert into GeneSet objects
	for(i in seq(along=symbols)) {
		genesetlist[[i]] <- GeneSet(symbols[[i]], geneIdType=geneIdType, setName=geneset.names[[i]])
		# library(org.Mm.eg.db)
		# genesetlist[[i]] <-mget(symbols[[i]], org.Mm.egENSEMBL2EG, ifnotfound=NA)	
	}
 	
 	# return(genesetlist)

	# convert into GeneSetCollection
	res <- GeneSetCollection(genesetlist)
	
	res
	}
