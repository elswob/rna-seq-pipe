require(GSEABase)

setMethod("collectionType", 
          signature=signature(object="GeneSetCollection"),
          function(object) {
			unlist(lapply(as.list(object), collectionType))
		  }
)
# collectionType(modules)


setMethod("contributor", 
          signature=signature(object="GeneSetCollection"),
          function(object) {
			unlist(lapply(as.list(object), contributor))
		  }
)
# contributor(modules)


setMethod("creationDate", 
          signature=signature(object="GeneSetCollection"),
          function(object) {
			unlist(lapply(as.list(object), creationDate))
		  }
)
# creationDate(modules)



setMethod("description", 
          signature=signature(object="GeneSetCollection"),
          function(object) {
			unlist(lapply(as.list(object), description))
		  }
)
# description(modules)



setMethod("geneIds", 
          signature=signature(object="GeneSetCollection"),
          function(object) {
			res <- lapply(as.list(object), geneIds)
			names(res) <- setName(object)
			res
		  }
)
# geneIds(modules)



setMethod("longDescription", 
          signature=signature(object="GeneSetCollection"),
          function(object) {
			unlist(lapply(as.list(object), longDescription))
		  }
)
# longDescription(modules)


setMethod("organism", 
          signature=signature(object="GeneSetCollection"),
          function(object) {
			unlist(lapply(as.list(object), organism))
		  }
)
# organism(modules)

setReplaceMethod("collectionType",
		signature=signature(object="GeneSetCollection",
							value="CollectionType"),
		function(object, value) {
			res <- GeneSetCollectionapply(object, function(gs) {collectionType(gs) <- value; gs})
			res
		}
)



setMethod("pubMedIds", 
          signature=signature(object="GeneSetCollection"),
          function(object) {
			lapply(as.list(object), pubMedIds)
		  }
)
# pubMedIds(modules)


setMethod("setIdentifier", 
          signature=signature(object="GeneSetCollection"),
          function(object) {
			unlist(lapply(as.list(object), setIdentifier))
		  }
)
# setIdentifier(modules)


setMethod("geneIdType", 
          signature=signature(object="GeneSetCollection"),
          function(object) {
			unlist(lapply(as.list(object), geneIdType))
		  }
)
# geneIdType(modules)

setReplaceMethod("geneIdType",
		signature=signature(object="GeneSetCollection",
							value="GeneIdentifierType"),
		function(object, value) {
			res <- GeneSetCollectionapply(object, function(gs) {geneIdType(gs) <- value; gs})
			res
		}
)

setMethod("setName", 
          signature=signature(object="GeneSetCollection"),
          function(object) {
			res <- GeneSetCollectionapply(object, setName)
			res <- unname(unlist(res))
			res
			# unlist(lapply(as.list(object), setName))
		  }
)
# setName(modules)

setReplaceMethod("setName", 
          signature=signature(object="GeneSetCollection", value="character"),
          function(object, value) {
			length(value) == length(object) || stop("lengths of object and value must match")
			for(i in seq(along=object)) {
				setName(object[[i]]) <- value[i]
			}
			object
		  }
)


setMethod("setVersion", 
          signature=signature(object="GeneSetCollection"),
          function(object) {
			unlist(lapply(as.list(object), setVersion))
		  }
)
# setVersion(modules)



setMethod("urls", 
          signature=signature(object="GeneSetCollection"),
          function(object) {
			unlist(lapply(as.list(object), urls))
		  }
)
# urls(modules)



#####
# New functions not already specified in GeneSet

setMethod("annotation", 
          signature=signature(object="GeneSet"),
          function(object) {
			  annotation( geneIdType(object) )
		  }
)
setReplaceMethod("annotation", 
          signature=signature(object="GeneSet", value="character"),
		  function(object, value) {
			  class <- paste( geneIdType(object)@type, "Identifier", sep="")
			  geneIdType(res) <- new(class, annotation=value)
			  res
		  }
)


setMethod("annotation", 
          signature=signature(object="GeneSetCollection"),
          function(object) {
			  unlist(lapply(as.list(object), function(gs) annotation(geneIdType(gs))))
		  }
)
setReplaceMethod("annotation", 
          signature=signature(object="GeneSetCollection", value="character"),
		  function(object, value) {
			  res <- GeneSetCollectionapply(object, function(gs) {annotation(gs) <- value; gs})
			  res
		  }
)




#' Apply a Function over a GeneSetCollection
#' Return a GeneSetCollection of the same length as \code{X}, 
#' each element of which is the result of applying \code{FUN} to the corresponding 
#' element of \code{X}.
#'
#' \code{FUN} is found by a call to \code{\link{match.fun}} and typically is specified
#' as a function or a symbol (e.g. a backquoted name) or a character
#' string specifying a function to be searched for from the
#' environment of the call to \code{\link{GeneSetCollectionapply}}.
#' 
#' Function \code{FUN} must be able to accept as input a \code{GeneSet} object.
#' 
#' @param X a GeneSetCollection.
#' @param FUN the function to be applied to each element of \code{X}. see
#'        \code{Details}.
#' @param \dots arguments passed to \code{FUN}
#' @return a \code{list}, or a \code{GeneSetCollection} object, depending on the return value of \code{FUN}
#'   If FUN ends up destroying all GeneSet's (set's them to \code{NA}), then \code{NA} is returned.
#' @author Mark Cowley, 2011-08-11
#' @export
GeneSetCollectionapply <- function(X, FUN, ...) {
	is(X, "GeneSetCollection") || stop("x must be a GeneSetCollection")

	res <- lapply(as.list(X), FUN, ...)
	names(res) <- sapply(as.list(X), setName)
	suppressWarnings(res <- res[!is.na(res)])
	if( length(res) == 0 ) return(NA)
	# convert a list of GeneSet's back into a GeneSetCollection
	if( is.list(res) && is(res[[1]], "GeneSet") ) {
		res <- GeneSetCollection(res, geneIdType=geneIdType(X), collectionType=collectionType(X))
	}

	res
}
# CHANGELOG
# 2011-11-07:
# - added names to res straight away
# 

# setMethod("as.list",
# 	signature=signature("GeneSetCollection"),
# 	function(x, ...) {
# 		res <- list()
# 		for(i in seq(along=x)) {
# 			res[[i]] <- x[[i]]
# 		}
# 		names(res) <- setName(x)
# 		
# 		res
# 	}
# )
# Error: evaluation nested too deeply: infinite recursion / options(expressions=)?


#' Apply a function to a GeneSetCollection over a snow cluster
#'
#' \code{FUN} is found by a call to \code{\link{match.fun}} and typically is specified
#' as a function or a symbol (e.g. a backquoted name) or a character
#' string specifying a function to be searched for from the
#' environment of the call to \code{\link{GeneSetCollectionClusterApply}}.
#' 
#' Function \code{FUN} must be able to accept as input a \code{GeneSet} object.
#'  
#' @param cl a snow cluster object
#' @param X a \code{GeneSetCollection} object
#' @param FUN the function to be applied to each element of \code{X}. see
#'        \code{Details}.
#' @param \dots arguments passed to \code{FUN}
#' @return a \code{list}, or a \code{GeneSetCollection} object, depending on the return value of \code{FUN}
#' @return
#' @author Mark Cowley, 2011-08-11
#' @export
GeneSetCollectionClusterApply <- function(cl, X, FUN, ...) {
	is(X, "GeneSetCollection") || stop("X must be a GeneSetCollection")
	
	X2 <- as.list(X); names(X2) <- sapply(X2, setName)
	res <- clusterApply(cl, X2, FUN, ...)
	if( is.list(res) && is(res[[1]], "GeneSet") ) {
		res <- res[match(sapply(res, setName), sapply(X, setName))]
		res <- GeneSetCollection(res, geneIdType=geneIdType(X), collectionType=collectionType(X))
	}
	
	res
}


GeneSetCollectionClusterApplyLB <- function(cl, X, FUN, ...) {
	is(X, "GeneSetCollection") || stop("X must be a GeneSetCollection")
	
	X2 <- as.list(X); names(X2) <- sapply(X2, setName)
	res <- clusterApplyLB(cl, X2, FUN, ...)
	if( is.list(res) && is(res[[1]], "GeneSet") ) {
		res <- res[match(sapply(res, setName), sapply(X, setName))]
		res <- GeneSetCollection(res, geneIdType=geneIdType(X), collectionType=collectionType(X))
	}
	
	res
}

#' Convert a GeneSetCollection to an environment
#' Convert a GeneSetCollection (geneset to gene mapping) to an environment mapping
#' from genes to genesets (note the reversal of direction in here.)
#'
#' TODO: can this be made generic using as?
#' 
#' @param x a GeneSetCollection object
#' @return an environment, where genes are the search keys, and gensets are the values.
#' @author Mark Cowley, 2011-08-12
#' @examples
#' g1 <- c("1", "10", "100", "1000", "10000", "100008586", "10001", "10002", "10003", "100037264")
#' gs1 <- GeneSet(g1, geneIdType=EntrezIdentifier("org.Mm.eg.db"), setName="gs1")
#' g2 <- c("100037265", "10004", "100049542", "100049716", "10005", "10006", "10007", "10008", "10009", "100093630")
#' gs2 <- GeneSet(g2, geneIdType=EntrezIdentifier("org.Mm.eg.db"), setName="gs2")
#' gsc <- GeneSetCollection(list(gs1, gs2))
#' env <- GeneSetCollection2environment(gsc)
#' env
#' names(as.list(env))
#' @export
GeneSetCollection2environment <- function(x) {
	is(x, "GeneSetCollection") || stop("x must be a GeneSetCollection")
	require(Biobase) || stop("required package 'Biobase' is not installed")
	l <- as.list(x)
	names(l) <- sapply(x, setName)
	geneset2gene <- lapply(l, geneIds)
	gene2geneset <- reverseSplit(geneset2gene)
	gene2geneset.env <- list2env(gene2geneset)
	
	return(gene2geneset.env)
}




# #' How many genes in a GeneSet or each element of a GeneSetCollection
# #'
# #' @param x a GeneSet or GeneSetCollection object
# #' @param na.rm logical: if \code{TRUE}, then only count the non-\code{NA} elements
# #' @return if \code{x} is a \code{GeneSet}, a \code{numeric(1)} indicating 
# #' the number of genes in the \code{GeneSet}.\cr
# #' if \code{x} is a \code{GeneSetCollection}, a \code{numeric(length(x))} indicating the number of 
# #'   genes in each GeneSet within \code{x}
# #' @author Mark Cowley, 2011-08-12
# #' @export
# numGenes <- function(object, na.rm=TRUE) {
# 	if(is(object, "GeneSetCollection")) sapply(object, numGenes)
# 	else if(is(object, "GeneSet")) {
# 		gi <- geneIds(object)
# 		sum( !is.na(gi) )
# 	}
# 	else stop("object is not a GeneSet or GeneSetCollection")
# }

setGeneric("numGenes",
  function(object, na.rm=TRUE)
    standardGeneric("numGenes")
)

setMethod("numGenes", 
          signature=signature(object="GeneSetCollection"),
          function(object, na.rm=TRUE) {
			sapply(object, numGenes, na.rm=na.rm)
		  }
)

setMethod("numGenes", 
          signature=signature(object="GeneSet"),
          function(object, na.rm=TRUE) {
				gi <- geneIds(object)
				if( na.rm ) sum( !is.na(gi) ) else length(gi)
		  }
)
