#' Generate (perhaps approximate) dendrogram for heirarchical clustering of the rows of mat.
#'
#' mat must be a numeric matrix with row labels.
#' mat will be column normalized before clustering. So no column can have constant values.
#'
#' Approximate methods will be employed if:
#'
#' - the number of columns exceeds pcLimit (default 512):
#'   Clustering will be performed using the first pcSize (default 128) principal components.
#'   If the number of rows exceeds pcMaxRows (default 2048), the principal components will be
#'   approximated using pcMaxRows randomly selected rows.
#'
#' - the number of rows exceeds hclustMax (default 25000), the rows will partitioned using
#'   kmeans 2 as many times as required until each partition has at most hclustMax rows.
#'
#'
#' Rclusterpp::Rclusterpp.hclust will be used for (subsets of) the matrix with at most hclustMax rows.
#' The parameters method, distance, and p are passed to it without modification.
#' @export
#' @param mat A numeric matrix with row labels
#' @param method Agglomeration method passed to hclust
#' @param distance Distance metric passed to hclust
#' @param p Minkowski distance parameter passed to hclust
#' @param pcLimit Maximum number of columns before using principal components
#' @param pcSize Number of principal components to use
#' @param pcMaxRows Maximum number of rows to use for principle components
#' @param hclustMax Maximum number to rows to hclust at a time
#' @return A dendrogram for the rows of mat
#'
#' @importFrom methods is
#' @importFrom stats as.dendrogram kmeans prcomp sd
aclust <- function (mat,
		   method = "ward",
		   distance = "euclidean",
		   p = 2,
		   pcLimit = 512,
		   pcSize = 128,
		   pcMaxRows = 2048,
		   hclustMax = 25000
		   ) {
  # Validate input matrix.
  stopifnot (is(mat,"matrix") && mode(mat) == "numeric");
  columnStandardDeviations <- apply (mat, 2, function(x)sd(x,na.rm=TRUE));
  stopifnot (all(columnStandardDeviations > 0));
  # Column normalize input matrix.
  cc <- apply (mat, 2, function(x)mean(x,na.rm=TRUE));
  mat <- t((t(mat) - cc)/columnStandardDeviations);

  # Generate perhaps approximate matrix to cluster.
  if (ncol(mat) <= pcLimit) {
    pcMat <- mat;
  } else {
    if (nrow(mat) <= pcMaxRows) {
      pc <- prcomp (mat, center = FALSE, rank = pcSize);
      pcMat <- pc$x;
    }
    else {
      rr <- sample.int(nrow(mat), pcMaxRows);
      pc <- prcomp (mat[rr,], center = FALSE, rank = pcSize);
      pcMat <- mat %*% pc$rotation;
    }
  }

  # "Recursively" subdivide matrix into subsets with at most hclustMax rows in each.
  # Perform hierarchical clustering on each fragment, then reassemble into an overall
  # clustering.
  #
  # To avoid exhausting R's stack on very large inputs, we use an iterative algorithm
  # that implements the tree traversal using a todo list.
  #
  # On the way down, we subdivide the rows into limited size clusters using kmeans 2.
  # At the leaves we perform hierarchical clustering.
  # On the way up, we recombine the dendrograms from the node's children.
  #
  todo <- list(list(rows=rownames(pcMat), down=TRUE, child=0, ddgs=list()));
  idx <- 1;
  while (idx > 0) {
    rows <- todo[[idx]]$rows;
    if (todo[[idx]]$down) {
      if (length(rows) > hclustMax) {
	# Too big for hclust. Split once using kmeans
        todo[[idx]]$down <- FALSE;  # Next time will be on the way up.
	km <- kmeans (pcMat[rows,], 2);
	todo[[idx+1]] <- list (rows=names(km$cluster)[which(km$cluster==1)], down=TRUE, child=1, ddgs=list());
	todo[[idx+2]] <- list (rows=names(km$cluster)[which(km$cluster==2)], down=TRUE, child=2, ddgs=list());
        idx <- idx + 2;
      } else {
	# Small enough for hclust.
	if (length(rows) == 1) {
	  # Just in case we get a really bad split from kmeans
	  ddg <- which (rownames(pcMat) == rows);
	  attr(ddg, "members") <- as.integer(1);
	  attr(ddg, "height") <- 0;
	  attr(ddg, "label") <- rows;
	  attr(ddg, "leaf") <- TRUE;
	  class(ddg) <- "dendrogram";
	} else {
	  hc <- Rclusterpp::Rclusterpp.hclust(pcMat[rows,], method = method, distance = distance, p = p);
	  ddg <- fixvalues (as.dendrogram(hc), rownames(pcMat));
	}
        # Update ddg in the todo entry for the parent node (if any).
        child <- todo[[idx]]$child;
        if (child > 0) {
          # idx-child is the index of the parent node.
          todo[[idx-child]]$ddgs[[child]] <- ddg;
        }
        # Pop fully processed leaf from the to-do list.
        todo[[idx]] <- NULL;
        idx <- idx - 1;
      }
    } else {
      # By the time a node is processed on the way up, its children will have saved their
      # dendrograms into the ddgs list in the node's todo entry.
      # The dendrogram value for a node is just a list of its children.  Convenient!
      ddg <- todo[[idx]]$ddgs;
      # Set dendrogram specific attributes and class.
      # Sum total number of members.  Determine maximum height of child nodes.
      nt <- attr (ddg[[1]], "members");
      h <- attr (ddg[[1]], "height");
      for (ii in 2:length(ddg)) {
	      nt <- nt + attr(ddg[[ii]], "members");
	      h <- max(h,attr(ddg[[ii]],"height"));
      }
      attr(ddg, "members") <- nt;
      attr(ddg, "height") <- h + 1;  # Will be set to proper value later.
      class(ddg) <- "dendrogram";
      # Update ddg in the todo entry for the parent node (if any).
      child <- todo[[idx]]$child;
      if (child > 0) {
	# idx-child is the index of the parent node.
	todo[[idx-child]]$ddgs[[child]] <- ddg;
      }
      # Pop fully processed node from the to-do list.
      todo[[idx]] <- NULL;
      idx <- idx - 1;
    }
  }
  # The final dendrogram will be left in ddg.

  # Scale dendrogram heights using original (column-centered) data before returning.
  rescale (ddg, mat)
}

# Fix the values of the leaves of ddg so that they equal the index of each leaf label in labels.
# Needed because hclust may only process a subset of rows, in which case the leaves will be assigned
# the wrong values.
#
# Use an iterative algorithm using a to-do list to traverse the tree so that we don't run out of
# stack space.
fixvalues <- function (ddg, labels) {
  todo <- list(list(ddg=ddg, down=TRUE, child=0));
  idx <- 1;
  while (idx > 0) {
    ddg <- todo[[idx]]$ddg;
    if (todo[[idx]]$down) {
      if (attr(ddg,"members")==1) {
        # Leaf node.
        # Set the numeric leaf value to the index of its label.
	# FYI: ddg-ddg zeroes the value, but keeps the dendrogram class and properties.
        ddg <- ddg-ddg+which(attr(ddg,"label") == labels);
        # Update ddg in the todo entry for the parent node (if any).
        child <- todo[[idx]]$child;
        if (child > 0) {
          # idx-child is the index of the parent node.
          todo[[idx-child]]$ddg[[child]] <- ddg;
        }
        # Pop fully processed leaf from the to-do list.
        todo[[idx]] <- NULL;
        idx <- idx - 1;
      } else {
        # On the way down through a non-leaf node:
        #   Push each child onto the to-do list.
        # FYI: Main loop will visit children from right to left.
        todo[[idx]]$down <- FALSE;  # Next time will be on the way up.
        for (m in 1:length(ddg)) {
          todo[[idx+m]] <- list (ddg=ddg[[m]], down=TRUE, child=m);
        }
        idx <- idx + length(ddg);
      }
    } else {
      # Process updated dendrogram node on the way back up.
      # Nothing to do but update ddg in the todo entry for the parent node (if any).
      child <- todo[[idx]]$child;
      if (child > 0) {
        # idx-child is the index of the parent node.
        todo[[idx-child]]$ddg[[child]] <- ddg;
      }
      # Pop this completed node from the to-do list.
      todo[[idx]] <- NULL;
      idx <- idx - 1;
    }
  }
  # ddg will be the updated root of the dendrogram.
  ddg
}
