#' Rescale a dendrogram based on a data matrix.
#'
#' Recalculate the height of the dendrogram branches based on the within-
#' cluster sum of differences of the data within the numeric matrix mat.
#' All labels of ddg must correspond to rows in mat.
#'
#' @export
#' @param ddg A dendrogram
#' @param mat A numeric matrix with row labels that include all the labels in ddg.
#' @return A dendrogram rescaled so that node heights equal the total within-cluster differences in mat.
#'
rescale <- function (ddg, mat) {
  # Recalculate the height of the dendrogram branches based on the within-
  # cluster sum of differences of the data within the numeric matrix mat.
  # All leaves of ddg must correspond to rows in mat; i.e.
  stopifnot (is (ddg, "dendrogram"));
  stopifnot (is (mat, "matrix") && mode(mat) == "numeric");
  stopifnot (length(intersect(labels(ddg),rownames(mat))) == attr(ddg,"members"));

  # Basic strategy is to evaluate recursively* the total variance across all
  # rows in the current node (for all columns in the matrix) from the branches
  # of the current node, based on the sample means and variances of the node's
  # branches using Result 1 from O'Neill (2014):
  # s_n^2 = 1/(n_1+n_2-1)[(n_1-1)s_1^2 + (n_2-1)s_2^2 + (n_1*n_2)/(n_1+n_2)(mean(x_1)-mean(x_2))^2]
  # We need to accumulate the means separately for each column, but once summarized
  # as variance, all columns can be summed into a single total.  To enable means to be accummulated
  # we actually accumulate sums and divide by the group size as needed.

  # * To prevent large, lop-sided dendrograms from blowing out the R stack, we use an iterative algorithm
  # that places nodes on a to-do stack. Each stack entry has four fields:
  # - The dendrogram for the node concerned
  # - Whether we are traversing down the dendrogram or going back up
  # - The index of this dendrogram in its parent.  It is zero for the root dendrogram node.
  # - A list of the data sums for each child (set when processing child on the way up).
  todo <- list(list(ddg=ddg, down=TRUE, child=0, sums=list()));
  idx <- 1;
  while (idx > 0) {
    ddg <- todo[[idx]]$ddg;
    if (todo[[idx]]$down) {
      if (attr(ddg,"members")==1) {
        # Leaf node.
        attr (ddg, "height") <- 0;  # Just to be sure. Should be unnecessary.
        # Update ddg and set sums in the todo entry for the parent node (if any).
        child <- todo[[idx]]$child;
        if (child > 0) {
          # idx-child is the index of the parent node.
          todo[[idx-child]]$ddg[[child]] <- ddg;
          todo[[idx-child]]$sums[[child]] <- mat[attr(ddg,"label"),];
        }
        # Pop fully processed leaf from the to-do list.
        todo[[idx]] <- NULL;
        idx <- idx - 1;
      } else {
        # On the way down through a non-leaf node.
        # Push each child onto the to-do list on the way down.
        # FYI: Main loop will visit children from right to left.
        todo[[idx]]$down <- FALSE;  # Next time will be on the way up.
        for (m in 1:length(ddg)) {
          todo[[idx+m]] <- list (ddg=ddg[[m]], down=TRUE, child=m, sums=list());
        }
        idx <- idx + length(ddg);
      }
    } else {
      # Compute updated dendrogram node on the way back up.  Note:
      # - All child nodes have been fully processed and updated.
      # - The updated child nodes have been saved in the dendrogram stored in this todo item.
      # - The data sums for each child node have been set in the sums list in this todo item.
      sums <- todo[[idx]]$sums;
      # Initialize the totals for this node from the leftmost child:
      nt <- attr(ddg[[1]], "members");
      st <- attr(ddg[[1]], "height") / nt;
      mt <- sums[[1]];
      # Iteratively add each successive child to the totals, updating:
      # - the number of elements in the total (nt),
      # - the total variance(st), and
      # - the column totals (mt) for the accumulated elements.
      for (rs in 2:length(ddg)) {
        nrs <- attr(ddg[[rs]], "members");
        srs <- attr(ddg[[rs]], "height") / nrs;
        mrs <- sums[[rs]];
        sdiff <- sum((mt/nt-mrs/nrs)**2);
        st <- ((nt-1)*st + (nrs-1)*srs + (nt*nrs)/(nt+nrs)*sdiff)/(nt+nrs-1);
        nt <- nt + nrs;
        mt <- mt + mrs;
      }
      # Update the height of the current dendrogram node:
      stopifnot (attr(ddg,"members") == nt);
      attr(ddg, "height") <- nt*st;
      # Update ddg and set sums in the todo entry for the parent node (if any).
      child <- todo[[idx]]$child;
      if (child > 0) {
        # idx-child is the index of the parent node.
        todo[[idx-child]]$ddg[[child]] <- ddg;
        todo[[idx-child]]$sums[[child]] <- mt;
      }
      # Pop this completed node from the to-do list.
      todo[[idx]] <- NULL;
      idx <- idx - 1;
    }
  }
  # ddg will be the updated initial dendrogram.
  return (ddg);
}
