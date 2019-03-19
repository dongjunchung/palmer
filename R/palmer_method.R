
#' @import pheatmap
#' @exportMethod predict plot

#' @title Summary of palmer results
#'
#' @description Summary of palmer results
#'
#' @aliases show.palmer
#'
#' @param object output of palmer function

setMethod(
  f="show",
  signature="palmer",
  definition=function( object ) {

    init <- object@init

    cat( "Summary: palmer (class: palmer)\n" )
    cat( "--------------------------------------------------\n" )
    cat( "Model settings \n")
    cat( "Number of samples to be analyzed: ", init$nsample, "\n", sep="" )
    cat( "Number of variables to be analyzed: ", init$nvariable, "\n", sep="" )
    cat( "Number of gene clusters: ", init$K, "\n", sep="" )
    cat( "Number of GO term clusters: ", init$L, "\n", sep="" )
    cat( "Number of bootstrapping: ", init$B, "\n", sep="" )
    cat( "--------------------------------------------------\n" )
  }
)




#' @title palmer result
#'
#' @description Return palmer result
#'
#' @aliases predict.palmer
#'#'
#' @param object output of palmer function


setMethod(
  f="predict",
  signature="palmer",
  definition=function( object ) {

    # extract objects

    result <- object@result
    gene <- cbind(result$Genecluster, result$Geneprob)
    go <- cbind(result$GOcluster,result$GOprob)

    return(list(
      Gene = gene,
      GO = go
    ))
  }
)



#' @title Plotting function
#'
#' @description Heatmap with clustered data
#'
#' @aliases plot.palmer
#'
#' @param x output of palmer function

setMethod(
  f="plot",
  signature=c("palmer","missing"),
  definition=function( x, y,... ) {

    data <- x@data
    result <- x@result
    KC <- result$Genecluster
    LC <- result$GOcluster
    geneset <- result$Geneset

    KC1 <- sort(KC); LC1 <- sort(LC)

    anno_col <- data.frame(GO=factor(LC1))
    anno_row=data.frame(Cluster=factor(KC1),Gene=factor(geneset[match(names(KC1),names(geneset))]))
    rownames(anno_row) <- names(KC1)

    row.index <- match(rownames(anno_row),rownames(data))
    col.index <- match(rownames(anno_col),colnames(data))

    color <- c("red","green","blue")
    gene <- color[unique(geneset)]; names(gene) <- unique(geneset)
    cluster <- color[unique(KC1)]; names(cluster) <- unique(KC1)
    go <- color[unique(LC1)]; names(go) <- unique(LC1)
    colcolors <- list(Gene=gene, Cluster=cluster, GO=go)

    heat <- pheatmap(data[row.index,col.index],color=c("white", "black"),cluster_row=FALSE, cluster_col=F,
                     annotation_col=anno_col, annotation_row=anno_row,
                     legend=T, annotation_colors=colcolors,annotation_legend=T,
                     show_rownames = T, show_colnames = T,...

    )

  }
)
