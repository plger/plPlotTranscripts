#' plPlotTranscripts
#'
#' Plots the exon-intron structure of a set of transcripts. See `vignette("plPlotTranscripts)`
#'
#' @param gr a GRanges object imported from a gff file (as done by `rtracklayer::import.gff`).
#' @param title the plot title (default none).
#' @param sortBy fields to sort the transcripts (default `c('transcript_type','structure')`).
#' If `NULL`, the transcripts are not sorted and the order is taken as is. In addition to the columns of `gr`'s
#' metadata, the following elements can be used: 'start', 'end', 'strand', and 'structure' (structure similarity).
#' @param txData the columns of `gr`'s metadata that should be displayed as hoverinfo.
#' In addition to the columns of `gr`, coordinates is also available (start-stop and strand).
#' Default transcript id, type, gene name and coordinates.
#' @param label.field the column of `gr` to use as transcript label (default `transcript_id`)
#' @param annotations an optional data.frame of annotations to add to the plot. The data.frame should contain the columns
#' `x` (coordinate, numeric), `y` (either an integer specifying the transcript, or a character vector specifying transcript ids),
#' and `text`
#' @param maxIntronSize the maximum size (in nt) of introns. Stretches of intron not overlapping any exon
#'  will be scaled down to this size if needed. Default 50. If set to NULL, actual genomic coordinates are used.
#' @param plotExonsSeparately logical; whether to plot each exon separately (allowing different coloring and hoverinfo)
#'  rather than whole transcripts. Default FALSE. the space between transcripts (default 10)
#' @param space the space between transcripts (default 10)
#' @param cdsW the line width of coding regions (default 15)
#' @param utrW the line width of untranslated regions (default 10)
#' @param intronW the line width of introns (default 2)
#' @param colorBy either the name of a metadata column in `gr` according to which the transcript colors should be
#' determined, or a named vector of values with transcript ids as names. Default (NULL) uses fixed colors.
#' @param intronColor the color of introns (default darkgrey). Ignored if `colorBy` is given.
#' @param exonColor the color of untranslated exonic regions (default darkgrey). Ignored if `colorBy` is given.
#' @param cdsColor the color of coding regions (default black). Ignored if `colorBy` is given.
#' @param ... further arguments passed directly to `add_annotations()` (ignored if `annotations=NULL`)
#'
#' @return a plotly widget.
#' @author Pierre-Luc Germain
#'
#' @export
plPlotTranscripts <- function( gr,
                               title='',
                               sortBy=c("transcript_type","structure"),
                               txData=c("transcript_id","transcript_type","gene_name","coordinates"),
                               label.field="transcript_id",
                               annotations=NULL,
                               maxIntronSize=50,
                               plotExonsSeparately=FALSE,
                               space=10,
                               cdsW=15,
                               utrW=10,
                               intronW=2,
                               intronColor="darkgrey",
                               colorBy=NULL,
                               exonColor="darkgrey",
                               cdsColor="black",
                               ... ){

  library(GenomicRanges)
  library(plotly)

  if(!is.null(annotations)){
    if(!is.data.frame(annotations) | !all(c("x","y","text") %in% colnames(annotations))){
      warning("`annotations` does not match requirements (should be a data.frame with at least the columns `x`, `y`, and `text`) and will be ignored. ")
      annotations <- NULL
    }
  }
  
  cn <- colnames(gr@elementMetadata)
  if( !("transcript_type" %in% cn) && "transcript_biotype" %in% cn ){
    gr$transcript_type <- gr$transcript_biotype
  }
  
  # we first identifiy the number of transcripts and eventually sort them
  gr <- gr[which(gr$type %in% c("transcript","exon","CDS")),]
  txs <- gr[which(gr$type=="transcript"),]
  if(!is.null(sortBy) && length(sortBy)>0){
    to <- .sortTranscripts(gr,sortBy)
    txs <- txs[to,,drop=F]
  }else{
    sortBy <- NULL
  }
  ntx <- length(txs)

  # color assignment
  cols <- NULL
  if(!is.null(colorBy)){
    if(length(colorBy)>1 || !(colorBy %in% names(txs@elementMetadata))){
      if(all(txs$transcript_id %in% names(colorBy))){
        cols <- .parse2Colors(colorBy[txs$transcript_id])
        names(cols) <- txs$transcript_id
        if(plotExonsSeparately){
          gr$color <- cols[gr$transcript_id]
          colorBy <- "color"
        }
      }else{
        warning("`colorBy` is neither a column of `gr` nor a named vector of transcript assignment, and will be ignored.")
        colorBy <- NULL
      }
    }else{
      if(plotExonsSeparately){
        gr$color <- as.character(.parse2Colors(gr@elementMetadata[,colorBy]))
        colorBy <- "color"
      }else{
        cols <- .parse2Colors(txs@elementMetadata[,colorBy])
        names(cols) <- txs$transcript_id
      }
    }
  }

  # preparing the hover-text
  if(plotExonsSeparately){
    # exons are plotted sperately, so we need to set hover text (including eventually coordinates) separately
    if("coordinates" %in% txData) gr$coordinates <- paste0(as.character(seqnames(txs)),":",start(gr),"-",end(gr)," ",strand(gr))
    gr$info <- apply(as.data.frame(gr@elementMetadata[,txData]),1,collapse="\n",FUN=paste)
    txs <- gr[which(gr$type=="transcript"),]
    if(!is.null(sortBy)) txs <- txs[to,,drop=F]
  }else{
    # while transcripts are plotted
    if("coordinates" %in% txData) txs$coordinates <- paste0(as.character(seqnames(txs)),":",start(txs),"-",end(txs)," ",strand(txs))
  }
  txtext <- apply(as.data.frame(txs@elementMetadata[,txData]),1,collapse="\n",FUN=paste)
  names(txtext) <- txs$transcript_id

  # y-coordinates are the index of the transcript times `space`
  y <- 1:ntx*space
  names(y) <- txs$transcript_id

  # defining an empty axis; will be used for y-axis, and if minIntronSize!=NULL also for x-axis
  noaxis <- list(
    zeroline = FALSE,
    showline = FALSE,
    showticklabels = FALSE,
    showgrid = FALSE)

  exons <- gr[which(gr$type=="exon"),]
  cds <- gr[which(gr$type=="CDS"),]


  if(!is.null(maxIntronSize)){
    # we scale everything to have small introns
    # first, we determine intronic regions that are not overlapping any exon
    exons <- sort(exons)
    introns <- cbind(end(exons)[-length(exons)],start(exons)[-1])
    introns <- introns[!duplicated(introns),,drop=F]
    introns <- introns[which(apply(introns,1,mis=maxIntronSize,FUN=function(x,mis){ (x[[2]]-x[[1]]) > mis})),,drop=F]
    # introns now contains the sections that will be cut
    reductions <- introns[,2]-introns[,1]-maxIntronSize
    for(i in 1:nrow(introns)){
      # for each section to be cut, we substract its length from
      # all coordinates that are >= to it
      w <- which(start(exons)>=introns[i,2])
      start(exons)[w] <- start(exons)[w]-reductions[i]
      end(exons)[w] <- end(exons)[w]-reductions[i]
      w <- which(start(cds)>=introns[i,2])
      start(cds)[w] <- start(cds)[w]-reductions[i]
      end(cds)[w] <- end(cds)[w]-reductions[i]
      w <- which(start(txs)>=introns[i,2])
      start(txs)[w] <- start(txs)[w]-reductions[i]
      w <- which(end(txs)>=introns[i,2])
      end(txs)[w] <- end(txs)[w]-reductions[i]
      if(!is.null(annotations)){
        w <- which(annotations$x>=introns[i,2])
        annotations$x[w] <- annotations$x[w]-reductions[i]
      }
      introns <- introns-reductions[i]
    }
    xaxis <- noaxis
  }else{
    xaxis <- list( title="Genomic position")
  }

  # plot initiation
  p <- plot_ly()

  # plot transcript lines (straight lines, without introns-exons)
  for(i in 1:ntx){
    p <- p %>% add_trace( type="scatter", mode = 'lines',
                          x=c(start(txs[i,]),end(txs[i,])),
                          y = y[i],
                          text=txtext[i],
                          hoverinfo='text',
                          name = txs$transcript_id[i],
                          line=list( color=ifelse(is.null(colorBy),intronColor,colorBy[txs$transcript_id[i]]),
                                     width=intronW))
  }

  # draw exons and CDS
  if(plotExonsSeparately){
    p <- .plPlotTranscripts.drawBlocks(p, exons, y, colorBy, "info", utrW, exonColor)
    if(length(cds)>0) p <- .plPlotTranscripts.drawBlocks(p, cds, y, colorBy, "info", cdsW, cdsColor)
  }else{
    p <- .plPlotTranscripts.drawTranscripts(p, exons, y, cols, txtext, utrW, exonColor)
    if(length(cds)>0) p <- .plPlotTranscripts.drawTranscripts(p, cds, y, cols, txtext, cdsW, cdsColor)
  }

  # transcript labels
  if(!is.null(label.field) && label.field %in% names(txs@elementMetadata)){
    p <- p %>% add_annotations( x=start(txs),
                                y=y,
                                text=paste0(txs@elementMetadata[,label.field],"  "),
                                xref = "x",
                                yref = "y",
                                xanchor='right',
                                showarrow=F)
  }

  # draw annotations
  if(!is.null(annotations)){
    if(is.factor(annotations$y)) annotations$y <- as.character(annotations$y)
    if(is.character(annotations$y)){
      names(y) <- txs@elementMetadata$transcript_id
      annotations$y <- y[annotations$y]
    }else{
      annotations$y <- annotations$y*space
    }
    p <- add_annotations(p, data=annotations, x=~x, y=~y, text=~text, ...)
  }
  p %>% layout(title=title, xaxis=xaxis, yaxis=noaxis, showlegend = FALSE)
}

.plPlotTranscripts.drawBlocks <- function(p, exons, y, colorField, textField, width, defaultColor){
  for(i in 1:length(exons)){
    ey <- y[exons$transcript_id[i]]
    p <- p %>% add_trace(
      x = c(start(exons)[i],end(exons)[i]),
      y = rep(ey,2),
      type = 'scatter',
      mode = 'lines',
      hoverinfo = 'text',
      text = exons@elementMetadata[[textField]][i],
      line = list( color=ifelse(is.null(colorField),defaultColor,exons@elementMetadata[i,colorField]),
                   width=width)
    )
  }
  p
}

.plPlotTranscripts.drawTranscripts <- function(p, exons, y, colors, info, width, defaultColor){
  exons <- split(exons,exons$transcript_id)
  for(i in 1:length(exons)){
    te <- sort(exons[[i]])
    ey <- y[te$transcript_id[1]]
    x <- as.numeric(t(cbind(start(te),end(te),rep(NA,length(te)))))
    p <- p %>% add_trace( type="scatter", mode = 'lines',
                          x = x,
                          y = rep(ey,length(x)),
                          text=info[te$transcript_id[1]],
                          hoverinfo='text',
                          line=list( color=ifelse(is.null(colors),defaultColor,colors[te$transcript_id[1]]),
                                     width=width)
    )
  }
  return(p)
}


.sortTranscripts <- function(gr, sortBy){
  nbtx <- length(unique(gr@elementMetadata[which(gr$type=="transcript"),"transcript_id"]))
  if(nbtx<=2) return(1:nbtx)
  if('structure' %in% sortBy){
    library(seriation)
    upos <- unique(c(start(gr),end(gr)))
    bm <- t(sapply(split(c(start(gr),end(gr)),rep(gr$transcript_id,2)), upos=upos, FUN=function(x, upos){ upos %in% x }))
    gr <- gr[which(gr$type=="transcript"),]
    bm <- bm[gr$transcript_id,]
    gr$structure <- get_order(seriate(dist(bm), method = "MDS_angle"))
  }else{
    gr <- gr[which(gr$type=="transcript"),]
  }
  if("start" %in% sortBy) gr$start <- start(gr)
  if("end" %in% sortBy) gr$end <- end(gr)
  if("strand" %in% sortBy) gr$strand <- strand(gr)
  tx <- as.data.frame(gr@elementMetadata, row.names = 1:length(gr))
  for(f in rev(sortBy)) tx <- tx[order(tx[[f]]),]
  return(as.numeric(row.names(tx)))
}
