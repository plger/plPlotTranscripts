#' is.color
#'
#' @param x a vector of values
#'
#' @return a logical vector of the same length as x indicating whether each element of x is a valid color
#' @export
is.color <- function(x){
  sapply(x, function(X){
    tryCatch(is.matrix(col2rgb(X)), error = function(e) FALSE)
  })
}


#' getQualitativePalette
#'
#' Returns a qualitative color palette of the given size. If less than 23 colors are required,
#' the colors are based on Paul Tol's palettes. If more, the `randomcoloR` package is used.
#'
#' @param nbcolors number of colors (from 1 to 22)
#'
#' @return A vector of colors
#'
#' @export
getQualitativePalette <- function(nbcolors){
  nbcolors <- round(nbcolors)
  switch(as.character(nbcolors),
         "1"=c("#4477AA"),
         "2"=c("#4477AA", "#CC6677"),
         "3"=c("#4477AA", "#DDCC77", "#CC6677"),
         "4"=c("#4477AA", "#117733", "#DDCC77", "#CC6677"),
         "5"=c("#332288", "#88CCEE", "#117733", "#DDCC77", "#CC6677"),
         "6"=c("#332288", "#88CCEE", "#117733", "#DDCC77", "#CC6677","#AA4499"),
         "7"=c("#332288", "#88CCEE", "#44AA99", "#117733", "#DDCC77", "#CC6677","#AA4499"),
         "8"=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677","#AA4499"),
         "9"=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677", "#882255", "#AA4499"),
         "10"=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#882255", "#AA4499"),
         "11"=c("#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#882255", "#AA4499"),
         "12"=c("#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#AA4466", "#882255", "#AA4499"),
         "13"=c("#882E72", "#B178A6", "#1965B0", "#5289C7", "#7BAFDE", "#4EB265", "#90C987", "#CAE0AB", "#F7EE55", "#F6C141", "#F1932D", "#E8601C", "#DC050C"),
         "14"=c("#882E72", "#B178A6", "#D6C1DE", "#1965B0", "#5289C7", "#7BAFDE", "#4EB265", "#90C987", "#CAE0AB", "#F7EE55", "#F6C141", "#F1932D", "#E8601C", "#DC050C"),
         "15"=c("#114477", "#4477AA", "#77AADD", "#117755", "#44AA88", "#99CCBB", "#777711", "#AAAA44", "#DDDD77", "#771111", "#AA4444", "#DD7777", "#771144", "#AA4477", "#DD77AA"),
         "16"=c("#114477", "#4477AA", "#77AADD", "#117755", "#44AA88", "#99CCBB", "#777711", "#AAAA44", "#DDDD77", "#771111", "#AA4444", "#DD7777", "#771144", "#AA4477", "#DD77AA", "black"),
         "17"=c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788"),
         "18"=c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788"),
         "19"=c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788", "black"),
         "20"= c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788"),
         "21"= c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788"),
         "22"= c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788", "black"),
         stop("Max 22 colors.")
  )
}


getLinearColMap <- function(x,colors=NULL,nbsteps=100){
  if(is.null(colors)){
    library(viridis)
    cols <- viridis(nbsteps)
  }else{
    cols <- colorRampPalette(colors)(nbsteps)
  }
  x <- (x-min(x))/(max(x)-min(x))
  cols[round(x*(nbsteps-1))+1]
}

.parse2Colors <- function(x,...){
  if(is.numeric(x)) return(getLinearColMap(x,...))
  if(all(is.color(x))) return(x)
  return( getQualitativePalette(length(unique(x)))[as.numeric(as.factor(x))] )
}
