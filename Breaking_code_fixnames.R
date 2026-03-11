# Exemple of function to correct the name in a seurat obje

fixnames <- function(x){
  if(class(x)=='character'){
    x[x=='488'] <- 'D4_repro'
    x[x=='509'] <- 'D0_repro'
    x[x=='489'] <- 'D7_repro'
    x[x=='547pos'] <- 'D42_transfo'
    x[x=='508'] <- 'poumon_total'
    x[x=='546pos'] <- 'D72_transfo'
    x[x=='467'] <- 'D16_transfo'
    x[x=='510'] <- 'D0_transfo'
  } else if(class(x)=='Seurat'){
    x$orig.ident[x$orig.ident=='488'] <- 'D4_repro'
    x$orig.ident[x$orig.ident=='509'] <- 'D0_repro'
    x$orig.ident[x$orig.ident=='489'] <- 'D7_repro'
    x$orig.ident[x$orig.ident=='547pos'] <- 'D42_transfo'
    x$orig.ident[x$orig.ident=='508'] <- 'poumon_total'
    x$orig.ident[x$orig.ident=='546pos'] <- 'D72_transfo'
    x$orig.ident[x$orig.ident=='467'] <- 'D16_transfo'
    x$orig.ident[x$orig.ident=='510'] <- 'D0_transfo'
  }
  else stop('x should be either a vector or a Seurat object')
  return(x)
}
