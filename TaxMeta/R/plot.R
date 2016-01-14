#' plot the contribution of taxonomic levels to variance in a taxonomic analysis
#' @name plot
#' @param res a result of class TaxMeta
#' @param taxonomy the taxonomy used to estimate res, must be entered as a named argument
#'
#' @author Philipp Neubauer
#' @references Neubauer.P. Minto, C, and Jensen, O.P. (in prep)
#' @seealso \code{\link{TaxMeta}}
#' @export
NULL

#' @method plot TaxMeta
#' @export
plot.TaxMeta <- function(...,taxonomy){

  if(missing(taxonomy)) stop('Please supply taxonomy as a named argument')
  if(nargs()==2){

    msd <- .summarise.res(...,taxonomy)

    ggplot(msd) +
      geom_point(aes(x=Factor, y=means), size=4) +
      geom_linerange(aes(x=Factor, y=means,ymin=q1,ymax=q3),size=1) +
      geom_linerange(aes(x=Factor, y=means,ymin=q11,ymax=q33),size=2) +
      theme_classic() +
      coord_flip() +
      #ylab(expression(Finite~population~variance~(log[10]~PPMR))) +
      ylab('Proportion of variance')

  } else {

    sds <- list(...)

    msd <- lapply(sds, .summarise.res, taxonomy)
    msd <- lapply(1:length(msd), function(x) {
      msd[[x]]$Model <- names(sds)[x]
      as.data.frame(msd[[x]])})

    msd <- do.call('rbind', msd)

    msd$Model <- factor(sapply(as.character(msd$Model),.simpleCap),
                       levels = sapply(names(sds),.simpleCap))

    cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#999999")

    dw=0.75
    ggplot(msd) +
      geom_point(aes(x=Factor, y=means, col=Model), size=4, position=position_dodge(width=dw)) +
      geom_linerange(aes(x=Factor, y=means,ymin=q1,ymax=q3,col=Model),size=1,position=position_dodge(width=dw)) +
      geom_linerange(aes(x=Factor, y=means,ymin=q11,ymax=q33,col=Model),size=2,position=position_dodge(width=dw)) +
      scale_colour_manual('Model',values=cbPalette) +
      theme_classic() +
      coord_flip() +
      #ylab(expression(Finite~population~variance~(log[10]~PPMR))) +
      ylab('Proportion of variance')

  }
}
