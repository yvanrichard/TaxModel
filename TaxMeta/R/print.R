#' print results a taxonomic analysis
#' @name print
#' @param res a result of class TaxMeta
#'
#' @author Philipp Neubauer
#' @references Neubauer.P. Minto, C, and Jensen, O.P. (in prep)
#' @seealso \code{\link{TaxMeta}}
#' @export
NULL

#' @method print TaxMeta
#' @export
print.TaxMeta <- function(res){

  if(all(!is(res$result)=="summary.mcmc")) {

    result <- summary(as.mcmc(lapply(res$result,function(x) as.mcmc(x[,!grepl('log_lik',colnames(x))]))))
  } else {
    result <- res$result
  }

  print(result)

  print(res$convergence)


  cat('\n','\n','waic',res$waic$waic,'waic SE', res$waic$se_waic,'\n')
  cat('\n','\n','loo',res$loo$loo,'loo SE', res$loo$se_loo,'\n')

}
