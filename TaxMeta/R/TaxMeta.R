#' Meta-analysis useing nested hierarchical models
#'
#' @param response the (column) name of the response variable
#' @param distribution a dsitribution that is supported here and implemented in JAGS (see the JAGS manual for reference). Should be one of
#' \itemize{
#'  \item \code{norm} for a normal distribution
#'  \item \code{gamma} for a gamma distribution
#'  \item \code{lnorm} for a log-normal distribution
#'}
#' @param study_epsilon the column name of individual study standard errors. This parameter is needed for a meta-analysis in the calssic sense, but is not needed if measurements are on the original scale, i.e., includes multiple datapoints for each level at the highest taxonomic resolution considered. If the latter is not the case, and study_epsilon is not available, variability at the highest taxonomic level cannot be estimated.
#' @param taxonomy the column names of taxonomic levels to be considered. To be given in the order from highest resolution (e.g., species) to lowest (e.g., order). For species names, please ensure that these are unique for unique species, or use a [genus species] in the species column.
#' @param cont_cov column names of continuous covariates, if any.
#' @param fixed_fx column names of fixed effects, if any.
#' @param random_fx column names of random effects, if any.
#' @param tree optional, a phylogenetic tree to be combined with taxonomic effects (see Hadfield & Nakagawa, 2010).
#' @param scale the scale parameter for the half-cauchy hierarchical priors. The scale See Gelman 2004 for the original rationale, and ?bayou::qhalfcauchy for an R implementation of this distribution.
#' @param type the type of model to run, one of
#' \itemize{
#'  \item \code{full} fully hierarchical model with taxonomic variances estimated within groups at each taxonomic level. The scale for each the taxonomic model as well as the other random effects is estimated hierachically, from a global hyper-scale parameter.
#'  \item \code{gamma} the scale for all (taxonomic and non-taxonomic) random effects is hierarchically estimated from a vague gamma prior distribution.
#'  \item \code{uniform} the scale for all t(axonomic and non-taxonomic) random effects is hierarchically estimated from a uniform prior distribution.
#'  \item \code{fixed} the scale for all t(axonomic and non-taxonomic) random effects is fixed to \code{scale}.
#'}
#' @param loo_waic pseudo leave-one-out cross validation (see Vehtari et al., 2015) and Wanatabe (or widely applicable) information criterion (Wanatabe, 2010). NOTE that this option leads to large output objects for large datasets, that may fill all your memory (and crash your machine)!
#' @param n.chains number of MCMC chains (default 3)
#' @param n.thin thinning interval of MCMC chains (default 10)
#' @param n.iter number of iterations (default 20e3)
#' @param n.burnin number of burnin iterations to discard (default 5e3)
#' @param return_MCMC return MCMC chains or just summary statistics (default).
#'
#' @return A (list) object of class TaxMeta
#' @author Philipp Neubauer
#'
#' @export
TaxMeta <- function(dataset,
                    response,
                    distribution,
                    study_epsilon = NULL,
                    taxonomy = c('species','genus','family','order'),
                    cont_cov = NULL,
                    fixed_fx = NULL,
                    random_fx = NULL,
                    tree = NULL,
                    scale = NULL,
                    type='full',
                    loo_waic = F,
                    n.chains = 3,
                    n.iter = 25e3,
                    n.thin=1e1,
                    n.burnin=5e3,
                    return_MCMC=F,
                    ...
){


  #if(is.null(scale)) {
  #  scalefunc <- function(scale, max_data) abs(bayou::phalfcauchy(max_data, scale = scale)-0.75)
  #  scale <- optimize(scalefunc,interval = c(0.001,1000),max(abs(dataset[,response])), maximum = F)$objective
  #}

  ##---- impute taxonomy
  if(!all(taxonomy %in% colnames(dataset))){
    stop("\n Full taxonomy not in dataset, please use the 'get_tax' function to get the taxonomy using the taxize package.", immediate. = T)

  }

  ##--- lower case taxonomy names

  tix <- which(colnames(dataset) %in% taxonomy)
  taxonomy <- tolower(taxonomy)
  colnames(dataset)[tix] <- taxonomy


  NOBS = nrow(dataset)
  resp <- dataset[,response]
  # take out predictions from
  not_pred <- which(!is.na(resp))
  nnot_pred <- length(not_pred)
  preds <- which(is.na(resp))
  npreds <- length(preds)

  afn <- function(x) as.numeric(factor(x))

  ##---- test for study epsilon and estimable effects
  if(is.null(study_epsilon)) {

    if(!any(table(dataset[,taxonomy[1]])>1)){
      warning("\n no study_epsilon given and no replication at highest taxonomic resolution; trying to estimate epsilon by dropping highest taxonomic level. If this is not what you want, abort now!",immediate. = T)

      taxonomy <- taxonomy[-1]

    } else {
      warning("\n no study_epsilon given; trying to estimate epsilon",immediate. = T)
    }

  }

  ##---- taxonomic effects
  model_tax <- list()
  ntax <- length(taxonomy)
  taxix <- matrix(NA, NOBS, ntax)
  #precaution
  dataset[,taxonomy] <- apply(dataset[,taxonomy],2,as.character)

  if(!all(grepl(' ',dataset[,'species']))) dataset[!grepl(' ',dataset[,'species']),'species'] <- apply(dataset[!grepl(' ',dataset[,'species']),c('genus','species')],1,paste,collapse=' ')
  taxonomy <- c(taxonomy,"grand")
  for (t in 1:ntax){c
    taxix[,t] <- afn(dataset[,taxonomy[t]])
    model_tax[[t]] <- .parse_tax(taxonomy[t],taxonomy[t+1],taxix[,t],type)
  }


  ##---- random effects
  if(!is.null(random_fx)){
    model_rfx <- list()
    nrfx <- length(random_fx)
    rfxix <- matrix(NA, NOBS, nrfx)
    for (r in 1:nrfx){
      rfxix[,r] <- afn(dataset[,random_fx[r]])
      model_rfx[[r]] <- .parse_rfx(random_fx[r],rfxix[,r])
    }
  }

  ##---- fixed and continuous effects
  if(!is.null(fixed_fx)) {covs <- model.matrix(~factor(dataset[,fixed_fx]))} else {covs <- data.frame(rep(1,NOBS))}
  if(!is.null(cont_cov)) {covs <- cbind(covs,dataset[,cont_cov])}
  ncovs = ncol(covs)

  # sub lognorml
  distribution_sub <- distribution
  if(distribution=='lnorm') distribution_sub = 'norm'

  Model <- cat(
    "model{ \n",
    "#data and prediction likelihood \n",
    "for(i in 1:NOBS){ \n",
    ifelse(loo_waic,paste0("log_lik[i] <- logdensity.",distribution_sub,"(resp[i],",.parse_dist(distribution, study_epsilon),') \n',''),''),
    paste0("resp[i] ~ d",distribution_sub,'(',.parse_dist(distribution, study_epsilon),')'), '\n',
    "mu[i] <-betas[1:ncovs]%*%covs[i,1:ncovs] +",
    ifelse(!is.null(random_fx),paste0(random_fx,'_mu[rfxix[i,',1:nrfx,']]',collapse='+'),''),
    paste0(taxonomy[1:ntax],'_mu[taxix[i,',1:ntax,']]',collapse='+'),"} \n",
    "# likelihood for loo and waic for observations only \n",
    '# save predictions \n',
    ifelse(any(preds),paste0("for(i in 1:npreds){ \n",
                             ifelse(!is.null(study_epsilon),
                                    "pred[i] <- resp[preds[i]] \n} \n",
                                    "pred[i] <- mu[preds[i]] \n} \n")),''),
    "# fixed and regression effects \n",
    .parse.fixed_fx(ncovs),
    "# random effects \n",
    ifelse(!is.null(random_fx),paste(paste(model_rfx),collapse='\n'),''),
    "# taxonomic effects \n",
    paste(paste(model_tax),collapse='\n'),
    "# scale model \n",
    .parse_scale_model(type),
    "# estimate epsilon if needed \n",
    ifelse(is.null(study_epsilon),"epsilon ~ dgamma(1e-6,1e-6) \n",''),
    "}", file='taxmodel.R')


  ###---- run jags model

  params <- .parse.params(type, preds, loo_waic,taxonomy[1:ntax],study_epsilon, random_fx)
  datas <- list(NOBS=NOBS,
                ncovs=ncovs,
                covs=covs,
                taxix=taxix,
                resp=resp
  )

  if(type!='full') datas$scale=scale
  if(any(preds)) {
    datas$preds=preds
    datas$npreds=npreds
  }
  if(!is.null(study_epsilon)) datas$epsilon = 1/dataset[,study_epsilon]^2

  #attach(data.frame(apply(dataset[,taxonomy[1:ntax]],2,afn)))
  if(type=='full') datas[taxonomy[2:ntax]] <- lapply(taxonomy[2:ntax],get,data.frame(apply(dataset[,taxonomy[2:ntax]],2,afn)))

  if(!is.null(random_fx)){
    #attach(data.frame(apply(dataset[,random_fx[1:nrfx]],2,afn)))
    #datas[random_fx[1:ntax]] <- lapply(random_fx[1:nrfx],get,2)
    datas$rfxix=rfxix
  }

  TM <- jags.model(file = 'taxmodel.R',
                   n.chains = n.chains,
                   data=datas)

  update(TM, n.burnin)
  res <- coda.samples(TM, variable.names = params,
                      n.iter = n.iter,thin = n.thin)

  res <- lapply(res, function(re) {
    colnames(re)[grepl('betas',colnames(re))] <- c('Intercept',unlist(apply(data.frame(dataset[,fixed_fx]),2,function(x) sort(unique(x)))),cont_cov)
    re
  })


  cres <- do.call('rbind',res)
  colnames(cres)[grepl('betas',colnames(cres))] <- c('Intercept',unlist(apply(data.frame(dataset[,fixed_fx]),2,function(x) sort(unique(x)))),cont_cov)

  try(diagn <- gelman.diag(as.mcmc(lapply(res,function(x) as.mcmc(x[,!grepl('log_lik',colnames(x))])))))

  if(!return_MCMC) {
    res <- summary(as.mcmc(lapply(res,function(x) as.mcmc(x[,!grepl('log_lik',colnames(x))]))))
  }

  if(loo_waic==T) {
    res_out <- list(
      result=res,
      convergence = diagn,
      waic = loo::waic(cres[,grepl('log_lik',colnames(cres))]),
      loo = loo::loo(cres[,grepl('log_lik',colnames(cres))])
    )
  } else {

    res_out <- list(
      result=res,
      convergence = diagn)

  }
  class(res_out) <- 'TaxMeta'
  res_out
}
