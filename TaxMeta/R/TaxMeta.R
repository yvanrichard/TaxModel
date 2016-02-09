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
#' @param regr_priors A data.frame of mean (1st column) and variances (2nd column) for the regression prior(s). If NULL (default) vague (Normal(0,1e6)) priors are set. Convergence can sometimes be improved with more constrained priors.
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
            save_tax = NULL,
            cont_cov = NULL,
            fixed_fx = NULL,
            random_fx = NULL,
            tree = NULL,
            scale = NULL,
            type='full',
            regr_priors = NULL,
            loo_waic = F,
            n.chains = 3,
            n.iter = 25e3,
            n.thin=1e1,
            n.burnin=5e3,
            return_MCMC=F,
            save_mod_data=F,
            ...
            ){

    timestamp <- format(Sys.time(), '%Y-%m-%d_%H-%M')
    modfile <- paste0(timestamp, "_taxmodel.R")
    
    ##   if(is.null(scale)) {
    ##    scalefunc <- function(scale, max_data) abs(bayou::phalfcauchy(max_data, scale = scale)-0.75)
    ##    scale <- optimize(scalefunc,interval = c(0.001,1000),max(abs(dataset[,response])), maximum = F)$objective
    ##   }

    ##---- impute taxonomy
    if(!all(taxonomy %in% colnames(dataset))){
        stop("\n Full taxonomy not in dataset, please use the 'get_tax' function to get the taxonomy using the taxize package.", immediate. = T)

    }

    ##--- lower case taxonomy names

    tix <- match(taxonomy, colnames(dataset))
    ## just to make sure its not factors!
    dataset[,taxonomy] <- apply(dataset[, taxonomy, drop=FALSE], 2, as.character)
    taxonomy <- tolower(taxonomy)
    colnames(dataset)[tix] <- taxonomy

    NOBS <- nrow(dataset)
    resp <- dataset[,response]
    ## take out predictions from
    not_pred <- which(!is.na(resp))
    nnot_pred <- length(not_pred)
    preds <- which(is.na(resp))
    npreds <- length(preds)

    ##---- test for study epsilon and estimable effects
    if(is.null(study_epsilon)) {

        if(!any(table(dataset[,taxonomy[1]])>1)){
            warning("\n No study_epsilon given and no replication at highest taxonomic resolution; trying to estimate epsilon by dropping highest taxonomic level. If this is not what you want, abort now!",immediate. = T)

            taxonomy <- taxonomy[-1]

        } else {
            warning("\n No study_epsilon given; trying to estimate epsilon",immediate. = T)
        }
    }

    ##---- taxonomic effects
    model_tax <- list()
    ntax <- length(taxonomy)
    taxix <- matrix(NA, NOBS, ntax)
    ##precaution
    ##dataset[,taxonomy] <- apply(dataset[,taxonomy],2,as.character)

    sp.as.genus.sp <- grepl(' ',dataset[,'species'])
    if(!all(sp.as.genus.sp))
        dataset[!sp.as.genus.sp, 'species'] <-
          apply(dataset[!sp.as.genus.sp, c('genus','species')], 1, paste, collapse=' ')
    taxonomy <- c(taxonomy, "grand")

    afn <- function(x) as.numeric(factor(x))

    for (t in 1:ntax){
        taxix[,t] <- afn(dataset[,taxonomy[t]])
        model_tax[[t]] <- .parse_tax(taxonomy[t], taxonomy[t+1], taxix[,t], type)
    }

    ##---- random effects
    if(!is.null(random_fx)){
        model_rfx <- list()
        nrfx <- length(random_fx)
        rfxix <- matrix(NA, NOBS, nrfx)
        for (r in 1:nrfx){
            rfxix[,r] <- afn(dataset[,random_fx[r]])
            model_rfx[[r]] <- .parse_rfx(random_fx[r],rfxix[,r],type)
        }
    }

    ##---- fixed and continuous effects
    if(!is.null(fixed_fx)) {
        if (length(fixed_fx)>1)
            stop('More than one fixed effects specified, which is not yet supported')
        covs <- model.matrix(~factor(dataset[,fixed_fx]))
    } else {
        covs <- data.frame(rep(1,NOBS))
    }
    if(!is.null(cont_cov)) {
        covs <- cbind(covs, dataset[,cont_cov])
    }
    ncovs <- ncol(covs)

    if(!is.null(regr_priors)) {
        regr_mu <- regr_priors[,1]
        regr_tau <- 1/regr_priors[,2]
    } else {
        regr_mu <- rep(0,ncovs)
        regr_tau <- rep(1e-6,ncovs)
    }

    ## sub lognorml
    distribution_sub <- distribution
    if (distribution == 'lnorm')
        distribution_sub <- 'norm'

    Model <- cat(
                "model{ \n\n",
                "# Data and prediction likelihood \n",
                "for(i in 1:NOBS){ \n\t",
                ifelse(loo_waic,
                       paste0("log_lik[i] <- logdensity.",distribution_sub,"(resp[i],",
                              .parse_dist(distribution, study_epsilon),') \n',''), ''),
                paste0("\tresp[i] ~ d",distribution_sub,'(',.parse_dist(distribution, study_epsilon),')'), '\n',
                "\tmu[i] <- betas[1:ncovs] %*% covs[i,1:ncovs] + ",
                ifelse(!is.null(random_fx),
                       paste0(paste0(random_fx,'_mu[rfxix[i,',1:nrfx,']]',collapse=' + '),' + '),''),
                paste0(taxonomy[1:ntax],'_mu[taxix[i,',1:ntax,']]',collapse=' + '),"} \n",
                "\n# Likelihood for loo and waic for observations only \n",
                '# Save predictions \n',
                ifelse(any(preds),paste0("for(i in 1:npreds){ \n",
                                         ifelse(!is.null(study_epsilon),
                                                "\tpred[i] <- resp[preds[i]] \n} \n",
                                                "\tpred[i] <- mu[preds[i]] \n} \n")),''),
                "\n# Fixed and regression effects \n",
                .parse.fixed_fx(ncovs),
                "\n# Random effects \n",
                ifelse(!is.null(random_fx),paste(paste(model_rfx),collapse='\n'),''),
                "\n# Taxonomic effects \n",
                paste(paste(model_tax),collapse='\n'),
                "\n\n# Scale model \n",
                .parse_scale_model(type),
                "\n# Estimate epsilon if needed \n",
                ifelse(is.null(study_epsilon),
                       "epsilon ~ dgamma(1e-6,1e-6) \nsd.epsilon <- sqrt(1/epsilon)",''),
           "}", sep='', file=modfile)


###---- run jags model

    params <- .parse.params(type, preds, loo_waic, taxonomy[1:ntax], study_epsilon, random_fx, save_tax)
    datas <- list(NOBS=NOBS,
                 ncovs=ncovs,
                 covs=covs,
                 taxix=taxix,
                 resp=resp,
                 regr_mu=regr_mu,
                 regr_tau=regr_tau
                 )

    if(type=='fixed')
        datas$scale <- scale
    if(any(preds)) {
        datas$preds <- preds
        datas$npreds <- npreds
    }
    if(!is.null(study_epsilon))
        datas$epsilon <- 1/dataset[,study_epsilon]^2

    ##attach(data.frame(apply(dataset[,taxonomy[1:ntax]],2,afn)))
    if(type=='full')
        datas[taxonomy[2:ntax]] <- lapply(taxonomy[2:ntax], get,
                                         data.frame(apply(dataset[,taxonomy[2:ntax]], 2, afn)))

    if(!is.null(random_fx)){
        ##attach(data.frame(apply(dataset[,random_fx[1:nrfx]],2,afn)))
        ##datas[random_fx[1:ntax]] <- lapply(random_fx[1:nrfx],get,2)
        datas$rfxix <- rfxix
    }

    if (save_mod_data)
        save(datas, file=paste0(timestamp, '_model_data.rdata'))

    ## Run model
    TM <- jags.model(file     = modfile,
                    n.chains = n.chains,
                    data     = datas)

    update(TM, n.burnin)
    res <- coda.samples(TM, variable.names = params,
                       n.iter = n.iter, thin = n.thin)

    res <- lapply(res, function(re) {
        colnames(re)[grepl('betas',colnames(re))] <-
          c('Intercept', unlist(apply(data.frame(dataset[,fixed_fx]), 2, function(x) sort(unique(x))))[-1], cont_cov)
        re
    })

    cres <- do.call('rbind',res)

    res <- coda::as.mcmc(lapply(res, function(x) as.mcmc(x[,!grepl('log_lik',colnames(x))])))

    diagn <- try(gelman.diag(coda::as.mcmc(lapply(res, function(x) coda::as.mcmc(x[, !grepl('pred',colnames(x))])))))

    if(!return_MCMC) {
        res <- summary(res)
    }

    res_out <- list(taxonomy = taxonomy, random_fx = random_fx, 
                   result = res, convergence = diagn)
    if (loo_waic == T) {
        res_out <- c(res_out, list(waic = loo::waic(cres[, grepl("log_lik", colnames(cres))]),
                                  loo = loo::loo(cres[, grepl("log_lik", colnames(cres))])))
    }
    class(res_out) <- 'TaxMeta'
    res_out
}
