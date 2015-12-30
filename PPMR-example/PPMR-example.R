require(dplyr)

ppmr.data <- ppmr.best

# need numbers
n <- function(x) length(unique(x))
nidx <- function(x) as.numeric(factor(x,levels=unique(x)))
ind = nidx(ppmr.data$Individual.ID)

## matching taxonomy for taxonomy indices
# match unique individual to species (Predator)
SPECIES <- nidx(ppmr.data$Predator)

# match unique species to family 
FAM <- nidx(ppmr.data$Family)

# match unique family to order
ORD <- nidx(ppmr.data$Order)

# match unique order to class
CL <- nidx(ppmr.data$Class)

ppmr.data <- ppmr.data %>% mutate(log_ten_Pred = log10(Predator.mass))

lpred <- 50
pred = sample.int(n = length(ppmr.data$log_ten_Pred),lpred)
PPMRdata_rfx_h <- list(log10Pmass = ppmr.data$log_ten_Pred,
                     log10PPMR = log10(ppmr.data$Predator.mass/ppmr.data$Prey.mass),
                     geoarea = nidx(ppmr.data$Geographic.location),
                     pred =pred ,
                     NOBS = nrow(ppmr.data),
                     NIND = n(ppmr.data$Individual.ID),
                     NPRED = lpred,
                     NSPECIES = n(ppmr.data$Predator),
                     NFAMILIES = n(ppmr.data$Family),
                     NORDERS = n(ppmr.data$Order),
                     NCLASSES = n(ppmr.data$Class),
                     GEOAREAS = n(ppmr.data$Geographic.location),
                     species = SPECIES,
                     class = CL,
                     family = FAM,
                     order = ORD,
                     ind = nidx(ppmr.data$Individual.ID),
                     species_pred = SPECIES[pred],
                     class_pred = CL[pred],
                     family_pred = FAM[pred],
                     order_pred = ORD[pred],
                     ind_pred = nidx(ppmr.data$Individual.ID)[pred]
)

require(R2jags)
source('jags.R')

PPMR_rfx <- jags.parallel(model.file = 'PPMR-model_rfx_r.R',
                        n.iter = 250000,
                        n.burnin = 50000,
                        DIC = T,
                        n.thin = 100,
                        data=PPMRdata_rfx_h,
                        n.chains = 2,
                        parameters.to.save = c('sd.species',
                                               'sd.family',
                                               'sd.order',
                                               'grandmu',
                                               'pred_i',
                                               'pred_f',
                                               'pred_o',
                                               'pred_s',
                                               'sigma.species',
                                               'sigma.genus',
                                               'sigma.family',
                                               'sigma.order',
                                               'sd.ind',
                                               'sd.ppmr',
                                               'sd.geoareafx',
                                               'beta',
                                               'hc_scale'))



#DD_rfx
#traceplot(DD_rfx)

PPMR_rfx_h <- jags.parallel(model.file = 'PPMR-model_rfx_rh.R',
                            n.iter = 250000,
                            n.burnin = 50000,
                            DIC = T,
                            n.thin = 100,
                            data=PPMRdata_rfx_h,
                            n.chains = 2,
                            parameters.to.save = c('sd.species',
                                                   'sd.ind',
                                                   'sd.family',
                                                   'sd.order',
                                                   'order.scale',
                                                   'family.scale',
                                                   'species.scale',
                                                   'grandmu',
                                                   'pred_i',
                                                   'pred_f',
                                                   'pred_o',
                                                   'pred_s',
                                                   'sd.sigma.species',
                                                   'sd.sigma.order',
                                                   'sd.sigma.family',
                                                   'sd.sigma.ind',
                                                   'sd.ppmr',
                                                   'sd.geoareafx',
                                                   'beta',
                                                   'hc_scale'))


as.mcmc.rjags <- function (x, subs=NULL, ...) 
{
  x <- x$BUGSoutput
  mclist <- vector("list", x$n.chains)
  mclis <- vector("list", x$n.chains)
  strt <- x$n.burnin + 1
  end <- x$n.iter
  if (is.null(subs)) {
    ord <- order(dimnames(x$sims.array)[[3]])
  } else {
    ord <- which(!grepl(subs,dimnames(x$sims.array)[[3]]))
  }
  for (i in 1:x$n.chains) {
    tmp1 <- x$sims.array[, i, ord]
    mclis[[i]] <- mcmc(tmp1, start = strt, end = end, thin = x$n.thin)
  }
  as.mcmc.list(mclis)
}

plot(as.mcmc(PPMR_rfx_h,subs='pred'))
crosscorr(as.mcmc(PPMR_rfx_h,subs='pred'))

save(PPMR_rfx,PPMR_rfx_h,PPMRdata_rfx_h,ppmr.best,file='PPMR.model.runs.Rdata')

#DD_rfx_h
#traceplot(DD_rfx_h)
