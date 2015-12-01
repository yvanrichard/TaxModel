model{

  for(p in 1:NPRED){

    pred_i[p] ~ dnorm(mu_pred_i[p],tau)
    mu_pred_i[p] <- grandmu + beta*log10Pmass[p]+speciesmu[species_pred[p]] + indmu[ind_pred[p]] + familymu[family_pred[p]] + ordermu[order_pred[p]] + geoareafx[geoarea[p]]
    
    pred_s[p] ~ dnorm(mu_pred_s[p],tau)
    mu_pred_s[p] <- grandmu + beta*log10Pmass[p]+speciesmu[species_pred[p]] + ind_predict[p] + familymu[family_pred[p]] + ordermu[order_pred[p]]+  geoareafx[geoarea[p]]
    
    pred_f[p] ~ dnorm(mu_pred_f[p],tau)
    mu_pred_f[p] <- grandmu + beta*log10Pmass[p]+species_predict[p] + ind_predict[p] + familymu[family_pred[p]] + ordermu[order_pred[p]] +  geoareafx[geoarea[p]]
    
    pred_o[p] ~ dnorm(mu_pred_o[p],tau)
    mu_pred_o[p] <- grandmu + beta*log10Pmass[p]+species_predict[p] + ind_predict[p] + family_predict[p] + ordermu[order_pred[p]] + geoareafx[geoarea[p]]
    
    species_predict[p] <- family.xi*species.eta_pred[p]
    species.eta_pred[p] ~ dnorm(0,family.prec)
    
    ind_predict[p] <- species.xi*ind.eta_pred[p]
    ind.eta_pred[p] ~ dnorm(0,species.prec)
 
    family_predict[p] <- order.xi*family.eta_pred[p]
    family.eta_pred[p] ~ dnorm(0,order.prec)
    
  }
  
  for(i in 1:NOBS){
    log10PPMR[i] ~ dnorm(mu[i],tau)
    mu[i] <- grandmu + beta*log10Pmass[i]+ speciesmu[species[i]] + indmu[ind[i]] + familymu[family[i]] + ordermu[order[i]]+ geoareafx[geoarea[i]]
    
  }
  
  beta ~ dnorm(0,1e-6)
  
  for(i in 1:NSPECIES){
    
    speciesmu[i] <- family.xi*species.eta[i]
    species.eta[i] ~ dnorm(0,family.prec)
    
  }
  
  # ind fx from order distributions
  for(g in 1:NIND){
    # ind mean drawn from family dist
    indmu[g] <- species.xi*ind.eta[g]
    # ind tau drawn from hierarchical distr over all families
    ind.eta[g] ~ dnorm(0,species.prec)
  }
  
  # family fx from order distribution
  for(f in 1:NFAMILIES){
    # family mean drawn from order dist
     familymu[f] <- order.xi*family.eta[f]
      family.eta[f] ~ dnorm(0,order.prec)
  }
  
  # order fx from class distribution
  for(o in 1:NORDERS){
    # order mean drawn from order dist
    ordermu[o] <- grand.xi*order.eta[o]
    order.eta[o] ~ dnorm(0,grand.prec)
    
  }
  
  # priors tax hierachy
  species.xi ~ dnorm(0,10)
  species.prec ~ dgamma(0.5,0.5)
  sigma.ind <- abs(species.xi)/sqrt(species.prec)
  
  family.xi ~ dnorm(0,10)
  family.prec ~ dgamma(0.5,0.5)
  sigma.species <- abs(family.xi)/sqrt(family.prec)

  order.xi ~ dnorm(0,10)
  order.prec ~ dgamma(0.5,0.5)
  sigma.family <- abs(order.xi)/sqrt(order.prec)
  
  grand.xi ~ dnorm(0,10)
  grand.prec ~ dgamma(0.5,0.5)
  sigma.order <- abs(grand.xi)/sqrt(grand.prec) 
  
  grandmu ~ dnorm(-3.6,0.2)
  
  # other rfx 
  geo.xi ~ dnorm(0,10)
  geo.prec ~ dgamma(0.5,0.5)
  
  for (d in 1:GEOAREAS){
    geoareafx[d] <- geo.xi*geo.eta[d] 
    geo.eta[d] ~ dnorm(0,geo.prec)
  }
  
  #   # other rfx
  #   
  #   for (h in 1:CONTINENTS){
  #     continent_eff[h] <- cont.xi*cont.eta[h] 
  #     cont.eta[h] ~ dnorm(0,cont.prec)
  #   }
  #   
  #   cont.xi ~ dnorm(0,0.1)
  #   cont.prec ~ dgamma(0.5,0.5)
  
  # finite populatio∂Ωn sds
  sd.ind <- sd(indmu)
  sd.species <- sd(speciesmu)
  sd.family  <- sd(familymu)
  sd.order   <- sd(ordermu)
  sd.ppmr <- 1/sqrt(tau)
  # sd.continent_eff <- sd(continent_eff)
  
  tau ~ dgamma(1,0.5)
  
}