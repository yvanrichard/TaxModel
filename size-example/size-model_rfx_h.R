model{
  
  for(p in 1:NPRED){
    pred_s[p] ~  dnorm(mu_pred_s[p],tau[p])
    mu_pred_s[p] <- grandmu + speciesmu[pred[p]] + genusmu[genus_pred[p]] + familymu[family_pred[p]] + ordermu[order_pred[p]]
    pred_g[p] ~  dnorm(mu_pred_g[p],tau[p])
    mu_pred_g[p] <- grandmu + species_predict[p] + genusmu[genus_pred[p]] + familymu[family_pred[p]] + ordermu[order_pred[p]]
    pred_f[p] ~  dnorm(mu_pred_f[p],tau[p])
    mu_pred_f[p] <- grandmu + species_predict_u[p] + genus_predict[p] + familymu[family_pred[p]] + ordermu[order_pred[p]]
    pred_o[p] ~  dnorm(mu_pred_o[p],tau[p])
    mu_pred_o[p] <- grandmu + species_predict_u[p] + genus_predict_u[p] + family_predict[p] + ordermu[order_pred[p]]
    
    species_predict[p] <- genus.xi[genus_pred[p]]*species.eta_pred[p]
    species_predict_u[p] <- genus.xi_pred*species.eta_pred[p]
    species.eta_pred[p] ~ dnorm(0,genus.prec)
    
    genus_predict[p] <- family.xi[family_pred[p]]*genus.eta_pred[p]
    genus_predict_u[p] <- family.xi_pred*genus.eta_pred[p]
    genus.eta_pred[p] ~ dnorm(0,family.prec)
    
    family_predict[p] <- order.xi[order_pred[p]]*family.eta_pred[p]
    family.eta_pred[p] ~ dnorm(0,order.prec)
    
  }
  
  for(s in 1:NSPECIES){
    lmass[s] ~ dnorm(mu[s],tau[s])
    mu[s] <- grandmu + speciesmu[s] + genusmu[genus[s]] + familymu[family[s]] + ordermu[order[s]]
    
    speciesmu[s] <- genus.xi[genus[s]]*species.eta[s]
    species.eta[s] ~ dnorm(0,genus.prec)
    
  }
  
  # genus fx from order distributions
  for(g in 1:NGENUS){
    # genus mean drawn from family dist
    genusmu[g] <- family.xi[family[g]]*genus.eta[g]
    # genus tau drawn from hierarchical distr over all families
    genus.eta[g] ~ dnorm(0,family.prec)
    genus.xi[g] ~ dnorm(0,genus.scale)
    sigma.species[g] <- abs(genus.xi[g])/sqrt(genus.prec)
  }
  
  # family fx from order distribution
  for(f in 1:NFAMILIES){
    # family mean drawn from order dist
     familymu[f] <- order.xi[order[f]]*family.eta[f]
      family.eta[f] ~ dnorm(0,order.prec)
    family.xi[f] ~ dnorm(0,family.scale)
    sigma.genus[f] <- abs(family.xi[f])/sqrt(family.prec)
  }
  
  # order fx from class distribution
  for(o in 1:NORDERS){
    # order mean drawn from order dist
    ordermu[o] <- grand.xi*order.eta[o]
    order.eta[o] ~ dnorm(0,grand.prec)
    order.xi[o] ~ dnorm(0,order.scale)
    sigma.family[o] <- abs(order.xi[o])/sqrt(order.prec)
  }
  
  # priors tax hierachy
  species.prec ~ dgamma(0.5,0.5)
  
  family.scale  ~ dgamma(2,hc_scale)
  family.prec ~ dgamma(0.5,0.5)
  family.xi_pred ~ dnorm(0,family.scale)
  
  genus.scale  ~ dgamma(2,hc_scale)
  genus.prec ~ dgamma(0.5,0.5)
  genus.xi_pred ~ dnorm(0,genus.scale)
  
  order.scale   ~ dgamma(2,hc_scale)
  order.prec ~ dgamma(0.5,0.5)
  
  grand.xi ~ dnorm(0,2/hc_scale)
  grand.prec ~ dgamma(0.5,0.5)
  sigma.order <- abs(grand.xi)/sqrt(grand.prec)
  
  grandmu ~ dnorm(0,1e-6)
  hc_scale ~ dunif(0.00001,100000)
  
  #   # other rfx
  #   
  #   for (h in 1:CONTINENTS){
  #     continent_eff[h] <- cont.xi*cont.eta[h] 
  #     cont.eta[h] ~ dnorm(0,cont.prec)
  #   }
  #   
  #   cont.xi ~ dnorm(0,0.1)
  #   cont.prec ~ dgamma(0.5,0.5)
  
  # finite population sds
  sd.sigma.species <- sd(sigma.species)
  sd.sigma.genus <- sd(sigma.genus)
  sd.sigma.family <- sd(sigma.family)
  
  sd.genus <- sd(genusmu)
  sd.species <- sd(speciesmu)
  sd.family  <- sd(familymu)
  sd.order   <- sd(ordermu)
  # sd.continent_eff <- sd(continent_eff)
  
  #tau ~ dgamma(0.001,0.001)
  
}