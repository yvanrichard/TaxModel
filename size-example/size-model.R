model{

  
  for(p in 1:NPRED){
    pred_s[p] ~  dnorm(mu_pred_s[p],tau[p])
    mu_pred_s[p] <- genusmu[genus_pred[p]] + genus.xi[genus_pred[p]]*species.eta[p]
    pred_g[p] ~  dnorm(mu_pred_g[p],tau[p])
    mu_pred_g[p] <- genusmu[genus_pred[p]] + genus.xi[genus_pred[p]]*species.eta_pred[p]
    pred_f[p] ~  dnorm(mu_pred_f[p],tau[p])
    mu_pred_f[p] <- genus_mu_predict[p] + genus.xi_pred*species.eta_pred[p]
    pred_o[p] ~  dnorm(mu_pred_o[p],tau[p])
    mu_pred_o[p] <- genus_mu_predict_u[p] + genus.xi_pred*species.eta_pred[p]
    
    species_predict[p] <- genus.xi[genus_pred[p]]*species.eta_pred[p]
    species.eta_pred[p] ~ dnorm(0,genus.prec)
    
    genus_mu_predict[p] <- familymu[family_pred[p]]+family.xi[family_pred[p]]*genus.eta_pred[p]
    genus_mu_predict_u[p] <- family_mu_predict[p] + family.xi_pred*genus.eta_pred[p]
    genus.eta_pred[p] ~ dnorm(0,family.prec)
    
    family_mu_predict[p] <- ordermu[order_pred[p]] + order.xi[order_pred[p]]*family.eta_pred[p]
    family.eta_pred[p] ~ dnorm(0,order.prec)
    
    
  }
  
    for(s in 1:NSPECIES){
    lmass[s] ~ dnorm(speciesmu[s],tau[s])
    # species mean drawn from family dist 
    speciesdev[s] <- speciesmu[s] - genusmu[genus[s]]
    speciesmu[s] <- genusmu[genus[s]] + genus.xi[genus[s]]*species.eta[s]
    
    #l_obs_sp[s] ~ dnorm(speciesmu[s]+species.xi[s]*sp_pred.eta[s],tau)
    # species tau drawn from hierarchical distr over all families
    #sp_pred.eta[s] ~ dnorm(0,species.prec)
    species.eta[s] ~ dnorm(0,genus.prec)
  }
  
  # genus fx from order distributions
  for(g in 1:NGENUS){
    # genus mean drawn from family dist
    genusdev[g] <- genusmu[g] - familymu[family[g]]
    genusmu[g] <- familymu[family[g]]+family.xi[family[g]]*genus.eta[g]
    #l_obs_gen[g] ~ dnorm(genusmu[g]+genus.xi[g]*gen_pred.eta[g],tau)
    # genus tau drawn from hierarchical distr over all families
    #gen_pred.eta[g] ~ dnorm(0,genus.prec)
    sigma.genus[g] <- abs(genus.xi[g])/sqrt(genus.prec)
    genus.eta[g] ~ dnorm(0,family.prec)
    genus.xi[g] ~ dnorm(0,genus.scale)
  }
  
  # family fx from order distribution
  for(f in 1:NFAMILIES){
    # family mean drawn from order dist
    familydev[f] <- familymu[f] - ordermu[order[f]]
    familymu[f] <- ordermu[order[f]]+order.xi[order[f]]*family.eta[f]
    #l_obs_fam[f] ~ dnorm(familymu[f]+family.xi[f]*fam_pred.eta[f],tau)
    # family tau drawn from hierarchical distr over all families
    #fam_pred.eta[f] ~ dnorm(0,family.prec)
    sigma.family[f] <- abs(family.xi[f])/sqrt(family.prec)
    family.eta[f] ~ dnorm(0,order.prec)
    family.xi[f] ~ dnorm(0,family.scale)
  }
  
  # order fx from class distribution
  for(o in 1:NORDERS){
    # order mean drawn from order dist
    ordermu[o] <- grandmu+grand.xi*order.eta[o]
    #l_obs_ord[o] ~ dnorm(ordermu[o]+order.xi[o]*ord_pred.eta[o],tau)
    # order tau drawn from hierarchical distr over all orders
    #ord_pred.eta[o] ~ dnorm(0,order.prec)
    order.eta[o] ~ dnorm(0,grand.prec)
    sigma.order[o] <- abs(order.xi[o])/sqrt(order.prec)
    
    order.xi[o] ~ dnorm(0,order.scale)
  }
  
  sigma.grand <- abs(grand.xi)/sqrt(grand.prec)
  
  # priors tax hierachy
  species.scale ~ dgamma(2,hc_scale)
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

  grandmu ~ dnorm(0,1e-6)
  hc_scale ~ dunif(0.00001,100000)
  
  # finite population sds
  sd.sigma.family <- sd(sigma.family)
  sd.sigma.genus <- sd(sigma.genus)
  sd.sigma.order <- sd(sigma.order)
  
  sd.genus <- sd(genusdev)
  sd.species <- sd(speciesdev)
  sd.family  <- sd(familydev)
  sd.order   <- sd(ordermu)
 
  }