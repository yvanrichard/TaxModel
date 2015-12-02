jags <- function (data, inits, parameters.to.save, model.file = "model.bug", 
          n.chains = 3, n.iter = 2000, n.burnin = floor(n.iter/2), 
          n.thin = max(1, floor((n.iter - n.burnin)/1000)), DIC = TRUE, 
          working.directory = NULL, jags.seed = 123, refresh = n.iter/50, 
          progress.bar = "text", digits = 5, RNGname = c("Wichmann-Hill", 
                                                         "Marsaglia-Multicarry", "Super-Duper", "Mersenne-Twister"), 
          jags.module = c("glm", "dic")) 
{
  if (!is.null(working.directory)) {
    working.directory <- path.expand(working.directory)
    savedWD <- getwd()
    setwd(working.directory)
    on.exit(setwd(savedWD))
  }
  else {
    savedWD <- getwd()
    working.directory <- savedWD
  }
  if (is.character(data) && length(data) == 1 && regexpr("\\.txt$", 
                                                         data) > 0) {
    if (all(basename(data) == data)) {
      fn2 <- file.path(working.directory, data)
      if (normalizePath(fn2) != normalizePath(data)) {
        try(file.copy(fn2, data, overwrite = TRUE))
      }
    }
    if (!file.exists(data)) {
      stop("File", data, "does not exist")
    }
    if (file.info(data)["size"] == 0) {
      stop("Empty data file ", data)
    }
    e <- new.env()
    eval(parse(data), e)
    data <- as.list(e)
  }
  else if (is.character(data) || (is.list(data) && all(sapply(data, 
                                                              is.character)))) {
    dlist <- lapply(as.list(data), get, envir = parent.frame(1))
    names(dlist) <- unlist(data)
    data <- dlist
  }
  else if (!is.list(data)) {
    stop("data must be a character vector of object names, a list of object names, or a list of objects")
  }
  if (is.function(model.file)) {
    temp <- tempfile("model")
    temp <- if (is.R() || .Platform$OS.type != "windows") {
      paste(temp, "txt", sep = ".")
    }
    else {
      gsub("\\.tmp$", ".txt", temp)
    }
    write.model(model.file, con = temp, digits = digits)
    model.file <- gsub("\\\\", "/", temp)
    if (!is.R()) 
      on.exit(file.remove(model.file), add = TRUE)
  }
  if (DIC) {
    parameters.to.save <- c(parameters.to.save, "deviance")
    rjags::load.module("dic", quiet = TRUE)
  }
  
  if (!missing(inits) && !is.function(inits) && !is.null(inits) && 
      (length(inits) != n.chains)) {
    stop("Number of initialized chains (length(inits)) != n.chains")
  }
  RNGname <- match.arg(RNGname)
  if (RNGname %in% c("Wichmann-Hill", "Marsaglia-Multicarry", 
                     "Super-Duper", "Mersenne-Twister")) {
    RNGname <- paste("base::", RNGname, sep = "")
  }
  else {
    stop("The name of the RNG is not correctly provided!")
  }
  if (!is.null(jags.module)) {
    n.module <- length(jags.module)
    for (m in 1:n.module) {
      rjags::load.module(jags.module[m])
    }
  }
  init.values <- vector("list", n.chains)
  if (missing(inits)) {
    for (i in 1:n.chains) {
      init.values[[i]]$.RNG.name <- RNGname
      init.values[[i]]$.RNG.seed <- runif(1, 0, 2^31)
    }
  }
  else if (is.null(inits)) {
    for (i in 1:n.chains) {
      init.values[[i]]$.RNG.name <- RNGname
      init.values[[i]]$.RNG.seed <- runif(1, 0, 2^31)
    }
  }
  else if (is.function(inits)) {
    if (any(names(formals(inits)) == "chain")) {
      for (i in 1:n.chains) {
        init.values[[i]] <- inits(chain = i)
        init.values[[i]]$.RNG.name <- RNGname
        init.values[[i]]$.RNG.seed <- runif(1, 0, 2^31)
      }
    }
    else {
      for (i in 1:n.chains) {
        init.values[[i]] <- inits()
        init.values[[i]]$.RNG.name <- RNGname
        init.values[[i]]$.RNG.seed <- runif(1, 0, 2^31)
      }
    }
  }
  else {
    if (!is.list(inits)) {
      stop("Invalid inits")
    }
    if (length(inits) != n.chains) {
      stop("Number of initialized chains (length(inits)) != n.chains")
    }
    for (i in 1:n.chains) {
      init.values[[i]] <- inits[[i]]
      init.values[[i]]$.RNG.name <- RNGname
      init.values[[i]]$.RNG.seed <- runif(1, 0, 2^31)
    }
  }
  m <- rjags::jags.model(model.file, data = data, inits = init.values, 
                  n.chains = n.chains, n.adapt = 0)
 
  rjags:::update.jags(m,n.burnin)
  samples <- rjags::coda.samples(model = m, variable.names = parameters.to.save, 
                          n.iter = (n.iter - n.burnin), thin = n.thin, by = refresh, 
                          progress.bar = progress.bar)
  fit <- R2jags:::mcmc2bugs(samples, model.file = model.file, program = "jags", 
                   DIC = DIC, DICOutput = NULL, n.iter = n.iter, n.burnin = n.burnin, 
                   n.thin = n.thin)
  out <- list(model = m, BUGSoutput = fit, parameters.to.save = parameters.to.save, 
              model.file = model.file, n.iter = n.iter, DIC = DIC)
  class(out) <- "rjags"
  return(out)
}

jags.parallel <- function (data, inits, parameters.to.save, model.file = "model.bug", 
          n.chains = 2, n.iter = 2000, n.burnin = floor(n.iter/2), 
          n.thin = max(1, floor((n.iter - n.burnin)/1000)), n.cluster = n.chains, 
          DIC = TRUE, working.directory = NULL, jags.seed = 123, digits = 5, 
          RNGname = c("Wichmann-Hill", "Marsaglia-Multicarry", "Super-Duper", 
                      "Mersenne-Twister"), jags.module = c("glm", "dic"), export_obj_names = NULL, 
          envir = .GlobalEnv) 
{
  jags.params <- parameters.to.save
  jags.inits <- if (missing(inits)) {
    NULL
  }
  else {
    inits
  }
  jags.model <- model.file
  if (is.character(data) && length(data) == 1 && regexpr("\\.txt$", 
                                                         data) > 0) {
    if (all(basename(data) == data)) {
      fn2 <- file.path(working.directory, data)
      if (normalizePath(fn2) != normalizePath(data)) {
        try(file.copy(fn2, data, overwrite = TRUE))
      }
    }
    if (!file.exists(data)) {
      stop("File", data, "does not exist")
    }
    if (file.info(data)["size"] == 0) {
      stop("Empty data file ", data)
    }
    e <- new.env()
    eval(parse(data), e)
    data <- as.list(e)
  }
  else if (is.character(data) || (is.list(data) && all(sapply(data, 
                                                              is.character)))) {
    dlist <- lapply(as.list(data), get, envir = parent.frame(1))
    names(dlist) <- unlist(data)
    data <- dlist
  }
  else if (!is.list(data)) {
    stop("data must be a character vector of object names, a list of object names, or a list of objects")
  }
  list2env(data, envir = envir)
  .runjags <- function() {
    jagsfit <- jags(data = eval(expression(data)), inits = jags.inits, 
                    parameters.to.save = eval(expression(jags.params)), 
                    model.file = eval(expression(jags.model)), n.chains = 1, 
                    n.iter = eval(expression(n.iter)), n.burnin = eval(expression(n.burnin)), 
                    n.thin = eval(expression(n.thin)), DIC = eval(expression(DIC)), 
                    working.directory = eval(expression(working.directory)), 
                    jags.seed = eval(expression(jags.seed)), progress.bar = "none", 
                    digits = eval(expression(digits)), RNGname = eval(expression(RNGname)), 
                    jags.module = eval(expression(jags.module)))
    return(jagsfit)
  }
  cl <- parallel::makeCluster(n.cluster, methods = FALSE)
  parallel::clusterExport(cl, c(names(data), "mcmc", "mcmc.list", export_obj_names,'jags'), 
                envir = envir)
  parallel::clusterSetRNGStream(cl, jags.seed)
  tryCatch(res <- parallel::clusterCall(cl, .runjags), finally = parallel::stopCluster(cl))
  result <- NULL
  model <- NULL
  for (ch in 1:n.chains) {
    result <- abind::abind(result, res[[ch]]$BUGSoutput$sims.array, 
                    along = 2)
    model[[ch]] <- res[[ch]]$model
  }
  if (is.function(model.file)) {
    model.file <- substitute(model.file)
  }
  result <- R2jags:::as.bugs.array2(result, model.file = model.file, 
                           program = "jags", DIC = DIC, n.iter = n.iter, n.burnin = n.burnin, 
                           n.thin = n.thin)
  out <- list(model = model, BUGSoutput = result, parameters.to.save = parameters.to.save, 
              model.file = model.file, n.iter = n.iter, DIC = DIC)
  class(out) <- c("rjags.parallel", "rjags")
  return(out)
}