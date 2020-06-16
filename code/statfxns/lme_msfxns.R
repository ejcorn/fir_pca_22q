
lme.ms <- function(cfg,df){
  # INPUTS:
  # cfg: list with the named elements below
  #   fixed, random: character names of independent variables in fixed and random effects
  #   response: character name of response (default to 'Score' for this script)
  #   id: character name of ids (Default to 'scanid' for this script)
  # df: data frame containing fixed.IVs and random.IVs and response
  # OUTPUTS:
  # lme call with input fixed and random effects
  # use this so I can do easy model selection without copy-pasting the same parameters
  fixed <- reformulate(termlabels = cfg$fixed,response = cfg$response,intercept = T) # make fixed formula
  random <- reformulate(termlabels = cfg$random,response = cfg$response,intercept = T) # make random formula
  random <- as.formula(paste(Reduce(paste,deparse(random)),'|',cfg$id)) # add id variable to random
  ctrl <- lmeControl(opt = 'optim') 
  # parse formulas as text in order to extract coefficients and compare models with anova
  return(        eval(bquote(   lme(fixed=.(fixed),random=.(random), 
                                    data = df,na.action=na.exclude,control=ctrl,method='ML')   )) )
}

lme4.ms <- function(cfg,df){
  # INPUTS:
  # cfg: list with the named elements below
  #   fixed, random: character names of independent variables in fixed and random effects
  #   response: character name of response (default to 'Score' for this script)
  #   id: character name of ids (Default to 'scanid' for this script)
  # df: data frame containing fixed.IVs and random.IVs and response
  # OUTPUTS:
  # lme call with input fixed and random effects
  # use this so I can do easy model selection without copy-pasting the same parameters
  fixed <- reformulate(termlabels = cfg$fixed,response = cfg$response,intercept = T) # make fixed formula
  f <- as.formula(paste(cfg$response,'~',as.character(fixed)[3],'+(',paste(cfg$random,collapse='+'),'|',cfg$id,')'))
  ctrl <- lmerControl(optCtrl = list(method='Nelder-Mead')) 
  # parse formulas as text in order to extract coefficients and compare models with anova
  return(        eval(bquote(   lmer(formula = .(f), 
                                    data = df,na.action=na.exclude,control=ctrl,REML=F)   )) )
}

lme.compare <- function(mdl.gold,mdl.test,coef.sig){
  # INPUT:
  # mdl.gold, mdl.test: compare new model (mdl.test) to current gold std model (mdl.gold)
  # using 3 criteria that all must be true. (c1) p < 0.05 for anova(m1,m2) AND (c2) AIC drops AND (c3) added coefficient is significant
  # coef.sig: name of coefficient to test significance for
  #
  # OUTPUT:
  # cond: true/false is model better
  c1 <- sort(anova(mdl.gold,mdl.test)$`p-value`) < 0.05 # likelihood ratio p < 0.05
  c2 <- summary(mdl.gold)$AIC > summary(mdl.test)$AIC # AIC decreases
  print(paste('checking p-val of',coef.sig))
  print(summary(mdl.test)$tTable[coef.sig,,drop=F])
  c3 <- summary(mdl.test)$tTable[coef.sig,'p-value'] < 0.05 # p-value of most recently added coefficient is < 0.05
  cond <- all(c(c1,c2,c3)) # all conditions must be true to advance
  return(cond)
}

lme.last.coef.sig <- function(coef.names){
  # INPUTS:
  # coef.names: character string of coefficient names
  #
  # OUTPUTS:
  # name or index of last coefficient in table from lme output
  last.coef.name <- rev(coef.names)[1] # get name of last coefficient
  if(grepl('*',last.coef.name)){ # if last coefficient is an interaction
    last.coef.name <- gsub('\\*',':',last.coef.name)
  }
  return(last.coef.name)
}

lme.stepup <- function(mdl.gold=NULL,test.cfg=NULL,lme.ms.fun=lme.ms){
  # INPUTS:
  # mdl.gold: current gold standard model
  # test.cfg: list of fixed and random effects to step through
  # lme.ms: model function that allows formulas to be evaluated as plain text
  #
  # OUTPUTS:
  # gold standard model and associated config
  
  # start with first element of test config and only advance if that model is GOOD
  # application is to start with simple models and only step up if each model is better than the next
  test.counter <- 1 # start at first and advance as long as : (c1) p < 0.05 for anova(m1,m2) AND (c2) AIC drops AND (c3) added coefficient is significant
  cond <- TRUE
  while(cond & test.counter < length(test.cfg)+1){ # while cond is still true keep adding stuff until you've tested all models
    print('testing:')
    cfg <- test.cfg[[test.counter]]
    print(rev(cfg$fixed)[1])
    fitSuccess <- tryCatch({mdl.test <- lme.ms.fun(cfg,df = mdl.gold$data)
        cond <- lme.compare(mdl.gold,mdl.test,coef.sig=lme.last.coef.sig(cfg$fixed))
            },error = function(err) {
              print(err)
              print('halting model advancement now')
              return(FALSE)
            })
    if(fitSuccess){
      if(cond){
        mdl.gold <- mdl.test
        cfg.gold <- cfg} # if test successful update gold standard}
      test.counter <- test.counter+1
    } else {cond <- FALSE}
  }
  return(list(mdl.gold=mdl.gold,cfg.gold=cfg.gold))
}

lme.stepdown <- function(mdl.gold=NULL,test.cfg=NULL,lme.ms.fun=lme.ms){
  # INPUTS:
  # mdl.gold: current gold standard model
  # test.cfg: list of fixed and random effects to step through
  # lme.ms: model function that allows formulas to be evaluated as plain text
  #
  # OUTPUTS:
  # gold standard model and associated config
  
  # start with first element of test config and only advance if that model is BAD
  # application is to start with complex models then step down until you find one better than gold
  test.counter <- 1 # start at first and advance as long as : NOT [ (c1) p < 0.05 for anova(m1,m2) AND (c2) AIC drops AND (c3) added coefficient is significant]
  adv <- TRUE # this tells you whether to keep advancing
  while(adv & test.counter < length(test.cfg)+1){ # while cond is still true keep adding stuff until you've tested all models
    print('testing:')
    cfg <- test.cfg[[test.counter]]
    print(rev(cfg$fixed)[1])
    fitSuccess <- tryCatch({mdl.test <- lme.ms.fun(cfg,df = mdl.gold$data)
    cond <- lme.compare(mdl.gold,mdl.test,coef.sig=lme.last.coef.sig(cfg$fixed)) # check if tested model is better than gold std
    },error = function(err) {
      print(err)
      print('halting model advancement now')
      return(FALSE)
    })
    if(fitSuccess){
      if(cond){ # if tested model IS better than gold, then store results and stop advancing
        mdl.gold <- mdl.test
        cfg.gold <- cfg # if test successful update gold standard
        adv <- FALSE
      } else if(!cond){adv <- TRUE} # if tested model is NOT better, then keep going
    } else {adv <- TRUE} # if you couldn't fit the first model, try again with the next less complex model
    test.counter <- test.counter+1
  }
  return(list(mdl.gold=mdl.gold,cfg.gold=cfg.gold))
}

lme.selectbest <- function(mdl.gold=NULL,cfg.gold=NULL,test.cfg=NULL,lme.ms.fun=lme.ms){
  # INPUTS:
  # mdl.gold: current gold standard model
  # cfg.gold: configuration of current gold standard model
  # test.cfg: list of fixed and random effects to step through
  # lme.ms: model function that allows formulas to be evaluated as plain text
  #
  # OUTPUTS:
  # gold standard model and associated config

  mdl.tests <- lapply(test.cfg, function(cfg) tryCatch({return(lme.ms(cfg,df = mdl.gold$data))},
                                          error=function(err){
                                            return('')
                                          }))
  if(sum(mdl.tests != '')>0){ # if any of the tested models fit successfully
    mdl.tests <- mdl.tests[mdl.tests != ''] # delete models that couldn't be fit
    # first check if any models are doing better than gold standard
    test.vs.gold <- sapply(1:length(mdl.tests), function(mdl.idx) lme.compare(mdl.gold,mdl.tests[[mdl.idx]],coef.sig=lme.last.coef.sig(test.cfg[[mdl.idx]]$fixed)))
    #mdl.all <- c(list(mdl.gold),mdl.tests) # these should be in order
    # first column is whether each model is better than mdl.gold
    if(!any(test.vs.gold)){ # if no test model is better than mdl.gold, return mdl.gold and its cfg list
      return(list(mdl.gold=mdl.gold,cfg.gold=cfg.gold))
    } else{ # otherwise, of the models better than gold, return the model that was better than the most other models tested 
      mdl.tests <- mdl.tests[test.vs.gold] # delete models not better than gold standard
      test.cfg <- test.cfg[test.vs.gold] # delete the configuration list elements associated with those models
      perf.comp <- matrix(NA,nrow=length(mdl.tests),ncol=length(mdl.tests)) # pairwise comparison of models that beat the gold standard
      for(mdl.idx1 in 1:length(mdl.tests)){
        for(mdl.idx2 in 1:length(mdl.tests)){
          perf.comp[mdl.idx1,mdl.idx2] <- lme.compare(mdl.tests[[mdl.idx2]],mdl.tests[[mdl.idx1]],coef.sig=lme.last.coef.sig(test.cfg[[mdl.idx1]]$fixed))
        }
      } # output mdl1 x mdl2 matrix of comparisons where rows indicate model being tested for coefficient significance and AIC drop and anova      
      mdl.idx <- which.max(rowSums(perf.comp))
      if(length(mdl.idx) > 1){print('ERROR: models tied'); return()}
      if(length(mdl.idx) == 0){print('ERROR: all tested models are equivalent'); return()}
      mdl.gold <- mdl.tests[[mdl.idx]]
      cfg.gold <- test.cfg[[mdl.idx]]
      return(list(mdl.gold=mdl.gold,cfg.gold=cfg.gold))
    }
  } else {return(list(mdl.gold=mdl.gold,cfg.gold=cfg.gold))}
}


lme.fixed.predict.ave <- function(mdl,fixed.terms){
  # INPUTS:
  # mdl: lme model object
  # fixed.terms: characters specifying which terms to get fixed effect overall line (no random effects)
  # using group average values for all other covariates
  #
  # OUTPUTS:
  # for every combination of levels in fixed.terms
  # get fitted values of model using fixed effects and average values for covariates 
  # return response variable averaged as well as the fitted values with .fixed appended
  
  df <- mdl$data[complete.cases(mdl$data),] # get complete cases
  df.agg <- aggregate(df,by=setNames(lapply(fixed.terms, function(term) df[,term]),fixed.terms),mean)
  df.agg[,'(Intercept)'] <- 1 # add intercept
  fixed.coefs <- mdl$coefficients$fixed
  fixed.coef.names <- names(fixed.coefs)
  if(any(grepl(':',fixed.coef.names))){ # if there are any interactions, make a new column in df.agg for the interaction
    for(int.idx in grep(':',fixed.coef.names)){
      coef.int <- fixed.coef.names[int.idx] # get name of interaction (var1:var2)
      int.components <- strsplit(coef.int,':',fixed = T)[[1]] # split at colon
      # create interaction predictor column by multiplying variables together
      # doing this with averaging only really makes sense with reasonably discrete variables
      df.agg[,coef.int] <- df.agg[,int.components[1]]*df.agg[,int.components[2]]
    }
  }
  resp.name <- attr(getResponse(mdl),'label')
  df.agg[,paste0(resp.name,'.fixed')] <- sapply(1:nrow(df.agg), function(j) sum(fixed.coefs * df.agg[j,fixed.coef.names]))
  return(df.agg[,c(resp.name,paste0(resp.name,'.fixed'),fixed.terms)])
}

lme.subj.coefs <- function(mdl,c.o.i){
  # INPUTS:
  # mdl: lme model object
  # c.o.i: character specifying which terms to get subject-specific parameters for
  #
  # OUTPUTS:
  # for coefficients in c.o.i
  # get coefficients for each subject by combining fixed and random effects, and interactions if present
  # return dataframe of subject-specific coefficients with subject identifiers
  # ** assumes you only have one id variable for random effects **
  
  fixed.coefs <- mdl$coefficients$fixed
  random.coefs <- mdl$coefficients$random[[1]] # ** assumes you only have one id variable for random effects **
  subj.coefs <- intersect(names(fixed.coefs),colnames(random.coefs)) 
  return(fixed.coefs[subj.coefs] + random.coefs[2,subj.coefs])
  
  if(any(grepl(':',fixed.coef.names))){ # if there are any interactions, make a new column in df.agg for the interaction
    for(int.idx in grep(':',fixed.coef.names)){
      coef.int <- fixed.coef.names[int.idx] # get name of interaction (var1:var2)
      int.components <- strsplit(coef.int,':',fixed = T)[[1]] # split at colon
      # create interaction predictor column by multiplying variables together
      # doing this with averaging only really makes sense with reasonably discrete variables
      df.agg[,coef.int] <- df.agg[,int.components[1]]*df.agg[,int.components[2]]
    }
  }
  resp.name <- attr(getResponse(mdl),'label')
  df.agg[,paste0(resp.name,'.fixed')] <- sapply(1:nrow(df.agg), function(j) sum(fixed.coefs * df.agg[j,fixed.coef.names]))
  return(df.agg[,c(resp.name,paste0(resp.name,'.fixed'),fixed.terms)])
  
  
}