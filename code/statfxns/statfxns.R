library(gsubfn)

pval.2tail.np <- function(test.val,dist){
	# test.val: individual value being compared to distribution
	# dist: vector,distribution of values under some null model, or otherwise
	# sig.fig: number of significant figures
	# compute 2-tailed p-value for test value occurring in distribution
	dist <- as.numeric(dist)
	pval.2tail <- 2*min(mean(test.val >= dist),mean(test.val <= dist))
	return(pval.2tail)
}


pval.label.np <- function(pvals,n,sig.fig=2){
	# pval: p-value obtained from non-parametric test
	# n: number of samples in distribution used to obtain p-value
	# sig.fig: number of significant figures to report
	# don't ever say p = 0; instead replace with p < 1/n
	
	# make p-value label, rounded to sig figs
	pval.txt <- paste('p =',signif(pvals,sig.fig))
	# replace p = 0 with p < 1/n
	pval.txt[pvals == 0] <- paste('p <',signif(1/n,sig.fig))
	return(pval.txt)
}

list.posthoc.correct <- function(X,method){
  # unlist a list, posthoc correct over all values according to "method"
  # relist the list in the same structure and return
  return(relist(flesh=p.adjust(unlist(X),method=method),skeleton=X))
}

get.ftest.pval <- function(m){
  # INPUTS:
  # given lm object m extract f-test pvalue
  x <- summary(m)
  return(pf(x$fstatistic[1],x$fstatistic[2],x$fstatistic[3],lower.tail=FALSE))
}

get.rsq <- function(m){
  # INPUTS:
  # given lm object m extract total r-squared
  return(summary(m)$r.squared)
}

get.coef <- function(m){
  # INPUTS:
  # given lm object m extract coefficient table
  return(summary(m)$coef)
}

outlier.mask <- function(x){
  # return mask of values in x that are outliers
  return(x %in% boxplot.stats(x)$out)
}

cohens.f2 <- function(full.model,xname){
  # from Selya et al. 2012
  # f2 = (R2.ab - R2.a) / (1 - R2.ab)
  y <- full.model$model[[1]]
  covariates.ind <- which(names(full.model$model) != xname)[-1] # remove y variable index
  x.c <- full.model$model[covariates.ind]
  R2.ab<-summary(full.model)$r.sq
  R2.a<-summary(lm(y~. , data= x.c))$r.sq
  f.2 <- (R2.ab - R2.a) / (1 - R2.ab)
  return(f.2)
}


# used in plot_subsample_cluster.R
maxprobassignment <- function(partitions){
	# INPUTS:
	# partitions: list whose elements contain cluster assignments for subsamples of data,
	# each element corresponding to a value of k. names of elements are value of k. 
	# 0's mean that the TR wasn't in that subsample

	# OUTPUTS:
	# x.max: list whose elements are:
	#	$maxprob: for each TR, what was the largest proportion of times that TR was assigned to the same cluster?
	#	$k: value of k

	k.rng <- names(partitions)
	x.max <- list()
	for(k in k.rng){
		x <- partitions[[k]]
		k.n <- as.numeric(k)
		x.count <- sapply(1:k.n, function(k) rowSums(x == k)) # count number assigned to cluster 1:k, excluding 0's (when TR wasn't in subsample)
		x.pdf <- sapply(1:k.n, function(k) x.count[,k] / rowSums(x.count)) # count fraction of times a TR was assigned to each cluster, out of all the times it was assigned to any cluster (either b/c no unique match or TR not in subsample)
		x.max[[k]] <- list()
		x.max[[k]]$maxprob <- apply(x.pdf,1,max) # get maximum value for all k values
		x.max[[k]]$k <- k # store k value		
	}
	return(x.max)
}

get.partial.resids <- function(m,xname){
  # INPUTS:
  # m: lm object
  # xname: character of x variable name
  #
  # OUTPUTS:
  # df: with columns x and y to make partial residual plot
  # partial residual wrt x_i is (y-y^) + b_i*x_i
  
  pr <- m$residuals + m$model[,xname]*m$coefficients[xname]
  df <- data.frame(x = m$model[,xname],y=pr)
  return(df)
  
}
