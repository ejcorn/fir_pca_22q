rename.column <- function(df,old_name,new_name){
	# rename column of data frame df that was called old_name to new_name
	colnames(df)[colnames(df) == old_name] <- new_name
	return(df)
}

remove.column <- function(df,cname){return(df[,-which(colnames(df) == cname)])}

collapse.columns <- function(df,cnames=colnames(df),groupby=NULL){
  # INPUTS:
  # df: dataframe
  # cnames: column names to perform operation on, default to all columns
  # groupby: column name to group variables by, treated separately from variables in cnames
  
  # OUTPUTS:
  # df.new: dataframe with 2 columns:
  # values: all columns in cnames vertically concatenated. 
  # names: elements of cnames corresponding to rows
  # groupby: groups of observations in df for each variable in cnames, these retain their names in df.new
  
  df.names <- do.call('cbind',lapply(cnames, function(n) rep(n,nrow(as.matrix(df[,cnames])))))  
  df.new <- data.frame(values = as.vector(as.matrix(df[,cnames])),names=as.vector(df.names),stringsAsFactors=F)
  if(!is.null(groupby)){
    for(group.var in groupby){
      df.grp <- do.call('cbind',lapply(cnames,function(n) df[,group.var]))
      df.new[,group.var] <- as.vector(df.grp)  
    }
    
  }
  return(df.new)
}

name <- function(x,x.names){
  # INPUTS:
  # x: vector or dataframe
  # x.names: names for elements or columns of x
  # OUTPUTS:
  # x with x.names as names
  names(x) <- x.names
  return(x)
}

paste0.combinations <- function(x,y){
  # INPUTS:
  # x: character vector
  # y: character vector
  #
  # OUTPUTS:
  # xy: every element of x pasted to every element of y, with x occuring first

  xy <- as.vector(sapply(x, function(i) 
    sapply(y, function(j) paste0(i,j))))
  return(xy)
}

tp.vec2mat <- function(v){
  # INPUTS:
  # v: vectorized transition probabilities
  # 
  # OUTPUTS:
  # m: tp matrix
  
  return(t(matrix(v,sqrt(length(v)),sqrt(length(v)))))
}

matrix.to.df <- function(m,dnames){
  # INPUTS:
  # m: matrix
  # dnames: 2 element list specifying dimnames for m
  #
  # OUTPUTS:
  # df:
  dimnames(m) <- dnames
  return(as.data.frame(m))

}

axis.sym <- function(x){
  # INPUTS:
  # x: vector
  # 
  # OUTPUTS:
  # c(mi,ma) for symmetric axis with extreme values determining bidirectional limits
  return(c(-max(abs(x)),max(abs(x))))
}