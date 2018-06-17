# FUNCTION
# stability selection using glmnet lasso with arbitrary subdivisions
# x = response vector
# m = predictor matrix
# PFER = per family error rate to control (if q is not set, it is chosen to achieve an error below this bound)
# L = number of subsample partitions
# t = threshold --> for the bound to hold, this must be chosen such that t > q/D!
# B = number of iterations
# q = number of covariates selected for each subsample (overrides PFER)
rs.stabsel = function( x, m, PFER=1, L=2, t=0.7, B=50, q=NULL, verbose=T ) {
  require(glmnet)
  
    D = ncol(m)
  
  # workaround since I use the colnames later -- assign them if there are none
  if(is.null(colnames(m))) { colnames(m) = paste0("V", 1:ncol(m)) }
  
  # glmnet needs a matrix
  if (class(m) == "data.frame") { m = as.matrix(m) }
  
  # calculate the error -- implement method to automatically select q for a given t and E!
  if (is.null(q)) { q = select_q(PFER, L, D, t); print(q) }
  
  
  error = expected_error_rate( L, q, D, t )
  if (verbose) { print(paste0("PFER controlled at ",error) ) }
  n = nrow(m)
  
  # check if condition for error bound holds
  if (t <= q/D ) { print( paste0("Set t to at least ", q/D, " for the error bounds to hold!") ) }
  
  # partitioning into L equal groups
  lvec = as.integer( seq(1, n, length.out = L+1) )
  
  # allocate accumulation vector and set start index
  acc = vector(mode="character", length=q*B*L)
  ai = 0
  
  # start iterations
  for (iter in 1:B) {
    # random order of observations
    part = sample( 1:n, n )
    
    # Lasso on L equal partitions of the observations, selecting the top q covariates by order of appearance in the regularization path
    for (i in 1:L) {
      idx = part[ (lvec[i]+(i!=1)):lvec[i+1] ]
      gl = glmnet( m[idx, ], x[ idx ], alpha=1 )$beta
      
      # collect q selected vars at each iteration
      acc[ (ai+1):(ai+q) ] = names( sort( apply(gl, 1, function(x) which(x != 0)[1] ) )[1:q] )
      ai = ai + q
    }
  }
  
  s = table(acc)
  s = (s / sum(s)) * q
  # note: this only returns frequencies for covariates that have been selected at least once
  
  return(s)
}




# FUNCTION
# stability selection using glmnet lasso with arbitrary subdivisions
# m = data matrix (rows = observations)
# PFER = per family error rate to control (if q is not set, it is chosen to achieve an error below this bound)
# L = number of subsample partitions
# t = threshold --> for the bound to hold, this must be chosen such that t > q/D!
# B = number of iterations
# q = number of covariates selected for each subsample (overrides PFER)
# n2 = distance transformation (TRUE/FALSE)
rs.stabsel.matrix = function( m, PFER=1, L=2, t=0.7, B=50, q=NULL, n2 = F, description = "", return_list=F ) {
  require(glmnet)
  require(foreach)
  require(doParallel)
  t1 = Sys.time()
  
  if (n2) {
    m = scale(m)
    m = apply(m, 2, convert2sq)
  }
  
  D = ncol(m)-1
  
  # workaround since I use the colnames later -- assign them if there are none
  if(is.null(colnames(m))) { colnames(m) = paste0("V", 1:ncol(m)) }
  
  # glmnet needs a matrix
  if (class(m) == "data.frame") { m = as.matrix(m) }
  
  # calculate the error -- implement method to automatically select q for a given t and E!
  if (is.null(q)) { q = select_q(PFER, L, D, t) }
  
  error = expected_error_rate( L, q, D, t )
  n = nrow(m)
  
  # check if condition for error bound holds
  if (t <= q/D ) { print( paste0("Set t to at least ", q/D, " for the error bounds to hold!") ) }
  
  # partitioning into L equal groups
  lvec = as.integer( seq(1, n, length.out = L+1) )
  
  # Generate list of subsample indices in blocks of size N/L
  ssl = list()
  for (b in 1:B) {
    part = sample( 1:n, n )
    for (i in 1:L) {
      ssl = c(ssl, list( part[ (lvec[i]+(i!=1)):lvec[i+1] ] ) )
    }
  }
  #return(ssl)
  
  # initialize accumulation matrix template
  acc = Matrix(0, nrow = ncol(m), ncol=ncol(m), dimnames = list(colnames(m), colnames(m)))
  
  # ITERATE through subsets
  
  if (!return_list) { combine_fun = add } else { combine_fun = c }
  
  accs = foreach(iter = 1:length(ssl), .combine = combine_fun ) %dopar% {
    require(glmnet)
    tmp = acc
    idx = ssl[[iter]]
    
    # ITERATE through covariates
    for (p in 1:ncol(tmp)) {
      gl = glmnet( m[idx, -p], m[ idx, p ], alpha=1 )$beta
      s = names( sort( apply(gl, 1, function(x) which(x != 0)[1] ) )[1:q] )
      s = s[!is.na(s)]
      tmp[s,p] = tmp[s,p] + 1
    }
    return(tmp)
  }
  
  
  # counts to frequencies
  if(!return_list) { accs = accs/(B*L) } else {
    accsnew = list()
    for (i in 1:length(accs)) {
      accsnew[[i]] = Reduce(add, accs[1:i])/i
    }
    accs = accsnew
  }
  
  return(list("adj"=accs,"time"= Sys.time() - t1, "cutoff"=t, "PFER"=PFER, "PFER_achieved"=error, "L"=L, "q"=q, "D"=D, "description"=description) )
}










# FUNCTION
# transform to distance vector
convert2sq = function( x ) {
  x = as.matrix( dist(x, method = "manhattan") )
  xmrow = apply(x,1,mean, na.rm=T)
  xmcol = apply(x,2,mean, na.rm=T)
  xm = mean(x)
  x = sweep( x, 1, xmrow )
  x = sweep( x, 2, xmcol ) 
  x = x + xm
  x = as.vector(x)
  return(x)
}


# FUNCTION
# Calculate PFER bound for a set of parameters (see rs.stabsel) according to Beinrucker et al. 2016, Corollary 2
expected_error_rate = function( L, q, D, t ) {
  kldist = function( p, q ) { p * log( p/q ) + (1-p) * log( (1-p)/(1-q) ) }  # kullback leibler distance
  out = (L * (1 - t) + 1 ) * exp( -L * kldist(t, q/D) )
  return( out )
}


# FUNCTION
# Select number of selected covars q per subsample to control for a given PFER, based on
# number of subsample partitions L, threshold t, and number of covars D
select_q = function( PFER = 1, L, D, t) {
  for (q in 1:D) {
    if ( expected_error_rate(L, q, D, t) > PFER ) { return( max(1, q-1) ) }
  }
}





# FUNCTION
# stability selection using distance partial correlation with arbitrary subdivisions
# m = data matrix (rows = observations)
# PFER = per family error rate to control (if q is not set, it is chosen to achieve an error below this bound)
# L = number of subsample partitions
# t = threshold --> for the bound to hold, this must be chosen such that t > q/D!
# B = number of iterations
# q = number of covariates selected for each subsample (overrides PFER)
rs.stabsel.dcor = function( m, PFER=1, L=2, t=0.7, B=50, q=NULL, max.n2.rows = 4000, distance.transformation = F, stabsel=T, description = "" ) {
  require(foreach)
  require(doParallel)
  require(corpcor)
  require(Matrix)
  t1 = Sys.time()
  
  D = ncol(m)-1
  
  
  if ( (nrow(m)/L)^2 > max.n2.rows && distance.transformation ) {
    L = ceiling( nrow(m) / sqrt(max.n2.rows) )
    warning( paste0("Parameter L overriden by max.n2.rows! Set to ", L ))
  }
  
  # workaround since I use the colnames later -- assign them if there are none
  if(is.null(colnames(m))) { colnames(m) = paste0("V", 1:ncol(m)) }
  
  # glmnet needs a matrix
  if (class(m) == "data.frame") { m = as.matrix(m) }
  
  # calculate the error -- implement method to automatically select q for a given t and E!
  if (is.null(q)) { q = select_q(PFER, L, D, t) }
  
  error = expected_error_rate( L, q, D, t )
  n = nrow(m)
  
  # check if condition for error bound holds
  if (t <= q/D ) { print( paste0("Set t to at least ", q/D, " for the error bounds to hold!") ) }
  
  # partitioning into L equal groups
  lvec = as.integer( seq(1, n, length.out = L+1) )
  
  # Generate list of subsample indices in blocks of size N/L
  ssl = list()
  for (b in 1:B) {
    part = sample( 1:n, n )
    for (i in 1:L) {
      ssl = c(ssl, list( part[ (lvec[i]+(i!=1)):lvec[i+1] ] ) )
    }
  }
  #return(ssl)
  
  # initialize accumulation matrix template
  acc = Matrix(0, nrow = ncol(m), ncol=ncol(m), dimnames = list(colnames(m), colnames(m)))
  
  # ITERATE through subsets
  accs = foreach(iter = 1:length(ssl), .combine = `+`) %dopar% {
    require(corpcor)
    tmp = acc
    idx = ssl[[iter]]
    
    if (distance.transformation) {
      tmp = abs( rs.distancePartialCorrelation( m[idx,] ) )
    } else {
      tmp = abs( pcor.shrink( m[idx,], verbose = F ) )
    }
    #tmp = abs( cor( m[idx,] ) )
    diag(tmp) = 0
    if( stabsel) { tmp = apply( tmp, 2, function(x) { x[ x < sort(x, decreasing = T)[q] ] = 0; return(x) } ) > 0 }
    
    return(tmp)
  }
  
  # counts to frequencies
  accs = Matrix( accs/(B*L), sparse = T )
  
  # do regular full parcor on top
  fullpc = matrix( pcor.shrink( m, verbose = F ), nrow=nrow(m) ) * (accs >= t)
  
  
  return(list("adj"=accs, "adj_pcor"=fullpc, "time"= Sys.time() - t1, "cutoff"=t, "PFER"=PFER, "PFER_achieved"=error, "L"=L, "q"=q, "D"=D, "B"=B, "description"=description) )
}



# FUNCTION
# computes distance partial correlations
rs.distancePartialCorrelation = function( m, n2_only = F ) {
  require(magrittr)
  
  convert2sq = function( x ) {
    x = as.matrix( dist(x, method = "manhattan") )
    xmrow = apply(x,1,mean, na.rm=T)
    xmcol = apply(x,2,mean, na.rm=T)
    xm = mean(x)
    x = sweep( x, 1, xmrow )
    x = sweep( x, 2, xmcol ) 
    x = x + xm
    x = as.vector(x)
    return(x)
  }
  
  
  m = scale(m)
  m2 = apply(m, 2, convert2sq)
  colnames(m2) = colnames(m)
  
  if(!n2_only) {
    require(corpcor)
    m2 = pcor.shrink(m2, verbose = F)
  }
  
  return(m2)
}















