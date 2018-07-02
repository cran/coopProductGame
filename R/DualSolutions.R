#-----------------------------------------------------------------------------#
# Auxiliary function to calculate multiple dual solutions                     #
# Cooperation in linear production programming problems                       #
#-----------------------------------------------------------------------------#

# getDualSolution -------------------------------------------------------------
# Imputation Set
getDualSolution <- function(a,A,b,init,basic, n1=N){
  
  # Total decision variables
  N <- ncol(A)
  # Total restrictions
  M <- nrow(A)
  # Position of nonbasic variables
  nonbasic <- (1:N)[-basic]
  # Simplex table
  simplexTable <- cbind(b, -A[, nonbasic, drop = FALSE])
  obfun <- apply(simplexTable[(M+n1-N+1):M,,drop=FALSE],2L,sum)
  simplexTable <- rbind(c(0,a[nonbasic]),simplexTable,obfun)
  obfun <- obfun[-1L]
  
  # Iterations 
  k <- 1
  while (!all(obfun > -1e-10) && (k <= (M+2*N))){
    pcol <- 1 + order(obfun)[1]
    neg <- 1 + (1:M)[simplexTable[2:(M+1),pcol] < -1e-10]
    ratios <- -simplexTable[neg,1L]/simplexTable[neg,pcol]
    prow <- neg[order(ratios)[1L]]
    simplexTable <- pivotSimplexTable(simplexTable,prow,pcol)
    temp <- basic[prow-1L]
    basic[prow-1L] <- nonbasic[pcol-1L]
    nonbasic[pcol-1L] <- temp
    obfun <- simplexTable[M+2L,-1L]
    k <- k+1
  }
  
  val.aux <- simplexTable[M+2,1L]
  
  if ((val.aux < 1e-10) && any(basic>n1)) {
    ar <- (1L:M)[basic>n1]
    for (j in seq_along(temp)) {
      prow <- 1+ar[j]
      pcol <- 1 + order(
        nonbasic[abs(simplexTable[prow,-1L])>1e-10])[1L]
      simplexTable <- pivotSimplexTable(simplexTable,prow,pcol)
      temp1 <- basic[prow-1L]
      basic[prow-1L] <- nonbasic[pcol-1L]
      nonbasic[pcol-1L] <- temp1
    }
  }
  
  soln <- rep(0,N)
  soln[basic] <- simplexTable[2:(M+1L),1L]
  solution <- data.frame()
  solution <- rbind(solution,soln[seq(1:2)])
  
  while (any(simplexTable[1,][-1][which(nonbasic<=n1)]<1e-10)){
    pcol <- 1+order(obfun)[1L]
    neg <- 1+ (1L:M)[simplexTable[2:(M+1),pcol] < -1e-10]
    ratios <- -simplexTable[neg,1L]/simplexTable[neg,pcol]
    prow <- neg[order(ratios)[1L]]
    if(is.na(prow)){break}
    simplexTable <- pivotSimplexTable(simplexTable,prow,pcol)
    temp <- basic[prow-1L]
    basic[prow-1L] <- nonbasic[pcol-1L]
    nonbasic[pcol-1L] <- temp
    obfun <- simplexTable[M+2L,-1L]
    
    val.aux <- simplexTable[M+2,1L]
    
    if ((val.aux < 1e-10) && any(basic>n1)) {
      ar <- (1L:M)[basic>n1]
      for (j in seq_along(temp)) {
        prow <- 1+ar[j]
        pcol <- 1 + order(
          nonbasic[abs(simplexTable[prow,-1L])>1e-10])[1L]
        simplexTable <- pivotSimplexTable(simplexTable,prow,pcol)
        temp1 <- basic[prow-1L]
        basic[prow-1L] <- nonbasic[pcol-1L]
        nonbasic[pcol-1L] <- temp1
      }
    }
    soln <- rep(0,N)
    soln[basic] <- simplexTable[2:(M+1L),1L]
    
    solution <- rbind(solution,soln[seq(1:2)])
    
    if(dim(solution)[1] > dim(unique(solution))[1]){
      solution <- unique(solution) 
      break
    }
  }

  return(solution)
  
}


# pivotSimplexTable -----------------------------------------------------------
# Pivot simplex table
pivotSimplexTable <- function(simplexTable, pivotRow, pivotCol) {
  
  pivot <- simplexTable[pivotRow, pivotCol]
  pcv <- simplexTable[,pivotCol]
  
  simplexTable[-pivotRow,] <- simplexTable[-pivotRow, ] - (simplexTable[-pivotRow, pivotCol] / pivot) %o% simplexTable[pivotRow, ]
  simplexTable[pivotRow, ] <- simplexTable[pivotRow, ] / (-pivot)
  simplexTable[pivotRow, pivotCol] <- 1 / pivot
  simplexTable[-pivotRow, pivotCol] <- pcv[-pivotRow] / pivot
  
  return(simplexTable)
}