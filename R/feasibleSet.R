#-----------------------------------------------------------------------------#
# Feasible Set ant Cut Points                                                 #
#-----------------------------------------------------------------------------#

# getRestrictions -------------------------------------------------------------
# Problem restrictions

getRestrictions <- function(A, b, box){
  
  m <- dim(A)[1]
  
  restrictions <- matrix(ncol=4,nrow=m)
  for (i in 1:m) {
    if (abs(A[i, 1]) < 1e-09) 
      restrictions[i,]<-c(box[1], b[i]/A[i, 2], box[2], b[i]/A[i,2])
    else if (abs(A[i, 2]) < 1e-09) 
      restrictions[i,] <- c(b[i]/A[i, 1], box[1], b[i]/A[i, 1], box[2])
    else {
      y.intercept <- b[i]/A[i, 2]
      if (y.intercept > box[2]) {
        x0 <- (b[i] - A[i, 2] * box[2])/A[i, 1]
        y0 <- box[2]
      }
      else if (y.intercept < box[1]) {
        x0 <- (b[i] - A[i, 2] * box[1])/A[i, 1]
        y0 <- box[1]
      }
      else {
        x0 <- box[1]
        y0 <- (b[i] - A[i, 1] * box[1])/A[i, 2]
      }
      y.intercept <- (b[i] - A[i, 1] * box[2])/A[i, 2]
      if (y.intercept > box[2]) {
        x1 <- (b[i] - A[i, 2] * box[2])/A[i, 1]
        y1 <- box[2]
      }
      else if (y.intercept < box[1]) {
        x1 <- (b[i] - A[i, 2] * box[1])/A[i, 1]
        y1 <- box[1]
      }
      else {
        x1 <- box[2]
        y1 <- (b[i] - A[i, 1] * box[2])/A[i, 2]
      }
      restrictions[i,] <- c(x0, y0, x1, y1)
    }
  }
  
  return(restrictions)
}


# getFeasibleSet --------------------------------------------------------------
# Feasible region of the problem
getFeasibleSet <- function(coeffMat, rhsVec,cutPointsAxis){
  
  feasible <- c()
  for (j in 1:dim(cutPointsAxis)[2]) {
    if (all(coeffMat %*% cutPointsAxis[, j, drop = FALSE] <= rhsVec) && 
        all(cutPointsAxis[, j] >= 0)){
      feasible <- c(feasible, j)
    }
  }
  feasibleSet <- t(cutPointsAxis[, feasible])
  feasibleSet <- feasibleSet[chull(feasibleSet), ]
  
  return(feasibleSet)
  
}

# getCutPointsAxis ------------------------------------------------------------
# Cut points with axis
getCutPointsAxis <- function(coeffMat, rhsVec){
  
  m <- nrow(coeffMat)
  
  aux1 <- rbind(coeffMat, diag(2))
  aux2 <- c(rhsVec, rep(0, 2))
  cutPointsAxis <- matrix(0, 2, 0)
  for (i in 1:(m + 1)) {
    for (j in (i + 1):(m + 2)) {
      H <- aux1[c(i, j), ]
      if (abs(det(H)) > 1e-09) {
        cutPointsAxis <- cbind(cutPointsAxis, solve(H, aux2[c(i, j)]))
      }
    }
  }
  
  return(cutPointsAxis)
}




