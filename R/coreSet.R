#-----------------------------------------------------------------------------#
# Auxiliary functions to calculate Core Set                                   #
# Cooperation in linear production programming problems                       #
#-----------------------------------------------------------------------------#

#' @import kappalab

# getImputationSet ------------------------------------------------------------
# Imputation Set
getImputationSet <- function(characteristicFunction) {
  
  v <- characteristicFunction
  N <- length(v)
  
  i1 <- c(v[N] - v[2] - v[3], v[2], v[3])
  i2 <- c(v[1], v[N] - v[1] - v[3], v[3])
  i3 <- c(v[1], v[2], v[N] - v[1] - v[2])
  
  imputationSet <- matrix(c(i1, i2, i3), ncol = 3, byrow = TRUE)
  
  return(imputationSet)
  
}


# point3Dto2D -----------------------------------------------------------------
# Point transformation to two dimensions
point3Dto2D <- function(x, alpha, beta) {
  
  coord_x <- x[2]*cos(beta*(pi/180)) - x[1]*cos(alpha*(pi/180))
  coord_y <- x[3] + x[1]*sin(alpha*(pi/180)) + x[2]*sin(beta*(pi/180))
  
  p <- c(coord_x, coord_y)
  names(p) <- c("V1", "V2")

  return(p)
  
}

# lineEq3D --------------------------------------------------------------------
# Equations for lines
lineEq3D <- function(A, B) {
  
  dirVec <- B-A
  
  rhsMat <- matrix(c(A[1], A[2], A[3]))
  
  return(list(coef = dirVec, rhs = rhsMat))
  
}

# getCoreSet ------------------------------------------------------------------
# Core Set for a given game
getCoreSet <- function(characteristicFunction) {
  
  v <- characteristicFunction  
  N <- length(v)
  
  # Imputation Set
  imputationSet <- getImputationSet(v)
  
  # Cut points
  x12 <- c(v[N] - v[6], v[2], v[6] - v[2])
  x13 <- c(v[N] - v[6], v[6] - v[3], v[3])
  x21 <- c(v[1], v[N] - v[5], v[5] - v[1])
  x23 <- c(v[5] - v[3], v[N] - v[5], v[3])
  x31 <- c(v[1], v[4] - v[1], v[N] - v[4])
  x32 <- c(v[4] - v[2], v[2], v[N] - v[4])
  
  
  coreSet <- matrix(c(x12, x13, x21, x23, x31, x32), ncol = 3, byrow = TRUE)
  
  # Intersection Points between lines
  x12 <- linesIntersec(coreSet, 1,2,3,4)
  x13 <- linesIntersec(coreSet, 1,2,5,6)
  x23 <- linesIntersec(coreSet, 3,4,5,6)
  interSet <- matrix(c(x12, x13, x23), ncol = 3, byrow = TRUE)
  
  # Core Set
  core_aux <- rbind(imputationSet, coreSet, interSet)
  core_aux <- unique(core_aux)
  
  coreSet <- matrix(ncol = 3)[-1, ]
  for (i in 1:nrow(core_aux)) {
    point_aux <- core_aux[i, ]
    # Comprobar que queda entre los limites
    l1 <- point_aux[1] >= v[1] & point_aux[1] <= (v[N] - v[6])
    l2 <- point_aux[2] >= v[2] & point_aux[2] <= (v[N] - v[5])
    l3 <- point_aux[3] >= v[3] & point_aux[3] <= (v[N] - v[4])

    if (all(l1, l2, l3)) {
      coreSet <- rbind(coreSet, point_aux)
    }
  }
  
  row.names(coreSet) <- NULL
  if (nrow(coreSet) == 0) {
    stop("Empty core")
  }
  
  return(coreSet)  
  
}

# linesIntersec ---------------------------------------------------------------
# Intersectio of lines
linesIntersec <- function(coreSet, l11, l12, l21, l22){
  
  A <- coreSet[l11,]
  AB <- coreSet[l12,]- coreSet[l11,]
  
  C <- coreSet[l21,]
  CD <- coreSet[l22,]- coreSet[l21,]
  
  if (!(all(AB == 0) || all(abs(CD) < 1e-10))){
    # Lines equations
    eqLine1 <- cbind(A,AB)
    eqLine2 <- cbind(C,CD)
    
    # Equation system
    eqSystem <- cbind(eqLine1[,2], -eqLine2[,2], eqLine2[,1]-eqLine1[,1])
    
    eqSol <- solve(eqSystem[1:2,1:2], eqSystem[1:2,3])
    
    point1 <- eqLine1[,1] + eqLine1[,2]*eqSol[1]
    point2 <- eqLine2[,1] + eqLine2[,2]*eqSol[2]
    
    if (all.equal(point1, point2)) {
      intersecPoint <- point1

    } else {
      cat("Do not intersect")
    }
  }else{
    intersecPoint <- NULL
  }
  
  return(intersecPoint)
  
}

# shapleyValue -------------------------------------------------------------
shapleyValueInt <- function (x){
  
  z <- as.matrix(x$Lex)
  z <- as.vector(z)
  coalitions <- kappalab::set.func(c(0, z))
  SV <- kappalab::Shapley.value(coalitions)
  SV <- as.matrix(SV)
  
  colnames(SV) <- "Shapley Value"
  Output <- list(SV = SV)
  class(Output) <- "ShapleyValue"
  
  return(Output)
}

