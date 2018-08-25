#-----------------------------------------------------------------------------#
# CoopProductGame package                                                     #
# Cooperation in linear production programming problems                       #
#-----------------------------------------------------------------------------#

# makeLP ----------------------------------------------------------------------
#' Make a linear production programming problem
#' 
#' Given a linear production problem \code{A \%*\% x <= b}, the
#' \code{makeLP} function creates a new \code{lpSolve} linear program model object.
#' 
#' @param c vector of benefits.
#' @param A production matrix.
#' @param b vector of resources.
#'
#' @return \code{makeLP} returns a \code{lpSolve} linear program model object. 
#' Specifically an R external pointer with class \code{lpExtPtr}.
#' 
#' @examples
#' # Vector of benefits
#' c <- c(68,52)
#' # Production matrix
#' A <- matrix(c(4, 5, 6, 2), ncol = 2, byrow = TRUE)
#' # Vector of resources
#' b <- c(4,33)
#' # Make the associated linear production problem 
#' prod <- makeLP(c, A, b)
#'   
#' @author D. Prieto
#'
#'   
#' @export
 

makeLP <- function(c, A, b){
  prod=make.lp(dim(A)[1],dim(A)[2])
  
  for (j in 1:dim(A)[1]){  
    set.row(prod,j,c(A[j,]))  
  }
  
  set.objfn(prod,c)
  set.constr.type(prod,rep("<=",dim(A)[1]))
  set.rhs(prod,b)
  lp.control(prod,sense="max")
  sol <- solve(prod)
  
  if (sol!=0){
    stop
    print("Error: optimal solution not found.")
  }else{
    return(prod)
  }
}

#-----------------------------------------------------------------------------#