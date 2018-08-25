#-----------------------------------------------------------------------------#
# CoopProductGame package                                                     #
# Cooperation in linear production programming problems                       #
#-----------------------------------------------------------------------------#

# productLinearProblem --------------------------------------------------------
#' Linear production programming problems
#' 
#' Given a linear production programming problem \code{A \%*\% x <= b}, the
#' \code{productLinearProblem} solves the problem by making use of
#'  \code{lpSolveAPI}.
#' 
#' @param c vector of benefits.
#' @param A production matrix.
#' @param b vector of resources.
#' @param plot logical value indicating  if the function displays graphical 
#' solution (\code{TRUE}) or not (\code{FALSE}). Note that this option only makes 
#' sense when we have a two-dimension problem.
#' @param show.data logical value indicating if the function displays the 
#' console output (\code{TRUE}) or not (\code{FALSE}). By default the 
#' value is \code{TRUE}.
#'
#' @return \code{productLinearProblem} returns and prints a list with the 
#' following components:
#' 
#' \code{ObjetiveValue} Value of the objetive function from a successfully solved 
#' linear production programming problem.
#' 
#' \code{OptimalSolution} Values of the variables from a successfully solved linear 
#' production programming problem.
#'  
#' @author D. Prieto
#' 
#' @examples
#' # Vector of benefits
#' c <- c(68,52)
#' # Production matrix
#' A <- matrix(c(4,5,6,2),ncol=2, byrow = TRUE)
#' # Matrix of resources. Each row is the vector of resources of each player
#' b <- c(4,33)
#' # Solution of the associated linear production game
#' productLinearProblem(c,A,b, show.data = TRUE)
#' 
#' # ------------------------------------------------------------------------
#' # Objetive value: 
#' # ------------------------------------------------------------------------
#' #   [1] "Z = 68"
#' # 
#' # ------------------------------------------------------------------------
#' # Optimal solution: 
#' # ------------------------------------------------------------------------
#' #   [1] 1 0
#' # ------------------------------------------------------------------------
#'   
#' @export




productLinearProblem <- function(c, A, b, plot=FALSE, show.data=FALSE){
  
  prod <- makeLP(c, A, b)
  valor_obj <- get.objective(prod)
  sol <- get.variables(prod)
  

  
  if(plot == TRUE){
    if(ncol(prod)!=2){
      stop("Error. The problem can not not be represented because not 
           have two decision variables.")
    }else{
      print(plotlm(prod, A, b, c))      
    }
  }
  
  if(show.data == TRUE){
    cat(" ------------------------------------------------------------------------\n")
    cat(" Objetive value: \n")
    cat(" ------------------------------------------------------------------------\n")
    print(paste0("Z = ", valor_obj))
    cat("\n")
    cat(" ------------------------------------------------------------------------\n")
    cat(" Optimal solution: \n")
    cat(" ------------------------------------------------------------------------\n")
    print(sol)
    cat(" ------------------------------------------------------------------------\n")
  }
  
  invisible(list(ObjetiveValue=valor_obj,OptimalSolution=sol))
}

