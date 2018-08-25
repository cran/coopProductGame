#-------------------------------------------------------------------------------------------------#
# CoopProductGame package                                                                         #
# Cooperation in linear production programming problems                                           #
#-------------------------------------------------------------------------------------------------#

# coopProductGame -------------------------------------------------------------------------------
#' Cooperative linear production games
#' 
#' Given a linear production problem \code{A\%*\%x <= B}, the
#' \code{coopProductGame} solves the problem by making use of \code{lpSolveAPI}
#' where each agent provides his own resources. 
#' 
#' @param c vector containing the benefits of the products.
#' @param A production matrix.
#' @param B matrix containing the amount of resources of the several players 
#' where each row is one player.
#' @param plot logical value indicating if the function displays graphical 
#' solution (\code{TRUE}) or not (\code{FALSE}). Note that this option only makes
#' sense when we have a two-dimension problem.
#' @param show.data logical value indicating if the function displays the 
#' console output (\code{TRUE}) or not (\code{FALSE}). By default the 
#' value is \code{TRUE}.
#'
#' @return \code{coopProductGame} returns a list with the solution of the problem,
#'  the objective value and a Owen allocation if it exists. If we have a two 
#'  dimension dual problem, the function returns all the Owen allocations 
#'  (if there are more than one we obtain the end points of the segment 
#'  that contains all possible allocations.)
#'  
#' @author D. Prieto
#' 
#' @examples
#' # Vector of benefits
#' c <- c(68, 52)
#' # Production matrix
#' A <- matrix(c(4, 5, 6, 2), ncol = 2, byrow = TRUE)
#' # Matrix of resources. Each row is the vector of resources of each player
#' B <- matrix(c(4, 6, 60, 33, 39, 0), ncol = 3, byrow = TRUE)
#' # Solution of the associated linear production game
#' coopProductGame(c, A, B, show.data = TRUE)
#' 
#'  # ------------------------------------------------------------------------
#'  # Optimal solution of the problem for each coalition: 
#'  # ------------------------------------------------------------------------
#'  # 
#'  # S={1}      1.00  0.00
#'  # S={2}      1.50  0.00
#'  # S={3}      0.00  0.00
#'  # S={1,2}    2.50  0.00
#'  # S={1,3}    1.68 11.45
#'  # S={2,3}    2.86 10.91
#'  # S={1,2,3} 10.00  6.00
#'  # 
#'  # ------------------------------------------------------------------------
#'  #   Cooperative production game: 
#'  # ------------------------------------------------------------------------
#'  #                S={0} S={1} S={2} S={3} S={1,2} S={1,3} S={2,3} S={1,2,3}
#'  # Associated game    0    68   102     0     170     710     762       992
#'  # ------------------------------------------------------------------------
#'  #   
#'  # ------------------------------------------------------------------------
#'  #   The game has a unique Owen's allocation:
#'  # ------------------------------------------------------------------------
#'  # [1] "(230, 282, 480)"
#'  # ------------------------------------------------------------------------
#' 
#'   
#' @import dplyr
#' @import GameTheory
#' @import ggplot2
#' @import lpSolveAPI
#' @import grid
#' @export


coopProductGame <- function(c, A, B, plot=FALSE, show.data=FALSE){
  
  
  if(plot==TRUE){
    aux <- print(linearProductionGame(c, A, B, plot=TRUE))
  }else{
    aux <- linearProductionGame(c, A, B)
  }
  
  
  game <- aux$game
  sol <- aux$Sol
  
  rownames(game) <- c("Associated game")
  
 
  shapleyValue <- shapleyValue(game[-1])
  nucleolo <- nucleolus(game[-1])
  Owen <- owenSet(c, A, B)
  
  
  if(show.data){
    
    auxOwen <- as.data.frame(Owen$Owen)
    colnames(auxOwen) <- c()
    rownames(auxOwen) <- c()
    
    cat(" ------------------------------------------------------------------------\n")
    cat(" Optimal solution of the problem for each coalition: \n")
    cat(" ------------------------------------------------------------------------\n")
    print(round(sol,2))
    cat("\n")
    cat(" ------------------------------------------------------------------------\n")
    cat(" Cooperative production game: \n")
    cat(" ------------------------------------------------------------------------\n")
    print(game)
    cat(" ------------------------------------------------------------------------\n")
    cat("\n")
    if(nrow(auxOwen) == 1 & Owen$multipleSolutions == 0){
      cat(" ------------------------------------------------------------------------\n")
      owenPrint <- paste0("(", paste(auxOwen, collapse = ", "), ")")
      cat("The linear production problem has a unique Owen's allocation:\n")
      cat(" ------------------------------------------------------------------------\n")
      print(owenPrint)
    }
    
    if(nrow(auxOwen) == 2 || Owen$multipleSolutions == 1){
      
      if(nrow(auxOwen) == 2){
        cat(" ------------------------------------------------------------------------\n")
        cat("The linear production problem has multiple Owen's allocations \n")
        cat("All the points between the following are allocations of the Owen Set: \n")
        cat(" ------------------------------------------------------------------------\n")
        owenPrint1 <- paste0("(", paste(auxOwen[1,], collapse = ", "), ")")
        owenPrint2 <- paste0("(", paste(auxOwen[2,], collapse = ", "), ")")
        print(owenPrint1)
        print(owenPrint2)
      }else{
        cat(" ------------------------------------------------------------------------\n")
        cat("The linear production problem has multiple Owen's allocations \n")
        cat("One of them is the following \n")
        cat(" ------------------------------------------------------------------------\n")
        owenPrint1 <- paste0("(", paste(auxOwen[1,], collapse = ", "), ")")
        print(owenPrint1)
      }
    }
    cat(" ------------------------------------------------------------------------\n")
    cat("\n")
    cat(" ------------------------------------------------------------------------\n")
    cat(" Shapley: \n")
    cat(" ------------------------------------------------------------------------\n")
    shapleyPrint <- paste0("(", paste(round(shapleyValue,3), collapse = ", "), ")")
    print(shapleyPrint)
    cat(" ------------------------------------------------------------------------\n")
    cat("\n")
    
    cat(" ------------------------------------------------------------------------\n")
    cat(" Nucleolus: \n")
    cat(" ------------------------------------------------------------------------\n")
    nucleoloPrint <- paste0("(", paste(round(nucleolo, 3), collapse = ", "), ")")
    print(nucleoloPrint)
    cat(" ------------------------------------------------------------------------\n")
    
  }
  
  colnames(Owen$Owen) <- seq(1:ncol(Owen$Owen))
  
  
  invisible(list(Sol = round(sol,2), game = as.numeric(game), Owen = Owen$Owen, 
                 Shapley = shapleyValue, Nucleolus = nucleolo))
  
  
}

#-------------------------------------------------------------------------------------------------#