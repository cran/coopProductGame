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
#' @param b matrix containing the amount of resources of the several players 
#' where each row is one player.
#' @param plot logical value indicating  if the function displays graphical 
#' solution (\code{TRUE}) or not (\code{FALSE}). Note that this option only make 
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
#' c <- c(68,52)
#' # Production matrix
#' A <- matrix(c(4,5,6,2),ncol=2, byrow = TRUE)
#' # Matrix of resources. Each row is the vector of resources of each player
#' b <- matrix(c(4,33,6,39,60,0),ncol=2, byrow = TRUE)
#' # Solution of the associated linear production game
#' coopProductGame(c,A,b, show.data = TRUE)
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


coopProductGame <- function(c, A, b, plot=FALSE, show.data=FALSE){
  
  coa <- coalitions(dim(b)[1])
  coa1 <- coa$Binary
  coa2 <- coa$Usual
  
  juego <- data.frame(0)
  colnames(juego) <- c("S={0}")
  barx <- data.frame()
  plots <- list()
  for (i in 2:dim(coa1)[1]){
    pos<-which(coa1[i,]!=0)
    
    if(length(pos)==1){
      bs <- b[pos,]
    }else{
      bs <- apply(b[pos,],2,sum)
    }
    
    prod <- makeLP(c, A, bs)
    juego[i] <- get.objective(prod)
    colnames(juego)[i] <- paste0("S={", paste0(strsplit(as.character(coa2[i]), split = "")[[1]], collapse = ","),"}")
    barx <- rbind(barx, get.variables(prod))
    rownames(barx)[i-1] <- paste0("S={", paste0(strsplit(as.character(coa2[i]), split = "")[[1]], collapse = ","),"}")
    
    if (plot==T){
      plots[[i-1]] <- plotlm(prod, A, bs, c, coa2[i])
    }
    
    if(length(pos)== dim(b)[1]){
      
      if(nrow(A)>2){
        dual <- get.dual.solution(prod)[2:(dim(b)[2]+1)]
        Owen <- t(b %*% dual)
      }else{
        dim_aux <- dim(A)[2]
        
        A_dual <- cbind(t(A), -diag(dim_aux))
        b_dual <- c
        c_dual <- c(bs, rep(0,dim_aux))
        
        m <- dim(A_dual)[1]
        n <- dim(A_dual)[2]
        
        bs <- b_dual
        c <- c(c_dual, rep(0,m))
        A <- cbind(A_dual, diag(m))
        init = c(rep(0, n), b_dual)
        basic = n + (1:m)
        
        n1 <- n
        
        out1 <- getDualSolution(a = c, A = A, b = bs, init = init, basic = basic, n1=n1)
        
        Owen <- data.frame()
        
        for(i in 1:nrow(out1)){
          Owen <- rbind(Owen, t(b %*% as.numeric(out1[i,])))
        }
      }
    }
  }
  
  if(plot==TRUE){
    numPlots = length(plots)
    layout <- matrix(seq(1, 2 * ceiling(numPlots/2)),
                     ncol = 2, nrow = ceiling(numPlots/2), byrow = TRUE)
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
  
  rownames(juego) <- c("Associated game")
  colnames(barx) <- rep("", ncol(barx))
  
  if(show.data){
    
    auxOwen <- as.data.frame(Owen)
    colnames(auxOwen) <- c()
    rownames(auxOwen) <- c()
    
    cat(" ------------------------------------------------------------------------\n")
    cat(" Optimal solution of the problem for each coalition: \n")
    cat(" ------------------------------------------------------------------------\n")
    print(round(barx,2))
    cat("\n")
    cat(" ------------------------------------------------------------------------\n")
    cat(" Cooperative production game: \n")
    cat(" ------------------------------------------------------------------------\n")
    print(juego)
    cat(" ------------------------------------------------------------------------\n")
    cat("\n")
    if(nrow(auxOwen) == 1){
      cat(" ------------------------------------------------------------------------\n")
      owenPrint <- paste0("(", auxOwen[1], ", ", auxOwen[2], ", ", auxOwen[3], ")")
      cat("The game has a unique Owen's allocation:\n")
      cat(" ------------------------------------------------------------------------\n")
      print(owenPrint)
    }
    
    if(nrow(auxOwen) == 2){
      cat(" ------------------------------------------------------------------------\n")
      cat("The game has multiple Owen's allocations \n")
      cat("All the points betwen the following are allocations of Owen Set: \n")
      cat(" ------------------------------------------------------------------------\n")
      owenPrint1 <- paste0("(", auxOwen[1,1], ", ", auxOwen[1,2], ", ", auxOwen[1,3], ")")
      owenPrint2 <- paste0("(", auxOwen[2,1], ", ", auxOwen[2,2], ", ", auxOwen[2,3], ")")
      print(owenPrint1)
      print(owenPrint2)
    }
    cat(" ------------------------------------------------------------------------\n")
  }
  
  
  
  invisible(list(Sol = round(barx,2), Valor_obj = as.numeric(juego), Owen = Owen))
  
  
}

#-------------------------------------------------------------------------------------------------#