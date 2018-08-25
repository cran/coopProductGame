#-------------------------------------------------------------------------------------------------#
# CoopProductGame package                                                                         #
# Cooperation in linear production programming problems                                           #
#-------------------------------------------------------------------------------------------------#

# linearProductionGame ----------------------------------------------------------------------------
#' Cooperative linear production games
#' 
#' Given a linear production problem, the
#' \code{linearProductionGame} function solves the problem by making use of \code{lpSolveAPI}
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
#' @return \code{linearProductionGame} returns a list with the solutions of 
#' the associated problem of each coalition and the objective value for coalition N.
#'  
#' @author D. Prieto
#' 
#' @examples
#' # Vector of benefits
#' c <- c(68,52)
#' # Production matrix
#' A <- matrix(c(4,5,6,2),ncol=2, byrow = TRUE)
#' # Matrix of resources. Each column is the vector of resources of each player
#' B <- matrix(c(4, 6, 60, 33, 39, 0),ncol = 3, byrow = TRUE)
#' # Solution of the associated linear production game
#' linearProductionGame(c, A, B, show.data = TRUE)
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
#' 
#'   
#' @import dplyr
#' @import GameTheory
#' @import ggplot2
#' @import lpSolveAPI
#' @import grid
#' @export


linearProductionGame <- function(c, A, B, plot=FALSE, show.data=FALSE){
  
  coa <- coalitions(dim(B)[2])
  coa1 <- coa$Binary
  coa2 <- coa$Usual
  
  juego <- data.frame(0)
  colnames(juego) <- c("S={0}")
  barx <- data.frame()
  plots <- list()
  for (i in 2:dim(coa1)[1]){
    pos<-which(coa1[i,]!=0)
    
    if(length(pos)==1){
      bs <- B[,pos]
    }else{
      bs <- apply(B[,pos],1,sum)
    }
    
    prod <- makeLP(c, A, bs)
    juego[i] <- get.objective(prod)
    colnames(juego)[i] <- paste0("S={", paste0(strsplit(as.character(coa2[i]), split = "")[[1]], collapse = ","),"}")
    barx <- rbind(barx, get.variables(prod))
    rownames(barx)[i-1] <- paste0("S={", paste0(strsplit(as.character(coa2[i]), split = "")[[1]], collapse = ","),"}")
    
    if (plot==T){
      plots[[i-1]] <- plotlm(prod, A, bs, c, coa2[i])
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
    
  }
  
  invisible(list(Sol = round(barx,2), game = juego))
  
  
}

#-------------------------------------------------------------------------------------------------#