#-----------------------------------------------------------------------------#
# CoopProductGame package                                                     #
# Cooperation in linear production programming problems                       #
#-----------------------------------------------------------------------------#

# nucleolo ----------------------------------------------------------------
#' Nucleolus solution
#' 
#' This function computes the nucleolus solution of a game with a maximum of 
#' 4 agents.
#' 
#' @param game a vector that represents the cooperative game.
#' @param show.data logical value indicating if the function displays the 
#' console output (\code{TRUE}) or not (\code{FALSE}). By default the 
#' value is \code{FALSE}.
#'
#' @return \code{nucleolus} returns and prints the Nucleolus Solution
#' of associated cooperative game.
#'  
#' @author D. Prieto
#' 
#' @examples
#' # Cooperative game
#' game <- c(68, 102, 0, 170, 710, 762, 992)
#' # Nucleolus solution
#' nucleolus(game, show.data = TRUE)
#' 
#'  # ------------------------
#'  # Nucleolus Solution
#'  # ------------------------
#'  # [1] "(149, 192, 651)"
#'   
#' @export
#' 
#' 

nucleolus <- function(game, show.data = FALSE){
  
  agents <- log(length(game)+1)/log(2)
  gameDef <- GameTheory::DefineGame(agents, game)
  
  invisible(capture.output(nucleolo <- suppressWarnings(GameTheory::Nucleolus(gameDef)$nucleolus)))
  nucleolo <- nucleolo[,"x(S)"]
  nucleolo <- t(as.matrix(nucleolo))
  colnames(nucleolo) <- seq(1:ncol(nucleolo))
  rownames(nucleolo) <- "Nucleolus solution"
  
  if(show.data == TRUE){
    cat(" ------------------------\n")
    cat(" Nucleolus solution: \n")
    cat(" ------------------------\n")
    print(paste0("(", paste(nucleolo, collapse = ", "), ")"))
  }
  
  invisible(nucleolo)
}

