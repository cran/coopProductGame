#-----------------------------------------------------------------------------#
# CoopProductGame package                                                     #
# Cooperation in linear production programming problems                       #
#-----------------------------------------------------------------------------#

# shapleyValue ----------------------------------------------------------------
#' Shapley Value Solution
#' 
#' Calculates the Shapley Value for a N-agent cooperative game.
#' 
#' @param game a vector that represents the cooperative game.
#' @param show.data logical value indicating if the function displays the 
#' console output (\code{TRUE}) or not (\code{FALSE}). By default the 
#' value is \code{FALSE}.
#'
#' @return \code{shapleyValue} returns and prints the Shapley Value
#' of associated cooperative game.
#'  
#' @author D. Prieto
#' 
#' @examples
#' # Cooperative game
#' game <- c(68, 102, 0, 170, 710, 762, 992)
#' # Shapley Value
#' shapleyValue(game, show.data = TRUE)
#' 
#'  # -----------------------------
#'  # Shapley Value Solution: 
#'  # -----------------------------
#'  # [1] "(229, 272, 491)"
#'   
#' @export

shapleyValue <- function(game, show.data = FALSE){
  
  agents <- log(length(game)+1)/log(2)
  gameDef <- GameTheory::DefineGame(agents, game)
  
  shapley <- shapleyValueInt(gameDef)$SV
  shapley <- t(shapley)
  
  if(show.data == TRUE){
    rownames(shapley) <- ""
    cat(" -----------------------------\n")
    cat(" Shapley Value Solution: \n")
    cat(" -----------------------------\n")
    print(paste0("(", paste(shapley, collapse = ", "), ")"))
  }
  
  invisible(shapley)
}
