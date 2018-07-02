#-----------------------------------------------------------------------------#
# CoopProductGame package                                                     #
# Cooperation in linear production programming problems                       #
#-----------------------------------------------------------------------------#

# coalitions ------------------------------------------------------------------
#' Coalitions for a given numbers of players \code{n}.
#' 
#' This functions gives all the coalitions, including the empty coalition,
#' for a number of players n.
#' 
#' @param n Number of players.
#'
#' @return A list with the following components:
#' 
#' \item{\code{Binary}}{Matrix where each row is a binary representation of the 
#' coalition.}
#'  
#' \item{\code{Usual}}{Vector with the usual representatios of the coalitions.}
#'     
#' @author D. Prieto
#' 
#' 
#' @examples
#' # Number of players:
#' n <- 3
#' # Associated coalitions:
#' coalitions(n)
#' 
#' # $Binary
#' #       [,1] [,2] [,3]
#' # [1,]    0    0    0
#' # [2,]    1    0    0
#' # [3,]    0    1    0
#' # [4,]    0    0    1
#' # [5,]    1    1    0
#' # [6,]    1    0    1
#' # [7,]    0    1    1
#' # [8,]    1    1    1
#' # 
#' # $Usual
#' # [1]   0   1   2   3  12  13  23 123
#' 
#' @export
#' @importFrom gtools permutations




coalitions <- function(n){

  coaBinary <- gtools::permutations(n = 2, r = n, v = c(0,1), repeats.allowed = TRUE)

  aux1<-c(0)
  for (i in 2:nrow(coaBinary)){
    aux<-which(coaBinary[i,]!=0)
    
    if(length(aux) == 1){
      aux1 <- c(aux1, aux)
    }else{
      au <- c()
      for(i in 1:length(aux)){
        au <- as.numeric(paste0(au, aux[i]))
      }
      aux1 <- c(aux1,au)
    }
  }
  
  coaBinary<-coaBinary[order(aux1),]
  
  coalitions<-sort(aux1)
  sol<-list(coaBinary,coalitions)
  names(sol)<-c("Binary","Usual")
  return(sol)
  
}