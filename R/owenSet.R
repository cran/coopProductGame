
#-----------------------------------------------------------------------------#
# CoopProductGame package                                                     #
# Cooperation in linear production programming problems                       #
#-----------------------------------------------------------------------------#

# owenSet ----------------------------------------------------------------------
#' Owen Set
#' 
#' This function computes the Owen Set of a linear production game
#' 
#' 
#' @param c vector containing the benefits of the products.
#' @param A production matrix.
#' @param B matrix containing the amount of resources of the several players 
#' where each row is one player.
#' @param show.data logical value indicating if the function displays the 
#' console output (\code{TRUE}) or not (\code{FALSE}). By default the 
#' value is \code{FALSE}.
#'
#' @return \code{owenSet} returns and prints the owen Set of associated 
#' linear production problem.
#'  
#' @author D. Prieto
#' 
#' @examples
#' # Vector of benefits
#' c <- c(68, 52)
#' # Production matrix
#' A <- matrix(c(4, 5, 6, 2), ncol=2, byrow = TRUE)
#' # Matrix of resources. Each row is the vector of resources of each player
#' B <- matrix(c(4, 6, 60, 33, 39, 0),ncol = 3, byrow = TRUE)
#' # Solution of the associated linear production game
#' owenSet(c, A, B, show.data = TRUE)
#' 
#'  # ------------------------------------------------------------------------
#'  # The linear production problem has a unique Owen's allocation:
#'  # ------------------------------------------------------------------------
#'  # [1] "(230, 282, 480)"
#'   
#' @export
#' 
#' 

owenSet <- function(c, A, B, show.data = FALSE){
  
  bN <- apply(B, 1, sum)
  prod <- makeLP(c, A, bN)
  multiple <- 0
  if(nrow(A)>2){
    
    c_dual <- bN
    A_dual <- t(A)
    B_dual <- c
    prod_dual=make.lp(dim(A_dual)[1],dim(A_dual)[2])
    for (j in 1:dim(A_dual)[1]){  
      set.row(prod_dual,j,c(A_dual[j,]))  
    }
    set.objfn(prod_dual,c_dual)
    set.constr.type(prod_dual,rep(">=",dim(A_dual)[1]))
    set.rhs(prod_dual,B_dual)
    sol <- solve(prod_dual)
    
    basics <- get.basis(prod_dual)
    aux <- get.basis(prod_dual, nonbasic = TRUE)
    nonbasics <- -aux[(length(basics)+1):length(aux)]
    
    
    duals <-get.sensitivity.rhs(prod_dual)$duals[nonbasics]
    if(any(duals==0)){
      multiple <- 1
    }
    
    dual <- get.dual.solution(prod)[2:(dim(B)[1]+1)]
    Owen <- round(t(t(B) %*% dual),3)
  }else{
    dim_aux <- dim(A)[2]
    
    A_dual <- cbind(t(A), -diag(dim_aux))
    b_dual <- c
    c_dual <- c(bN, rep(0,dim_aux))
    
    m <- dim(A_dual)[1]
    n <- dim(A_dual)[2]
    
    bN <- b_dual
    c <- c(c_dual, rep(0,m))
    A <- cbind(A_dual, diag(m))
    init = c(rep(0, n), b_dual)
    basic = n + (1:m)
    
    n1 <- n
    
    out1 <- getDualSolution(a = c, A = A, b = bN, init = init, basic = basic, n1=n1)
    
    Owen <- data.frame()
    
    for(i in 1:nrow(out1)){
      Owen <- rbind(Owen, t(t(B) %*% as.numeric(out1[i,])))
    }
    
    if(nrow(Owen) == 2){
      multiple <- 1
    }
    
  }
  
  auxOwen <- as.data.frame((Owen))
  colnames(auxOwen) <- seq(1:ncol(auxOwen))
  rownames(auxOwen) <- c()
  
  if(show.data){
  
    
    if(nrow(auxOwen) == 1 & multiple == 0){
      cat(" ------------------------------------------------------------------------\n")
      owenPrint <- paste0("(", paste(auxOwen, collapse = ", "), ")")
      cat("The linear production problem has a unique Owen's allocation:\n")
      cat(" ------------------------------------------------------------------------\n")
      print(owenPrint)
    }
    
    if(nrow(auxOwen) == 2 || multiple == 1){
      
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
  }
  
  
  invisible(list(Owen = auxOwen, multipleSolutions = multiple))
}
