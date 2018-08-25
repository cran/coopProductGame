#-----------------------------------------------------------------------------#
# CoopProductGame package                                                     #
# Cooperation in linear production programming problems                       #
#-----------------------------------------------------------------------------#

# plotCoreSet -----------------------------------------------------------------
#' Plot Core Set for cooperative production linear games.
#' 
#' Given a linear production game, the 
#' \code{plotCoreSet} function plots the imputation Set, Core Set and the most common
#' solutions (Nucleolus, Shapley Value and allocations 
#' of the Owen Set).
#' 
#' 
#' 
#' In most cases the Owen Set consists of a single allocation, but in some cases there are infinities.
#' In the case that there are infinite allocations, if the problem has two dimensions, 
#' they will be given by a line, which we will represent graphically.
#' If the problem has more than two dimensions, an allocation of all possible ones will be represented.
#' 
#' @param c vector containing the benefits of the products.
#' @param A production matrix.
#' @param B matrix containing the amount of resources of the several players 
#' where each row is one player.
#'
#' @return \code{plotCoreSet} returns a \code{ggplot} object with the imputation set
#' of the game, the core and the most common solutions.
#'  
#' @author D. Prieto
#'  
#' @seealso \code{\link{coopProductGame}}
#' 
#' 
#' @examples
#' # Vector of benefits
#' c <- c(68, 52)
#' # Production matrix
#' A <- matrix(c(4, 5, 6, 2), ncol = 2, byrow = TRUE)
#' # Matrix of resources. Each row is the vector of resources of each player
#' B <- matrix(c(4, 6, 60, 33, 39, 0), ncol = 3, byrow = TRUE)
#' # Solution of the associated linear production game
#' plotCoreSet(c, A, B)
#'   
#'   
#' @export
#' @importFrom stats dist
#' @importFrom utils capture.output
#' @importFrom grDevices chull

plotCoreSet <- function(c, A, B) {
  
  aux <- linearProductionGame(c, A, B)
  game <- round(as.numeric(aux$game)[-1],4)
  owenPoint <- owenSet(c, A, B)$Owen
  
  V1 = NULL
  V2= NULL
  imputationSet <- getImputationSet(game)
  coreSet <- getCoreSet(game)
  
  # Imputation Set
  imputationSet2D <- matrix(ncol = 2)[-1, ]
  for (i in 1:dim(imputationSet)[1]) {
    imputationSet2D <- rbind(imputationSet2D, point3Dto2D(imputationSet[i, ], 30, 30))
  }
  imputationSet2D <- as.data.frame(imputationSet2D)
  imputationSet2D <- imputationSet2D%>%
    dplyr::mutate(posh = 0.5,
                  posv = ifelse(V2 == max(V2), -1, 2))
  
  vertexLabels <- paste0( "(", round(imputationSet[,1],2), ",", round(imputationSet[,2],2), ",", round(imputationSet[,3],2), ")")
  rownames(imputationSet2D) <- vertexLabels
  
  gameDef <- GameTheory::DefineGame(3, game)
  invisible(capture.output(nucleolo <- suppressWarnings(GameTheory::Nucleolus(gameDef)$nucleolus)))
  nucleoulus <- round(nucleolo[,"x(S)"],2)
  nucleoulus2D <- data.frame(t(point3Dto2D(nucleoulus, 30, 30)),
                             textG = paste0( "(", nucleoulus[1], ",", nucleoulus[2], ",", nucleoulus[3], ")"))
  shapleyValue <- round(shapleyValueInt(gameDef)$SV,2)
  shapley2D <-  data.frame(t(point3Dto2D(shapleyValue, 30, 30)),
                           textG = paste0( "(", shapleyValue[1,], ",", shapleyValue[2,], ",", shapleyValue[3,], ")"))
  
  coreset2d <- matrix(ncol = 2)[-1, ]
  for (i in 1:nrow(coreSet)) {
    coreset2d <- rbind(coreset2d, point3Dto2D(coreSet[i, ], 30, 30))
  }
  
  
  if(nrow(owenPoint) == 1){
    owen2D <- data.frame(t(point3Dto2D(owenPoint, 30, 30)),
                         textG = paste0( "(", round(owenPoint[1],2), ",", round(owenPoint[2],2), ",", round(owenPoint[3],2), ")"))
  }else if(nrow(owenPoint) == 2){
    owen2D <- matrix(ncol=2, nrow = 2)
    namesOwen <- matrix(ncol = 1, nrow = 2)
    for(i in 1:nrow(owenPoint)){
      aux <- as.numeric(t(point3Dto2D(owenPoint[i,], 30, 30)))
      textO = paste0( "(", round(owenPoint[i,1],2), ",", round(owenPoint[i,2],2), ",", round(owenPoint[i,3],2), ")")
      owen2D[i,] <-  aux
      namesOwen[i,] <- textO
      
    }
    
    owen2D <- as.data.frame(cbind(owen2D, namesOwen))
  }
  
  colnames(owen2D) <- c("V1", "V2", "textG")
  
  solutions <- rbind(nucleoulus2D, shapley2D, owen2D)
  solutions <- as.data.frame(solutions, check.names = FALSE)
  solutions$names <- c("Nucleolus", "Shapley", rep("Owen Allocation", dim(owen2D)[1]))
  
  solutions2 <- solutions[,1:2]
  
  if(any(dist(solutions2)<10)){
    
    pos <- which(dist(solutions2)<10)
    
    if(any(pos == 1)){
      aj = c(-0.5 , 1.5, -0.5)
    }else{
      aj = c(-0.5,-0.5,1.5)
    }
    
  }else{
    aj = rep(-0.5, dim(solutions2)[1])
  }
  
  SolutionPoints <- solutions$names
  
  margin <- (max(imputationSet2D$V1) - min(imputationSet2D$V1))/10
  margin2 <- (max(imputationSet2D$V2) - min(imputationSet2D$V2))/10
  
  core <- ggplot2::ggplot()+
    ggplot2::geom_polygon(data=imputationSet2D, ggplot2::aes(x=V1,y=V2),fill="grey90", colour = "black", size = 1)+
    ggplot2::geom_point(data=imputationSet2D, ggplot2::aes(x=V1,y=V2),colour="black")+
    ggplot2::geom_text(data = imputationSet2D, ggplot2::aes(x=V1,y=V2), label = rownames(imputationSet2D),hjust=imputationSet2D$posh, vjust=imputationSet2D$posv, size = 4) +
    ggplot2::geom_polygon(data = as.data.frame(coreset2d[chull(coreset2d),]), ggplot2::aes(x=V1,y=V2), fill = "grey68", colour = "black", size = 1)+
    ggplot2::scale_x_continuous(limits = c(min(imputationSet2D$V1)-margin, max(imputationSet2D$V1)+margin))+
    ggplot2::scale_y_continuous(limits = c(min(imputationSet2D$V2)-margin2, max(imputationSet2D$V2)+margin2))
  
  if(nrow(solutions)>3){
    owenSegment <- solutions[which(solutions$names == "Owen Allocation"),]
    otherSolutions <- solutions[which(solutions$names != "Owen Allocation"),]
    SolutionPointsNames <- otherSolutions$names
    x1 = as.numeric(owenSegment$V1[1])
    x2 = as.numeric(owenSegment$V1[2])
    y1 = as.numeric(owenSegment$V2[1])
    y2 = as.numeric(owenSegment$V2[2])
    segment <- data.frame(x1, x2, y1, y2)
    core <- core +
      ggplot2::geom_segment(data = segment, aes(x = x1, y = y1, xend = x2, yend = y2, colour = "Owen Allocation"))+
      ggplot2::geom_point(data = owenSegment, ggplot2::aes(x=as.numeric(V1),y=as.numeric(V2), colour = "Owen Allocation"))+
      ggplot2::geom_text(data = owenSegment, ggplot2::aes(x=as.numeric(V1),y=as.numeric(V2), colour = "Owen Allocation"), label = owenSegment$textG, show.legend = FALSE, vjust = aj[1:2], size = 4)+
      ggplot2::geom_point(data = otherSolutions, ggplot2::aes(x=as.numeric(V1),y=as.numeric(V2), colour = SolutionPointsNames))+
      guides(color=guide_legend(override.aes=list(shape=c(16,NA,16),linetype=c(0,1,0)), title = ""))+
      ggplot2::geom_text(data = otherSolutions, ggplot2::aes(x=as.numeric(V1),y=as.numeric(V2), color = SolutionPointsNames), label = otherSolutions$textG, show.legend = FALSE, vjust = aj[1:2], size = 4)
    
    
  }else{
    
    core <- core +
      ggplot2::geom_point(data = solutions, ggplot2::aes(x=as.numeric(V1),y=as.numeric(V2), colour = SolutionPoints))+
      ggplot2::geom_text(data = solutions, ggplot2::aes(x=as.numeric(V1),y=as.numeric(V2), colour = SolutionPoints), label = solutions$textG,show.legend = FALSE, vjust = aj, size = 4)
  }
  
  
  core <- core +
    ggplot2::theme_bw()+
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), legend.position="bottom",
                   axis.title.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank())+
    ggplot2::ggtitle("Imputation Set and Core Solutions")
  
  return(core)
  
}
