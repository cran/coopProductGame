#-----------------------------------------------------------------------------#
# CoopProductGame package                                                     #
# Cooperation in linear production programming problems                       #
#-----------------------------------------------------------------------------#

# plotCoreSet -----------------------------------------------------------------
#' Plot Core Set for cooperative production linear games.
#' 
#' Given a linear production game \code{A \%*\% x <= b(S)}, the 
#' \code{plotCoreSet} plots the imputation Set, Core Set and the most common
#' solutions (Nucleoulus, Shapley Value and allocations 
#' of Owen Set.)
#' 
#' 
#' 
#' This function plots all allocations of Owen Set. If we have multiple 
#' allocations, we should pass the ends of the segment that represents all 
#' the possible allocations.
#' 
#' @param game vector that represents a cooperative linear production game.
#' @param owenPoint matrix with allocations of Owen Set.
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
#' # Game
#' game <- c(68,102,0,170,710,762,992)
#' # Owen allocation of the game
#' owenPoint <- matrix(c(230, 282, 480), nrow = 1)
#' # Plot Core set
#' plotCoreSet(game, owenPoint)
#'   
#'   
#' @export
#' @importFrom stats dist
#' @importFrom utils capture.output
#' @importFrom grDevices chull

plotCoreSet <- function(game, owenPoint) {
  
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
    dplyr::mutate(posh = ifelse(V1 == min(V1), 0.3, ifelse(V1 == max(V1), 0.6, 0.5)),
                  posv = ifelse(V2 == max(V2), -0.5, 1.2))
  
  vertexLabels <- paste0( "(", imputationSet[,1], ",", imputationSet[,2], ",", imputationSet[,3], ")")
  rownames(imputationSet2D) <- vertexLabels
  
  gameDef <- GameTheory::DefineGame(3, game)
  invisible(capture.output(nucleolo <- suppressWarnings(GameTheory::Nucleolus(gameDef)$nucleolus)))
  nucleoulus <- nucleolo[,"x(S)"]
  nucleoulus2D <- data.frame(t(point3Dto2D(nucleoulus, 30, 30)),
                             textG = paste0( "(", nucleoulus[1], ",", nucleoulus[2], ",", nucleoulus[3], ")"))
  shapleyValue <- shapleyValue(gameDef)$SV
  shapley2D <-  data.frame(t(point3Dto2D(shapleyValue, 30, 30)),
                           textG = paste0( "(", shapleyValue[1,], ",", shapleyValue[2,], ",", shapleyValue[3,], ")"))
  
  coreset2d <- matrix(ncol = 2)[-1, ]
  for (i in 1:nrow(coreSet)) {
    coreset2d <- rbind(coreset2d, point3Dto2D(coreSet[i, ], 30, 30))
  }
  
  
  if(nrow(owenPoint) == 1){
    owen2D <- data.frame(t(point3Dto2D(owenPoint, 30, 30)),
                         textG = paste0( "(", owenPoint[1], ",", owenPoint[2], ",", owenPoint[3], ")"))
  }else if(nrow(owenPoint) == 2){
    owen2D <- matrix(ncol=2, nrow = 2)
    namesOwen <- matrix(ncol = 1, nrow = 2)
    for(i in 1:nrow(owenPoint)){
      aux <- as.numeric(t(point3Dto2D(owenPoint[i,], 30, 30)))
      textO = paste0( "(", owenPoint[i,1], ",", owenPoint[i,2], ",", owenPoint[i,3], ")")
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
    
    if(pos == 1){
      aj = c(1.5 ,-0.5, -0.5)
    }else{
      aj = c(-0.5,-0.5,1.5)
    }
    
  }else{
    aj = rep(-0.5, dim(solutions2)[1])
  }
  
  SolutionPoints <- solutions$names
  
  core <- ggplot2::ggplot()+
    ggplot2::geom_polygon(data=imputationSet2D, ggplot2::aes(x=V1,y=V2),fill="grey90", colour = "black", size = 1)+
    ggplot2::geom_point(data=imputationSet2D, ggplot2::aes(x=V1,y=V2),colour="black")+
    ggplot2::geom_text(data = imputationSet2D, ggplot2::aes(x=V1,y=V2), label = rownames(imputationSet2D),hjust=imputationSet2D$posh, vjust=imputationSet2D$posv, size = 4) +
    ggplot2::geom_polygon(data = as.data.frame(coreset2d[chull(coreset2d),]), ggplot2::aes(x=V1,y=V2), fill = "grey68", colour = "black", size = 1)
  
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