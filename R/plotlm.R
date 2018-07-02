#-----------------------------------------------------------------------------#
# CoopProductGame package                                                     #
# Cooperation in linear production programming problems                       #
#-----------------------------------------------------------------------------#

# plotlm ----------------------------------------------------------------------
#' Plot method for linear production programming problems
#' 
#' This function plot the graphical solution of simple linear production 
#' programming problem \code{A \%*\% x <= b} with two decision 
#' variables. The decision variables must be real, nonnegative and must not 
#' have a finite upper bound. Only inequality constraints are supported.
#' 
#' @param prod a linear production programming problem of class \code{lpExtPtr}.
#' @param A production matrix.
#' @param b vector of resources.
#' @param c vector of benefits.
#' @param title title of plot. By default is \code{NULL}, so it returns a plot
#' without title.
#'
#' @return Returns and plot a \code{ggplot} object with graphical solution of 
#' the problem.
#' 
#' @author D. Prieto
#' 
#' @examples
#' # Vector of benefits
#' c <- c(68,52)
#' # Matrix of coefficients
#' A <- matrix(c(4,5,6,2), ncol = 2, byrow = TRUE)
#' # Vector of resources
#' b <- c(4,33)
#' # Make the associated linear program 
#' prod <- makeLP(c, A, b)
#' plotlm(prod, A, b, c)
#' 
#' @seealso \code{\link{makeLP}}.
#'   
#'   
#' @export
#' @importFrom grDevices contourLines
#' @importFrom graphics par

plotlm <- function(prod, A, b, c, title = NULL){
  
  V1 <- NULL
  V2 <- NULL
  
  # Problem dimensions
  m <- dim(prod)[1]
  n <- dim(prod)[2]
  
  # Test to verify that the problem can be represented
  if (n != 2 || 
      any(get.bounds(prod)$lower != 0) || 
      any(get.bounds(prod)$upper != Inf) || 
      any(get.type(prod) != "real")) 
    stop("The problem can not be represented")
  
  # Feasible Set
  cutPointsAxis <- getCutPointsAxis(A,b)
  feasibleSet <- getFeasibleSet(A, b, cutPointsAxis)
  limits <- round(max(feasibleSet) * c(-0.1, 1.4),3)

  box <- max(feasibleSet) * c(-0.05, 1.2)
  old.par <- par(mar = c(3, 2, 4, 2) + 0.1, pty = "s")
  on.exit(par(old.par))
  
  feasibleSet <- as.data.frame(feasibleSet)
  if(dim(feasibleSet)[2]==1){
      a <- ggplot(feasibleSet)+
        scale_x_continuous("",breaks = seq(0,limits[2]+1),limits = c(0,limits[2]+1) )+
        scale_y_continuous("",breaks = seq(0,limits[2]+1),limits = c(0,limits[2]+1) )+
        geom_segment(x=0,y=limits[1],xend=0,yend=limits[2] + 1,arrow = arrow(length = unit(0.5, "cm")))+
        geom_segment(x=limits[1],y=0,xend=limits[2] + 1,yend=0,arrow = arrow(length = unit(0.5, "cm")))
    }else{
      
      feasibleSet$names <- paste0("(", feasibleSet$V1, ", ", feasibleSet$V2, ")")

      a <- ggplot(feasibleSet)+
        geom_polygon(aes(x=feasibleSet$V1,y=feasibleSet$V2),fill="lightgreen")+
        scale_x_continuous("",breaks = unique(round(seq(0,limits[2], by = ifelse(round(limits[2]/5) ==0, 1, round(limits[2]/5))))),limits = c(limits[1],limits[2]) )+
        scale_y_continuous("",breaks = unique(round(seq(0,limits[2], by = ifelse(round(limits[2]/5) ==0, 1, round(limits[2]/5))))),limits = c(limits[1],limits[2]) )+
        geom_segment(x=0,y=limits[1],xend=0,yend=limits[2],arrow = arrow(length = unit(0.5, "cm")))+
        geom_segment(x=limits[1],y=0,xend=limits[2],yend=0,arrow = arrow(length = unit(0.5, "cm")))
    }
  
  if(!is.null(title)){
    if(!is.character(title)){
      title <- as.character(title)
    }
    if(title == "123"){
      title = "N"
    }
    a <- a + ggtitle(paste("S = {", paste(strsplit(title, split = "")[[1]], collapse = ", "),"}")) 
  }
  
  a <- a + ggplot2::theme_bw() +
    theme(plot.title = ggplot2::element_text(hjust = 0.5), panel.grid.major = ggplot2::element_blank(), 
                 panel.grid.minor = ggplot2::element_blank()) 
    
  
  # Restrictions of the problem
  restrictions <- getRestrictions(A, b, box)
    
  for(i in 1:m){
    a <- a + geom_segment(x = restrictions[i,1], y = restrictions[i,2], xend = restrictions[i,3], yend = restrictions[i,4])
  }

  # Solution of the problem
  solucion <- get.variables(prod)
  a <- a + geom_point(aes(x = solucion[1], y = solucion[2]), colour = "red", size = 3)+ 
    geom_text(aes(x = solucion[1], y = solucion[2],label=paste0("(", round(solucion[1],2),", ", round(solucion[2],2),")")),hjust=-0.2, vjust=-0.2, colour = "red")
  
  if(dim(feasibleSet)[2]!=1){
    
    sol <- get.objective(prod)
    xs<-seq(limits[1],limits[2],0.1)
    ys<-xs
    f<-function(x,y)c[1]*x+c[2]*y - sol
    zs<-outer(xs,ys,FUN=f)

    dat <- contourLines(xs, ys, zs, levels = c((-(limits[2]-limits[1])*2), 0, (limits[2]-limits[1])*2))
    
    solid.x <- c()
    solid.y <- c()
    for(i in 1:length(dat)){
      if(dat[[i]]$level == 0){
        solid.x = c(solid.x, dat[[i]]$x)
        solid.y = c(solid.y, dat[[i]]$y)
      }
    }
    solid.x <- unique(solid.x)
    solid.y <- unique(solid.y)
    
    dash.min <- data.frame()
    for(i in 1:length(dat)){
      if(dat[[i]]$level < 0){
        aux = cbind(dat[[i]]$x, dat[[i]]$y)
        dash.min <- unique(rbind(dash.min, aux))
      }
    }
    
    dash.max <- data.frame()
    for(i in 1:length(dat)){
      if(dat[[i]]$level > 0){
        aux = cbind(dat[[i]]$x, dat[[i]]$y)
        dash.max <- unique(rbind(dash.max, aux))
      }
    }
    
    a <- a + geom_line(data = dash.min, 
                       mapping = aes(x = V1, y = V2), colour = "red", linetype="dashed")+
      geom_line(data = data.frame(a = solid.x, b = solid.y), 
                mapping = aes(x = a, y = b), colour = "red", linetype="solid", size = 1)+
      geom_line(data = dash.max, 
                mapping = aes(x = V1, y = V2), colour = "red",linetype="dashed") 
  }
    
  return(a)
}

#-----------------------------------------------------------------------------#