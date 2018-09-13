#' Circular Manhattan plot                                                     
#'                                                                              
#' Some text here
#' 
#' 
#'                                                                              
#' @param pVal Vector of p-values 
#' @param pos Vector of genomic positions for each p-value
#' @param chr Corresponding vector of chromosomes for p-value
#' @param annot annotation for some of the top hits (not supported yet)
#' @param log True if p-values are on logarithmic form. False as default. 
#' @param filename Name if saving to file 
#' @param width Width when saving plot                                          
#' @param height Height when saving plot                                        
#'                                                                              
#' @return returns a ggplot2 object, piped to X11                               
#' @export                                                                      
#' @import ggplot2          
manhattan.circular <- function(pVal, pos, chr, annot, log = FALSE, filename = NULL, width = 10, height = 10){




gg <- ggplot(dat, aes(x = pos_cum, y = logP, colour = chr)) +                         
        geom_point(size = .4) +                                                 
        coord_polar() +                                                         
        ylim(min(dat$logP) - 3, 0) +                                            
        xlim(0,max(dat$pos_cum)) 

		return(gg)
}


#' Regular Manhattan plot
#' 
#' Some text here
#' 
#' 
#' @param pVal Vector of p-values 
#' @param pos Vector of genomic positions for each p-value
#' @param chr Corresponding vector of chromosomes for p-value
#' @param annot annotation for some of the top hits (not supported yet)
#' @param log True if p-values are on logarithmic form. False as default. 
#' @param filename Name if saving to file 
#' @param width Width when saving plot                                          
#' @param height Height when saving plot                                        
#'                                                                              
#' @return returns a ggplot2 object, piped to X11                               
#' @export                                                                      
#' @import ggplot2          
manhattan.regular <- function(pVal, pos, chr, annot, log = FALSE, filename = NULL, width = 100, height = 10){

gg <- ggplot(dat, aes(x = pos_cum, y = logP, colour = chr)) +                         
        geom_point(size = .4) +                                                 
        ylim(min(dat$logP) - 3, 0) +                                            
        xlim(0,max(dat$pos_cum)) 

		return(gg)
}



.manhattan <- function(pVal, pos, chr, annot = NULL, log = FALSE){
	
#	o		<- order(pos, chr) ## If list is not ordered!! Assume to be ordered for now! 
#	chr		<- chr[o]; pos <- pos[o]
	ll		<- 100 # Add space between the columns 
	pos2	<- split(pos, chr)
	pos_cum	<- numeric(length(pos))
	pos_cum[chr == chr[1]] <- pos2[[1]]









}




