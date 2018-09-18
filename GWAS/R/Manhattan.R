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
manhattan.circular <- function(pVal, pos, chr, annot, filename = NULL, width = 10, height = 10){


	sign.th	 <- c(log10(max(1.01,pVal[p.adjust(pVal, method = "fdr") < 0.05])),  log10(max(1.01,pVal[pVal < 0.05/length(pVal)]))) 
	pVal	 <- c(1,pVal); pos <- c(0,pos); chr  <- c(chr[1],chr)
	dat 	 <- .manhattan(pVal, pos, chr)
	dat$logP <- -1*dat$logP
	chr.cent <- sapply(split(dat$pos_cum, chr),median)
	xlim	 <- c(0, max(dat$pos_cum))


gg <- ggplot(dat, aes(x = pos_cum, y = logP, colour = chr)) +                         
        geom_point(size = .4) +                                                 
        coord_polar() +                                                         
        ylim(min(dat$logP) - 3, 0) +                                            
		scale_x_continuous(breaks = chr.cent,
        	labels=seq_along(chr.cent), 
			limits  = xlim) + 
		theme(axis.text.y=element_blank()) +
		geom_hline(yintercept = sign.th, colour = c("blue","red")) + 
		theme(legend.position="none") + 
		scale_fill_brewer(palette="Set1")


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
manhattan.regular <- function(pVal, pos, chr, annot, filename = NULL, width = 100, height = 10, wb = FALSE){


	sign.th	 <- c(-1*log10(max(1.01,pVal[p.adjust(pVal, method = "fdr") < 0.05])),  -1*log10(max(1.01,pVal[pVal < 0.05/length(pVal)]))) 
	dat 	 <- .manhattan(pVal, pos, chr)
	chr.cent <- sapply(split(dat$pos_cum, chr),median)
	xlim	 <- c(0, max(dat$pos_cum))
	nChr	 <- length(unique(chr))
	if(wb){
		col  <- rep(c("#999999","#000000"),nChr/2)
		if(nChr %% 2 == 1)col <- c(col, "#999999")
		dat$col <- as.character(rep(col, table(chr)))	
	}else{
		dat$col <- rep(rep(RColorBrewer::brewer.pal(12,"Set3"),3)[1:nChr], table(chr)) 
	}
	

gg <- ggplot(dat, aes(x = pos_cum, y = logP, group = chr)) +                         
        geom_point(size = .4, aes(colour = col)) +                                                 
        ylim(0,max(dat$logP)) +                                            
		scale_x_continuous(breaks = chr.cent,
        	labels=seq_along(chr.cent), 
			limits  = xlim) + 
		geom_hline(yintercept = sign.th, colour = c("blue","red")) + 
		theme(legend.position="none") 
		
		return(gg)
}

manhattan <- function(...) manhattan.regular(...)


.manhattan <- function(pVal, pos, chr, annot = NULL, log = FALSE){
	
#	o		<- order(pos, chr) ## If list is not ordered!! Assume to be ordered for now! 
#	chr		<- chr[o]; pos <- pos[o]
	
	pos2	<- split(pos, chr)
	pos_cum	<- numeric(length(pos))
	ll		<- median(pos2[[1]]/4) 
	pos_cum[chr == chr[1]] <- pos2[[1]] + ll - min(pos2[[1]])

	for(ch in seq_along(pos2)[-1]){
		pos_cum[chr == ch] <- pos2[[ch]] + max(pos_cum) + ll - min(pos2[[ch]]) 
	}

	return(data.frame(pVal = pVal, logP = -log10(pVal), 
				pos_cum = pos_cum, chr = as.factor(chr)))
				

}









