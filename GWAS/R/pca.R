#' PCA plot of SNP data  
#' 
#' Plot the PCA or SVA composition of a large data set of SNP. For very large
#' data sets, the SVA decomposition is strongly recommended, due to its 
#' computational efficiency over PCA. 
#'
#' @param x A matrix with SNP data, each column represents a SNP.
#' @param col (Optional) a group factor such as Batch or population sub-structure  
#' @param file Filename if QQ plot is to be save as a FPD/JPEG. NOTE file extension 
#' should be included in the filename.
#' @param width Width when saving plot
#' @param height Height when saving plot
#' @param svd If SVD composition is to be used instead of PCA, recommended for large data sets. Default is false. 
#' @return returns a ggplot2 object, piped to X11 
#' @export
#' @import ggplot2
#' @import irlba::irlba
pca <- function(x, col = "Population", main = "", file = NULL, width = NULL, height = NULL, svd = FALSE){
        
		if(nrow(x) > ncol(x)){
			cat("Number of SNPs less than number of individuals, will transpose the data matrix\n")
			x	<- t(x)
		}

		x       	<- na.omit(x)
		if(svd){
			L 		<- irlba(x, 2)
			comp	<- L$u * (nrow(x) -1)
		}else{
			pca		<- prcomp(x)
			comp	<- pca$x[,1:2]
		}
	

	dat	<- data.frame(PC1 = comp[,1], PC2 = comp[,2], col = col) 

	myColours <- color_pal(length(unique(col))); names(myColours) <- unique(col)

   gg <- ggplot(dat, aes(x = PC1, y = PC2)) +                        
         labs(title = main)
                                                                                
    if(length(unique(col)) > 1){                       
         gg <- gg +  geom_point(aes(colour = col)) +                       
         scale_colour_manual(name = "Population",values = myColours)
    }else{                                                                      
        gg  <- gg + geom_point(colour = "#383838")                              
    }                                                                           
    ## Add title of X-Y axis with SVA dicomposition      
		 
    if(!is.null(file)){
        if(length(grep("jpeg", file)) == 0) file <- paste(file,"jpeg", sep=".")
         qq <-    gg +
            ggsave(file= file, width = width, height = height, units = "mm", pointsize = 2, dpi = 900)
         }
    gg
}

color_pal <- function(n){

	pal <- c("#7FC97F","#BEAED4","#FDC086","#FFFF99","#386CB0","#F0027F","#BF5B17","#666666")
	if(n > 8)
		pal <- rep(pal, ceiling(n/8))

	return(pal[1:n])
}

