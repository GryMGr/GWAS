#' Volcano plot 
#'
#' Plot the effect sizes against the p-values for diagnostics of GWAS results, 
#' or QC when comparing multiple cohorts. 
#' 
#'
#' @param x A result object with 3+ columns. 1st col is coefficients, 2nd is 
#' standard error and 3rd is P-value.  
#' @param file Filename if plot is to be save as a JPEG file.  NOTE file 
#' extention  should not be included in the filename.
#' @param ncol Number of columns if the volcano plot spans a matrix layout if 
#' multiple exposure models is investigated. 
#' @param width Width when saving plot
#' @param height Height when saving plot
#' @param free If the axis should be en sync or not between the plots. 
#' 
#' @return returns a ggplot2 object, piped to X11 
#' @export
#' @import ggplot2
volcano <- function(x, main = "", file = NULL, ncol = 2, width = NULL, height = NULL, free = "free"){

    if(class(x) == "list"){
        x           <- lapply(x,na.omit)
        nn          <- names(x)
        dat			<- sapply(x function(x)data.frame(Effect = x[,1], LogP = -log10(x[,3]), Significant = .significant(x))) 
        ymax        <- sapply(dat, function(x)max(x$Pval, na.rm = TRUE))
        	if(free != "free") ymax <- max(ymax)
		xmax		<- sapply(dat,function(x)range(x$Effect)) 
		dat         <- do.call(rbind,dat)
			if(free != "free") xmax <- range(dat$Effect)
        dat$model   <- factor(rep(nn,sapply(x,nrow)),levels = nn)
        text        <- data.frame(model = nn,  x = 0, y = ymax - ymax/10, text = LETTERS[1:length(nn)])
        wrap       <- TRUE
    }else{
        x       <- na.omit(x)
        wrap    <- FALSE
    	dat		<- data.frame(Effect = x[,1], LogP = -log10(x[,3]), Significant = .significant(x[,3]))
		
	}
    myColours <- c("#D55E00","#E69F00","#56B4E9"); names(myColours) <- c("Bonferroni", "FDR", "NonSign")

    gg <- ggplot(dat, aes(x = Effect, y = LogP)) +
         scale_colour_manual(name = "Significant",values = myColours) +
         labs(title = main)

    if(any(dat$Significant %in% c("FDR", "Bonferroni"))){
         gg <- gg +  geom_point(aes(colour = Significant))
    }else{
        gg  <- gg + geom_point(colour = "#56B4E9")
    }

    if(wrap){
        gg <-   gg +
                facet_wrap(~model, ncol = ncol, scales = free) +
                geom_text(data = text, aes(x=x, y=y,label=text), hjust=0, size=10, show.legend=FALSE)
    }else if(!wrap & overlay){
        gg <-   gg + scale_x_continuous(limits = c(0, max(dat$Observed)))
    }

    if(!is.null(file)){
        if(length(grep("jpeg", file)) == 0) file <- paste(file,"jpeg", sep=".")
         qq <-    gg +
            ggsave(file= file, width = width, height = height, units = "mm", pointsize = 2, dpi = 900)
         }
    gg
}


.significant <- function(x){

	p.adj <- ifelse(p.adjust(x, method = "bonferroni") < 0.05, "Bonferroni", ifelse(p.adjust(x, method = "BH") < 0.05, "FDR", "NonSign"))
	return(factor(p.adj, levels = c("Bonferroni", "FDR", "NonSign")))
}

