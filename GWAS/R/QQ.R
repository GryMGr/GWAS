#' QQ plot 
#'
#' Print out the top CpGs from an EWAS model, with three collumns,
#' the effect size, the standard deviation and the p-value for each
#' CpG. The CpGs are sorted by p-value if sort = TRUE. 
#'
#' @param x A result object from the EWAS, of class 
#' @param file Filename if QQ plot is to be save as a FPD/JPEG. NOTE file extention 
#' should not be included in the filename.
#' @param ncol Number of colunms if the QQ plot spans a matrix layout
#' @param overlay If multiple QQ plots is to be overlayed instead of ploted in seperate windows
#' @param width Width when saving plot
#' @param height Height when saving plot
#' @param free If the sub-plot windows should be locked with the same axsis 
#' @return returns a ggplot2 object, piped to X11 
#' @export
#' @import ggplot2
QQ <- function(x, main = "", file = NULL, ncol = 2, width = NULL, height = NULL, overlay = FALSE, free = "free"){

    if(class(x) == "list"){
        x           <- lapply(x,na.omit)
        nn          <- names(x)
        dat         <- lapply(x, .QQdata)
        ymax        <- sapply(dat, function(x)max(x$Observed, na.rm = TRUE))
                if(free != "free") ymax <- max(ymax)
        dat         <- do.call(rbind,dat)
        dat$model   <- factor(rep(nn,sapply(x,nrow)),levels = nn)
        text        <- data.frame(model = nn,  x = 0, y = ymax - ymax/10, text = LETTERS[1:length(nn)])
        wrap       <- TRUE
    }else{
        x       <- na.omit(x)
        dat     <- .QQdata(x)
        wrap    <- FALSE
    }
        dat$Significant <- factor(dat$Significant, levels = c("Bonferroni", "FDR", "NonSign"))
        myColours <- c("#D55E00","#E69F00","#56B4E9"); names(myColours) <- c("Bonferroni", "FDR", "NonSign")

    gg <- ggplot(dat, aes(x = Expected, y = Observed)) +
         scale_colour_manual(name = "Significant",values = myColours) +
         labs(title = main) +
         geom_abline(intercept = 0, slope = 1, size = .3, colour = "darkgray")

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

.QQdata <- function(x){

            if(any(x[,3] == 0)){
                logP           <- ifelse(x[,3] == 0, yes = 2*pnorm(abs(x[,1])/x[,2],log.p = TRUE, lower.tail = FALSE)/log(10), no = log10(x[,3]))
                Observed       <- -sort(logP)
            }else{
                Observed    <- -log10(sort(x[,3]))
            }

            Expected    <- -log10(1:nrow(x)/nrow(x))
            fdr         <-  sort(p.adjust(x[,3],method = "BH"))
            bonferroni  <-  sort(p.adjust(x[,3],method = "bonferroni"))
            Significant <-  factor(ifelse(bonferroni < 0.05, yes = "Bonferroni", no = ifelse(fdr < 0.05, yes = "FDR", no = "NonSign")))

            out <- data.frame(Observed, Expected, fdr, bonferroni, Significant)
            return(out)
}
