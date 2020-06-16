plot.quantiles <- function(fit,z,...){
  require("ggplot2")
  plot.options <- list(...)
  if(is.null(plot.options$title)){plot.options$title <- "Quantiles plot"}
  if(is.null(plot.options$xlabelaxis)){plot.options$xlabelaxis <- "Years"}
  if(is.null(plot.options$ylabelaxis)){plot.options$ylabelaxis <- "Returns"}
  if(is.null(plot.options$pdftitle)){plot.options$pdftitle <- "QuantilesPlot"}
  if(is.null(plot.options$wd)){plot.options$wd <- getwd()}
  color <- NULL
  color[1] <- "blue"
  color[2] <- "red"
  color[3] <- "green"
  num.quant <- dim(fit$quant.mat)[2]
  plotdata <- data.frame(y=z$y,date=z$date)
  quant.mat <- fit$quant.mat
  if(num.quant>7){stop('too many quantiles to plot (functions handles 7 thetas only). lower the number of thetas')}
  if(num.quant%%2==0){warning('number of quantiles is even. check and set to odd (including median, theta = 0.5)')}
  p <- ggplot(plotdata, aes(x=as.Date(plotdata$date), y=plotdata$y) ) + #x must have class "date"
    geom_line(color='#c0392b') + xlab(plot.options$xlabelaxis) + ylab(plot.options$ylabelaxis) + ggtitle(plot.options$title)
  p + theme(axis.text.x=element_text(angle=60, hjust=1)) + 
    scale_x_date(date_labels = "%Y %m %d", breaks = function(x) seq.Date(from = min(x), to = max(x),by = "1 years" ),
                 minor_breaks = function(x) seq.Date(from = min(x), to = max(x), by = "2 years")) + 
    for(i in 1:floor(num.quant/2)){
      p <- p + geom_line(aes_string(y = quant.mat[,i]))
      p <- p + geom_line(aes_string(y = quant.mat[,dim(quant.mat)[2]+1-i])) 
      p <- p + geom_ribbon(aes_string(ymin=quant.mat[,i],ymax=quant.mat[,dim(quant.mat)[2]+1-i]), fill=color[i], alpha="0.5", show.legend = TRUE) 
    }
  p <- p + geom_line(aes_string(y = quant.mat[,ceiling(num.quant/2)]), color='grey', alpha="0.5", show.legend = TRUE)
  setwd(plot.options$wd)
  pdf(file=plot.options$pdftitle)
  print(p)
  dev.off()
  return(p)
}