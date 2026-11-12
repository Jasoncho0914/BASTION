# Plot against numeric coordinates internally; preserve calendar axis labels.
.bastion_plot_panel <- function(time,mean,lower=NULL,upper=NULL,main="",
    xlab="Time",ylab="Value",show_x=TRUE,col="black",band_col=grDevices::adjustcolor("#0072B2",.50),...) {
  xx<-as.numeric(time)
  order<-order(xx)
  xx<-xx[order];mean<-mean[order]
  limits<-range(c(mean,lower,upper),finite=TRUE)
  if(!length(xx) || !all(is.finite(limits)))stop("No finite values to plot")
  date_axis<-inherits(time,c("POSIXct","Date"))
  graphics::plot(xx,mean,type="n",ylim=limits,main=main,xlab=xlab,ylab=ylab,
    xaxt=if(date_axis || !show_x)"n" else "s",...)
  if(date_axis && show_x) {
    if(inherits(time,"Date"))graphics::axis.Date(1,at=pretty(time)) else
      graphics::axis.POSIXct(1,x=time,at=pretty(time))
  }
  if(!is.null(lower) && !is.null(upper)) {
    lower<-lower[order];upper<-upper[order]
    good<-is.finite(xx)&is.finite(lower)&is.finite(upper)
    indices<-which(good)
    groups<-if(length(indices))split(indices,cumsum(c(TRUE,diff(indices)!=1))) else list()
    for(ii in groups)if(length(ii))graphics::polygon(c(xx[ii],rev(xx[ii])),
      c(lower[ii],rev(upper[ii])),col=band_col,border=NA)
  }
  graphics::lines(xx,mean,col=col)
}

#' Plot a BASTION decomposition
#' @param x A `bastion_fit` object.
#' @param ... Additional arguments passed to each base graphics plot.
#' @return The fitted object, invisibly.
#' @export
plot.bastion_fit <- function(x,...) {
  old<-graphics::par(no.readonly=TRUE);on.exit(graphics::par(old))
  panels<-ncol(x$components)+1L
  graphics::par(mfrow=c(panels,1L),mar=c(2.5,5.5,1.8,1),oma=c(3,0,2,0),
    mgp=c(3.5,.8,0),cex.axis=.8,cex.lab=.85,cex.main=.9,las=1)
  .bastion_plot_panel(x$time,x$observed,main="Observed",xlab="",show_x=FALSE,...)
  for(name in names(x$components))
    .bastion_plot_panel(x$time,x$components[[name]],
      lower=x$intervals$lower[[name]],upper=x$intervals$upper[[name]],
      main=name,xlab="",show_x=identical(name,tail(names(x$components),1L)),...)
  cl<-x$metadata$fitting$cl
  label<-if(is.null(cl))"Shaded bands: credible intervals" else
    paste0("Shaded bands: ",format(100*cl,trim=TRUE),"% credible intervals")
  graphics::mtext("Time",side=1,outer=TRUE,line=1,cex=.8)
  graphics::mtext(label,side=3,outer=TRUE,line=.4,cex=.8)
  invisible(x)
}

#' Plot a BASTION component with credible intervals
#' @param fit A `bastion_fit` object.
#' @param component Exact column name in `fit$components`.
#' @param part Seasonal view: `all` (default), `total`, `recurring`, or `deviation`.
#' @param col Posterior mean line color.
#' @param band_col Credible band fill color.
#' @param ... Additional arguments passed to each base graphics plot.
#' @return The fitted object, invisibly.
#' @export
plot_component <- function(fit,component,part=c("all","total","recurring","deviation"),
    col="black",band_col=grDevices::adjustcolor("#0072B2",.50),...) {
  if(!inherits(fit,"bastion_fit"))stop("fit must be a bastion_fit object")
  if(length(component)!=1L || !component %in% names(fit$components))
    stop("component must be one of: ",paste(names(fit$components),collapse=", "))
  part<-match.arg(part)
  seasonal<-component %in% names(fit$seasonality)
  if(!seasonal && !part %in% c("all","total"))
    stop("recurring and deviation apply only to seasonal components")
  views<-if(seasonal && part=="all")c("total","recurring","deviation") else
    if(part=="all")"total" else part
  old<-graphics::par(no.readonly=TRUE);on.exit(graphics::par(old))
  graphics::par(mfrow=c(length(views),1L),mar=c(4.5,5.5,2.5,1),mgp=c(3,.8,0),
    cex.axis=.8,cex.lab=.85,cex.main=.9,las=1)
  cl<-fit$metadata$fitting$cl
  level<-if(is.null(cl))"credible interval" else paste0(format(100*cl,trim=TRUE),"% credible interval")
  for(view in views) {
    if(view=="total") {
      tt<-fit$time;mu<-fit$components[[component]]
      lo<-fit$intervals$lower[[component]];hi<-fit$intervals$upper[[component]]
      axis_label<-"Time"
    } else {
      block<-fit$seasonality[[component]][[view]]
      mu<-block$mean;lo<-block$lower;hi<-block$upper
      tt<-if(view=="recurring")fit$metadata$mapping$recurring_positions[[component]] else fit$time
      axis_label<-if(view=="recurring") {
        if(is.null(fit$metadata$mapping$seasonality))"Position within cycle" else "Position within cycle [0, 1)"
      } else "Time"
    }
    .bastion_plot_panel(tt,mu,lo,hi,
      main=paste0(component,if(seasonal)paste0(" - ",view) else ""," (",level,")"),
      xlab=axis_label,col=col,band_col=band_col,...)
  }
  invisible(fit)
}
