detr <- function(x,smoothing) {

spline.p.for.smooth.Pspline<-function(year)  {
p <- .5/((6*(cos(2*pi/year)-1)^2/(cos(2*pi/year)+2)))
return(p)
}


ifelse ((is.null(ncol(x))),x2 <- ts.union(x,x), x2 <- x)

fyarray <- start(x2)[1]
lyarray <- end(x2)[1]
begin <- fy(x2)
end <- ly(x2)

smoothedarray <- x2

spline.p <- spline.p.for.smooth.Pspline(smoothing)


for (i in 1:ncol(x2)) {
temp <- smooth.Pspline((begin[i]:end[i]),x2[(begin[i]-fyarray+1):(end[i]-fyarray+1),i],spar=spline.p,method=1)
smoothedarray[(begin[i]-fyarray+1):(end[i]-fyarray+1),i]  <- temp$ysmth
}

if (is.null(ncol(x))) {smoothedarray <- smoothedarray[,-1]}

residualarray <- x-smoothedarray
ratioarray <- x/smoothedarray

return(ratioarray)
}

