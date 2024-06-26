## Scripts to plot forecast for the multi-stock SAM model (until the package is updated)
## Author: Christoffer Moesgaard Albertsen (cmoe@aqua.dtu.dk)
## Year: 2023
## Copyright: Technical University of Denmark

axisInside <- function(side){
    i1 <- if(side %in% c(2,4)){3:4}else{1:2}
    i2 <- if(side == 2){1}else if(side == 4){2}else if(side == 1){3}else{4}
    i3 <- if(side %in% c(2,4)){1:2}else{3:4}
    ticks <- axisTicks(par("usr")[i1],FALSE)
    sgn <- if(side %in% c(1,2)){ 1 }else{ -1}
    pos <- ((1:4 + 1) %% 4 + 1)[side]
    tA <- par("usr")[i2]+ sgn * 0.01 * diff(par("usr")[i3])
    tB <- ticks
    tx <- if(side %in% c(2,4)){ tA }else{ tB }
    ty <- if(side %in% c(1,3)){ tA }else{ tB }
    text(tx,ty,ticks,pos = pos)
    sx0 <- if(side %in% c(2,4)){ par("usr")[i2] }else{ ticks }
    sx1 <- if(side %in% c(2,4)){ par("usr")[i2]+ sgn * 0.01 * diff(par("usr")[i3]) }else{ ticks }
    sy0 <- if(side %in% c(1,3)){ par("usr")[i2] }else{ ticks }
    sy1 <- if(side %in% c(1,3)){ par("usr")[i2]+ sgn * 0.01 * diff(par("usr")[i3]) }else{ ticks }
    segments(sx0,sy0,sx1,sy1)
}


collapse <- function(x,
                     collap=c(rep(", ",length(x)-2),ifelse(length(x)>2,", and "," and "))){
    if(length(x) == 1)
        return(x)
    if(length(collap) < length(x)-1)
        collap <- rep(collap,len=length(x)-1)
    str <- x[1]
    for(i in 2:length(x))
        str <- paste(str,x[i],sep=collap[i-1])
    str
}
tofrac <- Vectorize(function(x, dollar = TRUE){
    if(dollar){
    paste0("$\\frac{",
           paste0(strsplit(as.character(MASS::as.fractions(x)),"/")[[1]],collapse = "}{"),
           "}$")
    }else{
    paste0("\\frac{",
           paste0(strsplit(as.character(MASS::as.fractions(x)),"/")[[1]],collapse = "}{"),
           "}")
    }
})
displayNum <- Vectorize(function(x,capitalize = FALSE){
    y <- switch(as.character(x),
                "0"="Zero",
                "1"="One",
                "2"="Two",
                "3"="Three",
                "4"="Four",
                "5"="Five",
                "6"="Six",
                "7"="Seven",
                "8"="Eight",
                "9"="Nine",
                formatC(as.numeric(x),big.mark=",",small.mark=".",digits=0,format="f"))
    if(!capitalize & x < 10)
        return(tolower(y))
    return(y)                
},"x")

formatDate <- function(x, format="",locale = Sys.getlocale("LC_TIME"), ...){
    ol <- Sys.getlocale("LC_TIME");
    Sys.setlocale("LC_TIME",locale)
    v <- strftime(x,format=format,...)
    ## invisible(apply(matrix(strsplit(ol,"(;|=)")[[1]],nrow=2),2,function(xx){
    ##         tryCatch({Sys.setlocale(xx[1],xx[2])},error=function(x)"")
    ##         }))
    Sys.setlocale("LC_TIME",ol)
    v    
}

fd <- function(x,digits = 1) formatC(x, digits = digits, format="f")


dtucols <- function(x){
    ## DTU color palette
    names <- c("corporate", "purple", "orange", "blue", "grey", "yellow", "black", "navy", "green", "red", "pink", "brightgreen")
    if(missing(x))
        x <- names
    cols <- c("#990000", # Red
              "#79238E", #Purple
              "#FC7634", #Orange
              "#2F3EEA", #Blue
              "#DADADA", #Grey
              "#F6D04D", #Yellow
              "#000000", #Black
              "#030F4F", #Navy
              "#008835", #Green
              "#E83F48", #Red2
              "#F7BBB1", #Pink
              "#1FD082" #Bright green          
              )
    names(cols) <- names
    if(is.numeric(x)){
        out <- cols[x]
    }else if(is.character(x)){
        out <- cols[match(x, names)]
    }else{
        stop("Wrong input")
    }
    out    
}
palette(dtucols())



makePlot <- function(fit, F, stock = seq_len(length(fit)+1),cols = c(seq_along(fit),"black")[stock], ylim = NULL, xlim = NULL, trans = 0.3, lwd = 5, xlab = "Year", ylab = F, legend = c(names(fit),"Total")[stock], legend.pos = "top", legend.ncol = length(stock), add = FALSE, ci = TRUE, axes=TRUE){
    ## plotfn <- eval(parse(text = paste0(F,"plot")))
    tabfn <- eval(parse(text = paste0(F,"table")))
    tabL <- tabfn(fit, addTotal=TRUE, returnList = TRUE)[stock]
    if(F%in%c("catch","land") && any(stock == length(fit)+1))
        tabL[[which(stock == length(fit)+1)]] <- head(tabL[[which(stock == length(fit)+1)]],-1)
    if(F=="fbar")
      tabL<-lapply(tabL, function(x) x[-nrow(x),])
    if(is.null(ylim))
        ylim <- range(tabL,finite=TRUE,na.rm=TRUE)
    if(is.null(xlim))
        xlim <- range(sapply(tabL,function(x)as.numeric(rownames(x))))
    if(!add){
        plot(0,0,type="n", ylim = ylim, xlim = xlim, xlab=ifelse(axes,xlab,NA),ylab=ifelse(axes,ylab,NA), las = 1,axes=axes)
        grid()
    }
    if(ci){
        for(i in seq_along(stock)){
            xx <- as.numeric(rownames(tabL[[i]]))
            polygon(c(xx,rev(xx)), c(tabL[[i]][,2],rev(tabL[[i]][,3])), col = caMisc::addTrans(cols[i],trans), border = NA)
        }
    }
    for(i in seq_along(stock)){
        xx <- as.numeric(rownames(tabL[[i]]))
        lines(xx,tabL[[i]][,1], col = cols[i], lwd = lwd)
    }
    if(!add)
        legend(legend.pos, legend = legend, col = cols, lwd = lwd, ncol = legend.ncol, bty="n")
    
}

makeRetroPlot <- function(fl, F, stock = seq_len(length(basefit)+1),cols = c(seq_along(fit),"black")[stock], peel_col = NULL,ylim=NULL,xlim=NULL, ...){
    tabfn <- eval(parse(text = paste0(F,"table")))
    tabL0 <- tabfn(attr(fl,"fit"), addTotal=TRUE, returnList = TRUE)[stock]
    tabL <- lapply(lapply(fl, tabfn, addTotal=TRUE, returnList = TRUE),"[",stock)
    if(is.null(ylim))
        ylim <- range(0,c(tabL0,lapply(tabL,function(x)lapply(x,function(y)y[,1]))),finite=TRUE,na.rm=TRUE)
    if(is.null(xlim))
        xlim <- range(lapply(tabL0,function(x)as.numeric(rownames(x))),sapply(tabL,function(x)range(sapply(x,function(y)as.numeric(rownames(y))))))
    makePlot(attr(fl,"fit"),F=F, stock=stock,xlim=xlim,ylim=ylim, cols=cols,...)
    if(is.null(peel_col))
        peel_col <- function(i)caMisc::addTint(cols,seq(0.7,0.2,len=length(fl)+1)[i+1])
    a <- lapply(rev(seq_along(fl)), function(i)makePlot(fl[[i]],F=F, stock=stock,add=TRUE,ci=FALSE, cols = peel_col(i)))  
}

makeRetroPanel <- function(fl,F,ylab=F, addRho = TRUE, peel_col = NULL){
    par(mfrow = c(2,2), mar=c(0.1,0.1,0.1,0.1), oma = c(6,6,6,6))
    bf <- attr(fl,"fit")
    ns <- length(bf)
    if(addRho){
        ## mm10 <- mohn(fl,addTotal=TRUE)
        mm5 <- mohn(tail(fl,5), addTotal=TRUE)
    }
    for(ss in seq_len(ns+1)){
        if(is.null(peel_col)){
            peel_col0 <- function(i)addTint(ifelse(ss>ns,"black",ss),seq(0.5,0,len=length(fl))[i])
        }else{
            peel_col0 <- peel_col
        }

        makeRetroPlot(fl,F,stock=ss,cols="grey",peel_col=peel_col0,legend.pos=NA,axes=FALSE)
        if(is.null(names(fl)))
            names(fl) <- attr(bf,"m_data")$maxYearAll - rev(seq_along(fl))
        legend("topright",legend=c("Reference",rev(names(fl))),col=c("grey",sapply(rev(seq_along(fl)),peel_col0)),lwd=3,bty="n")
        mtext(c(names(bf),"Total")[ss],line=-2)
        if(addRho){
            ## vv10 <- mm10[match(F,c("rec","ssb","fbar","catch")),ss]
            vv5 <- mm5[match(F,c("rec","ssb","fbar","catch")),ss]
            ## mtext(paste("Mohn's",expression(rho),"(10 peels) =",sprintf("%.3f",vv10)),line=-3,cex=0.9)
            mtext(paste("Mohn's",expression(rho),"(5 peels) =",sprintf("%.3f",vv5)),line=-4,cex=0.9)
        }
        axis(ifelse(ss %in% c(1,3),2,4),las=1)
        axis(ifelse(ss %in% c(1,2),3,1),las=1)
        mtext("Year",outer=TRUE,side=c(1,1,3,3),line=3,at=c(0.25,0.75))
        mtext(ylab,outer=TRUE,side=c(2,2,4,4),line=4,at=c(0.25,0.75))
        box()
    }
}


makeComparePanel <- function(fl,F,ylab=F, addRho = TRUE){
    par(mfrow = c(2,2), mar=c(0.1,0.1,0.1,0.1), oma = c(6,6,6,6))
    bf <- attr(fl,"fit")
    ns <- length(bf)
    ## if(addRho)
    ##     mm <- mohn(fl,addTotal=TRUE)
    for(ss in seq_len(ns+1)){
        makeRetroPlot(fl,F,stock=ss,cols="black",peel_col=function(i)addTrans(ifelse(ss>ns,"black",ss),0.1),legend.pos=NA,axes=FALSE)
        #legend("topright",legend=c("Reference",2021:2017),col=c("grey",sapply(rev(seq(0.5,0,-0.1)),function(i)addTint(ifelse(ss>ns,"black",ss),i))),lwd=3,bty="n")
        mtext(c(names(bf),"Total")[ss],line=-2)
        ## if(addRho){
        ##     vv <- mm[match(F,c("rec","ssb","fbar","catch")),ss]
        ##     mtext(paste("Mohn's",expression(rho),"=",sprintf("%.3f",vv)),line=-3,cex=0.9)
        ## }
        axis(ifelse(ss %in% c(1,3),2,4),las=1)
        axis(ifelse(ss %in% c(1,2),3,1),las=1)
        mtext("Year",outer=TRUE,side=c(1,1,3,3),line=3,at=c(0.25,0.75))
        mtext(ylab,outer=TRUE,side=c(2,2,4,4),line=4,at=c(0.25,0.75))
        box()
    }
}

bgtext <- function(x, y, labels, ...,
                   cex = 1, font = NULL,
                   bg = "white", bgex = 1, border = NA){
    w <- strwidth(labels, cex = cex, font = font)
    h <- strheight(labels, cex = cex, font = font)
    rect(x - 0.5 * bgex * w,
         y - 0.5 * bgex * h,
         x + 0.5 * bgex * w,
         y + 0.5 * bgex * h,
         col = bg, border = NA)
    text(x,y,labels, cex=cex, font=font, ...)
}


observedCatch <- function(fit){
    aux <- attr(fit,"m_data")$sharedObs$aux
    logobs <- attr(fit,"m_data")$sharedObs$logobs
    ii <- aux[,"fleet"] %in% which(attr(fit,"m_data")$sharedObs$fleetTypes==0)
    CMW <- Reduce("+",lapply(attr(fit,"m_data")$sam,function(x) x$catchMeanWeight)) / length(fit)
    CAY <- xtabs(exp(logobs[ii]) ~ aux[ii,"year"] + aux[ii,"age"] + aux[ii,"fleet"]) * CMW
    apply(CAY,1,sum)
}

grabFit <- function(f){
    e <- new.env()
    load(f,e)   
    ## attr(e$fit,"withBestConf") <- e$withBestConf
    ## attr(e$fit,"validResid") <- !any(is.nan(e$RESID$residual))
    e$fit
}

library(stockassessment)
library(multiStockassessment)
library(caMisc) # install.packages('caMisc', repos = c('https://calbertsen.r-universe.dev', 'https://cloud.r-project.org'))
suppressPackageStartupMessages(library(sf))
suppressPackageStartupMessages(library(rworldxtra))

grabSR <- function(f){
    e <- new.env()
    load(f,e)
    fmsy <- ifelse(!match("MSY",rownames(e$rpD$Northwest$tables$F)),NA,do.call(sprintf,c(fmt="%.4f (%.4f - %.4f)",as.list(e$rpD$Northwest$tables$F["MSY",]))))
    data.frame(Name = stockassessment:::recruitmentProperties(e$fit[[1]])$name,
               AIC = AIC(e$fit),
               "F~MSY~" = fmsy)
}


`[.msam_retro` <- function(x,i){
    x2 <- unclass(x)
    v <- x2[i]
    attr(v,"fit") <- attr(x,"fit")
    class(v) <- class(x)
    names(v) <- names(x)[i]
    v
}


makeSensitivityPanel <- function(fl, F, ylab) makeRetroPanel(fl,F,ylab=ylab,addRho=FALSE, peel_col = function(i) addTrans(i+3 + (i>=(2)),ifelse(attr(fl[[i]],"m_sdrep")$pdHess && attr(fl[[i]],"m_opt")$convergence==0,1,0.1)))



landtable <- function (fit, obs.show = FALSE, returnList = FALSE, ...) 
{
    CW <- lapply(attr(fit, "m_data")$sam, function(x) x$landMeanWeight)
    xx <- lapply(CW, function(x) as.integer(rownames(x[, , 1])))
    ret <- multiStockassessment:::tableit.msam(fit, x = xx, "logLand", trans = exp, returnList = TRUE, 
        ...)
    if (obs.show) {
        sopStock <- lapply(fit, function(ff) catchtable(ff, obs.show = TRUE)[, 
            "sop.catch"])
        for (i in 1:length(ret)) ret[[i]] <- cbind(ret[[i]], 
            sop.catch = sopStock[[i]])
    }
    if (returnList) 
        return(ret)
    return(cbindYearTables(ret))
}


getFCTab <- function(ff, F, stock){
    if(stock == length(ff)+1){
        ## if(F=="fbar"){
        ##     fulltab <- attr(ff$result[[1]],"tab")
        ##     return(fulltab[,grep(sprintf("^%s:",F),colnames(fulltab))]*NA)
        ## }
        vv <- Reduce("+",lapply(ff$result, function(r) sapply(r,function(x) x[[F]])))
        if(F=="fbar")
            vv <- vv / length(ff$result)
        tab <- t(apply(vv,2, quantile, c(0.5,0.025,0.975)))
        colnames(tab) <- paste0(F,":",c("median","low","high"))
        rownames(tab) <- sapply(ff$result[[1]],function(x) x$year)
    }else{
        fulltab <- attr(ff$result[[stock]],"tab")
        tab <- fulltab[,grep(sprintf("^%s:",F),tolower(colnames(fulltab)))]
    }
    if(F != "ssb")
        tab <- head(tab,-1)
    tab
}


addFC <- function(ii, rL, F, stock){
    ff <- rL[[ii]]
    tab <- getFCTab(ff,F,stock)
    yy <- as.numeric(rownames(tab))
    offset <- seq(-0.2,0.2,len=length(rL))[ii]
    arrows(yy+offset,tab[,2],yy+offset,tab[,3], angle=90,code=3,length=0.1,lwd=3,col=addTrans(ii+3+(ii>=2),0.5))
    points(yy+offset,tab[,1], cex=1.5,pch=19, col = ii+3+(ii>=2))
}


## fitMS is the fit
## rL is a list of forecasts (i.e., FC from the forecast script on TAF)
## F is what should be plotted (e.g., F="ssb" for SSB) Can be anything with a *table function.
makeFCPanel <- function(fitMS, rL,F,ylab=F){
    par(mfrow = c(2,2), mar=c(0.1,0.1,0.1,0.1), oma = c(6,6,6,6))
    ns <- length(fitMS)
    tabfn <- eval(parse(text = paste0(F,"table")))
 
    for(ss in seq_len(ns+1)){
        ylm <- range(0,tabfn(fitMS,returnList=TRUE,addTotal=TRUE)[[ss]],lapply(rL,FUN=getFCTab, F=F,stock=ss), na.rm=TRUE)
        xlm <- range(as.numeric(rownames(tabfn(fitMS,returnList=TRUE,addTotal=TRUE)[[ss]])),
                     lapply(lapply(rL,FUN=getFCTab, F=F,stock=ss),function(x) as.numeric(rownames(x))), na.rm=TRUE)
        makePlot(fitMS,F,ss,xlim=xlm, ylim = ylm,cols="grey", axes=FALSE,legend.pos="bottom")
        if(ss <= length(fitMS)){
            if(F == "fbar"){
                abline(h=Fmsy[ss],lty=4,col="darkgrey",lwd=3)
                bgtext(par("usr")[1]+strwidth(expression(F[MSY]))/2 + 0.01 * diff(par("usr")[1:2]),Fmsy[ss],expression(F[MSY]), bg = rgb(1,1,1,0.3))
            }else if(F == "ssb"){
                abline(h=Blim[ss],lty=4,col="darkgrey",lwd=3)
                bgtext(par("usr")[1]+strwidth(expression(B[lim]))/2 + 0.01 * diff(par("usr")[1:2]),Blim[ss],expression(B[lim]), bg = rgb(1,1,1,0.3))
                abline(h=MSYBtrigger[ss],lty=4,col="darkgrey",lwd=3)
                bgtext(par("usr")[1]+strwidth(expression(MSY~B[trigger]))/2 + 0.01 * diff(par("usr")[1:2]),MSYBtrigger[ss],expression(MSY~B[trigger]), bg = rgb(1,1,1,0.3))
            }
        }
        invisible(lapply(seq_along(rL),FUN=addFC,rL=rL,F=F,stock=ss))
        legend("topleft",sapply(rL,function(x)x$name),col=seq_along(rL)+3+(seq_along(rL)>=2),pch=19,bty="n")
        axis(ifelse(ss %in% c(1,3),2,4),las=1)
        axis(ifelse(ss %in% c(1,2),3,1),las=1)
        mtext("Year",outer=TRUE,side=c(1,1,3,3),line=3,at=c(0.25,0.75))
        mtext(ylab,outer=TRUE,side=c(2,2,4,4),line=4,at=c(0.25,0.75))
        box()
    }
}


## ff should be e.g. FC[[1]] from the forecast script on TAF
## Can be used as: do.call("rbind",lapply(FC, MakeScenarioLine))
## To give all the lines
MakeScenarioLine <- function(ff){
    tabs <- lapply(ff$result,attr,which="tab")
    CC <- lapply(tabs, function(x) sprintf("%d\\newline\\color{gray}(%d - %d)\\color{black}",x["2023","catch:median"],x["2023","catch:low"],x["2023","catch:high"]))
    LL <- lapply(tabs, function(x) sprintf("%d\\newline\\color{gray}(%d - %d)\\color{black}",x["2023","Land:median"],x["2023","Land:low"],x["2023","Land:high"]))
    DD <- lapply(tabs, function(x) sprintf("%d\\newline\\color{gray}(%d - %d)\\color{black}",x["2023","Discard:median"],x["2023","Discard:low"],x["2023","Discard:high"]))
    FF <- lapply(tabs, function(x) sprintf("%.3f\\newline\\color{gray}(%.3f - %.3f)\\color{black}",x["2023","fbar:median"],x["2023","fbar:low"],x["2023","fbar:high"]))
    FFL <- lapply(tabs, function(x) sprintf("%.3f\\newline\\color{gray}(%.3f - %.3f)\\color{black}",x["2023","fbarL:median"],x["2023","fbarL:low"],x["2023","fbarL:high"]))
    FFD <- lapply(tabs, function(x) sprintf("%.3f\\newline\\color{gray}(%.3f - %.3f)\\color{black}",x["2023","fbarD:median"],x["2023","fbarD:low"],x["2023","fbarD:high"]))
    SS <- lapply(tabs, function(x) sprintf("%d\\newline\\color{gray}(%d - %d)\\color{black}",x["2024","ssb:median"],x["2024","ssb:low"],x["2024","ssb:high"]))
    pctSSB <- lapply(ff$result, function(x){
        vv <- (x[[4]]$ssb / x[[3]]$ssb - 1)*100
        sprintf("%.1f %%\\newline\\color{gray}(%.1f %% - %.1f %%)\\color{black}",median(vv),quantile(vv,0.025),quantile(vv,0.975))
    })
    pctTAC <- lapply(seq_along(tabs), function(ii) sprintf("%.1f %%\\newline\\color{gray}(%.1f %% - %.1f %%)\\color{black}",(tabs[[ii]]["2023","catch:median"] / IntermediateYearCatch[[ii]] - 1)*100,(tabs[[ii]]["2023","catch:low"]/IntermediateYearCatch[[ii]]-1)*100,(tabs[[ii]]["2023","catch:high"] - IntermediateYearCatch[[ii]] - 1)*100))
    pctBelowBlim <- lapply(seq_along(ff$result), function(i){
        vv <- (ff$result[[i]][[4]]$ssb < Blim[i])
        sprintf("%.1f %%",mean(vv)*100)
    })
    res <- lapply(seq_along(tabs), function(ii){
        yy <- rownames(tabs[[ii]])
        ll <- list(gsub("TAC","Catch",ff$name),
                   unname(CC[[ii]]),
                   unname(LL[[ii]]),
                   unname(DD[[ii]]),
                   unname(FF[[ii]]),
                   unname(FFL[[ii]]),
                   unname(FFD[[ii]]),
                   unname(SS[[ii]]),
                   unname(pctSSB[[ii]]),
                   unname(pctTAC[[ii]]),
                   unname(pctBelowBlim[[ii]]))
        nn <- c("Basis",
                       sprintf("Total\ catch\ (%s)",head(tail(yy,2),1)),
                       sprintf("Projected\ landings\ (%s)",head(tail(yy,2),1)),
                       sprintf("Projected\ discards\ (%s)",head(tail(yy,2),1)),
                       sprintf("F\\textsubscript{total}\ (ages\ 2\ -\ 4)\ (%s)",head(tail(yy,2),1)),
                       sprintf("F\\textsubscript{landings}\ (ages\ 2\ -\ 4)\ (%s)",head(tail(yy,2),1)),
                       sprintf("F\\textsubscript{discard}\ (ages\ 2\ -\ 4)\ (%s)",head(tail(yy,2),1)),
                       sprintf("SSB\ (%s)",tail(yy,1)),
                       "%\ SSB\ change",
                       "%\ Catch\ change",
                       sprintf("%%\ probability\ of\ falling\ below\ Blim\ in\ %s",tail(yy,1)))
        as.data.frame(ll,check.names=FALSE,cut.names=FALSE,col.names=nn)
    })
    names(res) <- names(tabs)
    res
}

