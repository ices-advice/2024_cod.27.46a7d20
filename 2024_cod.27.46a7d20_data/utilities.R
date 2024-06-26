boot.weight.at.age <-function(d,minAge,maxAge,WALmods,Nboot=5){
  
  dkeep <-d
  
  area <-unlist(strsplit(as.character(dkeep$YearsubArea[1]),":"))[2]
  print(as.character(dkeep$YearsubArea[1]))
  WALmodel <-WALmods[[area]]
  
  
  hauls <-dkeep[[2]]$haul.id
  bwts <-weight.at.age.fixed.lw(dkeep,minAge=minAge,maxAge=maxAge,WALmodel=WALmodel)
  # rownames(bwts) <-"Value"
  
  set.seed(12345)
  
  one.boot <-function(i,minAge,maxAge,WALmodel){
    #    boothauls <-sample(hauls,replace=TRUE)
    #    d[[2]] <-do.call("rbind",lapply(boothauls,function(x){subset(dkeep[[2]],haul.id==x)}))
    #    d[[1]] <-do.call("rbind",lapply(boothauls,function(x){subset(dkeep[[1]],haul.id==x)}))
    #    tmp <-weight.at.age.fixed.lw(d,minAge=minAge,maxAge=maxAge,WALmodel=WALmodel)
    
    b.idx <-sample(1:length(hauls),replace=TRUE)
    boothauls <-hauls[b.idx]
    bhcount <-table(boothauls)
    
    d[[2]] <-dkeep[[2]][b.idx,]
    
    d[[1]] <-dkeep[[1]]
    d[[1]]$nr <-0
    d[[1]]$nr <-bhcount[match(d[[1]]$haul.id,names(bhcount))]
    row.idx <-which(d[[1]]$nr!=0)
    ca.rows <-rep(row.idx,times=d[[1]]$nr[row.idx])
    d[[1]] <-d[[1]][ca.rows,]
    
    tmp <-weight.at.age.fixed.lw(d,minAge=minAge,maxAge=maxAge,WALmodel=WALmodel)
    
    
    return(tmp)
  }
  
  out <-lapply(1:Nboot,one.boot,minAge,maxAge,WALmodel)
  
  bwts <-bind_rows(bwts,out)
  
  bwts <-as.matrix(bwts)
  #  dimnames(bwts)[[2]] <-paste0("Age",seq(minAge,maxAge,1))
  dimnames(bwts)[[2]] <-paste0("Age",colnames(bwts))
  dimnames(bwts)[[1]] <-c("Value",paste0("b",1:Nboot))
  
  return(bwts)
}


weight.at.age.fixed.lw<-function(d,minAge,maxAge,WALmodel){
  # Uses externally fitted Wt-Len model
  checkSpectrum(d);
  
  d[[1]] = subset(d[[1]],Age>=minAge);
  d[[1]]$Age[ d[[1]]$Age > maxAge ] = maxAge;
  cm.b=attr(d,"cm.breaks")
  
  WL = exp(predict(WALmodel,newdata=data.frame(LngtCm=cm.b[-length(cm.b)],IndWgt=NA)))
  
  pla = length.given.age(d,minAge,maxAge);
  if(!is.null(pla)) { 
    ## find E(w(l) | a ) = sum( w(l)*p(l|a) )
    Ewa = numeric(ncol(pla))
    for(a in 1:ncol(pla)){
      Ewa[a] = sum( WL*pla[,a],na.rm=TRUE)
    }
    #  
    Ewa[Ewa==0] <-NA
    names(Ewa) <-colnames(pla)
  }else{
    Ewa <-rep(NA,maxAge-minAge+1)
    names(Ewa) <-as.character(minAge:maxAge)
  }
  Ewa  
}


length.given.age<-function(d,minAge,maxAge,ALK=NULL){
  d[[1]] = subset(d[[1]],Age>=minAge);
  d[[1]]$Age[ d[[1]]$Age > maxAge ] = maxAge;
  cm.b=attr(d,"cm.breaks") 
  ## find length distr. given age p( l | a)
  ## p( l_i | a ) = p( a | l_i ) * p(l_i) / sum( p( a | l) p(l) )
  pl = colSums(d$N)/sum(d$N)
  d[[1]]$sizeGroup <- cut(d[[1]]$LngtCm, breaks = cm.b, 
                          right = FALSE)
  if(is.null(ALK)){
    if(any( minAge:maxAge %in% d[[1]]$Age) ){
      tab = xtabs(NoAtALK ~ sizeGroup + Age, data = d[[1]])
      pal = tab/rowSums(tab)
    }  else return(NULL);
  } else pal=ALK;
  
  pla = pal;
  for(i in 1:ncol(pla)){
    pa = sum(pal[,i]*pl,na.rm=TRUE)
    pla[,i] = pal[,i]*pl / pa; 
  }
  
  pla
}


weight.len.mod <-function(indat){
  dat <-indat[[1]]
  lwmod <-lm(log(IndWgt)~I(log(LngtCm)),data=subset(dat,!is.na(IndWgt) & IndWgt>0))
  return(lwmod)
}


weight.len <-function(indat){
  dat <-indat[[1]]
  lwmod <-lm(log(IndWgt)~I(log(LngtCm)),data=subset(dat,!is.na(IndWgt) & IndWgt>0))
  coef <-list(loga=coef(lwmod)[[1]],b=coef(lwmod)[[2]],loga.sd=summary(lwmod)$coefficients[1,2],b.sd=summary(lwmod)$coefficients[2,2])
  return(coef)
}




na2zero <- function(x){
  x[is.na(x)] <- 0
  x
}


####' Function to cut years and ages from data by summing
####'
####' Years outside the range (minYear,maxYear) are removed.
####' Ages are summed to form a new plus group at age maxAge
####' @param x matrix with data
####' @return matrix with new data
####' @author Christoffer Moesgaard Albertsen
cutSum <- function(x){
  a <- as.numeric(colnames(x))
  y <- as.numeric(rownames(x))
  if(all(a <= maxAge) && all(y >= minYear) && all(y <= maxYear)) return(x)
  at <- attributes(x)
  v <- cbind(x[y>=minYear & y<=maxYear,a<maxAge,drop=FALSE],rowSums(x[y>=minYear & y<=maxYear,a>=maxAge,drop=FALSE]))
  colnames(v) <- a[a<=maxAge]
  nms <- setdiff(names(at),names(attributes(v)))
  attributes(v)[nms] <- at[nms]
  v
}

####' Function to cut years and ages from data by dropping
####'
####' Years outside the range (minYear,maxYear) are removed.
####' Ages are removed to form a new plus group at age maxAge
####' @param x matrix with data
####' @return matrix with new data
####' @author Christoffer Moesgaard Albertsen
cutDrop <- function(x, maxA = maxAge){
  a <- as.numeric(colnames(x))
  y <- as.numeric(rownames(x))
  if(all(a <= maxA) && all(y >= minYear) && all(y <= maxYear)) return(x)
  at <- attributes(x)
  v <- x[y>=minYear & y<=maxYear,a<=maxA,drop=FALSE]
  nms <- setdiff(names(at),names(attributes(v)))
  attributes(v)[nms] <- at[nms]
  v
}

####' Function to cut years and ages from data by averaging
####'
####' Years outside the range (minYear,maxYear) are removed.
####' Ages are averages to form a new plus group at age maxAge
####' @param x matrix with data
####' @return matrix with new data
####' @author Christoffer Moesgaard Albertsen
cutMean <- function(x, w = matrix(1,nrow(x),ncol(x))){
  a <- as.numeric(colnames(x))
  y <- as.numeric(rownames(x))
  if(all(a <= maxAge) && all(y >= minYear) && all(y <= maxYear)) return(x)
  at <- attributes(x)
  v1 <- x[y>=minYear & y<=maxYear,a<maxAge,drop=FALSE]
  x2 <- x[y>=minYear & y<=maxYear,a>=maxAge,drop=FALSE]
  w2 <- w[y>=minYear & y<=maxYear,a>=maxAge,drop=FALSE]
  rsw2 <- rowSums(w2)
  v2 <- rowSums(x2*w2) / rsw2
  v2[is.na(rsw2) | rsw2 == 0] <- NA
  v <- cbind(v1,v2)
  colnames(v) <- a[a<=maxAge]
  rownames(v) <- y[y>=minYear & y<=maxYear]
  nms <- setdiff(names(at),names(attributes(v)))
  attributes(v)[nms] <- at[nms]
  v
}


##------------catch correlations--------------------------
empty <-
  function(x,y,groups,subscripts, panel.number, packet.number) {
    # do nufink
  }

diag.panel.cm <-
  function(draw,axis.line.col,...)
  {
    diag.panel.splom(draw=F, axis.line.col=axis.line.col,...)
  }

main.panel.cm <-
  function(x,y,groups,subscripts, panel.number, packet.number, .tol=0.05) {
    # do sumfink
    panel.xyplot(x,y,pch=19, col=grey(0.35), cex=0.2)
    
    # fit a linear model to the data
    lm1 <- lm(y ~ x)
    x1<-0:20/20
    fit <- suppressWarnings(predict.lm(lm1, newdata=data.frame(x=x1), se.fit=T))
    y1 <- fit$fit
    yu <- y1 + 2*fit$se
    yl <- y1 - 2*fit$se
    
    sig <- identical(anova(lm1)$"Pr(>F)"[1]<.tol,T)
    
    if (sig) {
      line.f <- list(lwd=3, lty=1, col="#000000")
      line.ci <- list(lwd=2, lty=1, col="red4")
    } else {
      line.f <- list(lwd=1, lty=1, col="#0000FF")
      line.ci <- list(lwd=1, lty=2, col="#0000FF")
    }
    
    panel.lines(x1,y1, lwd=line.f$lwd, lty=line.f$lty, col=line.f$col)
    panel.lines(x1,yu, lwd=line.ci$lwd, lty=line.ci$lty, col=line.ci$col)
    panel.lines(x1,yl, lwd=line.ci$lwd, lty=line.ci$lty, col=line.ci$col)
    
    # draw in axes
    grid.draw(linesGrob(x=c(-0.1, 1, 1),
                        y=c( 0  , 0, 1.1),
                        gp = gpar(lwd = 2, col=grey(0.5)),
                        default.units = "npc"))
  }

panel.pairs.cm <-
  function(z, subscripts, panel.subscripts) {
    panel.pairs(z,
                panel.subscripts = panel.subscripts,
                subscripts       = subscripts,
                lower.panel      = empty,
                diag.panel       = diag.panel.cm,
                upper.panel      = main.panel.cm,
                axis.line.col    = "#FFFFFF",
    )
  }

centre.log <-
  function(mat) {
    mat[mat<=0] <- NA
    mat <- log(mat)
    apply(mat,2,function(x) (x-min(x, na.rm=T))/diff(range(x,na.rm=T)))
  }

plot.index.corr <-
  function(object, wndows=T) {
    par(bty="n")
    trellis.par.set(box.rectangle=list(col="white"))
    for (i in seq(length(object))) {
      #select one tuning fleet
      tune.mat <- t(object[[i]]@catch.n@.Data[,,1,1,1,1]) #NW/GL edit: add sixth dimension
      # make cohort matrix
      n <- dim(tune.mat)[2]
      cohort.mat <- matrix(NA, ncol=n, nrow=dim(tune.mat)[1]+n-1)
      colnames(cohort.mat) <- colnames(tune.mat)
      for (j in 1:n) {
        cohort.mat[,j] <- c(rep(NA,n-j),tune.mat[,j],rep(NA,j-1))
      }
      main <- object[[i]]@name
      if (wndows) windows()
      print(splom(~centre.log(cohort.mat), superpanel=panel.pairs.cm, xlab=main, col="white"))
    }
  }


##------------catchcurve function--------------------------
catchcurve <- function(stock, agerange=c(NULL, NULL)) {
  #Taken without change from rds
  #def.par <- par(no.readonly = TRUE)# save default, for resetting...
  if(!inherits(stock, "FLStock"))
    stop("First input must be an 'FLStock' object!")
  #windows(width=7, height=6)
  #par(mfrow=c(1,1))
  dims   <- dimnames(stock@catch.n)
  arange <- stock@range["min"] :stock@range["max"]
  nages  <- length(as.vector(dims$age))
  
  if(nages > 2) {
    yrange <- stock@range["minyear"]:stock@range["maxyear"]
    catch  <- as.vector(stock@catch.n)
    
    t.     <- cbind(rep(yrange,each=max(arange)-min(arange)+1), arange, catch)
    tdf    <- as.data.frame(cbind(t., t.[,1]-t.[,2]))
    
    suppressWarnings(plot(tdf[,1],log(tdf[,3]),type="null",xlab="year of catch",ylab="log-catch"))
    uc     <- sort(unique(tdf[tdf[,2]==arange,4]))
    for(j in 1:length(uc)) {
      if(!is.null(agerange))
        arange <- agerange[1]:agerange[2]
      v. <- tdf[tdf[,4]==uc[j] ,]
      v. <- v.[is.element(v.[,2], arange),]
      if(length(v.[,1])>1) {
        
        lines(v.[,1],log(v.[,3]),col="red")
        text(v.[,1],log(v.[,3]),as.character(v.[,2]),col="black",cex=0.6)
        
      }
    }
  }
  #par(def.par)
}

#Cod34a    <- no.discards(read.FLStock(paste(path, "Cod347.idx", sep="")))
#catchcurve(Cod34a,c(0,7))


##------------catchcurvegrad function--------------------------
catchcurvegrad <- function(stock, agerange=c(NULL, NULL)) {
  #Adapted from rds by jdo [23/08/05] with slight modification (age range, lm and plotting)
  #def.par <- par(no.readonly = TRUE)# save default, for resetting...
  if(!inherits(stock, "FLStock"))
    stop("First input must be an 'FLStock' object!")
  
  #par(mfrow=c(1,1))
  dims   <- dimnames(stock@catch.n)
  arange <- stock@range["min"] :stock@range["max"]
  nages  <- length(as.vector(dims$age))
  
  if(nages > 2) {
    yrange <- stock@range["minyear"]:stock@range["maxyear"]
    catch  <- as.vector(stock@catch.n)
    
    t.     <- cbind(rep(yrange,each=max(arange)-min(arange)+1), arange, catch)
    tdf    <- as.data.frame(cbind(t., t.[,1]-t.[,2]))
    tdg    <- tdf[tdf$arange>=agerange[1] & tdf$arange<=agerange[2],]
    colnames(tdg) <- list("year", "age", "catch", "cohort")
    tdg <- na.omit(tdg) #NW/GL edit: remove intermediate year
    
    res <- c(0,0)
    for(cht in sort(unique(tdg$cohort))) {
      tdc <- tdg[tdg$cohort==cht,]
      colnames(tdc) <- list("year", "age", "catch", "cohort")
      if(length(tdc$year)==agerange[2]-agerange[1]+1) {
        grad <- lm(log(catch)~age, data  = tdc)
        rss  <- c(cht,-as.numeric(grad$coefficients[2]))
        res  <- rbind(res, rss)
      }
    }
    res<-res[-1,]
    colnames(res) <- c("cohort", "negative gradient")
    plot(res)
    lines(res)
    title(paste("Ages",agerange[1],"to",agerange[2]))
  }
  #par(def.par)
  return(res)
}


##------------cpue_plot function--------------------------------
single_cpue_plot_cohort <- function(cpue) {
  ##	program written by rds, but slightly adapted by jdo
  #def.par <- par(no.readonly = TRUE)# save default, for resetting...
  if(!inherits(cpue, "FLIndices"))
    stop("FLIndex must be an 'FLIndex' object!")
  
  #par(mfrow=c(2,2))
  counter <- 0
  
  for(i in 1:length(cpue)) {
    
    counter <- counter + 1
    if(counter == 5 || counter == 9) {
      ##windows(width=7, height=6)
      ##par(mfrow=c(2,2))
    }
    
    yrange  <- cpue[[i]]@range[4]:cpue[[i]]@range[5]
    arange  <- cpue[[i]]@range[1] :cpue[[i]]@range[2]
    index   <- as.vector(cpue[[i]]@index/rowMeans(cpue[[i]]@index))
    index <- log(index)
    
    t.      <- cbind(rep(yrange,each=length(arange)),arange, index)
    t.      <- cbind(t., t.[,1]-t.[,2])
    
    suppressWarnings(plot(t.[,4],t.[,3],type="null",xlab="Cohort",ylab="Log-mean-standardised index"))
    for(j in sort(unique(t.[,2]))) {
      lines(t.[t.[,2]==j,4],t.[t.[,2]==j,3],col="red")
      text (t.[t.[,2]==j,4],t.[t.[,2]==j,3],as.character(t.[t.[,2]==j,2]),col="black", cex=0.6)
    }
    title(main=cpue[[i]]@name)
  }
  #par(def.par)
}
#################################################################################################

##------------cpue_plot function--------------------------------
single_cpue_plot_year <- function(cpue) {
  ##	program written by rds, but slightly adapted by jdo
  #def.par <- par(no.readonly = TRUE)# save default, for resetting...
  if(!inherits(cpue, "FLIndices"))
    stop("FLIndex must be an 'FLIndex' object!")
  
  ##par(mfrow=c(2,2))
  counter <- 0
  
  for(i in 1:length(cpue)) {
    
    counter <- counter + 1
    if(counter == 5 || counter == 9) {
      ##windows(width=7, height=6)
      ##par(mfrow=c(2,2))
    }
    
    yrange  <- cpue[[i]]@range[4]:cpue[[i]]@range[5]
    arange  <- cpue[[i]]@range[1] :cpue[[i]]@range[2]
    index   <- as.vector(cpue[[i]]@index/rowMeans(cpue[[i]]@index))
    index <- log(index)
    
    t.      <- cbind(rep(yrange,each=length(arange)),arange, index)
    t.      <- cbind(t., t.[,1])
    
    suppressWarnings(plot(t.[,4],t.[,3],type="null",xlab="Year",ylab="Log-mean-standardised index"))
    for(j in sort(unique(t.[,2]))) {
      lines(t.[t.[,2]==j,4],t.[t.[,2]==j,3],col="red")
      text (t.[t.[,2]==j,4],t.[t.[,2]==j,3],as.character(t.[t.[,2]==j,2]),col="black", cex=0.6)
    }
    title(main=cpue[[i]]@name)
  }
  #par(def.par)
}

##------------catchcurve_index function--------------------------
single_catchcurve_index <- function(cpue) {
  ##	plots catch curves for survey indices
  ##	rds 27/01/04
  ##	modified jdo to work with FLCore 16/08/05
  #def.par <- par(no.readonly = TRUE)# save default, for resetting...
  if(!inherits(cpue, "FLIndices"))
    stop("Input must be an 'FLIndices' object!")
  
  ##par(mfrow=c(2,2))
  counter <- 0
  
  for(i in 1:length(cpue)) {
    
    dims  <- dimnames(cpue[[i]]@index)
    nages <- length(as.vector(dims$age))
    if(nages > 2) {
      counter <- counter + 1
      if(counter == 5 || counter == 9) {
        ##windows(width=7, height=6)
        ##par(mfrow=c(2,2))
      }
      yrange <- cpue[[i]]@range["minyear"]:cpue[[i]]@range["maxyear"]
      arange <- cpue[[i]]@range["min"] :cpue[[i]]@range["max"]
      catch  <- as.vector(cpue[[i]]@catch.n)
      effort <- as.vector(cpue[[i]]@effort)
      
      t.     <- cbind(rep(yrange,each=max(arange)-min(arange)+1), arange, catch, rep(effort,each=max(arange)-min(arange)+1))
      cohort <- t.[,1]-t.[,2]
      t.     <- cbind(t.,cohort)
      tdf    <- as.data.frame(t.)
      
      #suppressWarnings(plot(tdf[,1],log(tdf[,3]/tdf[,4]),type="n",xlab="year of survey",ylab="log-abundance index"))
      plot(tdf[,1],log(tdf[,3]/tdf[,4]),type="n",xlab="Year",ylab="log-abundance index")
      uc     <- sort(unique(tdf[,5]))
      for(j in 1:length(uc)) {
        v. <- tdf[tdf[,5]==uc[j],]
        if(length(v.)>1) {
          lines(v.[,1],log(v.[,3]/v.[,4]),col="red")
          text(v.[,1],log(v.[,3]/v.[,4]),as.character(v.[,2]),col="black",cex=0.6)
        }
      }
      title(main=cpue[[i]]@name)
    }
  }
  #par(def.par)
}

##------------catchcurvegrad_index function--------------------------
single_catchcurvegrad_index <- function(cpue,agerange) {
  ##	returns a vector of catch curve gradients for each cohort
  ##	rds 26/07/04
  ##	modified jdo to work with FLCore & FLIndex (instead of FLStock) 16/08/05
  #def.par <- par(no.readonly = TRUE)# save default, for resetting...
  if(!inherits(cpue, "FLIndices"))
    stop("First input must be an 'FLIndices' object!")
  if(length(cpue) != length(agerange) || !is.list(agerange))
    stop("Second input must be list of same length as first input")
  
  #par(mfrow=c(2,2))
  counter <- 0
  reslist <- list(1)
  if(length(cpue)>1) reslist <- c(reslist,c(1:(length(cpue)-1)))
  mina<-agerange
  
  for(i in 1:length(cpue)) {
    
    dims  <- dimnames(cpue[[i]]@index)
    nages <- length(as.vector(dims$age))
    mina[[i]] <- agerange[[i]][2]-agerange[[i]][1]+1
    if(nages > 2) {
      counter <- counter + 1
      if(counter == 5 || counter == 9) {
        #windows(width=7, height=6)
        #par(mfrow=c(2,2))
      }
      yrange    <- cpue[[i]]@range["minyear"]:cpue[[i]]@range["maxyear"]
      arange    <- cpue[[i]]@range["min"]:cpue[[i]]@range["max"]
      cpu       <- cpue[[i]]@index[as.character(arange),]
      logindrat <- cpu
      logindrat[1:nrow(cpu)-1,1:ncol(cpu)-1] <- -log(cpu[1:nrow(cpu)-1,1:ncol(cpu)-1]/cpu[-1,-1])
      logindrat[nrow(cpu),] <- logindrat[nrow(cpu)-1,]
      logindrat[,ncol(cpu)] <- logindrat[,ncol(cpu)-1]
      cpu       <- as.vector(cpu)
      cpu       <- ifelse(cpu==0,NA,cpu)
      logindrat <- as.vector(logindrat)
      
      t.     <- cbind(rep(yrange,each=max(arange)-min(arange)+1), arange, cpu)
      t..    <- cbind(rep(yrange,each=max(arange)-min(arange)+1), arange, logindrat)
      tdf    <- as.data.frame(cbind(t., t.[,1]-t.[,2]))
      tdf.   <- as.data.frame(cbind(t.., t..[,1]-t..[,2]))
      colnames(tdf) <- list("year", "age", "cpue", "cohort")
      colnames(tdf.) <- list("year", "age", "cpue", "cohort")
      
      tdg    <- tdf[tdf$age>=agerange[[i]][1] & tdf$age<=agerange[[i]][2],]
      tdg.    <- tdf.[tdf.$age>=agerange[[i]][1] & tdf.$age<=agerange[[i]][2],]
      
      uc     <- sort(unique(tdg[,4]))
      res    <- c(0,0)
      suppressWarnings(plot(tdg.[,4],-tdg.[,3],type="null",xlab="cohort",ylab="negative gradient"))
      for(j in 1:length(uc)) {
        dat  <- tdg[tdg$cohort==uc[j],]
        if(length(dat$year)==mina[[i]]) {
          grad <- lm(log(cpue)~year, data  = dat)
          rss  <- c(uc[j],-as.numeric(grad$coefficients[2]))
          res  <- rbind(res, rss)
        }
      }
      res <- res[-1,]
      colnames(res) <- c("cohort", "negative gradient")
      lines(res)
      points(res)
      title(main=paste(cpue[[i]]@name,"- ages",agerange[[i]][1],"to",agerange[[i]][2],sep=" "))
      reslist[[i]]<-res
    }
  }
  #par(def.par)
  #return(reslist)
}

##------------corplotswithin_index function--------------------------
corplotswithin_index <- function(cpue,fl) {
  #Written by jdo [22/08/05], based loosely on functions by rds
  def.par <- par(no.readonly = TRUE)# save default, for resetting...
  if(!inherits(cpue, "FLIndices"))
    stop("First input must be an 'FLIndices' object!")
  #win.graph() # windows(width=7, height=6)
  par(mfrow=c(2,3))
  counter <- 0
  dims  <- dimnames(cpue[[fl]]@index)
  nages <- length(as.vector(dims$age))
  yrange    <- cpue[[fl]]@range["minyear"]:cpue[[fl]]@range["maxyear"]
  arange    <- cpue[[fl]]@range["min"]:cpue[[fl]]@range["max"]
  cpu       <- cpue[[fl]]@index[as.character(arange),]
  #print(cpu)
  cpu <- log(cpu)
  #print(cpu)
  t.     <- cbind(rep(yrange,each=max(arange)-min(arange)+1), arange, cpu)
  tdf    <- as.data.frame(cbind(t., t.[,1]-t.[,2]))
  colnames(tdf) <- list("year", "age", "cpue", "cohort")
  tdf$cohort<-substr(as.character(tdf$cohort),3,4)
  #print(tdf)
  if(nages > 1) {
    for(i in arange[1:(length(arange)-1)]) {
      counter <- counter + 1
      # 			if(counter == 5) {
      #         #win.graph() # windows(width=7, height=6)
      # 				par(mfrow=c(2,2))
      # 				counter <- 0
      # 			}
      d1<-tdf[tdf$age==i,]$cpue[1:(length(tdf[tdf$age==i,]$cpue)-1)]
      d1t<-tdf[tdf$age==i,]$cpue[1:(length(tdf[tdf$age==i,]$cpue)-1)]
      #print("d1")
      #print(d1)
      d2<-tdf[tdf$age==i+1,]$cpue[-1]
      d2t<-tdf[tdf$age==i+1,]$cpue[-1]
      #print("d2")
      #print(d2)
      d3<-tdf[tdf$age==i,]$cohort[1:(length(tdf[tdf$age==i,]$cohort)-1)]
      #print("d3")
      #print(d3)
      d1.<-d1[is.finite(d1)&is.finite(d2)]
      d2.<-d2[is.finite(d1)&is.finite(d2)]
      d3.<-d3[is.finite(d1)&is.finite(d2)]
      d4.<-min(d1.)+0.2*(max(d1.)-min(d1.))
      d<-matrix(c(d1.,d2.),length(d1.),2)
      dtt<-matrix(c(d1t,d2t),length(d1t),2)
      #print(d)
      #print(dtt)
      newf <- data.frame(x = seq(min(d2.), max(d2.), (max(d2.) - min(d2.))/(length(d2.)-1)))
      predict(lm(d2.~d1.), newf, se.fit = TRUE)
      pred.w.plim <- predict(rlm(d2.~d1.), newf, interval="prediction")
      pred.w.clim <- predict(rlm(d2.~d1.), newf, interval="confidence")
      #matplot(new$x,cbind(pred.w.clim, pred.w.plim[,-1]),
      #       lty=c(1,2,2,3,3), type="l", ,xlab=paste("Log-numbers at age",i),
      #       ylab=paste("Log-numbers at age",i+1))			
      xl <- c(min(d[,1])-1,max(d[,1])+1)
      yl <- c(min(d[,2])-1,max(d[,2])+1) #c(min(pred.w.clim[,2],d2.)-1,max(pred.w.clim[,3],d2.)+1)
      plot(d,type="n", xlim=xl, ylim=yl, 
           xlab=paste("Log-numbers at age",i),
           ylab=paste("Log-numbers at age",i+1))
      title(main=cpue[[fl]]@name)
      text(d1.,d2.,d3.,col="black",cex=0.6)
      abline(lm(d2.~d1.), xlim=xl, ylim=yl)
      abline(rlm(d2.~d1.,maxit=50),lty=3, xlim=xl, ylim=yl)
      #print(length(d[,1]))
      #print()
      if(is.finite(dtt[length(dtt[,1]),1])& is.finite(dtt[length(dtt[,2]),2]))
      {
        text(d[length(d[,1]),1],d[length(d[,2]),2],label="[  ]",font=2,col="red",cex=1.0)
        # text(d[length(d[,1])-1,1],d[length(d[,2])-1,2],label="[  ]",font=2,col="orange",cex=1.0)
        # text(d[length(d[,1])-2,1],d[length(d[,2])-2,2],label="[  ]",font=2,col="blue",cex=1.0)
      }
      #			print(d)
      #			print(pred.w.clim[order(pred.w.clim$fit,pred.w.clim$lwr,pred.w.clim$upr),])
      #lines(sort(d[,1]),sort(pred.w.clim[,3]), xlim=xl, ylim=yl,lty=2)
      #lines(sort(d[,1]),sort(pred.w.clim[,2]), xlim=xl, ylim=yl,lty=2)
      
      lines(sort(d[,1]),sort(pred.w.plim[,3]), xlim=xl, ylim=yl,lty=3)
      lines(sort(d[,1]),sort(pred.w.plim[,2]), xlim=xl, ylim=yl,lty=3)
      text(d4.,y=max(yl)-0.5,label=paste("cor =",substr(as.character(corr(d)),1,5)),col=1)
    }
  }
  par(def.par)
}

corplotsbetween_index1 <- function(cpue,fl1,fl2) {
  #Written by jdo [22/08/05], based loosely on functions by rds
  def.par <- par(no.readonly = TRUE)# save default, for resetting...
  if(!inherits(cpue, "FLIndices"))
    stop("First input must be an 'FLIndices' object!")
  #win.graph() # windows(width=7, height=6)
  #par(mfrow=c(1,1))
  counter <- 0
  
  minage<-max(cpue[[fl1]]@range["min"],cpue[[fl2]]@range["min"])
  maxage<-min(cpue[[fl1]]@range["max"],cpue[[fl2]]@range["max"])
  minyr<-max(cpue[[fl1]]@range["minyear"],cpue[[fl2]]@range["minyear"])
  maxyr<-min(cpue[[fl1]]@range["maxyear"],cpue[[fl2]]@range["maxyear"])
  
  arange <- minage:maxage
  yrange <- minyr:maxyr
  cpu1 <- cpue[[fl1]]@index[as.character(arange),as.character(yrange)]
  cpu1 <- log(cpu1)
  cpu2 <- cpue[[fl2]]@index[as.character(arange),as.character(yrange)]
  cpu2 <- log(cpu2)
  t1     <- cbind(rep(yrange,each=max(arange)-min(arange)+1), arange, cpu1)
  t2     <- cbind(rep(yrange,each=max(arange)-min(arange)+1), arange, cpu2)
  tdf1    <- as.data.frame(cbind(t1, t1[,1]-t1[,2]))
  tdf2    <- as.data.frame(cbind(t2, t2[,1]-t2[,2]))
  colnames(tdf1) <- list("year", "age", "cpue", "cohort")
  colnames(tdf2) <- list("year", "age", "cpue", "cohort")
  tdf1$cohort<-substr(as.character(tdf1$cohort),3,4)
  tdf2$cohort<-substr(as.character(tdf2$cohort),3,4)
  #print("data")
  
  for(i in arange[1:length(arange)]) {
    counter <- counter + 1
    # 		if(counter == 5) {
    # 		  counter <- 0
    #       #win.graph() # windows(width=7, height=6)
    # 	    par(mfrow=c(2,2))
    # 		}
    d1<-tdf1[tdf1$age==i,]$cpue
    d2<-tdf2[tdf2$age==i,]$cpue
    d3<-tdf1[tdf1$age==i,]$cohort
    d1.<-d1[is.finite(d1)&is.finite(d2)]
    d2.<-d2[is.finite(d1)&is.finite(d2)]
    d3.<-d3[is.finite(d1)&is.finite(d2)]
    d4.<-min(d1.)+0.2*(max(d1.)-min(d1.))
    d<-cbind(d1.,d2.)
    
    #plot(d,type="n",xlab=paste("Log-numbers:",cpue[[fl1]]@name),ylab=paste("Log-numbers:",cpue[[fl2]]@name))
    #title(main=paste("Age",i))
    #text(d1.,d2.,d3.,col="black",cex=0.6)
    #abline(lm(d2.~d1.))
    #abline(rlm(d2.~d1.,maxit=50),lty=3)
    #text(d4.,y=max(d2.),label=paste("cor =",substr(as.character(corr(d)),1,5)),col=1)
    
    newf <- data.frame(x = seq(min(d2.), max(d2.), (max(d2.) - min(d2.))/(length(d2.)-1)))
    predict(lm(d2.~d1.), newf, se.fit = TRUE)
    pred.w.plim <- predict(lm(d2.~d1.), newf, interval="prediction")
    pred.w.clim <- predict(lm(d2.~d1.), newf, interval="confidence")
    #matplot(new$x,cbind(pred.w.clim, pred.w.plim[,-1]),
    #       lty=c(1,2,2,3,3), type="l", ,xlab=paste("Log-numbers at age",i),
    #       ylab=paste("Log-numbers at age",i+1))			
    xl <- c(min(d[,1])-1,max(d[,1])+1)
    yl <- c(min(d[,2])-1,max(d[,2])+1) #c(min(pred.w.clim[,2],d2.)-1,max(pred.w.clim[,3],d2.)+1)
    plot(d,type="n", xlim=xl, ylim=yl, 
         xlab=paste("Log-numbers:",cpue[[fl1]]@name),
         ylab=paste("Log-numbers:",cpue[[fl2]]@name))
    title(main=paste("Age",i))       
    text(d1.,d2.,d3.,col="black",cex=0.6)
    abline(lm(d2.~d1.), xlim=xl, ylim=yl)
    abline(rlm(d2.~d1.,maxit=50),lty=3, xlim=xl, ylim=yl)
    text(d[length(d[,1]),1],d[length(d[,2]),2],label="[  ]",font=2,col="red",cex=1.0)
    # text(d[length(d[,1])-1,1],d[length(d[,2])-1,2],label="[  ]",font=2,col="orange",cex=1.0)
    # text(d[length(d[,1])-2,1],d[length(d[,2])-2,2],label="[  ]",font=2,col="blue",cex=1.0)
    
    #		print(d)
    #		print(pred.w.clim[order(pred.w.clim$fit,pred.w.clim$lwr,pred.w.clim$upr),])
    #   lines(sort(d[,1]),sort(pred.w.clim[,3]), xlim=xl, ylim=yl,lty=2)
    #   lines(sort(d[,1]),sort(pred.w.clim[,2]), xlim=xl, ylim=yl,lty=2)
    lines(sort(d[,1]),sort(pred.w.plim[,3]), xlim=xl, ylim=yl,lty=3)
    lines(sort(d[,1]),sort(pred.w.plim[,2]), xlim=xl, ylim=yl,lty=3)
    text(d4.,y=max(yl),label=paste("cor =",substr(as.character(corr(d)),1,5)),col=1)
    
    
  }
  par(def.par)
}