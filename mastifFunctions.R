
getEcoReg <- function(tdata, 
                      rpath = "/Users/jimclark/makeMastOnJimClark/ecoRegionFiles/WWFecoRegions" ){
  
  
  #https://www.worldwildlife.org/publications/terrestrial-ecoregions-of-the-world
  # tdata has plot, lon, lat
  
  require(maptools)
  require(ggmap)
  require(rgeos)
  require(sp)
  require(rgdal)
  
  tdata$plot <- as.character(tdata$plot)
  
  pnow <- getwd()
  
  plots <- unique( tdata$plot )
  mm    <- match(plots, tdata$plot)
  
  dat   <- data.frame(Longitude = tdata$lon[mm], Latitude = tdata$lat[mm],
                      names = tdata$plot[mm])
  wkeep <- which(is.finite(dat$Longitude))
  dat <- dat[wkeep,]
  coordinates(dat) <- ~ Longitude + Latitude
  
  setwd( rpath )
  wwf <- readOGR(".", "wwf_terr_ecos")
  
  proj4string(dat) <- CRS("+proj=longlat")
  
  wwf3 <- spTransform(dat, proj4string(wwf))
  dat$ecoRegWWF <- over(wwf3, wwf)$ECO_NAME 
  dat$ecoCode <- over(wwf3, wwf)$eco_code 
  
  dat$ecoRegWWF <- shortEcoReg(dat$ecoRegWWF)
  ecoReg <- as.character(dat$ecoRegWWF)
  
  setwd(pnow)
  
  kk <- match(tdata$plot, plots)
  ecoReg[kk]
}


shortEcoReg <- function(ecoRegWWF){
  
  ecoRegWWF <- as.character(ecoRegWWF)
  ecoRegWWF <- .replaceString(ecoRegWWF, ' and ', '/')
  ecoRegWWF <- .replaceString(ecoRegWWF, 'Eastern ', 'E ')
  ecoRegWWF <- .replaceString(ecoRegWWF, 'Northeastern ', 'NE ')
  ecoRegWWF <- .replaceString(ecoRegWWF, 'Northwestern ', 'NW ')
  ecoRegWWF <- .replaceString(ecoRegWWF, 'Northwest ', 'NW ')
  ecoRegWWF <- .replaceString(ecoRegWWF, 'Western ', 'W ')
  ecoRegWWF <- .replaceString(ecoRegWWF, 'Southeastern ', 'SE ')
  ecoRegWWF <- .replaceString(ecoRegWWF, 'South ', 'S ')
  ecoRegWWF <- .replaceString(ecoRegWWF, 'Northern ', 'N ')
  ecoRegWWF <- .replaceString(ecoRegWWF, 'North ', 'N ')
  ecoRegWWF <- .replaceString(ecoRegWWF, 'Southern ', 'S ')
  ecoRegWWF <- .replaceString(ecoRegWWF, 'Middle ', 'Mid ')
  ecoRegWWF <- .replaceString(ecoRegWWF, 'Southwestern ', 'SW ')
  ecoRegWWF <- .replaceString(ecoRegWWF, 'Southwest ', 'SW ')
  #ecoRegWWF <- .replaceString(dat$ecoRegWWF, 'Central ', 'C ')
  
  ecoRegWWF
  
}

treesNearTraps <- function(tdata, xytree, xytrap, meters = 60){
  
  # retain only trees within meters of the closest seed trap
  
  tmp  <- nn2( xytrap[,c('x','y')], xytree[,c('x','y')], k = 1)
  wrow <- which(tmp[[2]] < meters)
  
  tr <- columnPaste(tdata$plot,tdata$tree)
  xy <- columnPaste(xytree$plot,xytree$tree)[wrow]
  
  tdata <- tdata[tr %in% xy,]
  
  list(treeData = tdata, xytree = xytree[wrow,])
}
  
  
gen4code <- function(xx, nn=4){   # needed for combineds runs
  
  #shorten genus name in genusSpecies string to nn characters
  
  FAC  <- FALSE
  if(is.factor(xx))FAC <- TRUE
  xx   <- as.character(xx)
  fc   <- sapply(gregexpr("[A-Z]",xx),'[',1)
  wf   <- which(fc > 5)
  if(length(wf) > 0){
    gen <- substr(xx[wf],1,4)
    spp <- substr(xx[wf],fc[wf],1000)
    xx[wf] <- columnPaste(gen,spp,'')
  }
  if(FAC)xx <- as.factor(xx)
  xx
}

.appendMatrix <- function(m1,m2,fill=NA,SORT = FALSE,asNumbers = FALSE){  
  
  # matches matrices by column names
  # asNumbers: if column heads are numbers and SORT, then sort numerically
  
  if(length(m1) == 0){
    if(is.matrix(m2)){
      m3 <- m2
    } else {
      m3 <- matrix(m2,nrow=1)
    }
    if( !is.null(names(m2)) )colnames(m3) <- names(m2)
    return(m3)
  }
  if(length(m2) == 0){
    if(!is.matrix(m1))m1 <- matrix(m1,nrow=1)
    return(m1)
  }
  if( is.vector(m1) | (length(m1) > 0 & !is.matrix(m1)) ){
    nn <- names(m1)
    if(is.null(nn))message('cannot append matrix without names')
    m1 <- matrix(m1,1)
    colnames(m1) <- nn
  }  
  if( is.vector(m2) | (length(m2) > 0 & !is.matrix(m2)) ){
    nn <- names(m2)
    if(is.null(nn))message('cannot append matrix without names')
    m2 <- matrix(m2,1)
    colnames(m2) <- nn
  }
  
  c1 <- colnames(m1)
  c2 <- colnames(m2)
  r1 <- rownames(m1)
  r2 <- rownames(m2)
  n1 <- nrow(m1)
  n2 <- nrow(m2)
  
  allc <-  unique( c(c1,c2) ) 
  if(SORT & !asNumbers)allc <- sort(allc)
  if(SORT & asNumbers){
    ac <- as.numeric(allc)
    allc <- as.character( sort(ac) )
  }
  
  nr <- n1 + n2
  nc <- length(allc)
  
  if(is.null(r1))r1 <- paste('r',c(1:n1),sep='-')
  if(is.null(r2))r2 <- paste('r',c((n1+1):nr),sep='-')
  new <- c(r1,r2)
  
  mat1 <- match(c1,allc)
  mat2 <- match(c2,allc)
  
  out <- matrix(fill,nr,nc)
  colnames(out) <- allc
  rownames(out) <- new
  
  out[1:n1,mat1] <- m1
  out[(n1+1):nr,mat2] <- m2
  out
}

buildSpecByPlot <- function(cnames, mat, plot){
  
  if(is.matrix(mat)){
    wg <- which(!cnames %in% rownames(mat))
    if(length(wg) > 0){
      mm <- matrix(0, length(wg), ncol(mat) )
      rownames(mm) <- cnames[wg]
      mat <- rbind(mat, mm)
    }
    
    if(!plot %in% colnames(mat)){
      mat <- cbind(mat, 0)
      colnames(mat)[ncol(mat)] <- plot
    }
    mat[cnames, plot] <- 1
    
  }else{
    mat <- matrix(1, length(cnames), ncol=1)
    rownames(mat) <- cnames
    colnames(mat)  <- plot
  }
  mat
}

dbetaBinom <- function(y, n, p, sd, log = FALSE){
  
  # sd is sqrt(var) of beta(p|a, b), not of betaBinomial
  # not normalized for n
  
  
  n[ n < y ] <- y[n < y]
  bb <- n*0
  
  ww <- which(sd == 0)
  if(length(ww) > 0)bb[ww] <- dbinom(y[ww], n[ww], p[ww], log=T)
  
  ww <- which(sd > 0)
  tiny <- 1e-5
  a    <- p[ww]^2/sd[ww]/sd[ww]*(1 - p[ww]) - p[ww]
  a[a < tiny] <- tiny
  b  <- a*(1/p[ww] - 1)
  bb[ww] <- lchoose(n[ww], y[ww]) + lbeta(y[ww] + a, n[ww] - y[ww] + b) - lbeta(a, b)
  
  if(log)return(bb)
  
  exp(bb)
}

myrmultinom <- function(size,p,ASVECTOR = FALSE){  
  
  # n multinomial r.v. for a n by ncol(p) matrix of probs
  # each row of p is a probability vector
  # size is one integer or a length-n vector of integers
  # if ASVECTOR  = TRUE all size == 1, returns a vector of columns, otherwise a matrix
  
 # p <- row2Mat(p)
  
  n     <- nrow(p)
  J     <- ncol(p)
  
  if(length(size) == 1)size <- rep(size,n)
  
  jord  <- sample(J,J)    #randomize order
  
  rs <- rowSums(p)
  ws <- which(rs != 1)
  if(length(ws) > 0){
    p[ws,] <- p[ws,]/rs[ws]
  }
  
  p <- p[,jord,drop = FALSE]
  
  sizej <- size
  sumj  <- rep(0,n)
  dpj   <- rep(1,n)
  pj    <- p
  wj    <- c(1:n)
  
  if(ASVECTOR){        #  only if all size == 1
    
    yy <- size*0
    
    for(j in 1:(J-1)){
      a     <- round(pj[wj,1],10)
      tmp  <- rbinom(length(wj),sizej[wj],a)
      yy[wj[tmp == 1]] <- j
      sumj[wj]  <- sumj[wj] + tmp
      sizej <- size - sumj                       # no. remaining to choose
      dpj   <- dpj - p[,j]                       # Pr for remainder
      pj    <- matrix(p[,c((j+1):J)]/dpj,nrow(p))
      wj    <- which(sumj < size,arr.ind  = TRUE) 
    }
    
    yy[yy == 0] <- J
    yy <- jord[yy]
    
    return(yy)
  }
  
  yy  <- matrix(0,n,J)
  
  for(j in 1:(J-1)){
    a     <- round(pj[wj,1],10)
    yy[wj,j] <- rbinom(length(wj),sizej[wj],a)
    sumj  <- sumj + yy[,j]
    sizej <- size - sumj                       # no. remaining to choose
    dpj   <- dpj - p[,j]                       # Pr for remainder
    pj    <- matrix(p[,c((j+1):J)]/dpj,nrow(p))
    wj    <- which(sumj < size,arr.ind  = TRUE) 
  }
  
  if(n == 1)yy[,J] <- size - sum(yy)
  if(n > 1) yy[,J] <- size - rowSums(yy)
  
  yy[,jord] <- yy
  yy
  
}

mastIDmatrix <- function(treeData, seedData, genus, 
                         specNames = NULL, seedNames = NULL, 
                         censMin = NULL, verbose, ngen = 4){
  
  # possible seed ID errors: seedNames counted where seedNames of species 
  #    is missing. A '1' in R matrix indicates a 
  # must supply either 'genus' or both 'specNames' and 'seedNames'
  # ngen: no. of characters to match in genus
  
  CENS <- FALSE
  
  specIn <- specNames
  seedIn <- seedNames
  
  plotT <- treeData$plot <- as.character(treeData$plot)
  plotS <- seedData$plot <- as.character(seedData$plot)
  treeData$species <- as.character(treeData$species)
  
  if(is.null(specNames)){
    trows  <- which( startsWith( treeData$species, substr(genus,1,4) ) )
    specNames <- sort(unique(treeData$species[trows]))
  }else{
    trows <- which(treeData$species %in% specNames)
  }
  if(length(trows) == 0)stop('\nspecNames not found in treeData\n')
  
  if(is.null(seedNames)){
    scols  <- which( startsWith( colnames(seedData), substr(genus,1,4) ) )
  }else{
    scols <- which(colnames(seedData) %in% seedNames)
  }
  if(length(scols) == 0)stop('\nseedNames not found in seedData\n')
  
  snames <- colnames(seedData)[scols]
  
  tplots <- plotT[trows]
  ws     <- which( seedData[, scols, drop = FALSE] > 0, arr.ind  = TRUE )
  splots <- plotS[seedData$plot[ws[,1]]]
  
  allPlots <- sort(unique( plotT ) )
  
  tdata <- treeData[trows,]
  sdata <- seedData[plotS %in% allPlots,c('plot','trap',snames)]
  
  # possible errors where seed type is missing
  spec <- sort(unique(tdata$species))
  specBySeed <- matrix(0, length(spec), length(snames))
  rownames(specBySeed) <- spec
  colnames(specBySeed) <- snames
  
  if(!is.null(censMin)){
    CENS <- TRUE
    scens <- colnames(censMin)[colnames(censMin) %in% seedNames]
    pcens <- columnSplit(rownames(censMin),'-')[,1]
  }
  
  UN <- FALSE
  gu <- grep('UNKN', snames)
  if(length(gu) > 0)UN <- TRUE
  uncol <- rep(0, length(spec))
  names(uncol) <- spec
  
  for(j in 1:length(allPlots)){
    
    ws <- which( sdata$plot == allPlots[j])
    wt <- which( tdata$plot == allPlots[j])
    sm <- suppressWarnings( 
      apply( sdata[ws, snames, drop = FALSE], 2, max, na.rm  = TRUE)
    )
    sj <- names( sm[sm > 0] )
      
    tj  <- names( table(tdata$species[wt]) )
    win <- which(!tj %in% sj)                          #species not in seedNames
    if(length(win) > 0)specBySeed[ tj[win],sj ] <- 1
    if(CENS){
      wc <- which( pcens == allPlots[j])
      cj <- colSums( censMin[wc, scens, drop = FALSE], na.rm  = TRUE )
      cj <- names( cj[cj > 0] )
      wic <- which(!tj %in% cj)       
      if(length(wic) > 0)specBySeed[ tj[wic],cj ] <- 1
    }
    if(UN)uncol[tj] <- 1
  }
  
  wj <- match(rownames(specBySeed), colnames(specBySeed))
  wf <- which(is.finite(wj))
  specBySeed[ cbind(wf,wj[wf]) ] <- 1
  if(UN)specBySeed[,gu] <- uncol
  
  # missing seedNames
  wc <- colnames(specBySeed)[colSums(specBySeed) == 0]
  if(length(wc) > 0){
    wmiss <- colnames(specBySeed)[wc]
    specBySeed <- specBySeed[,!colnames(specBySeed) %in% wc, drop = FALSE]
    sdata <- sdata[,!colnames(sdata) %in% wc, drop = FALSE]
    snames <- colnames(specBySeed)
  }
  
  # missing specNames, seedNames
  wc <- rownames(specBySeed)[rowSums(specBySeed) == 0]
  if(length(wc) > 0)specBySeed[wc,] <- 1 
  
  specNames <- rownames(specBySeed)
  seedNames <- colnames(specBySeed)
  
  wm <- which(!seedIn %in% seedNames)
  if(length(wm) > 0){
    sa <- paste0(seedIn[wm], collapse = ', ')
    if(verbose)cat( paste('\nseedNames could not be used:', sa, '\n') )
    seedNames <- seedNames[!seedNames %in% sa]
  }
  wm <- which(!specIn %in% specNames)
  if(length(wm) > 0){
    sa <- paste0(specIn[wm], collapse = ', ')
    if(verbose)cat( paste('\nspecNames could not be used:', sa, '\n') )
    specNames <- specNames[!specNames %in% sa]
  }
  
  specBySeed <- specBySeed[drop = FALSE,specNames,seedNames]
  
  list(R = specBySeed, specNames = specNames,
       seedNames = seedNames, seedData = sdata)
}

mastPriors <- function(file, specNames, code, genus = 'NULL'){
  
  # if not found in 'code' column can use 'genus', if supplied
  # if only 'genus', then a single vector returned
  # code columns in priorParameters.txt are 'code...'
  
  specNames <- as.character(specNames)
  
  ns <- length(specNames)
  
  if( endsWith(file, '.txt') )
     priorVals <- read.table( file, header  = TRUE, stringsAsFactors = F)
  if( endsWith(file, '.csv') )
    priorVals <- read.csv( file, header  = TRUE, stringsAsFactors = F)
  
  ccols <- grep('code',colnames(priorVals))
  wcol  <- which(colnames(priorVals) == code)
  mcols <- ccols[!ccols == wcol]
  if(length(mcols) > 0)priorVals <- priorVals[,-mcols]
  
  wr <- which(priorVals[,'genus'] == genus & is.na(priorVals[,code]))
  genRow <- priorVals[drop = FALSE,wr,]
  
  if( code == 'genus'){
    rownames(genRow) <- prr$genus
    return(genRow[,-ccols, drop = FALSE])
  }
  prr <- priorVals[ 1:ns, ]   # dummy values
  
  wrow <- match(specNames, priorVals[,code])
  wf   <- which(is.finite(wrow))
  if( length(wf) > 0){
    prr[wf,]   <- priorVals[wrow[wf],]
  }
  
  wn <- which(!is.finite(wrow))   # use genus for missing

  if(length(wn) > 0){
    wr <- which(priorVals[,'genus'] == genus & is.na(priorVals[,code]))
    if(length(wr) > 0){
      prr[wn,] <-  priorVals[drop = FALSE, wr,] 
      prr[wn, code] <- specNames[wn]
    }else{
      wd <- which(priorVals$genus == 'default')
      prr[wn,] <- priorVals[drop = FALSE, wd,] 
      prr[wn, code] <- specNames[wn]
      prr[, 'genus'] <- genus
    }
  }
  rownames(prr) <- specNames
  prr
}

buildSeedByPlot <- function(sdata, snames){
  
  plot  <- as.character(sdata$plot)
  plots <- sort(unique(plot))
  
  snames <- snames[snames %in% colnames(sdata)]
  ntype  <- length(snames)
  nseed  <- nrow(sdata)
  
  if(ntype == 0)return(numeric(0))
  
  seedCount <- as.matrix(sdata[,snames])
  #  colnames(seedCount) <- snames
  
  wm <- match(as.character(sdata$plot),plots)

  ii <- rep(wm, ntype)
  jj <- rep(1:ntype, each=nseed)
  totalSeed <- t( .myBy(as.vector(seedCount),ii,jj,fun='sum') )
  colnames(totalSeed) <- plots
  rownames(totalSeed) <- paste('seeds_',colnames(seedCount), sep='')
  
  totalSeed
}

mastClimate <- function( file, plots, years, months, FUN = 'mean', 
                         vname = character(0), lastYear = 2021){
  
  # return covariate for a vector of plots and years
  # vname variable name; special treatment of 'degDays'
  # plots and years are vectors of the same length
  # months is a vector in (1, 12)
  
  fform <- integer(0)
  
  tform <- grep('.csv',file)
  if(length(tform) == 0)ffrom <- grep('.txt', file)
  
  if(length(fform) == 0 & length(tform) == 0)
    stop('\nfile must be .csv or .txt format\n')

  if(length(tform) > 0)data <- read.csv(file, header  = TRUE, row.names=1)
  if(length(fform) > 0)data <- read.table(file, header  = TRUE, row.names=1)
  
  colnames(data) <- .replaceString(colnames(data),'X','')
  rownames(data) <- .fixNames(rownames(data), all  = TRUE, MODE='character')$fixed
  
  # precip must be positive
  PREC <- FALSE
  dmin <- min(data, na.rm  = TRUE)
  if(dmin > 0)PREC <- TRUE
  
  plots <- .fixNames(plots, all  = TRUE, MODE='character')$fixed
  
  # small number of plots, but many points
  
  allPlots <- sort(unique(plots))
  
  tmp <- columnSplit(colnames(data), '_')
  yr  <- as.numeric( tmp[,1] )
  mo  <- as.numeric( tmp[,2] )
  yd  <- sort(unique(yr))
  
  my <- max(years)
  n.ahead <- my - max(yd)
  nyr     <- length(yd) + max( c(0, n.ahead) )
  
  if(n.ahead > 0){
    
    lastyr <- max(yr)
    lastmo <- max(mo)
    yseq   <- (lastyr + 1):(lastyr + n.ahead)
    mseq   <- (lastmo + 1):(lastmo + 12)
    ymore  <- rep( yseq, each = 12 )
    mmore  <- rep( mseq, n.ahead )
    
    wk <- which(mseq > 12)             # data do not end in december
    if(length(wk) > 0){
      mseq[wk] <- mseq[wk] - 12
      mmore <- rep( mseq, n.ahead )
      ytmp <- ymore*0
      ytmp[1] <- ymore[1] - 1
      ytmp[mseq == 1] <- ymore[mseq == 1] 
      for(j in 2:length(ytmp)) if(ytmp[j] == 0)ytmp[j] <- ytmp[j-1] 
      ymore <- ytmp
    }

    yr <- c(yr, ymore)
    mo <- c(mo, mmore)
    
    dmore <- matrix(NA, nrow(data),length(ymore))
    colnames(dmore) <- paste(ymore,'_',mmore,sep='')
    data <- cbind(data, dmore)
  }
  
  yd <- sort(unique(yr))
  xx <- aa <- numeric(0)  # hold mean and monthly anomaly
  wy <- match(yr, yd)
  wm <- match(mo, 1:12)
  
  missing <- character(0)
  
  for(j in 1:length(allPlots)){
    
    wj <- which(rownames(data) == allPlots[j])
    
    if(length(wj) == 0){
      xx <- cbind(xx, rep(NA, nyr))
      missing <- c(missing, allPlots[j])
      stop( paste('plot',allPlots[j],'is missing from',file) )
      next
    }
    
    dj   <- data[wj,] 
    if(is.list(dj))dj <- unlist(dj)
    
    wf     <- which(is.finite(dj))
    nahead <- length(dj) - max(wf)
    dj     <- dj[1:max(wf)]
    
    if(nahead > 0){  # predict forward
      fitAR <-  arima(dj, order = c(2, 0, 0), 
                      seasonal = list(order = c(1, 0, 0), period = 12))
      ptmp  <- predict(fitAR, n.ahead = nahead)
      ptmp  <- ptmp$pred + rnorm(nahead, 0, ptmp$se) 
      dj    <- c(dj, signif(ptmp, 3))
    }
    
    mmat <- matrix(NA, nyr, 12)
    mmat[ cbind(wy, wm) ] <- dj
    
    monthMu   <- colMeans(mmat, na.rm=T)  # monthly means over all years
    monthAnom <- t(t(mmat) - monthMu)
    
    if( vname == 'degDays' ){
      mvec <- mmat
      mvec[ mvec < 0 ] <- 0
      mvec <- 30*rowSums(mvec, na.rm=T)
      avec <- mvec*0
    }else{
      mvec <- suppressWarnings(
        apply(mmat[,months, drop = FALSE], 1, FUN, na.rm  = TRUE)
      )
      avec <- suppressWarnings(
        apply(monthAnom[,months, drop = FALSE], 1, FUN, na.rm  = TRUE)
      )
    }
    xx   <- cbind(xx, mvec)
    aa   <- cbind(aa, avec)
  }
  
  
  if(length(missing) > 0){
    pmiss <- paste0( missing, collapse=', ')
    warning( paste('\nMissing plots in covariate file:\n', pmiss) )
  }
  
  colnames(xx) <- colnames(aa) <- allPlots
  xx[!is.finite(xx)] <- NA
  aa[!is.finite(aa)] <- NA
  
  if(PREC)xx[xx < 0] <- 0
  
  iy <- match(years, yd)
  ip <- match(plots, colnames(xx))
  
  yy   <- xx[ cbind(iy, ip) ]
  anom <- aa[ cbind(iy, ip) ]
  
  site <- signif( apply(xx, 2, mean, na.rm=T)[ip], 4)
  tc   <- paste0(num2Month(months), collapse='')
  
  xm <- signif(cbind(yy, site, anom), 4)
  colnames(xm) <- paste(vname, tc, c('','SiteMean','Anom'), sep='')
  
  xm[,3] <- xm[,1] - xm[,2]
  
  if(length(missing) > 0){
    miss <- unique( columnPaste(plots[-wfy], years[-wfy], '_') )
    pmiss <- paste0( miss, collapse=', ')
    warning( paste('\nMissing plot_years in covariate file:\n', pmiss) )
  }

  list(x = xm, missingPlots = missing)
}

num2Month <- function(monthNum){
  
  mNames   <- c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug',
                'Sep','Oct','Nov','Dec')
  mNames[monthNum]
}

month2Num <- function(monthName){
  
  match(monthName, num2Month(1:12) )
}

lowerFirstLetter <- function(xx){
  s <- unlist(strsplit(xx, " "))
  s <- paste(tolower(substring(s, 1, 1)), substring(s, 2),
             sep = "", collapse = " ")
  unlist(strsplit(s, " "))
}

upperFirstLetter <- function(xx, FIRSTONLY = F){
  
  # FIRSTONLY - only first letter of first word when a string has multiple words
  
  if(FIRSTONLY){
    return( paste(toupper(substring(xx, 1, 1)), substring(xx, 2),
            sep = "") )
  }
  s <- unlist(strsplit(xx, " "))
  s <- paste(toupper(substring(s, 1, 1)), substring(s, 2),
             sep = "", collapse = " ")
  unlist(strsplit(s, " "))
}

getSigFig <- function(x){
  length(gregexpr("[[:digit:]]", as.character(x))[[1]])
}

combineSpecies <- function(specVec, specNames, combineSpecs){
  
  if(is.null(combineSpecs))
    return( list(specNames = specNames, species = specVec) )
  
  if(!is.character(specVec))specVec <- as.character(specVec)
  
  if(!is.null(combineSpecs)){
    wf <- which(combineSpecs[,'from'] %in% as.character(specVec))
    if(length(wf) > 0){
      specs <- as.character(specVec)
      for(m in wf){
        ws <- which(specs == combineSpecs[m,'from'])
        specs[ws] <- combineSpecs[m,'to']
      }
  #    specVec <- as.factor(specs)
      specVec <- specs
    }
    specNames <- sort(unique(as.character(specVec)))
  }
  list(specNames = specNames, species = specVec)
}

combineSeedNames <- function(sdata, snames, cseed=NULL){
  
  if(is.null(cseed))
    return( list(seedNames = snames, seedData = sdata) )
  
  from <- cseed[,1]
  to   <- cseed[,2]
  
  c1 <- matrix( paste(from,'_min',sep=''), ncol=1)
  c2 <- matrix( paste(from,'_max',sep=''), ncol=1)
  
  cnew <- rbind( cbind(c1, to), cbind(c2, to) )
  
  cmb <- rbind( cbind(from, to) , cnew )
  
  wmm <- grep('_min_min',cmb[,1])
  if(length(wmm) > 0){
    cmb[wmm,1] <- .replaceString(cmb[wmm,1], '_min_min','_min')
  }
  wmm <- grep('_min_max',cmb[,1])
  if(length(wmm) > 0){
    cmb[wmm,1] <- .replaceString(cmb[wmm,1], '_min_max','_min')
  }
  cmb <- cmb[!duplicated(cmb[,1]),]
  cmb <- cmb[order(cmb[,1]),]
  
  smb1 <- paste(snames,'_min',sep='')
  smb2 <- paste(snames,'_max',sep='')
  smb  <- c(snames, smb1, smb2)
  
  wmm <- grep('_min_min',smb)
  if(length(wmm) > 0){
    smb[wmm] <- .replaceString(smb[wmm], '_min_min','_min')
  }
  wmm <- grep('_min_max',smb)
  if(length(wmm) > 0){
    smb[wmm] <- .replaceString(smb[wmm], '_min_max','_min')
  }
  smb <- smb[!duplicated(smb)]
  
  ww <- which(smb %in% cmb[,'from'])
  if(length(ww) == 0 | is.null(cmb))
    return(list(seedNames = snames, seedData = sdata))
  
  mm <- match(smb[ww], cmb[,'from'])
  
 # sall <- smb[-mm]
  
  for(k in 1:length(mm)){
    fromk <- cmb[mm[k],'from']
    tok   <- cmb[mm[k],'to']
    smb[ww[k]] <- tok
    if( !tok %in% colnames(sdata) & fromk %in% colnames(sdata) ){
      colnames(sdata)[colnames(sdata) == fromk] <- tok
      next
    }
    if( tok %in% colnames(sdata) & fromk %in% colnames(sdata) ){
      sdata[ ,tok ] <- sdata[ ,tok ] +
        sdata[ ,fromk ]
    }else{
      colnames(sdata)[colnames(sdata) == fromk] <- tok
    }
  }
  
  wm <- which(colnames(sdata) %in% cmb[,'from'])
  if(length(wm) > 0)sdata <- sdata[,-wm]
 
  seedNames <- snames[!duplicated(snames)]
  seedNames <- seedNames[seedNames %in% colnames(sdata)]

  list(seedNames = seedNames, seedData = sdata)
}

trimCharVec <- function(cvec, string=NULL, good=NULL, bad=NULL){
  
  #delete elements from character vector cvec:
  #    containing string
  #    in bad
  #    not in good
  
  if( !is.null(string) ){
    wm <- grep(string, cvec)
    if(length(wm) > 0)cvec <- cvec[-wm]
  }
  if( !is.null(good) ){
    wm <- which(!cvec %in% good)
    if(length(wm) > 0)cvec <- cvec[-wm]
  }
  if( !is.null(bad) ){
    wm <- which(cvec %in% bad)
    if(length(wm) > 0)cvec <- cvec[-wm]
  }
  
  cvec
}

.trimRows <- function(xmat1, xmat2, xcol, STOP = FALSE){
  
  # xmat1 will have rows trimmed to include xmat2, no row matching
  
  kword <- character(0)
  
  x1 <- xmat1[,xcol]
  x2 <- xmat2[,xcol]
  if(is.factor(x1) | is.factor(x2)){
    x1 <- as.character(x1)
    x2 <- as.character(x2)
  }
  
  wm <- match(x1, x2)
  wf <- which(is.finite(wm))
  if(length(wf) < length(wm)){
    if(STOP){
      stop(paste('\nduplicates in', xcol))
    }else{
      kword <- paste('Values have been trimmed in ', xcol,'.',sep='')
      xmat1 <- xmat1[wf,]
    }
  }
  list(mat1 = xmat1, words = kword)
}

cleanSeedData <- function(sdata, xytrap, seedNames, verbose){
  
  words <- character(0)
  sdata$plot <- as.character(sdata$plot)
  sdata$trap <- as.character(sdata$trap)
  
  if(is.factor(sdata$year))sdata$year <- factor2integer(sdata$year)
  
  sdata$trapID <- columnPaste(sdata$plot,sdata$trap)
  xytrap$trapID <- columnPaste(xytrap$plot,xytrap$trap)
  
  if(!'area' %in% colnames(sdata)){
    sdata$area <- 1
    kword <- 'area is missing from seedData.'
    warning(kword)
    words <- paste(words, kword )
  }
  if( is.character(sdata$area) ){
    kword <- 'seedData$area was coerced to numeric.'
    warning( paste('\nNote: ', kword, '\n', sep='') )
    words <- paste(words, kword)
    sdata$area <- as.numeric(sdata$area)
  }
  if(!'active' %in% colnames(sdata)){
    kword <- 'An $active column was added to seedData.'
    if(verbose)cat( paste("\nNote: ", kword, sep="") )
    sdata$active <- 1
  }
  if( is.character(sdata$active) ){
    kword <- 'seedData$active coerced to numeric.'
    warning( paste('\nNote: ', kword, '\n', sep='') )
    words <- paste(words, kword)
    sdata$active <- as.numeric(sdata$active)
  }
  if(max(sdata$active, na.rm  = TRUE) > 1 | min(sdata$active, na.rm  = TRUE) < 0){
    kword <- 'Some seedData$active are outside (0, 1).'
    warning( paste('\nNote: ', kword, '\n', sep='') )
    words <- paste(words, kword)
  }
  
  wm <- match(colnames(sdata), seedNames)
  wf <- which(is.finite(wm))
  colnames(sdata)[wf] <- .fixNames(colnames(sdata)[wf], all  = TRUE, MODE='character')$fixed
  seedNames <- .fixNames(seedNames, all  = TRUE, MODE='character')$fixed
  
  #censored columns
  seedNamesAll <- seedNames
  m1 <- grep('_min',colnames(sdata))
  if(length(m1) > 0)seedNamesAll <- c( seedNamesAll, colnames(sdata)[m1] )
  m2 <- grep('_max',colnames(sdata))
  if(length(m2) > 0)seedNamesAll <- c( seedNamesAll, colnames(sdata)[m2] )
  
  nseed <- length(seedNamesAll)
  
  #fill missing with zero count, zero active
  sdata  <- sdata[order(sdata$plot, sdata$trap, sdata$year),]
  trapID <- columnPaste(sdata$plot, sdata$trap)
  id     <- sort(unique(trapID))
  snew   <- numeric(0)
  
  for(k in 1:length(id)){
    
    wk <- which(trapID == id[k])
    pk <- sdata$plot[wk[1]]
    ak <- range(sdata$year[sdata$plot == pk])
    ak <- ak[1]:max(ak)
    yk <- sdata$year[wk]
    if( (max(ak) - min(ak) + 1) == length(yk) ){
      smat <- sdata[wk,c('plot','trap','year','area','active',seedNamesAll)]
    }else{
  
      sk <- min(ak):max(ak)
      nk <- length(sk)
      
      kmat <- data.frame(plot = rep(sdata$plot[wk[1]],nk),
                         trap = rep(sdata$trap[wk[1]],nk),
                         year = sk, area = rep(sdata$area[wk[1]],nk))
      cmat <- matrix(0, nk, 1+nseed)
      colnames(cmat) <- c('active', seedNamesAll)
      smat <- as.matrix(sdata[wk,c('active', seedNamesAll)])
      cmat[match(yk,sk),] <- smat
      smat <- cbind(kmat, cmat)
    }
    snew <- rbind(snew, smat)
  }
  sid <- columnPaste(snew$plot,snew$trap)
  sid <- columnPaste(sid,snew$year,'_')
  rownames(snew) <- sid
  
  sdata <- snew
  
  ccols <- c('plot','trap','year','area','active',seedNamesAll)
  wm <- which(!ccols %in% colnames(sdata))
  if(length(wm) > 0)
    stop(paste0('\nNMissing columns from seedData:\n', ccols[wm], collapse=', '))
  
  ccols <- c('plot','trap','x','y')
  wm <- which(!ccols %in% colnames(xytrap))
  if(length(wm) > 0)
    stop(paste0('\nMissing columns from xytrap:\n', ccols[wm], collapse=', '))
  
  wna <- which(is.na(sdata$active) )
  if(length(wna) > 0){
    kword <- 'Some values are undefined in seedData$active.'
    if(verbose)cat( paste('\nNote: ', kword, '\n', sep='') )
    words <- paste(words, kword)
    sdata$active[wna] <- 0
  }
  wna <- which(is.na(sdata$area))
  if(length(wna) > 0)
    stop('\nSome values undefined or zero in seedData$area\n')
  
  sdata  <- .fixNamesVector( c('plot','trap'), sdata, MODE='character')
  xytrap <- .fixNamesVector( c('plot','trap'), xytrap, MODE='character')
  if(is.character(sdata$year))sdata$year <- as.numeric(sdata$year)
  sdata$trapID  <- columnPaste(sdata$plot,sdata$trap)
  sdata$plotYr  <- columnPaste(sdata$plot,sdata$year,'_')
  xytrap$trapID <- columnPaste(xytrap$plot,xytrap$trap)
  
  ww <- which(!sdata$trapID %in% xytrap$trapID)
  
  if(length(ww) > 0){
    wp <- paste0( sort(unique(sdata$trapID[ww])), collapse=', ')
    if(verbose)cat( paste('\nRemoved seedData missing from xytrap: ',wp,'\n',sep='') )
    sdata <- sdata[-ww,]
  }
   
  sdata <- sdata[,c('plot','trap','trapID','year','plotYr',
                    'area','active',seedNamesAll)]
  
  xytrap <- .cleanRows(xytrap, 'trapID')
  xytrap <- xytrap[xytrap$trapID %in% sdata$trapID,]
  
  list(sdata = sdata, xytrap = xytrap, words = words)
}

cleanTreeData <- function(tdata, xytree, specNames){
  
  words <- character(0)
  
  ccols <- c('plot','tree','species','year','diam')
  wm <- which(!ccols %in% colnames(tdata))
  if(length(wm) > 0)
    stop(paste0('\nMissing columns from treeData:\n', ccols[wm], collapse=', '))
  
  tdata  <- .fixNamesVector( c('plot','tree','species'), tdata, MODE='character')
  if(is.character(tdata$year))tdata$year <- as.numeric(tdata$year)
  
  if(is.null(xytree)){       # bogus tree locations if missing
    kword <- ' xytree is missing from inputs.'
    warning(kword)
    words <- paste(words, kword )
    
    tid <- columnPaste(tdata$plot, tdata$tree)
    tid <- sort(unique(tid))
    pt  <- columnSplit(tid, '-')
    xytree <- data.frame(plot = pt[,1], tree = pt[,2], x = 0, y = 0)
  }
  xytree  <- .fixNamesVector( c('plot','tree'), xytree, MODE='character')
  
  tdata <- tdata[as.character(tdata$species) %in% specNames,]
  tdata$treeID  <- columnPaste(tdata$plot, tdata$tree)
  xytree$treeID <- columnPaste(xytree$plot, xytree$tree)
  xytree <- xytree[as.character(xytree$treeID) %in% as.character(tdata$treeID),]
  
  mm <- match( as.character(xytree$treeID), as.character(tdata$treeID) )
  xytree$species <- tdata$species[mm]
  
  ccols <- c('plot','tree','x','y')
  wm <- which(!ccols %in% colnames(xytree))
  if(length(wm) > 0)
    stop(paste0('\nMissing columns from xytree:\n', ccols[wm], collapse=', '))
  
  tdata$plotYr  <- columnPaste(tdata$plot, tdata$year,'_')
  
  tmp <- .trimRows(xytree, tdata, 'treeID')
  xytree <-  tmp$mat1
  words <- paste(words, tmp$words)
  
  if(!'species' %in% colnames(xytree)){
    wm <- match(xytree$treeID, tdata$treeID)
    xytree$species <- tdata$species[wm]
  }
  
  if( is.character(tdata$diam) ){
    kword <- ' treeData$diam coerced to numeric.'
    warning( paste('\nNote: ', kword, '\n', sep='') )
    words <- paste(words, kword)
    tdata$diam <- as.numeric(tdata$diam)
  }
  
  if(is.factor(tdata$year))tdata$year <- factor2integer(tdata$year)
  
  xytree <- .cleanRows(xytree, 'treeID')
  xytree <- .trimRows(xytree, tdata, 'treeID')[[1]]
 # tdata  <- .trimRows(tdata, xytree, 'treeID')[[1]]
  
  list(tdata = tdata, xytree = xytree, specNames = specNames,
       words = words)
}

cleanInputs <- function(inputs, beforeFirst = 4,
                        afterLast = 4, p = 0, verbose = FALSE){
  
  xytree <- xytrap <- sdata <- NULL
  
  SEEDDATA <- TRUE
  if( !'seedData' %in% names(inputs) )SEEDDATA <- FALSE
  
  SEEDCENSOR <- TREESONLY <- FALSE   # censored counts with 'seedNames_min', 'seedNames_max'
  censMin <- censMax <- NULL
  #minDiam <- 0
  words <- character(0)
  cropCols <- c('cropCount','cropMin','cropMax')
  
  tdata        <- inputs$treeData
  specNames    <- inputs$specNames
  combineSpecs <- inputs$combineSpecs
  tdata$plot   <- .fixNames(tdata$plot, all  = TRUE, MODE='character')$fixed
  tdata$year   <- as.numeric(tdata$year)
  tdata$year   <- factor2integer(tdata$year)
  tdata        <- tdata[as.character(tdata$species) %in% specNames,]
  tdata$tree   <- .fixNames(tdata$tree, all  = TRUE, MODE='character')$fixed
  tdata$treeID <- columnPaste(tdata$plot, tdata$tree)
  priorTable   <- inputs$priorTable
  
  if(!is.null(priorTable)){
    minDiam <- priorTable[tdata$species, 'minDiam']
    maxDiam <- priorTable[tdata$species, 'maxDiam']
    maxFec  <- priorTable[tdata$species, 'maxFec']
  }else{
    if(verbose){
      cat( '\nNote: missing priorTable, used defaults for minDiam, maxDiam, maxFec' )
    }
    minDiam <- rep(10, nrow(tdata))
    maxDiam <- rep(40, nrow(tdata))
    maxFec  <- rep(1e+8, nrow(tdata))
  }
  
  ww <- which(is.na(tdata$year))
  if(length(ww) > 0)
    stop('treeData$year has NAs')
  if(min(tdata$year) < 1700)
    stop('treeData$year < 1700')
  
  tdata  <- tdata[order(tdata$plot, tdata$tree, tdata$year),]
  
  plotInput <- sort(unique(tdata$plot))
  
  if(SEEDDATA){
    xytree$plot   <- .fixNames(xytree$plot, all  = TRUE, MODE='character')$fixed
    xytree        <- inputs$xytree
    xytree        <- xytree[order(xytree$plot, xytree$tree),]
    xytree$plot   <- .fixNames(xytree$plot, all  = TRUE, MODE='character')$fixed
    xytree$treeID <- columnPaste(xytree$plot, xytree$tree)
    
    sdata  <- inputs$seedData
    xytrap <- inputs$xytrap
    combineSeeds <- inputs$combineSeeds
    seedNames    <- inputs$seedNames
    censMin      <- inputs$censMin
    censMax      <- inputs$censMax
    sdata$year   <- as.numeric(sdata$year)
    sdata$year   <- factor2integer(sdata$year)
    plotInput    <- sort(unique(c(tdata$plot,sdata$plot)))
    xytrap$plot  <- .fixNames(xytrap$plot, all  = TRUE, MODE='character')$fixed
    sdata$plot   <- .fixNames(sdata$plot, all  = TRUE, MODE='character')$fixed
    
    ww <- which(is.na(xytree$x) | is.na(xytree$y))
    if(length(ww) > 0)xytree <- xytree[-ww,]
    
    ww <- which(is.na(xytrap$x) | is.na(xytrap$y))
    if(length(ww) > 0)xytrap <- xytrap[-ww,]
    
    ww <- which(!xytrap$plot %in% sdata$plot)
    if(length(ww) > 0){
      nop   <- unique( xytrap$plot[ww] )
      sdata <- sdata[!sdata$plot %in% nop,]
    }
  }
  
  pname <- .fixNames(plotInput, all  = TRUE, MODE='character')$fixed
  names(plotInput) <- pname
  
  specNames <- .fixNames(inputs$specNames, all  = TRUE, MODE='character')$fixed
  
  if(length(which(duplicated(c(specNames)))) > 0)stop('duplicate specNames')
  
  tmp <- combineSpecies(tdata$species, specNames, combineSpecs)
  tdata$species <- tmp$species
  specNames     <- tmp$specNames
  
  
  # must have traps or crop count
  ll <- !tdata$treeID %in% xytree$treeID
  ww <- which(ll)
  
  
  if( length(ww) > 0 ){
    
    wc <- tdata$treeID %in% xytree$treeID
    
    if( 'cropCount' %in% colnames(tdata) ){
      wc <- wc | is.finite(tdata$cropCount) 
   }
    if('cropMin' %in% colnames(tdata) ){
      wc <- wc | is.finite(tdata$cropMin)
    }
    
    wk <- which(!wc)
  
    if(length(wk) > 0){
      tdata <- tdata[-wk,]
    }
    wk <- which(!tdata$treeID %in% xytree$treeID)
    if(length(wk) > 0){
      tdata <- rbind(tdata[-wk,],tdata[wk,])   # TREESONLY moved to end
    }
  }
  
  if(SEEDDATA){
    
    sdata$year  <- as.numeric(sdata$year)
    sdata$plot  <- .fixNames(sdata$plot, all  = TRUE, MODE='character')$fixed
    xytrap$plot <- .fixNames(xytrap$plot, all  = TRUE, MODE='character')$fixed
    
    sdata  <- sdata[order(sdata$plot, sdata$trap, sdata$year),]
    xytrap <- xytrap[order(xytrap$plot, xytrap$trap),]
    
    # columns in seedNames or seedNames_min/max
    scols     <- c("plot","trap","year","area","active")
    countCols <- colnames(sdata)[!colnames(sdata) %in% scols]
    
    ck <- numeric(0)
    for(k in 1:length(seedNames)){
      wk <- grep(seedNames[k], countCols)
      ck <- c(ck, wk)
    }
    countCols <- countCols[sort(unique(ck))]
    
    smin <- paste(seedNames,'_min',sep='')
    smax <- paste(seedNames,'_max',sep='')
    
    countCols <- countCols[countCols %in% c(seedNames, smin, smax)]
    
    tmp <- combineSeedNames(sdata, countCols, combineSeeds)
    sdata     <- tmp$seedData
    seedNames <- tmp$seedNames
    nseed     <- length(seedNames)
    
    sdata <- sdata[,c(scols,seedNames)]
    
    seedPlots <- sort(unique(sdata$plot))
    treePlots <- sort(unique(tdata$plot))
    
    sdata$trapID <- columnPaste(sdata$plot,sdata$trap)
    rownames(sdata) <- NULL
    
    plotTrapYr <- columnPaste(sdata$trapID, sdata$year, '_')
    
    rnames <- columnPaste(sdata$trapID, sdata$year, '_')
    
    rownames(sdata) <- rnames
    
    if( is.null(censMin) ){   # if censMin not built yet
      sall <- c(seedNames, paste(seedNames, '_min',sep=''),
                paste(seedNames,'_max',sep=''))
      countCols <- countCols[countCols %in% sall] 
      sdata <- sdata[,c(scols, countCols)]
    }
    
    xytrap$plot   <- as.character(xytrap$plot)
    sdata$plot    <- as.character(sdata$plot)
    sdata$trapID  <- columnPaste(sdata$plot,sdata$trap)
    xytrap$trapID <- columnPaste(xytrap$plot,xytrap$trap)
    
    #only seedData on plots with trees
    psave     <- intersect(seedPlots, treePlots)
    if(length(psave) == 0)stop('tree and seed data not from same plots')
    sdata     <- sdata[sdata$plot %in% psave,]
    xytrap    <- xytrap[xytrap$plot %in% psave,]
    
    # treesOnly
    ww <- which(!tdata$plot %in% sdata$plot)    # not on a seed trap plot
    
    if(length(ww) > 0){
      
      # if no crop counts, retain only trapped plots
      
      wcrop <- which(cropCols %in% colnames(tdata))
      
      if( length(wcrop) == 0 ){
        tdata     <- tdata[tdata$plot %in% psave,]
      }else{                                  # there are crop counts
        
        # not on trapped plots, but finite crop counts: keep but move to bottom
        w1 <- !tdata$plot %in% sdata$plot
        w2 <- rep(FALSE, length(w1))
        for(m in wcrop){
          w2[ is.finite( tdata[,cropCols[m]] ) ] <- TRUE
        }
        wk <- which(w1 & w2)
        
    #    wk <- which(!tdata$plot %in% sdata$plot &
    #                (is.finite(tdata$cropCount) | is.finite(tdata$cropMin) | is.finite(tdata$cropMax)))
    #    if(length(wk) == 0){
    #      tdata     <- tdata[tdata$plot %in% psave,]
    #    }else{
    #      wm <- which(!tdata$plot %in% sdata$plot &
    #                    !is.finite(tdata$cropCount) & 
    #                    !is.finite(tdata$cropMin) & !is.finite(tdata$cropMax))
    #      if(length(wm) > 0){
    #        tdata <- tdata[-wm,]
    #      }else{
            if(length(wk) > 0)tdata <- rbind(tdata[-wk,],tdata[wk,])
            TREESONLY <- TRUE
    #      }
    #    }
      }
    }
    
    seedPlots <- sort(unique(sdata$plot))
    treePlots <- sort(unique(tdata$plot))
    
    #tdata cannot start more than p yr before or after sdata
    
    pb <- p + beforeFirst
    pl <- p + afterLast
    
    womit <- numeric(0)
    for(k in 1:length(seedPlots)){
      kp <- range(sdata$year[sdata$plot == seedPlots[k]])
      wt <- which(tdata$plot == seedPlots[k])
      wk <- which(tdata$year[wt] < (kp[1] - pb))
      if(length(wk) > 0)womit <- c(womit, wt[wk])
      wk <- which(tdata$year[wt] > (kp[2] + pl))
      if(length(wk) > 0)womit <- c(womit, wt[wk])
    }
    
    # trees with fecundity estimates
    if(TREESONLY){
      wf <- which(is.finite(tdata$cropCount))
      womit <- womit[!womit %in% wf]
    }
    
    if(length(womit) > 0)tdata <- tdata[-womit,]
    
    # censor inactive traps
    seedNames <- .fixNames(inputs$seedNames, all  = TRUE, MODE='character')$fixed
    
    ww <- which(!seedNames %in% colnames(sdata))
    if(length(ww) > 0)seedNames <- seedNames[-ww]
    
    #fix seed names in censored columns
    ntype <- length(seedNames)
    scens <- .multivarChainNames( c('min', 'max'), seedNames )
    wcens <- which(scens %in% colnames(sdata))   # data input as censored
    mcols <- integer(0)
    
    if(length(wcens) > 0){  # repair names, but leave '_min', '_max'
      m1cols <- grep('_min',colnames(sdata))
      ncc <- .replaceString(colnames(sdata)[m1cols],'_min','')
      mc <- .fixNames(ncc, all  = TRUE, MODE='character')$fixed
      colnames(sdata)[m1cols] <- paste(mc,'_min',sep='')
      
      m2cols <- grep('_max',colnames(sdata))
      ncc <- .replaceString(colnames(sdata)[m2cols],'_max','')
      mc <- .fixNames(ncc, all  = TRUE, MODE='character')$fixed
      colnames(sdata)[m2cols] <- paste(mc,'_max',sep='')
      mcols <- c(m1cols, m2cols)
      
      seedNames <- .fixNames(seedNames, all  = TRUE, MODE='character')$fixed
    }
    
    wcc <- c(1:ncol(sdata))
    if(length(mcols) > 0)wcc <- wcc[-mcols]
    colnames(sdata)[wcc] <- .fixNames(colnames(sdata)[wcc], all  = TRUE, 
                                      MODE='character')$fixed
    
    tmp <- cleanSeedData(sdata, xytrap, seedNames, verbose)
    sdata  <- tmp$sdata
    xytrap <- tmp$xytrap
    words  <- paste(words, tmp$words)
    
    scols <- c("plot","trap", "trapID","year","plotYr","area","active")
    countCols <- colnames(sdata)[!colnames(sdata) %in% scols]
    
    sdata <- sdata[,c(scols, countCols)]
    
    if(length(wcens) > 0){  # only if censMin/censMax not built yet
      
      # move _min, _max columns out of sdata to censMin, censMax, 
      # only for rows
      
      SEEDCENSOR <- TRUE
      
      scens <- paste(seedNames, '_min', sep='')
      tcens <- paste(seedNames, '_max', sep='')
      wcens <- which(scens %in% colnames(sdata))   # data input as censored
      scens <- scens[wcens]
      wcens <- which(tcens %in% colnames(sdata))   # data input as censored
      tcens <- tcens[wcens]
      
      stypes <- columnSplit(scens,'_min')
      ttypes <- columnSplit(tcens,'_max')
      stypes <- sort(unique(c(stypes, ttypes, seedNames))) #seedNames from censored columns
      
      slo <- matrix(0, nrow(sdata), length(stypes))  
      colnames(slo) <- stypes
      snew <- shi <- slo 
      
      # censored rows:
      ms <- match(stypes, seedNames)
      wf <- which(is.finite(ms))
      
      snew[,stypes[ms[wf]]] <- as.matrix(sdata[,stypes[ms[wf]]])
      
      # uncensored rows
      scensName <- .replaceString(scens, '_min', '')
      slo[,scensName] <- as.matrix(sdata[,scens])
      shi[,scensName] <- as.matrix(sdata[,tcens])
      
      wnot  <- sort(unique(which(is.na(slo),arr.ind  = TRUE)[,1]))
      wcens <- sort(unique(which(is.na(snew),arr.ind  = TRUE)[,1]))
      
      sinit <- round( (slo[wcens,] + shi[wcens,])/2 )
      sinit[!is.finite(sinit)] <- slo[wcens][!is.finite(sinit)]
      sinit[!is.finite(sinit)] <- 0
      
      snew[wcens,] <- sinit
      
      snew[is.na(snew)] <- 0
      
      censMin <- cbind(wcens, slo[wcens,])
      censMax <- cbind(wcens, shi[wcens,])
      colnames(censMin)[1] <- colnames(censMax)[1] <- 'srow'
      censMin[is.na(censMin)] <- 0
      censMax[is.na(censMax)] <- censMin[is.na(censMax)]
      
      colnames(censMin) <- colnames(censMax) <- 
        .fixNames(colnames(censMin), all  = TRUE, MODE='character')$fixed
      
      sdata <- cbind(sdata[,scols], snew)
      
      sdata$plot <- .fixNames(sdata$plot, all  = TRUE, MODE='character')$fixed
      srr <- columnPaste(sdata$plot, sdata$trap)
      srr <- columnPaste(srr,sdata$year)
      rownames(sdata) <- srr
      rownames(censMin) <- rownames(censMax) <- rownames(sdata)[wcens]
    }
    
    srow <- which(sdata$active < 1)
    rrow <- which(is.na(sdata[,seedNames]),arr.ind  = TRUE)
    if(is.matrix(rrow))rrow <- rrow[,1]
    
    srow <- sort(unique(c(srow,rrow)))
    
    if(length(censMin) > 0){
      mm   <- match(rownames(censMin),rownames(sdata))
      srow <- srow[!srow %in% mm]
    }
    
    nseed <- length(seedNames)
    
    if(length(srow) > 0){ # only if not already in censMin/censMax
      
      # seedNames identified in plotYr
      
      cmin <- sdata[srow,seedNames, drop = FALSE]
      cmin[is.na(cmin)] <- 0
      
      sdata[srow,seedNames] <- cmin
      
      if(length(censMin) == 0){
        censMin <- cmin
        censMax <- cmin*0 + Inf
      }else{
        wnew <- which(!colnames(cmin) %in% colnames(censMin))
        if(length(wnew) > 0){
          newmat <- matrix(0, nrow(censMin), length(wnew))
          colnames(newmat) <- colnames(cmin)[wnew]
          censMin <- cbind(censMin, newmat)
          newmat <- newmat + Inf
          censMax <- cbind(censMax, newmat)
        }
        wnew <- which(!colnames(censMin) %in% colnames(cmin))  
        snn <- seedNames[seedNames %in% colnames(censMin)]
        cmax <- cmin
        cmax <- cmax + Inf
        
        censMin <- rbind(censMin[, snn, drop = FALSE], cmin[,snn, drop = FALSE])
        censMax <- rbind(censMax[, snn], drop = FALSE, cmax[,snn, drop = FALSE])
      }
    }
    sdata$trapID <- columnPaste(sdata$plot, sdata$trap)
  }
  
  plots     <- sort(unique(tdata$plot))
  
  if(SEEDDATA){
    
    tmp <- cleanTreeData(tdata, xytree, specNames)
    tdata     <- tmp$tdata
    xytree    <- tmp$xytree
    specNames <- tmp$specNames
    words     <- paste(words, tmp$words)
    plots     <- sort(unique(tdata$plot))
    
    wseed <- which(sdata$plot %in% plots)
    sdata  <- sdata[wseed,]
    xytrap <- xytrap[xytrap$plot %in% plots,]
    if(!is.null(censMin)){
      tmp <- trimCens(sdata, censMin, censMax)
      censMin <- tmp$censMin
      censMax <- tmp$censMax
    }
    
    # check coordinates
    ntt    <- table(xytree$plot)
    utt    <- tapply(xytree$x, xytree$plot, range)
    tplots <- .fixNames(names(utt))$fixed
    metersX <- matrix( round(unlist(utt)), ncol=2, byrow  = TRUE )
    colnames(metersX) <- c('minX', 'maxX')
    dx <- apply(metersX,1,diff)
    
    metersY <- matrix( round(unlist(tapply(xytree$y, xytree$plot, range))), ncol=2,
                       byrow  = TRUE )
    colnames(metersY) <- c('minY', 'maxY')
    dy <- apply(metersY,1,diff)
    
    ha <- round(apply(metersX,1,diff)*apply(metersY,1,diff)/10000, 2)
    metersX <- cbind(metersX, dx)
    metersY <- cbind(metersY, dy)
    
    
    mmm <- cbind(metersX, metersY)
    rownames(mmm) <- tplots
    mmm <- cbind(ntt[tplots],mmm)
    colnames(mmm)[1] <- 'trees'
    
    if(verbose){
      if(max(ha) > 100)cat('\nNote: plot area from xytree > 100 ha? See below:')
      cat('\n\nSpatial range for trees (dx, dy): \n')
      print(cbind(mmm))
    }
    
  }else{
    
    tdata$plotYr  <- columnPaste(tdata$plot, tdata$year,'_')
  }
  
  if(TREESONLY){
    wp <- which(!plots %in% rownames(mmm))
    if(length(wp) > 0){
      pw <- paste0(plots[wp], collapse=', ')
      pw <- paste('Trees not on trapped plots: ',pw,'.',sep='')
      words <- paste(words, pw)
      if(verbose)cat(paste('\n',pw,'\n',sep=''))
    }
  }
  
  if(SEEDDATA){
    
    ntt     <- table(xytrap$plot)
    splots  <- names(ntt)
    metersX <- matrix( round(unlist(tapply(xytrap$x, xytrap$plot, range))), ncol=2,
                       byrow  = TRUE )
    colnames(metersX) <- c('minX', 'maxX')
    dx <- apply(metersX,1,diff)
    
    metersY <- matrix( round(unlist(tapply(xytrap$y, xytrap$plot, range))), ncol=2,
                       byrow  = TRUE )
    colnames(metersY) <- c('minY', 'maxY')
    dy <- apply(metersY,1,diff)
    
    ha <- round( apply(metersX,1,diff)*apply(metersY,1,diff)/10000, 2 )
    metersX <- cbind(metersX, dx)
    metersY <- cbind(metersY, dy)
    
    mmm <- cbind(metersX, metersY)
    
    rownames(mmm) <- splots
    mmm <- cbind(ntt[splots],mmm)
    colnames(mmm)[1] <- 'traps'
    
    if(verbose){
      if(max(ha) > 100)cat('\nNote: plot area from xytrap > 100 ha? See below:')
      if(min(ha) < .01)cat('\nNote: plot area from xytrap < 1/100 ha? See below:')
      cat('\n\nSpatial range for trap data in meters: \n')
      print(mmm)
      
      if(!is.null(censMin)){
        stab <- table(sdata$plot)
        ttab <- matrix(0, 3, length(stab) )
        colnames(ttab) <- names(stab)
        ttab[2, ] <- stab
        mm   <- match(rownames(censMin),rownames(sdata))
        ctab <- table(sdata$plot[mm])
        ttab[1,names(ctab)] <- ctab
        ttab[3,] <- round(ttab[1,]/ttab[2,],3)
        rownames(ttab) <- c('censored','total','fraction')
        cat('\n\nCensored seed collections: \n')
        print(t(ttab))
      }
    }
  }
  
  # diameter
  rdiam <- range(tdata$diam,na.rm  = TRUE)
  if(rdiam[1] <= 0){
    ww <- which(tdata$diam <= 0)
    if(length(ww) > 0)tdata <- tdata[-ww,]
    kword <- ' Removed diameters <= 0.'
    if(verbose)cat( paste('\nNote: ', kword, '\n', sep='') )
    words <- paste(words, kword)
  }
  if(rdiam[2] > 500 & verbose)cat('\nNote: some diameters > 5 m')
  
  wna <- which(is.na(tdata$diam))
  if(length(wna) > 0){
    mm <- paste0( unique(tdata$treeID[wna]), collapse=', ')
    tdata <- tdata[-wna,]
    kword <- paste(' Removed treeData with missing diam:\n ', mm,sep='')
    if(verbose)cat('\nNote: ', kword,'\n')
    words <- paste(words, kword)
  }
  
  minDiam <- priorTable[tdata$species, 'minDiam']
  maxDiam <- priorTable[tdata$species, 'maxDiam']
  maxFec  <- priorTable[tdata$species, 'maxFec']
  
  if(is.null(minDiam))minDiam <- 5
  if(is.null(maxDiam))maxDiam <- 50
  if(is.null(maxFec))maxFec <- 1e+6
  
  
  wna <- which(tdata$diam > minDiam)
  if(length(wna) == 0){
    stop(paste('\nno trees > minDiam:\n ',minDiam,sep=''))
  }
  
  if(SEEDDATA)xytree <- .trimRows(xytree, tdata, 'treeID')[[1]]
  
  if('repr' %in% names(tdata)){
    rr <- suppressWarnings( range(tdata$repr,na.rm  = TRUE) )
    
    if(rr[1] == 1 & sum(rr, na.rm  = TRUE) == nrow(tdata)){
      kword <- ' All trees declared to be reproductive in treeData$repr.'
      warning(paste('\nNote: ', kword, '\n', sep='') )
      words <- paste( words, kword )
    }
    if(rr[2] == 0 & sum(rr, na.rm  = TRUE) == nrow(tdata)){
      kword <- ' All trees declared to be immature in treeData$repr.'
      warning(paste('\nNote: ', kword, '\n', sep='') )
      words <- paste( words, kword )
      tdata$repr[tdata$diam < minDiam] <- NA
    }
  }
  
  plots <- sort(unique(as.character(tdata$plot)))
  
  EXTEND <- FALSE
  if(EXTEND){
    for(j in 1:length(plots)){
      
      ts <- unique(tdata$year[as.character(tdata$plot) == plots[j]])
      rs <- ts
      if(SEEDDATA){
        ss <- unique(sdata$year[as.character(sdata$plot) == plots[j]])
        rs <- ss
      }
      
      rs <- rs[rs <= (max(ts) + afterLast)]
      rs <- rs[rs >= (min(ts) - beforeFirst)]
      ts <- ts[ts %in% rs]     # change for tree obs outside seed obs yrs
      
      pyr <- (min(rs)-p):(max(rs)+p)
      
      if(SEEDDATA){
        wj <- which(as.character(sdata$plot) == plots[j] & 
                      !sdata$year %in% pyr)
        if(length(wj) > 0)sdata <- sdata[-wj,]
      }
      wi <- which(as.character(tdata$plot) == plots[j])
      wj <- which(as.character(tdata$plot) == plots[j] & !tdata$year %in% pyr)
      
      if(length(wj) > 0){
        if(length(wi) == length(wj)){   #census does not fall in seed trap interval
          yi <- sort(unique(tdata$year[wj]))
          dd <- outer(X = yi, Y = pyr, function(X,Y) (X - Y)^2) #closest trap yr
          dy <- pyr[ apply(dd, 1, which.min) ]
          di <- match(tdata$year[wj], yi)
          tdata$year[wj] <- dy[di]
          wj <- which(as.character(tdata$plot) == plots[j] & !tdata$year %in% pyr)
        }
        if(TREESONLY){
          wc <- which(is.finite(tdata$cropFraction))
          wj <- wj[!wj %in% wc]
        }
        if(length(wj) > 0)tdata <- tdata[-wj,]
      }
    }
  }
  
  specNames <- sort(unique(tdata$species))
  plots     <- sort(unique(as.character(tdata$plot)))
  
  # duplicated tree years
  ty <- with(tdata, table(treeID, year) )
  treeYr <- columnPaste(tdata$treeID, tdata$year,'_')
  
  if(max(ty) > 1){
    wm <- which(ty > 1,arr.ind  = TRUE)
    cw <- unique(rownames(wm))
    cy <- paste0( cw , collapse=', ')
    
    wn <- which(!duplicated(treeYr))
    tdata <- tdata[wn,]
    
    kword <- paste( ' Removed treeData$tree with duplicate years:', cy )
    if(verbose)cat(paste('\nNote: ', kword, '\n', sep='') )
    words <- paste( words, kword )
    words <- paste(words, kword)
 #   tdata <- tdata[!as.character(tdata$treeID) %in% cw,]
  }
  
  
  # after tdata$year extended by p, beforeFirst, afterLast, 
  # remove traps beyond range(tdata$year) by plot
  # don't remove tree years between range(tdata$year), they will be interpolated
  # sdata <- trimPlotYr(tdata, sdata, beforeFirst, afterLast, p)
  
  years <- range( tdata$year ) 
  if(SEEDDATA){
    tmp <- .trimRows(xytree, tdata, 'treeID')
    xytree <- tmp$mat1
    words  <- paste(words, tmp$words)
    years <- range(c(years, sdata$year))
  }
  years <- min(years):max(years)
  
  # too rare
  tid     <- tdata[!duplicated(tdata$treeID),]
  specTab <- table(tid$species)
  wna     <- which(specTab < 5)
  
  if(length(wna) == length(specNames))
    stop( paste('All species too rare: ', paste0(names(wna),collapse=', '), sep='') )
  
  if(length(wna) > 0){               # species too rare
    
    bad   <- names(specTab)[wna]
    ww    <- which(!tdata$species %in% bad)
    tdata <- tdata[ww,]
    
    tmp <- .trimRows(xytree, tdata, 'treeID')
    xytree <- tmp[[1]]
    words  <- paste(words, tmp[[2]])
    
    rn     <- paste0(names(specTab)[wna],collapse=', ')
    specNames <- specNames[!specNames %in% names(specTab)[wna]]
    nspec  <- length(specNames)
    
    if(verbose){
      cat('\nNote: too rare, removed:\n')
      print( rn )
    }
    
    if(SEEDDATA){
      
      wukn <- grep('UNKN',bad)   # do not remove UNKN type
      if(length(wukn) > 0)bad <- bad[-wukn]
      
      badFruit <- which(colnames(sdata) %in% paste(bad, 'fruit', sep=''))
      badCones <- which(colnames(sdata) %in% paste(bad, 'cones', sep=''))
      
      bad <- c(bad, colnames(sdata)[c(badFruit, badCones)])
      
      wz <- which(bad %in% colnames(sdata))
      if(length(wz) > 0){
        for(b in wz){
          wb <- which(colnames(sdata) == bad[b])
          sdata <- sdata[,-wb]
          if(SEEDCENSOR){
            wb <- which(colnames(censMin) == bad[b])
            censMin <- censMin[,-wb, drop = FALSE]
            censMax <- censMax[,-wb, drop = FALSE]
          }
        }
        seedNames <- seedNames[!seedNames %in% bad]
      }
    }
    
    ptab <- table(tdata$plot)
    w0   <- which(ptab == 0)
    
    if(length(w0) > 0){
      pmiss <- names(ptab)[w0]
      wm <- which(as.character(xytree$plot) %in% pmiss)
      if(length(wm) > 0)xytree <- xytree[-wm,]
      
      if(SEEDDATA){
        wm <- which(as.character(sdata$plot) %in% pmiss)
        if(length(wm) > 0)sdata <- sdata[-wm,]
        
        wm <- which(as.character(xytrap$plot) %in% pmiss)
        if(length(wm) > 0)xytrap <- xytrap[-wm,]
      }
      plots  <- plots[!plots %in% pmiss]
    }
  }
  
  # remove plots where there are no trees
  plotRm <- character(0)
  
  for(j in plots){
    t2 <- which(as.character(tdata$plot) == j)
    if(length(t2) == 0){
      plotRm <- c(plotRm,j)
      next
    }
  }
  
  if(length(plotRm) > 0){
    plots  <- plots[!plots %in% plotRm]
    wr     <- which(tdata$plot %in% plotRm)
    if(length(wr) > 0)tdata  <- tdata[-wr,]
    wr     <- which(xytree$plot %in% plotRm)
    if(length(wr) > 0)xytree  <- xytree[-wr,]
    if(SEEDDATA){
      wr     <- which(sdata$plot %in% plotRm)
      if(length(wr) > 0)sdata  <- sdata[-wr,]
      wr     <- which(xytrap$plot %in% plotRm)
      if(length(wr) > 0)xytrap <- xytrap[-wr,]
    }
  }
  
  # retain trees in plot-years that have seed traps
  
  plotYears <- sort(unique(as.character(tdata$plotYr)))
  
  if(SEEDDATA){
    sdata        <- sdata[sdata$plot %in% plots,]
    sdata$plotYr <- columnPaste(sdata$plot, sdata$year,'_')
    xytrap       <- xytrap[xytrap$plot %in% plots,]
    plotYears <- sort(unique(c(as.character(sdata$plotYr), 
                               as.character(tdata$plotYr))))
    sdata$plotyr <- match( sdata$plotYr, plotYears )
    
    if(!is.null(censMin)){
      tmp <- trimCens(sdata, censMin, censMax)
      censMin <- tmp$censMin
      censMax <- tmp$censMax
    }
  }
  
  
  tdata$plotyr <- match( tdata$plotYr, plotYears )
  specNames    <- sort(unique(as.character(tdata$species)))
  plots        <- sort(unique(as.character(tdata$plot)))
  
  # seedNames, specNames
  
  if(SEEDDATA){
    
    xytree    <- xytree[xytree$treeID %in% tdata$treeID,]
    
    gg <- grep('UNKN',seedNames)
    
    if(length(gg) > 0){
      kword <- paste0( 'The unknown seed type is ',seedNames[gg],'.',sep=' ')
      if(verbose)cat( paste('\nNote: ', kword, '\n', sep='' ) )
      words <- paste(words, kword)
    }
    if(length(gg) > 1)
      stop( 'Only one seedName can have "UNKN" in name')
    
    ww <- which(!specNames %in% seedNames)
    if(length(ww) > 0){
      if(length(gg) == 0){
        kword <- ' There is no "UNKN" in seedNames'
        if(verbose)cat( '\nNote: ', kword, '\n' )
        words <- paste(words, kword)
      }
    }
    
    ww <- which(!seedNames %in% specNames)
    
    if(length(ww) > 0){
      
      wg <- which(!ww %in% gg)
      
      if(length(wg) > 0){
        
        missName <- seedNames[ww[wg]]
        mname    <- paste0(missName, collapse=', ')
        kword <- paste(' seedNames that are not in specNames and not "UNKN":\n', mname)
        
        if(verbose)cat( paste('\n', kword, '\n', sep='') )
        words <- paste(words, kword, '.', sep='')
        
        # appended specName
        ispec <- character(0)
        for(i in 1:length(specNames)){
          iss <- grep(specNames[i], missName)
          if(length(iss) > 0)ispec <- c(ispec, specNames[i])
        }
        
        # if not, add to UNKN class
        mcol <- grep(missName[1], colnames(sdata))
        
        if(length(gg) > 0){   # there is an UNKN type
          if(length(mcol) > 0 & length(ispec) == 0){
            if(length(missName) > 1){
              for(k in 1:length(missName)){
                mcol <- c(mcol, grep(missName[k], colnames(sdata)))
              }
            }
            mcol <- unique(mcol)
            
            smc <- sdata[,mcol]
            if(length(mcol) > 1)smc <- rowSums(smc, na.rm  = TRUE)
            
            sdata[,seedNames[gg]] <- sdata[,seedNames[gg]] + smc
            sdata <- sdata[,-mcol]
            kword <- paste0(' Moved ', missName, ' to "UNKN" class',
                            collapse=', ')
            if(verbose)cat( paste('\n', kword, '\n', sep='') )
            words <- paste(words, kword)
          }
        }else{
          sdata <- sdata[,-mcol]  #there is no UNKN type
          #     seedNames <- seedNames[!seedNames %in% missName]
        }
        
        if(length(ispec) == 0)seedNames <- seedNames[-ww[wg]]
      }
    }
    sdata  <- sdata[order(sdata$plot,sdata$trap,sdata$year),]
    xytrap <- xytrap[order(xytrap$plot, xytrap$trap),]
    
    if(SEEDCENSOR | length(censMin) > 0){
      
      tmp <- trimCens(sdata, censMin, censMax)
      censMin <- tmp$censMin
      censMax <- tmp$censMax
      
      wc <- which( colnames(censMin) %in% colnames(sdata) )
      censMin <- censMin[,wc, drop = FALSE]
      censMax <- censMax[,wc, drop = FALSE]
    }
  }
  
  if(is.null(tdata$obs))tdata$obs <- 1
  

  inputs$treeData  <- tdata
  inputs$specNames <- specNames
  inputs$plotInput <- plotInput
  inputs$inwords   <- words
  inputs$TREESONLY <- TREESONLY
  
  if(SEEDDATA){
    inputs$xytrap    <- xytrap
    inputs$seedNames <- seedNames
    inputs$seedData  <- sdata
    inputs$xytree    <- xytree
    inputs$censMin <- censMin
    inputs$censMax <- censMax
  }

  inputs
}

trimCens <- function(sdata, censMin, censMax){
  
  mm <- match(rownames(censMin),rownames(sdata))
  wf <- which(is.finite(mm))
  censMin <- censMin[drop = FALSE, wf,]
  censMax <- censMax[drop = FALSE, wf,]
  list(censMin = censMin, censMax = censMax)
}

trimPlotYr <- function(tdata, sdata, beforeFirst, afterLast, p){
  
  tmp  <- table(tdata$plot, tdata$year)
  
  tmp[tmp > 1] <- 1
  tmat <- tmp*matrix(as.numeric(colnames(tmp)), nrow(tmp), ncol(tmp), byrow  = TRUE)
  tmat[tmat == 0] <- NA
  mm   <- apply(tmp*tmat, 1, range, na.rm  = TRUE)
  
  mm[1,] <- mm[1,] - max(c(beforeFirst, p))
  mm[2,] <- mm[2,] + max(c(afterLast, p))
  ty   <- character(0)
  for(j in 1:ncol(mm)){
    jy <- paste(colnames(mm)[j],mm[1,j]:mm[2,j],sep='_')
    ty <- c(ty, jy)
  }
  ww <- which(!sdata$plotYr %in% ty)
  if(length(ww) > 0){
    sdata <- sdata[-ww,]
  }
  sdata
}


addObsTrap <- function(tdata, sdata){
  
  plots <- sort(unique(tdata$plot))
  nplot <- length(plots)
  obsTrap <- rep(0, nrow(tdata))
  
  for(j in 1:nplot){
    wj  <- which(as.character(sdata$plot) == plots[j])
    if(length(wj) == 0)next
    sjj <- sdata[wj,]
    pyr <- min(sjj$year):max(sjj$year)
    wtt <- which(as.character(tdata$plot) == plots[j] & 
                   tdata$year %in% pyr)
    obsTrap[wtt] <- 1
  }
  obsTrap
}

mastFillCensus <- function(inputs, beforeFirst = 15,
                           afterLast = 15, p = 0, verbose = FALSE){
  
  # fill census with seed trap years, interpolated diameter
  # beforeFirst - assumed present for no. of years before tree first observed
  # afterLast   - no. yr after last observed
  # p       - if AR(p > 0) model, fills p yr before and after tree observed
  
  words <- character(0)
  SEEDCENSOR <- TREESONLY <- FILLED <- FALSE
  SEEDDATA   <- TRUE
  if(!'seedData' %in% names(inputs))SEEDDATA <- FALSE
  
  if('FILLED' %in% names(inputs))return(inputs)
  
  AR <- FALSE
  if(p > 0){
    AR <- TRUE
    kword <- paste(' This is an AR(', p, ') model. ', sep='')
    words <- paste(words, kword)
    if(verbose)cat( '\n', kword, '\n' )
  }
  
  inputs$treeData$year <- factor2integer(inputs$treeData$year)
  inputs$treeData$plot <- .fixNames(inputs$treeData$plot, all  = TRUE, 'character')$fixed
  
  ss <- tapply(inputs$treeData$year, 
               list(plot = inputs$treeData$plot), range, na.rm  = TRUE)
  sn <- names(ss)
  censusYr <- matrix( unlist(ss), ncol=2, byrow  = TRUE )
  rownames(censusYr) <- sn
  
  if(SEEDDATA){
    inputs$seedData$year <- factor2integer(inputs$seedData$year)
    inputs$seedData$plot <- .fixNames(inputs$seedData$plot, all  = TRUE, 'character')$fixed
    ss <- tapply(inputs$seedData$year, 
                 list(plot = inputs$seedData$plot),range)
    sn <- names(ss)
    trapYr <- matrix( unlist(ss), ncol=2, byrow  = TRUE )
    rownames(trapYr) <- sn
  }
  
  inputs    <- cleanInputs(inputs, beforeFirst, afterLast, p, verbose = verbose)
  specNames <- inputs$specNames
  seedNames <- inputs$seedNames  
  tdata     <- inputs$treeData   
  sdata     <- inputs$seedData  
  xytree    <- inputs$xytree     
  xytrap    <- inputs$xytrap
  censMin   <- inputs$censMin
  censMax   <- inputs$censMax
  words     <- c(words, inputs$inwords)
  TREESONLY <- inputs$TREESONLY
  rownames(tdata) <- columnPaste(tdata$treeID,tdata$year,'_')
  
  # rare species removed in cleanInputs, 
  #    now determine if there are still cropCounts/cropMin
  if('cropCount' %in% colnames(tdata)){
    rc <- suppressWarnings( max(tdata$cropCount, na.rm=T) )
    if(rc < 0){
      CONES <- FALSE
      tdata <- tdata[,!colnames(tdata) %in% 
                               c('cropCount','cropFraction','cropFractionSd')]
    }
    rc <- suppressWarnings( max(tdata$cropMin, na.rm=T) )
    if(rc < 0){
      CONES <- FALSE
      tdata <- tdata[,!colnames(tdata) %in% 
                       c('cropMin','cropMax')]
    }
  }
  
  if(!is.null(censMin))SEEDCENSOR <- TRUE
  

  years   <- range( tdata$year[tdata$obs == 1] )  
 # treeIDs <- sort( unique(as.character(tdata$treeID)) )
  plots   <- sort( unique(as.character(tdata$plot)) )
  nplot   <- length(plots)
  nyr     <- length(years)
  
  if(SEEDDATA){
    sdata$obs     <- 1
    tdata$obsTrap <- addObsTrap(tdata, sdata)
    years   <- range( c(years, sdata$year) )
    trapIDs <- sort( unique(as.character(sdata$trapID)) )
  }
  
  years <- years[1]:years[2]
  allYears <- (min(years) - p):(max(years) + p)
  
  rownames(tdata) <- columnPaste(tdata$treeID,tdata$year,sep='_')
  
  vtypes <- getVarType(colnames(tdata), tdata, i=tdata$treeID, j = tdata$year)
  
  if('repr' %in% names(tdata))vtypes$repr = 'ij'
  if('repMu' %in% names(tdata))vtypes$repMu = 'ij'
  if('repSd' %in% names(tdata))vtypes$repSd = 'ij'
  if('province' %in% names(tdata))vtypes$province = 'i'
  
  
  fecMin <- fecMax <- cropMin <- cropMax <- NULL
  
  
  if( 'cropMin' %in% names(tdata) ){ # if censored classes, then fecMin, fecMax already provided
    vtypes$cropMin = 'ij'
    cropMin <- tdata$cropMin
    names(cropMin) <- rownames(tdata)
    
    wcrop <- which( is.finite(cropMin) )
    if(!'fecMin' %in% colnames(tdata)){
      tdata$fecMin <- NA
    }
    mcrop <- which( !is.finite(tdata$fecMin[wcrop]) )
    
    if( 'seedTraits' %in% names(inputs) & length(mcrop) > 0 ){  
      fecMin <- tdata$cropMin[wcrop]*inputs$seedTraits[ tdata$species[wcrop], 'seedsPerFruit' ]
      tdata$fecMin[wcrop[mcrop]] <- fecMin[mcrop]
      names(fecMin) <- rownames(tdata)
      if(verbose)cat(paste('\nNote: cropMin values without fecMin--used seedsPerFruit\n') )
    }
  }
  if('cropMax' %in% names(tdata)){
    vtypes$cropMax = 'ij'
    cropMax <- tdata$cropMax
    names(cropMax) <- rownames(tdata)
    
    wcrop <- which( is.finite(cropMax) )
    if(!'fecMax' %in% colnames(tdata)){
      tdata$fecMax <- NA
    }
    mcrop <- which( !is.finite(tdata$fecMax[wcrop]) )
    
    if( 'seedTraits' %in% names(inputs) & length(mcrop) > 0 ){
      fecMax <- tdata$cropMax*inputs$seedTraits[ tdata$species, 'seedsPerFruit' ]
      fecMax[fecMax < 1] <- 1
      tdata$fecMax[wcrop[mcrop]] <- fecMax[mcrop]
      tdata$fecMax <- fecMax
      names(fecMax) <- rownames(tdata)
      if(verbose)cat(paste('\nNote: cropMax values without fecMax--used seedsPerFruit\n') )
    }
  }
  
  if('priorTable' %in% names(inputs)){
    if('fecMax' %in% colnames(priorTable)){
      mf <- inputs$priorTable[tdata$species,'maxFec']
      tdata$fecMax[ tdata$fecMax > mf ] <- mf[ tdata$fecMax > mf ]
    }
  }
  if('fecMax' %in% colnames(tdata)){
    vtypes$fecMax = 'ij'
    fecMax <- tdata$fecMax
    names(fecMax) <- rownames(tdata)
  }
  
  wnull <- which(is.null(vtypes))
  if(length(wnull) > 0){
    for(k in 1:length(wnull)) vtypes[[k]] <- 'ij'
  }
  
  # gap years for individuals
  
  ww   <- numeric(0)
  if(length(table(tdata$year)) > 1){
    ttab <- table(tdata$treeID, tdata$year)
    wmin <- apply(ttab, 1, which.max)
    wmax <- apply(ttab[,ncol(ttab):1], 1, which.max)
    wmax <- c(ncol(ttab):1)[wmax]
    wtot <- rowSums(ttab)
    ww   <- which(wmax > wmin & wtot < (wmax - wmin))
  }
  
  if(SEEDDATA){                                 # traps without tree data
    kk <- which(!as.character(sdata$plotYr) %in% 
                  as.character(tdata$plotYr) )
    ww <- sort( unique(c(ww,kk) ))
  }
    
  treeID  <- tdata$treeID
  treeIDs <- treeID[ !duplicated(treeID) ]
  
  if( length(ww) > 0 | AR ){  # fill missing tree years
    
    yrSeq <- years
    if(AR)yrSeq <- allYears
    ijFull <- numeric(0)
    nyr    <- length(yrSeq)
    
    mm <- match(as.character(tdata$treeID), treeIDs)
    ijIndex <- cbind( mm, match(tdata$year, yrSeq) )
    
    matFull <- matrix(0 , length(treeIDs), nyr)
    rownames(matFull) <- treeIDs
    matFull[ ijIndex ] <- 1
    
    # seed years by plot
    if(SEEDDATA){  
      
      kjIndex <- numeric(0)
      splot <- unique(sdata$plot)
      tplot <- columnSplit(treeIDs,'-')[,1]
      for(m in 1:length(splot)){
        wm <- which(tplot == splot[m])
        my <- range( sdata$year[ sdata$plot == splot[m] ] )
        my <- match(my, yrSeq)
        my <- my[1]:my[length(my)]
        km <- as.matrix( expand.grid(wm,my) )
        kjIndex <- rbind( kjIndex, km )
      }
      matFull[ kjIndex ] <- 1
    }
      
    # years for p lags
    ijp <- ijIndex
    ijp[,2] <- ijp[,2] + p
    ijp <- ijp[ ijp[,2] <= nyr, ]
    matFull[ ijp ] <- 1
    ijp <- kjIndex
    ijp[,2] <- ijp[,2] + p
    ijp <- ijp[ ijp[,2] <= nyr, ]
    matFull[ ijp ] <- 1
    
    ijp <- ijIndex
    ijp[,2] <- ijp[,2] - p
    ijp <- ijp[ ijp[,2] > 0, ]
    matFull[ ijp ] <- 1
    ijp <- kjIndex
    ijp[,2] <- ijp[,2] - p
    ijp <- ijp[ ijp[,2] > 0, ]
    matFull[ ijp ] <- 1
    
    ymat  <- matrix( yrSeq, length(treeIDs), nyr, byrow=T)
    first <- apply( matFull, 1, which.max )
    last  <- apply( ymat*matFull, 1, which.max )
    
    for(k in 1:nyr){
      wk <- which( k >= first & k <= last )
      matFull[wk,k] <- 1
    }
    
    ijFull <- which(matFull == 1, arr.ind=T)
    
    tdata$times <- match(tdata$year, allYears)
    if(SEEDDATA)sdata$times <- match(sdata$year, allYears)
    
    vtypes    <- getVarType(colnames(tdata), tdata, i=tdata$treeID, j = tdata$year) 
    vtypes$obs   <- 'ij'
  #  vtypes$times <- 'j'
    vtypes$diam  <- 'ij'
    vtypes$repr  <- 'ij'
    vtypes$obs   <- 'ij'
    vtypes$year  <- 'time'
    if('cropCount' %in% colnames(data))vtypes$cropCount <- 'none'
    if('cropFraction' %in% colnames(data))vtypes$cropFraction <- 'none'
    if('cropFractionSd' %in% colnames(data))vtypes$cropFractionSd <- 'none'
    if('fecMin' %in% colnames(data))vtypes$fecMin <- 'none'
    if('fecMax' %in% colnames(data))vtypes$fecMax <- 'none'
    
   # jtimes <- 1:nyr
    
    tmp <- fillMissing(variables = vtypes, data = tdata, icol = 'treeID', jcol = 'year', 
                       jtimes = allYears, ijFull = ijFull)
    data <- tmp$data
    
   # data$year <- allYears[ data$times ]
    data$plot      <- as.character(data$plot)
    data$tree      <- as.character(data$tree)
    rnames         <- columnPaste(data$treeID,data$year, '_')
    rownames(data) <- rnames
    data$diam      <- round(data$diam,1)
    
    if('province' %in% colnames(data)){
      ptab <- table(data$province,data$plot)
      wtab <- which(ptab > 0, arr.ind  = TRUE)
      wtab <- cbind(rownames(ptab)[wtab[,1]], colnames(ptab)[wtab[,2]])
      mm   <- match(data$plot,wtab[,2])
      data$province <- wtab[mm,1]
    }
    
    wm <- match(rownames(data), rownames(tdata))
    wf <- which(is.finite(wm))
    
    data$repr <- NA
    if('repr' %in% colnames(tdata))data$repr[wf] <- tdata$repr[wm[wf]]
    
    if('cropCount' %in% colnames(data)){
      
      data$cropCount <- data$cropFraction <- data$cropFractionSd <- NA
      
      data$cropCount[wf]    <- ceiling( tdata$cropCount[wm[wf]] )
      if('cropFraction' %in% colnames(tdata))data$cropFraction[wf] <- tdata$cropFraction[wm[wf]]

      if(!'cropFractionSd' %in% colnames(data))data$cropFractionSd <- NA
      
      if('cropFractionSd' %in% colnames(tdata)){
        data$cropFractionSd[wf] <- tdata$cropFractionSd[wm[wf]]
      }
      wmm <- which(is.na(data$cropFractionSd) & !is.na(data$cropFraction))
      if(length(wmm) > 0){
        data$cropFractionSd[wmm] <- .2*dbeta(data$cropFraction[wmm], .1,2) + 1e-3
      }
      data$repr[data$cropCount > 0] <- 1
    }
    
    tdata        <- data
    tdata$obs    <- as.numeric(as.character(tdata$obs))
    tdata$obs[is.na(tdata$obs)] <- 0
    tdata$plotYr <- columnPaste(tdata$plot,tdata$year,'_') 
    tdata$treeID <- as.character(tdata$treeID)
    
    if(SEEDDATA){
      xytree <- .trimRows(xytree, tdata, 'treeID')[[1]]   # trees having too few years
      
      ss <- tapply(sdata$year, list(plot = sdata$plot),range)    
      sn <- names(ss)
      trapYr <- matrix( unlist(ss), ncol=2, byrow  = TRUE )
      rownames(trapYr) <- sn
      
      tdata$obsTrap <- tdata$obs <- 0
      for(j in 1:nplot){
        kj <- as.character(tdata$plot) == plots[j]
        wj <- which(kj)
        if(!plots[j] %in% rownames(trapYr))next
        yr <- (trapYr[plots[j],1]):(trapYr[plots[j],2])
        tdata$obsTrap[wj] <- 0
        tdata$obsTrap[kj & (tdata$year %in% yr)] <- 1
        tdata$obs[kj & (tdata$year %in% yr)] <- 1
      }
      sdata <- trimPlotYr(tdata, sdata, beforeFirst, afterLast, p)
    }
    
    tdata$obs[is.finite(tdata$cropFraction)] <- 1
    tdata$obs[is.finite(tdata$cropMin)] <- 1
    
    if(SEEDCENSOR){
      tmp <- trimCens(sdata, censMin, censMax)
      censMin <- tmp$censMin
      censMax <- tmp$censMax
    }
    
    plotYears <- sort(unique(c(as.character(tdata$plotYr), 
                               as.character(sdata$plotYr))))
    tdata$plotyr <- match( as.character(tdata$plotYr), plotYears )
    if(SEEDDATA)sdata$plotyr <- match( as.character(sdata$plotYr), plotYears )
    
    rownames(tdata) <- columnPaste(tdata$treeID, tdata$year, '_')
    
    if(!'repr' %in% colnames(tdata))tdata$repr <- NA
    
    if('fecMin' %in% names(tdata) & !is.null(fecMin)){
      wm <- match(names(fecMin), rownames(tdata))
      wf <- which(is.finite(wm))
      tdata$fecMin <- NA
      tdata$fecMin[wm[wf]] <- fecMin[wf]
    }
    if('fecMax' %in% names(tdata) & !is.null(fecMax)){
      wm <- match(names(fecMax), rownames(tdata))
      wf <- which(is.finite(wm))
      tdata$fecMax <- NA
      tdata$fecMax[wm[wf]] <- fecMax[wf]
    }
    if('cropMin' %in% names(tdata)){
      wm <- match(names(cropMin), rownames(tdata))
      wf <- which(is.finite(wm))
      tdata$cropMin <- NA
      tdata$cropMin[wm[wf]] <- cropMin[wf]
      tdata$repr[tdata$cropMin > 0] <- 1
    }
    if('cropMax' %in% names(tdata)){
      wm <- match(names(cropMax), rownames(tdata))
      wf <- which(is.finite(wm))
      tdata$cropMax <- NA
      tdata$cropMax[wm[wf]] <- cropMax[wf]
      tdata$repr[tdata$cropMax > 0] <- 1
    }
  }
  
  
  allYears <- sort(unique(tdata$year))
  
  tdata$times <- match(tdata$year, allYears)
  sdata$times <- match(sdata$year, allYears)
  
  if(TREESONLY & SEEDDATA){
    ww <- which(!tdata$plot %in% sdata$plot)
    if(length(ww) > 0){
      tdata <- rbind(tdata[-ww,],tdata[ww,])
    }
  }
  
  attr(tdata, 'plag') <- p
  
  inputs$treeData  <- tdata
  inputs$xytree    <- xytree
  inputs$specNames <- specNames
  inputs$inwords   <- words
  inputs$TREESONLY <- TREESONLY
  
  if(SEEDDATA){
    inputs$seedData  <- sdata
    inputs$xytrap    <- xytrap
    inputs$seedNames <- seedNames
    inputs$censMin   <- censMin
    inputs$censMax   <- censMax
  }
  
  inputs$FILLED <- T
  
  inputs
}

vec2mat <- function(xx, ROW = FALSE){
  
  #if(ROW) make row vector
  
  if(is.matrix(xx))return(xx)
  
  cc <- names(xx)
  xx <- matrix(xx)
  rownames(xx) <- cc
  if(!ROW)xx <- t(xx)
  xx
}

.cleanRows <- function(xmat, xcol, STOP = FALSE, verbose = FALSE){
  
  ww <- which(duplicated(xmat[,xcol]))
  if(length(ww) == 0)return(xmat)
  
  if(STOP)stop(paste('duplicates in', xcol))
  
  tvec <- xmat[ww,xcol]
  if(length(ww) > 1)tvec <- paste0(tvec,collapse=', ')
  
  if(verbose)cat(paste('\nNote: removed duplicates in', xcol, ':\n', tvec) )
  
  xmat[-ww,]
}

HMC <- function (ff, fMin, fMax, ep, L, tree, sdat, ug, 
                 mu, sg, zz, R, SAMPR, distance, 
                 obsTrapRows, obsYr, seedNames, USPEC){
  
  #Hamiltonian Markov chain update
  
  getU <- function(q, U  = TRUE){   # yq = log(fec)
    
    # for Hamiltonian
    nseed <- ncol(R)
    fq <- exp(q)
    
    ww <- which(fq > fMax)
    vv <- which(fq > fMax)
    
    fq[ww] <- fMax[ww]
    fq[vv] <- fMin[vv]
    
    
    if( SAMPR | nseed > 1 ){
      fq <- matrix(fq,length(fq),ncol=ncol(R))*R[drop = FALSE,tree$specPlot,]
    }else{
      fq <- matrix(fq,ncol=1)
    }
    
    uvec <- ug[1]
    
    if(USPEC){
      uvec <- matrix( ug[attr(distance, 'species')], nrow(distance), ncol(distance) )
    }
    
    dmat <- uvec/pi/(uvec + distance^2)^2
    dmat[dmat < 1e-8] <- 0
    dmat[is.na(dmat)] <- 0
    
    
    plotyrs <- unique(sdat$plotyr)
    
    lambda <- kernYrRcpp(dmat, fq*zz, seedrow = sdat$drow,
                         treecol = tree$dcol, plotyrs, 
                         treeplotYr = tree[,'plotyr'], seedplotYr = sdat[,'plotyr'])
    ss <- as.matrix(sdat[,seedNames])
    lambda[lambda < 1e-6] <- 1e-6
    
    if(U){
      mmat  <- matrix(0,max(sdat$plotyr),1)
      sprob <- -ss*log(lambda) + activeArea*lambda
      ii    <- rep(sdat$plotyr, nseed)
      tmp   <- .myBy(as.vector(sprob), ii, ii*0+1, summat=mmat, fun='sum')
      tprob <- 1/sg*(q - mu)^2 + tmp[tree$plotyr]
      
      return( tprob ) 
    }
    
    kmat <- dmat[sdat$drow,tree$dcol]
    smat <- -ss/lambda + activeArea
    
    if(nseed == 1){
      svec <- ff*colSums( kmat*as.vector(smat) )
    }else{
      svec <- rep(0,length(q))
      for(m in 1:nseed){
        sv <-  colSums(kmat*as.vector(smat[,m]) )*fq[,m]
        svec <- svec + sv
      }
    }
    svec + (q - mu)/sg
  }
  
  q <- log( ff )
  p <- currentP <- rnorm(length(q)) 
  
  activeArea <- sdat$area
  
  # half step for momentum 
  
  p <- p - ep*getU(q, U = FALSE)/2
  
  
  # Alternate full steps for position and momentum
  wall <- 1:length(q)
  
  for (i in 1:L){
    
    q[wall] <- q[wall] + ep[wall]*p[wall]  
    
    if(i < L)p[wall] <- p[wall] - ep[wall]*getU(q, U = FALSE)[wall] 
    wall <- which(q < log(fMax))
  }
  
  # half step for momentum end
  p <- p - ep*getU(q, U = FALSE)/2
  
  # Negate momentum at end of trajectory to make proposal symmetric
  p <- -p
  
  # Evaluate potential and kinetic energies at start and end of trajectory
  currentU  <- getU(log( ff ), U  = TRUE) 
  proposedU <- getU(q, U   = TRUE)
  currentK  <- currentP^2/2 
  proposedK <- p^2/2
  
  # Accept or reject the state at end of trajectory, returning either 
  # the position at the end of the trajectory or the initial position
  sp <- currentU - proposedU + currentK - proposedK
  ii <- tree$plotyr
  
  pnow <- .myBy(sp*zz, ii, ii*0 + 1, fun='sum')
  
  a <- exp(pnow)
  wa <- which( runif(length(a)) < a)
  
  if(length(wa) > 0){
    wt <- which(tree$plotyr %in% wa & zz == 1)
    ff[wt] <- exp( q[wt] )
    ep[wt] <- ep[wt]*1.1
    ep[-wt] <- ep[-wt]*.9
  }else{
    ep <- ep*.9
  }
  
  list(fg = ff, epsilon = ep, accept = length(wa))
}

msarLagTemplate <- function(plag, data, icol, jcol, gcol, ocol, yeGr, verbose = FALSE){
  
  # index for sampling betaYr
  # data - data.frame created by msarSetup
  # icol - individual column (integer or name)
  # jcol - time column (integer)
  # gcol - group column (integer or name)
  # ocol - indicator for observed (1) or interpolated (0)
  
  ifac  <- gfac <- FALSE
  idata <- data[,icol]
  jdata <- data[,jcol]
  gdata <- data[,gcol]
  odata <- data[,ocol]
  
  jindex   <- sort( unique(jdata) )
  groups   <- yeGr
  times    <- sort( unique(jdata) )
  ngroup   <- length(groups)
  nyr      <- length(times)
  lagGroup <- vector('list',ngroup)
  names(lagGroup) <- groups
  
  lagMatrix <-  numeric(0)
  lagGroup  <- numeric(0)
  
  nall <- 0
  
  for(m in 1:ngroup){
    
    wm   <- which(gdata == m & odata == 1) 
    lagm <- numeric(0)
    
    im <- idata[wm]
    tall <- unique(im)
    
    jm <- jdata[wm]
    
    orm <- order(im,jdata[wm])
    
    imk <- match(im,tall)
    jmk <- match(jm,jindex)
    tmat <- matrix(NA,length(tall),nyr)
    tmat[ cbind(imk,jmk) ] <- wm
    rownames(tmat) <- tall
    colnames(tmat) <- jindex
    
    for(j in (plag+1):nyr){
      
      jname <- paste(tall,times[j],sep='-')
      tj    <- tmat[,j:(j-plag),drop = FALSE]
      rownames(tj) <- jname
      lagm <- rbind(lagm,tj)
    }
    wm <- unique( which(is.na(lagm),arr.ind  = TRUE)[,1] )
    if(length(wm) > 0)lagm <- lagm[drop = FALSE, -wm,]
    
    if(nrow(lagm) == 0)next
  
    nall <- nall + nrow(lagm)
    
    colnames(lagm) <- paste('lag',c(0:plag),sep='-')
    
    ttt <- columnSplit(rownames(lagm),'-')
    ord <- order(ttt[,1],ttt[,2],ttt[,3])
    
    lagm <- lagm[drop = FALSE, ord,]
    lagMatrix <- rbind(lagMatrix,lagm)
    
    lagGroup  <- c(lagGroup,rep(m,nrow(lagm)))
  }
  
  if(verbose){
    cat( paste('\nNumber of full observations with AR(', plag, ') model\ is: ', sep='') )
    print( nall ) 
    if(nall < 10)stop('not enough observations for AR(p), try reducing p')
  }
  
  list(matrix = lagMatrix, group = lagGroup)
}

getVarType <- function(vnames, data, i, j){
  
  # 'i'  - individual variable
  # 'j'  - time variable
  # 'ij' - individual/time
  
  if(is.factor(i))i <- as.character(i)
  
  id <- sort(unique(i))
  ni <- length(id)
  yr <- sort(unique(j))
  ny <- length(yr)
  i <- match(i, id)
  j <- match(j, yr)
  
  ij <- cbind(i, j)
  
  
  vnew <- vector('list',length(vnames))
  names(vnew) <- vnames
  
  for(k in 1:length(vnames)){
    
    cj <- data[,vnames[k]]
    if(is.factor(cj))cj <- as.character(cj)
    mj <- matrix(NA, ni, ny)
    mj[ ij ] <- cj
    
    rr <- suppressWarnings( apply(mj,1,range,na.rm  = TRUE) )
    rc <- suppressWarnings( apply(mj,2,range,na.rm  = TRUE) )
    rowSame <- all(rr[1,] == rr[2,])
    colSame <- all(rc[1,] == rc[2,])
    if(is.na(rowSame))rowSame <- FALSE
    if(is.na(colSame))colSame <- FALSE
    if(!rowSame & !colSame)vnew[[k]] <- 'ij'
    if( rowSame &  colSame)vnew[[k]] <- 'i'
    if(!rowSame &  colSame)vnew[[k]] <- 'j'
    if( rowSame & !colSame)vnew[[k]] <- 'i'
  }
  vnew
}

msarSetup <- function(data, plag, icol, jcol, gcol = NULL, yeGr,
                      minGroup=10, verbose = FALSE){
  
  # icol - column names in data for individual index; integer or factor
  # jcol - column name in data for time index; integer
  # gcol - column name for group index; integer or factor
  # pcol - column name for group names (character)
  # gmat - matrix of individual columns to retain in output
  # plag - AR lag 
  # minGroup - minimum group size
  
  if(!is.data.frame(data))stop('data must be a data.frame')
  
  huge <- 1e+10
  
 # data$obs <- 1
  
  if(is.factor(data[,icol]))data[,icol] <- droplevels(data[,icol])
  if(is.factor(data[,jcol]))data[,jcol] <- droplevels(data[,jcol])
  if(is.null(gcol)){
    gcol <- 'group'
    data$group <- rep(1,nrow(data))
  }else{
    if(is.factor(data[,gcol]))data[,gcol] <- droplevels(data[,gcol])
  }
  
  i <- as.character(data[,icol])
  j <- data[,jcol]
  g <- data[,gcol]
  if(is.factor(g))g <- as.character(g)
  
  iall <- unique(as.character(i))                  # preserve order
  
  jobs <- range(data[data$obsTrap == 1 | data$obs == 1,jcol])     
  jall <- range(range(j))
  if(jall[1] < (jobs[1] - plag))jall[1] <- jall[1] - plag
  if(jall[2] > (jobs[2] + plag))jall[2] <- jall[2] + plag
  
  jall <- c(jall[1]:jall[2])
  
  if( !is.null(yeGr)){
    groups <- yeGr
  }else{
    groups <- unique(g)
  }
  
  ii <- match(i,iall)
  jj <- match(j,jall)  # original times
  ni <- length(iall)
  nj <- length(jall)
  
  io <- i[data$obsTrap == 1 | data$obs == 1]
  jo <- j[data$obsTrap == 1 | data$obs == 1]
  ijIndex <- cbind( match(io,iall),match(jo,jall))
  ijFull  <- cbind( ii,jj )
  
  
 # baseMat <- matrix(0, length(iall), length(jall))
 # baseMat[ cbind(ii, jj) ] <- 1
 # rownames(baseMat) <- iall
 # colnames(baseMat) <- jall
 # base <- baseMat
 # base[base == 0] <- NA
  
  #groups with sufficient times
  tmp <- table(g,ii)
  tmp[tmp < plag] <- 0
  tmp[tmp > 1] <- 1
  tsum <- rowSums(tmp)
  
  if(verbose){
    cat('\nNote: no. trees with > plag years, by group:\n')
    print(tsum)
  }
  
  wlow <- which(tsum < minGroup)
  if(length(wlow) > 0){
    pll <- paste0(names(tsum)[wlow], collapse = ', ')
    if(verbose)cat( paste('\nsmall group(s): ',pll,'\n',sep='') )
  }
  
  vtypes    <- getVarType(colnames(data), data, ii, jj) 
  vtypes$obs <- 'ij'
  if('cropCount' %in% colnames(data))vtypes$cropCount <- 'none'
  if('cropFraction' %in% colnames(data))vtypes$cropFraction <- 'none'
  if('cropFractionSd' %in% colnames(data))vtypes$cropFractionSd <- 'none'
  if('fecMin' %in% colnames(data))vtypes$fecMin <- 'none'
  if('fecMax' %in% colnames(data))vtypes$fecMax <- 'none'
  
  tpy <- as.character(columnPaste(data$treeID,data$year))
  
  tmp <- fillMissing(vtypes, data, icol, jcol, ijFull)
  data   <- tmp$data
  naVars <- tmp$naVars
  tpn    <- as.character(columnPaste(data$treeID,data$year))
  data$obs[ which(!tpn %in% tpy) ] <- 0
  rm(tpy,tpn)
  
  inew <- numeric(0)
  
  for(i in 1:length(iall)){           # speed up needed
    wi <- which(data[,icol] == iall[i])
    inew <- c(inew,wi)
  }
    
  data <- data[inew,]
  
  
  i <- match(data[,icol], iall)
  g <- match(as.character(data[,gcol]),as.character(groups) )
  
  wideGroup <- apply(table(i,g),1,which.max) 
  
  ngroup <- length(groups)
  betaYr <- round( matrix(rnorm(ngroup*plag,0,.1),ngroup,plag), 3)
  rownames(betaYr) <- groups
  colnames(betaYr) <- paste('lag',c(1:plag),sep='-')
  
  list(xdata = data, times = jall, yeGr = groups,
       groupByInd = wideGroup, betaYr = betaYr, plag = plag)
}

fillMissing <- function(variables, data, icol, jcol, 
                        jtimes, ijFull = NULL){
  
  # ijFull is location in i by j matrix that includes added obs
  
  if(!is.data.frame(data))data <- as.data.frame(data, stringsAsFactors = F)
  
  if(!'obs' %in% colnames(data))data$obs <- 1
  
  id <- data[ ,icol]
  it <- data[ ,jcol]

  
  ids <- unique(id)
  
  rt <- range( c(jtimes, it) )
  jtimes <- rt[1]:rt[2]
  
  ijIndex <- cbind( match(id, ids), match(it, jtimes) )
  ijIndex <- ijIndex[ order(ijIndex[,1], ijIndex[,2]), ]
  
  if(is.null(ijFull))ijFull <- ijIndex
  ijFull <- ijFull[ order(ijFull[,1], ijFull[,2]), ]
  
  ny <- length(jtimes)
  ni <- length(ids)

  newData <- vector('list',ncol(data))
  names(newData) <- names(data)
  ffact <- which( sapply(data, is.factor) )
  naVars <- character(0)
  
  for(k in 1:ncol(data)){           #expand to pre- and post-data
    
    jmat  <- matrix(NA, ni, ny)
    
    vtype <- variables[names(data)[k]]
    tinySlope <- .00001
    
    kvar <- data[ ,k]
    sigFig <- getSigFig(kvar[1])
    
    if(k %in% ffact)kvar <- as.character(kvar)
    
    
    jmat[ ijIndex ]  <- kvar
    w0 <- which( is.na(jmat),arr.ind   = TRUE)
    w1 <- which(!is.na(jmat),arr.ind   = TRUE )
    
    if(vtype == 'time'){
      jmat <- matrix( jtimes, ni, ny, byrow=TRUE)
    }
    
    if(vtype == 'i'){
      w1 <- w1[!duplicated(w1[,1]),]
      w2 <- w0
      w2[,2] <- w1[match(w2[,1],w1[,1]),2]
      jmat[ w0 ] <- jmat[ w2 ]
    }
    if(vtype == 'j'){
      if(colnames(data)[k] == jcol){
        jmat <- matrix(jtimes,ni,ny, byrow   = TRUE)
      }else{
        w1 <- w1[!duplicated(w1[,2]),]
        w2 <- w0
        w2[,1] <- w1[match(w2[,2],w1[,2]),1]
        jmat[ w0 ] <- jmat[ w2 ]
      }
    }
    if( vtype == 'ij' ){            # use trend
      if( is.numeric(kvar) & !all(is.na(kvar)) ){
        minVal <- suppressWarnings( apply(jmat, 1, min, na.rm=T) )
        maxVal <- suppressWarnings( apply(jmat, 1, max, na.rm=T) )
        
        INCREASING <- FALSE
        
        if( colnames(data)[k] == 'diam'){
          
          tinySlope <- .005
          minVal[ minVal < .25 ] <- .25
          
          zz <- which(minVal >= maxVal)
          if(length(zz) > 0){
            minVal[zz] <- minVal[zz] - .1*ny
            maxVal[zz] <- maxVal[zz] + .1*ny
          }
          minVal[ minVal < .1 ] <- .1
            
          INCREASING <- TRUE
        }
        if( colnames(data)[k] == 'repMu'){
          INCREASING <- TRUE
          minVal <- 0
          maxVal <- 1
        }
        jmat <- .interpRows(jmat, INCREASING=INCREASING, minVal=minVal, maxVal=maxVal,
                            defaultValue=NULL,tinySlope=tinySlope) 
      }else{
        naVars <- c(naVars, colnames(data)[k])
      }
    }
    
    ktmp <- jmat[ijFull]
    if( k %in% ffact ) ktmp <- as.factor(ktmp)

    newData[[k]] <- ktmp
  }
  
  xdata <- data.frame(newData, stringsAsFactors = F)
  rownames(xdata) <- NULL
  
  
  list(data = xdata, naVars = naVars)
}

.propZ <- function(znow, last0first1, matYr){
  
  # repr - known repr from tdata
  # random walk proposal
  
  new <- matYr + sample( c(-1:1), nrow(znow), replace   = TRUE)
  
  ww  <- which(new < last0first1[,'last0'])
  new[ww] <- last0first1[ww,'last0']
  
  ww  <- which(new > last0first1[,'first1'])
  new[ww] <- last0first1[ww,'first1']
  
  new[last0first1[,'all0'] == 1] <- ncol(znow) + 1
  new[last0first1[,'all1'] == 1] <- 1
  
  down <- which(new < matYr & new > 0)
  znow[ cbind(down,new[down]) ] <- 1   # advance 1 year
  
  up <- which(new > matYr & new < ncol(znow))
  znow[ cbind(up,matYr[up]) ] <- 0     # delay 1 year
  
  znow[last0first1[,'all0'] == 1, ] <- 0
  znow[last0first1[,'all1'] == 1, ] <- 1
  
  wna <- which(is.na(znow))
  if(length(wna) > 0){
    znow[wna] <- 0
    znow <- t( apply( znow, 1, cumsum ) )
    znow[znow > 1] <- 1
  }
  
  list(zmat = znow, matYr = new)
}

.boxCoeffsLabs <- function( boxPars, labels, colLabs = NULL, cex=1, 
                            xadj = 0){
  
  ncols <- length(labels)
  if(is.null(colLabs))colLabs <- rep('black',ncols)
  
  at <- boxPars$xtick
  if(is.matrix(at))at <- colMeans( at )
  
  yfig <- par('usr')[3:4]
  dy   <- diff(yfig)
  yloc <- boxPars$stats[1,]
  pos  <- 2
  
  if(length(at) > 1){
    yends <- boxPars$stats[c(1,nrow(boxPars$stats)),]
    dends <- rbind(yfig[1] - yends[1,], yfig[2] - yends[2,])
    sides <- apply( abs(dends), 2, which.max)
    wt  <- which(sides == 1)
    yloc <- yends[1,]
    pos  <- 2
  }
  
  text( at, yloc - dy/20,labels, 
        offset = -.1,
        col=colLabs, pos = pos, srt=90, cex = cex)
  
  
  #if(length(wt) > 0)text(at[wt], yends[1,wt] - dy/20,labels[wt], 
  #                       offset = -.1,
  #                       col=colLabs[wt], pos=2, srt=90, cex = cex)
 # wt  <- which(sides == 2)
 # if(length(wt) > 0)text(at[wt], yends[2,wt] + dy/20,labels[wt], 
 #                        offset = -.1,
 #                        col=colLabs[wt], pos=4, srt=90, cex = cex)
}

getBins <- function(xx, nbin=15, pow=1){
  
  xrange <- range(xx, na.rm   = TRUE)
  bins <- unique( quantile(as.vector(xx[xx > 0]),seq(0,1,length=nbin)^pow) )
  bins <- c(0,bins)
  dbin <- diff(bins)
  
  bins[-1] <- bins[-1] - diff(bins)/2
  bins <- c(bins,max(bins)+1)
  db   <- diff(bins)
  w0   <- which(db < 1)
  
  if(length(w0) > 0)bins <- bins[-(w0+1)]
  sort(unique(bins))
}

plotCoeffs <- function(beta, bchain, burnin, ng, specNames, specCol, cex=.8,
                       xlab=''){
  
  fnames  <- .coeffNames(rownames(beta))
  nspec   <- length(specNames)
  
  nint <- grep(':',colnames(bchain))
  intt <- c(1:ncol(bchain))
  intt <- c(1:ncol(bchain))
  if(length(nint) > 0)intt <- intt[-nint]
  
  bc <- bchain[burnin:ng,intt,drop = FALSE]
  colnames(bc) <- .replaceString(colnames(bc),'species','')
  
  ylim <- quantile(bc, c(.1, .99))
  ylim[1] <- ylim[1] - .7*diff(ylim)
  ylim[2] <- ylim[2] + .7*diff(ylim)
  
  boxPars <- .boxCoeffs( chain = bc, snames = specNames, xaxt = 'n',
                        xlab = xlab, ylab='', ylim=ylim,
                        cols = specCol[specNames], addSpec = '')
  
  .boxCoeffsLabs( boxPars, specNames, specCol[specNames], cex=cex, xadj=1/nspec/2 )
  .plotLabel('intercept','bottomleft', cex = cex)
  
  bnint <- numeric(0)
  
  if(length(nint) > 0){
    
    bc <- bchain[burnin:ng,nint,drop = FALSE]
    colnames(bc) <- .replaceString(colnames(bc),'species','')
    scol <- rep(character(0), ncol(bc))
    
    ylim <- quantile(bc, c(.1, .99))
    ylim[1] <- ylim[1] - .5*diff(ylim)
    ylim[2] <- ylim[2] + .5*diff(ylim)
    
    bcc <- .boxCoeffsMultiSpec(bc, specNames, xlab = '', xaxt='n', ylim=ylim,
                               ylab=' ', cols = specCol[specNames], addSpec = '')
    bnint <- append(bnint, list(bcc))
  }
  invisible(list(intercept = boxPars, slopes = bnint))
}


.mastPlot2File <- function( output, plotPars ){
  
  parfile <- 'plotPars'
  if(is.null(plotPars))plotPars <- list()
  
  RMD <- plotPars$RMD
  if(RMD == 'pdf'){
    RMD <- 'pandoc'
    plotPars$MAPS <- FALSE
  }
  
  if(!RMD %in% c('html','pandoc','pdf'))plotPars$RMD  <- 'html'
  
  rfile   <- 'tmp.r'
  file.create(rfile)
  fileConn <- file(rfile)
  rtext <- '.mastPlot(output, plotPars)'
  
  writeLines(rtext, con = fileConn)
  close(fileConn)
  
  knitr::spin(rfile)
  file.remove( rfile )
  
  
  ttab <- table(output$inputs$treeData$plot)
  t2   <- names(ttab)[which.max(ttab)]
  t1   <- output$inputs$specNames[1]
  
  fname <- paste(t1,t2,'report.Rmd',sep='_')
  
  tmp <- readLines('tmp.md')
  rmdOut <- file(fname, "w")
  
  tmp <- .replaceString(tmp,".mastPlot(output, plotPars)","")
  
  capline <- grep('## ',tmp)
  newCaps <- .replaceString(tmp[capline],'## ','')

  caption <- grep('plot of chunk', tmp)
  cwords  <- tmp[caption - 3]
  cwords  <- .replaceString(cwords, ' ','')
  
  fnum <- c(1:length(caption))
  
  for(m in 1:length(caption)){
    tmp[caption[m]]  <- .replaceString(tmp[caption[m]], 'plot of chunk unnamed-chunk-1',
                                       newCaps[m])
  }
  comment <- grep('##',tmp)
  tmp <- tmp[-comment]
  tick <- grep('`',tmp)
  tmp <- tmp[-tick]
  
  day <- paste('date: ', date(), '\n', sep='')
  
  doc <- "output: 'html_document'\n"
  if(RMD == 'pandoc')doc <- "output: pdf_document\n"
  
  ptab <- paste("\n   print( knitr::kable(tj, knitr.kable.NA = '', caption = caps[wj, 2], format = '",
                RMD,"') )\n", sep='')
  
  top <- paste("title: '",t1, " at ", t2,"'\n", sep='')
  
  #new lines to rmdOut file
  cat("---\n", file = rmdOut)
  cat(top, file = rmdOut)
  cat(day, file = rmdOut)
  cat(doc, file = rmdOut)
  cat("---\n", file = rmdOut)
  
  cat("\nFor more background on this summary see http://rpubs.com/jimclark/281413\n", 
      file = rmdOut)
  
  cat("\n```{r, results = 'asis', echo = F}\n", file = rmdOut)
  cat("\ntfiles <- list.files('tables', full.names   = TRUE)",file = rmdOut)
  cat("\ntfiles <- tfiles[!tfiles %in% c('tables/words.txt', 'tables/captions.txt')]",
      file = rmdOut)
  cat("\ncaps  <- read.table('tables/captions.txt', as.is   = TRUE, header = FALSE)", 
      file = rmdOut)
  
  cat("\nwords <- read.table('tables/words.txt', header = FALSE)", file = rmdOut)
  cat("\ncolnames(words) <- NULL", file = rmdOut)
  cat("\ncat( paste0(unlist(words), collapse=' ') )\n", file=rmdOut)
  cat("\nfor(j in 1:length(tfiles)){\n",file = rmdOut)
  cat("\n   wj <- which(caps[,1] == tfiles[j])", file = rmdOut)
  cat("\n   tj <- read.table(tfiles[j], header   = TRUE)",file = rmdOut)
  cat(ptab,file = rmdOut)
  cat("\n}\n",file = rmdOut)
  cat("```\n", file = rmdOut)
  
  #all lines from rmdOut file
  for(i in 1:length(tmp)){
    cat(tmp[i], "\n", file = rmdOut, sep = "\t")
  }
  close(rmdOut)
}

mastPlot <- function(output, plotPars = NULL){
  
  # RMD - generate Rmarkdown
  
  verbose <- FALSE
  
  if(!is.null(plotPars)){
    if('RMD' %in% names(plotPars)){
      .mastPlot2File( output, plotPars )
      return()
    }
  }
  .mastPlot( output, plotPars, verbose )
}
  
.mastPlot <- function(output, plotPars, verbose){
    
  SAVEPLOTS <- SEEDCENSOR <- FALSE
  YR <- AR <- RANDOM <- TV <- SAMPR <- PREDICT <- SPACETIME <- FALSE
  SEEDDATA <- TRUE
  
  RMD <- NULL
  CONSOLE <- POINTS <- RMAT <- TRUE
  CONES <- FALSE
  caption <- character(0)
  
  data <- treeData <- priorUgibbs <- fecPred <- plotDims <- seedTraits <-
    plotHaByYr <- plotArea <- keepIter <- aUgibbs <- alphaMu <- NULL
  
  seedNames <- specNames <- R <- formulaFec <- formulaRep <- xfec <- xrep <-
    tdata <- seedData <- xytree <- xytrap <- distall <- ng <- burnin <- nplot <-
    ntree <- ntrap <- nyr <- maxFec <- bfec <- brep <- upar <- rgibbs <- 
    betaFec <- betaRep <- rMu <- rSe <- usigma <- fecMu <- fecSe <- matrMu <- 
    seedPred <- inputs <- chains <- parameters <- predictions <- 
    upars <- dpars <- trueValues <- betaYrMu <- betaYrSe <-  
    sgibbs <- ugibbs <- omegaE <- predPlots <- betaYrRand <- priorTable <-
    betaYrRandSE <- prediction <- eigenMu <- facLevels <- specPlots <- NULL
  randGroups <- formulaRan <- rnGroups <- reIndex <- xrandCols <- NULL  
  specGroups <- plotGroups <- yrIndex <- randomEffect <- yearEffect <- NULL
  pacfMat <- pacfSe <- pacsMat <- pacsSe <- obsRows <- NULL
  modelYears <- seedPredGrid <- treePredGrid <- NULL
  notFit <- character(0)
  censMin <- censMax <- NULL
  
  ugibbs <- matrix(0)
  
  outFolder <- 'mastPlots'
  yeGr <- NULL
  plotsPerPage <- 4
  MAPS <- TRUE
  
  for(k in 1:length(output))assign( names(output)[k], output[[k]] )
  for(k in 1:length(inputs))assign( names(inputs)[k], inputs[[k]] )
  for(k in 1:length(chains))assign( names(chains)[k], chains[[k]] )
  for(k in 1:length(parameters))assign( names(parameters)[k], parameters[[k]] )
  for(k in 1:length(prediction))assign( names(prediction)[k], prediction[[k]] )
  if('arList' %in% names(data)){
    for(k in 1:length(data$arList))
      assign( names(data$arList)[k], data$arList[[k]] )
    AR <- TRUE
    plag <- ncol(data$arList$betaYr)
  }
  
  if(!'seedData' %in% names(inputs)){
    SEEDDATA <- FALSE
  }else{
    if(length(seedData) == 0)SEEDDATA <- FALSE
  }
  
  if(!is.null(plotPars)){
    for(k in 1:length(plotPars))assign( names(plotPars)[k], plotPars[[k]] )
    if( 'trueValues' %in% names(plotPars)  ){
      TV <- TRUE
    }
  }
  if('trueValues' %in% names(inputs)){
    TV <- TRUE
    for(k in 1:length(inputs$trueValues))
      assign( names(inputs$trueValues)[k], inputs$trueValues[[k]] )
  }
  if('cropCount' %in% colnames(treeData))CONES <- TRUE
  notFit    <- output$data$setupData$notFit
  maxFec    <- output$inputs$maxFec
  specNames <- output$inputs$specNames
  seedNames <- output$inputs$seedNames
  if('rgibbs' %in% names(output$chains))SAMPR <- TRUE
  if(!RMAT)SAMPR <- FALSE

  
  if('randomEffect' %in% names(output$inputs)){
    RANDOM <- TRUE
    for(k in 1:length(randomEffect))
      assign( names(randomEffect)[k], randomEffect[[k]] )
    agibbs <- .orderChain(agibbs, specNames)
  }
  if('yearEffect' %in% names(output$inputs)){
    YR <- TRUE
    yrIndex <- output$data$setupYear$yrIndex
    if(is.null(yrIndex))yrIndex <- output$inputs$yrIndex
    if('p' %in% names(output$inputs$yearEffect)){
      if(output$inputs$yearEffect$p > 0){
        AR <- TRUE
        YR <- FALSE
      }
    }
    yeGr <- as.character(output$data$setupYear$yeGr)
    if(is.null(yeGr))yeGr <- output$data$setupData$yeGr
 #   yeGr <- .replaceString(yeGr, ' ', '')
    bygibbsR <- .orderChain(bygibbsR, yeGr)
    if(ncol(ugibbs) > 1)ugibbs <- .orderChain(ugibbs, specNames)
  }
  if(!YR)yeGr <- as.character(output$data$setupData$yeGr)
  ngroup <- length(yeGr)
  
  
  ngLab <- ng
  burnLab <- burnin
  if( 'keepIter' %in% names(inputs) ){
    burnin <- ceiling( burnin*keepIter/ng )
    ng <- keepIter
  }
  nspec <- length(specNames)
  
  if(!is.null(RMD))SAVEPLOTS <- CONSOLE <- FALSE
  if(SAVEPLOTS)CONSOLE <- FALSE
  
  if(SAVEPLOTS & verbose){
    tt <- paste('\nPlots saved to ',outFolder,'/\n',sep='')
    cat(tt)
  }
  
  if(AR)YR <- FALSE
  if(!is.null(censMin) & !is.null(censMax))SEEDCENSOR <- TRUE
  
  xmean <- output$data$setupData$xmean  # to unstandardize xfec, xrep
  xsd   <- output$data$setupData$xsd
  
  tdata <- output$inputs$treeData
  sdata <- output$inputs$seedData
  
  xfecu2s <- output$data$setupData$xfecu2s
  xrepu2s <- output$data$setupData$xrepu2s
  
  xfecs2u <- output$data$setupData$xfecs2u
  xreps2u <- output$data$setupData$xreps2u
  
  fnames  <- .coeffNames(rownames(betaFec))
  rnames  <- .coeffNames(rownames(betaRep))
  
  xfec <- output$data$setupData$xfec
  xrep <- output$data$setupData$xrep
  Qf   <- ncol(xfec)/nspec
  Qr   <- ncol(xrep)/nspec
  
  if( !is.matrix(pacfMat) )pacfMat <- t( as.matrix( pacfMat ) )
 # if( !is.matrix(pacsMat) )pacsMat <- t( as.matrix( pacsMat ) )
  if( !is.matrix(pacfSe) )pacfSe <- t( as.matrix( pacfSe ) )
  
  
  
  ###############
  rm(treeData)
  ##############
  
  outTables <- 'tables/table_'
  if(file.exists('tables')){
    fi <- list.files('tables', full.names   = TRUE)
    for(j in 1:length(fi)){
      file.remove(fi[j])
    }
  }else{
    dir.create('tables')
  }
  
  captions <- character(0)
  
  sumOut <- summary.mastif( output, verbose = FALSE )
  
  if(!is.null(RMD)){
    ww <- which(!names(sumOut) == 'words')
    for(k in ww){
      xk <- sumOut[[k]]
      ok <- paste( outTables, names(sumOut)[k], '.txt', sep='')
      
      cap <- attr(xk, 'caption')
      if(is.null(cap))cap <- 'tcap'
      write.table(xk, quote = FALSE, file = ok, row.names   = TRUE) 
      
      captions <- rbind(captions, cbind(ok, cap) )
    }
  }
  words <- paste0( sumOut$words, collapse = '. ')
  words <- .replaceString(words, '\n', '')
  
  if(!is.null(RMD)){
    write.table(captions, file='tables/captions.txt', row.names = FALSE, col.names = FALSE)
    write.table(sumOut$words, file='tables/words.txt', row.names = FALSE, col.names = FALSE)
  }
  
  if(is.null(yeGr))yeGr <- 'all'
  
  nspec  <- length(specNames)
  years  <- sort(unique(tdata$year))
  nyr    <- length(years)
  ngroup <- length(yeGr)
  plots  <- sort(unique(as.character(tdata$plot)))
  nplot  <- length(plots)
  
  tmp1 <- as.character(tdata$plot)
  tmp2 <- tdata$year
  
  if(SEEDDATA){
    ntype  <- length(seedNames)
    if(!is.null(seedPredGrid)){
      PREDICT <- TRUE
      tmp1 <- c(tmp1, as.character(seedPredGrid$plot))
      tmp2 <- c(tmp2, seedPredGrid$year)
    }
  }
  
  plotYrTable <- table(plot = tmp1, year = tmp2)
  rm(tmp1)
  rm(tmp2)
  
  cfun <- colorRampPalette( c('#66c2a5','#fc8d62','#8da0cb') )
  specCol <- cfun(nspec)
  names(specCol) <- specNames
  cols <- specCol
  
  gfun <- colorRampPalette( c("#8DD3C7", "#BEBADA", "#FB8072",
                              "#80B1D3", "#FDB462") )
  groupCol <- gfun(ngroup)
  names(groupCol) <- yeGr
  
  gfun <- colorRampPalette( c("forestgreen","#8DD3C7", "#BEBADA", "#FB8072",
                              "#80B1D3", "#FDB462", "brown") )
  plotCol <- gfun(nplot)
  names(plotCol) <- plots
  
  if(SAVEPLOTS){
    ff <- file.exists(outFolder)
    if(!ff)dir.create(outFolder)
  }
  
  
  ########### MCMC chains
  
  if(is.null(RMD)) graphics.off()
  
  refVals <- NULL

  words <- .chainPlot(chains$brep, burnin, 'maturation', ngLab, burnLab,
             refVals = refVals, CONSOLE, RMD, SAVEPLOTS, outFolder)
  
  bfecFit   <- chains$bfec[,!colnames(chains$bfec) %in% notFit, drop = FALSE]
  words <- .chainPlot(bfecFit, burnin, 'fecundity', ngLab, burnLab, 
             refVals = refVals, CONSOLE, RMD, SAVEPLOTS, outFolder)
  
  # species with estimated seed dispersal
  trapSpec <- tapply(tdata$fit, tdata$species, sum)
  
  trapSpec <- names(trapSpec)[ trapSpec > 10 ]
  trapSpec <- trapSpec[ trapSpec %in% colnames(ugibbs) ]
  if(length(trapSpec) == 0)trapSpec <- colnames(ugibbs)
  
  if(SEEDDATA){
    
    USPEC <- TRUE
    if(ncol(ugibbs) == 1)USPEC <- FALSE
    
    ALLONE <- FALSE
    if(USPEC)ALLONE <- TRUE
    
  #  ugibbs  <- chains$ugibbs[burnin:ng,trapSpec,drop = FALSE]
    
    sord <- order( colMeans(ugibbs), decreasing   = TRUE) 
    sord <- colnames(ugibbs)[sord]
    ugibbs <- ugibbs[,sord,drop=FALSE]
    
    if(TV){
      refVals <- inputs$trueValues$upar[sord]
      if(is.na(refVals))refVals <- inputs$trueValues$upar
    }
    
    intval <- priorTable[drop = FALSE, sord, c('priorU','priorVU','minU','maxU')]
    colnames(intval) <- c('mean','var','min','max')
    ylim <- range(ugibbs)
    
    ql <- min( qnorm(.05,intval[,'mean'],sqrt(intval[,'var'])) )
    qh <- max( qnorm(.95,intval[,'mean'],sqrt(intval[,'var'])) )
    
    if(ylim[1] > ql)ylim[1] <- ql
    if(ylim[2] < qh)ylim[2] <- qh
    
    words <- .chainPlot(ugibbs, burnin,
                        label = 'dispersal parameter u (with prior)',
                        ngLab = ngLab, burnLab = burnLab, refVals = refVals, CONSOLE, 
                        RMD, SAVEPLOTS, outFolder, ALLONE  = TRUE, cols = specCol, 
                        ylim = ylim, intval=intval)
    if(USPEC){
      words <- .chainPlot(chains$priorUgibbs, burnin, 'dispersal mean and variance', 
                          ngLab, burnLab,
                          refVals = NULL, CONSOLE, RMD, SAVEPLOTS, outFolder)
    }
    
    words <- .chainPlot(chains$sgibbs, 
                        burnin, 'variance sigma', ngLab, burnLab,
                        refVals = NULL, CONSOLE, RMD, 
                        SAVEPLOTS, outFolder)
  }else{
    words <- .chainPlot(chains$sgibbs[,1,drop = FALSE], 
                        burnin, 'variance sigma', ngLab, burnLab,
                        refVals = NULL, CONSOLE, RMD, 
                        SAVEPLOTS, outFolder)
  }
  if(ncol(ugibbs) > 1){
    sord <- order( colMeans(ugibbs), decreasing   = FALSE) 
    sord <- colnames(ugibbs)[sord]
    ugibbs <- ugibbs[,sord,drop=FALSE]
  }

############ R matrix
  
  if(SAMPR){
    
    mg   <- chains$rgibbs
    posR <- attr(rMu,'posR')
    rff  <- NULL
    
    tmp <- columnSplit(colnames(mg),'_')
    seedCols <- tmp[,1]
    
    wdash <- grep('-',tmp[1,2])
    if(length(wdash) > 0){
      tmp  <- columnSplit(tmp[,2],'-')
      specCols <- tmp[,1]
      plotCols <- tmp[,2]
      np <- length(plots)
    }else{
      specCols <- tmp[,2]
      plotCols <- NULL
      np <- 1
    }
    
    for(j in 1:np){
      if(!is.null(plotCols)){
        wj <- which(plotCols == plots[j])
        if(length(wj) == 0)next
        rj <- range(mg[,wj])
        if(rj[1] == 1 | rj[2] == 0)next
      }else{
        wj <- 1:ncol(mg)
      }
      
      label <- paste('M matrix', plots[j])
      
      if(TV){
        tt <- columnSplit(colnames(mg)[wj])
        rff <- inputs$trueValues$R[ cbind(tt[,2],tt[,1]) ]
      }
      
      words <- .chainPlot(mg[,wj,drop = FALSE], burnin, label, 
                          ngLab, burnLab,ylim=c(0,1),
                 refVals = rff, CONSOLE, RMD, SAVEPLOTS, outFolder)
    }
    
    if(is.null(RMD)) graphics.off()
    
    if(SAVEPLOTS)pdf( file=.outFile(outFolder,'Mpars.pdf') )
    
    espec <- sort(unique(specCols))
    
    par(mfrow=c(length(espec),1), bty='n', mar=c(.5,4,.5,3), oma=c(4,3,2,4))
    tlab <- ''
    
    plotColors <- plotCol[plotCols]
    plpt <- character(0)
    
    for(j in 1:nspec){
      
      cj <- which(specCols  == specNames[j])
      if(length(cj) == 0)next
      
      unkn <- grep('UNKN',colnames(mg)[cj])
      if(length(unkn) > 0){
        cj <- c(cj[-unkn],cj[unkn])
      }
      
      cols <- plotColors[cj]
      plpt <- c(plpt, unique(names(cols)))
      
      boxPars <- .boxCoeffs(mg[burnin:ng,cj,drop = FALSE], specNames[j], xlab = tlab,
                            ylab=specNames[j], addSpec='', ylim=c(0,1),
                            cols = cols, yaxt='n')
      wt <- c(match(seedCols[cj], seedNames), 1000)
      wt[length(wt)] <- wt[length(wt)] + 1
      wt <- which(diff(wt) != 0)
      abline(v=wt+.5,col='grey')
      text(wt-.1,.5,seedCols[cj[wt]], srt=90, cex=.8)
      
      abline(h=1,col='grey', lwd=1)
      axis(2, at=c(0,1), las=2)
      tlab <- ''
    }
    
    mtext('Fraction from these species...', side=2, line=.4, outer  = TRUE, cex=1.2)
    mtext('...counted as these seed types', side=3, line=.4, outer  = TRUE, cex=1.2)
    
    ncol <- round(length(plpt)/4) + 1
    
    plpt <- unique(plpt)
    
    cornerLegend('bottomright', plpt, text.col = plotCol[plpt],
                 cex=.9, bty='n', ncol=ncol)
    
    if(CONSOLE)
      readline('species -> seed type -- return to continue ')
    if(SAVEPLOTS)dev.off( )
    
    if(is.null(RMD)){
      graphics.off()
    }else{
      words <- paste( 'Contribution of each species to seedNames' )
      message( words )
      caption <- c(caption, words)
    }
    
  ########################
  
    #inverse probability species h| seedtype m
    
    fec <- output$prediction$fecPred
    
    nsim <- 40
    ksim <- sample(c(burnin:ng), nsim, replace  = TRUE)
    
    npairs <- columnSplit(colnames(rgibbs),'_')[,c(2,1)]
    
    rmat <- sprob <- sprob2 <- parameters$rMu*0
    
    for(m in 1:length(plots)){
      
      wj <- grep( paste('-',plots[m],sep=''), colnames(rgibbs))
      
      if(length(wj) <= 1)next
      mm <- rgibbs[,wj,drop = FALSE]
      mp <- npairs[wj,]
      
      wm <- which(fec$plot == plots[m])
      ff <- fec$fecEstMu[wm]
      ss <- fec$fecEstSe[wm] + .00001
      mt <- fec$matrEst[wm]
      
      rrow <- grep(plots[m],rownames(rmat))
      
      for(k in 1:nsim){
        
        rmat <- rmat*0
        rmat[npairs[wj,]] <- mm[ksim[k],]
        rmm  <- rmat[drop = FALSE,rrow,]
        
        fk <- .tnorm(length(ff), 0, maxFec, ff, ss)*mt
        tf <- tapply(ff, list(species = fec$species[wm]), FUN=sum)
        tf <- tf/sum(tf)
        names(tf) <- paste(names(tf),plots[m],sep='-')
        tmat <- rmm*0
        tmat[names(tf),] <- rep(tf, ncol(tmat))
        
        sf <- rmm*tmat/matrix( colSums(rmm*tmat),nrow(rmm), ncol(rmm), byrow  = TRUE )
        sf[is.na(sf)] <- 0
        sprob[rrow,] <- sprob[rrow,] + sf
        sprob2[rrow,] <- sprob2[rrow,] + sf^2
      }
    }
    seed2SpecMu <- sprob/nsim
    sse <- sprob2/nsim - seed2SpecMu^2
    sse[sse < 1e-30] <- 0
    seed2SpecSe <- sqrt( sse )
    
    if(max(seed2SpecSe, na.rm  = TRUE) > 1e-20){
      
      if(SAVEPLOTS)pdf( file=.outFile(outFolder,'undiffSeed.pdf') )
      
      kplots <- character(0)
      
      for(m in 1:length(plots)){
        rrow <- grep(plots[m],rownames(seed2SpecMu))
        smm  <- seed2SpecMu[rrow,]
        wk   <-  which(smm > 1e-20 & smm < 1, arr.ind  = TRUE)
        if(length(wk) == 0)next
        kplots <- c(kplots,plots[m])
      }
      
      tt <- .getPlotLayout(length(kplots))
      par(mfrow=tt$mfrow, bty='n', mar=c(3,3,1,1), oma=c(3,3,3,4))
      aspec <- character(0)
      
      ucol <- grep('UNKN', colnames(seed2SpecMu))
      
      ISPLOT <- FALSE
      
      for(m in 1:length(kplots)){
        
        rrow <- grep(kplots[m],rownames(seed2SpecMu))
        smm  <- seed2SpecMu[rrow,]
        emm  <- seed2SpecSe[rrow,]
        wk   <-  which(smm > 1e-20 & smm < 1, arr.ind  = TRUE)
        if(length(wk) == 0)next
        
        wk <- wk[wk[,2] == ucol,]
        if(length(wk) == 0)next
        
        ISPLOT <- TRUE
        
        cspec <- columnSplit(rownames(wk),'-')[,1]
        aspec <- c(aspec,cspec)
        colm <- specCol[cspec]
        tmp  <- barplot(smm[wk], beside  = TRUE, col=.getColor(colm,.5), 
                        border=colm, ylim=c(0,1),yaxt='n', lwd=2)
        labels <- FALSE
        if(m %in% tt$left)labels=c(0,1)
        axis(2, c(0,1), labels=labels)
        segments(tmp, smm[wk], tmp, smm[wk] + 1.96*emm[wk], lwd = 1.5, col=colm)
        segments(tmp-.1, smm[wk] + 1.96*emm[wk], tmp+.1, smm[wk] + 1.96*emm[wk],
                 lwd=1.5, col=colm)
        .plotLabel(kplots[m], 'topright', above  = TRUE)
      }
      
      if(ISPLOT){
        mtext('...from these species', side=1, line=.4, outer  = TRUE, cex=1.2)
        mtext('Fraction of unknown seed type...', side=2, line=.4, outer  = TRUE, cex=1.2)
        
        aspec <- sort(unique(aspec))
        
        cornerLegend('bottomright', aspec, text.col = specCol[aspec],
                     cex=1.1, bty='n', ncol=1)
        
        if(CONSOLE)
          readline('Species to undiff seed -- return to continue ')
        if(SAVEPLOTS)dev.off( )
        
        if(is.null(RMD)){
          graphics.off()
        }else{
          words <- paste( 'Contribution of each species to undifferentiated type' )
          message( words )
          caption <- c(caption, words)
        }
      }else{
        dev.off()
      }
    }
  }
  
  ##############################
  
  if(RANDOM){
    
    aMu <- parameters$aMu
    
    att <- aMu*0
    wc <- 1
    if(length(aMu) > 1){
      diag(att) <- 1
      wc <- which(att == 1)
    }
    
    vaa <- agibbs[,wc,drop = FALSE]
    vrr <- apply(vaa,2,sd)
    if( max(vrr) > 1e-5 ){
      
      words <- .chainPlot(aUgibbs[,wc,drop = FALSE], burnin, 
                          'random effects variance', 
                          ngLab, burnLab,
                 refVals = NULL, CONSOLE, RMD, SAVEPLOTS, outFolder)
      caption <- c(caption, words)
    }
    if(is.null(RMD)) graphics.off()
  }
  
  
  ########### coefficients
  
  if(SAVEPLOTS)pdf( file=.outFile(outFolder,'fecundityCoeffs.pdf') )

  fitCols <- colnames(xfecu2s)
  
  # standardized betas
  bfecStnd <- bfec*0
  bfecStnd[,fitCols] <- bfec[,fitCols]%*%t(xfecu2s)
  brepStnd <- brep%*%t(xrepu2s)
  colnames(brepStnd) <- colnames(brep)
  
  wss <- grep('species', colnames(bfecStnd))
  wff <- grep('species', rownames(betaFec))
  
  if(length(wss) > 0 & length(wff) == 0)
    rownames(betaFec) <- colnames(bfecStnd)
  
  if(nspec == 1){
              
    par(mfrow=c(1,2), mar=c(5,4,2,1), bty='n')
    tmp <- .boxplotQuant( brepStnd[drop = FALSE,burnin:ng,], add = F, xaxt = 'n',
                          xlim = NULL, ylab='Standard deviations', 
                          outline = FALSE, col=.getColor('black', .2), 
                          border='black', lty=1, boxfill=NULL )
    
    axis(1, at = c(1:nrow(betaRep)), labels=rnames, las = 2)
    abline(h=0, col='grey', lwd=2, cex.axis=.6)
    title('a) Maturation')
    if(TV)abline(h=inputs$trueValues$betaRep, lwd=2, lty=2, col='grey')
    
    tmp <- .boxplotQuant( bfecStnd[drop = FALSE,burnin:ng,!colnames(bfecStnd) %in% notFit], 
                          add = F, xaxt = 'n',
                          xlim = NULL, ylab='', 
                          outline = FALSE, col=.getColor('black', .2), 
                          border='black', lty=1, boxfill=NULL )
    axis(1, at = c(1:nrow(betaFec)), labels=fnames, las = 2)
    abline(h=0, col='grey', lty=2)
    title('b) Fecundity')
    if(TV){
      abline(h=inputs$trueValues$betaFec, lwd=2, lty=2, col='grey')
      text(2,inputs$trueValues$betaFec[1] + 2,'true values')
    }
    
  }else{
    
    par(mfrow=c(2,2), mar=c(1,3,1,.2), bty='n', oma=c(2,2,1,1))
    
    rr <- rownames(betaFec)
    if(length(notFit) > 0)rr <- rr[!rr %in% notFit]

    plotCoeffs(beta = betaRep, bchain = brepStnd, burnin, ng, 
               specNames, specCol, cex=.85,
               xlab = 'Maturation')
    
    plotCoeffs(betaFec[rr,], bfecStnd[,rr], burnin, ng, specNames, specCol, cex=.85,
               xlab = 'Fecundity')
    mtext('Standardardized estimates',side=2, outer  = TRUE)
  }
  
  if(CONSOLE)readline('fecundity, maturation -- return to continue ')
  if(SAVEPLOTS)dev.off( )
  
  if(is.null(RMD)){
    graphics.off()
  }else{
    words <- paste( 'Posterior estimates for maturation and fecundity parameters' )
    message( words )
    caption <- c(caption, words)
  }
  
  ############ fecundity credible interval by diameter, individual species
  
  diam90 <- NULL  # if SAVEPLOTS, this list holds 90% CI for diameter effect
  
  obsRows  <- which(tdata$obs == 1)
  
  if(SAVEPLOTS){
    
    nsim <- 500
    
    dseq <- seq(1, 500, by=1)
    dmat <- cbind(1, dseq, dseq^2)
    
    mrList <- vector('list',nspec)
    names(mrList) <- specNames
    dfList <- srList <- mrList
    
    for(j in 1:nspec){
      
      fspec  <- fecPred[obsRows,]
      wrow   <- which(fspec$species == specNames[j])
      fspec  <- fspec[wrow,]
      dspec  <- fspec$diam
      
      # plot up to last 5 trees
      qd    <- 5/length(dspec)
      dtops <- quantile(dspec, 1 - qd)
      
      xmat <- dmat[dmat[,2] < dtops,]
      ntt  <- nrow(xmat)
      
      #  ntt  <- length(ttt)                          # tree yr
      ksamp <- sample(burnin:ng,nsim,replace  = TRUE)   # MCMC row
      bcols <- colnames(brep)
      fcols <- colnames(bfec)
      if(nspec > 1){
        bcols <- bcols[ grep(specNames[j], bcols) ]
        fcols <- fcols[ grep(specNames[j], fcols) ]
        fcols <- fcols[ c(1, which( endsWith(fcols, ':diam') | endsWith(fcols, ':I(diam^2)'))) ]
      }
      fcols <- fcols[ c(1, which( endsWith(fcols, 'diam') | endsWith(fcols, 'I(diam^2)'))) ]
      
      mk <- sk  <- matrix(NA, ntt, nsim)
      df <- matrix(NA, 2, nsim)
      rownames(df) <- c('dstar','fstar')
      
      for(k in 1:nsim){  # predict maturation, fecundity
        
        bf <- t( bfec[ksamp[k],fcols, drop = FALSE] )
        br <- t( brep[ksamp[k],bcols, drop = FALSE] )
        
        dstar <- -bf[2]/2/bf[3]
        fstar <- exp(  bf[1] + bf[2]*dstar + bf[3]*dstar^2 )
        df[,k] <- c(dstar,fstar)
        
        mk[,k] <- pnorm( xmat[,1:2]%*%br )
        sk[,k] <- exp( xmat%*%bf )*mk[,k]
      }
      
      ww <- which(df[1,] < dtops & df[1,] > 0)
      df <- df[,ww]
      ll <- length( ww )/nsim
      
      mm1 <- ss1 <- dd1 <- numeric(0)
      
      if( length(ww)/nsim > .5 ){
        mm1 <- apply( mk[,ww], 1, quantile, c(.5, .05, .95) )  # 90% credible interval for maturation
        ss1 <- apply( sk[,ww], 1, quantile, c(.5, .05, .95) )  # for 90% credible interval on fecundity
        dd1 <- apply( df, 1, quantile, c(.5, .05, .95) )
        attr(dd1, 'fraction') <- ll
      }
      
      mrList[[j]] <- mm1  # 90% credible interval for maturation
      srList[[j]] <- ss1  # for 90% credible interval on fecundity
      dff <- dd1
      dfList[[j]] <- dff
      
      #   ymax <- max( c(ymax, srList[[j]][1,]) )
      #   xmax <- max( c(xmax, dtops) )
    }
    
    for(j in 1:nspec){
      
      mr <- mrList[[j]]
      sr <- srList[[j]]
      df <- dfList[[j]]
      
      if(length(sr) == 0)next
      
      ds <- 1:ncol(sr)
      fjname <- paste('diam90_',specNames[j], '.pdf', sep='')
      
      pdf( file=.outFile(outFolder,fjname) )
      
      par(bty='n')
      
      fspec  <- fecPred[obsRows,]
      wrow   <- which(fspec$species == specNames[j])
      fspec  <- fspec[wrow,]
      
      y1 <-  max(sr)
      y2 <- max(fspec$fecEstMu, na.rm=T)
      y3 <- df[1,2]*4
      
      ymax <- min( c(y1, y2, y3), na.rm=T )
      ymax <- max( c(ymax, 1.2*max( sr[1,])) )
      
      plot(ds, sr[1,], type='l' , ylim = c(0, ymax), lwd=2, col = 'brown',
           xlab = 'Diameter (cm)', ylab = 'Seeds (90% CI)')
      .shadeInterval( ds, t(sr[2:3,]) , col=.getColor('tan', .1) )
      .plotLabel(specNames[j])
      
      if(SAVEPLOTS)dev.off( )
    }
    diam90 <- list(maturation = mrList, fecundity = srList, 
                   inflection = dfList)
  }
   
  
  
  ############ maturation, fecundity
  
  if(SAVEPLOTS)pdf( file=.outFile(outFolder,'maturation.pdf') )
  
  mfrow <- c(nspec,2)
  par(mfrow=mfrow, bty='n', mar=c(1,3,.5,3),  oma=c(3,2,1.2,2))
  
  nsim <- 500
  lseq <- seq(0, 1, length=1000)
  lseq <- cumsum(lseq*(1 - lseq))
  lseq <- lseq/max(lseq)
  lseq <- lseq[!duplicated(lseq)]
  
  fss   <- rep(specNames[1], length(obsRows))
  
  if(nspec > 1){
    scols <- paste('species',specNames,sep='')
    xftmp <- xfec[,scols]
    xrtmp <- xrep[,scols]
    fss   <- specNames[ apply(xftmp[obsRows,], 1, which.max) ]
  }
  
  # simulate first to evaluate scale
  
  mrList <- vector('list',nspec)
  names(mrList) <- specNames
  drList <- srList <- arList <- mrList
  
  xmax <- ymax <- 0
  
  for(j in 1:nspec){
    
    fspec  <- fecPred[obsRows,]
    wrow   <- which(fss == specNames[j])
    fspec  <- fspec[wrow,]
    pspec  <- fspec$plot
    mspec  <- fspec$matrEst
    dspec  <- fspec$diam
    
    # plot up to last 5 trees
    qd <- 5/length(dspec)
    dtops <- quantile(dspec, 1 - qd)
    
    xrspec <- xrep[rownames(fspec),]
    xfspec <- xfec[rownames(fspec),]
    
    pspec <- as.character(pspec)
    pall  <- sort(unique(pspec))
    
    dq   <- round( quantile(dspec, lseq), 1)
    dq   <- unique(dq)
    
    ttt  <- nn2(dspec, dq, k = 1)[[1]]
    www  <- which(duplicated(ttt))
    if(length(www) > 0)ttt <- ttt[-www]
    
    ntt  <- length(ttt)                          # tree yr
    ksamp <- sample(burnin:ng,nsim,replace  = TRUE)   # MCMC row
    bcols <- colnames(brep)
    fcols <- colnames(bfec)
    if(nspec > 1){
      bcols <- bcols[ grep(specNames[j], bcols) ]
      fcols <- fcols[ grep(specNames[j], fcols) ]
    }
    
    ddcol <- grep('diam',bcols)
    fdcol <- grep('diam',fcols)
    
    xtt <- xrspec[ttt,bcols]            # standardized
    ftt <- xfspec[ttt,fcols]
    dtt <- dspec[ttt]
    
    fcol <- c(1, grep('diam',colnames(ftt)))  # isolate diameter effect
    rcol <- c(1, grep('diam',colnames(xtt)))
    
    fcols <- fcols[fcol]
    bcols <- bcols[rcol]
    
    grf <- grep(':diam:',fcols)
    if(length(grf) > 0)fcols <- fcols[-grf]
    grf <- grep('diam:',fcols)
    if(length(grf) > 0)fcols <- fcols[-grf]
    
    grb <- grep(':diam:',bcols)
    if(length(grf) > 0)bcols <- bcols[-grf]
    grf <- grep('diam:',bcols)
    if(length(grf) > 0)bcols <- bcols[-grf]
    
    xtmu <- matrix(colMeans(ftt), nrow(ftt), ncol(ftt), byrow  = TRUE)
    gcol <- grep('diam',colnames(ftt))
    
    rk <- mk <- ak <- sk   <- matrix(NA, ntt, nsim)
    
    for(k in 1:nsim){  # predict maturation, fecundity
      
      ss <- sqrt( sgibbs[ksamp[k],1] )
      bf <- bfecStnd[ksamp[k],fcols, drop = FALSE] 
      bint <- 0
      
      if(RANDOM){                             # check standardized
        bfr <- rmvnormRcpp(1, aMu[1,]*0, aMu)
        names(bfr) <- colnames(aMu)
        bfr <- bfr[colnames(aMu) %in% fcols]
        bf[,names(bfr)] <- bf[,names(bfr)] + bfr
      }
      
      mk[,k] <- pnorm( xtt[,bcols]%*%t( brepStnd[ksamp[k],bcols, drop = FALSE] ) )
      muk <- ftt[,fcols]%*%t( bf )
      #   sk[,k] <- exp( muk + rnorm(ntt, 0, ss) )*mk[,k]
      
      sk[,k] <- exp( muk )*mk[,k]
      ak[,k] <- exp( muk + rnorm(ntt, 0, ss) )*mk[,k]
    }
    
    wd <- which(dtt <= dtops)
    
    mrList[[j]] <- apply( mk[wd,], 1, quantile, c(.5, .05, .95) )        # 90% credible interval for maturation
    srList[[j]] <- sqrt(apply( sk[wd,], 1, quantile, c(.5, .05, .95) ))  # for 90% credible interval on fecundity
    arList[[j]] <- sqrt(apply( ak[wd,], 1, quantile, c(.5, .05, .95) ))  # for 90% predictive interval
    drList[[j]] <- dtt[wd]
    ymax <- max( c(ymax, srList[[j]][1,]) )
    xmax <- max( c(xmax, dtops) )
  }
  

  ym   <- max( c(1.2*ymax, 20) )  # sets vertical scale
  tt   <- sqrtSeq(ym)
  at   <- tt$at
  labs <- tt$labs
  xlim <- c(0, max(tdata$diam, na.rm  = TRUE))
  npat <- 5
  if(nspec > 4)npat <- 4
  if(nspec > 7)npat <- 3
  
  # sqrt scale
  for(j in 1:nspec){
    
    fspec  <- fecPred[obsRows,]
    wrow   <- which(fss == specNames[j])
    fspec  <- fspec[wrow,]
    pspec  <- fspec$plot
    mspec  <- fspec$matrEst
    dspec  <- fspec$diam
    
    mr <- mrList[[j]]
    sr <- srList[[j]]
    ar <- arList[[j]]
    dtt <- drList[[j]]
    
    xaxt <- 'n'
    if(j == nspec)xaxt <- 's'
    
    plot(NULL, xlim=xlim,ylim=c(0,1),xlab='',
         ylab='', xaxt = xaxt)
    if(j < nspec)axis(1, labels = FALSE)
    wp <- match(as.character(pspec),plots)
    
    fj <- mspec
    wj <- which(fj < .001 | fj > .999)
    fj[wj] <- jitter(fj[wj],.15)
    if(POINTS)points(dspec,fj, pch=16, col=.getColor(plotCol[pspec],.3),
                     cex=.3)
    
    .shadeInterval( dtt, t(mr[2:3,]) , col=
                      .getColor('tan', .1) )
    lines(dtt, mr[1,], col='white', lwd=5)
    lines(dtt, mr[1,], lwd=2)
    
    if(j == 1)title('Maturation')
    
    fpp <- fspec$fecEstMu*fspec$matrEst
    
    xaxt <- 'n'
    if(j == nspec)xaxt <- 's'
    
    # histogram of diameter polygon getpoly
    
    plot(NULL, xlim=xlim, ylim=range(at), xlab='', ylab='', 
         xaxt = xaxt, yaxt='n')
    if(j < nspec)axis(1, labels = FALSE)
    axis(2, at = at, labels = labs)
    
    
    tov  <- hist(fspec$diam, nclass=20, plot = FALSE)
    ovec <- tov$count
    
    tmp <- .getPoly(tov$breaks[-1],ovec)
    cnt <- tmp[2,]
    cnt <- log10(cnt)
    cnt[cnt < 0] <- 0
    
    pmax   <- .4
    pscale <- pmax*max(at)/(1 + max(cnt))
    polygon( tmp[1,], pscale*cnt, col=.getColor(specCol[specNames[j]],.3), lwd=2, 
             border=specCol[specNames[j]])
    
    pat <- max(cnt/pmax) + 1
    mpat <- floor( pat )
    
    pat <- seq(0, mpat, length= npat )
    lpat <- round(10^pat,-2)
    if(lpat[2] == 0)lpat <- round(10^pat,-1)
    if(lpat[2] == 0)lpat <- round(10^pat,0)
    apat <- pat*max(at)/mpat
    axis(side=4, at=apat, labels = lpat)
    abline(h=apat, col=.getColor(specCol[specNames[j]],.3),lty=2,
           lwd=2)
    
    if(POINTS)points(dspec,sqrt(fpp), pch=16, col=.getColor(plotCol[pspec],.3),
                     cex=.5)
    
    .shadeInterval( dtt, t(ar[2:3,]) , col=
                      .getColor('wheat', .1) )
    .shadeInterval( dtt, t(sr[2:3,]) , col=
                      .getColor('tan', .01) )
    
    lines(dtt, sr[1,], col='white', lwd=6)
    lines(dtt, sr[1,], lwd=2)
    
    if(nspec > 1).plotLabel(specNames[j])
    if(j == 1)title('Fecundity')
  }
  mtext('Diameter (cm)', side=1, line=1, outer  = TRUE)
  mtext('Probability', side=2, line=0, outer  = TRUE)
  mtext('log trees', side=4, line=0, outer  = TRUE)
  
  if(CONSOLE)
    readline('maturation, fecundity by diameter -- return to continue ')
  if(SAVEPLOTS)dev.off( )
  
  if(is.null(RMD)){
    graphics.off()
  }else{
    words <- paste( 'Predictive distribution for maturation, fecundity from posterior (shaded intervals) and latent states (dots) ' )
    message( words )
    caption <- c(caption, words)
  }
  
  ########### random coefficients
  
  if(RANDOM){
    
    xrandCols  <- output$data$setupRandom$xrandCols
    formulaRan <- output$data$setupRandom$formulaRan
    
    betaFec <- output$parameters$betaFec[,1, drop = FALSE]
    names(xrandCols) <- rownames(output$parameters$betaFec)[xrandCols]
    alphaUMu <- parameters$alphaUMu
    intercept <- paste('species',specNames,sep='')
    if(nspec == 1)intercept <- '(Intercept)'
    slopes    <- names(xrandCols)[!names(xrandCols) %in% intercept]
    
    rlim <- range(alphaUMu)
    
    if(rlim[2] > 0){
      
      npp <- length(slopes)/length(intercept)
      
      npp <- 1
      
      mfrow <- c(1,1)
      if(npp > 1)mfrow <- c(2,2)
      
      if(SAVEPLOTS)pdf( file=.outFile(outFolder,'randomCoeffs.pdf') )
      par(mfrow=mfrow, bty='n', mar=c(4,4,1,1), cex=1.2)
      
      icol <- match(intercept, colnames(alphaMu))
      icol <- icol[is.finite(icol)]
      xlim <- range(betaFec[icol,] + alphaMu[,icol])
      intercept <- intercept[is.finite(icol)]
      
      if(length(slopes) == 0){
        
        breaks <- seq(rlim[1]-.1, rlim[2]+.1,length=20)
        hvals <- xvals <- numeric(0)
        for(s in 1:length(specNames)){
          ra <- alphaUMu[,icol[s]]
          ra <- ra[ra != 0]
          tmp <- hist(ra, plot = FALSE, breaks = breaks)
          hvals <- cbind(hvals,tmp$density)
          xvals <- cbind(xvals,tmp$mids)
        }
        ylim <- range(hvals, na.rm  = TRUE)
        xl <- range(betaFec[colnames(alphaUMu),1],na.rm  = TRUE) + range(rlim, na.rm  = TRUE)
        xl[1] <- xl[1] - 1
        xl[2] <- xl[2] + 1
        
        plot(NA,xlim=xl,ylim = ylim, xlab='log fecundity (fixed plus random)', ylab='frequency')
        for(s in 1:length(specNames)){
          
          hp <- which(hvals[,s] > 0)
          if(length(hp) == 0)next
          
          ws <- suppressWarnings( range( hp ) )
          
          
          
          ss <- ws[1]:ws[2]
          xs <- betaFec[icol[s],1] + xvals[ss,1]
          xs <- c(xs[1],xs,xs[length(xs)])
          ys <- c(0, hvals[ss,s], 0)
          lines( xs, ys,type='s',col=specCol[s], lwd=2)
          lines( xs, xs*0,col=specCol[s], lwd=2)
        }
      }else{
        
        ww <- grep(':', slopes)
        
        if(length(ww) > 0){
          slopeLab <- columnSplit(slopes, ':')[,2]
        }else{
          slopeLab <- slopes
        }
        jcol <- match(slopes, colnames(alphaMu))
        ylim <- range(betaFec[jcol,1] + alphaMu[,jcol])
        
        if(length(icol) == 1)jcol <- jcol[1]
        
        slab <- unique(slopeLab)
        jcol <- jcol[slopeLab == slab[1]]
        
        plot(betaFec[icol,1],betaFec[jcol,1],col=.getColor(specCol,.7), 
             xlim=xlim, ylim=ylim, cex=.8, xlab='Intercept',ylab=slopeLab[1])
        abline(h = 0, lwd=2, lty=2, col='grey')
        abline(v = 0, lwd=2, lty=2, col='grey')
        
        for(s in 1:length(specNames)){
          points(betaFec[icol[s],1] + alphaMu[,icol[s]],
                 betaFec[jcol[s],1] + alphaMu[,jcol[s]],
                 col=.getColor(specCol[s],.3), pch=16, cex=1)
        }
      }
      legend('topright', specNames, text.col = specCol,bty='n')
      
      if(CONSOLE)
        readline('fixed plus random effects -- return to continue ')
      if(SAVEPLOTS)dev.off( )
      
      if(is.null(RMD)){
        graphics.off()
      }else{
        words <- 'Posterior estimates of fixed plus random effects'
        message( words )
        caption <- c(caption, words)
      }
    }
  }
  
  ############ dispersal
  
  if(SEEDDATA){
    
    cols  <- specCol[colnames(ugibbs)]
    upars <- parameters$upars[drop=FALSE, colnames(ugibbs),]
    
    if(ncol(ugibbs) > 1){
      
      if( is.null(rownames(upars)) )rownames(upars) <- 'mean'
      
      if(SAVEPLOTS)pdf( file=.outFile(outFolder,'dispersalCoeffs.pdf') )
      
      ncc <- ncol(ugibbs)
      sideMar <- min( c(3, 3 + 10/ncc) )
      
      par(mfrow=c(1,1), mar=c(3,5,5,4), bty='n', oma=c(2,sideMar,1,sideMar))
      
      ylim <- range(ugibbs)
      ylim[1] <- .75*ylim[1]
      ylim[2] <- 1.25*ylim[2]
      
      ylab1 <- expression( paste(hat(u),' (', plain(m)^2, ')') )
      ylab2 <- expression( paste('Kernel mean ', bar(d), ' (m)') )
      
      # uord <- order( colMeans(ugibbs, na.rm=T) )
      # colu <- cols[uord]
      
      atvals <- c(1:ncc)/(ncc + 1)
      atvals <- atvals - mean(atvals)
      
      tmp <- .boxplotQuant( ugibbs, xaxt='n', add = F, xlim = NULL, 
                            ylab=ylab1, ylim = ylim, boxwex=.6,
                            at = atvals,
                            outline = FALSE, col=.getColor(cols,.3), 
                            border=cols, lty=1, boxfill=NULL )
      tmp$xtick <- atvals
      .boxCoeffsLabs( boxPars = tmp, labels = names(cols), 
                      colLabs = cols, cex=.8)
      
      rr <- range( pi*sqrt(ugibbs)/2 )
      rm <- -1
      by <- 10
      if(diff(rr) < 20){
        by <- 5
        rm <- 0
      }
      if(diff(rr) < 10){
        by <- 2
        rr[1] <- rr[1] - 1
        rr[2] <- rr[2] + 1
      }
      
      rr <- round(rr, rm)
      rr <- seq(rr[1],rr[2], by=by)
      
      axis(4, at = (2*rr/pi)^2, labels = rr )
      mtext(ylab2, 4, line=0, outer  = TRUE)
      abline(v=par('usr')[1:2], lwd=2)
      
      if(CONSOLE)
        readline('dispersal by group -- return to continue ')
      if(SAVEPLOTS)dev.off( )
      if(is.null(RMD)){
        graphics.off()
      }else{
        words <- paste( 'Posterior estimates for dispersal parameters' )
        message( words )
        caption <- c(caption, words)
      }
    }
    
    if(is.null(RMD)) graphics.off()
    
    if(SAVEPLOTS)pdf( file=.outFile(outFolder,'seedShadow.pdf') )
    
    par(mfrow=c(1,1), bty='n')
    
    xfec  <- data$setupData$xfec
    
    qd     <- .75
    dcol  <- grep(':diam',colnames(xfec))
    dtcol <- grep('diam',colnames(tdata))
    lab60 <- tdata[,dtcol]
    lab60 <- lab60[lab60 != 0]
    lab60 <- quantile(lab60, qd)
    lab60 <- signif(lab60,1)
    if(lab60 < 10)lab60 <- 10
    ord <- order(tdata[,dtcol])
    sdd <- tdata[ord,dtcol]
    dd <- findInterval(lab60, sdd)
    lab60 <- sdd[dd]
    diam60 <- xfec[ord[dd],dcol]
    diam60 <- diam60[diam60 != 0]
    
    # each group at mean for other predictors
    xbar   <- numeric(0)
    rnames <- character(0)
    ucol   <- numeric(0)
    
    snames <- colnames(ugibbs)
    
    for(j in 1:ncol(ugibbs)){
      
      wj   <- which(tdata$species == snames[j])
      xmu  <- colMeans(xfec[drop = FALSE,wj,])
      
      w0 <- which(xmu != 0)       #insert diameter value
      w0 <- intersect(w0, dcol)
      xmu[w0] <- diam60
      
      xbar <- rbind(xbar, xmu)
      rnames <- c(rnames, snames[j])
    }
    
    rownames(xbar) <- rnames
    if(!USPEC)xbar <- xbar[drop = FALSE, 1, ]
    
    nsim <- 2000
    ns  <- 100
    buffer <- 15
    
    dpars <- parameters$dpars[drop=FALSE, colnames(ugibbs),]
    
    mm <- 6*max(dpars)
    if(USPEC)mm <- 6*max(dpars[,'estimate'])
    
    dseq <- 10^seq(0,log10(mm),length=ns) - 1
    dseq <- matrix( dseq, ns, nsim)
    
    ij <- sample(1:nrow(ugibbs), nsim, replace  = TRUE)
    
    ssList <- numeric(0)
    maxy   <- 0
    
    bf <- bfecStnd[ij,colnames(xbar)]
    
    if(RANDOM){
      br <- rmvnormRcpp(nrow(bf), aMu[1,]*0, aMu)
      colnames(br) <- colnames(aMu)
      bf[,colnames(br)] <- bf[,colnames(br)] + br
    }
    ff <-  bf%*%t(xbar)
    ff <- ff[,colnames(ugibbs),drop=FALSE]
    xbar <- xbar[drop=F, colnames(ugibbs),]
    
    maxy <- numeric(0)
    
    kss <- 1:(ns-buffer)
    keepSeq <- c(rev(kss), kss) 
    dss  <- c(-rev(dseq[kss,1]),dseq[kss,1])
    
    for(k in 1:ncol(ugibbs)){
      
      uj <- ugibbs[ij,k]
      sj <- sgibbs[ij,'sigma']
      
      kj <- uj/pi/(uj + dseq^2)^2
      kj <- kj*matrix(exp(ff[,k] + sj/2), ns, nsim, byrow  = TRUE)
      kj[!is.finite(kj)] <- NA
      qj <- t( apply(kj, 1, quantile, c(.5, .05, .95), na.rm  = TRUE ) )
      
      for(m in 1:3)qj[,m] <- runmed(qj[,m], k = buffer, endrule='constant')
      maxy <- c(maxy, max( qj[,1] ))
      
      ssList <- append(ssList, list( qj[keepSeq,] ) )
    }
    
    maxy[maxy > 1e+10] <- 1e+10
    maxy[maxy < 1e-10] <- 1e-10
    names(ssList) <- rownames(xbar)
    
    par(bty='n', mar=c(5,5,1,1))
    
    labSeeds <- expression( paste('Seeds (', plain(m)^-2, ')') )
    
    if(USPEC){
      rmax <- diff( range(log10(maxy)) )
    }else{
      rmax <- log10(maxy)
    }
    SQRT <- FALSE
    if( rmax > 3 ){
      SQRT <- TRUE
      smax <- sqrt(max(maxy))
      plot(NULL, xlim=range(dss), ylim=c(0,smax),
           xlab='Distance (m)', ylab = labSeeds, yaxt='n')
      tt   <- sqrtSeq(1.2*smax)
      at   <- tt$at
      labs <- tt$labs
      axis(2, at = at, labels = labs)
      
    }else{
      plot(NULL, xlim=range(dss), ylim=c(0,1.2*max(maxy)),
           xlab='Distance (m)', ylab = labSeeds)
    }
    
    for(k in 1:length(ssList)){
      if(maxy[k] <= 1e-10 | maxy[k] >= 1e+10)next
      ss <- ssList[[k]][,2:3]
      
      if(SQRT)ss <- sqrt(ss)
      .shadeInterval( dss, ss , col=
                        .getColor(cols[names(cols)[k]], .2) )
    }
    mc <- numeric(0)
    for(k in 1:length(ssList)){
      ss <- ssList[[k]][,1]
      if( max(ss,na.rm  = TRUE) >= 1e+10)next
      if(SQRT)ss <- sqrt(ss)
      lines(dss, ss, col='white',lwd=6)
      lines(dss, ss, col=cols[names(cols)[k]], lwd=2)
      mc <- c(mc, max(ssList[[k]][,1]))
    }
    
    if(USPEC){
      ord <- order(mc, decreasing  = TRUE)
      legend('topright',names(cols)[ord],
             text.col=cols[names(cols)[ord]], bty='n')
    }
    legend('topleft', paste(round(lab60),' cm diameter tree',sep=''), bty='n')
    
    if(CONSOLE)
      readline('seed shadow -- return to continue ')
    if(SAVEPLOTS)dev.off( )
    
    if(is.null(RMD)){
      graphics.off()
    }else{
      words <- paste( 'Seed shadow predictions for a ', round(lab60),' cm diameter tree',sep='')
      message( words )
      caption <- c(caption, words)
    }
    
    
    ########### predicted seed counts, true fecundity
    
    xs <- rowSums(sdata[sdata$obs == 1,seedNames,drop = FALSE])
    pcols <- grep('predMean',colnames(seedPred))      # by seed trap
    ecols <- grep('estMean',colnames(seedPred))
    
    ys <- seedPred[,pcols]
    zs <- seedPred[,ecols]
    
    ww <- which(is.finite(xs) & xs > 0)
    
    if(length(ww) > 10){
      
      if(SAVEPLOTS)pdf( file=.outFile(outFolder,'seedPrediction.pdf') )
      
      xs <- xs[ww]
      ys <- ys[ww]
      zs <- zs[ww]
      
      ylim <- quantile(c(ys, zs), c(0,.99),na.rm  = TRUE)
      ylim[2] <- ylim[2]*1.5
      xlim <- range(xs, na.rm  = TRUE) + 1
      
      mfrow <- c(1,2)
      title <- 'a) From posterior distribution'
      
      if(TV)mfrow <- c(2,2)
      if(SEEDCENSOR)mfrow <- c(2,2)
      
      par(mfrow=mfrow, mar=c(4,4,2,1), bty='l')
      
      bins <- getBins(xs, nbin=10, pow=.4)
      nbin <- length(bins)
      
      opt <- list(log = FALSE, xlabel='Observed', bins = bins,
                  nbin=nbin, ylabel='Predicted', col='forestgreen', 
                  ylimit = ylim, xlimit = xlim, SQRT  = TRUE)
      tmp <- .plotObsPred(xs, ys, opt = opt)
      .plotLabel(title, above  = TRUE, cex=1)
      abline(0, 1, lwd=2, col='white')
      abline(0,1,lty=2)
      abline(h = mean(ys, na.rm=T),lty=2)
      if(SEEDCENSOR){
        cs <- sqrt(seedPred[rownames(censMin),pcols])
        xl <- rowSums(censMin[,-1, drop = FALSE])
        xh <- rowSums(censMax[,-1, drop = FALSE])
        xh[xh > xlim[2]] <- xlim[2]
        opt <- list(log = FALSE, xlabel='Censored interval', bins = bins,
                    nbin=nbin, ylabel=' ', col='white', ptcol='white',
                    ylimit = ylim, xlimit = xlim, SQRT  = TRUE)
        tmp <- .plotObsPred(xs, ys, opt = opt)
        segments(xl, cs, xh, cs, col=.getColor('black',.1))
        abline(0, 1, lwd=2, col='white')
        abline(0,1,lty=2)
        .plotLabel('b) Censored observations', above  = TRUE, cex=1)
      }
      
      lab <- 'b) From fecundity estimate'
      if(SEEDCENSOR)lab <- 'c) From fecundity estimate'
      
      opt <- list(log = FALSE, xlabel='Observed', bins = bins, atx = tmp$atx, labx = tmp$labx,
                  aty = tmp$aty, laby = tmp$laby,
                  ylabel='', col='forestgreen', #ylimit=ylim, xlimit = xlim,
                  nbin=nbin, SQRT  = TRUE)
      tmp <- .plotObsPred(xs, zs, opt = opt)
      .plotLabel(lab, above  = TRUE, cex=.8)
      abline(0, 1, lwd=2, col='white')
      abline(0,1,lty=2)
      abline(h = mean(zs, na.rm=T),lty=2)
      
      cwords <- 'prediction'
      
      if(SEEDCENSOR){
        cs <- sqrt(seedPred[rownames(censMin),ecols])
        xl <- rowSums(censMin[,-1, drop = FALSE])
        xh <- rowSums(censMax[,-1, drop = FALSE])
        xh[xh > xlim[2]] <- xlim[2]
        opt <- list(log = FALSE, xlabel='Censored interval', bins = bins,
                    nbin=nbin, ylabel=' ', col='white', ptcol='white',
                    ylimit = ylim, xlimit = xlim, SQRT  = TRUE)
        tmp <- .plotObsPred(xs, ys, opt = opt)
        segments(xl, cs, xh, cs, col=.getColor('black',.1))
        .plotLabel('d) Censored from fecundity estimate', above  = TRUE, cex=.9)
        abline(0, 1, lwd=2, col='white')
        abline(0,1,lty=2)
      }
      
      if(TV){
        
        xs <- inputs$trueValues$fec
        ys <- prediction$fecPred$fecEstMu
        names(ys) <- rownames(prediction$fecPred)
        
        anames <- intersect(names(xs),names(ys))
        anames <- anames[anames %in% names(xs) & anames %in% names(ys)]
        ys <- ys[anames]
        xs <- xs[anames]
        
        ylim <- quantile(ys[ys > 1],c(.02,.99))
        ylim[1] <- max(c(ylim[1],1))
        
        ws <- which(xs > 0 & ys > 0)
        xlim <- quantile(xs[ws],c(.02,.99))
        
        bins <- getBins(xs, pow = .2)
        nbin <- length(bins)
        
        opt <- list(xlimit = xlim, ylimit = ylim, bins = bins,
                    nbin=nbin, log = FALSE, xlabel='True values', 
                    ylabel='Estimates', col='darkgreen', SQRT  = TRUE)
        .plotObsPred(xs, ys, opt = opt)
        .plotLabel('c) Fecundity prediction', above  = TRUE, cex=.8)
        abline(0, 1, lwd=2, col='white')
        abline(0,1,lty=2)
        
        xs <- inputs$trueValues$repr
        ys <- prediction$fecPred$matrEst
        names(ys) <- rownames(prediction$fecPred)
        
        anames <- intersect(names(xs),names(ys))
        anames <- anames[anames %in% names(xs) & anames %in% names(ys)]
        ys <- ys[anames]
        xs <- xs[anames]
        
        tmp <- boxplot( ys ~ xs, plot = FALSE)
        
        ss  <- apply( prediction$fecPred[,c('matrEst','matrPred')], 2, quantile, 
                      pnorm(c(-1.96,-1,0,1,1.96)), na.rm  = TRUE ) 
        ss  <- ss[names(xs)]
        
        ps <- c(.025, .05, .1, .157)
        ps <- c(ps, .5, 1 - ps)
        qs <- unlist( by(ys, INDICES=xs, FUN=quantile, ps ) )
        
        
        tmp$stats <- matrix(qs, length(ps), 2)
        bxp(tmp, add = F, yaxt = 'n',  varwidth  = TRUE, xlab='True values',
            ylim=c(0,1), ylab='', boxwex=.5,
            outline = FALSE, col=.getColor('black', .2), 
            border='black', lty=1, boxfill=NULL)
        axis(2, at=c(0,1))
        .plotLabel('d) Maturation prediction', above  = TRUE, cex=.8)
      }
      
      if(SEEDDATA | TV){
        
        if(CONSOLE)
          readline( paste(cwords, '-- return to continue') )
        if(SAVEPLOTS)dev.off( )
        
        if(is.null(RMD)){
          graphics.off()
        }else{
          words <- 'Prediction from the posterior distribition'
          message( words )
          caption <- c(caption, words)
        }
      }
    } 
  }###################### end plot obs/pred seedData
  
  
  wc    <- which(is.finite(tdata$cropCount))
  if(length(wc) == 0)CONES <- FALSE
  
  if( CONES | 'cropMin' %in% colnames(tdata) ){
    
    if(SAVEPLOTS)pdf( file=.outFile(outFolder,'cropTrees.pdf') )
    
    par(mfrow=c(1,1),bty='n', mar=c(4,4,2,1), omi=c(.8,.2,0,.5))
    
    cobs  <- seedTraits[tdata$species,'seedsPerFruit']*tdata$cropCount/tdata$cropFraction
    
    if('cropMin' %in% colnames(tdata)){
      wcrop <- which(is.finite(tdata$cropMin))
      cobs[wcrop] <- (tdata$fecMin[wcrop] + tdata$fecMax[wcrop])/2 #  mean
      cobs[wcrop][!is.finite(cobs[wcrop])] <- tdata$fecMin[wcrop][!is.finite(cobs[wcrop])]
    }
    
    wc    <- which(is.finite(cobs))
    cobs  <- sqrt(cobs[wc])                    #sqrt scale
    cest  <- prediction$fecPred$fecEstMu[wc]
    cese  <- prediction$fecPred$fecEstSe[wc]
    cpre  <- prediction$fecPred$fecPredMu[wc]
    cpse  <- prediction$fecPred$fecPredSe[wc]
    
    cplo <- cpre - cpse
    cplo[cplo < 0] <- 0
    cplo <- sqrt(cplo)
    cphi <- sqrt(cpre + cpse)
    
    celo <- cest - cese 
    celo[celo < 0] <- 0
    celo <- sqrt(celo)
    cehi <- sqrt(cest + cese)
    
    xlim <- c(0, max(cobs))      # sqrt scale
    ylim <- sqrt( c(0, max(cpre + cpse)) )
    if(ylim[2] < max(cobs))ylim[2] <- max(cobs)
    
    ccens <- rep(0, nrow(tdata))
    ws    <- numeric(0)
    
    ccenlo <- ccenhi <- NA
    if('cropMin' %in% colnames(inputs$treeData)){
      ws <- which(is.finite(tdata$cropMin))
      ccenlo <- tdata$fecMin[ws]
      ccens  <- ccenlo
      sest  <- prediction$fecPred$fecEstMu[ws]
      sese  <- prediction$fecPred$fecEstSe[ws]
      selo  <- sest - sese
      selo[selo < 0] <- 0
      selo  <- sqrt(selo)
      selo  <- sqrt(selo)
    }
    if('cropMax' %in% colnames(inputs$treeData)){
      ws <- which(is.finite(tdata$cropMax))
      ccenhi <- tdata$fecMax[ws]
      ccens  <- (ccens + ccenhi)/2
      sehi  <- sqrt(sest + sese)
    }
    if('cropMin' %in% colnames(inputs$treeData) |
       'cropMax' %in% colnames(inputs$treeData)){
      sobs  <- sqrt(ccens)
      slo   <- sqrt(ccenlo)
      shi   <- sqrt(ccenhi)
      shi[ shi > max(sehi) ] <- max(sehi)
    }
    
    xlab <- expression( frac('count', 'crop fraction') %*% 'seeds per cone' )
      
  #    expression( paste('Kernel mean ', bar(d), ' (m)') )
    tt   <- sqrtSeq(1.2*xlim[2])
    atx   <- tt$at
    labsx <- tt$labs
    
    tt   <- sqrtSeq(1.2*ylim[2])
    aty   <- tt$at
    labsy <- tt$labs
 
    plot(cobs, sqrt(cpre), xlim = xlim, ylim = ylim, xaxt='n', yaxt='n',
         pch = 3, xlab='',ylab='Seeds per tree', col='thistle4')
    axis(1, at = atx, labels = labsx)
    axis(2, at = aty, labels = labsy)
    
    suppressWarnings(
      arrows( cobs, cplo, cobs, cphi, lwd=1,
              angle=90, length=.05, col='thistle4', code=3)
    )
    
    if(length(ws) > 0){
      suppressWarnings(
      arrows( sobs, selo, sobs, sehi, lwd=1,
              angle=90, length=.05, col=.getColor('mediumaquamarine',.8), code=3) )
     suppressWarnings(
      arrows( slo, sqrt(sest), shi, sqrt(sest), lwd=1,
              angle=90, length=.05, col=.getColor('mediumaquamarine',.8), code=3) )
    }
    
    points(cobs,sqrt( cest ), pch = 3)

    suppressWarnings(
      arrows( cobs, celo, cobs, sqrt(cest + cese) + .1, lwd=1,
              angle=90, length=.05, col=1, code=3)
    )
    abline(0, 1, lwd=2, col='grey',lty=2)
    leg <- c('Estimated','Predicted')
    text.col <- c('black','grey')
    if(length(ws) > 0){
      leg <- c('Estimated censored counts','Estimated crop counts','Predicted')
      text.col <- c(.getColor('mediumaquamarine',.8), 'black', 'thistle4')
    }
    legend('topleft', legend = leg, text.col = text.col,
           box.col = 'grey' )
    mtext(xlab, side=1, outer  = TRUE, line=0)
    
    if(CONSOLE)
      readline('Fecundity for cropCount trees -- return to continue ')
    if(SAVEPLOTS)dev.off( )
    if(is.null(RMD)){
      graphics.off()
    }else{
      words <- 'Estimated/predicted fecundity for trees with cone counts'
      message( words )
      caption <- c(caption, words)
    }
  }
  
  # true parameter values
  
  if(TV){
    
    brepTrue <- inputs$trueValues$betaRep
    bfecTrue <- inputs$trueValues$betaFec
    
    betaFec <- as.matrix(output$parameters$betaFec[,1:4])
    betaRep <- as.matrix(output$parameters$betaRep[,1:4])
    
    if(SAVEPLOTS)pdf( file=.outFile(outFolder,'trueParameters.pdf') )
    
    par(mfrow=c(1,2),bty='n', mar=c(3,3,3,4), oma=c(2,2,1,1))

    sc  <- grep('diam',rownames(betaRep))
    col <- rep('black', nrow(betaRep))
    col[-sc] <- 'gray'
    slp <- c(mean(inputs$trueValues$betaRep[sc,1]), mean(betaRep[sc,1]) )
    
    xlim <- range( c(brepTrue, betaRep[,3:4]) )
    
    plot( brepTrue, betaRep[,1], xlim = xlim, ylim = xlim, xlab='', 
          ylab=' ', pch=3, col = col)
    abline(0, 1, lwd=2, lty=2, col=.getColor('black',.4))
    suppressWarnings(
      arrows( inputs$trueValues$betaRep, betaRep[,3], 
              inputs$trueValues$betaRep, betaRep[,4], lwd=2,
              angle=90, length=.05, col=col, code=3)
    )
    text(slp[1],slp[2], 'Slopes', pos=2)
    .plotLabel('a) Maturation parameters', above  = TRUE, cex=1)
    
    sc  <- grep('diam',rownames(betaFec))
    col <- rep('black', length(betaFec))
    col[-sc] <- 'gray'
    slp <- c(mean(inputs$trueValues$betaFec[sc,1]), mean(betaFec[sc,1]) )
    
    bt <- bfecTrue
    bf <- betaFec
    
    ylim <- range( c(bt, bf[,-2]) )
    
    plot( bt, bf[,1],  ylim = ylim,
          xlab='', ylab='', pch=3, col=col)
    abline(0,1, lwd=2, lty=2, col=.getColor('black',.4))
    suppressWarnings(
      arrows( bt, bf[,3], 
              bt, bf[,4], lwd=2,
              angle=90, length=.05, code=3, col=col)
    )
    text(slp[1],slp[2], 'Slopes', pos=2)
    .plotLabel('b) Fecundity parameters', above  = TRUE, cex=1)
    mtext('True values', side=1, line=0, outer  = TRUE)
    mtext('Estimates', side=2, line=0, outer  = TRUE)
    
    if(CONSOLE)
      readline('parameter recovery -- return to continue ')
    if(SAVEPLOTS)dev.off( )
    if(is.null(RMD)){
      graphics.off()
    }else{
      words <- 'Posterior estimates compared with true values'
      message( words )
      caption <- c(caption, words)
    }
  }
  
  ################# year effects
  
  YE <- 'betaYr' %in% names(parameters)
  
  if(is.null(parameters$betaYr))YE <- YR <- FALSE
  
  if(YE)YE <- nrow(parameters$betaYr) > 1 
  if(YE){
    YEE <- 'betaYrRand' %in% names(parameters)
    if(YEE){
      if(sum(parameters$betaYrRand) == 0)YE <- YR <- FALSE
    }
  }
  
  if( YE ){
    
    if(is.null(RMD)) graphics.off()
    
    file <- paste('yearEffect.pdf',sep='')
    
    if(SAVEPLOTS)pdf( file=.outFile(outFolder,file) )
    
    if('lagGroup' %in% names(inputs))AR <- TRUE
    
    betaYrMu <- parameters$betaYr[,1, drop = FALSE]
    betaYrSe <- parameters$betaYr[,2, drop = FALSE]
    
    if(ncol(betaYrMu) == 1){
      betaYrMu <- t(betaYrMu)
      betaYrSe <- t(betaYrSe)
    }
    
    RANDYR <- FALSE
    
    if( 'betaYrRand' %in% names(parameters) ){  #combine fixed and random
      
      betaYrRand   <- parameters$betaYrRand
      betaYrRandSE <- parameters$betaYrRandSE
      betaYrRandSE[is.na(betaYrRandSE)] <- 0
      bmu <- bsd <- betaYrRand*0
      RANDYR <- TRUE
      
      if(!AR){
        ttab <- table( tdata$groupName, tdata$year )
      }else{
        ttab <- betaYrRand
      }
      ttab[ttab != 0] <- 1
      ctab <- betaYrRand*0
      ctab[,colnames(ttab)] <- ttab
      
      for(k in 1:nrow(betaYrRand)){
        
        bmu[k,] <- betaYrRand[k,]
        bsd[k,] <- betaYrRandSE[k,]
        bmu[k,] <- bmu[k,]*ctab[rownames(betaYrRand)[k],]
        bsd[k,] <- bsd[k,]*ctab[rownames(betaYrRand)[k],]
      }
      
      betaYrMu <- bmu
      betaYrSe <- bsd
      betaYrSe[is.na(betaYrSe)] <- 0
    }
    
    if(AR){
      par(mfrow=c(1,2),bty='n', mar=c(5,4,2,1))
      xlab <- 'lag (yr)'
      yr <- c(1:plag)
      mr <- .5
    }
    
    xtt <- seq(1900,2100,by=5)                #reference xtick
    
    par(mfrow=c(1,1),bty='n', mar=c(4,4,2,5), mai=c(1,1,1,1.1))
    if(AR){
      yr <- xtick <- 1:plag
      xlab <- 'lag (yr)'
    }
    if(YR & !AR){
      yr <- which(colSums(betaYrMu) != 0) 
      betaYrMu <- betaYrMu[,yr]
      betaYrSe <- betaYrSe[,yr]
      names(yr) <- .replaceString(names(yr),'yr-','')
      yr <- as.numeric(names(yr))
      xtick <- min(yr):max(yr)
      xlab <- ''
      if(length(yr) > 15)xtick <- xtick[xtick %in% xtt]
    }
    mr  <- max(betaYrMu + betaYrSe, na.rm  = TRUE)
    mm  <- max(betaYrMu - betaYrSe, na.rm  = TRUE)
    mr  <- max(c(abs(mr),abs(mm)))
    
    ylim = mr*c(-1.5, 1.5)
    
    plot(NULL, xlim = range(yr), ylim = ylim, xlab = xlab, 
         ylab = 'log fecundity', xaxt='n')
    axis(1, xtick)
    abline(h = 0, lty=2, col='grey', lwd=2)
    leg <- character(0); col <- numeric(0)
    
    if( (YR & !RANDYR) | (!AR & !YR) ){  # there is one group
      loHi <- cbind( betaYrMu - 1.96*betaYrSe,
                     betaYrMu + 1.96*betaYrSe)
      .shadeInterval(yr,loHi,col='black',PLOT  = TRUE,add  = TRUE, trans = .3)
      lines(yr, betaYrMu, col=.getColor('white',.7), lwd=5)
      lines(yr, betaYrMu, col='grey', lwd=2)
    }
    col <- numeric(0)
    
    if(!is.null(betaYrRand)){
      betaYrRand   <- betaYrRand[drop = FALSE,yeGr,]
      betaYrRandSE <- betaYrRandSE[drop = FALSE,yeGr,]
    }
    
    if(ngroup == 1){
      wj <- which(is.finite(betaYrMu) & betaYrMu != 0) 
      loHi <- cbind( betaYrMu[wj] - 1.96*betaYrSe[wj],
                     betaYrMu[wj] + 1.96*betaYrSe[wj])
      .shadeInterval(yr[wj],loHi,col=groupCol[1],PLOT  = TRUE,add  = TRUE, trans = .3)
      lines(yr[wj], betaYrMu[wj], col=.getColor('white',.7), lwd=5)
      lines(yr[wj], betaYrMu[wj], col=groupCol[1], lwd=2)
    }else{
      
      for(j in 1:ngroup){
        
        nj <- yeGr[j]
        wj <- which(is.finite(betaYrMu[nj,]) & betaYrMu[nj,] != 0) 
        if(length(wj) < 2)next
        kj <- wj[ which(diff(wj) == 1) ]
        dj <- c(1, which(diff(wj) > 1)) #, max(wj))
        sj <- kj[dj]
        ej <- c(wj[dj-1], max(wj))
        ej <- ej[is.finite(sj)]
        sj <- sj[is.finite(sj)]
        
        sj <- sj[ is.finite(ej) ]
        
        for(m in 1:length(sj)){
          mm   <- sj[m]:ej[m]
          if(length(mm) < 2)next
          loHi <- cbind( betaYrMu[nj,mm] - 1.96*betaYrSe[nj,mm],
                         betaYrMu[nj,mm] + 1.96*betaYrSe[nj,mm])
          .shadeInterval(yr[mm],loHi,col=groupCol[nj],PLOT  = TRUE,add  = TRUE, trans = .3)
        }
      }
      for(j in 1:ngroup){
        
        nj <- yeGr[j]
        wj <- which(is.finite(betaYrMu[nj,]) & betaYrMu[nj,] != 0) 
        if(length(wj) < 2)next
        kj <- wj[ which(diff(wj) == 1) ]
        dj <- c(1, which(diff(wj) > 1)) #, max(wj))
        sj <- kj[dj]
        ej <- c(wj[dj-1], max(wj))
        ej <- ej[is.finite(sj)]
        sj <- sj[is.finite(sj)]
        
        sj <- sj[ is.finite(ej) ]
        
        for(m in 1:length(sj)){
          mm   <- sj[m]:ej[m]
          lines(yr[mm], betaYrMu[nj,mm], col=.getColor('white',.7), lwd=5)
          lines(yr[mm], betaYrMu[nj,mm], col=groupCol[nj], lwd=2)
        }
      }
      cornerLegend('topright', yeGr, text.col = groupCol[yeGr],
                   cex=.8, bty='n')
    }
  }
  
  if(YE & AR){
    title('AR coefficients',adj=0, font.main=1, font.lab=1,
          cex.main=.9)
    
    if(CONSOLE)
      readline('lag effect groups -- return to continue ')
    if(SAVEPLOTS)dev.off()
    if(is.null(RMD)){
      graphics.off()
    }else{
      words <- 'Posterior estimates of lag effects, by random group'
      message( words )
      caption <- c(caption, words)
    }
  }
  
  if(YE)title('Year effects, +/- 1 se',adj=1, font.main=1, font.lab=1,
              cex.main=.9)
  
  if(YE | AR){
    if(CONSOLE)
      readline('year effect all groups -- return to continue ')
    if(SAVEPLOTS)dev.off( )
    if(is.null(RMD)){
      graphics.off()
    }else{
      words <- 'Posterior estimate of year effects, by random group'
      message( words )
      caption <- c(caption, words)
    }
  }
  
  if( (YR | AR) & 'groups' %in% names(yearEffect) ){  # by region
    
    if(is.null(RMD)) graphics.off()
    
    file <- 'yearEffectByGroup.pdf'
    if(AR)file <- 'lagEffectByGroup.pdf'
    
    if(SAVEPLOTS)pdf( file=.outFile(outFolder,file) )
    
    yeGr <- output$data$setupYear$yeGr
    if(is.null(yeGr))yeGr <- as.character(output$data$setupData$yeGr)
    region <- yearEffect$plotGroups
    spec   <- yearEffect$specGroups
    
    wd <- grep('-',yeGr)
    if(length(wd) > 0){
      spec   <- columnSplit(yeGr,sep='-')[,1]
      region <- columnSplit(yeGr,sep='-')[,2]
    }else{
      if(length(spec) > 0)spec <- yeGr
      if(length(region) > 0)region <- yeGr
    }
    regs <- unique(region)
    nreg <- length(regs)
    
    if(is.null(nreg) | nreg == 0){
      regs <- spec
      nreg <- length(spec)
    }
    
    if(nreg > 0){
      
      nr <- ceiling(nreg/2)
      
      par(mfrow=c(nr,2), bty='n', mar=c(2,2,.1,2), oma=c(2,3,1,1))  
      
      if(YR)yval <- 4
      if(AR)yval <- 1
      
      for(k in 1:nreg){
        
        wk <- which(region == regs[k])
        
        plot(NULL, xlim = range(yr), ylim = c(-yval,yval), xlab = xlab, 
             ylab = '', xaxt='n',yaxt='n')
        if(k == nreg)axis(1, xtick)
        axis(2,c(-yval,0,yval))
        abline(h = 0, lty=2, col='grey', lwd=2)
        
        wll <- numeric(0)
        
        for(j in wk){
          
          nj <- yeGr[j]
          wj <- which(is.finite(betaYrMu[nj,]) & betaYrMu[nj,] != 0) 
          if(length(wj) < 2)next
          wll  <- c(wll,j)
          wj <- which(is.finite(betaYrMu[nj,]) & betaYrMu[nj,] != 0) 
          kj <- wj[ which(diff(wj) == 1) ]
          dj <- c(1, which(diff(wj) > 1)) #, max(wj))
          sj <- kj[dj]
          ej <- c(wj[dj-1], max(wj))
          ej <- ej[is.finite(sj)]
          sj <- sj[is.finite(sj) & is.finite(ej)]
          
          for(m in 1:length(sj)){
            
            mm   <- sj[m]:ej[m]
            if(length(mm) < 2)next
            loHi <- cbind( betaYrMu[nj,mm] - 1.96*betaYrSe[nj,mm],
                           betaYrMu[nj,mm] + 1.96*betaYrSe[nj,mm])
            .shadeInterval(yr[mm],loHi,col=groupCol[nj],PLOT  = TRUE,add  = TRUE, trans = .3)
          }
        }
        for(j in wk){
          
          nj <- yeGr[j]
          
          ww <- which(is.na(betaYrMu[nj,]) | betaYrMu[nj,] == 0)
          yj <- yr
          bj <- betaYrMu[nj,]
          bj[ww] <- NA
          yj[ww] <- NA
          
          lines(yj, bj, col=.getColor('white',.7), lwd=5)
          lines(yj, bj, col=groupCol[nj], lwd=2)
          
        }
        
        if(length(wll) == 0){
          par(new  = TRUE)
          next
        }
        
        legend('bottomleft', yeGr[wll], text.col = groupCol[yeGr[wll]],
               cex=1., bty='n')
        .plotLabel(region[j],'topleft')
        
        par(new = F)
      }
      
      if(AR){
        mtext('Lag',side=1,line=1,outer  = TRUE)
      }else{
        mtext('Year',side=1,line=1,outer  = TRUE)
      }
      mtext('log fecundity',side=2,line=1,outer  = TRUE)
    }
    
    
    if(CONSOLE)
      readline('year effect by group -- return to continue ')
    if(SAVEPLOTS)dev.off( )
    if(is.null(RMD)){
      graphics.off()
    }else{
      words <- 'Posterior estimate of year effects, by group'
      message( words )
      caption <- c(caption, words)
    }
  }
  
  ############# pacf
  
 # pacfMat <- output$parameters$pacfMat
 # pacsMat <- output$parameters$pacsMat
  
  
  nyy <- ceiling(nyr*.6)
  if(nyy > 10)nyy <- 10
  
  if( nyr > 3 & sum(pacfMat,na.rm  = TRUE) != 0 ){
    
    if(is.null(RMD)) graphics.off()
    
    file <- paste('pacf.pdf',sep='')
    if(SAVEPLOTS)pdf( file=.outFile(outFolder,file) )
    
    par(mfrow=c(1,2),bty='n', mar=c(3,2,1,.5), oma=c(2,3,1,1))
    
    mr  <- .9
    ylim <- range(pacfMat[,-1,drop=F],na.rm  = TRUE)
    if(SEEDDATA)ylim <- range( c(ylim, pacsMat[-1]), na.rm  = TRUE )
    
    plot(NULL, xlim = c(1,nyy), ylim = ylim, xaxt = 'n',
         xlab = '', ylab = '')
    axis(1, at=c(1:nyy))
    abline(h = 0, lty=2, col='grey', lwd=2)
    leg <- character(0); col <- numeric(0)
    lag <- c(0:(nyr-1))
    
    if(yeGr[1] %in% specNames & specNames[1] %in% rownames(pacfMat)){
      pacfMat <- pacfMat[drop = FALSE,specNames,]
      cols <- specCol[specNames]
      leg  <- specNames
    }else{
      cols <- groupCol
      leg <- rownames(pacfMat)
    }
    
    pacCol <- gfun(nrow(pacfMat))
    names(pacCol) <- rownames(pacfMat)
    
    for(j in 1:nrow(pacfMat)){
      
      wj <- which( is.finite(pacfMat[j,]))
      wj <- wj[-1]                   # omit lag zero
      nj <- length(wj)
      if(nj < 2)next
      
      ac <- pacfMat[j,wj]
      
      lines(lag[wj], ac, col=pacCol[j], lwd=2)
      loHi <- cbind( ac - 1.96*pacfSe[j,wj],
                     ac + 1.96*pacfSe[j,wj])
      .shadeInterval(lag[wj],loHi,col=pacCol[j],PLOT  = TRUE,add  = TRUE, 
                     trans = .2)
      
      up <- which(loHi[,1] > 0)
      up <- up[up > 1]
      points(lag[wj[up]],ac[up],cex=1.3,pch=16,col=.getColor(pacCol[j],.5))
      
      up <- up[up < 10]
      up <- paste0( up[up > 2], collapse=', ')
      ll <- rownames(pacfMat)[j]
      leg <- c(leg,ll)
      col <- c(col,j)
    }
    if(length(leg) > 1)legend('topright', leg, text.col = pacCol, 
                              bty='n', cex=.6)
    .plotLabel('a) log Fecundity', location='topleft',above  = TRUE, cex=.9)
    
    if(SEEDDATA){
      xlim <-  c(1,nyy)
      plot(NULL, xlim = xlim, ylim = ylim, xaxt = 'n', yaxt='n',
           xlab = '', ylab = '')
      axis(1, at=c(1:nyy))
      axis(2, labels = FALSE)
      abline(h=0, lty=2,lwd=2, col='grey')
      leg <- character(0); col <- numeric(0)
      lag <- c(0:nyy)
      
      wj <- which(pacsMat != 0)
      wj <- wj[wj != 1]
      nj <- length(wj)
      
      if(nj > 2){
        
        ac <- pacsMat[wj]
        lines(lag[wj], ac, col=1, lwd=2)
        segments( lag[wj], ac*0, lag[wj], ac)
      }
      .plotLabel('b) Seed counts', location='topleft', above  = TRUE, cex=.9)
      mtext('Lag (yr)', side=1, outer  = TRUE)
      mtext('PACF', side=2, outer  = TRUE, line=1)
    }
    
    if(CONSOLE)
      readline('partial ACF -- return to continue ')
    if(SAVEPLOTS)dev.off( )
    if(is.null(RMD)){
      graphics.off()
    }else{
      words <- 'Partial ACF'
      message( words )
      caption <- c(caption, words)
    }
  }
  
  if( 'groups' %in% names(yearEffect) & sum(pacfMat,na.rm  = TRUE) != 0 ){
    
    file <- paste('pacfByGroup.pdf',sep='')
    
    if(is.null(RMD)) graphics.off()
    
    if(SAVEPLOTS)pdf( file=.outFile(outFolder,file) )
    
    mfrow <- .getPlotLayout(nspec)
    par(mfrow=mfrow$mfrow, bty='n', mar=c(1,1,1,1), oma=c(3,3,1,3))  
    
    preg <- columnSplit(rownames(pacfMat),'-')
    
    for(k in 1:nspec){
      
      wk <- which(preg[,1] == specNames[k])
      
      ylab <- xlab <- FALSE
      if(k %in% mfrow$left)ylab  = TRUE
      if(k %in% mfrow$bottom)xlab  = TRUE
      
      plot(NULL, xlim = c(1,nyy), ylim = ylim, xaxt = 'n', yaxt='n',
           xlab = '', ylab = '')
      axis(1, at=c(1:nyy), labels=xlab)
      axis(2,labels=ylab)
      abline(h = 0, lty=2, col='grey', lwd=2)
      leg <- character(0); col <- numeric(0)
      lag <- c(0:(nyr-1))
      
      for(j in wk){
        wj <- which( is.finite(pacfMat[j,]))
        wj <- wj[-1]                   # omit lag zero
        nj <- length(wj)
        if(nj < 2)next
        
        ac <- pacfMat[j,wj]
        
        lines(lag[wj], ac, col=plotCol[preg[j,2]], lwd=2)
        loHi <- cbind( ac - 1.96*pacfSe[j,wj],
                       ac + 1.96*pacfSe[j,wj])
        .shadeInterval(lag[wj],loHi,col=plotCol[preg[j,2]],PLOT  = TRUE,add  = TRUE, 
                       trans = .4)
        
        up <- which(loHi[,1] > 0)
        up <- up[up > 1]
        points(lag[wj[up]],ac[up],cex=1.3,pch=16,
               col=.getColor(plotCol[preg[j,2]],.5))
        
        up <- up[up < 10]
        up <- paste0( up[up > 2], collapse=', ')
        ll <- rownames(pacfMat)[j]
        #    leg <- c(leg,ll)
        #    col <- c(col,j)
      }
      .plotLabel(specNames[k],'topright') 
    }
    if( nspec == prod(mfrow$mfrow) ){
      cornerLegend('bottomright',plots,bty='n',text.col=plotCol[plots], cex=.9)
    }else{
      plot(NULL, xlim=xlim, ylim=ylim, xlab='', ylab='', axes = FALSE)
      legend('topleft',plots,bty='n',text.col=plotCol[plots], cex=1.2)
    }
    mtext('Lag (yr)',side=1,line=1,outer  = TRUE)
    mtext('PACF',side=2,line=1,outer  = TRUE)
    
    if(CONSOLE)
      readline('partial ACF -- return to continue ')
    if(SAVEPLOTS)dev.off( )
    if(is.null(RMD)){
      graphics.off()
    }else{
      words <- 'Partial ACF'
      message( words )
      caption <- c(caption, words)
    }
  }
  
  ############# eigenvalues AR
  
  if(AR){
    
    if(is.null(RMD)) graphics.off()
    
    file <- paste('eigenAR.pdf',sep='')
    if(SAVEPLOTS)pdf( file=.outFile(outFolder,file) )
    
    par(bty='n', mai=c(1,1,1,1),mar=c(5,5,1,1), mfrow=c(1,1))
    
    ename <-rownames(eigenMu)
    ename <- rep(ename,each=ncol(eigenMu))
    text.col <- rep(groupCol[yeGr],each=ncol(eigenMu))
    
    xlab <- expression( paste( plain(Re),' ',lambda ))
    ylab <- expression( paste( plain(Im),' ',lambda ))
    
    xseq <- seq(-1,1,length=100)
    yseq <- sqrt(1 - xseq^2)
    plot(eigenMu,xlim=c(-1.2,1.2),ylim=c(-1.1,1.1), 
         cex=.1,xlab=xlab,ylab=ylab)
    lines(xseq,yseq,lwd=2,col='grey',lt=2)
    lines(xseq,-yseq,lwd=2,col='grey',lt=2)
    lines(c(0,0),c(-1,1),col='grey',lt=2)
    lines(c(-1,1),c(0,0),col='grey',lt=2)
    text(Re(eigenMu),Im(eigenMu),ename, cex=.9, col=text.col)
    
    if(CONSOLE)
      readline('ACF eigenvalues -- return to continue ')
    if(SAVEPLOTS)dev.off( )
    if(is.null(RMD)){
      graphics.off()
    }else{
      words <- 'ACF eigenvalues'
      message( words )
      caption <- c(caption, words)
    }
  }
  
  ############# fecundity and seed prediction
  
  if(SEEDDATA){
    
    tyears <- years  <- sort(unique(sdata$year))
    tplots <- pplots <- sort(unique(as.character(sdata$plot)))
    if(AR)tyears <- sort(unique(tdata$year))
    
    if(PREDICT & length(inputs$predList$years) > 1){
      
      file <- paste('prediction.pdf',sep='')
      if(SAVEPLOTS)pdf( file=.outFile(outFolder,file) )
      
      pyears <- sort(unique(c(fecPred$year,seedPredGrid$year)))
      pplots <- sort(unique(as.character(seedPredGrid$plot)))
      tyears <- sort(unique(c(years,pyears,tyears)))
      tplots <- sort(unique(c(plots,pplots)))
      
      # prediction grid closest to traps
      
      pcol <- grep('meanM2',colnames(seedPredGrid))
      scol <- grep('seM2',colnames(seedPredGrid))
      
      kcol <- c("trapID","plot","year","trap","plotYr","plotyr",
                "drow","area","active")
      spred <- numeric(0)
      
      for(j in 1:length(pplots)){
        
        xyt <- xytrap[xytrap$plot == pplots[j],]
        
        sp <- seedPredGrid[seedPredGrid$plot == pplots[j],]
        so <- sdata[sdata$plot == pplots[j],]
        if(nrow(so) == 0)next
        sxy <- xyt[match(so$trapID,xyt$trapID),c('x','y')]
        
        spi <- as.character( columnPaste(sp$trapID,sp$year) )
        soo <- as.character( columnPaste(so$trapID,so$year) )
        spi <- match(soo,spi)
        wf  <- which(is.finite(spi))
        
        spp <- sp[spi[wf],pcol,drop = FALSE]
        countPerM2 <- rowSums(so[wf,seedNames,drop = FALSE])/so$area[wf]
        predPerM2  <- rowSums(spp)
        
        sall  <- cbind(so[wf,],spp,countPerM2,predPerM2)
        spred <- rbind(spred,sall)
      }
      
      #error by year
      rmse <- sqrt( (spred$countPerM2 - spred$predPerM2)^2 )
      aerr <- spred$countPerM2 - spred$predPerM2
      
      xlim <- range(spred$year,na.rm  = TRUE)
      
      maxMod <- NULL
      if(!is.null(modelYears))maxMod <- max(modelYears)
      
      pplots <- as.character(sort(unique(spred$plot)))
      
      mfrow <- .getPlotLayout(length(pplots))$mfrow
      opt   <- list(log = FALSE, xlabel='Year', POINTS = FALSE, 
                    ylabel='Residual', col='brown', add  = TRUE)
      
      par(mfrow=mfrow,bty='n', mar=c(3,3,1,1), oma=c(2,2,0,2))
      
      xseq <- c(0,2^c(0:15))[-2]
      
      ylabel <- expression( paste('Count (', plain(m)^-2,')') )
      zlabel <- expression( bar(y) %+-% plain(sd) )
      
      npl <- 0
      
      for(j in 1:length(pplots)){
        
        wj    <- which(spred$plot == pplots[j])
        obs   <- spred$year[wj]
        yMean <- spred$predPerM2[wj]
        yObs  <- spred$countPerM2[wj]
        
        if( max(yMean, na.rm  = TRUE) == 0 | length(yMean) < 3 )next
        
        tj <- by(yMean, obs, quantile, probs=pnorm(c(0,-1,1)), na.rm  = TRUE)
        cj <- names(tj)
        tj <- matrix( unlist(tj), ncol=3, byrow  = TRUE )
        rownames(tj) <- cj
        tj <- sqrt(tj)     
        ww <- which(is.finite(tj[,1]))
        
        if(length(ww) < 2)next
        
        yj <- as.numeric(rownames(tj))
        
        omu <- by(yObs, obs, quantile, probs=pnorm(c(0,-1,1)), na.rm  = TRUE)
        cj <- names(omu)
        oj <- matrix( unlist(omu), ncol=3, byrow  = TRUE )
        rownames(oj) <- cj
        oj <- sqrt(oj)     # not for residuals
        
        smax <- max( c(tj,oj,5) )
        tt   <- sqrtSeq(smax)
        at   <- tt$at
        labs <- tt$labs
        
        #  xlim <- range(obs,na.rm  = TRUE)
        ylim <- range(at)
        
        plot(NULL,xlim=xlim,ylim=ylim, ylab='', xlab='',yaxt='n')
        axis(2,at=at, labels = labs)
        
        if(!is.null(maxMod)){
          rect(maxMod+.5,ylim[1],xlim[2]+.5,ylim[2],col='wheat', border='wheat')
        }
        
        .shadeInterval(yj,loHi=tj[drop = FALSE,ww,2:3], col=.getColor('grey', .8))
        abline(h=0,lty=2,lwd=4,col='white')
        
        .shadeInterval(yj,loHi=oj[drop = FALSE,ww,2:3], col=.getColor('green', .3))
        
        lines(yj, tj[ww,1], col='white', lwd=8)
        lines(yj, tj[ww,1], lwd=3)
        lines(yj, oj[ww,1], col=.getColor('white',.5), lwd=8)
        lines(yj, oj[ww,1], col=.getColor('forestgreen',.7),lwd=3)
        
        points(jitter(obs),sqrt(yObs),pch=16,col=.getColor('forestgreen',.2))
        
        .plotLabel(pplots[j],'topleft')
        
        npl <- npl + 1
      }
      
      if(npl > 0){
        mtext('Year',side=1,line=0,outer  = TRUE, cex=1.4)
        mtext(ylabel,side=2,line=0,outer  = TRUE, cex=1.4)
        mtext(zlabel,side=4,line=0,outer  = TRUE, cex=1.4)
        
        if(CONSOLE)
          readline('observed (green), predicted (black), shaded forecast (if modelYears) -- return to continue ')
        if(SAVEPLOTS)dev.off( )
        if(is.null(RMD)){
          graphics.off()
        }else{
          words <- 'Observed (green), predicted (black), shaded forecast (if modelYears)'
          message( words )
          caption <- c(caption, words)
        }
      }else{
        dev.off()
      }
    }
    
    yfun    <- colorRampPalette( c('tan', 'brown','turquoise','steelblue') )
    yearCol <- yfun(nyr)
    names(yearCol) <- tyears
  }
  
  ########## tree correlations over years
  
  if(length(years) > 3){
    
    if(is.null(RMD)) graphics.off()
    
    file <- paste('treeCor.pdf',sep='')
    
    if(SAVEPLOTS)pdf( file=.outFile(outFolder,file) )
    
    breaks <- seq(-1.1,1.1,by=.1)
    ylim <- c(0,5)
    nplot <- length(plots)
    
    ppt <- character(0)
    
    for(j in 1:nplot){
      
      wjk <- tdata$tnum[ tdata$plot == plots[j] & 
                           prediction$fecPred$matrPred > .5]
      njk <- length(unique(wjk))
      if(njk < 2)next
      
      if(njk > 1000){
        wjk <- njk[sample(1000)]
        njk <- length(wjk)
      }
      ojk <- omegaE[wjk,wjk]
      ojk[is.na(ojk)] <- 0
      if(max(ojk) == 0)next
      
      ppt <- c(ppt,plots[j])
    }
    
    npp <- length(ppt)
    
    mfrow <- .getPlotLayout(npp)
    par(mfrow=mfrow$mfrow, mar=c(1,1,1,2), oma=c(3,3,1,1), bty='n')
    
    for(j in 1:npp){
      
      jk <- 0
      sk <- character(0)
      ek <- numeric(0)
      
      for(k in 1:nspec){
        
        wjk <- tdata$tnum[ tdata$species == specNames[k] &
                             tdata$plot == ppt[j] ]
        njk <- length(unique(wjk))
        if(njk < 2)next
        
        wjk <- sort(unique(wjk))
        ojk <- omegaE[wjk,wjk]
        
        ojk[is.na(ojk)] <- 0
        oj <- ojk
        diag(oj) <- 0
        rs <- which( rowSums(oj) == 0 )
        diag(oj) <- diag(ojk)
        if(length(rs) > 0)oj <- oj[-rs,-rs]
        
        if(length(oj) < 2)next
        
        oj[oj > .95] <- .95
        oj[oj < -.95] <- -.95
        
        diag(oj) <- 1
        
        jk <- jk + 1
        sk <- c(sk,specNames[k])
        
        oj <- oj[lower.tri(ojk)]
        ovec <- hist(oj, breaks = breaks, plot = FALSE)$density
        
        tmp <- .getPoly(breaks[-1],ovec)
        if(jk == 1){
          plot(tmp[1,], tmp[2,],type='s',lwd=2,
               col=.getColor(specCol[specNames[k]],.3),
               xlab='', ylab='', ylim=ylim, xaxt='n', yaxt='n')
          axis(1, at=c(-1,0,1), labels=c(-1,0,1))
          axis(2, labels  = TRUE)
        }
        
        tmp <- .getPoly(breaks[-1],ovec)
        polygon( tmp[1,], tmp[2,], col=.getColor(specCol[specNames[k]],.3), lwd=2, 
                 border=specCol[specNames[k]])
      }
      if(length(sk) == 0)next
      .plotLabel(ppt[j],'topleft')
      legend('topright',sk,text.col=specCol[sk],bty='n')
    }
    
    if(jk > 0){
      mtext('Correlation', side=1, outer  = TRUE, line=1)
      mtext('Density', side=2, outer  = TRUE, line=1)
    }
    
    if(CONSOLE)
      readline('tree correlation in time -- return to continue ')
    if(SAVEPLOTS)dev.off( )
    if(is.null(RMD)){
      graphics.off()
    }else{
      words <- 'Correlations between trees, over time'
      message( words )
      caption <- c(caption, words)
    }
  }
  
  ############# predicted maps
  
  if(!SEEDDATA)return( invisible(caption) )
  
  if(MAPS){
    
    treeSymbol <- output$prediction$fecPred$fecEstMu
    if(is.null(treeSymbol))treeSymbol <- treeData$diam
    
    if(is.null(RMD)) graphics.off()
    
    mpp <- 1:length(plots)
    
    plotDims <- as.matrix(plotDims)
    
    for(m in mpp){
      
      xlim <- plotDims[plots[m],c('xmin','xmax')]
      if(is.na(xlim[1]))next
      
      ylim <- plotDims[plots[m],c('ymin','ymax')]
      
      dx <- diff(xlim)
      dy <- diff(ylim)
      ratio <- dx/dy
      
      pyr <- plotYrTable[plots[m],]
      pyr <- as.numeric( colnames(plotYrTable)[pyr > 0] )
      
      mfrow <- c(2, 2)
      if( max(c(dx,dy)) > 100){
        if(ratio > 2) mfrow <- c(2,1)
        if(ratio < .5)mfrow <- c(1,2)
      }
      
      if(!is.null(RMD))mfrow <- c(1,1)
      
      nperPage <- prod(mfrow)
      
      yrm <- years[years %in% pyr]
      ny  <- length(yrm)
      
      k   <- 0
      add <- FALSE
      o   <- 1:nperPage
      o   <- o[o <= nyr]
      
      seedMax <- as.matrix( sdata[sdata$plot == plots[m] & sdata$year %in% yrm ,seedNames] )
      if(is.matrix(seedMax)){
        if(nrow(seedMax) == 0){
          seedMax <- rbind(seedMax,0)
        }
        seedMax <- rowSums(seedMax, na.rm  = TRUE)
      }
      
      seedMax <- max(seedMax, na.rm  = TRUE) + 1
      
      while( max(o) <=  nyr & length(o) > 0 ){
        
        if(length(o) < 5 & is.null(RMD))mfrow <- c(2,2)
        
        yr <- yrm[o]
        yr <- yr[is.finite(yr)]
        if(length(yr) == 0)break
        
        if(is.null(RMD)) graphics.off()
        
        file <- paste('map_',plots[m],'_',yr[1],'.pdf',sep='')
        
        if(SAVEPLOTS)pdf( file=.outFile(outFolder,file) )
        
        mapList <- output
        mapList$treeSymbol <- treeSymbol
        mapList$mapPlot <- plots[m]
        mapList$xlim <- xlim
        mapList$ylim <- ylim
        mapList$PREDICT <- TRUE
        mapList$mapYears <- yr
        mapList$treeScale <- .9
        mapList$trapScale <- 1
        mapList$mfrow <- mfrow
        mapList$seedMax <- seedMax
        mapList$plotScale <- 1
        mapList$COLORSCALE <- FALSE
        mapList$LEGEND <- TRUE
        mapList$RMD <- RMD
        
        add <- mastMap(mapList)
        
        if(add)scaleBar('m', value = 20, yadj=.07)
        
        if(CONSOLE)
          readline('predicted fecundity, seed data maps -- return to continue ')
        if(SAVEPLOTS)dev.off( )
        if(is.null(RMD)){
          graphics.off()
        }else{
          words <- paste('Prediction map for plot', plots[m],
                         'in year', yr)
          message( words )
          caption <- c(caption, words)
        }
        
        o <- o + nperPage
        o <- o[o <= nyr]
        
        if(length(o) == 0){
          break
        }
        
        if(!add)next
      }  
    }
  }
  
  ################# spatio-temporal correlation
  
  if(SPACETIME){
    
    # trees/sites ordered by similarity at zero lag
    
    mvs <- suppressWarnings(
      meanVarianceScore(output, Q = pnorm(c(0, -1, 1)), nsim=30, LAGMAT  = TRUE,
                        cyr = 8, ktree = 20, maxArea = 30^2, CLOSE  = TRUE)
    )
    
    treeCov <- mvs$lagCanopy
    trapCov <- mvs$lagGround
    nkk <- length(treeCov)
    plotk <- names(treeCov)
    gridArea <- mvs$gridArea
    
    if(is.null(RMD)) graphics.off()
    
    if(nkk > 0){
      
      col2 <- colorRampPalette(c("#67001F", "#67001F", "#B2182B", "#D6604D", "#F4A582",
                                 "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                                 "#4393C3", "#2166AC", "#053061", "#053061"))
      for(k in 1:nkk){
        
        wpt <- which(names(treeCov) == plotk[k])
        wpc <- which(names(trapCov) == plotk[k])
        
        if(length(wpt) == 0 | length(wpc) == 0)next
        
        tvar <- treeCov[[wpt]]
        cvar <- trapCov[[wpc]]
        
        if(length(tvar) < 2 | length(cvar) < 2) next
        
        tmp  <- columnSplit(colnames(tvar),'_')
        
        klag <- as.numeric(tmp[,ncol(tmp)])
       
        if(nrow(tvar) < 2)next
        
        km <- max(klag)
        if(km > 5)km <- 5
        
        file <- paste('spaceTime',plotk[k],'.pdf',sep='')
        if(SAVEPLOTS)pdf( file=.outFile(outFolder,file) )
        
        par(mfrow=c(2,km+1), mar=c(2,1,1,1), oma=c(2,1,2,1), bty='n',xpd  = TRUE)
        
        order  <- 'hclust'
        
        for(i in 0:km){
          
          wl <- which(klag == i)
          stree <- tvar[,wl]
          rnames <- rownames(stree)
          
          if(i > 0){
            order <- 'original'
            colnames(stree) <- rownames(stree)
            stree <- stree[rnames,rnames]
          }
          
          tmp <- corrplot(stree, is.corr  = TRUE, method='color', col=rev(col2(200)),
                          tl.pos='n', cl.length=3, cl.lim=c(-1,1), type='lower',
                          order=order, cl.pos='n')
          rnames <- rownames(tmp)
          tlab <- paste('lag',i)
          title(main = list(tlab, cex = 1.5,
                            font = 1), line = -2, adj=1)
          if(i == 0)title(main = list("Canopy", cex = 1.5,
                                      font = 3))
        }
        
        tmp  <- columnSplit(colnames(cvar),'_')
        klag <- as.numeric(tmp[,ncol(tmp)])
        cvar[cvar < -1] <- 0
        cvar[cvar > 1] <- 0
        
        order  <- 'hclust'
        
        for(i in 0:km){
          wl <- which(klag == i)
          if(length(wl) == 0)next
          stree <- cvar[,wl]
          
          if(ncol(stree) < nrow(stree)){
            sc <- columnSplit(colnames(stree),'_')[,1]
            w0 <- which(!rownames(stree) %in% sc)
            c1 <- c(1:(w0-1))
            c2 <- c((w0+1):ncol(stree))
            c1 <- c1[c1 > 0]
            c2 <- c2[c2 < nrow(stree) & c2 > w0]
            stree <- cbind( stree[,c1], 0, stree[,c2] )
            colnames(stree)[w0] <- paste(rownames(stree)[w0],i,sep='_')
          }
          
          if(i > 0){
            order <- 'original'
            colnames(stree) <- rownames(stree)
            stree <- stree[rnames,rnames]
          }
          
          tmp <- corrplot(stree, is.corr  = TRUE, method='color', col=rev(col2(200)),
                          tl.pos='n', cl.length=3, cl.lim=c(-1,1), type='lower',
                          cl.pos='n')
          rnames <- rownames(tmp)
          if(i == 0)title(main = list("Forest floor", cex = 1.5,
                                      font = 3))
        }
        mtext(plotk[k], 1, outer  = TRUE, line=0)
        
        
        if(CONSOLE)
          readline('tree-time (above), space-time (below) -- return to continue ')
        if(SAVEPLOTS)dev.off()
        if(is.null(RMD)){
          graphics.off()
        }else{
          words <- 'Tree-time (above), space-time (below)'
          message( words )
          caption <- c(caption, words)
        }
      }
    }
    
    ############## score by scale
    
   # darea <- 100
    
    mvs <- suppressWarnings(
      meanVarianceScore(output, Q = pnorm(c(0, -1, 1)), ktree = 20, 
                        nsim=20, LAGMAT  = TRUE, CLOSE  = TRUE)
    )
    gridArea <- mvs$gridArea
    
    if(is.null(RMD)) graphics.off()
    
    scoreT <- scoreS <- scoreTse <- scoreSse <- numeric(0)
    pname  <- character(0)
    
    for(k in 1:length(plots)){
      
      wt <- which(names(mvs$scoreTree) == plots[k])
      ws <- which(names(mvs$scoreSeed) == plots[k])
      
      if(length(wt) == 0 | length(ws) == 0)next
      
      dtree <- mvs$scoreTree[[wt]]
      dseed <- mvs$scoreSeed[[ws]]
      dtreeSE <- mvs$scoreTreeSe[[wt]]
      dseedSE <- mvs$scoreSeedSe[[ws]]
      
      if(length(dtree) > 2 & length(dseed) > 2){
        scoreT <- append(scoreT, list(dtree))
        scoreS <- append(scoreS, list(dseed))
        scoreTse <- append(scoreTse, list(dtree))
        scoreSse <- append(scoreSse, list(dseed))
        
        pname  <- c(pname,plots[k])
      }
    }
    
    names(scoreT) <- names(scoreS) <- names(scoreTse) <- 
      names(scoreSse) <- pname
  
    
    pk <- names(scoreT)[k]
    dss <- scoreT
    ylab  <- 'Number of hosts'
    title <- 'Canopy'
    file  <- 'resourceScoreCanopy.pdf'
    carea <- 1
    q <- seq(0, 1, length=15)^1
    cols <- .getColor('black',q)
    
    npp <- length(scoreT)
    
    if(npp > 0){
      
      for(j in 1:2){
        
        if(j == 2){
          dss <- scoreS
          file <- 'resourceScoreGround.pdf'
          yy <- as.numeric(rownames(dk))
          ylab <- 'Distance (m)'
          title <- 'Forest floor'
        }
        
        zscale <- range(sapply(dss, range, na.rm  = TRUE))
        cseq   <- seq(zscale[1], zscale[2], length=30)
        
        xscale <- max(sapply(dss, ncol))
        yscale <- max(sapply(dss, nrow))
        if(yscale < 100)yscale <- 100
        if(j == 2 & yscale < 20)yscale <- 20
        
        xlim <- log(1 + c(0, xscale))
        ylim <- log(1 + c(0, yscale+10))
        
        if(SAVEPLOTS)pdf( file=.outFile(outFolder,file) )
        
        mff <- .getPlotLayout(length(dss))
        
        par(mfrow=mff$mfrow, bty='l',mar=c(1,1,.1,.1), oma=c(3,3,1,1))
        
        for(k in 1:npp){
          
          dk <- dss[[k]]
          xx <- columnSplit(colnames(dk),'_')[,2]
          xx <- as.numeric(xx)
          yy <- c(1:nrow(dk))
          
          if(j == 2){
            aa <- as.numeric(rownames(dk))
            yy <- (pi*aa^2)/10000
          }
          
          levels <- quantile(dk, q )
          levels <- sort(unique(levels))
          
          ytick <- c(1, 5, c(1:5)*10 )
          lx <- log(xx + 1)
          ly <- log(yy + 1)
          ltick <- log(ytick + 1) 
          
          if(j == 1)ylim[1] <- ly[1]/2
          if(j == 2)ylim[1] <- diff(ly[1:2])/2
          
          plot(NA, axes = F, xlim=xlim, ylim=ylim, xlab='', ylab='')
          contour(lx, ly, t(dk), levels=levels, col=cols, labels='', add  = TRUE,
                  axes = FALSE)
          .filled.contour(lx, ly, t(dk), levels=levels, col=cols)
          
          if(length(yy) > 10){
            #     yy <- c(0,yy)
            #     ly <- c(0,ly)
            bb <- ceiling(length(yy)/10)
            ss <- seq(1, length(yy), by=bb)
            ytick <- ytick[ss]
            ltick <- ltick[ss]
          }
          
          xlabs <- ylabs <- FALSE
          
          if(k %in% mff$bottom)xlabs <- xx
          
          if(j == 2){
            ytick <- c(100, 1000, 10000, 50000, 100000)/10000
            ltick <- log(ytick + 1)
          }
          
          if(k %in% mff$left){
            ylabs <- ytick
            if(j == 2)ylabs <- c('100 m2', '1000 m2', '1 ha', '5 ha','10 ha')
          }
          
          axis(1, at=lx, labels=xlabs)
          axis(2, at=ltick, labels=ylabs)
          
          .plotLabel(names(dss)[k],'topright')
        }
        mtext('Years', 1, outer  = TRUE, line=1)
        mtext(ylab, 2, outer  = TRUE, line=1)
        
        endLabs <- signif(range(dk),1)
        
        clist <- list( kplot=1, ytick=NULL, text.col = 'black',
                       cols=cols, labside='right', text.col=col,
                       bg='grey', endLabels=endLabs) 
        cornerScale( clist )
        
        if(CONSOLE)
          readline('score by scale -- return to continue ')
        if(SAVEPLOTS)dev.off()
        if(is.null(RMD)){
          graphics.off()
        }else{
          words <- 'Score by scale'
          message( words )
          caption <- c(caption, words)
        }
      }
      
      # plot comparison
      
      nhost <- 5
      nhy   <- 2
      tmat  <- matrix(NA, npp, 3)
      colnames(tmat) <- c('mu','lo','hi')
      rownames(tmat) <- names(scoreT)
      smat <- tmat
      
      for(k in 1:npp){
        
        nn  <- nhost
        dkk <- scoreT[[k]]
        if(nn > nrow(dkk))nn <- nrow(dkk)
        skk <- scoreTse[[k]]
        mu  <- dkk[nn,nhy-1]
        ss  <- skk[nn,nhy-1]
        tmat[k,] <- c(mu, mu + ss*c(-1,1))
        
        nn  <- nhost
        dkk <- scoreS[[k]]
        if(nn > nrow(dkk))nn <- nrow(dkk)
        skk <- scoreSse[[k]]
        mu  <- dkk[nn,nhy-1]
        ss  <- skk[nn,nhy-1]
        smat[k,] <- c(mu, mu + ss*c(-1,1))
      }
      
      file <- 'scoreByPlot.pdf'
      if(SAVEPLOTS)pdf( file=.outFile(outFolder,file) )
      
      par(mfrow=c(1,1), mar=c(4,4,2,2), bty='n',xpd = FALSE)
      xlim <- range(tmat) + c(-.1,.2)
      ylim <- range(smat) + c(-.1,.2)
      
      ylab <- paste('Forest floor, ', nhost*gridArea,'m2')
      xlab <- paste('Canopy score,',nhost,'host trees')
      
      plot(NA, xlim=xlim, ylim=ylim, xlab=xlab,
           ylab = ylab)
      points(tmat[,1],smat[,1])
      segments(tmat[,1], smat[,2], tmat[,1], smat[,3])
      segments(tmat[,2], smat[,1], tmat[,3], smat[,1])
      text(tmat[,1]+.1, smat[,1]+.1,names(mvs$scoreTree), pos=4)
      abline(0,1,col='grey',lwd=2,lty=2)
      .plotLabel(paste(nhy,' yr'), 'bottomright')
      
      if(CONSOLE)
        readline('score by plot -- return to continue ')
      if(SAVEPLOTS)dev.off()
      if(is.null(RMD)){
        graphics.off()
      }else{
        words <- 'Score by plot'
        message( words )
        caption <- c(caption, words)
      }
    }
    
    if(is.null(RMD)) graphics.off()
    
    file <- 'scoreComponents.pdf'
    if(SAVEPLOTS)pdf( file=.outFile(outFolder,file) )
    
    score4plot <- function(totalScore, cname){
      
      mcol <- paste(cname,'_mu',sep='')
      vcol <- paste(cname,'_stdDev',sep='')
      scol <- paste(cname,'_score', sep='')
      
      mrow <- grep('_mu',rownames(totalScore))
      cmu <- totalScore[mrow,mcol]
      chi <- sqrt(cmu + totalScore[mrow+1,mcol])
      clo <- sqrt(cmu - totalScore[mrow+1,mcol])
      
      cse <- totalScore[mrow+1,mcol]
      cse[!is.finite(cse)] <- cmu[!is.finite(cse)]
      cup <- cmu + cse
      
      vmu <- totalScore[mrow,vcol]
      vse <- totalScore[mrow,vcol]
      vse[!is.finite(vse)] <- vmu[!is.finite(vse)]
      vup <- vmu + vse
      
      scoreMu <- totalScore[mrow, scol]
      scoreSe <- totalScore[mrow+1, scol]
      
      vmat <- rbind(cmu, vmu)
      ylim <- range(c(cmu, vmu, cup, vup))
      list(vmat = vmat, mup = cup, vup = vup, ylim = ylim,
           score = cbind(scoreMu, scoreSe))
    }
    
    par(mfcol=c(2,2), mar=c(1,2,1,4), oma=c(1,2,1,1)) # canopy mean vs variance
    cols <- c("darkgreen", "brown")
    
    totalScore <- mvs$totalScore
    
    ##################
    
    splot      <- columnSplit(rownames(totalScore))[,1]
    spp        <- splot[!duplicated(splot)]
    
    ttt  <- score4plot(totalScore, 'canopy')
    vmat <- ttt$vmat
    
    ylim[2] <- (ttt$ylim[2] - ttt$ylim[1])*10
    ylim[1] <- ttt$ylim[1]*.5
    
    xl <- barplot( vmat, width=.5, beside  = TRUE, ylim=ylim,
                   col = cols, xaxt='n', log='y')
    smu <- as.vector( vmat )
    upm <- rbind(ttt$mup, ttt$vup)
    sup <- as.vector( upm )
    stt <- apply(upm, 2, max, na.rm  = TRUE)
    
    segments(xl, smu, xl, sup)
    segments(xl-.2, sup, xl+.2, sup)
    mtext('n host trees',2,line=3) 
    text(colMeans(xl)-.5, 1.3*stt, spp, srt=90, pos=4, cex=.85)
      
    title('Canopy')
    
    score <- ttt$score
    sup   <- score[,1] + score[,2]
    sdn   <- score[,1] - score[,2]
    
    ylim <- range(c(sdn, sup))
    
    wn <- which(score[,1] < 0)
    cc <- rep(cols[1],length(sup))
    cc[wn] <- cols[2]
    
    xl <- barplot( score[,1], width=.5, beside  = TRUE, ylim=ylim,
                   col = cc, xaxt='n')
    mtext('Score',2,line=3) 
    sb <- sup
    sb[wn] <- sdn[wn]

    segments(xl, score[,1], xl, sb)
    segments(xl-.1, sb, xl+.1, sb)
    
    ttt  <- score4plot(totalScore, 'ground')
    vmat <- ttt$vmat
    
    ylim[2] <- (ttt$ylim[2] - ttt$ylim[1])*10
    ylim[1] <- ttt$ylim[1]*.5
    
    xl <- barplot( vmat, width=.5, beside  = TRUE, ylim=ylim,
                   col = cols, xaxt='n', log='y')
    mtext('area',2,line=3) 
    smu <- as.vector( vmat )
    sup <- as.vector( rbind(ttt$mup, ttt$vup) )
    segments(xl, smu, xl, sup)
    segments(xl-.2, sup, xl+.2, sup)
    
    title('Forest floor')
    
    legend('topright', legend=c('mean benefit','variance cost'), 
           text.col=cols, bty='n')
    
    score <- ttt$score
    sup   <- score[,1] + score[,2]
    sdn   <- score[,1] - score[,2]
    
    ylim <- range(c(sdn, sup))
    
    wn <- which(score[,1] < 0)
    cc <- rep(cols[1],length(sup))
    cc[wn] <- cols[2]
    
    xl <- barplot( score[,1], width=.5, beside  = TRUE, ylim=ylim,
                   col = cc, xaxt='n')
    sb <- sup
    sb[wn] <- sdn[wn]
    
    segments(xl, score[,1], xl, sb)
    segments(xl-.2, sb, xl+.2, sb)
    
    if(CONSOLE)
      readline('score components -- return to continue ')
    if(SAVEPLOTS)dev.off()
    if(is.null(RMD)){
      graphics.off()
    }else{
      words <- 'Score by plot'
      message( words )
      caption <- c(caption, words)
    }
    
    ################ entropy
    
    entropy <- mvs$entropy
    
    if(!is.null(entropy)){
      
      if(nrow(entropy) > 4){
        
        if(is.null(RMD)) graphics.off()
        
        file <- 'entropy.pdf'
        if(SAVEPLOTS)pdf( file=.outFile(outFolder,file) )
        
        par(mfrow=c(1,2), mar=c(4,4,1,.5), oma=c(1,1,1,1), bty='n',xpd = FALSE)
        
        entropy[!is.finite(entropy)] <- NA
        
        xl <- range(entropy[,1], na.rm  = TRUE )
        dx <- .3*diff(xl)
        xl[2] <- xl[2] + dx
        yl <- range(entropy[,1])
        
        we <- grep('tree-tree',rownames(entropy))
        wr <- grep('site-site',rownames(entropy))
        rnames <- unlist( strsplit(rownames(entropy)[we],'_tree-tree') )
        
        xl <- range(entropy[we,1], na.rm  = TRUE ) + c(-1,1)
        dx <- .3*diff(xl)
        xl[2] <- xl[2] + dx
        yl <- range(entropy[wr,1], na.rm  = TRUE) + c(-1,1)
         
        plot(entropy[we,1],entropy[wr,1], xlim=xl, ylim=yl, xlab='', 
             ylab='Forest floor entropy', cex=.01, yaxt='n')
        axis(2, line=1)
        abline(0,1,lty=2)
        
        par(new = FALSE,xpd  = TRUE)
        text(entropy[we,1],entropy[wr,1],rnames)
        mtext('Canopy entropy',1,outer  = TRUE,line=-1)
        .plotLabel('a) Spatial', above  = TRUE)
        
        we <- grep('tree-lag',rownames(entropy))
        wr <- grep('site-lag',rownames(entropy))
        rnames <- unlist( strsplit(rownames(entropy)[we],'_tree-lag') )
        
        xl <- range(entropy[we,1], na.rm  = TRUE ) + c(-1,1)
        dx <- .3*diff(xl)
        xl[2] <- xl[2] + dx
        yl <- range(entropy[wr,1], na.rm  = TRUE) + c(-1,1)
        
        plot(entropy[we,1],entropy[wr,1], xlim=xl, ylim=yl, xlab='', ylab='', 
             cex=.01, yaxt='n')
        axis(2, line=1)
        # abline(0,1,lty=2)
        
        par(new = FALSE,xpd  = TRUE)
        text(entropy[we,1],entropy[wr,1],rnames)
        .plotLabel('b) Temporal', above  = TRUE)
        
        par(new  = TRUE,xpd = FALSE)
        
        if(CONSOLE)
          readline('entropy -- return to continue ')
        if(SAVEPLOTS)dev.off()
        if(is.null(RMD)){
          graphics.off()
        }else{
          words <- 'Entropy'
          message( words )
          caption <- c(caption, words)
        }
      }
    }
  }
   
  invisible( list(caption = caption, diam90 = diam90) )
}
  
.getPoly <- function(x,y){
  
  dx <- diff(x)
  xx <- c(x[1] - dx[1]/2, x[-1] - dx/2)
  xx <- rep(xx,each=2)
  yy <- rep(y,each=2)
  yy <- c(0,yy,0)
  xx <- c(xx,xx[length(xx)],xx[length(xx)])
  rbind(xx,yy)
}

.chainPlot <- function(mat, burnin, label, ngLab = NULL, burnLab = NULL,
                       refVals = NULL, CONSOLE, RMD,
                       SAVEPLOTS = FALSE, outFolder='', ALLONE = F, 
                       cols = NULL, ylim = NULL, intval=NULL){
  
  words <- character(0)
  
  if(is.null(ngLab))ngLab <- ng
  if(is.null(burnLab))burnLab <- burnin
  
  if(!is.null(refVals)){
    if(length(refVals) == 1 & ncol(mat) > 1)refVals <- rep(refVals, ncol(mat))
  }
  
  cseq <- 1:nrow(mat)
  if(length(cseq) > 2000)cseq <- round(seq(1,length(cseq),length=1000))
  
  if(SAVEPLOTS){
    fileName <- .replaceString(label,',','')
    fileName <- .replaceString(label,' ','')
    fileName <- paste(fileName,'.pdf',sep='')
    pdf( file=.outFile(outFolder,fileName) )
  }
  
  cnames <- .coeffNames( colnames(mat) )
  colnames(mat) <- cnames
  
  npp <- length(cnames)
  if(npp > 25)npp <- 25
  
  if(ALLONE)npp <- 1
  
  mfrow <- .getPlotLayout(npp)
  par(mfrow=mfrow$mfrow, bty='n', mar=c(2,2,2,2), oma=c(2,3,1,1)) 
  
  cex <- 1/(1 + mfrow$mfrow[2])^.1
  
  ng <- nrow(mat)
  
  cseq <- 1:ng
  burnline <- burnin
  ss   <- burnin:ng
  if(nrow(mat) > 1000){
    cseq <- seq(1,ng,length=1000)
    burnline <- burnin/ng*1000
    ss <- cseq[ cseq > burnin ]
  }
  
  if(is.null(ylim) & ALLONE)ylim <- range(mat)
  
  NEWY <- FALSE
  if(is.null(ylim))NEWY <- TRUE
  
  naa <- 0
  
  if(is.null(cols))cols <- rep('black', ncol(mat))
  
  xmax <- max(cseq)
  
  for(j in 1:ncol(mat)){
    
    if(j %in% c(25, 50)){
      naa <- naa + 1
      if(CONSOLE){
        lab <- paste(label, ' -- return to continue')
        readline(lab)
      }
      if(SAVEPLOTS) {
        dev.off( )
      }
      if(!is.null(RMD)){
        words <- paste( 'MCMC chains for', label, 'with 95% coverage' )
        message( words )
      }
      if(SAVEPLOTS){
        fileName <- .replaceString(label,',','')
        fileName <- .replaceString(label,' ','')
        fileName <- paste(fileName,'_',letters[naa],'.pdf',sep='')
        pdf( file=.outFile(outFolder,fileName) )
      }
      npp <- ncol(mat) - j
      if(npp > 25)npp <- 25
      mfrow <- .getPlotLayout(npp)
      par(mfrow=mfrow$mfrow, bty='n', mar=c(2,2,2,2), oma=c(2,3,1,1)) 
    }
    
    xlabels <- FALSE
    
    if(j %in% mfrow$bottom)xlabels <- TRUE
    
    mj <- mat[,j]
    
    if(NEWY){
      ylim <- range(mj)
      if(!is.null(refVals) & NEWY){
        ylim <- range(c(refVals[j], mj))
        expd <- diff(ylim)
        ylim <- ylim + c(-1, 1)*.5*expd
      }
    }
    if(cnames[j] %in% c('sigma','mspe'))ylim <- c(0, 1.2*ylim[2])
      
    
    if( j == 1 | !ALLONE ){
      plot(mj[cseq], type='l', ylim = ylim, xaxt='n', xlab='', ylab='',
           col=cols[j])
      if(xlabels){
        axis(1, at = c(0, burnline, 1000), labels = c(0, burnLab, ngLab))
      }else{
        axis(1, at = c(0, burnline, 1000), labels = F)
      }
    }else{
      lines(mj[cseq], col=cols[j])
    }
    if(!is.null(intval)){
      int <- intval[colnames(mat)[j],]
      if(length(int) > 0){
        if('min' %in% colnames(intval)){
          lines(c(0,0),intval[j,c('min','max')], col=cols[j])
        }
        if('mean' %in% colnames(intval) &
           'var' %in% colnames(intval)){
          dseq <- seq(ylim[1],ylim[2],length=100)
          pline <- dnorm(dseq, intval[j,'mean'], sqrt(intval[j,'var']))
          xx <- c(xmax*.8*pline, dseq*0)
          yy <- c(dseq, rev(dseq))
          polygon(xx, yy, col=.getColor(cols[j], .2))
          lines(xx, yy, col=cols[j])
        }
      }
    }
    
    q <- quantile( mj[ss], c(.025, .5, .975) )
    for(k in 1:3){
      segments(burnline, q[k], 1000, q[k], col='white',lwd=1)
      segments(burnline, q[k], 1000, q[k], col=cols[j], lty=2)
    }
    segments(burnline, q[1], burnline, q[3], col='white',lwd=1)
    segments(burnline, q[1], burnline, q[3], col=cols[j], lty=2)
    if(!is.null(refVals))abline(h = refVals[j], col='blue')
    if(!ALLONE){
      .plotLabel(cnames[j], above  = TRUE, cex=cex)
    }
  }
  mtext('Iteration', outer  = TRUE, side=1, line=1)
  mtext('Parameter value', outer  = TRUE, side=2, line=1)
  if(ALLONE){
    legend('topright', colnames(mat), text.col=cols, bty='n')
  }
 
  if(CONSOLE){
    lab <- paste(label, ' -- return to continue')
    readline(lab)
  }
  if(SAVEPLOTS){
    dev.off( )
  }
  if(!is.null(RMD)){
    words <- paste( 'MCMC chains for', label, 'with 95% coverage' )
    message( words )
  }
  invisible(words)
}
 
.outFile <- function(outFolder=character(0),file){
  paste(outFolder,file,sep='/')
}

.plotLabel <- function(label,location='topleft',cex=1.3,font=1,
                       above = FALSE, below = FALSE, bg=NULL, wrap = 1000){
  
  # wrap - no. of characters to wrap label to two lines
  
  if(above){
    adj <- 0
    if(location == 'topright')adj=1
    title(label,adj=adj, font.main = font, font.lab =font, cex.main=cex)
    return()
  }
  if(below){
    adj <- 0
    if(location == 'bottomright')adj=1
    mtext(label,side=1,adj=adj, outer = FALSE,font.main = font, font.lab =font,cex=cex)
    return()
  }
  
  if(is.null(bg)){
    tmp <- legend(location,legend=' ',bty='n')
  } else {
    tmp <- legend(location,legend=label,bg=bg,border=bg,text.col=bg,bty='o')
  }
  
  xt <- tmp$rect$left # + tmp$rect$w
  yt <- tmp$text$y
  
  pos <- 4
  tmp <- grep('right',location)
  if(length(tmp) > 0)pos <- 2
  
  XX <- par()$xlog
  YY <- par()$ylog
  
  if(XX)xt <- 10^xt
  if(YY)yt <- 10^yt
  
  if( wrap > nchar(label) ){
    text(xt,yt,label,cex=cex,font=font,pos=pos)
    return()
  }
  
  # split at spaces
  words <- unlist( strsplit(label, ' ') )
  if(length(words) == 1){
    text(xt,yt,label,cex=cex,font=font,pos=pos)
    return()
  }
  
  nw <- cumsum( nchar(words) + 1 )
  br <- which(nw <= wrap)
  if(length(br) == 0){
    br <- 1
    wrap <- nchar(words[1]) + 1
  }
  n2 <- which(nw > wrap)
  lab1 <- words[br]
  lab2 <- words[n2]
  if(length(br) > 1) lab1 <- paste0(words[br], collapse = ' ')
  if(length(n2) > 1) lab2 <- paste0(words[n2], collapse = ' ')
  
  yaxp <- par('yaxp')
  da <- diff(yaxp[1:2])
  yz <- yt + .2*c(-da,da)
  if(YY){
    da <- diff(log(yaxp[1:2]))
    yz <- exp( log(yt) + .15*c(-da,0) )
  }
  
  text(xt,yz[2],lab1,cex=cex,font=font,pos=pos)
  text(xt,yz[1],lab2,cex=cex,font=font,pos=pos)

}

commas4numbers <- function(x){
  
  x  <- as.character( round(x) )
  lx <- l1 <- nchar(x)
  nr <- floor( lx/3 ) + 1
  
  l0 <- lx - 2
  xn <- character(0)
  for(k in 1:nr){
    xk <- substr(x, l0, l1 )
    xn <- paste(',', xk,xn, sep='')
    l0 <- l0 - 3
    l1 <- l1 - 3
    if(l0 < 1)l0 <- 1
  }
  if( startsWith(xn, ',') ) xn <- substr(xn, 2, 1000)
  if( startsWith(xn, ',') ) xn <- substr(xn, 2, 1000)
  if( startsWith(xn, ',') ) xn <- substr(xn, 2, 1000)
  
  xn
}
.boxCoeffs <- function(chain, snames, xlab = "", ylab='Coefficient',
                       addSpec = 'species', ylim=NULL, cols = NULL,
                       xaxt = 's', yaxt = 's'){
  
  nspec  <- length(snames)
  cnames <- colnames(chain)
  xn     <- character(0)
  vnames <- numeric(0)
  iname  <- character(0)
  gnames <- paste(addSpec,snames,sep='')

  for(j in 1:nspec){
    ij     <- which(cnames == gnames[j])
    if(length(ij) > 0)iname <- 'intercept'
    wj     <- grep(gnames[j],cnames)
    if(length(wj) > 0)vnames <- rbind(vnames, wj)
    wk <- grep(':',cnames[wj])
    if(length(wk) > 0){
      xn <- matrix( unlist( strsplit(cnames[wj[wk]],':')),
                    ncol=2, byrow  = TRUE)[,2]
    }
  }
    
  rownames(vnames) <- snames[vnames[,1]]
  colnames(vnames) <- c(iname,xn)
  nv <- ncol(vnames)
  
  nss <- nrow(vnames)
  
  atvals <- c(1:nss)/(nss + 1)
  atvals <- atvals - mean(atvals)
  sseq   <- c(1:nv)
  xlim   <- c(.5, nv +.5)
  
  if(is.null(ylim)){
    ylim <- range(chain)
    ylim[1] <- ylim[1] - diff(ylim)*.25
  }
  
  add <- FALSE
  if(is.null(cols))cols <- seq(1:nss)
  
  xlabel <- ''
  
  stats <- xtick <- numeric(0)
  
  for(j in 1:nv){
    
    jcol <- vnames[,j]
    if(j > 1)add <- TRUE
    chainj <- chain[,jcol, drop = FALSE]
    mj     <- mean(chainj)
    sj     <- sd(chainj)
    chainj <- chainj #/sj
    
    if(j == nv)xlabel <- xlab
    
    colj <- cols[j]
    
    if(any(colnames(chainj) %in% snames))colj <- cols[colnames(chainj)]
    
    wi <- grep(':',colnames(chainj))
    if(length(wi) > 0){
      tt <- columnSplit(colnames(chainj)[wi],':')
      tk <- tt[,1]
      colj <- cols[tk]
      colnames(chainj) <- tt[,2]
    }
    
    boxPars <- .boxplotQuant( chainj, xaxt=xaxt, add = add,
                          at = atvals + j, xlim = xlim,
                          outline = FALSE, ylim=ylim,
                          col=.getColor(colj, .5), 
                          border=colj, lty=1,
                          ylab = ylab, yaxt = yaxt)
    stats <- cbind(stats,boxPars$stats)
    xtick <- rbind(xtick, atvals+j)
    
  }
  .plotLabel(xlab,'topleft',above  = TRUE)
 # legend('topright',snames, text.col=1:nspec, bty='n')
  abline(h = c(0), lwd=1, col=.getColor('grey',.6))
  boxPars$stats <- stats
  boxPars$xtick <- xtick
  
  
  invisible(boxPars)
}

.boxCoeffsMultiSpec <- function(chain, snames, xlab = "", ylab='Coefficient',
                       addSpec = 'species', ylim = NULL, cols = NULL,
                       xaxt = 's', yaxt = 's', cex = .85){
  
  nspec  <- length(snames)
  cnames <- colnames(chain)
  
  snn <- xnn <- cnames
  
  for(k in 1:length(snames)){
    sk <- paste(snames[k],':',sep='')
    wk <- grep( sk, cnames )
    if(length(wk) == 0)next
    xnn[wk] <- .replaceString(xnn[wk], sk, '')
    snn[wk] <- snames[k]
  }
  
  xnames   <- unique(xnn)
  allnames <- cbind(snn, xnn)
  
  if(is.null(ylim)){
    ylim <- range(chain)
    ylim[1] <- ylim[1] - diff(ylim)*.25
  }
  
  nv <- length(xnames)
  
  atvals <- c(1:nspec)/(nspec + 1)
  atvals <- atvals - mean(atvals)
  sseq   <- c(1:nv)
  xlim   <- c(.5, nv +.5)
  stats <- numeric(0)
  add   <- FALSE
  
  for(j in 1:nv){
    
    jcol <- which(allnames[,2] == xnames[j])
    scol <- which(snames %in% allnames[jcol,1])
    colj <- cols[scol]
    
    if(j > 1)add <- TRUE
    chainj <- chain[,jcol, drop = FALSE]
    at <- atvals + j
    
    boxPars <- .boxplotQuant( chainj, xaxt=xaxt, add = add,
                              at = at[scol], xlim = xlim,
                              outline = FALSE, ylim=ylim,
                              col=.getColor(colj, .5), 
                              border=colj, lty=1,
                              ylab = ylab, yaxt = yaxt)
    if(j == 1)abline(h=0, col='grey')
    
    if(nv == 1){
      .plotLabel(xnames[j], 'bottomleft', cex = cex)
      next
    }
    bstat <- boxPars$stats
    brange <- range(bstat)
    if(brange[1] > 0)brange[1] <- 0
    if(brange[2] < 0)brange[2] <- 0
    wside  <- which.max( abs(ylim - brange) )
    pos <- 2
    if(wside == 2)pos <- 4
    text(mean(at), brange[wside], xnames[j], pos=pos, srt = 90,
         cex = cex)
    if(j < nv)abline(v = (max(at) + diff(at[1:2])), lty=2, col='grey')
  }
  at
}

.coeffNames <- function(cvec){
  
  # clean names  for coefficients
  
  fnames  <- .replaceString(cvec, '(Intercept)', 'intercept' )
  fnames  <- .replaceString(fnames, 'I(', '' )
  fnames  <- .replaceString(fnames, '))', ')' )
  fnames  <- .replaceString(fnames, 'species','')
  fnames  <- .replaceString(fnames, '^2)','^2')
  fnames
}

.fixNamesVector <- function(vnames, data, MODE='keep'){
  
  wgg <- match(vnames,names(data))
  
  for(k in wgg){
    data[[k]] <- .fixNames(data[[k]], all  = TRUE, MODE)$fixed
  }
  data
}

.fixNames <- function(cvec, all = FALSE, MODE='keep'){
  
  # MODE == 'character', 'factor', or 'keep' (return same mode)
  
  cdup <- numeric(0)
  
  FACT <- FALSE
  if(is.factor(cvec)){
    FACT <- TRUE
    cvec <- as.character(cvec)
  }
  if(is.numeric(cvec)){
    cvec <- as.character(cvec)
    wdot <- grep('.',cvec)
    if(length(wdot) > 0)cvec <- .replaceString(cvec,'.','dot')
  }
  if(all) cvec <- .replaceString(cvec, '_','')
  cvec <- .replaceString(cvec, '-','')
  cvec <- .replaceString(cvec, ' ','')
  cvec <- .replaceString(cvec, "'","")
 # cvec <- .replaceString(cvec, ".","dot")
  cvec <- .replaceString(cvec, '"','')
  if( (FACT | MODE == 'factor') & MODE != 'character' ){
    cvec <- as.factor(cvec)
    droplevels(cvec)
  }
  
  cvec <- .replaceString(cvec, 'acerPenn', 'acerPens')
  
  wd <- which(duplicated(cvec))
  if(length(wd) > 0)cdup <- wd
    
  list(fixed = cvec, dups = cdup)
}

.setupR <- function(sdata, tdata, seedNames, specNames, verbose, unknown='UNKN'){
  
  SAMPR <- TRUE
  UCOLS <- FALSE
  
  wun <- grep('UNKN', seedNames)
  if(length(wun) > 1){
    stop("\ncan only have one seedNames with 'UNKN' class\n")
  }
  if(length(wun) == 1)UCOLS <- TRUE
    
  if(!'specPlot' %in% colnames(tdata))
    tdata$specPlot <- columnPaste(tdata$species,tdata$plot)
  
  if( length(seedNames) == 1 )SAMPR <- FALSE
  
  plots  <- sort(unique(as.character(tdata$plot)))
  priorR <- numeric(0)
  
  for(j in 1:length(plots)){
    
    wj   <- which(tdata$plot == plots[j])
    jtab <- table(tdata$species[wj])
    jtab <- jtab[jtab > 0]
    
    ws   <- which(sdata$plot == plots[j])
    stab <- colSums(sdata[drop = FALSE,ws,seedNames], na.rm  = TRUE)
    sname <- names(stab)
    unkn <- grep('UNKN', sname)
    UNKN <- sname[unkn]
    
    jvec <- matrix(jtab)
    JJ   <- crossprod(t(jvec)) + diag(.01, length(jvec))
    rr <- crossprod(matrix(stab, nrow=1),t(jvec))%*%solve( JJ )
    rownames(rr) <- names(stab)
    colnames(rr) <- names(jtab)
    rr <- t(rr)
    
    if(length(rr) == 1){
      rr[1] <- 1
      rownames(rr) <- paste(names(jtab), plots[j],sep='-')
      priorR <- rbind(priorR,rr)
    }else{
      
      for(m in 1:nrow(rr)){
        
        # no seeds counted
        if(sum(rr[m,]) == 0){ 
          if(length(UNKN) > 0){
            rr[m,UNKN] <- 1
          }
          next
        }
        
        # counted as a different species
        OTHER <- FALSE
        wm <- which( rr[m,] > 0 )
        wk <- which( !names(wm) %in% c(rownames(rr), UNKN) )
        keep <- character(0)
        if(length(wk) > 0){
          keep  <- names(wm[wk])
          OTHER <- TRUE
        }
        
        # to unknown class
        if(UCOLS){
          wm <- which(!colnames(rr) %in% c(rownames(rr)[m],UNKN, keep) )
          if(length(wm) > 0)rr[m,wm] <- 0
        }
        
        # to same class
        wm <- which(colnames(rr) == rownames(rr)[m])
        if(length(wm) > 0){
          if(OTHER)rr[m,wm] <- rr[m,wm]*10
          if(UCOLS & OTHER)rr[m,unkn] <- rr[m,unkn]*10
          if(!UCOLS){
            rr[m,wm] <- 1
            rr[m,-wm] <- 0
          }
        }
      }
      rownames(rr) <- paste(names(jtab), plots[j],sep='-')
      rr <- sweep(rr, 1, rowSums(rr), '/')
      priorR <- rbind(priorR,rr)
    }
  }
  
  priorR[!is.finite(priorR)] <- 0
  
  seedCount <- as.matrix(sdata[,seedNames,drop = FALSE])

  tt <- columnSplit(rownames(priorR),'-')
  
  attr(priorR,'species') <- tt[,1]
  attr(priorR,'plot')    <- tt[,2]
  
  ws <- which(rowSums(priorR) == 0)
  if(length(ws) > 0){
    rr <- mastIDmatrix(tdata, sdata, 
                       specNames = specNames, seedNames = seedNames, 
                       verbose = verbose)$R
    if(ncol(rr) > 1)rr <- sweep(rr, 1, rowSums(rr), '/')
    priorR[ws,colnames(rr)] <- rr[attr(priorR,'species')[ws],]
  }
  
  priorRwt <- priorR*10
  
  posR <- which(!priorR %in% c(0, 1))
  
  if( all(priorR %in% c(0, 1)) )SAMPR <- FALSE
  
  return( list(SAMPR = SAMPR, R = priorR, priorR = priorR, priorRwt = priorRwt,
               seedCount = seedCount, posR = posR, tdata = tdata) )
}

setupZ <- function(tdata, xytree, specNames, years, minD, maxD, maxFec, CONES, 
                   seedTraits=NULL, verbose){
  
  SEEDDATA <- TRUE
  if(is.null(xytree))SEEDDATA <- FALSE
  
  years <- min(tdata$year,na.rm  = TRUE):max(tdata$year,na.rm  = TRUE)
  maxF  <- specPriorVector(maxFec, tdata)
  nspec <- length(specNames)
  
  tdata$treeID <- columnPaste(tdata$plot,tdata$tree)
  
 # tid  <- columnPaste(tdata$plot,tdata$tree)
  tids <- unique(tdata$treeID)
  dcol <- match(tdata$treeID,tids) #do again after reorder
  ntree <- length(tids)
  
  nyr <- length(years)
  
  # initialize repro
  if(!'repr' %in% colnames(tdata)){
    tdata$repr <- NA
  }else{
    tdata$repr[tdata$repr < .5] <- 0
    tdata$repr[tdata$repr >= .5] <- 1
  }
  
  if( 'cropMax' %in% colnames(tdata) ){
    tdata$repr[ tdata$cropMax > 0 ]  <- 1
 #   tdata$repr[ tdata$cropMax == 0 ] <- 0
  }
  if( 'cropMin' %in% colnames(tdata) ){
    tdata$repr[ tdata$cropMin > 0 ] <- 1
  }
  
  if( !'fecMin' %in% colnames(tdata) ){
    tdata$fecMin <- 1e-4
  }else{
    wm <- which(is.na(tdata$fecMin))
    if(length(wm) > 0){
      tdata$fecMin[wm] <- 1e-4
    }
  }
  if(!'fecMax' %in% colnames(tdata)){
    tdata$fecMax <- maxF
  }else{
    wm <- which( is.na(tdata$fecMax) | tdata$fecMax == 0 )
    if(length(wm) > 0){
      tdata$fecMax[wm] <- maxF[wm]
    }
  }
  tdata$fecMin[ tdata$fecMin < 1e-4 ] <- 1e-4
  
  fstart <- rep(NA, nrow(tdata))
  if('lastFec' %in% names(tdata))fstart <- tdata$lastFec
  
  if(CONES){
    
    tdata$repr[tdata$cropCount > 0] <- 1
    
    if( is.null(seedTraits) ){
      warning('cannot use treeData$cropCount without inputs$seedTraits, 
              assumed = 1')
      seedTraits <- matrix(1, nspec, 1)
      rownames(seedTraits) <- .fixNames(specNames, MODE='character')$fixed
      colnames(seedTraits) <- 'seedsPerFruit'
    }
    seedTraits[,'seedsPerFruit'] <- ceiling(seedTraits[,'seedsPerFruit'])
    
    if( !'cropFraction' %in% colnames(tdata) ){
      kwords <- "\nNote: missing column _cropFraction_, assumed = 1\n"
      words <- c(words, kwords)
      tdata$cropFraction <- .95
    }
    tdata$cropFraction[tdata$cropFraction > .99] <- .99
    tdata$cropCount <- ceiling(tdata$cropCount)
    
    if(!'cropFractionSd' %in% colnames(tdata)){
      tdata$cropFractionSd <- NA
    }
    
    wm <- which( is.finite(tdata$cropFraction) & 
                   !is.finite(tdata$cropFractionSd) )
    if(length(wm) > 0){
      fs <- .1*dbeta(tdata$cropFraction[wm], .1,2) + 1e-3
      tdata$cropFractionSd[wm] <- signif(fs, 3)
    }
    
    ww <- which( tdata$cropCount == 0 & tdata$cropFraction == 0 )
    if(length(ww) > 0){
      if(verbose){
        cat('\nNote: deleted finite cropCount with cropFraction = 0:\n')
        print(tdata$treeID[ww])
        tdata$cropFraction[ww] <- tdata$cropCount[ww] <- NA
      }
    }
      
    ww <- which( is.finite(tdata$cropCount) & tdata$cropFraction == 0)
    if(length(ww) > 0){
      if(verbose)cat('\nNote: deleted finite cropCount with cropFraction = 0:\n')
      print(tdata$treeID[ww])
      tdata$cropFraction[ww] <- tdata$cropCount[ww] <- NA
    }               
  }
  

  if('serotinous' %in% colnames(tdata)){
    
    ww <- which(tdata$serotinous == 1)
    if(length(ww) > 0){
      snames <- unique(tdata$treeID[ww])
      tdata$repr[tdata$treeID %in% snames] <- 1
    }
  }
  
  iy  <- match(tdata$year, years)
  nyr <- length(years)
  
  zknown <- matrix(NA, ntree, nyr)
  rownames(zknown) <- tids
  colnames(zknown) <- years
  
  dmat <- dminMat <- dmaxMat <- zknown
  zknown[ cbind(dcol, iy) ]  <- tdata$repr 
  dmat[ cbind(dcol, iy) ]    <- tdata$diam
  dminMat[ cbind(dcol, iy) ] <- minD
  dmaxMat[ cbind(dcol, iy) ] <- maxD
  
  mdc <- apply(dminMat,1,max, na.rm  = TRUE)
  dminMat[,1:nyr] <- mdc
  mdc <- apply(dmaxMat,1,max, na.rm  = TRUE)
  dmaxMat[,1:nyr] <- mdc
  
  mmin <- apply(dmat,1,min, na.rm  = TRUE)
  mmax <- apply(dmat,1,max, na.rm  = TRUE)
  
  mtmp <- ntmp <- dminMat*0   
  mtmp[,1:nyr] <- mmin - 1.5  # could be replaced with growth trend
  ntmp[,1:nyr] <- mmax + 1.5
  
  dtmp <- dmat
  dtmp[is.na(dtmp)] <- 0
  mcum <- t( apply(dtmp, 1, cumsum))
  dmat[mcum == 0] <- mtmp[mcum == 0]
  dmat[is.na(dmat)] <- ntmp[is.na(dmat)]
  
  #after first observed mature
  zyr <- zknown
  zyr[is.na(zyr)] <- 0
  zyr <- t( apply(zyr, 1, cumsum))
  zyr[zyr > 1] <- 1
  zknown[zyr == 1] <- 1
  
  last0 <- rep(0, nrow(zknown))
  names(last0) <- tids
  first1 <- last0 + ncol(dmat) + 1
  
  znew <- zknown
  
  for(k in 1:nyr){   # obs repr
    zk <- zknown[,k]
    ww <- which(zk == 0)
  #  if(k > 1)ww <- which(zk == 0 & zknown[,k-1] == 0)
    last0[ww] <- k
    ww <- which(zk == 1 & first1 > k)
    first1[ww] <- k                 # known mature
    w0 <- which( dmat[,k] < dminMat[,k] & first1 > k) # small and not yet obs mature
    if(length(w0) > 0)last0[w0] <- k
    
    zk[w0] <- 0
    znew[,k] <- zk
  }
  
  for(k in nyr:1){                                                # note reverse
    w1 <- which(dmat[,k] > dmaxMat[,k] & last0 < k & first1 > k)  # large and after last obs immature
    if(length(w1) > 0)first1[w1] <- k
    znew[w1,k] <- 1
  }
  
  zknown <- znew
  
  all0 <- all1 <- last0*0
  all0[last0 >= nyr] <- 1
  all1[first1 == 1] <- 1
  
  # mature first yr in data?
  fyr <- tapply(tdata$year, tdata$treeID, min)  # 1st yr in data
  fyr <- fyr[names(last0)]
  f1  <- first1
  f1[f1 > length(years)] <- length(years) - 1
  zyr <- years[f1]                          # 1st yr mature
  
  ww <- which(fyr == zyr)   
  if(length(ww) > 0)all1[ww] <- 1
  
  # is last yr in a data a zero?
  lyr <- tapply(tdata$year, tdata$treeID, max)  # 1st yr in data
  lyr <- lyr[names(last0)]
  l0  <- last0
  l0[l0 == 0] <- 1
  zyr <- years[l0]
  
  ww <- which(lyr == zyr)   
  if(length(ww) > 0)all0[ww] <- 1
  
  last0[all0 == 1] <-  ncol(zknown)
  
  last0first1 <- cbind(last0, first1, all0, all1)
  
  #initial values
  zmat <- zknown
  zmat[all0 == 1,] <- 0
  zmat[all1 == 1,] <- 1
  
  matYr <- round( (last0 + first1)/2 )
  matYr[matYr == 0] <- 1
  matYr[all0 == 1] <- ncol(zmat) + 1
  
  for(k in 1:length(years)){
    zk <- zknown[,k]
    zk[which(k <= last0)] <- 0
    zk[which(k >= first1)] <- 1
    zknown[,k] <- zk
    
    zk <- zmat[,k]
    wna <- which(is.na(zk))
    zk[k < matYr] <- 0
    zk[k >= matYr] <- 1
    zmat[,k] <- zk
  }
  
  tids <- unique(tdata$treeID)
  dcol <- match(tdata$treeID, tids)
     
  tyindex <- cbind(dcol, iy)  #tree-yr index
  z    <- zmat[ tyindex ]
  
  zknownVec <- zknown[tyindex]
  
  tdata$repr <- zknownVec
  tdata$fecMin[ is.na(zknownVec) ] <- 1e-4
  tdata$fecMax[ is.na(zknownVec) ] <- maxF[ is.na(zknownVec) ]
  ww <- which(zknownVec == 1)
  if(length(ww) > 0)tdata$fecMax[ ww ] <- maxF[ ww ]
  tdata$fecMin[ zknownVec == 0 ] <- 1e-4
  
  ww <- which(tdata$fecMin < 1 & tdata$repr == 1)
  if(length(ww) > 0)tdata$fecMin[ww] <- 1
  
  ww <- which(is.na(tdata$fecMin) & tdata$repr == 1)
  if(length(ww) > 0)tdata$fecMin[ww] <- 1
  
  tdata$fecMax[tdata$fecMax < 1] <- 1
  tdata$fecMin[tdata$fecMin < 1e-4] <- 1e-4
  
  last   <- which(last0first1[,'all0'] == 1)
  snames <- rownames(last0first1)[last]  # always immature
  
  scc <- numeric(0)
  fstart[ is.na(fstart) & tdata$repr == 0 ] <- .01
  
  
  if(CONES){
    
    ww <- which(is.finite(tdata$cropCount) &
                  !is.finite(tdata$cropFraction))
    if(length(ww) > 0)tdata$cropFraction[ww] <- .95
    
    ww <- which(!is.finite(tdata$cropCount) &
                  is.finite(tdata$cropFraction))
    if(length(ww) > 0)tdata$cropFraction[ww] <- NA
    
    ww <- which(!is.finite(tdata$cropFractionSd) &
                  is.finite(tdata$cropFraction))
    if(length(ww) > 0)tdata$cropFractionSd[ww] <- 0
    
    ww <- which(is.finite(tdata$cropFraction))
     
    cll <- seedTraits[tdata$species[ww],'seedsPerFruit']*tdata$cropCount[ww] # minimum at least those counted
    scc <- cll/tdata$cropFraction[ww]
    
    chi <- scc*2
    clo <- scc*.5
    clo[clo < cll] <- cll[clo < cll]
    clo[clo < 1 & cll >= 1] <- 1
    chi[ chi < 10 ] <- 10
    
    tdata$fecMin[ww] <- clo
    tdata$fecMax[ww] <- chi
    
    ww <- which(tdata$fecMax <= tdata$fecMin)
    if(length(ww) > 0){
      print('fecMax < fecMin')
      print(tdata[ww[1:10],])
      stop()
    }
    
    fstart[ww] <- (clo + 2*scc + chi)/4
    
    # non-zero cones cannot be always immature
    specMatr <- unique(tdata$treeID[which(tdata$cropCount > 0)])
    wm <- which(last0first1[,'all0'] == 1 &
                  rownames(last0first1) %in% specMatr)
    if(length(wm) > 0)last0first1[wm,'all0'] <- 0
  }
  
  
  if('cropMin' %in% colnames(tdata)){
    ww <- which(is.finite(tdata$cropMin))
    fs <- ( tdata$fecMin[ww] +  tdata$fecMax[ww] )/2
    fs[ !is.finite(fs) ] <- tdata$fecMin[ww][ !is.finite(fs) ]
    fs[ fs < 1e-4 ] <- 1e-4
    fstart[ww] <- fs
  }
 
  ww <- which(!is.finite(fstart))
  if(length(ww) > 0){
      zw <- z[ww]
      fw <- zw
      tl <- tdata$fecMin[ww] + zw
      tl[ zw == 0 ] <- .01
      th <- tdata$fecMax[ww] + zw
      tl[ zw == 1 & tl < 1 ] <- 1
      th[ zw == 0] <- .99
      fm <- sqrt( tl * th )
      fw <- .tnorm( length(zw), tl, th, fm, 10 )
      fstart[ww] <- fw
  }
  
  fstart[ fstart < 1e-4 ] <- 1e-4
  ww <- which(fstart > tdata$fecMax | fstart < tdata$fecMin)
  if(length(ww) > 0){
    fstart[ww] <- .tnorm(length(ww), tdata$fecMin[ww], tdata$fecMax[ww], fstart[ww], 10)
  }
  
  # fit means inclusion in distall: excludes known immature, serotinous,
  #                                 trees not in xytree (cropCount only)
  
  if(SEEDDATA){
    snames <- c(snames, tdata$treeID[!tdata$treeID %in% xytree$treeID]) # no location
  }else{
    snames <- rownames(last0first1)  # none
  }
  ww <- which(tdata$serotinous == 1)
  if(length(ww) > 0)snames <- c(snames, tdata$treeID[ww])
  snames <- unique(snames)
  
  fit <- rep(0, nrow(last0first1))
  
  last  <- which(rownames(last0first1) %in% snames)
  first <- which(!rownames(last0first1) %in% snames)
  fit[first] <- 1
  last0first1 <- cbind(last0first1, fit)
  
  wnew   <- c(first, last)
  zmat   <- zmat[wnew,]
  zknown <- zknown[wnew,]
  matYr  <- matYr[wnew]
  last0first1 <- last0first1[wnew,]
  
  tdata$fit <- 1
    
  if(length(snames) > 0){  # there are some known immature, put immature at end
    mf <- which(!tdata$treeID %in% snames)
    ml <- which(tdata$treeID %in% snames)
    tdata$fit[ml] <- 0
    mm <- c(mf, ml)
    tdata <- tdata[mm,]
    fstart <- fstart[mm]
    z <- z[mm]
  }
  
  tdata$tnum <- match(tdata$treeID, rownames(last0first1))
  fstart[fstart < 1e-4] <- 1e-4
  fecMaxCurrent <- tdata$fecMax
  fecMinCurrent <- tdata$fecMin
  
  
  
  list(z = z, zmat = zmat, zknown = zknown, matYr = matYr, seedTraits = seedTraits,
       last0first1 = last0first1, tdata = tdata, fstart = fstart,
       fecMinCurrent = fecMinCurrent, fecMaxCurrent = fecMaxCurrent)
}
  
getPredGrid <- function(predList, tdata, sdata, xytree, xytrap, group, 
                        specNames, plotDims, trapRows, verbose = FALSE){
  
  mapMeters  <- predList$mapMeters
  mapPlot    <- predList$plots
  mapYear    <- predList$years
  
  if(is.null(mapMeters)){
    mapMeters <- 5
    if(verbose)cat('\nMissing mapMeters for prediction grid set to 5 m\n')
  }
  
  plotYrComb <- table( tdata$plot[trapRows], tdata$year[trapRows] )
  plotYrComb <- plotYrComb[drop = FALSE,rownames(plotYrComb) %in% mapPlot,]
  plotYrComb <- plotYrComb[,colnames(plotYrComb) %in% mapYear, drop = FALSE ]
  
  predList$years <- colnames(plotYrComb)[colSums(plotYrComb) > 0]
  
  npred      <- nrow(plotYrComb)
  predList$plots <- predList$plots[predList$plots %in% rownames(plotYrComb)]
  
  if(sum(plotYrComb) == 0){
    if(verbose){
      cat('\n\nPlot-years in predList missing from data:\n')
      print( paste(predList$plot,': ', predList$year, sep='') )
      cat('\n\n')
    }
    return( list(seedPred = NULL, distPred = NULL) )
  }
  
  if(length(mapMeters) == 1 & npred > 1)mapMeters <- rep( mapMeters, npred )
  
  seedPred <- numeric(0)
  drowj <- drowTot <- 0
  
  distPred <- grp <- spp <- numeric(0)
  treeid   <- trapid <- numeric(0)
  
  gridSize <- rep(0, npred)
  names(gridSize) <- rownames(plotYrComb)
  
  for(j in 1:npred){
    
    wj <- which(plotYrComb[j,] > 0)
    pj <- rownames(plotYrComb)[j]
    wy <- as.numeric(colnames(plotYrComb)[wj])
  
    jplot <- as.matrix( plotDims[rownames(plotDims) == pj,] )
    
    jx <- c(jplot['xmin',1] - mapMeters[j]/2,
            jplot['xmax',1] + mapMeters[j]/2)
    jy <- c(jplot['ymin',1] - mapMeters[j]/2,
            jplot['ymax',1] + mapMeters[j]/2)
    
    sx <- seq(jx[1],jx[2],by=mapMeters[j])
    sy <- seq(jy[1],jy[2],by=mapMeters[j])
    
    jgrid <- expand.grid(x = sx, y = sy)
    gridSize[j] <- nrow(jgrid)
    
    yrj   <- rep( wy, nrow(jgrid) )
    
    jseq  <- rep(1:nrow(jgrid), each=length(wy) )
    jgrid <- jgrid[ jseq ,]
    
    dgrid   <- drowTot + jseq
    drowTot <- max(dgrid)
    trapID  <- paste(pj,'-g',dgrid,sep='')
    
    dj  <- data.frame(trapID = trapID, year = yrj, jgrid, drow=0, dgrid = dgrid)
    dj$plot  <- pj
    
    #includes trap years from data
    wmatch <- which(sdata$plot %in% pj)
    id     <- as.character(unique( sdata$trapID[wmatch] ) )
    sdd    <- sdata[match(id,sdata$trapID),]
    
    xy  <- xytrap[match(sdd$trapID,xytrap$trapID),c('x','y')]
    jj  <- rep( c(1:nrow(sdd)), each=length(wy))
    yrj <- rep( wy, nrow(sdd) )
    drAll  <- sort(unique(sdd$drow))
    dgrid  <- drowTot + c(1:length(drAll))
    dgrid  <- dgrid[ match(sdd$drow[jj],drAll) ]
    drowTot <- max(dgrid)
    
    tj <- data.frame(trapID = sdd$trapID[jj], year = yrj, 
                     x = xy[jj,'x'], y = xy[jj,'y'],
                     drow = sdd$drow[jj], dgrid=dgrid, plot = sdd$plot[jj])
    dj <- rbind(dj,tj)
    
    seedPred <- rbind(seedPred, dj)
  }
  
  tdat <- tdata[tdata$plot %in% predList$plots &
                  tdata$year %in% predList$years,]
  
  xyt <- xytree[xytree$plot %in% predList$plots,]
  
  tmp <- setupDistMat(tdat, seedPred, xyt, seedPred, verbose)
  distPred <- tmp$distall
  seedPred <- tmp$sdata
  treePred <- tmp$tdata
  
  
  plotYrComb <- cbind(plotYrComb, mapMeters, gridSize)
  seedPred$active  <- seedPred$area <- 1 # note for 1 m2
  attr(distPred,'species') <- spp
  
 # distPred[distPred == 0] <- 100000
  distPred <- round(distPred,1)
  
  seedPred$drow <- match(as.character(seedPred$trapID), rownames(distPred))
  
  rownames(seedPred) <- NULL
  
  if(verbose){
    cat("\nPrediction grid size: ")
    cat("If too large, increase predList$mapMeters:\n")
    print( plotYrComb[, c('mapMeters','gridSize')] )
  }
  
  list(seedPred = seedPred, distPred = distPred, 
       treePred = treePred, predList = predList)
}

cleanFactors <- function(x){
  
  #fix factor levels
  
  scode <- names(x[ which(sapply( x, is.factor )) ])
  if(length(scode) > 0){
    for(j in 1:length(scode)) {
      x[,scode[j]] <- droplevels(x[,scode[j]])
    }
  }
  x
}
  
.setupRandom <- function(randomEffect, tdata, xfec, xFecNames, specNames){
  
  tdata$species <- as.factor(tdata$species)
  
  nspec      <- length(specNames)
  formulaRan <- randomEffect$formulaRan
  if(nspec > 1)formulaRan <- .specFormula(randomEffect$formulaRan)
  xx        <- .getDesign(formulaRan, tdata)$x
  if(nspec > 1)xx <- xx[,grep('species',colnames(xx)),drop = FALSE]  # CHECK for 1 spp
  
  xrandCols  <- match(colnames(xx),colnames(xfec))
  
  if( !is.finite(min(xrandCols)) )
    stop('\nthere are variables in formulaRan that are missing from formulaFec\n')
  
  Qrand      <- length(xrandCols)
  reI        <- as.character(tdata[,randomEffect$randGroups])
  rnGroups   <- unique(reI)
  reIndex    <- match(reI, rnGroups)
  names(reIndex) <- rnGroups[reIndex]
  reGroups   <- unique(reIndex)
  names(reGroups) <- names(reIndex)[match(reGroups,reIndex)]

  nRand      <- length(reGroups)
  Arand      <- priorVA <- diag(1, Qrand)
  dfA        <- ceiling( Qrand + 1  + nRand/2)
  alphaRand  <- matrix(0, nRand, Qrand)
  colnames(alphaRand) <- xFecNames[xrandCols]
  rownames(alphaRand) <- rnGroups
  
  XX <- crossprod(xx)
  diag(XX) <- diag(XX) + .00000001
  xrands2u <- solve(XX)%*%crossprod(xx,xfec[,xrandCols]) 
  xrands2u[abs(xrands2u) < 1e-8] <- 0
  
  list(formulaRan = formulaRan, xrandCols = xrandCols, Qrand = Qrand,
       rnGroups = rnGroups, reIndex = reIndex, reGroups = reGroups, 
       Arand = Arand, dfA = dfA, alphaRand = alphaRand, priorVA = priorVA,
       xrands2u = xrands2u)
}

getPlotDims <- function(xytree, xytrap){
  
  plots <- sort( unique(as.character(xytree$plot) ) )
  npp   <- length(plots)
  
  pdims <- numeric(0)
  
  for(j in 1:npp){
    
    wt <- which(xytree$plot == plots[j])
    ws <- which(xytrap$plot == plots[j])
    
    jx <- range( c(xytree$x[wt], xytrap$x[ws]) )
    jy <- range( c(xytree$y[wt], xytrap$y[ws]) )
    jx[1] <- floor( jx[1] - 1 )
    jx[2] <- ceiling( jx[2] + 1 )
    jy[1] <- floor( jy[1] - 1 )
    jy[2] <- ceiling( jy[2] + 1 )
    
    area <- diff(jx)*diff(jy)/10000
    
    if(area > 200){
      cat( paste('\nPlot area > 200 ha:', plots[j], 'is', area, 'ha\n') )
      stop('check coordinates for xytree, xytrap')
    }
    
    pdims <- rbind( pdims, c(jx, jy, area) )
  }
  colnames(pdims) <- c('xmin','xmax','ymin','ymax','area')
  rownames(pdims) <- plots
  pdims
}
  
.orderChain <- function(xchain, snames){
  
  if(!snames[1] %in% colnames(xchain))return(xchain)
  
  ns <- length(snames)
  
  mnames <- .coeffNames(colnames(xchain))
  first  <- mnames[1:ns]
  tmp    <- grep('_',first)
  
  if(length(tmp) > 0){
    first <- matrix( unlist(strsplit(first,'_')),ncol=2,byrow  = TRUE)[,2]
  }
  
  orr    <- match(snames,first)
  if(is.na(orr[1]))return(xchain)
  
  newChain <- xchain*0
  
  k <- orr
  m <- 1:ns
  while(max(k) <= ncol(xchain)){
    
    newChain[,m] <- xchain[,k]
    colnames(newChain)[m] <- colnames(xchain)[k]
    
    m <- m + ns
    k <- k + ns
  }
  newChain
}

factor2integer <- function(fvec){
  as.numeric(as.character(fvec))
}

formit <- function(form, nspec){
  
  ff   <- as.character(form)
  ff   <- .replaceString(ff, ':', '*')
  
  form <- as.formula( paste(ff, collapse=' ') )
  
  if(nspec > 1){
    fc   <- .replaceString( as.character(form), 'species *','')
    fc   <- as.formula( paste( fc, collapse = ' ') )
    form <- .specFormula(fc)
  }
  .fixFormula(form)
}

.fixFormula <- function(form){
  
  # remove I(log()) from formula
  
  fchar <- as.character(form)[2]
  tmp   <- gregexpr('I(log(', fchar, fixed  = TRUE)[[1]]
  
  if(tmp[1] < 0)return(form)
  
  if(length(tmp) > 0){
    
    while(tmp[1] > 0){
      tmp   <- gregexpr('I(log(', fchar, fixed  = TRUE)[[1]]
      end <- gregexpr(')', fchar, fixed  = TRUE)[[1]]
      we <- end[ min(which(end > tmp[1])) ]
      substr(fchar, we, we) <- " "
      substr(fchar, tmp[1], (tmp[1] + attr(tmp,'match.length')[1]) ) <- "  log("
      fchar <- .replaceString(fchar,'  ',' ')
      fchar <- .replaceString(fchar,' )',')')
      tmp   <- gregexpr('I(log(', fchar, fixed  = TRUE)[[1]]
    }
  }
  as.formula( paste('~ ',  fchar, collapse = ' ') )
}


setupPriors <- function(treeData, specNames, priorTable, priorList = NULL, 
                        priorDist = NULL, priorVDist = NULL, maxDist = NULL, minDist = NULL,
                        minDiam = NULL, maxDiam = NULL, sigmaMu = NULL, maxF = NULL, maxFec = NULL,
                        ug = NULL, priorTauWt = NULL, priorVU = NULL, ARSETUP = F, USPEC = F){
  
  if(!is.null(priorList)){
    
    if(length(priorList) > 1){ #priors by species
      for(k in 1:length(priorList)){
        wk <- which(names(priorList[[k]]) == 'priorDist')
        if(length(wk) == 1)priorDist <- unlist(priorList[[k]][wk]) 
        wk <- which(names(priorList[[k]]) == 'priorVDist')
        if(length(wk) == 1)priorVDist <- unlist(priorList[[k]][wk])
        wk <- which(names(priorList[[k]]) == 'minDist')
        if(length(wk) == 1)minDist <- unlist(priorList[[k]][wk])
        wk <- which(names(priorList[[k]]) == 'maxDist')
        if(length(wk) == 1)maxDist <- unlist(priorList[[k]][wk])
        wk <- which(names(priorList[[k]]) == 'minDiam')
        if(length(wk) == 1)minDiam <- unlist(priorList[[k]][wk])
        wk <- which(names(priorList[[k]]) == 'maxDiam')
        if(length(wk) == 1)maxDiam <- unlist(priorList[[k]][wk])
        wk <- which(names(priorList[[k]]) == 'maxF')
        if(length(wk) == 1)maxFec <- unlist(priorList[[k]][wk])
      }
    }
  }else{
    priorList <- list(priorDist = priorDist, priorVDist = priorVDist, 
                      maxDist = maxDist, minDist = minDist,
                      minDiam = minDiam, maxDiam = maxDiam, 
                      sigmaMu = sigmaMu, maxF = maxF)
  }
  
  pcols <- c("priorDist","priorVDist","minDist","maxDist","minDiam",
             "maxDiam","maxFec")
  
  nspec <- length(specNames)
  
  if(is.null(priorTable)){
    priorTable <- matrix(NA, nspec, length(pcols))
    colnames(priorTable) <- pcols
    rownames(priorTable) <- specNames
    for(k in 1:length(pcols)){
      priorTable[,pcols[k]] <- get( pcols[k][1] )
    }
  }else{
    if(nrow(priorTable) > 1 & !ARSETUP ){
      if(var(priorTable[,'priorDist']) > 0)USPEC <- TRUE
    }
    pm <- which(!pcols %in% colnames(priorTable))
    if(length(pm) > 0){
      for(k in pm){
        priorTable <- cbind(priorTable, get( pcols[k][1] ))
        colnames(priorTable)[k] <- pcols[k]
      }
    }
  }
  
  priorTable <- priorTable[drop = FALSE,specNames,]
  if(nrow(priorTable) > 1){
    w1 <- w2 <- w3 <- TRUE
    if(!is.null(priorTable[,'priorDist']))
      w1 <- var(priorTable[,'priorDist'], na.rm  = TRUE) == 0 
    if(!is.null(priorTable[,'minDist']))
      w2 <- var(priorTable[,'minDist'], na.rm  = TRUE) == 0 
    if(!is.null(priorTable[,'maxDist']))
      w3 <- var(priorTable[,'maxDist'], na.rm  = TRUE) == 0 
    if(w1 & w2 & w3)USPEC <- FALSE
  }
  if(nrow(priorTable) == 1)USPEC <- FALSE
  
  priorU  <- round( (2*priorTable[,'priorDist']/pi)^2 )
  priorVU <- round( (2/pi)^2*priorTable[,'priorVDist']^2 )
  maxU    <- round( (2*priorTable[,'maxDist']/pi)^2 )
  minU    <- round( (2*priorTable[,'minDist']/pi)^2 )
  
  priorTable <- cbind(priorTable,priorU, priorVU, minU, maxU)
  
  umean <- mean(priorTable[,'priorU'])
  propU <- mean(priorTable[,'priorU'])/100
  uvar  <- mean(priorTable[,'priorVU'])
  if(is.null(ug))ug <- mean(priorTable[,'priorU'])
  
  for(k in 1:ncol(priorTable)){
    pk <- priorTable[,k]
    names(pk) <- rownames(priorTable)
    assign(colnames(priorTable)[k], pk)
  }
  priorTable <- priorTable[,!duplicated(colnames(priorTable)),drop = FALSE]
  ug   <- priorU
  maxF <- max(maxFec)
  
  npt   <- 1:nspec
  if(!USPEC)npt <- 1
  
  if(is.null(priorTauWt))priorTauWt <- ceiling(nrow(treeData)/nspec/10)
  tau1 <- priorTauWt
  tau2 <- priorVU*(tau1 - 1)
  
  priorU  <- mean( priorU )
  priorVU <- mean( priorVU )
  
  if(nspec == 1){
    maxU    <- max( maxU )
    minU    <- min( minU )
    names(minU) <- names(maxU) <- specNames
    priorDist  <- priorDist[1]
    priorVDist <- priorVDist[1]
    maxDist    <- maxDist[1]
    minDist    <- minDist[1]
    maxFec     <- maxFec[1]
  }else{
    ug <- ug[specNames]
    minU <- minU[specNames]
    maxU <- maxU[specNames]
    minDiam <- minDiam[specNames]
    maxDiam <- maxDiam[specNames]
    maxFec  <- maxFec[specNames]
  }
  
  if(is.null(sigmaMu))sigmaMu <- 5
  sigmaWt <- sqrt(nrow(treeData))
  
  list(priorTable = priorTable, priorList = priorList, priorDist = priorDist,
       priorVDist = priorVDist, maxDist = maxDist, minDist = minDist,
       minDiam = minDiam, maxDiam = maxDiam, maxFec = maxFec,
       sigmaMu = sigmaMu, sigmaWt = sigmaWt,
       maxF = maxF, umean = umean, priorU = priorU, priorVU = priorVU, 
       minU = minU, maxU = maxU, propU = propU, uvar = uvar, ug = ug, 
       USPEC = USPEC, npt = npt, tau1 = tau1, tau2 = tau2)
}

mastif <- function(inputs, formulaFec=NULL, formulaRep=NULL, 
                   ng = NULL, burnin = NULL, predList = NULL, 
                   yearEffect = NULL, randomEffect = NULL, 
                   modelYears = NULL, plotDims = NULL){   
  data  <-  NULL
  
  if( class(inputs) == 'mastif' ){
    
    inputs$inputs$ng <- ng
    inputs$inputs$burnin <- burnin
    
    parameters <- inputs$parameters
    priorTable <- inputs$inputs$priorTable
    
    data   <- inputs$data
    if( !is.null(modelYears) ){
      inputs$inputs$tdataOut <- inputs$prediction$tdataOut
      inputs$inputs$sdataOut <- inputs$prediction$sdataOut
    }
    inputs <- inputs$inputs
    inputs$parameters <- parameters
    inputs$priorTable <- priorTable
    class(inputs) <- 'mastif'
  }

  if(is.null(ng))stop("\nsupply no. MCMC steps, 'ng'\n")
  if(is.null(burnin))stop("\nsupply 'burnin'\n")
  
  .mast(inputs, data, formulaFec, formulaRep, predList, yearEffect, 
        randomEffect, modelYears, ng, burnin) 
}
   
.mast <- function(inputs, data, formulaFec, formulaRep, predList, yearEffect,
                  randomEffect, modelYears, ng, burnin){
   
  upar <- xytree <- xytrap <- specNames <- treeData <- seedData <-
    seedNames <- arList <- times <- xmean <- xfecCols <- xrepCols <-
    groupByInd <- dfA <- xrands2u <- lagGroup <- lagMatrix <- xfecs2u <-
    xreps2u <- Qrand <- xfecU <- xrepU <- seedTraits <- 
    plotDims <- plotArea <- tdataOut <- sdataOut <- specPlots <- 
    plotNames <- distall <- trapRows <- fstart <- 
    fecMinCurrent <- fecMaxCurrent <- NULL
  
  notFit <- NULL
  maxU <- minU <- npt <- priorU <- priorVU <- tau1 <- tau2 <- NULL
  censMin <- censMax <- NULL
  SEEDCENSOR <- CONES <- RANDOM <- YR <- AR <- ARSETUP <- USPEC <- 
    TREESONLY <- FECWT <- verbose <- SEEDDATA <- SAMPR <- FALSE
  
  words <- inwords <- character(0)
  
  priorList <- priorTable <- NULL

  priorDist <- 25; priorVDist <- 40; maxDist <- 70; minDist <- 4
  minDiam <- 10; maxDiam <- 40; maxF = 1e+8; sigmaMu <- 1
 
  plag  <- p  <- 0; maxFec   <- 1e+8; priorVtau <- 6
  ug <- priorTauWt <- NULL
  alphaRand <- Arand <- priorB <- priorIVB <- betaPrior <- NULL
  
  if('seedData' %in% names(inputs))SEEDDATA <- TRUE
  
  PREDSEED <- TRUE
  if(is.null(predList))PREDSEED <- FALSE
  
  betaYr <- betaLag <- yeGr <- plots <- years <- NULL
  facLevels <- character(0)
  ngroup <- 1
  nng    <- ng
  
  if( class(inputs) == 'mastif' ){
    
    ARSETUP <- TRUE
    
    ww <- which(!names(inputs) %in% c('inputs','chains','fit',
                                      'burnin', 'ng', 'predList'))
    
    for(k in ww)assign( names(inputs)[k], inputs[[k]] )
 #   tdata <- treeData
    
    for(k in 1:length(data$setupData)){
      assign( names(data$setupData)[k], data$setupData[[k]] )
    }
    if('arList' %in% names(data)){
      for(k in 1:length(data$arList))
        assign( names(data$arList)[k], data$arList[[k]] )
      AR <- TRUE
    }
    if('setupRandom' %in% names(data)){
      for(k in 1:length(data$setupRandom))
        assign( names(data$setupRandom)[k], data$setupRandom[[k]] )
      RANDOM <- TRUE
    }
    if('setupYear' %in% names(data)){
      for(k in 1:length(data$setupYear)){
        if(names(data$setupYear)[k] == 'yrIndex')next
        assign( names(data$setupYear)[k], data$setupYear[[k]] )
      }
      YR <- TRUE
    }
    ug <- inputs$parameters$upars[specNames,1]
    ug[is.na(ug)] <- ug[1] 
    
    if(length(ug) > 1){
      names(ug) <- specNames
      if( !( diff(range(ug)) == 0) )USPEC <- TRUE
    }
    
    upar <- ug
    years <- sort(unique(treeData$year))
    years <- min(years):max(years)
    nyr <- length(years)
    if(!is.null(predList)){
      predList$years <- predList$years[predList$years %in% years]
      if('plots' %in% names(predList))
        predList$plots <- .fixNames(predList$plots, all  = TRUE)$fixed
    }
    yrIndex <- yrIndex[,!duplicated(colnames(yrIndex))]
    
    R <- parameters$rMu
    
    seedTable   <- inputs$seedByPlot
    matYr       <- inputs$matYr 
    last0first1 <- inputs$last0first1
    
    zmat <- matrix(0, nrow(last0first1), nyr)
    zmat[ cbind(1:nrow(zmat), matYr) ] <- 1
    zmat <- t( apply(zmat, 1, cumsum) )
    zmat[zmat > 1] <- 1
    rownames(zmat) <- rownames(last0first1)
    
    ij <- cbind( match(treeData$treeID, rownames(zmat)), match(treeData$year, years) )
    z <- zmat[ij]
    
  }else{ # not class(inputs) == 'mastif'
    
    if(!is.null(yearEffect)){
      if('p' %in% names(yearEffect))plag <- yearEffect$p
      
      if('groups' %in% names(yearEffect)){
        
        ygg <- inputs$treeData[,yearEffect$groups,drop = FALSE]
        if(ncol(ygg) > 1)ygg <- columnPaste(ygg[,1],ygg[,2])
        yee <- table(ygg)
        yee <- yee[yee > 10]
        if(length(yee) < 2){
          if(verbose){
            cat('\nCannot use random groups specified in yearEffect:\n')
            print(yearEffect$groups)
            print(yee)
          }
          yearEffect <- yearEffect[!names(yearEffect) == 'groups']
        }
      }
    }
    
    if(!'FILLED' %in% names(inputs)){
      inputs <- mastFillCensus(inputs, p = plag)  
    }
    
    for(k in 1:length(inputs))assign( names(inputs)[k], inputs[[k]] )
    years <- unique( range( treeData$year ) )
    
    
    if(length(years) == 1)stop( '\nmust have > 1 year of data\n' )
    
    if(SEEDDATA){
      priorR <- mastIDmatrix( inputs$treeData, inputs$seedData, 
                              specNames = inputs$specNames,
                              seedNames = inputs$seedNames,
                              censMin = inputs$censMin, verbose = verbose)$R
      if(is.matrix(priorR)){
        inputs$specNames <- specNames <- rownames(priorR)
        inputs$seedNames <- seedNames <- colnames(priorR)
      }
      
      if(!is.null(censMin))SEEDCENSOR <- TRUE
      
      years <- range(c(treeData$year,seedData$year))
    }
    
    words <- c(words, inwords)
    years <- years[1]:years[2]
  }
  
  keepIter <- 4000
  
  plots <- .fixNames( sort(unique(as.character(treeData$plot))), all  = TRUE)$fixed
  nspec <- length(specNames)
  
  if(!is.null(randomEffect)){
    randomEffect$formulaRan <- .fixFormula(randomEffect$formulaRan)
    if('randGroups' %in% names(randomEffect)){
      if(randomEffect$randGroups == 'tree')randomEffect$randGroups <- 'treeID'
      randGroups <- randomEffect$randGroups
    }
  }
  
  if(!SEEDDATA){
    PREDSEED <- FALSE
    predList <- NULL
  }
  
  if(!is.null(predList)){
    
    PREDSEED <- TRUE
    
    if(!'plots' %in% names(predList))stop('\npredList must include plots\n')
    
    predList$plots <- .fixNames(predList$plots, all  = TRUE)$fixed
    predList$plots <- predList$plots[predList$plots %in% plots]
    if(length(predList$plots) == 0)
      stop('\nPrediction plots do not occur in treeData\n')
    if(!'mapGrid' %in% names(predList))predList$mapGrid <- 5
  }
  
  if(!is.null(yearEffect)){
    YR <- TRUE
    if('p' %in% names(yearEffect)){
      plag <- yearEffect$p
    }
  }
  
  if(plag > 0){
    AR <- TRUE
    YR <- FALSE
  }
  
  ng <- nng
  plots <- sort(unique(as.character(treeData$plot)))
  
  if(SEEDDATA){
    tmp <- checkPlotDims(plots, years, xytree, xytrap, plotDims, plotArea)
    plotDims <- tmp$plotDims
    plotArea <- tmp$plotArea
  }
  
 # sigmaWt <- sqrt(nrow(treeData))
  
  plist <- setupPriors(treeData, specNames, priorTable, priorList, priorDist, 
                          priorVDist, maxDist, minDist,
                          minDiam, maxDiam, sigmaMu, maxF, maxFec,
                          ug, priorTauWt, priorVU, ARSETUP, USPEC)
  for(k in 1:length(plist))assign( names(plist)[k], plist[[k]] )
  
  if(ARSETUP)ug <- upar
  inputs$USPEC <- USPEC
  
  if(verbose){
    cat('\nPrior parameter values:\n')
    print(priorTable[npt,])
  }
  
  if(!is.null(yearEffect))plag <- yearEffect$p
  
  
  formulaFec <- formit(formulaFec, nspec)
  formulaRep <- formit(formulaRep, nspec)
  
  if( !is.null(modelYears) ){
    
    inputs$modelYears <- modelYears
    
    wtree <- which(treeData$year %in% modelYears)
    wtrap <- which(seedData$year %in% modelYears)
    
    tdataOut <- treeData[-wtree,]
    sdataOut <- seedData[-wtrap,]
    sdataOut$seedM2  <- round(rowSums( as.matrix(sdataOut[,seedNames]) )/
                                sdataOut$area,1)
    
    xy <- xytrap[ match(sdataOut$trapID, xytrap$trapID),c('x','y')]
    sdataOut <- cbind(xy, sdataOut)
  }
  
  tdata <- treeData
  sdata <- seedData
  
  ccone <- grep('cropCount', colnames(tdata))
  if(length(ccone) != 0){
    CONES <- TRUE
    
    if(!'cropFraction' %in% colnames(tdata)){
      warning('no treeData$cropFraction provided for treeData$cropCount, assumed = .95')
      tdata$cropFraction <- NA
      tdata$cropFraction[is.finite(tdata$cropCount)] <- .95
    }
    tdata$cropCount <- ceiling(tdata$cropCount)
  }
  
  if(!'seedTraits' %in% names(inputs)){
    if(CONES)warning("inputs$seedTraits matrix not found for treeData$cropCount")
    seedTraits <- matrix(1, nspec, 2)
    colnames(seedTraits) <- c('gmPerSeed','seedsPerFruit')
    rownames(seedTraits) <- specNames
  }
  seedTraits[,'seedsPerFruit'] <- ceiling(seedTraits[,'seedsPerFruit'])
  inputs$seedTraits <- seedTraits
  
  ##################
  rm(treeData)
  rm(seedData)
  ##################
  
  if( !ARSETUP ){ 
    
    notFit <- betaPrior$notFit
    
    tmp <- .setupData(formulaFec, formulaRep, mainEffects, tdata, sdata,
                      xytree, xytrap, specNames, seedNames, AR, YR, yearEffect, 
                      minDiam, maxDiam, TREESONLY, maxFec, CONES, 
                      notFit, seedTraits = seedTraits, verbose)
    for(k in 1:length(tmp))assign( names(tmp)[k], tmp[[k]] ) 
    yeGr   <- as.character(yeGr)
    ngroup <- length(yeGr)
    tdata$treeID <- as.character(tdata$treeID)
    
    
    notFit  <- notFit[ notFit %in% colnames(xfec) ]
    notCols <- match(notFit, colnames(xfec))
    
    tid  <- unique(tdata$treeID)
    tnum <- match(tdata$treeID, tid)
    tdata$tnum <- tnum
    if('tnum' %in% colnames(yrIndex)){
      yrIndex[,'tnum'] <- tnum
    }else{
      yrIndex <- cbind(yrIndex, tnum)
    }
    ntree <- length(tid)
    
    
    setupData <- tmp[ !names(tmp) %in% 
                        c("tdata","sdata","seedNames","specNames",
                          "xytree","xytrap") ]
    
    data      <- append(data, list(setupData = setupData))
    
    seedTable <- buildSeedByPlot(sdata, seedNames)
    
    if(length(censMin) > 0){
      
      tmp <- trimCens(sdata, censMin, censMax)
      censMin <- as.data.frame(tmp$censMin)
      censMax <- as.data.frame(tmp$censMax)
      
      ctmp <- censMin
      ctmp$plot <- sdata[rownames(censMin),'plot']
      
      censTable <- buildSeedByPlot(ctmp, seedNames)
      
      rownames(censTable) <- .replaceString(rownames(censTable),'seeds','min')
      
      wk <- which(!colnames(seedTable) %in% colnames(censTable))
      if(length(wk) > 0){
        moreCols <- matrix(0, nrow(censTable), length(wk))
        colnames(moreCols) <- colnames(seedTable)[wk]
        censTable <- cbind(censTable, moreCols)
        seedTable <- rbind(seedTable, censTable[drop = FALSE,,colnames(seedTable)])
        attr(seedTable,'caption') <- 
          'min_seedNames rows give sum of minimum values for censored traps'
      }
    }
    
    if(verbose & SEEDDATA){
      cat('\nSeed count by plot:\n')
      print(seedTable)
  #    if(sum(seedTable) < 10)stop('\nnot enough seeds\n')
    }
    
    if(AR){
      data      <- append(data, list(arList = arList))
      for(k in 1:length(arList))assign( names(arList)[k], arList[[k]] ) 
      nyrAR    <- length(times)
      
      groupYr  <- yrIndex[,'group'] # for AR, group and groupYr identical
      yrIndex  <- cbind(yrIndex,groupYr)
      preYr    <- years[-c(1:plag)]
    }
    
    years <- range(tdata$year)
    years <- years[1]:years[2]

    if(YR){
      setupYear <- list(yeGr = yeGr, yrIndex = yrIndex)
      
      data     <- append(data, list(setupYear = setupYear))
      
      betaYr <- matrix(0, ngroup, length(years))
      rownames(betaYr) <- yeGr
      colnames(betaYr) <- years
      
      if(!'groupName' %in% colnames(tdata)){
        tdata$groupName <- yeGr[1]
        
        group <- match( as.character(tdata$groupName), yeGr )
        
        tdata$group <- yrIndex[,'group'] <- group
      }
      
      betaYrR  <- betaYr
      betaYrF  <- betaYrR[drop = FALSE,1,]
    }
  } ################
   
  notFit  <- notFit[ notFit %in% colnames(xfec) ]
  notCols <- match(notFit, colnames(xfec))
  
  tid   <- unique(tdata$treeID)
  ntree <- length(tid)
  
  years <- range(tdata$year)
  years <- years[1]:years[2]
  
  if(!'dcol' %in% colnames(yrIndex)){
    yrIndex <- cbind(yrIndex,tdata$dcol)
    colnames(yrIndex)[ncol(yrIndex)] <- 'dcol'
  }
  
  yrIndex[,'dcol'] <- tdata$dcol
  yrIndex[,'year'] <- match(tdata$year,years)
  if( is.list(yrIndex) )yrIndex <- as.matrix(yrIndex)
  
  obsYr <- obsTimes <- NULL
  
  if(SEEDDATA){
    if(is.null(upar))upar <- ug
    
    nseed    <- nrow(sdata)
    nsobs    <- table(sdata$plot)
    ntrap    <- nrow(xytrap)
    obsYr    <- sort(unique(tdata$year[tdata$obsTrap == 1])) # there are sdata!
    obsTimes <- match(obsYr,years)        # when there are trap years   
  }
  
  nplot <- length(plots)
  n     <- nrow(xfec)
  ntobs <- table(tdata$plot)
 
  ttab  <- table(tdata$plot, tdata$year)
  wtab  <- which(ttab > 0, arr.ind  = TRUE) 
 
  nyr    <- length(years)
  
  RANDYR <- FALSE
  tdata$species <- as.character(tdata$species)
  spp <- match( tdata$species, specNames )
  yrIndex <- cbind(yrIndex, spp)

  tdata     <- cleanFactors(tdata)
  plotYears <- sort(unique(tdata$plotYr))
  
  nacc   <- length(years)
  nspec  <- length(specNames)
  
  yrIndex <- yrIndex[,!duplicated(colnames(yrIndex))]
  
  UNSTAND <- TRUE
  if(length(xmean) == 0)UNSTAND <- FALSE
  
  Qfec <- ncol(xfec)
  Qrep <- ncol(xrep)
  xFecNames <- colnames(xfec)
  xRepNames <- colnames(xrep)
  
  nSpecPlot <- max(yrIndex[,'specPlot'])
  
  npacf   <- ceiling(nacc/2)
  pacfMat <- matrix(0, nSpecPlot, npacf)
  rownames(pacfMat) <- as.character(specPlots)
  colnames(pacfMat) <- paste('lag',0:(npacf-1),sep='-')

  pacf2  <- pacfMat
  acfMat <- pacfMat
  pacN   <- pacfMat
  
  if(YR | AR){
    ngroup <- length(yeGr)
    muyr <- rep(0,nrow(tdata))
    if(ngroup > 1)RANDYR <- TRUE
  }
  
  if( !ARSETUP ){        
    # random effects
    rnGroups <- reIndex <- reGroups <- priorVA <- NULL
    alphaRand <- Arand <- NULL
    if( !is.null(randomEffect) ){
      RANDOM     <- TRUE
      if('tree' %in% randomEffect$randGroups)
        randomEffect$randGroups[randomEffect$randGroups == 'tree'] <- 'treeID'
      
      tmp <- .setupRandom(randomEffect, tdata, xfec, xFecNames, specNames) 
      
      we <- which( table(tmp$reIndex) > 2)   # replication for random groups?
      if(length(we) < 2){                    # less than 2 groups
        RANDOM <- FALSE
        randomEffect <- NULL
        if(verbose)print('too few groups for randomEffect')
      }else{
        for(k in 1:length(tmp))assign( names(tmp)[k], tmp[[k]] ) 
        data <- append(data, list(setupRandom = tmp))
        if(length(Arand) == 1)ONEA <- TRUE
      }
    }
  }        
  ##############################################
  if(is.null(yeGr))yeGr <- specNames[1]
  
  ONEF <- ONER <- ONEA <- FALSE     
  if(ncol(xfec) == 1)ONEF <- TRUE  # intercept only
  if(ncol(xrep) == 1)ONER <- TRUE  # intercept only
  rownames(yrIndex) <- NULL
  
  
  obsRows <- which(tdata$obs == 1) ############check
  
  if(SEEDDATA){
    
    tdata$obs[tdata$plotYr %in% sdata$plotYr] <- 1  # note: more than just tdata$fit == 1
    
    seedPredGrid <- distPred <- NULL
    nseedPred <- 0
    
    # species to seed type
    tmp <- .setupR(sdata, tdata, seedNames, specNames, verbose = verbose)
    R         <- tmp$R
    priorR    <- tmp$priorR
    priorRwt  <- tmp$priorRwt
    SAMPR     <- tmp$SAMPR
    seedCount <- tmp$seedCount
    posR      <- tmp$posR
    
    obsRowSeed  <- which(sdata$year %in% obsYr)
    obsTrapRows <- sort(intersect( obsRows, trapRows ) )
    
    if( is.null(names(ug)) )names(ug) <- specNames[1:length(ug)]
    
    if(PREDSEED){
      
      if( is.character(predList$years) ) predList$years <- 
          as.numeric(predList$years)
      
      tmp <- getPredGrid(predList, tdata[obsTrapRows,], sdata, xytree, xytrap, 
                         group = yrIndex[,'group'], specNames, plotDims,
                         trapRows = trapRows)
      predList  <- tmp$predList
      sdatPred  <- tmp$seedPred
      distPred  <- tmp$distPred
      tdatPred  <- tmp$treePred
      nseedPred <- nrow(sdatPred)
      
      rownames(sdatPred) <- columnPaste(sdatPred$trapID, sdatPred$year, '_')
      
      if( !is.null(nseedPred) ){ 
        
        tdatPred <- tdatPred[,c('plot','treeID','dcol','species','specPlot','year')]
        tdatPred$row <- match( rownames(tdatPred), rownames(tdata) )
        
        plotYrs <- sort(unique(tdata$plotYr))
        
        tdatPred$plotYr <- columnPaste(tdatPred$plot,tdatPred$year,'_')
        tdatPred$plotyr <- match(tdatPred$plotYr,plotYrs)
        sdatPred$plotYr <- columnPaste(sdatPred$plot,sdatPred$year,'_')
        sdatPred$plotyr <- match(sdatPred$plotYr,plotYrs)
      }else{
        PREDSEED <- FALSE
      }
    }
  }
  
  
  if(YR){
    cnames  <- paste('yr',1:nyr,sep='-')
    sygibbs <- matrix(0,ng, nyr)
    colnames(sygibbs) <- cnames        # variance, random years
  }
  
  if(AR)Alag   <- diag(.1, plag, plag)
  
  ##################### fecundity subset
  # tdata$fit == 1:    could be mature, on seed-trap plot, not serotinous
  # trapRows:          tdata rows to include in kernel (have seedData)
  #                    exclude tdata$fit == 0 and tdata$serotinous == 1
  # obsTrapRows:       trapRows in years with seedData
  ####################
  
  yrIndex[,'tnum'] <- tdata$tnum
  yrIndex[,'year'] <- match(tdata$year, years)
  
  
  if( 'lastFec' %in% names(tdata) ){
    
    fs <- tdata$lastFec
    fs[ !is.na(fs) & fstart > 1 ] <- fstart[ !is.na(fs) & fstart > 1 ] # already have lastFec
    
    fs <- tdata$lastFec
    fs[fs == 0] <- 1e-4
    
    wu <- which(fs < 1 & z == 1)
    if(length(wu) > 0)fs[wu] <- .tnorm(length(wu), 1, 5, 1.3, 10)
    
    wu <- which(fs > 1 & z == 0)
    if(length(wu) > 0)fs[wu] <- NA
    
    fstart <- fs
    fg <- fstart
  }
  
  wna <- which(is.na(fstart))
  
  if(length(wna) > 0 & SEEDDATA){
    
    fg <- .initEM(last0first1, yeGr, distall, ug[1], tdata, sdata, 
                  specNames, seedNames, R, SAMPR, USPEC, years, trapRows, 
                  plotYears, z, xfec, fstart, verbose)
  }else{
    fg <- fstart
  }
  tdata$fecMin[ tdata$fecMin < 1e-4 ] <- 1e-4
  
  propF <- fg/10
  propF[propF < .0001] <- .0001
  
  if('groupName' %in% colnames(tdata) & verbose){
    cat('\nTree-years by plot and random group name\n')
    print(table(tdata$plot, tdata$groupName))
  }
  
  # reg variance
  sg <- sigmaMu
  s1 <- sigmaWt
  s2 <- sigmaMu*(s1 - 1)
  
  bgFec <- matrix(0, Qfec, 1)
  bgRep <- matrix(0, Qrep, 1)
  rownames(bgFec) <- xFecNames
  rownames(bgRep) <- xRepNames
  
  if(!is.null(betaPrior)){
    betaPrior <- .getBetaPrior(betaPrior, bgFec, bgRep, specNames)
  }
  
  if(FECWT){
    freq <- rep(1, nrow(tdata))
    dseq <- c(seq(0, 100, by = 10), 1000)
    for(j in 1:nspec){
      wj <- which(tdata$species == specNames[j])
      jtab <- table( cut(tdata$diam[wj], dseq) )
      jtab <- jtab/sum(jtab)
      dtab <- findInterval(tdata$diam[wj], dseq)
      freq[wj] <- jtab[dtab]
    }
    freq[ is.na(freq) ] <- .1
    freq[freq < .02] <- .02
 #   tdata$fecWt <- sqrt(1/freq)      # should be = freq
    tdata$fecWt <- 1/freq
  }
  
  ngroup <- length(yeGr)
  
  fitCols <- 1:Qfec
  if(length(notCols) > 0)fitCols <- fitCols[-notCols]
  
  
  #prior bgRep
  rVPI <- diag(10, nrow(bgRep))
  rvp  <- bgRep*0
  rvp[1:nspec] <- -.5
  rvp[ endsWith(rownames(bgRep),':diam') ] <- 1
  
  
  #######################
  
  .updateBeta <- .wrapperBeta(rvp, rVPI, priorB, priorIVB, SAMPR, obsRows,
                              tdata, xfecCols, xrepCols, last0first1, ntree, nyr, 
                              betaPrior, years, YR, AR, yrIndex,
                              RANDOM, reIndex, xrandCols, RANDYR,
                              fitCols, specNames, FECWT)
  
  if(SEEDDATA){
    .updateU <- .wrapperU(distall, tdata, minU, maxU, priorU, priorVU,
                          seedNames, nspec, trapRows, obsRowSeed, obsYr, 
                          tau1, tau2, SAMPR, RANDYR, USPEC)
  }
  
  predYr  <- sort( unique(tdata$year) )
  
  tcols <- c('specPlot','species','dcol','year','plotYr','plotyr','obs',
             'fecMin','fecMax')
  if(AR)tcols <- c(tcols,'times')
  if(CONES)tcols <- c(tcols,'cropCount', 'cropFraction', 'cropFractionSd')
  if('cropMin' %in% colnames(tdata))tcols <- c(tcols,'cropMin')
  if('cropMax' %in% colnames(tdata))tcols <- c(tcols,'cropMax')

  updateProp <- c( 1:1000, seq(1001, 10000, by=100) )
  updateProp <- updateProp[updateProp < .9*ng]
  
  pHMC <- .03
  if(nrow(tdata) > 50000) pHMC <- 0
  
  .updateFecundity <- .wrapperStates( SAMPR, USPEC, RANDOM, SEEDDATA, obsTimes, 
                                      plotYears, sdata, tdat = tdata[,tcols], seedNames,
                                      last0first1, distall, YR, AR, trapRows,
                                      obsRows, obsTrapRows, obsYr, predYr, obsRowSeed,
                                      ntree, years, nyr, xrandCols, reIndex, 
                                      yrIndex, plag, groupByInd, RANDYR, updateProp,
                                      seedTraits, pHMC)
 
  ######################## already done
  #tmp <- .unstandBeta(formula = formulaFec, xdata = tdata, xnow = xfec, 
  #                    xmean = xmean, notCols = notCols, specNames = specNames)
  #xfecU   <- tmp$x          
  #xfecs2u <- tmp$s2u      # multiply by beta_s to get beta_u
  #xfecu2s <- tmp$u2s      # multiply by beta_u to get beta_s
  
  
  ########################
  
  ikeep <- 1:ng
  if(ng < keepIter)keepIter <- ng
  if(keepIter < ng){
    ikeep <- round(seq(1, ng, length.out=keepIter))
    ikeep <- ikeep[!duplicated(ikeep)]
  }
  nkeep <- length(ikeep)
  
  bfgibbs  <- matrix(0, nkeep, Qfec); colnames(bfgibbs) <- xFecNames #unstandardized
  brgibbs  <- matrix(0, nkeep, Qrep); colnames(brgibbs) <- xRepNames
  bygibbsF <- bygibbsR <- NULL
  bsgibbs  <- bfgibbs #standardized--prediction from unstandardized may not work
  
  if(SEEDDATA){ 
    minU <- minU[specNames]
    maxU <- maxU[specNames]
    ug <- ug[specNames]
    
    
    ugibbs <- matrix(0, nkeep, nspec)
    colnames(ugibbs) <- specNames
    if(!USPEC) ugibbs <- ugibbs[,1,drop = FALSE]
    
    if(USPEC){
      priorUgibbs <- matrix(0,nkeep,2)
      colnames(priorUgibbs) <- c('mean','var')
    }
    ug[1:nspec] <- .tnorm( nspec, minU, maxU, ug, 5) 
    if(!USPEC)ug[1:nspec] <- ug[1]
  }
  
  sgibbs <- matrix(0, nkeep, 3)
  colnames(sgibbs) <- c('sigma','rmspe','deviance')
  
  ncols <- nyr
  if(AR){
    ncols <- plag
    cnames <- paste('lag',c(1:plag),sep='-')
  }
  
  betaYrF  <- betaYrR <- matrix(0, 1, ncols)
  if( YR ) sgYr <- rep(1, ncols)
  if(YR | AR){
    betaYrF  <- matrix(0, 1, ncols)
    betaYrR  <- matrix(0, ngroup, ncols)
    rownames(betaYrR) <- yeGr
    colnames(betaYrF) <- colnames(betaYrR) <- cnames
    
    bygibbsF <- matrix(NA, nkeep, length(betaYrF))
    bygibbsR <- matrix(0, nkeep, length(betaYrR))
    colnames(bygibbsF) <- colnames(betaYrF)
    colnames(bygibbsR) <- .multivarChainNames(yeGr, colnames(betaYrR))
    bygibbsN <- bygibbsR
  }

  if(AR){
    if(plag == 1){
      Gmat <- matrix(0,1,1)
    }else{
      Gmat  <- rbind(0, cbind( diag(plag-1), 0 ) )
    }
    eigenMat <- eigen1 <- eigen2 <- betaYrR*0
  }
  
  if(RANDOM){
    agibbs <- matrix(NA, nkeep, length(Arand))
    colnames(agibbs) <- .multivarChainNames(xFecNames[xrandCols],
                                            xFecNames[xrandCols])
    aUgibbs <- agibbs
    asum <- asum2 <- aUsum <- aUsum2 <- alphaRand*0
  }
  
  colnames(brgibbs) <- xRepNames
  
  if(SAMPR){
    rgibbs <- matrix(0, nkeep, length(R))
    colnames(rgibbs) <- .multivarChainNames(rownames(R), colnames(R) )
    rgibbs <- rgibbs[,posR]
  }
  
  accept <- rep(0, length(plotYears))
  
  
  pars  <- list(fg = fg,# fecMinCurrent = tdata$fecMin, 
              #  fecMaxCurrent = tdata$fecMax, 
                sg = sg, bgFec = bgFec, bgRep = bgRep,
                betaYrR = betaYrR*0, betaYrF = betaYrF, alphaRand = alphaRand, 
                Arand = Arand)
  if(SEEDDATA){
    pars$ug    <- ug
    pars$umean <- umean
    pars$uvar  <- uvar
    pars$R     <- R
  }
  
  mufec <- xfec%*%bgFec
  muyr  <- muran <- mufec*0
  
  # draw from probit
  tdata$repMu[!is.finite(tdata$repMu)] <- .5
  wlo <- rep(-Inf, length(z))
  whi <- rep(Inf, length(z))
  whi[z == 0] <- 0
  wlo[z == 1] <- 0
  w <- .tnorm(length(z), wlo, whi, tdata$repMu, 1)
  
  tmp <- .updateBeta(pars, xfec, xrep, w, z, zmat, matYr, muyr)
  bgFec <- pars$bgFec <- tmp$bgFec
  bgRep <- pars$bgRep <- tmp$bgRep
  
  if(ARSETUP){
    bgFec[rownames(parameters$betaFec),1] <- parameters$betaFec[,1]
    bgRep[rownames(parameters$betaRep),1] <- parameters$betaRep[,1]
    sg <- parameters$sigma[1,1]
    
    pars$bgFec <- bgFec
    pars$bgRep <- bgRep
    pars$sg <- sg
  }
  
  if(SEEDDATA){
    tmp <- .updateU(pars, z, propU, sdata)
    ug    <- pars$ug    <- tmp$ug
    umean <- pars$umean <- tmp$umean
    uvar  <- pars$uvar  <- tmp$uvar
    propU <- tmp$propU
    if(!USPEC)ug[1:nspec] <- pars$ug <- ug[1]
  }
  
  # tree correlation
  nSpecPlot <- max(yrIndex[,'specPlot'])
  
  fmat <- matrix(0,ntree,nyr)
  treeID <- unique(tdata$treeID)
  ntree  <- length(treeID)
  fmat <- matrix(0,ntree,nyr)
  yrIndex[,'tnum'] <- match(tdata$treeID, treeID)
  rownames(fmat) <- treeID
  colnames(fmat) <- years
  
  nspec  <- length(specNames)
  omegaE <- omegaN <- matrix(0,ntree,ntree)
  
  if(PREDSEED){  #predict species, not seedNames
    specPredSum <- specPredSum2 <- matrix(0, nseedPred, nspec)
    colnames(specPredSum) <- specNames
  }
  
  gupdate <- c(40, 80, 160, 240, 600, 800, 1200, 2000)  
  g1      <- 1
  yupdate <- sort( sample(burnin:ng, 50, replace  = TRUE) )
  yupdate <- unique(yupdate)
  pupdate <- burnin:ng
  if(length(pupdate) > 100)pupdate <- 
    unique( round( seq(burnin, ng, length=100) ))
  
  if(SEEDDATA){
    svarEst <- rep(0, nrow(sdata))
    names(svarEst) <- rownames(sdata)
    svarEst <- svarPred <- svarEst2 <- svarPred2 <- svarEst
    specSum <- specSum2  <- matrix(0,nrow(sdata),nspec)
    colnames(specSum)    <- colnames(specSum2) <- specNames
    rownames(specSum)    <- rownames(specSum2) <- rownames(sdata)
    activeArea <- sdata$area
    
    if(SEEDCENSOR){ #locations of censored seed counts
      censMin <- as.matrix(censMin)
      censMax <- as.matrix(censMax)
      censIndex <- which(censMax > censMin, arr.ind  = TRUE)
      censIndex <- sort(unique(censIndex[,1]))
      areaCens <- sdata[rownames(censMin),'area']
      cens2sdata <- match(rownames(censMin),rownames(sdata))
    }
  }
  
  
  ntoty  <- ntotyy <- rmspe <- deviance <- 0
  ntot   <- 0
  zest   <- zpred <- fest <- fest2 <- fpred <- fpred2 <- fg*0 # fecundity 
  sumDev <- ndev <- 0   #for DIC
  
  nPlotYr    <- max(tdata$plotyr)
  acceptRate <- nPlotYr/5
  
  ########### growth rate for maturation transition
  
  grow <-  smat <- rmat <- matrix( NA, ntree, nyr)
  grow[ yrIndex[,c('tnum','year')] ] <- tdata$diam
  gdev <- sweep(grow, 1, rowMeans(grow, na.rm  = TRUE),'-')
  tdev <- c(1:nyr) - (1 + nyr)/2
  tdev <- matrix(tdev, ntree, nyr, byrow  = TRUE)
  cgt  <- rowMeans( gdev * tdev, na.rm  = TRUE )   # covariance growth/time
  tvr  <- apply(tdev, 1, var, na.rm  = TRUE)
  slp  <- cgt/tvr                            # mean growth ratre
  slp[slp < .005] <- .005
  growVec <- slp[ tdata$tnum ]
  
  pbar <- txtProgressBar(min=1,max=ng,style=1)
  
  if(SEEDDATA)logScoreStates <- logScoreFull <- seedCount*0
  
  nprob <- 0
  
  fg[fg > tdata$fecMax] <- tdata$fecMax[fg > tdata$fecMax]
  fg[fg < tdata$fecMin] <- tdata$fecMin[fg < tdata$fecMin]
  pars$fg <- fg
  
  # trees with multiple years
  w2 <- table(tdata$treeID)
  multYear <- names(w2)[w2 > 1]
  MULTYR   <- tdata$treeID %in% multYear
  
  epsilon <- (log(fg) + log(tdata$fecMax + 1.01 - fg))/100000
  
  gk <- 0
  
  ##########################
#  fg    <- pars$fg <- trueValues$fec[rownames(tdata)]
 # z     <- trueValues$repr[rownames(tdata)]
  ########################
  
  if(RANDOM){
    minmax <- 4
    amu    <- rep(0, nspec)
    names(amu) <- colnames(xfec)[xrandCols]
    if(nspec == 1 & length(xrandCols) == 1)names(amu) <- specNames
    ispec  <- tdata$species[ match(treeID, tdata$treeID) ] # species column, assumes only random intercepts
    ispec  <- match( ispec, specNames )
    imat   <- matrix(0, ntree, nspec)
    imat[ cbind(1:ntree, ispec) ] <- 1
    rownames(imat) <- treeID
    colnames(imat) <- specNames
  }
  
  # individual trends
  sfmat  <- matrix(NA, ntree, nyr)
  rownames(sfmat) <- treeID
  symat <- matrix(1:ncol(sfmat), nrow(sfmat),ncol(sfmat), byrow=T)
  #sgmat  <- matrix(NA, ntree, nyr)
  rownames(symat) <- treeID
  slopesRate <- slopesNyr <- rep(0, ntree)
  
  
  for(g in 1:ng){ ####################
    
    
    if(g %in% ikeep)gk <- gk + 1    
    if(gk > ng)break
    
    pars$fg <- fg
    yg      <- log(fg)      
    mufec   <- xfec%*%bgFec
    
    if(RANDOM){
      
      yg <- yg - mufec
      if(YR){
        if(RANDYR){
          yg <- yg - betaYrR[yrIndex[,c('group','year')]] 
        }else{
          yg <- yg - betaYrF[yrIndex[,'year']] 
        }
      }
      if(AR)yg <- yg - muyr
      
      tt <- table(reIndex[z == 1])                # only mature individuals
      reGroups <- as.numeric(names(tt)[tt > 2])   # multiyear
      wrow     <- which(reIndex %in% reGroups)
      krow     <- which(!reIndex %in% reGroups)

      tmp <- .updateAlphaRand(ntree, yA = yg[wrow,], 
                              xfecA = xfec[wrow,xrandCols,drop = FALSE], 
                              sg, reIndexA = reIndex[wrow], reGroups,
                              Arand, priorVA, dfA, specNames, minmax = 1.5)
      Arand       <- pars$Arand <- tmp$Arand
      alphaRand   <- pars$alphaRand <- tmp$alphaRand
      meanRand    <- tmp$meanRand                      # mean RE by species
      
      if(g >= burnin)alphaRandU  <- tmp$alphaRand%*%t(xrands2u)
      
      
      if(g %in% ikeep){
        ArandU       <- xrands2u%*%tmp$Arand%*%t(xrands2u)
        agibbs[gk,]  <- as.vector(Arand)   
        aUgibbs[gk,] <- as.vector(ArandU)
      }

      alphaRand  <- alphaRand*imat
      if(g >= burnin)alphaRandU <- alphaRandU*imat
      
      muran <- alphaRand[reIndex,]
      if(length(Arand) > 1)muran <- rowSums( muran )
      
      
      # single year - use if only random intercepts, change if there are random slopes
      
      if( length(krow) > 0 ){
        amu[ names(meanRand) ] <- meanRand
        imu <- amu[yrIndex[krow,'spp']] 
        asd <- sqrt( diag(Arand)[yrIndex[krow,'spp']] )
        muran[krow] <- .tnorm( length(krow), -minmax, minmax, imu, asd)
      }
    }
    
    if(YR){
      yg <- log(fg) - mufec
      if(RANDOM)yg <- yg - muran
      
      tmp      <- .updateBetaYr(yg, z, sg, sgYr, betaYrF, betaYrR, yrIndex, yeGr,
                                RANDYR, obs = tdata$obs)
      betaYrF  <- pars$betaYrF <- tmp$betaYrF
      betaYrR  <- pars$betaYrR <- tmp$betaYrR
      sgYr     <- pars$sgYr    <- tmp$sgYr
      wfinite  <- tmp$wfinite
      
      if(g %in% ikeep){
        bygibbsF[gk,] <- betaYrF
        bygibbsR[gk,] <- betaYrR
        bygibbsN[gk,wfinite] <-  1
        sygibbs[gk,]  <- sgYr
      }
      muyr <- betaYrF[yrIndex[,'year']] 
      if(RANDYR) muyr <- muyr + betaYrR[ yrIndex[,c('group','year')] ]
    }
    
    if(AR){
      yg <- log(fg) 
      mu <- mufec
      if(RANDOM)mu <- mu + muran
      
      tmp <- .updateBetaAR_RE(betaYrF, betaYrR, Alag, yg, mu, z, 
                              lagGroup, lagMatrix, plag, ngroup, sg)
      betaYrF <- pars$betaYrF <- tmp$betaYrF
      betaYrR <- pars$betaYrR <- tmp$betaYrR
      wfinite <- which(betaYrR != 0)
      Alag    <- tmp$Alag
      muyr <- tmp$ylag
      if(g %in% ikeep){
        bygibbsF[gk,] <- betaYrF
        bygibbsR[gk,] <- betaYrR
        bygibbsN[gk,wfinite] <-  1
      }
    }
    
    wlo   <- 10*(z - 1)
    whi   <- 10*z
    w     <- .tnorm(length(z), wlo, whi, xrep%*%bgRep, 1)
    tmp   <- .updateBeta(pars, xfec, xrep, w, z, zmat, matYr, muyr)
    bgFec <- pars$bgFec <- tmp$bgFec
    bgRep <- pars$bgRep <- tmp$bgRep
    
    if(SEEDDATA){
      tmp   <- .updateU(pars, z, propU, sdata)
      ug    <- pars$ug    <- tmp$ug
      umean <- pars$umean <- tmp$umean
      uvar  <- pars$uvar  <- tmp$uvar
      propU <- tmp$propU
    }
    
    tmp <- .updateFecundity(g, pars, xfec, xrep, propF, z, zmat, matYr, muyr,
                            epsilon = epsilon)
    fg    <- pars$fg <- tmp$fg
    fecMinCurrent <- tmp$fecMinCurrent
    fecMaxCurrent <- tmp$fecMaxCurrent
    z     <- tmp$z
    zmat  <- tmp$zmat
    matYr <- tmp$matYr
    propF <- tmp$propF
    epsilon <- tmp$epsilon
    
    
    muf <- xfec%*%bgFec
    if(YR | AR)muf <- muf + muyr
    if(RANDOM)muf  <- muf + muran
    
    wrow   <- obsRows[z[obsRows] == 1]
    
    sg <- pars$sg <- .updateVariance(log(fg[wrow]), muf[wrow], s1, s2)
    
    if(sg > 30)sg <- pars$sg <- 30
    
    if(SAMPR){
      tmp  <- .updateR(ug, fz = fg[obsTrapRows]*z[obsTrapRows], SAMPR, USPEC, distall, 
                       sdata, seedNames, 
                       tdat = tdata[obsTrapRows,c('specPlot','year','plotyr','dcol')], 
                       R, priorR, priorRwt, obsYr, posR, plots)
      pars$R <- R <- tmp
      if(g %in% ikeep)rgibbs[gk,]  <- R[posR]
    }
     
    
    
    
    ################## predicted state
    
    brep <- (xreps2u%*%bgRep)[-c(1:nspec),]  #unstandardize: d log(S)/dD = beta_diam
    names(brep) <- .replaceString(names(brep),'species','')
    names(brep) <- .replaceString(names(brep),':diam','')
    brep[brep < 0] <- 0
    brep <- brep[tdata$species]
    
    dSdt <- brep*growVec   # d log(S)/dt = d log(S)/dD dD/dt
    
    tvec <- pnorm(xrep%*%bgRep)
    tii  <- tapply(tvec, tdata$tnum, mean, na.rm  = TRUE)
    tii  <- matrix(tii, ntree, nyr)
    
    smat[ yrIndex[,c('tnum','year')] ] <- tvec
    smat[ is.na(smat) ] <- tii[ is.na(smat) ]
    smat[,-1] <- 1 - smat[,1]
    
    rmat[ yrIndex[,c('tnum','year')] ] <- dSdt
    rmat[,1] <- 0
    rmat[is.na(rmat)] <- 0
    rcum <- t( apply(1 - rmat, 1, cumprod) )
    rmat <- rmat*rcum*smat
    rmat[,1] <- smat[,1]
    
    rmat <- cbind(rmat, 1 - rowSums(rmat))
    
    rmat <- myrmultinom(1, rmat)[,-ncol(rmat)]
    rmat[ !is.finite(rmat) ] <- 0
    rmat <- t( apply(rmat,1,cumsum) )
    zw   <- rmat[ yrIndex[,c('tnum','year')] ]
    
    
    flo <- fhi <- zw*0
    flo[zw == 0] <- -9.21034
    fhi[zw == 1] <- log(fecMaxCurrent[zw == 1] + .1)
    
    ymu <- .tnorm( length(muf), flo, fhi, muf, sqrt(sg) )
    fmu <- exp( ymu )                          # mean prediction
 
    # save unstandardized
    bfSave <- bgFec
    brSave <- bgRep
    
    if(UNSTAND){
      if(length(xfecs2u) > 0)bfSave[fitCols,1] <- xfecs2u%*%bgFec[fitCols,1]  
      if(length(xreps2u) > 0)brSave <- xreps2u%*%bgRep
      if(length(notFit) > 0)bfSave[notFit,1] <- 0
    }
    
    if(g %in% ikeep){
      bfgibbs[gk,] <- bfSave
      brgibbs[gk,] <- brSave
      bsgibbs[gk,] <- bgFec
      if(SEEDDATA){
        if(USPEC){
          ugibbs[gk,]  <- ug
          priorUgibbs[gk,] <- c(umean,uvar)
        }else{
          ugibbs[gk,1] <- ug[1]
        }
      }
    }
    
    if(g %in% gupdate & gk > 10 & SEEDDATA){
      gi <- (gk - 20):gk
      uu <- ug/4
      gi <- gi[gi > 0]
      if(USPEC){
        propU  <- apply(ugibbs[gi,],2,sd) + .01
        propU[propU > uu] <- uu[propU > uu]
      }else{
        propU <- sd(ugibbs[gi,]) + .01
        propU <- min(c(propU,uu))
      }
    }
    
    
    
    
    if(SEEDCENSOR){                # missing censored traps
      lf <- .getLambda(tdata[obsTrapRows,c('specPlot','year','plotyr','dcol')],
                       sdata[rownames(censMin),c('year','plotyr','drow')],
                       areaCens, ug, fg[obsTrapRows]*z[obsTrapRows], R, 
                       SAMPR, USPEC, distall, obsYr, PERAREA = FALSE)  # per trap
      lf[lf < 1e-9] <- 1e-9
      
      ttt   <- rtpois( lo = censMin[,seedNames,drop=F], 
                       hi = censMax[,seedNames,drop=F], mu = lf )
      sdata[rownames(censMin), colnames(ttt) ] <- ttt
    }
    
    if(g %in% pupdate & SEEDDATA){
        
      nprob <- nprob + 1
      
      # estimated fecundity per m2
      fz <- fg[obsTrapRows]*z[obsTrapRows]  
      lm <- .getLambda(tdata[obsTrapRows,c('specPlot','year','plotyr','dcol')],
                       sdata[,c('year','plotyr','drow')],
                       AA=1, ug, fz, R, 
                       SAMPR, USPEC, distall, obsYr, PERAREA  = TRUE, SPECPRED  = TRUE)  
      lm[lm < 1e-9] <- 1e-9
      pm    <- matrix( rpois(length(lm), lm), nrow(lm), ncol(lm) )
      specSum  <- specSum + pm
      specSum2 <- specSum2 + pm^2
      
      # estimated seeds per trap
      ls <- lm*sdata$area                      
      pf <- matrix( rpois(length(ls), ls), nrow(ls), ncol(ls) )
  #    SvarEst[nprob,] <- rowSums(pf)
      spf      <- rowSums(pf)
      svarEst  <- svarEst + spf
      svarEst2 <- svarEst2 + spf^2
      

      # predicted per trap
      la    <- .getLambda(tdata[obsTrapRows,c('specPlot','year','plotyr','dcol')], 
                          sdata[,c('year','plotyr','drow')],
                          AA = activeArea, ug, fmu[obsTrapRows]*zw[obsTrapRows], 
                          R, SAMPR, USPEC, distall, obsYr, PERAREA = FALSE) 
      la[la < 1e-9] <- 1e-9
      pg    <- matrix( rpois(length(la), la), nrow(la), ncol(la) )
      
      resid <- (seedCount - pg)^2
      if(SEEDCENSOR)resid[cens2sdata[censIndex],] <- 0
      
      rmspe <- sqrt( mean( resid, na.rm  = TRUE ) )
      
  #    SvarPred[nprob,] <- rowSums(pg)
      spf      <- rowSums(pg)
      svarPred  <- svarPred + spf
      svarPred2 <- svarPred2 + spf^2
      
      # deviance from predicted fecundity
      dev   <- dpois(seedCount, la, log  = TRUE)
      
      

      if(SEEDCENSOR){
        mm <- match(rownames(censMin),rownames(sdata)[obsRowSeed])
        
        pr <- dtpois(censMin[,seedNames,drop=F], censMax[,seedNames,drop=F], 
                     la[cens2sdata,], 
                     index = censIndex) 
        dev[cens2sdata,] <- log( pr )
      }
      deviance <- sum(dev)
      sumDev <- sumDev - 2*deviance
      ndev <- ndev + 1
      
      # score from predicted fecundity
      logScoreStates <- logScoreStates - dpois(seedCount, la, log  = TRUE)
      logScoreFull   <- logScoreFull - dev
      
      if(PREDSEED){  #NOTE: from estimated, not predicted fg
        ls <- .getLambda(tdatPred[,c('specPlot','year','plotyr','dcol')], 
                         sdatPred[,c('year','plotyr','drow')], AA = 1,
                         ug, fz[tdatPred[,'row']], R, 
                         SAMPR, USPEC, distPred, predList$years, PERAREA  = TRUE,
                         SPECPRED  = TRUE)   # per m^2
        ls <- ls + 1e-9
        ps <- matrix( rpois(nseedPred*ncol(ls), ls), nseedPred, nspec )
        specPredSum  <- specPredSum + ps
        specPredSum2 <- specPredSum2 + ps^2
      }
    }
    
    if(g %in% ikeep)sgibbs[gk,]  <- c(sg, rmspe, deviance)
    
    if(g %in% yupdate){
      
      ntoty  <- ntoty + 1
      
      # individual rate of change
      
      fzm <- fg*z
      
      sfmat[ yrIndex[,c('tnum','year')] ] <- log(fzm)
      sfmat[ sfmat < 2 ] <- NA
      sgmat <- sfmat*0 + 1
      
      # for rates
      wmat <- which( rowSums(sgmat, na.rm=T) > 1 )
      
      if(length(wmat) > 1){
        
        ntotyy <- ntotyy + 1
        
        symat[is.na(sfmat)] <- NA
        fmm  <- rowMeans(sfmat[wmat,], na.rm=T)
        ymm  <- rowMeans(symat[wmat,], na.rm=T)
        sfmat[wmat,] <- sweep(sfmat[wmat,], 1, fmm,'-')
        symat[wmat,] <- sweep(symat[wmat,], 1, ymm,'-')
        snmat <- rowSums(symat[wmat,]*0 + 1, na.rm=T)
        cvv   <- rowSums(sfmat[wmat,]*symat[wmat,],na.rm=T)
        vvv   <- rowSums(symat[wmat,]^2, na.rm=T)
        slope <- cvv/vvv
        wf    <- which(is.finite(slope))
        
        slopesRate[ wmat[wf] ] <- slopesRate[ wmat[wf] ] + slope[wf]
        slopesNyr[ wmat[wf] ]  <- slopesNyr[ wmat[wf] ] + snmat[wf]
      }
      
      
      ########################
      
      # by group
      if(AR){
        for(j in 1:ngroup){
          Gmat[1,] <- betaYrF + betaYrR[j,]
          eigenMat[j,] <- eigen(Gmat)$values
        }
        eigen1 <- eigen1 + eigenMat
        eigen2 <- eigen2 + eigenMat^2
      }
      
      fecRes <- fg
      yRes <- yg
      yRes <- yRes - mean(yRes)           #tree autocorrelation
      
      for(m in 1:nSpecPlot){############### MATCH TREEID ROWS FOR FMAT
        
        fmat <- fmat*0
        
        wm   <- which(yrIndex[,'specPlot'] == m & tdata$obs == 1 & fg > 1)
        dm   <- tdata$tnum[wm]
        ym   <- yrIndex[wm,'year']
        dr   <- unique(dm)
        if(length(dr) < 2)next
        
        fmat[ cbind(dm,ym) ] <- fecRes[wm]
        fm  <- fmat[unique(dm),unique(ym)]
        
        # between trees
        ff <- suppressWarnings( cor(t(fm), use="pairwise.complete.obs") ) 
        wf <- which(is.finite(ff))
        
        omegaE[dr,dr][wf] <- omegaE[dr,dr][wf] + ff[wf]
        omegaN[dr,dr][wf] <- omegaN[dr,dr][wf] + 1
        
        acm <- acfEmp(yRes[wm], tdata$tnum[wm], yrIndex[wm,'year'])
        wfin <- which(is.finite(acm))
        wfin <- wfin[wfin <= npacf]
        if(length(wfin) == 0)next
        
        pa <- try( pacf(acm[wfin]), T )
        
        if( !inherits(pa,'try-error') ){
          pacfMat[m,wfin] <- pacfMat[m,wfin] + pa
          pacf2[m,wfin]   <- pacf2[m,wfin] + pa^2
          acfMat[m,wfin] <- acm[wfin]
          pacN[m,wfin] <- pacN[m,wfin] + 1
        }
      }
    }
     
    if(g >= burnin){
      
      ntot <- ntot + 1
      
      zest  <- zest + z
      fz    <- fg*z
      fest  <- fest + fz  # only add to mature 
      fest2 <- fest2 + fz^2
      
      zpred  <- zpred + zw
      fz     <- fmu*zw
      fpred  <- fpred + fz
      fpred2 <- fpred2 + fz^2
      
      if(RANDOM){
        asum   <- asum + alphaRand
        asum2  <- asum2 + alphaRand^2
        aUsum  <- aUsum + alphaRandU
        aUsum2 <- aUsum2 + alphaRandU^2
      }
    }
    setTxtProgressBar(pbar,g)
  } ###########################################################
       
   # to re-initialize
  tdata$lastFec  <- fg
  tdata$lastRepr <- z
  
  # fecundity
  matrEst  <- zest/ntot
  matrPred <- zpred/ntot
  
  fecEstMu <- fest/ntot     #  same order at tdata
  fecEstSe <- fest2/ntot - fecEstMu^2
  fecEstSe[fecEstSe < 0] <- 0
  fecEstSe <- sqrt(fecEstSe)
  
  fecPredMu <- fpred/ntot     #  pred|z = 1
  fecPredSe <- fpred2/ntot - fecPredMu^2
  fecPredSe[fecPredSe < 0] <- 0
  fecPredSe <- sqrt(fecPredSe)
  
  fecPred <- tdata[,c('plot','treeID','species','year','diam','dcol')]
  if(YR | AR){
    ygr <- yrIndex
    colnames(ygr) <- paste('ind_',colnames(ygr),sep='')
    fecPred <- cbind(fecPred, ygr)
  }
  if(CONES){
    cropCount <- seedTraits[tdata$species,'seedsPerFruit']*tdata$cropCount/tdata$cropFraction
    fecPred <- cbind(fecPred, cropCount)
  }
  mpm <- round( cbind(matrEst, matrPred), 3)
  colnames(mpm) <- c('matrEst','matrPred')
  
  fpm <- cbind(fecEstMu, fecEstSe, fecPredMu, fecPredSe)
  fpm[is.na(fpm)] <- 0
  fpm <- round(fpm,1)
  colnames(fpm) <- c('fecEstMu', 'fecEstSe', 'fecPredMu', 'fecPredSe')
  
  fecPred <- cbind(fecPred, mpm, fpm)
  
  if(CONES){
    rmspeCrop <- sqrt( mean( (fecPredMu - cropCount)^2, na.rm=T ) )
  }
  
  if(SEEDDATA){
    
    scols <- c('plot','year','trapID','drow','area')
    countPerTrap <- rowSums(sdata[,seedNames,drop = FALSE])
    
 #   seedEst <- cbind( colMeans(SvarEst), apply(SvarEst, 2, sd) )
    seedEst <- svarEst/nprob
    s2      <- svarEst2/nprob - seedEst^2
    seedEst <- cbind(seedEst, sqrt(s2))
    colnames(seedEst) <- c('estMeanTrap', 'estSeTrap')
    
  #  seedPred <- cbind( colMeans(SvarPred), apply(SvarPred, 2, sd) )
    seedPred <- svarPred/nprob
    s2       <- svarPred2/nprob - seedPred^2
    seedPred <- cbind(seedPred, sqrt(s2))
    colnames(seedPred) <- c('predMeanTrap', 'predSeTrap')
    
    svv <- (countPerTrap - seedPred[,'predMeanTrap'])^2 
    predSeError  <- signif( sqrt(svv) , 3)
    
    seedSpecMu <- specSum/nprob
    seedSpecSe <- specSum2/nprob - seedSpecMu^2
    seedSpecSe[seedSpecSe < 0] <- 0
    seedSpecSe <- sqrt(seedSpecSe)
    
    seedPred <- data.frame( cbind( sdata[,scols], countPerTrap,
                                   signif(seedEst, 3),
                                   signif(seedPred, 3),
                                   predSeError),
                            stringsAsFactors = F)
    m1 <- paste(specNames,'meanM2',sep='_')
    m2 <- paste(specNames,'sdM2',sep='_')
    m2Mu <- seedSpecMu[,specNames,drop = FALSE]
    m2Se <- sqrt( seedSpecSe[,specNames,drop = FALSE]^2 ) 
    colnames(m2Mu) <- m1
    colnames(m2Se) <- m2
    seedPred <- cbind(seedPred, signif(m2Mu, 3), signif(m2Se,3))
    
    inflation <- signif( predSeError/(.1 + seedPred$predSeTrap), 3)
    
    pvv <- seedPred$predSeTrap^2  # predictive variance
    
    seedPred <- cbind(seedPred, inflation)
    
    #entire data set
    
    meanPredErrSd <- mean(predSeError, na.rm  = TRUE)
    meanPredSd    <- mean(seedPred$predSeTrap, na.rm  = TRUE)
    meanInflation <- mean(inflation, na.rm  = TRUE)
    
    # fit  
    
  #  nss <- length(obsRowSeed)
    
    MM <- FALSE
    if(!all(seedTraits[,'gmPerSeed'] == 1))MM <- TRUE
    
    if( MM ){
      mss <- matrix(seedTraits[specNames,'gmPerSeed'],nrow(sdata),nspec,byrow  = TRUE)
      massMu <- seedSpecMu[,specNames,drop = FALSE]*mss
      massSe <- sqrt( seedSpecSe[,specNames,drop = FALSE]^2*mss^2 ) 
      m1 <- paste(colnames(massMu),'meanGmM2',sep='_')
      m2 <- paste(colnames(massSe),'sdGmM2',sep='_')
      colnames(massMu) <- m1
      colnames(massSe) <- m2
      seedPred <- cbind(seedPred, signif(massMu, 3), signif(massSe,3))
    }
    
    if(PREDSEED){
      scols <- c('plot','trapID','year','x','y','drow','dgrid','area','active')
      specMu <- specPredSum/nprob
      sse <- specPredSum2/nprob - specMu^2
      specSe <- sqrt(sse)
      colnames(specMu) <- paste(colnames(specMu),'_meanM2',sep='')
      colnames(specSe) <- paste(colnames(specSe),'_sdM2',sep='')
      
      nss <- nrow(sdatPred)
      mss <- matrix(seedTraits[specNames,'gmPerSeed'],nss,nspec,byrow  = TRUE)
      massMu <- specMu*mss
      massSe <- sqrt( specSe^2*mss^2 ) 
      m1 <- paste(specNames,'meanGmM2',sep='_')
      m2 <- paste(specNames,'sdGmM2',sep='_')
      colnames(massMu) <- m1
      colnames(massSe) <- m2
      
      preds <- signif(cbind(specMu, specSe, massMu, massSe), 3)
      
      seedPredGrid <- data.frame( cbind( sdatPred[,scols], preds ) )
      treePredGrid <- cbind(tdatPred, fecPred[tdatPred$row,])
      
      # out-of-sample
      if(!is.null(modelYears)){
        
        sdataOut$plotTrapYr <- columnPaste(sdataOut$trapID,sdataOut$year)
        sdataOut$fore <- sdataOut$year - max(modelYears)
        
        tdataOut$plotTreeYr <- columnPaste(tdataOut$treeID,tdataOut$year)
        treePredGrid$plotTreeYr <- columnPaste(treePredGrid$treeID,treePredGrid$year)
        
        fec <- matrix(NA,nrow(tdataOut),4)
        colnames(fec) <- colnames(fpm)
        ww <- which(tdataOut$plotTreeYr %in% treePredGrid$plotTreeYr)
        qq <- match(tdataOut$plotTreeYr[ww], treePredGrid$plotTreeYr)
        fec[ww,] <- fpm[qq,]
        tdataOut <- cbind(tdataOut,fec)
      }
    }
    
    # seed, fecundity acf
    seedRes <- t( t(seedCount) - colMeans(seedCount, na.rm  = TRUE) )
    ww <- colSums(seedCount, na.rm  = TRUE)
    
    pacsMat <- NULL
    
    if(sum(ww) > 0){
      seedRes <- seedRes[,ww > 0,drop = FALSE]
      
      ii <- rep(sdata$drow,ncol(seedRes))
      yy <- match(sdata$year,years)
      jj <- rep(yy, ncol(seedRes))
      kk <- is.finite(ii)
      ii <- ii[kk]
      jj <- jj[kk]
      sr <- as.vector(seedRes)[kk]
      
      acm  <- acfEmp(sr, ii, jj)
      wfin <- which(is.finite(acm))
      pacsMat <- acm*0
      pacsMat[wfin] <- pacf(acm[wfin])
    }
  }
  
  kg <- which(ikeep >= burnin & ikeep <= ng)
  
  betaFec <- .chain2tab(bfgibbs[drop = F, kg,]) #unstandardized
  betaRep <- .chain2tab(brgibbs[drop = F, kg,])
  betaStd <- .chain2tab(bsgibbs[drop = F, kg,]) #standardized
  
  acfMat  <- acfMat/ntoty
  pacfMat <- pacfMat/pacN
  pacfSe  <- pacf2/pacN - pacfMat^2
  pacfMat[!is.finite(pacfMat)] <- 0
  pacfSe[!is.finite(pacfSe)] <- 0
  
  wk <- which(rowSums(pacfSe) > 0)
  pacfMat <- pacfMat[wk,,drop=F]
  pacfSe  <- pacfSe[wk,,drop=F]
  
  omegaE <- omegaE/omegaN
  
  ################ individual slopes
  
  trendRate <- slopesRate/ntotyy
  trendNyr  <- slopesNyr/ntotyy
  trendTree <- signif( cbind( trendRate, trendNyr ), 3 )
  colnames(trendTree) <- c('logf/yr', 'years')
  
  # by species-plot
  si <- tdata$species[ match(treeID, tdata$treeID) ]
  ip <- tdata$plot[ match(treeID, tdata$treeID) ]
  ll <- list(plot = ip, species = si)
  
  # weighted by series length
  ws <- tapply( trendRate*trendNyr, ll, sum, na.rm=T)
  wt <- tapply( trendNyr, ll, sum, na.rm=T)
  trendMu <- ws/wt
  
  vt <- tapply( trendRate^2*trendNyr, ll, sum, na.rm=T)
  trendSe <- sqrt( vt/wt - trendMu^2 )
  
  colnames(trendMu) <- paste(colnames(trendMu),'Mean')
  colnames(trendSe) <- paste(colnames(trendSe),'SE')
  
  trendPlotSpec <- signif( cbind(trendMu, trendSe), 4 )
  
  # by species
 
  # weighted by series length
  ws <- tapply( trendRate*trendNyr, si, sum, na.rm=T)
  wt <- tapply( trendNyr, si, sum, na.rm=T)
  trendMu <- ws/wt
  
  vt <- tapply( trendRate^2*trendNyr, si, sum, na.rm=T)
  trendSe <- sqrt( vt/wt - trendMu^2 )
  
  names(trendMu) <- paste(names(trendMu),'Mean')
  names(trendSe) <- paste(names(trendSe),'SE')
  
  trendSpec <- signif( c(trendMu, trendSe), 4 )
  
  trendEst <- list( trendSpec = trendSpec, trendPlotSpec = trendPlotSpec, 
                    trendTree = trendTree )
  
  
  ################ years/lags
  
  if(YR | AR){
    ncol <- nyr
    ccol <- years
    if(AR){
      ncol <- plag
      ccol <- colnames(bygibbsF)
    }
    
    betaYrMu <- betaYrSe <- matrix(nrow(betaYrF), ncol=ncol )
    wg <- which(rowSums(abs(bygibbsF)) != 0)
    wg <- wg[wg %in% kg]
    
   
    if(length(wg) > 0){
      
      betaYr <- .chain2tab(bygibbsF[drop = FALSE, wg,])
      betaYrMu <- matrix( colMeans(bygibbsF[drop = FALSE, wg,], na.rm  = TRUE), 
                          nrow(betaYrF), ncol=ncol )
      betaYrSe <- matrix( apply(bygibbsF[drop = FALSE, wg,], 2, sd, na.rm  = TRUE), 
                          nrow(betaYrF), ncol=ncol )
    }
 
    betaYrRand <- betaYrRandSE <- betaYrMu*0
    
    if(RANDYR){
      betaYrRand <- betaYrRandSE <- betaYrR*0
      wg <- which(rowSums(abs(bygibbsR)) != 0)
      wg <- wg[wg %in% kg]
  
      if(length(wg) > 0){   ###colSums(x, na.rm = na.rm, dims = dims, ...) :  'x' must be an array of at least two dimensions
        brsum <- matrix( colSums(bygibbsR[wg,], na.rm  = TRUE), 
                         nrow(betaYrR), ncol )
        brn <- matrix( colSums(bygibbsN[wg,], na.rm  = TRUE), 
                       nrow(betaYrR), ncol )
        betaYrRand <- brsum/brn
        brn2 <- matrix( colSums(bygibbsR[wg,]^2, na.rm  = TRUE), 
                        nrow(betaYrR), ncol )
        ser <- sqrt( brn2/brn - betaYrRand^2 )
        betaYrRandSE <- ser
        betaYrRand[!is.finite(betaYrRand)] <- 0
      }
      rownames(betaYrRand) <- rownames(betaYrRandSE) <- yeGr
      colnames(betaYrRand) <- colnames(betaYrRandSE) <- ccol
    }
  }
  
  ################ REs
  
  if(RANDOM){
    alphaMu <- asum/ntot
    av    <- asum2/ntot - alphaMu^2
    av[av < 0] <- 0
    alphaSe <- sqrt(av)
    
    alphaUMu <- aUsum/ntot
    av    <- aUsum2/ntot - alphaUMu^2
    av[av < 0] <- 0
    alphaUSe <- sqrt(av)
    
    if(ONEA){
      aMu <- mean(agibbs[kg,])
      aSe <- sd(agibbs[kg,])
      aUMu <- mean(agibbs[kg,])
      aUSe <- sd(agibbs[kg,])
    }else{
      colnames(alphaMu) <- colnames(alphaSe) <- 
        colnames(alphaUMu) <- colnames(alphaUSe) <- xFecNames[xrandCols]
      rownames(alphaMu) <- rownames(alphaSe) <- 
        rownames(alphaUMu) <- rownames(alphaUSe) <- as.character()
      aMu <- matrix( apply(agibbs[drop = FALSE,kg,], 2, mean), Qrand, Qrand )
      aSe <- matrix( apply(agibbs[drop = FALSE,kg,], 2, sd), Qrand, Qrand )
      aUMu <- matrix( apply(aUgibbs[drop = FALSE,kg,], 2, mean), Qrand, Qrand )
      aUSe <- matrix( apply(aUgibbs[drop = FALSE,kg,], 2, sd), Qrand, Qrand )
      colnames(aMu) <- rownames(aMu) <- colnames(aSe) <- rownames(aSe) <-
        colnames(aUMu) <- rownames(aUMu) <- colnames(aUSe) <- rownames(aUSe) <-
        colnames(alphaMu)
      names(xrandCols) <- xFecNames[xrandCols]
    }
  }
  
  mmu <- 1
  if(SAMPR){
    
    tmu <- apply(rgibbs[kg,], 2, mean)
    tse <- apply(rgibbs[kg,], 2, sd)
    
    mmu <- mse <- priorR*0
    mmu[posR] <- tmu
    mse[posR] <- tse
    
    mmu[mmu == 0 & priorR == 1] <- 1
    
    attr(mmu,'posR') <- attr(mse,'posR') <- posR
    
    colnames(mmu) <- colnames(mse) <- seedNames
    rownames(mmu) <- rownames(mse) <- rownames(R)
  }
  
  #################### dispersal
  
  if(SEEDDATA){
    
    dgibbs <- pi*sqrt(ugibbs)/2
    upars <- .chain2tab(ugibbs[kg,,drop = FALSE])[,c(1:4)]
    dpars <- .chain2tab(dgibbs[kg,,drop = FALSE])[,c(1:4)]
    
    if( USPEC ){
      
      uByGroup <- colMeans(ugibbs[drop = FALSE, kg,])
      dByGroup <- colMeans(dgibbs[drop = FALSE, kg,])
      priorDgibbs <- pi*sqrt(priorUgibbs[drop = FALSE, kg,])/2
      
      uall <- .chain2tab(priorUgibbs[kg,,drop = FALSE])[,c(1:4)]
      dall <- .chain2tab(priorDgibbs[,,drop = FALSE])[,c(1:4)]
      
      upars <- rbind(upars,uall)
      dpars <- rbind(dpars,dall)
      
    }else{
      upars <- .chain2tab(ugibbs[kg,,drop = FALSE])[,c(1:4)]
      uByGroup <- mean(ugibbs[kg,])
      
      dpars <- .chain2tab(dgibbs[kg,,drop = FALSE])[,c(1:4)]
      dByGroup <- mean(dgibbs[kg,])
    }
  }
  
  su <- .chain2tab(sgibbs[kg,])[,c(1:4)]
  
  # coefficients are saved unstandardized 
  
  if(UNSTAND){
    xfec <- xfecU
    xrep <- xrepU
  }
  
  zp <- pnorm(xrep%*%betaRep[,1])
  fp <- xfec%*%matrix(betaStd[,1], ncol=1) 
  
  if(YR){
    byr <- betaYrMu[ yrIndex[,'year'] ] 
    if(RANDYR)byr <- byr + betaYrRand[ yrIndex[,c('group','year')] ]
    byr[is.na(byr)] <- 0
    fp <- fp + byr
  }
  
  ############### AR
  
  if(AR){

    eigenMu <- eigen1/ntoty
    eigenSe <- sqrt( eigen2/ntoty - eigenMu^2 )
    
    yg   <- log(fecPred$fecPredMu)      
    yg[!is.finite(yg)] <- 0
    ylag <- yg*0
    mu   <- fp
    
    nl <- nrow(lagMatrix)
    zp <- as.vector(zp)
    zlag <- matrix(zp[lagMatrix[,-1]],nl,plag)
    zlag[zlag < .5] <- 0
    zlag[zlag > 0] <- 1
    xm <- matrix( yg[ lagMatrix[,-1]],nl,plag)*zlag
    ylag[lagMatrix[,1]] <- xm%*%t(betaYrMu)
    
    # random effects
    if(RANDYR){
      for(m in 1:ngroup){
        tg <- which(lagGroup == m)
        if(length(tg) == 0)next
        ylag[lagMatrix[tg,1]] <- xm[drop = FALSE,tg,]%*%t(betaYrRand[drop = FALSE,m,])
      }
    }
    ylag[!is.finite(ylag)] <- 0
    fp <- mu + ylag
  }
  
  ################ fit
  
  fit <- NULL
  
  if(SEEDDATA){
    
    fp <- exp(fp)
    
    meanDev <- sumDev/ndev
    
    la <- .getLambda(tdata[obsTrapRows,c('specPlot','year','plotyr','dcol')], 
                     sdata[,c('year','plotyr','drow')], activeArea,
                     uByGroup, as.vector(fp[obsTrapRows])*as.vector(zp[obsTrapRows]), 
                     mmu, SAMPR, USPEC, distall, obsYr, PERAREA = FALSE)
    la <- la + 1e-9
    pd  <- meanDev - 2*sum(dpois(seedCount, la, log  = TRUE))
    DIC <- pd + meanDev
    
    RMSPE <- mean(sgibbs[kg,'rmspe'])
    
    logScoreStates <- logScoreStates/nprob
    logScoreFull   <- logScoreFull/nprob
    
    fit <- list( DIC = round(DIC), scoreStates = mean(logScoreStates),
                 predictScore = round( mean(logScoreFull) ), 
                 RMSPE = signif(RMSPE,3),
                 meanPredErrSd = signif(meanPredErrSd), 
                 meanPredSd = signif(meanPredSd, 3),
                 meanInflation = signif(mean(inflation), 3))#,TGstats  = TGstats)
  }
  if(CONES)fit$rmspeCrop <- rmspeCrop
  
  inputs$treeData   <- tdata
  inputs$seedData   <- sdata
  inputs$formulaFec <- formulaFec
  inputs$formulaRep <- formulaRep
 
  inputs$ng         <- ng
  inputs$burnin     <- burnin
  inputs$keepIter   <- keepIter
  inputs$plotDims   <- plotDims
  inputs$plotArea   <- plotArea
  inputs$specNames  <- specNames
 
  inputs$xytree     <- xytree
  inputs$xytrap     <- xytrap
  inputs$yrIndex    <- yrIndex
  inputs$obsRows    <- obsRows

  inputs$maxFec     <- maxFec
  inputs$summary    <- words
  inputs$priorTable <- priorTable
 
  inputs$matYr      <- matYr 
  inputs$last0first1 <- last0first1
  if(SEEDCENSOR)inputs$censIndex <- censIndex

  if(!is.null(plotDims))inputs$plotDims <- plotDims
  if(!is.null(seedTraits))inputs$seedTraits <- seedTraits
  if(!is.null(yearEffect) & !'yearEffect' %in% names(inputs))
               inputs$yearEffect <- yearEffect

  chains <- list(bfec = .orderChain(bfgibbs, specNames), 
                 brep = .orderChain(brgibbs, specNames), 
                 sgibbs = sgibbs)
  parameters <- list( betaFec = betaFec,  
                      betaStd = betaStd,
                      betaRep = betaRep,
                      sigma = su, acfMat = acfMat,
                      pacfMat = pacfMat,
                      pacfSe = pacfSe, omegaE = omegaE,
                      omegaN = omegaN,
                      trendEst = trendEst)
  
  if(SAMPR){
    chains <- append(chains, list(rgibbs = rgibbs))
    parameters <- append(parameters, list(rMu = signif(mmu,3), 
                                          rSe = signif(mse,3)))
    inputs$priorR     <- priorR
  }
  fecPred$obs <- tdata$obs
  prediction  <-  list( fecPred = fecPred )

  if(YR){
    chains <- append(chains, list(sygibbs = sygibbs))
  }
  if(AR){
    parameters <- append(parameters, 
                         list(eigenMu = eigenMu, eigenSe = eigenSe))
    prediction <- append(prediction, list(tdataOut = tdataOut))
  }
  
  if(AR | YR){
    if(yeGr[1] %in% specNames){
      tmp <- .orderChain(bygibbsR, specNames)
    }
    
    chains     <- append(chains, list(bygibbsF = bygibbsF, bygibbsR = bygibbsR) )
    parameters <- append(parameters, list(betaYr = betaYr))
    
    if(RANDYR)parameters <- append(parameters, 
                                   list(betaYrRand = signif(betaYrRand,3),
                                        betaYrRandSE = signif(betaYrRandSE,3)))
  }
  if(SEEDDATA){
    inputs$upar       <- ug
    inputs$obsRowSeed <- obsRowSeed
    inputs$obsTrapRows <- obsTrapRows
    inputs$seedByPlot <- seedTable
    inputs$seedNames  <- seedNames
    
    parameters$upars <- upars
    parameters$dpars <- dpars
    parameters$pacsMat <- pacsMat
    chains$ugibbs <- ugibbs
    prediction <- append(prediction, list(seedPred = seedPred, sdataOut = sdataOut))
    if(USPEC)chains <- append(chains, list(priorUgibbs = priorUgibbs))
    
    if(PREDSEED) {
      prediction <- append(prediction, 
                           list(seedPredGrid = seedPredGrid,
                                treePredGrid = treePredGrid) )
      if(!'predList' %in% names(inputs))inputs <- append(inputs, list(predList = predList))
    }
  }
  if(RANDOM){
    inputs$randomEffect <- randomEffect
    chains     <- append(chains, list(agibbs = .orderChain(agibbs, specNames),
                                      aUgibbs = .orderChain(aUgibbs, specNames)) )
    parameters <- append(parameters, 
                         list(alphaMu = alphaMu, alphaSe = alphaSe,
                              aMu = aMu, aSe = aSe, 
                              alphaUMu = alphaUMu, alphaUSe = alphaUSe,
                              aUMu = aUMu, aUSe = aUSe) )
  }
  data$setupData$distall <- distall
  
  chains     <- chains[ sort( names(chains) )]
  inputs     <- inputs[ sort( names(inputs) )]
  data       <- data[ sort(names(data)) ]
  inputs$ng  <- ng
  inputs$burnin <- burnin
  inputs$obsYr  <- table(tdata$plot[tdata$obsTrap == 1], 
                         tdata$year[tdata$obsTrap == 1])
  parameters <- parameters[ sort( names(parameters) )]
  
  out <- list(inputs = inputs, chains = chains, data = data, 
              parameters = parameters, prediction = prediction )
  if(length(fit) > 0)out$fit <- fit
  
  class(out) <- 'mastif'
  out
} 
      
.chain2tab <- function(chain){
  
  if(!is.matrix(chain))chain <- matrix(chain, ncol=1)
  
  mu <- colMeans(chain)    #standardized
  
  SE <- apply(chain,2,sd)
  CI <- apply(chain,2,quantile,c(.025,.975))
  splus <- rep(' ', length=length(SE))
  splus[CI[1,] > 0 | CI[2,] < 0] <- '*'
  
  tab <- cbind( mu, SE, t(CI))
  tab <- signif(tab, 3)
  colnames(tab) <- c('estimate','SE','CI_025','CI_975')
  tab <- as.data.frame(tab)
  tab$sig95 <- splus
  attr(tab, 'note') <- '* indicates that zero is outside the 95% CI'
  
  tab
}


mergeSeedGrid <- function(gnow, gnew, scols=NULL){
  
  #  gnow - current seedPredGrid
  #  gnew - output$prediction$seedPredGrid to merge
  
  if(length(gnow) == 0)return(gnew)
  
  if(is.null(scols))scols <- c("plot","trapID","year","x","y",
                               "drow","dgrid","area","active")
  
  cnow <- colnames(gnow)[!colnames(gnow) %in% scols]
  cnew <- colnames(gnew)[!colnames(gnew) %in% scols]
  snow <- gnow[,scols]
  snew <- gnew[,scols]
  pnow <- gnow[,cnow]
  pnew <- gnew[,cnew]
  
  snow$trapID <- as.character(snow$trapID)
  snew$trapID <- as.character(snew$trapID)
  
  idnow <- columnPaste(snow$trapID, snow$year)
  idnew <- columnPaste(snew$trapID, snew$year)
  
  mm <- match(idnew, idnow)
  wf <- which(is.finite(mm))
  
  idfull <- sort(unique(c(idnow, idnew)))
  cfull  <- sort(unique(c(cnow, cnew))) 
  
  sfull <- matrix(0, length(idfull), length(cfull))
  colnames(sfull) <- cfull
  
  mm <- match(idnow, idfull)
  nn <- match(idnew, idfull)
  
  ifull <- vector( 'list', length(scols))
  names(ifull) <- scols
  
  for(k in 1:length(scols)){
    knew <- rep(NA, length(idfull))
    knew[mm] <- snow[,k]
    knew[nn] <- snew[,k]
    ifull[[k]] <- knew
  }
  ifull <- as.data.frame(ifull, stringsAsFactors = FALSE)
  
  sfull[mm, cnow] <- as.matrix( gnow[,cnow] )
  sfull[nn, cnew] <- as.matrix( gnew[,cnew] )
  ifull <- cbind(ifull, sfull)
  rownames(ifull) <- idfull
  ifull
}




meanVarianceScore <- function(output, ktree = 30, #maxSite = 100, 
                              maxArea = 20^2, cyr = 5,
                              LAGMAT = FALSE, Q = c(.5, .025, .975), nsim = 1,
                              CLOSE = FALSE){
  
  # ktree   - no. trees visited
  # maxSite - no. m2 visited
  # cyr     - no. years
  # CLOSE  = TRUE: sample small neighborhoods (selected to be close)
  # CLOSE = F: random neighborhoods
  
  lagCanopy <- lagGround <- NULL
  
  xytrap    <- output$inputs$xytrap
  xytree    <- output$inputs$xytree
  specNames <- output$inputs$specNames
  seedNames <- output$inputs$seedNames
  seedTraits  <- output$inputs$seedTraits
  
  
  ntrap     <- nrow(xytrap)
  
  fecPred      <- output$prediction$fecPred
  seedPred     <- output$prediction$seedPred
  
  seedPredGrid <- output$prediction$seedPredGrid
  if(is.null(seedPredGrid))seedPredGrid <- output$prediction$seedPred
  yr <- range(seedPredGrid$year)
  years <- yr[1]:yr[2]
  nyr   <- length(years)
  
  plots     <- sort(unique(as.character(seedPredGrid$plot)))
  nplot     <- length(plots)
  
  if(!'x' %in% colnames(seedPredGrid)){      # check this use of seedPred for mm
    xytrap <- output$inputs$xytrap
    mm <- match( as.character(seedPred$trapID), as.character(xytrap$trapID) )
    seedPredGrid$x <- xytrap$x[mm]
    seedPredGrid$y <- xytrap$y[mm]
  }
  
  treeID <- sort(unique(as.character(fecPred$treeID)))
  ntree  <- length(treeID)
  fecPred$dcol <- match(as.character(fecPred$treeID),treeID) #DIFFERENT FROM TDATA$DCOL
  
  fmat <- matrix(0, ntree, nyr)
  colnames(fmat) <- years
  rownames(fmat) <- treeID
  
  if(is.null(seedTraits)){
    seedTraits <- matrix(1,length(specNames),2)
    colnames(seedTraits) <- c('gmPerSeed','seedsPerFruit')
    rownames(seedTraits) <- specNames
  }
  
  seedTraits <- seedTraits[unique(rownames(seedTraits)),,drop = FALSE]
  
  scoreT <- scoreS <- deltaT <- deltaS <- numeric(0)
  scoreTse <- scoreSse <- deltaTse <- deltaSse <- numeric(0)
  treeCor <- trapCor <- numeric(0)
  entropy <- domain <- numeric(0)
  resourceScore <- resourceMean <- totalVar <- numeric(0)
  win <- floor(nyr/2)
  if(win > 10)win <- 10
  
  GRID <- SMASS <- FALSE
  meanCols <- grep('_meanGmM2', colnames(seedPredGrid))
  if(length(meanCols) == 0){
    meanCols <- grep('_meanM2', colnames(seedPredGrid))
    rs       <- columnSplit(colnames(seedPredGrid)[meanCols],'_meanM2')[,1]
    ws       <- which(rs %in% specNames)
    rs       <- rs[ws]
    meanCols <- meanCols[ws]
    rSeedMat <- seedPredGrid[,meanCols,drop = FALSE]*
      matrix(seedTraits[rs,'gmPerSeed'],nrow(seedPredGrid),length(rs),byrow  = TRUE)
    colnames(rSeedMat) <- paste(specNames,'_meanGmM2',sep='')
    seedPredGrid <- cbind(seedPredGrid, rSeedMat)
  }
  
  sgrid <- seedPredGrid
  GRID  <- FALSE
  
  
  meanCols <- grep('_meanGmM2', colnames(sgrid))
  
  trapID  <- as.character(sgrid$trapID)
  allTrap <- sort(unique(trapID))
  ntrap   <- length(allTrap)
  smat <- matrix(0, ntrap, nyr)
  rownames(smat) <- allTrap
  
  colnames(smat) <- years
  meanNames <- c('trees_PerTree','sites_PerSite','trees_PerYr', 'sites_PerYr')
  eNames <- c('tree-tree','site-site','tree-lag','site-lag')
  scoreNames <- c('gmTree','gmM2',eNames)
  
  kk   <- 0
  if(nsim > 1){
    pbar <- txtProgressBar(min=1,max=nplot*nsim,style=1)
    cat('\nScore\n')
  }
  
  fecAll <- seedAll <- numeric(0)   #mass basis
  totalScore <- numeric(0)
  
  if(win > 1){
    
    rjtree <- rjtrap <- rjall <- character(0)
    
    for(m in 1:nplot){
      
      wp <- fecPred$plot == plots[m]
      wo <- fecPred$obs == 1
      wm <- fecPred$matrEst > .5
      wy <- fecPred$year %in% years
      
      wk <- wp&wo&wm&wy
      wc <- which(wk)   #observed, mature
      
      if(length(wc) == 0)next
      
      yrm <- sort(unique(c(fecPred$year[wp],sgrid$year[sgrid$plot == plots[m]])))
      yrm <- yrm[yrm %in% years]
      
      #   dm   <- tdata$dcol[wc]
      dm <- fecPred$dcol[wc]
      
      dr   <- unique(dm)             # unique trees
      
      if(length(dr) <= 1)next
      
      # grid area for plot
      
      splot <- sgrid[sgrid$plot == plots[m],]
      splot <- splot[!duplicated(splot$trapID),]
      
      dx <- round( abs(diff(splot$x)))
      dy <- round( abs(diff(splot$y)))
      dx <- table( dx[dx > 1] )
      dy <- table( dy[dy > 1] )
      dx <- as.numeric(names(dx)[which.max(dx)])
      dy <- as.numeric(names(dy)[which.max(dy)])
      darea <- dx*dy
      ddist <- round( sqrt(darea), 0)
      rseq  <- c(0, seq(10, 2000, by = ddist))
      maxSite <- 1 + maxArea/darea
      
      ym   <- fecPred$year[wc]
      ym   <- match(ym,years)
      
      ytab  <- table(fecPred$year[wc])
      ykeep <- as.numeric(names(ytab))[ytab > 0]
      ckeep <- match(ykeep, years)
      
      emat <- rmat <- vmat <- matrix(0, nsim, 4)
      cmat <- matrix(0, nsim, 6)
      size <- matrix(0, nsim, 5)
      
      tseq  <- c(1:length(dr))
      
      scols <- paste('0_',c(0:win),sep='')
      
      ssmat <- matrix(0, length(rseq), length(scols))
      colnames(ssmat) <- scols
      rownames(ssmat) <- rseq
      
      tdmat <- matrix(0, length(tseq),length(scols))
      colnames(tdmat) <- scols
      rownames(tdmat) <- tseq
      
      sdmat <- nsdmat <- ssmat <- nstmat <- ssmat*0
      tdmat <- ntdmat <- ttmat <- nttmat <- tdmat*0
      
      sdmat2 <- ssmat2 <- ssmat*0
      tdmat2 <- ttmat2 <- tdmat*0
      
      rjtree <- c(rjtree,plots[m])
      
      fmat <- fmat*0
      
      yall <- sort(unique(fecPred$year[wc]))
      nyrk <- length(yall)
      
      tkeep <- ckeep
      
      fmat[ cbind(dm,ym) ] <- fecPred$fecEstMu[wc]
      rkeep <- which(rowSums(fmat) > .99)
      fm    <- fmat[drop = FALSE,rkeep,tkeep]
      sm    <- fecPred$species[match(rownames(fm),
                                     as.character(fecPred$treeID)) ]
      sm   <- as.character( sm )
      rvec <- seedTraits[sm, 'gmPerSeed']
      
      fr <- fm*rvec
      rownames(fr) <- columnPaste(rownames(fr),names(rvec),'_')
      
      
      fecAll <- append(fecAll, list(fr))    # all fecund trees, mass basis
      names(fecAll)[length(fecAll)] <- plots[m]
      
      totTT  <- totTT2 <- totSS <- totSS2 <- rep(0, 3)
      ntt <- nss <- rep(0,3)
      
      for(k in 1:nsim){
        
        kk <- kk + 1
        if(nsim > 1)setTxtProgressBar(pbar, kk)
        
        entTtree <- entTseed <- entYtree <-  entYseed <- rtot <- rTtree <- 
          rTseed <- rYtree <- rYseed <- Tk <- Ts <- Yk <- Ys <- carea <- 
          varTk <- varTs <- varYk <- varYs <- NA
        
        #   fcor <- NULL
        
        rkeep <- 1:nrow(fm)
        fk <- fm*rvec       #CHANGED TO MASS
        
        # must have cyr years
        if(length(rkeep) > ktree)rkeep <- sample(rkeep, ktree)
        ykeep <- ckeep
        
        if(length(ckeep) > cyr){
          ykeep <- sample(1:(length(ckeep) - cyr + 1), 1)
          ykeep <- ykeep:(ykeep + cyr - 1)
        }
        
        if(length(dr) > 1){   # canopy 
          
          dd    <- runif(length(rkeep)*length(ykeep), 0, 1)
          fk    <- fk[drop = FALSE,rkeep,ykeep] + dd # 0's are less than 1
          
          #    fm <- fmat
          
          if(length(fk) > 2 & length(tkeep) > 2){
            
            nyrk <- ncol(fk)
            ntrk <- nrow(fk)
            
            ff <- sample(nrow(fk))
            
            fcor <- makeCrossCov( tmat = fk[drop = FALSE,ff,], win = win, MATRIX  = TRUE, COR  = TRUE )[[1]]
            fcor[!is.finite(fcor)] <- 0
            
            if(LAGMAT){
              
              ff <- sample(nrow(fk))
              
              tt <- makeCrossCov( tmat = fk[drop = FALSE,ff,], win = win, MATRIX  = TRUE, COR = FALSE )
              mcov <- tt$lagCov
              
              if(length(mcov) > 1){
                
                if(k == 1){
                  tcovTot <- fcor
                }else{
                  tcovTot <- tcovTot + fcor
                }
                
                tmu  <- tt$lagMean
                
                ############################# total space-time covariance
                
                ml     <- as.numeric( columnSplit(colnames(mcov),'_')[,2] )   
                maxLag <- 1 + max(ml)
                
                hh <- grep('_-', colnames(mcov))
                totCov <- mcov[drop = FALSE,,-hh]
                nn <- cbind(rownames(totCov), paste(rownames(totCov),'_0',sep=''))
                totVar <- totCov[nn]*maxLag # each variance repeated for each lag (2nd-order stat.)
                totCov[nn] <- 0
                scov   <- sum(totCov*2)       # each covariance twice
                totVar <- sum(totVar) + scov
             
                totMu    <- sum(fk[ff,])     # total yield over ktree trees, lag years
                totScore <- log(totMu) - 1/2*log(totVar)
                
                tot <- cbind(totMu, totVar, totScore)
                wf  <- which(is.finite(tot))
                
                totTT[wf]  <- totTT[wf] + tot[wf]
                totTT2[wf] <- totTT2[wf] + tot[wf]^2
                
                ntt[wf] <- ntt[wf] + 1
                
                gg <- grep( paste(rownames(mcov)[1],'_-',sep=''), 
                            colnames(mcov) )                      #based on first site
                gg <- c(grep(paste(rownames(mcov)[1],'_0', sep=''), 
                             colnames(mcov) ), gg)
                tmu <- tmu[,gg,drop = FALSE]
                
                scov <- mcov[,gg,drop = FALSE]
                vars <- mcov[ cbind(rownames(mcov), paste(rownames(mcov),'_0',sep='')) ]
                vars <- matrix(vars, nrow(scov), ncol(scov))
                
                scov[-1] <- scov[-1]*2   #covariances count twice
                
                wtrow <- matrix(1:nrow(scov),nrow(scov),ncol(scov)) 
                wtcol <- matrix(1:ncol(scov),nrow(scov),ncol(scov), byrow  = TRUE) 
                
                v1 <- wtcol*vars[1]   # count diagonal elements
                v2 <- wtrow*vars
                
                scov[1] <- 0
                scov <- scov + v1 + v2
                
                scum0 <- t(apply(scov, 1, cumsum))
                scum  <- apply(t(scum0), 1, cumsum)
                if(!is.matrix(scum))scum <- scum0*0 + scum
                
                
                tcum0 <- t(apply(tmu, 1, cumsum))
                tcum <- apply(t(tcum0), 1, cumsum)
                if(!is.matrix(scum))tcum <- tcum0*0 + tcum
                
                
                tscore <- log(tcum) - 1/2*log(scum)
                rscore <- log( sum(tmu) ) - 1/2*log( length(tmu)*var(as.vector(fm)) )
                delta  <- tscore - rscore
              }
            }
            
            rmm   <- fk               #ALREADY MASS UNITS (SEE ABOVE)
            wm    <- which(rmm > 0)
            rtree <- mean( rmm[wm] ) # mean per reproductive tree
            
            Tk  <- cov(t(rmm))       # tree cov
            Yk  <- cov(rmm)          # year cov
            
            varTk <- sum(Tk)
            varYk <- sum(Yk)
            
            rTtree <- rYtree <- sum(rmm)
            
            if(!is.na(max(Tk))){
              tmp        <- var2score(rTtree, varTk, Tk) # var between trees 
              entTtree   <- tmp$entropy
            }
            if(!is.na(max(Yk))){
              tmp        <- var2score(rYtree, varYk, Yk) # var between years
              entYtree   <- tmp$entropy
            }
          } #end canopy
        }
        
        ############# seed traps or seed prediction grid
        
        smat <- smat*0   # ntrap X nyr: ykeep is index of yrs to keep
        
        wy <- unique(sgrid$year[sgrid$plot == plots[m]])
        wy <- wy[wy %in% years]
        wm <- which(sgrid$plot == plots[m] & 
                      sgrid$year == wy[1])
        ix <- sample(wm, 1)
        wo <- wm[wm != ix]
        
        dist <- NULL
        scor <- NULL
        
        if(CLOSE){
          distSite <- .distmat( sgrid$x[ix], sgrid$y[ix], 
                                sgrid$x[wo], sgrid$y[wo] )
          distSite <- distSite[distSite > 0]
          oo <- order(distSite)
          if(length(oo) > maxSite)oo <- oo[1:maxSite]
          wm <- c(ix,wo[oo])
          
          dist <- c(0,distSite[oo])
        }
        
        trapIDS <- as.character( sgrid$trapID[wm] )
        if(!is.null(dist))names(dist) <- trapIDS
        wm <- which( as.character(sgrid$trapID) %in% trapIDS &
                       sgrid$year %in% yrm &          
                       sgrid$plot == plots[m])
        j <- match(sgrid$year[wm], years)
        i <- match(as.character(sgrid$trapID[wm]),rownames(smat))
        
        smat[ cbind(i,j) ] <- rowSums(sgrid[wm, meanCols,drop = FALSE],
                                      na.rm  = TRUE)
        zmat <- smat[unique(i),]*darea            # scale by grid area
        zmat <- zmat + runif(length(zmat), 0, .1)
        rmm <- zmat[drop = FALSE, , ykeep]
        
        if(k == 1){
          seedAll <- append(seedAll, list(rmm))
          names(seedAll)[length(seedAll)] <- plots[m]
        }
        
        crr <- ncol(rmm)
        
        if(LAGMAT & crr > 2){
       
          tmp <- rmm

          tt   <- makeCrossCov( tmat = tmp, win = win, MATRIX  = TRUE, COR = FALSE )
          mcov <- tt$lagCov
          tmu  <- tt$lagMean
          
          if(length(mcov) > 1){
            
            ml     <- as.numeric( columnSplit(colnames(mcov),'_')[,2] )   
            maxLag <- 1 + max(ml)
            
            hh <- grep('_-', colnames(mcov))
            totCov <- mcov[drop = FALSE,,-hh]
            nn <- cbind(rownames(totCov), paste(rownames(totCov),'_0',sep=''))
            totVar <- totCov[nn]*maxLag # each variance repeated for each lag (2nd-order stat.)
            totCov[nn] <- 0
            scov   <- sum(totCov*2)       # covariance twice
            totVar <- sum(totVar) + scov
            #      totMu  <- sum(tmu[,1])
            totMu    <- sum(tmp)
            totScore <- log(totMu) - 1/2*log(totVar)
            
            tot <- cbind(totMu, totVar, totScore)
            wf  <- which(is.finite(tot))
            
            totSS[wf] <- totSS[wf] + tot[wf]
            totSS2[wf] <- totSS2[wf] + tot[wf]^2
            nss[wf] <- nss[wf] + 1
          }
          
          zrows <- rownames(mcov)
          zcols <- paste(rownames(mcov),'_0',sep='')
          
          if(!is.null(dist)){
            
            p1 <- paste('^',rownames(mcov)[1],'_',sep='')
            
            gg <- grep(p1,colnames(mcov))
            scov <- mcov
            #      if(length(gg) > 0){
            scov <- mcov[,gg, drop = FALSE]
            tmu <- tmu[,gg, drop = FALSE]
            #      }
            gg <- grep('-',colnames(scov))
            hh <- grep('_0',colnames(scov))
            gg <- sort(c(gg,hh))
            scov <- scov[,gg, drop = FALSE]
            tmu <- tmu[,gg, drop = FALSE]
            
          }else{
            gg <- grep( paste(rownames(mcov)[1],'_-',sep=''), 
                        colnames(mcov) )
            gg <- c(grep(paste(rownames(mcov)[1],'_0', sep=''), 
                         colnames(mcov) ), gg)
            scov <- mcov[,gg]
            tmu <- tmu[,gg]
            
            wz    <- which(!zcols %in% colnames(mcov))
            if(is.matrix(scov) & length(wz) > 0){
              zcols <- zcols[-wz]
              zrows <- zrows[-wz]
              scov  <- scov[-wz,]
              tmu   <- tmu[-wz,]
            }
          }
          
          
          if(length(mcov) > 2)vars <- mcov[ cbind(zrows, zcols) ]
          
          if(is.matrix(scov)){
            
            vars <- matrix(vars, nrow(scov), ncol(scov))
            
            scov[-1] <- scov[-1]*2   #covariances count twice
            
            wtrow <- matrix(1:nrow(scov),nrow(scov),ncol(scov)) 
            wtcol <- matrix(1:ncol(scov),nrow(scov),ncol(scov), byrow  = TRUE) 
            
            v1 <- wtcol*vars[1]   # count diagonal elements
            v2 <- wtrow*vars
            
            scov[1] <- 0
            scov    <- scov + v1 + v2
            
            scum <- t(apply(scov, 1, cumsum))
            scum <- apply(t(scum), 1, cumsum)
            
            tcum <- t(apply(tmu, 1, cumsum))
            tcum <- apply(t(tcum), 1, cumsum)
            
            tscore <- log(tcum) - 1/2*log(scum)
            rscore <- log( sum(tmu) ) - 1/2*log( length(tmu)*var(as.vector(tmp)) )
            delta  <- tscore - rscore
            
            delta  <- vec2mat(delta)
            tscore <- vec2mat(tscore)
            
            if(k == 1){
              sdmat  <- delta
              sdmat2 <- delta^2
              ssmat  <- tscore
              ssmat2 <- tscore^2
            }else{
              sdmat  <- sdmat + delta
              sdmat2 <- sdmat2 + delta^2
              ssmat  <- ssmat + tscore
              ssmat2 <- ssmat2 + tscore^2
            }
            
          }
          
          scor <- makeCrossCov( tmat = tmp, win = win, MATRIX  = TRUE, COR  = TRUE )[[1]]
          scor[!is.finite(scor)] <- 0
          if(k == 1){
            scovTot <- scor
          }else{
            scovTot <- scovTot + scor
          }
          
          if(k == 1){
            
            rjtrap   <- c(rjtrap, plots[m])
          }
        }
 
        if(is.matrix(rmm) & length(rmm) > 1){
          
          rcol <- rmm
          if(ncol(rcol) < nrow(rcol)){
            ss <- sample(nrow(rcol), ncol(rcol)-1)
            rcol <- rcol[ss,]
            rcol <- rcol + .tnorm(length(rcol), 0, .00001, .00001/2, .00001)
          }
          
          rtot <- sum(rmm)
          
          Ts <- cov(t(rcol))  # site cov
          Ys <- cov(rmm)     # year cov
          
          varTs <- sum(Ts)
          varYs <- sum(Ys)
          
          rTseed <- rYseed <- sum(rmm)
          
          tmp        <- var2score(rTseed, varTs, Ts)
          entTseed   <- tmp$entropy
          tmp        <- var2score(rYseed, varYs, Ys)
          entYseed   <- tmp$entropy
          
          n1 <- n2 <- n3 <- n4 <- NA
          if(is.matrix(Tk))n1 <- nrow(Tk)
          if(is.matrix(Ts))n2 <- nrow(Ts)
          if(is.matrix(Yk))n3 <- nrow(Yk)
          if(is.matrix(Ys))n4 <- nrow(Ys)
          
          emat[k,] <- c(entTtree, entTseed, entYtree, entYseed)
          cmat[k,] <- c(rtree, rtot)
          rmat[k,] <- c(rTtree, rTseed, rYtree, rYseed)
          size[k,] <- c(n1,n2,n3, n4, round(carea)) 
          vmat[k,] <- c(varTk, varTs, varYk, varYs)
        }
        
        
      }#######################
      
      lagCanopy <- append(lagCanopy, list(tcovTot/nsim))
      lagGround <- append(lagGround, list(scovTot/nsim))
      
      tcm <- totTT/ntt
      tcs <- totTT2/ntt - tcm^2 
      tcs[tcs < 0] <- 0
      tcm[2] <- sqrt(tcm[2])
      tcs[2] <- sqrt(tcs[2])
      tcs <- sqrt( tcs )
      
      tgm <- totSS/nss
      tgs <- totSS2/nss - tgm^2 
      tgs[tgs < 0] <- 0
      tgm[2] <- sqrt(tgm[2])
      tgs[2] <- sqrt(tgs[2])
      tgs <- sqrt( tgs )
      
      t1 <- rbind(tcm,tcs)
      t2 <- rbind(tgm,tgs)
      
      totScore <- cbind(t1, t2)
      
      ctt <- c('mu','stdDev','score')
      
      colnames(totScore) <- paste( c(rep('canopy',3),rep('ground',3)),
                                   ctt,sep='_' )
      rownames(totScore) <- paste(plots[m],c('mu','se'),sep='_')
      
      totalScore <- rbind(totalScore,totScore)
      
      
      if( LAGMAT){
        TREE <- GROUND <- TRUE
        deltaTree <- scoreTree <- deltaSeed <- scoreSeed <- 
          deltaTrSe <- scoreTrSe <- deltaSdSe <- scoreSdSe <- NULL
        
        if(sum(ntdmat) == 0)TREE <- FALSE
        if(sum(nttmat) == 0)GROUND <- FALSE
        
        if(TREE){
          deltaTree <- tdmat/nsim
          scoreTree <- ttmat/nsim
          deltaTrSe <- tdmat2/nsim - deltaTree^2
          scoreTrSe <- ttmat2/nsim - scoreTree^2
          wr <- which(rowSums(deltaTree,na.rm  = TRUE) != 0)
          wc <- which(colSums(deltaTree,na.rm  = TRUE) != 0)
          deltaTree <- deltaTree[drop = FALSE,wr,wc]
          deltaTrSe <- sqrt(deltaTrSe[drop = FALSE,wr,wc])
          scoreTree <- scoreTree[drop = FALSE,wr,wc]
          scoreTrSe <- sqrt(scoreTrSe[drop = FALSE,wr,wc])
        }
        if(GROUND){
          deltaSeed <- sdmat/nsim
          scoreSeed <- ssmat/nsim
          deltaSdSe <- sdmat2/nsim - deltaSeed^2
          scoreSdSe <- ssmat2/nsim - scoreSeed^2
          wr <- which(rowSums(deltaSeed,na.rm  = TRUE) != 0)
          wc <- which(colSums(deltaSeed,na.rm  = TRUE) != 0)
          deltaSeed <- deltaSeed[drop = FALSE,wr,wc]
          scoreSeed <- scoreSeed[drop = FALSE,wr,wc]
          deltaSdSe <- sqrt(deltaSdSe[drop = FALSE,wr,wc])
          scoreSdSe <- sqrt(scoreSdSe[drop = FALSE,wr,wc])
        }
        
        deltaT <- append(deltaT, list(deltaTree))
        deltaS <- append(deltaS, list(deltaSeed))
        scoreT <- append(scoreT, list(scoreTree))
        scoreS <- append(scoreS, list(scoreSeed))
        deltaTse <- append(deltaTse, list(deltaTrSe))
        deltaSse <- append(deltaSse, list(deltaSdSe))
        scoreTse <- append(scoreTse, list(scoreTrSe))
        scoreSse <- append(scoreSse, list(scoreSdSe))
        
      }
      
      emat[!is.finite(emat)] <- NA
      
      ee <- signif(t( apply(emat, 2, quantile, Q,na.rm  = TRUE )),3)
      rownames(ee) <- paste(plots[m],eNames,sep='_')
      vv <- t( apply(vmat, 2, quantile, Q,na.rm  = TRUE ))
      rownames(vv) <- paste(plots[m],eNames,sep='_')
      
      rjall <- c(rjall,plots[m])
      
      ii     <- apply(size,2,mean)
      domain <- rbind(domain, ii)
      
      entropy  <- rbind( entropy, ee )
      totalVar <- rbind( totalVar, vv)
      
    } #########end plot loop
    
    if(LAGMAT){
      if(length(lagCanopy) > 0)names(lagCanopy) <- rjtree
      if(length(lagGround) > 0)names(lagGround) <- rjtrap
      names(scoreT) <- rjtree
      names(scoreS) <- rjtrap
      names(deltaT) <- rjtree
      names(deltaS) <-  rjtrap
    }
    
    if(length(domain) > 0){
      rownames(domain) <- rjall
      colnames(domain) <- c('trees','sites','treeYr','siteYr','areaM2')
    }
  }
  
  if(LAGMAT){
    tmp$lagCanopy <- lagCanopy
    tmp$lagGround <- lagGround
    tmp$scoreSeed <- scoreS
    tmp$scoreTree <- scoreT
    tmp$deltaSeed <- deltaS
    tmp$deltaTree <- deltaT
    tmp$scoreSeedSe <- scoreSse
    tmp$scoreTreeSe <- scoreTse
    tmp$deltaSeedSe <- deltaSse
    tmp$deltaTreeSe <- deltaTse
  }
  for(k in 1:length(tmp)){
    if(is.null(tmp[[k]]) | is.numeric(0))next
    kcol <- which(sapply(tmp[[k]],is.numeric))
    jcol <- which(sapply(tmp[[k]],is.factor))
    kcol <- intersect(kcol,!jcol)
    for(j in kcol){
      tmp[[k]][,j][ is.finite( tmp[[k]][,j] ) ] <- NA
    }
  }
  tmp$fecAll  <- fecAll
  tmp$seedAll <- seedAll
  tmp$totalScore <- totalScore
  tmp$gridArea <- sqrt( darea )
  
  tmp
  
}

var2score <- function(rmean, totVr, rvar){
  
  # rmean - mean over sites or years
  # rvar  - covariance matrix
  # totVr - total variance
  # ndim  - dimension of covariance matrix
  
  if(length(rvar) < 4)return( list(score = NA, entropy = NA) )
  if(nrow(rvar) > 100){
    ss <- sample(nrow(rvar),100)
    rvar <- rvar[ss,ss]
  }
  
  ndim <- nrow(rvar)
  
  score <- log(rmean) - 1/2*suppressWarnings(log(totVr))
  
  dt    <- determinant(rvar)$modulus
  if(!is.finite(dt)){
    ev <- eigen(rvar)$values
    dt <- sum( log(ev[ev > 0]) )
  }
  entr  <- ndim/2*(1 + log(2*pi)) + dt/2
  entr  <- entr/ndim
  list(score = score, entropy = entr)
}

crossCovSetup <- function(tmat, win){
  
  nt <- ncol(tmat)
  if(win > nt/2)win <- floor(nt/2)
  
  ni    <- nrow(tmat)
  lead  <- -c(-win:win)
  ntt   <- length(lead)
  mgrid <- as.matrix( expand.grid(1:ntt,1:ni,1:ni ) )
  colnames(mgrid) <- c('t','i1','i2')
  
  ld    <- lead[mgrid[,'t']]
  mgrid <- cbind(ld,mgrid)
  colnames(mgrid)[1] <- 'lead'
  
  cc    <- columnPaste(mgrid[,'i1'],mgrid[,'i2'])
  mdex  <- columnPaste(mgrid[,'t'],cc)
  keep  <- which(mgrid[,'i2'] >= mgrid[,'i1'])
  mgrid <- mgrid[keep,]
  mdex  <- mdex[keep]
  
  mgrid <- as.data.frame(mgrid)
  
  rn <- rownames(tmat)
  if(!is.null(rn)){
    mgrid$ID1 <- rn[mgrid[,'i1']]
    mgrid$ID2 <- rn[mgrid[,'i2']]
  }
    
  list(win = win, ntt = ntt, ni = ni, mgrid = mgrid, 
       mdex = mdex, lead = lead)
}
  
makeCrossCov <- function(tmat, win = 5, MATRIX = FALSE, COR = FALSE ){
  
  # tmat - responses by time matrix
  # cross covariance of each row against population
  # MATRIX  - n by n*lag matrix
  # !MATRIX - n*n*lag vector
  # COR - correlation matrix
  
  tiny <- 1e-8
  
  nt  <- ncol(tmat)
  mid <- round(nt/2)
  rt  <- mid + c(-win,win)
  
  if(rt[2] <= rt[1])return(NULL)
  
  imat <- sweep(tmat, 1, rowMeans(tmat), '-')
  imat[imat == 0] <- tiny                        # no variation
  
  if(win > (ncol(tmat)-1))win <- ncol(tmat) - 1
  
  tmp <- crossCovSetup(tmat, win)
  ntt <- tmp$ntt
  ni  <- tmp$ni
  mgrid <- tmp$mgrid
  mdex  <- tmp$mdex
  lead  <- tmp$lead
  win   <- tmp$win
  
  crossCov <- rep(0,length(mdex))
  names(crossCov) <- mdex
  totSeed <- crossCov
  
  if(COR)ivar <- apply(imat,1,var)*(nt-1)/nt
  
  for(i in 1:(win+1)){
    
    ii <- 1:(nt - win + i - 1)
    pp <- (win - i + 2):nt
    
    ii  <- ii[ii > 0]
    pp  <- pp[pp <= nt]
    ldd <- pp[1] - ii[1]
    
    for(m in 1:ni){
      wm   <- which(mgrid[,'i1'] == m & mgrid[,'lead'] == ldd )
      mdx  <- mgrid[drop = FALSE,wm,]
      
      tres <- rowMeans( imat[mdx[,'i1'],ii,drop = FALSE]*imat[mdx[,'i2'],pp,drop = FALSE] )
      
      tm2 <- tmat[mdx[,'i2'],pp,drop = FALSE]
      if(ldd == 0)tm2[ rownames(tm2) == rownames(tmat[mdx[,'i1'],]) ] <- 0
      tsum <- rowMeans( tmat[mdx[,'i1'],ii,drop = FALSE] + tm2 )
      
      if(COR)tres <- tres/sqrt(ivar[mdx[1,'i1']]*ivar[mdx[,'i2']])
      crossCov[ wm ] <- tres
      totSeed[ wm ]  <- tsum
        
      if(ldd != 0){
        wm   <- which(mgrid[,'i1'] == m & mgrid[,'lead'] == -ldd )
        mdx  <- mgrid[drop = FALSE,wm,]
        tres <- rowMeans( imat[mdx[,'i1'], pp,drop = FALSE]*imat[mdx[,'i2'], ii,drop = FALSE] )
        if( identical(mdx[,'i1'], mdx[,'i2']) ){
          tres <- rowMeans( tmat[mdx[,'i1'],ii,drop = FALSE] )
        }else{
          tsum <- rowMeans( tmat[mdx[,'i1'], pp,drop = FALSE] + tmat[mdx[,'i2'], ii,drop = FALSE] )
        }
        
        if(COR)tres <- tres/sqrt(ivar[mdx[1,'i1']]*ivar[mdx[,'i2']])
        crossCov[ wm ] <- tres
        totSeed[ wm ] <- tsum
      }
    }
  }
  
  totSeed <- round(totSeed, 2)
  
  covMu <- cbind(mgrid, crossCov, totSeed)
  
  if(!MATRIX)return( lagCov = covMu, lagMean = NULL )
  
  t2 <- covMu$t[covMu$lead == 0][1]
  
  rn <- columnPaste(covMu[,'ID1'],covMu[,'lead'],'_')
  cn <- columnPaste(covMu[,'ID2'],t2,'_')
  rt <- rn[!duplicated(rn)]
  ct <- cn[!duplicated(cn)]
  
  stmat <- ttmat <- matrix(0,length(rt),length(ct))
  rownames(stmat) <- rownames(ttmat) <- rt
  colnames(stmat) <- colnames(ttmat) <- ct
  cii <- cbind(rn,cn)
  stmat[ cii ] <- covMu[,'crossCov']
  stmat <- t(stmat)
  
  ttmat <- stmat*0
  tmean <- rowMeans(tmat)
  rm    <- columnSplit(ct,'_')[,1]
  ttmat[1:nrow(ttmat),] <- tmean[rm]
  
  rnn <- columnSplit(rownames(stmat),'_')[,1]
  rownames(stmat) <- rnn
  list(lagCov = stmat, lagMean = ttmat)
}
  
.updateBetaAR <- function(betaYr, yg, mu, z,lagGroup, plag, ngroup, sg)  {
  
  # AR(p), fixed groups
  
  ylag <- yg*0
  
  for(m in 1:ngroup){
    
    lmm <- lagGroup[[m]]
    lmm <- lmm[drop = FALSE,z[lmm[,plag+1]] == 1,]
    if(nrow(lmm) <= (ncol(lmm)+5))next
    
    ym  <- yg[lmm[,1]] - mu[lmm[,1]]
    xm  <- matrix( yg[ lmm[,-1] ], nrow(lmm) )
    V   <- solve( crossprod(xm) )*sg
    v   <- crossprod(xm,ym)/sg
    tmp <- rmvnormRcpp(1,V%*%v,V) 
    whi <- which( abs(tmp) > 1)
    if(length(whi) > 0){
      tmp <- .tnormMVNmatrix(tmp,tmp,V,tmp*0 - 1,tmp*0 + 1, 
                             whichSample=c(1:length(tmp)) )
    }
    betaYr[m,] <- tmp
    
    lmm <- lagGroup[[m]]
    xm  <- matrix( yg[ lmm[,-1] ], nrow(lmm) )
    ylag[lmm[,1]] <- xm%*%t(tmp)
  }
  
  wfinite <- which(!betaYr == 0)
  list( betaYr = betaYr, ylag = ylag, wfinite = wfinite)
}

.updateBetaAR_RE <- function(betaYrF, betaYrR, Alag,
                             yg, mu, z, lagGroup, lagMatrix, plag, ngroup, sg){
  
  # AR(p), random groups
  # betaYrF - 1 by plag fixed effects
  # betaYrR - ngroup by plag random effects
  # Alag    - random effects covariance
  
  ylag <- yg*0
  
  if(ngroup == 1){
    
    # fixed effects
    cg   <- which(z[lagMatrix[,plag+1]] == 1)  # mature plag yr ago
    
    if(length(cg) <= plag){
      warning('too few mature individuals in AR groups--fewer groups or smaller p')
      return( list( betaYrF = betaYrF, betaYrR = betaYrR, ylag = ylag, Alag = Alag) )
    }
    lmm  <- lagMatrix[drop = FALSE,cg,]
    xm   <- matrix( yg[ lmm[,-1] ], nrow(lmm))
    mvec <- yg[ lmm[,1] ] - mu[ lmm[,1] ] - 
      rowSums( betaYrR[ lagGroup[cg],]*xm )  # remove random AR effects
    
    v <- crossprod(xm,mvec)/sg
    V <- solve(crossprod(xm)/sg + diag(1,plag))
    
    tmp <- rmvnormRcpp(1,V%*%v,V) 
    whi <- which( abs(tmp) > 1)                               # for stability (may not be desirable)
    if(length(whi) > 0){
      tmp <- .tnormMVNmatrix(tmp,tmp,V,tmp*0 - 1,tmp*0 + 1, 
                             whichSample=c(1:length(tmp)) )
    }
    
    betaYrF <- tmp
    
    ylag[lmm[,1]] <- ylag[lmm[,1]] + xm%*%t(betaYrF)
    return( list( betaYrF = betaYrF, betaYrR = betaYrR, ylag = ylag, Alag = Alag) )
  }
  
  # random effects
  
  AIlag <- solve(Alag)
  
  for(m in 1:ngroup){
    
    tg  <- lagGroup == m
    cg  <- z[lagMatrix[,plag+1]] == 1  # mature plag yr ago
    wm  <- which(tg & cg)
    
    if(length(wm) < 2){      #
      betaYrR[m,] <- 0
      next
    }
      
    if(length(wm) < plag){
      V <- Alag
      v <- t(betaYrF)*0
    }else{
      lmm <- lagMatrix[drop = FALSE, wm, ]
      xm  <- matrix( yg[ lmm[,-1] ], nrow(lmm) )
      mvec <- yg[lmm[,1]] - mu[lmm[,1]] # - xm%*%t(betaYrF)  #
      v    <- crossprod(xm, mvec)/sg
      V    <- solve(crossprod(xm)/sg + AIlag)
    }
    tmp <- rmvnormRcpp(1,V%*%v,V) 
    whi <- which( abs(tmp) > 1)
    if(length(whi) > 0){
      tmp <- .tnormMVNmatrix(tmp,tmp,V,tmp*0 - 1,tmp*0 + 1, 
                             whichSample=c(1:length(tmp)) )
    }
    betaYrR[m,] <- tmp 
    if(length(wm) > plag)ylag[lmm[,1]] <- ylag[lmm[,1]] + xm%*%betaYrR[m,]
  }
  
  # random effects covariance
  wr   <- which(rowSums(betaYrR) != 0)
  LL   <- crossprod(betaYrR[drop = FALSE,wr,]) 
  Alag <- .updateCovariance(LL, diag(1,plag), length(wr), plag+1)
  
  list( betaYrF = betaYrF*0, betaYrR = betaYrR, ylag = ylag, Alag = Alag)
}

acfEmp <- function(res, irow, time){
  
  # empirical (for model-based use .updateBetaAc)
  # assumes res values have unique irow and time
  # detrend
  mt <- max(time)
  st <- c(0:(mt-1))
  ni <- max(irow)
  
  resMat <- matrix(0,ni,mt)
  resMat[ cbind(irow,time) ] <- res
  xy <- t( matrix(st*t(resMat),mt) ) 
  xx <-  matrix(st^2,ni,mt,byrow  = TRUE) 
  bs <- rowSums(xy,na.rm  = TRUE)/rowSums(xx,na.rm  = TRUE)
  bi <- rowMeans(resMat,na.rm  = TRUE) - bs*mean(st)
  yy <- resMat - bi + matrix(bs,ncol=1)%*%st 
  res <- yy[ cbind(irow, time) ]
  
  FF  <- .myBy(res, irow, time, matrix(0,max(irow),max(time)),fun='sum')
  FF  <- crossprod(FF)
  z   <- FF[lower.tri(FF, diag  = TRUE)]
  ii  <- row(FF)-col(FF)
  ii  <- ii[lower.tri(ii, diag  = TRUE)] + 1
  tmp <- as.vector(.myBy(z,ii,ii*0+1,fun='sum'))
  tmp  <- tmp/tmp[1]

  names(tmp) <- paste('lag',st,sep='-')
  
  tmp
}

pacf <- function(xacf){
  
  nj <- length(xacf)
  st <- c(0:(nj-1))
  rmat <- diag(1,nj)
  rmat <- matrix( xacf[1 + abs( row(rmat)-col(rmat) )], nj, nj)
  tmp <- c(1, solve( toeplitz(xacf[1:(nj-1)]), xacf[2:nj] ) )
  names(tmp) <- paste('lag',st,sep='-')
  tmp
}

.mapSetup <- function(xlim,ylim,scale=NULL,widex=10.5,widey=6.5){  
  
  #scale is x per inch
  #new means not a new plot
  
  if(is.null(scale))scale <- 1
  
  px   <- diff(xlim)/scale
  py   <- diff(ylim)/scale
  
  if(px > widex){
    dx <- widex/px
    px <- widex
    py <- py*dx
  }
  if(py > widey){
    dx <- widey/py
    py <- widey
    px <- px*dx
  }
  
  par(pin=c(px,py))
  invisible( c(px,py) )
}

checkPlotDims <- function(plots, years, xytree = NULL, xytrap = NULL, 
                          plotDims, plotArea, verbose = FALSE){
  
  if(is.null(plotDims)){
    plotDims <- getPlotDims(xytree, xytrap)
  }else{
    rownames(plotDims) <- .fixNames(rownames(plotDims), all  = TRUE)$fixed
    if(ncol(plotDims) != 5)
      stop('\nplotDim must have 5 columns: xmin, xmax, ymin, ymax, area\n')
    
    wc <- which(!plots %in% rownames(plotDims))
    if(length(wc) > 0){
      xx <- paste('\nNote',plots[wc],' missing from plotDims\n ')
      if(verbose)cat(xx)
      moreRows <- matrix(NA, length(wc), 5)
      rownames(moreRows) <- plots[wc]
      plotDims <- rbind(plotDims,moreRows)
    }
  }
  ww <- which(!plots %in% rownames(plotDims))
  if(length(ww) > 0){
    pp <- matrix(NA, length(ww), ncol(plotDims))
    rownames(pp) <- plots[ww]
    plotDims <- rbind(plotDims, pp)
  }
  plotDims <- plotDims[drop = FALSE,plots,]
  
  if(is.null(plotArea)){
    plotArea <- matrix(plotDims[,'area'],nrow(plotDims),length(years))
    rownames(plotArea) <- plots
    colnames(plotArea) <- years
  }else{
    rownames(plotArea) <- .fixNames(rownames(plotArea), all  = TRUE)$fixed
  }
  wc <- which(!plots %in% rownames(plotArea))
  if(length(wc) > 0){
    moreRows <- matrix(NA, length(wc), ncol(plotArea))
    rownames(moreRows) <- plots[wc]
    colnames(moreRows) <- colnames(plotArea)
    plotArea <- rbind(plotArea,moreRows)
    plotArea <- plotArea[drop = FALSE,plots,]
  }
    
  list(plotDims = plotDims, plotArea = plotArea)
}

mastMap <- function(mapList){
  
  # if PREDICT, needs seedPredGrid, treePredGrid
  
  mapPlot <- mapYears <- treeSymbol <- xlim <- ylim <- NULL
  seedMax <- fecMax <- fecPred <- seedPred <- RMD <- 
    seedPredGrid <- treePredGrid <- treeData <- seedData <- 
    specNames <- seedNames <- NULL
  
  PREDICT <- LEGEND <- SCALEBAR  <- verbose <- FALSE
  MAPTRAPS <- MAPTREES <- COLORSCALE <- TRUE  
  
  treeScale <- trapScale <- plotScale  <- 1
  scaleValue <- 20
  cex <- .9
  mfrow <- c(1,1)
  
  if( 'chains' %in% names(mapList) )class(mapList) <- 'mastif'
  
  if( class(mapList) == 'mastif' ){
    indat <- c('treeData', 'seedData', 'specNames', 'seedNames',
               'xytree', 'xytrap')
    mi    <- match(indat, names(mapList$inputs))
    for(i in mi){
      mapList <- append( mapList, mapList$inputs[i] )
    }
    fecPred  <- mapList$prediction$fecPred
    seedPred <- mapList$prediction$seedPred
    treePredGrid  <- mapList$prediction$treePredGrid
    seedPredGrid  <- mapList$prediction$seedPredGrid
  }
  
  minPars <- c('specNames','seedNames','treeData','seedData','xytree','xytrap',
               'mapPlot','mapYears')
  
  wk <- which(!minPars %in% names(mapList))
  if(length(wk) > 0){
    mp <- paste0(minPars[wk], collapse=', ')
    stop('\nmissing from mapList: ',mp)
  }
  minPars <- c(minPars, 'treeSymbol')
  for(k in 1:length(minPars)){
    wk <- which(names(mapList) == minPars[k])
    assign( minPars[k], mapList[[wk[1]]] )
  }
  
  treeData <- treeData[as.character(treeData$species) %in% specNames,]
  
  
  tdat <- treeData
  sdat <- seedData

  rm(treeData)
  rm(seedData)
  
  plotVars <- c('fecPred','seedPred','PREDICT','treeScale',
                'trapScale','xlim','ylim','MAPTRAPS','MAPTREES','seedMax','NULL',
                'fecMax','mfrow','LEGEND','plotScale ','SCALEBAR','scaleValue',
                'COLORSCALE','RMD')
  
  for(k in 1:length(plotVars)){
    wk <- which(names(mapList) == plotVars[k])
    if(length(wk) == 0)next
    assign( plotVars[k], mapList[[wk]] )
  }
  
  
  mapPlot     <- .fixNames(mapPlot, all  = TRUE)$fixed
  tdat$plot  <- .fixNames(tdat$plot, all  = TRUE)$fixed
  sdat$plot  <- .fixNames(sdat$plot, all  = TRUE)$fixed
  xytree$plot <- .fixNames(xytree$plot, all  = TRUE)$fixed
  xytrap$plot <- .fixNames(xytrap$plot, all  = TRUE)$fixed
  
  if(length(specNames) == 1)LEGEND <- FALSE
  
  xytrap$trapID <- columnPaste(xytrap$plot,xytrap$trap)
  sdat$trapID <- columnPaste(sdat$plot,sdat$trap)
  xytree$treeID <- columnPaste(xytree$plot,xytree$tree)
  tdat$treeID <- columnPaste(tdat$plot,tdat$tree)
  
  
  mapTreeYr <- mapYears
  
  wp <- which(tdat$plot == mapPlot)
  if(length(wp) == 0 & MAPTREES)stop('\nno trees on this plot\n')
  
  censYr <- table(tdat$plot,tdat$year)[drop = FALSE,mapPlot,]
  cmiss  <- which(censYr[as.character(mapYears)] == 0)
  cmiss  <- as.numeric(names(cmiss))
  cmiss  <- sort(c(cmiss, mapYears[!mapYears %in% names(censYr)]))
  
  if( MAPTREES & length(cmiss) > 0 ){
    
    yy <- as.numeric(colnames(censYr))
    yy[censYr == 0] <- Inf
    
    dy <- abs( outer( mapYears, yy, '-' ) )
    cy <- yy[ apply(dy, 1, which.min) ]
    mapTreeYr <- cy
    
    mpp <- paste0(mapTreeYr, collapse=', ')
  }
  
  wtdata <- which(tdat$year %in% mapTreeYr & 
                  as.character(tdat$plot) %in% mapPlot)
  tree <- tdat[wtdata,]
  
  wt <- match(as.character(tree$treeID),
              as.character(xytree$treeID) )
  tree$x <- xytree$x[wt]
  tree$y <- xytree$y[wt]
  
  fmat <- tree
  if(is.null(treeSymbol)){
    treeSymbol <- tree$diam
  }else{
    treeSymbol <- treeSymbol[wtdata]
  }
  
  wsdata <- which(sdat$year %in% mapYears & 
                    as.character(sdat$plot) %in% mapPlot)
  seed <- sdat[wsdata,]
  ws   <- match(as.character(seed$trapID),
              as.character(xytrap$trapID) )
  seed$x <- xytrap$x[ws]
  seed$y <- xytrap$y[ws]
  
  if(is.null(seedPredGrid))PREDICT <- FALSE
  
  wtpred <- numeric(0)
  
  if(PREDICT){
    wtpred <- which(treePredGrid$year %in% mapTreeYr & 
                      treePredGrid$plot %in% mapPlot)
    if(length(wtpred) == 0)PREDICT <- FALSE

    treePredGrid <- treePredGrid[wtpred,]

    wt <- match(as.character(treePredGrid$treeID),
                as.character(xytree$treeID) )
    treePredGrid$x <- xytree$x[wt]
    treePredGrid$y <- xytree$y[wt]
    
    fmat <- treePredGrid
    
    if(is.null(treeSymbol)){
      treeSymbol <- treePredGrid$fecEstMu
      fmat <- tree
    }
    
    wsdata <- which(seedPred$year %in% mapYears & 
                      as.character(seedPred$plot) %in% mapPlot)
    seedPred <- seedPred[wsdata,]
    wt <- match(as.character(seedPred$trapID),
                        as.character(xytrap$trapID) )
    seedPred$x <- xytrap$x[wt]
    seedPred$y <- xytrap$y[wt]
    
    wspred <- which(seedPredGrid$year %in% mapYears & 
                      seedPredGrid$plot %in% mapPlot)
    seedPredGrid <- seedPredGrid[wspred,]
    
    scols <- paste(specNames, '_meanM2', sep='')
    scols <- c('plot','trapID','year','x','y',scols)
    
    seedPredGrid <- rbind(seedPred[,scols],seedPredGrid[,scols])
    
    
  }else{
    if(is.null(treeSymbol))treeSymbol <- tree$diam
  }
  
  if(length(wtdata) == 0 & length(wtpred) == 0 & MAPTREES){
    if(verbose)cat('\nNo obs or preds for mapPlot, mapYears\n')
    return(add = FALSE)
  }
  
  SEED <- TRUE
  snames <- paste(seedNames,'_meanM2',sep='')
  seedCount <- as.matrix(seed[,seedNames, drop = FALSE])
  if(nrow(seed) == 0)SEED <- FALSE
  
  if(is.null(seedMax) & length(seedCount) > 0)
    seedMax <- max( seedCount, na.rm  = TRUE )
  
  if(is.null(fecMax)){
    if(MAPTREES){
      fecMax <- max( treeSymbol, na.rm  = TRUE )
    }else{
      fecMax <- 0
    }
  }
  fecMax <- max(fecMax, na.rm  = TRUE)
  
  nspec <- length(specNames)
  
  cfun <- colorRampPalette( c('#66c2a5','#fc8d62','#8da0cb') )
  specCol <- cfun(nspec) 
  names(specCol) <- specNames
  
  xlimk <- ylimk <- numeric(0)
  dx <- dy <- numeric(0)
  
  npp  <- length(mapPlot)
  
  for(j in 1:npp){
    
    if( is.null(xlim) ){
      wxy1 <- which(xytree$plot == mapPlot[j])
      wxy2 <- which(xytrap$plot == mapPlot[j])
      
      xlimj <- range( c(xytree[wxy1,'x'], xytrap[wxy2,'x']) )
      ylimj <- range( c(xytree[wxy1,'y'], xytrap[wxy2,'y']) )
      if(PREDICT){
        wj <- which(seedPredGrid$plot == mapPlot[j])
        xlimj <- range( c(xlimj, seedPredGrid$x[wj]) )
        ylimj <- range( c(ylimj, seedPredGrid$y[wj]) )
      }

    }else{
      xlimj <-  xlim
      ylimj <-  ylim
    }
    dxj <- diff(xlimj)
    dyj <- diff(ylimj)
    
    dx <- c(dx, dxj)
    dy <- c(dy, dyj)
    
    xlimk <- rbind(xlimk, xlimj)
    ylimk <- rbind(ylimk, ylimj)
  }
  xlimit <- matrix(xlimk,npp,2,byrow = FALSE)
  ylimit <- matrix(ylimk,npp,2,byrow = FALSE)
  rownames(xlimit) <- rownames(ylimit) <- mapPlot
  
  rr  <- apply( rbind(xlimit, ylimit), 2, range)
  sc  <- max( apply(rr,1,diff) )/20
  
  opin <- par()$pin
  
  obs <- oyr <- numeric(0)
  if(SEED){
    stab <- with( seed, table(plot, year) )
    stab <- stab[drop = FALSE, mapPlot, ]
    obs  <- stab[, colnames(stab) %in% mapYears, drop = FALSE]
    oyr  <- as.numeric( colnames(obs)[ colSums(obs) > 0 ] )
  }
  
  pyr  <- numeric(0)
  pred <- NULL
  
  if(PREDICT){
    ptab <- with( seedPredGrid, table(plot, year) )
    if(!mapPlot %in% rownames(ptab)){
      if(verbose)cat('\nNo prediction for this plot\n')
      PREDICT <- FALSE
    }else{
      ptab <- ptab[drop = FALSE, mapPlot, ]
      wss  <- which(colnames(ptab) %in% mapYears)
      if(length(wss) == 0){
        if(verbose)cat('\nNo prediction for this plot-year\n')
        PREDICT <- FALSE
      }else{
        pred <- ptab[, colnames(ptab) %in% mapYears, drop = FALSE]
        pyr  <- colnames(pred)[ colSums(pred) > 0 ]
      }
    }
    pyr <- as.numeric(pyr)
  }
  yr  <- sort( unique( c(oyr, pyr) ) )
  if( length(oyr) > 0 ){
    oyr <- oyr[oyr %in% yr]
    obs <- obs[, as.character(oyr), drop = FALSE]
  }else{
    obs <- NULL
  }

  if(PREDICT){
    if( length(pyr) > 0 ){
      pyr <- pyr[pyr %in% yr]
      pred <- pred[, as.character(pyr), drop = FALSE]
    }
  }
  
  specAll  <- table(tree$species)
  specAll  <- specAll[specAll > 0]
  specPred <- table(tree$species[tree$plot %in% rownames(pred) &
                                   tree$year %in% pyr])
  specPred  <- names(specPred)[specPred > 0]
  colList <- numeric(0)
  
  for(j in 1:npp){
    
    WOJ <- WPJ <- FALSE
    
    jobs <- jpred <- jyr <- numeric(0)
    
    if(!is.null(obs)){
      jobs  <- obs[drop = FALSE, mapPlot[j],]
      WOJ <- TRUE
    }
    if(!is.null(pred)){
      jpred <- pred[drop = FALSE, mapPlot[j],]
      WPJ <- TRUE
    }
    jyr   <- sort(unique(as.numeric(c(colnames(jobs),colnames(jpred)))))
    njyr  <- length(jyr)
    
    if(is.null(mfrow))mfrow <- c(1,1)
    
    suppressWarnings(
      par(bty='o',mar=c(1,.4,2,.4), oma=c(3,3,1,1))
    )
    if(LEGEND)par(oma=c(3,3,1,3))
    
    par(mfrow=mfrow)
    
    if(!is.null(RMD)){
      mfrow <- c(1,1)
      njyr  <- 1
      par( mfrow=mfrow, mar=c(2,2,2,1), bty='o' )
    }
    
    if(njyr == 1){
      scale <- max( c(dx[j],dy[j]) )/2
      if(scale < 50)scale <- 50
    }else{
      if(prod(mfrow) == 1)mfrow <- .getPlotLayout(njyr)$mfrow
      par( mfrow=mfrow, mar=c(2,2,2,1), bty='o' )
      scale <- max(c(dx[j],dy[j]))/plotScale/3*max(mfrow)
    }
    
    if(!is.null(RMD))scale <- scale*1.1
    .mapSetup(xlimit[j,], ylimit[j,], scale =  scale)
    
    cyr <- as.character(yr)
    
    for(k in 1:njyr){
      
      add <- WO <- WP <- WT <- FALSE
      
      if(WOJ){
        if( cyr[k] %in% colnames(obs) &
            MAPTRAPS )WO <- obs[1,colnames(obs) == cyr[k]] > 0
      }
      if(WPJ){
        if( cyr[k] %in% colnames(pred) )WP <- pred[,colnames(pred) == cyr[k]] > 0
      }
      
      if(!WO & !WP & !MAPTREES)next
      
      if(WP){  #predicted surface, fecundity
        
        tmp <- .pmap(specNames=specNames, # xytree=xytree, 
                     plot=mapPlot[j], MAPTREES,
                     year=jyr[k], seedPredGrid=seedPredGrid, 
                     treePredGrid = treePredGrid,
                     xlim=xlimit[j,], ylim=ylimit[j,], treeScale, trapScale,
                     sCol = specCol[specNames], add=add)
        add <- tmp$add
        tmp <- tmp[names(tmp) != 'add']
        
        if(add)colList <- append(colList, list(tmp))
        names(colList)[length(colList)] <- mapPlot[j]
      }
      
      if(WO){  #observed seed
        
        seedk <- seed[seed$year == jyr[k],]
        sx <- seedk$x
        sy <- seedk$y
        z  <- as.matrix( seedk[,seedNames] )
        z  <- rowSums(z, na.rm  = TRUE)
        
        w1 <- which(z > 0)
        w0 <- which(z == 0)
    
        
        if(length(w1) > 0){
          z <- z/seedMax
          z <- 5*sqrt(sc*z)*trapScale
          symbols(sx[w1], sy[w1], squares=z[w1], inches = F,
                  xlab='', ylab='',bg = .getColor('black',.3), 
                  fg = .getColor('black',.5), add=add,
                  xaxt='n', yaxt='n', xlim=xlimit[j,], ylim=ylimit[j,])
          for(i in 1:4)axis(i, labels = FALSE, tck=.01)
          add <- TRUE
        }
        if(add == F){
          plot(NULL, xlim=xlimit[j,], ylim=ylimit[j,], xlab='', ylab='',
               axes = F)
          for(i in 1:4)axis(i, labels = FALSE, tck=.01)
          add <- TRUE
        }
        if(length(w0) > 0)points(sx[w0], sy[w0], pch=3, cex=.3,
                                 col=.getColor('black',.8))
        add <- TRUE
      }
      
      if(MAPTREES & !PREDICT){
        
        treek <-  fmat[drop = FALSE, fmat$year == jyr[k],]
        z     <- treeSymbol[drop = FALSE, fmat$year == jyr[k]]
        if(nrow(treek) == 0 & length(mapTreeYr) >= k){
          treek <-  fmat[drop = FALSE, fmat$year == mapTreeYr[k],]
          z     <- treeSymbol[drop = FALSE, fmat$year == mapTreeYr[k]]
        }
        
        if(nrow(treek) == 0)next
        
        sx <- treek$x
        sy <- treek$y
        
        if(!all(is.na(z))){
          mmm <- max(z, na.rm  = TRUE)
          if(mmm == 0)mmm <- 1
          z <- 1*sc*z/mmm*treeScale
          ic <- match(treek$species, specNames)
          
          symbols(sx, sy, circles=z*1.3, inches = F, add=add,
                  xlim = xlim, ylim = ylim,
                  fg = .getColor('white',.5), xlab = '', ylab = '',
                  bg =.getColor('white',.5), xaxt='n', yaxt='n')
          
          symbols(sx, sy, circles=z, inches = F, add  = TRUE,
                  xlab = '', ylab = '',
                  fg = .getColor(specCol[ specNames[ic] ],.6), 
                  bg=.getColor(specCol[ specNames[ic] ],.3), xaxt='n', yaxt='n')
          if(!add)for(i in 1:4)axis(i, labels = FALSE, tck=.01)
        }
      }
      
      if(!add)next
      
      cex <- sum(mfrow)^(-.1)
      .plotLabel(jyr[k], 'topright',cex = cex, above  = TRUE) 
      .plotLabel(mapPlot[j],'topleft',cex = cex,above  = TRUE)
    }
  }
  if(PREDICT)colList <- colList[ !duplicated(names(colList)) ]
  
  if(SCALEBAR)scaleBar('m', value = scaleValue, yadj=.07, cex=.8)
  if(LEGEND){
    cornerLegend('bottomright', names(specAll),
           text.col = specCol[match(names(specAll),specNames)],
           cex=.7, bty='n')
  }
  
  if(COLORSCALE & PREDICT){
    
    # use last plot (bottom or right side)
    cols <- colList[[length(colList)]]
    nss   <- length(cols$species)
    endLabels <- NULL
    
    for(k in 1:nss){
      
      if(k == nss)endLabels <- c(0, signif(seedMax,1) )
      
      ck <- .getColor(specCol[cols$species[k]], cols$colorLevels)
      
      clist <- list( kplot=k, ytick=NULL, text.col = 'black',
                     cols=ck, labside='right', text.col=col,
                     bg='grey', endLabels=endLabels) 
      cornerScale( clist )
  #    kk <- kk + 1
    }
  }
  invisible(add)
}
 
cornerLegend <- function(...) {
  suppressWarnings(
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
              mar=c(0, 0, 0, 0), new  = TRUE)
  )
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}

cornerScale <- function( clist ) {
  
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
              mar=c(0, 0, 0, 0), new  = TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  
  .cornerLegendScale( clist )
}

.cornerLegendScale <- function( clist ){  
  
  # left and right corners: xx = (x1,x2), y = (y1,y2)
  # bg is color of border
  # cols  - matching color sequence
   kplot <- 1
  
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
              mar=c(0, 0, 0, 0), new  = TRUE)
  on.exit(par(opar))
  
  xx <- yy <- cols <- NULL
  ytick <- scale <- text.col <- text.col <- bg <- endLabels <- NULL
  labside <- 'right'
  
  for(k in 1:length(clist))assign( names(clist)[k], clist[[k]] ) 
  
  xx <- -1.08 + .1*(kplot - 1) + c(.1, .2)
  yy <- c(-1, -.75)
  
  nn <- length(cols)
  ys <- seq(yy[1],yy[2],length=nn)
  if(is.null(scale))scale <- ys
  
  for(j in 1:(length(scale)-1)){
    
    rect(xx[1],ys[j],xx[2],ys[j+1],col=cols[j],border=NA)
  }
  if(!is.null(bg))rect(xx[1],yy[1],xx[2],yy[2],border=bg,lwd=1)
  if(!is.null(ytick)){
    
    ys <- diff(yy)/diff(range(ytick))*ytick
    yt <- ys - min(ys) + yy[1]
    
    for(j in 1:length(yt)){
      lines(xx,yt[c(j,j)])
    }
  }
  if(!is.null(endLabels)){ 
    if(labside == 'right')text(diff(xx)+c(xx[2],xx[2]),yy,endLabels)
    if(labside == 'left')text(c(xx[1],xx[1]),yy,endLabels,pos=2)
  }
}

values2grid <- function(x, y, z, nx=NULL, ny=NULL, dx=NULL, dy=NULL,
                        ksearch = 4, MATFORMAT  = TRUE){
  
  xl <- range(x)
  yl <- range(y)
  
  xs <- seq(xl[1],xl[2],length=nx)
  ys <- seq(yl[1],yl[2],length=ny)
  
  grid <- as.matrix(expand.grid(xs,ys))
  
  tmp <- nn2(cbind(x,y),grid,k=ksearch)
  nn  <- tmp[[1]]
  wt  <- tmp[[2]]
  mn  <- min(wt[wt > 0])
  wt  <- 1/(wt + mn)
  
  zz  <- matrix(z[nn],nrow(nn),ncol(nn))
  zz  <- rowSums(zz*wt)/rowSums(wt)
  
  if(!MATFORMAT)return(  cbind(grid,zz) )
  
  zmat <- matrix(NA, nx, ny)
  ix  <- match(grid[,1],xs)
  iy  <- match(grid[,2],ys)
  zmat[ cbind(ix,iy) ] <- zz
  
  rownames(zmat) <- xs
  colnames(zmat) <- ys
  
  list( x = xs, y = ys, z = zmat )
}

.pmap <- function(specNames=NULL, plot=NULL, MAPTREES  = TRUE,
                  year=NULL, seedPredGrid=NULL,
                  treePredGrid, xlim, ylim, treeScale, trapScale,
                  sCol = 'blue', add = FALSE){
  
  #multiple species for single plot-year
  
  ADD <- FALSE
  if(add)ADD <- TRUE
  
  pnames   <- paste(specNames,'_meanM2',sep='')
  nspec    <- length(specNames)
  predCols <- c('x','y',pnames)
  
  if(length(sCol) == 1 & nspec > 1)sCol <- rep(sCol, nspec)
  
  fec  <- treePredGrid[,'fecEstMu']
  matr <- treePredGrid[,'matrPred']
  smat <- as.matrix(seedPredGrid[,predCols])
  
  nx <- ceiling(diff(xlim)/3)
  ny <- ceiling(diff(ylim)/3)
  
  fecMax  <- max(fec,na.rm  = TRUE)
  seedMax <- max(smat[,pnames],na.rm  = TRUE)
  
  rr  <- apply( rbind(xlim, ylim), 2, range)
  sc  <- max( apply(rr,1,diff) )/20
  
  wspec <- which(colSums(smat[,pnames, drop = FALSE]) > 0)
  
  sn <- 4
  if(seedMax > 1)sn <- 3
  if(seedMax > 10)sn <- 2
  
  q <- seq(0, 1, length=10)^.3
  q[1] <- .3
  
  levels <- signif( quantile(smat[,pnames], q ) , sn)
  levels <- c(levels, signif(max(seedMax)*1.3, sn) )
  levels <- sort(unique(levels))
  colorLevels <- seq(.001, .99, length=length(levels))^3
  
  for(k in wspec){
    
    tmp <- values2grid(x=smat[,'x'],y=smat[,'y'],z=smat[,pnames[k]],
                       nx=nx, ny=ny)
    xseq <- tmp$x
    yseq <- tmp$y
    zmat <- tmp$z
    
    col <- .getColor(sCol[specNames[k]],colorLevels)
    contour(xseq, yseq, zmat, levels=levels, add=add, 
            col=col, labcex = 1, frame.plot = FALSE,
            drawlabels = FALSE, axes = F)
    .filled.contour(xseq, yseq, zmat, levels=levels, col=col)
    if(!add){
      for(m in 1:4)axis(m, labels = FALSE, tck=.03)
    }
    if(!ADD)add <- TRUE
    
  }
  
  if(MAPTREES){
    for(k in wspec){
      wk <- which( treePredGrid$species == specNames[k] )
      if(length(wk) == 0) next
      
      z   <- 5*sc*fec[wk]*treeScale
      z[z > 0] <- z[z > 0]/fecMax
      .mapSpec(x = treePredGrid[wk,'x'],y = treePredGrid[wk,'y'], z*1.5,  
               add = add, 
               mapx = xlim, mapy = ylim,
               colVec='white', fill = 'white')
      .mapSpec(x = treePredGrid[wk,'x'],y = treePredGrid[wk,'y'], z,  add  = TRUE, 
               mapx = xlim, mapy = ylim,
               colVec=sCol[specNames[k]], fill = sCol[specNames[k]])
    }
  }
  invisible( list(species = wspec, colorLevels = colorLevels, add = add) )
}

.updateCovariance <- function(SS, priorSS, n, df){
  
  SI   <- solveRcpp(SS + df*priorSS)
  sinv <- .rwish(n + df, SI)
  solveRcpp(sinv)
}

.updateAlphaRand <- function(ntree, yA, xfecA, sg, reIndexA, reGroups,
                             Arand, priorVA, dfA, specNames, minmax = 3){
  
  # any individual that is mature only one year cannot have random effects
  
  ONEA <- FALSE
  if(length(Arand) == 1)ONEA <- TRUE
  
  arand <- matrix(0, ntree, ncol(Arand))
  
  alphaRand <- randEffectRcpp(gindex=reIndexA, groups = reGroups, 
                              xfecA, yA, sg, solve(Arand))
  alphaRand[alphaRand < -minmax] <- -minmax
  alphaRand[alphaRand > minmax]  <- minmax
  
  if(ONEA){
    mrand     <- mean(alphaRand)
    names(mrand) <- specNames
    mrand[ !is.finite(mrand) ] <- 0
    alphaRand <- alphaRand - mrand
    arand[reGroups,] <- alphaRand
    Arand <- matrix(1/rgamma(1,1 + length(reGroups)/2, 1 + 1/2*sum(alphaRand^2)),1)
  }else{
    mrand <- colSums(alphaRand)/colSums(xfecA)
    mrand[ !is.finite(mrand) ] <- 0
    alphaRand <- alphaRand - matrix(mrand, nrow(alphaRand), ncol(alphaRand), byrow  = TRUE)
    arand[reGroups,] <- alphaRand
    AA    <- crossprod(alphaRand)
    Arand <- .updateCovariance(AA, priorVA, length(reGroups), dfA)
  }
  arand[arand < -minmax] <- -minmax
  arand[arand > minmax]  <- minmax
  
  list(alphaRand = arand, Arand = Arand, meanRand = mrand)
}

.rwish <- function(df,SS){
  z  <- matrix(rnorm(df*nrow(SS)),df,nrow(SS))%*%chol(SS)
  crossprod(z)
}

.riwish <- function(df,S){
  solveRcpp(.rwish(df,solveRcpp(S)))
}

print.mastif <- function(x, ...){
  
  rMu <- rSe <- usigma <- betaYrMu <- betaYrSe <- NULL
  
  cat("\nDIC:\n")
  print( round(x$fit$DIC,0) )
  
  cat("\nFecundity coefficients:\n")
  print( signif(x$parameters$betaFec, 3) )
  
  cat("\nMaturation coefficients:\n")
  print( signif(x$parameters$betaRep, 3) )
  
  if('betaYrMu' %in% names(x$parameters)){
    cat("\nYear effects:\n")
    print( signif(betaYrMu, 3) )
    print( signif(betaYrSe, 3) )
  }
  
  if('rgibbs' %in% names(x$chains)){
    cat("\nSpecies to seed type matrix R:\n")
    print( signif(rMu, 3) )
    print( signif(rSe, 3) )
  }
  
  cat("\nSigma, RMSPE:\n")
  usigma <- x$parameters$sigma
  print( signif(usigma, 4) )
  
  cat("\nDispersal parameter u (m^2):\n")
  print( x$parameters$upars )
  
  cat("\nKernel mean distance (m):\n")
  print( x$parameters$dpars )
}

summary.mastif <- function(object, verbose  = TRUE, latex = FALSE, ...){ 
  
  SEEDDATA <- TRUE
  
  betaFec   <- object$parameters$betaFec
  betaRep   <- object$parameters$betaRep
  priorTable <- object$inputs$priorTable
  seedNames <- object$inputs$seedNames
  specNames <- object$inputs$specNames
  tdata     <- object$inputs$treeData
  sdata     <- object$inputs$seedData
  plots     <- sort(unique(as.character(tdata$plot)))
  ntype     <- length(seedNames)
  nseed     <- nrow(sdata)
  nplot     <- length(plots)
  words     <- character(0)
  

  
  if(is.null(sdata))SEEDDATA <- FALSE
  
  if(latex)verbose <- FALSE
  
  out <- list()

  AR <- YR <- RANDOM <- SAMPR <- FALSE
  
  trapRMSPE <- object$fit$RMSPE
  trapDIC   <- object$fit$DIC
  cropRMSPE <- object$fit$rmspeCrop
  
  if("arList" %in% names(object$data)){
    AR <- TRUE
    plag <- object$data$arList$p
  }
  
  #model
  model <- rbind( as.character(object$inputs$formulaRep), 
                  as.character(object$inputs$formulaFec ))[, 2, drop = FALSE]
  model <- .replaceString(model, 'I(', '')
  model <- .replaceString(model, '^2)', '^2')
  rownames(model) <- c('Maturation:', 'Fecundity:')
  attr(model, 'caption') <- 'Model terms '
  re <- object$inputs$randomEffect
  if(!is.null(re)){
    ref <- paste( re$randGroups, re$formulaRan, sep='_' )[2]
    ref <- paste( ' + r(', ref, ')', sep='')
    model[2,1] <- paste(model[2,1], ref)
    attr(model, 'caption') <- paste(attr(model,'caption'), 
                                    'r(g_m) indicates a random effect for groups g and model m. ', sep='')
  }
  ye <- object$inputs$yearEffect
  if(!is.null(ye)){
    yef <- unlist(ye)
    if(length(yef) > 1)yef <- paste0(yef, collapse = '_')
    if(!AR){
      yef <- paste( ' + y(', yef, ')', sep='')
      model[2,1] <- paste(model[2,1], yef)
      attr(model, 'caption') <- 
        paste(attr(model,'caption'), 
              'y(h) indicates a year effect for random h. ', sep='')
    }
    if(AR){
      yef <- paste( ' + AR(',plag,', ', yef, ')', sep='')
      model[2,1] <- paste(model[2,1], yef)
      attr(model, 'caption') <- 
        paste(attr(model,'caption'), 
              'AR(p, h) indicates an AR(p) model for random group h. ', sep='')
    }
  }
  model <- .replaceString(model, ' * ', 'X')
  model <- .replaceString(model, '*', 'X')
  model <- .replaceString(model, '  ', ' ')
  model[1,1] <- paste("'", model[1,1], "'", sep='')
  model[2,1] <- paste("'", model[2,1], "'", sep='')
  colnames(model) <- 'model'
  out$amodel <- model
  
  
  prNames <- c("priorU", "priorVU", "minU", "maxU", "minDiam","maxDiam","maxFec")
  
  out$aprior <- priorTable[,prNames]
  attr(out$aprior, 'caption') <- 'Prior parameter values'
  
  
  #data summary
  wd <- which(!duplicated(tdata$treeID))
  trees <- table(tdata$species[wd],tdata$plot[wd])[,plots,drop = FALSE]
  rownames(trees) <- paste('trees',rownames(trees), sep='_')
  
  treeYears <- table(tdata$species,tdata$plot)[,plots,drop = FALSE]
  rownames(treeYears) <- paste('tree-yrs',rownames(treeYears), sep='_')
  dataTab <- rbind(trees,treeYears)
  
  if(SEEDDATA){
    wd <- which(!duplicated(sdata$trapID))
    traps <- table(sdata$plot[wd])[plots]
    ntr   <- names(traps)
    traps <- matrix(traps,1)
    colnames(traps) <- ntr
    rownames(traps) <- 'traps'
    
    trapYears <- table(sdata$plot)[plots]
    dataTab <- rbind(dataTab,traps,trapYears)
    rownames(dataTab)[nrow(dataTab)] <- 'trap-yrs'
    
    totalSeed <- buildSeedByPlot(sdata, seedNames)
    censSeed  <- buildSeedByPlot(sdata, paste(seedNames,'_min',sep=''))
    
    ww <- which(!colnames(dataTab) %in% colnames(totalSeed))
    if(length(ww) > 0){
      mm <- matrix(NA, nrow(totalSeed),ncol=length(ww))
      colnames(mm) <- colnames(dataTab)[ww]
      totalSeed <- cbind(totalSeed, mm)
    }
    
    total   <- totalSeed[,colnames(dataTab),drop = FALSE]
    dataTab <- rbind(dataTab,total)
    
    if(nplot > 1){
      total   <- rowSums(dataTab)
      dataTab <- cbind(dataTab,total)
    }
  }
  
  rownames(betaFec) <- .replaceString(rownames(betaFec), 'species','')
  rownames(betaRep) <- .replaceString(rownames(betaRep), 'species','')
  
  
  attr(dataTab, 'caption') <- 'Summary of observations by plot'
  attr(betaFec, 'caption') <- 'Fecundity coefficients (betaFec, log scale)'
  attr(betaRep, 'caption') <- 'Maturation coefficients (betaRep, probit)'
  
  out$adata   <- dataTab
  out$betaFec <- betaFec[,1:4]
  out$betaRep <- betaRep[,1:4]
  attr(out$betaFec, 'caption') <- 'Fecundity coefficients (betaFec, log scale)'
  attr(out$betaRep, 'caption') <- 'Maturation coefficients (betaRep, probit)'
  
  if(verbose){
    cat('\nData summary:\n')
    print(dataTab)
    
    cat("\nFecundity parameters (log scale):\n")
    print( betaFec )
    
    cat("\nMaturation parameters (probit):\n")
    print( betaRep )
  }
  
  if('betaYrMu' %in% names(object$parameters)){
    YR <- TRUE
    byr <- object$parameters$betaYr[,1:3]
    attr(byr, 'caption') <- 'Year effects for fecundity (log scale)'
    
    if(verbose){
      cat("\nYear effects, only mature individuals:\n")
      print( signif(byr, 4) )
    }
    out$betaYr <- signif(byr, 4)
  }
  
  if('betaYrRand' %in% names(object$parameters)){
    
    byrRand <- signif(object$parameters$betaYrRand, 4)
    byrRSE  <- signif(object$parameters$betaYrRandSE, 4)
    attr(byrRand, 'caption') <- 'Year random effects for fecundity by random group'
    attr(byrRSE, 'caption')  <- 'Year standard error for fecundity by random group'
    
    if(sum(byrRand) != 0){
      
      if(verbose){
        cat("\nYear effects, random group means:\n")
        print( byrRand )
        cat("\nYear effects, standard deviation between groups:\n")
        print( byrRSE )
      }
      
      out$betaYrRand   <- byrRand
      out$betaYrRandSE <- byrRSE
    }
  }
  
  if( !is.matrix(object$parameters$pacfMat) )object$parameters$pacfMat <- 
    t( as.matrix( object$parameters$pacfMat ) )
 # if( !is.matrix(object$parameters$pacsMat) )object$parameters$pacsMat <- 
 #   t( as.matrix( object$parameters$pacsMat ) )
  if( !is.matrix(object$parameters$pacfSe) )object$parameters$pacfSe <- 
    t( as.matrix( object$parameters$pacfSe ) )
  
  pacfmu <- signif( object$parameters$pacfMat[,-1,drop = FALSE], 4 )
  pacfse <- signif( object$parameters$pacfSe[,-1,drop = FALSE], 4 )
  
  if(ncol(pacfmu) > 8)pacfmu <- pacfmu[,1:8]
  if(ncol(pacfse) > 8)pacfse <- pacfse[,1:8]
  
  attr(pacfmu, 'caption') <- 'Partial autocorrelation in fecundity (log scale)'
  attr(pacfse, 'caption') <- 'Partial autocorrelation standard error in fecundity'
 
  out$pacfmu <- pacfmu
  out$pacfse <- pacfse
  
  
  if(SEEDDATA){
    pmat <- object$parameters$pacsMat
    if(!is.null(pmat)){
      pacfss <- signif( pmat[-1,drop = FALSE], 4 )
      if(length(pacfss) > 8)pacfss <- pacfss[1:8]
      attr(pacfss, 'caption') <- 'Partial autocorrelation in seed rain, all plots'
      out$pacfss <- pacfss
    }
  }
  
  
  if('aMu' %in% names(object$parameters)){
    RANDOM <- TRUE
    amu <- diag( object$parameters$aMu)
    ase <- diag( object$parameters$aSe )
    
    arand <- signif(cbind(amu,ase),3)
    rownames(arand) <- .replaceString(rownames(arand),'species','')
    colnames(arand) <- c('  Estimate','  Standard_error')
    
    if(verbose){
      cat("\nDiagonal elements of random effects covariance matrix (aMu):\n")
      print( arand )
    }
    attr(arand, 'caption') <- 'Diagonal elements of random effects covariance matrix (A)'
    out$arand <- arand
  }
  
  if('rgibbs' %in% names(object$chains)){
    SAMPR <- TRUE
    
    rMu <- signif(object$parameters$rMu, 4)
    rSe <- signif(object$parameters$rSe, 4)
    
    attr(rMu,'species') <- attr(rSe,'species') <- 
      attr(rMu,'posR') <- attr(rSe,'posR') <- 
      attr(rMu,'plot') <- attr(rSe,'plot') <- NULL
    
    if(verbose){
      cat("\nSpecies to seed type matrix R:\n")
      print( rMu )
      
      cat("\nStandard errors for R:\n")
      print( rSe )
    }
    
    attr(rMu, 'caption') <- 'ID error matrix mean estimate (R)'
    attr(rSe, 'caption') <- 'ID error matrix standard error (R)'
    
    out$rMu <- rMu
    out$rSe <- rSe
  }
  
  sigma <- object$parameters$sigma
  if(!SEEDDATA)sigma <- sigma[drop=FALSE,'sigma',]
  
  if(verbose){
    cat("\nSigma, RMSPE:\n")
    print( signif(sigma, 3) )
  }
  attr(sigma, 'caption') <- 'Sigma^2, RMSPE, and deviance'
  out$sigma <- signif(sigma, 3)
  
  if(SEEDDATA){
    utab <- object$parameters$upars
    
    upars <- utab[grep('u',rownames(utab)),]
    pu    <- rownames(utab)
    rownames(utab) <- NULL
    utab  <- data.frame(parameter = pu, utab)
    
    rownames(utab) <- NULL
    
    if(verbose){
      cat("\nKernel estimates:\n")
      print(utab)
    }
    attr(utab, 'caption') <- 'Kernel estimates, parameter u in 2Dt kernel and kernel mean d'
  }
  
  if(AR){
    out$eigenMu <- signif(object$parameters$eigenMu, 3)
    out$eigenSe <- signif(object$parameters$eigenSe, 3)
    attr(out$eigenMu, 'caption') <- 'Eigenvalues for autoregressive terms'
    attr(out$eigenSe, 'caption') <- 'Eigenvalue standard errors for autoregressive terms'
  }
  
  ty <- dataTab[grep('tree-yrs',rownames(dataTab)),]
  
  ntreeYr  <- sum(ty)
  years    <- object$data$setupData$years
  ntree    <- nrow( object$data$setupData$zmat )
  plots    <- object$data$setupData$plots
  nplot    <- length(plots)
  nyr      <- length(years)
  
  words <- paste("The sample contains ", ntreeYr, " tree-years on ", ntree, 
                 " individuals. There are ", nplot, " plots sampled over ", nyr, 
                 " years. ", sep="")
  if(SEEDDATA){
    ntrapYr  <- sum(dataTab['trap-yrs',])
    ntrap    <- nrow(object$data$setupData$distall)
    words <- paste(words, 'There are ', ntrapYr, ' trap-years on', ntrap, 
                 ', seed traps. The RMSPE is ',
                 signif(object$fit$RMSPE,3),
                 ", and the DIC is ",round(object$fit$DIC),". ",
                   sep='')
  }
  if('rmspeCrop' %in% names(object$fit)){
    ww <- which(is.finite(object$inputs$treeData$cropCount))
    ctab <- table(object$inputs$treeData$species[ww])
    cnames <- paste0( names(ctab), collapse=', ')
    ctot   <- sum(ctab)
    words <- paste(words, 'There are ',ctot,' cropCounts on ',cnames,'. ',
                   sep='')
  }
  
  if('trueValues' %in% names(object$inputs)){
    words <- paste('The data are generated by mastSim. ', words, sep='')
  }
  
  if('censMin' %in% names(object$inputs)){
    cmin <- object$inputs$censMin
    cmax <- object$inputs$censMax
    if(length(cmin) > 0){
      ncens <- nrow(cmin)
      words <- paste(words, paste("There are", ncens, 
                              "censored seed trap observations") )
    }
  }
  
  if(!is.null(object$inputs$inwords)){
    ww <- paste0(object$inputs$inwords, collapse='. ')
    ww <- .replaceString(ww, '\n','. ')
    ww <- .replaceString(ww, '\n','')
    ww <- .replaceString(ww, '\"','')
    ww <- .replaceString(ww, '. .','. ')
    ww <- .replaceString(ww, ':.',': ')
    words <- paste( words, ww, sep='. ')
    words <- .replaceString( words, '. .', '.')
  }
  
  if(RANDOM){
    
    morewords <- paste("Random effects were fitted on ",
                            length(object$data$setupRandom$rnGroups), 
                            " individuals.", sep='')
    words <- paste(words, morewords)
  }
  
  if(verbose){
    cat("\n",words)
  }
  
  fit <- numeric(0)
  if('rmspeCrop' %in% names(object$fit)){
    fit <- matrix(object$fit$rmspeCrop, 1, 1)
    colnames(fit) <- 'cropRMSPE'
  }
  if(SEEDDATA){
    sfit <- cbind(round(object$fit$DIC), object$fit$RMSPE)
    colnames(sfit) <- c('  DIC','  RMSPE')
    fit <- cbind(fit,sfit)
  }
  attr(fit, 'caption') <- 'Model fit'
  
  out$fit   <- fit
  out$words <- words
  
  out <- out[order(names(out))]
  
  if(latex){
    
    tnames <- names(out)
    for(k in 1:length(tnames)){
      xk <- out[tnames[k]][[1]]
      if(tnames[k] %in% c('fit', 'words'))next
      print('')
      print(tnames[k])
      xp <- as.data.frame(xk)
      ww <- which(sapply(xp, is.complex))
      if(length(ww) > 0)for(w in ww)xp[[w]] <- as.character(xp[[w]])
      print( xtable( xp, caption=attr(xk,'caption') ) )
    }
  }
 
  class(out) <- "summary.mastif"
  invisible(out) 
}

mastSim <- function(sim){       # setup and simulate fecundity data
  
  if(!'seedNames' %in% names(sim))
    stop('\nsim must include seedNames\n')
  if(!'specNames' %in% names(sim))
    stop('\nsim must include specNames\n')
  
  .dsim( sim ) 
}
 
.dsim <- function(sim){
  
  nyr <- 5; ntree  <-  10; ntrap <- 20; plotWide  <-  100; nplot  <-  3
  meanDist  <-  25; Q  <-  2
  minDiam <- 10
  maxDiam <- 40
  yearEffect <- NULL
  facLevels <- character(0)
  seedNames <- specNames <- c('piceaGlauca','piceaMariana','piceaUNKN')
  SPECS <- SAMPR <- AR <- USPEC <- FALSE
  maxFec <- 1e+7
  
  for(k in 1:length(sim))assign( names(sim)[k], sim[[k]] )
  specNames <- sort(specNames)
  S         <- length(specNames)
  ntreePlot <- rpois(nplot, ntree) + 1  # trees per plot
  ntrapPlot <- rpois(nplot, ntrap) + 4  # traps per plot
  nyrPlot   <- rpois(nplot, nyr) + 1
  plotNames <- paste('p',1:nplot,sep='')
  yearNames <- 2017 - c(max(nyrPlot):1)
  if(length(specNames) > 1)SPECS <- TRUE
  if(length(seedNames) > 1)SAMPR <- TRUE
  nyr <- max(nyrPlot)
  
  upar <- (2*meanDist/pi)^2
  
  year <- plot <- tree <- trap <- yrsd <- plsd <- numeric(0)
  
  for(j in 1:nplot){
    
    tree <- c( tree, rep(1:ntreePlot[j], each=nyrPlot[j]) )
    year <- c( year, rep(1:nyrPlot[j], ntreePlot[j]) ) 
    plot <- c( plot, rep(j, (nyrPlot[j]*ntreePlot[j]) ) )
    trap <- c( trap, rep(1:ntrapPlot[j], each=nyrPlot[j]) )
    yrsd <- c( yrsd, rep(1:nyrPlot[j], ntrapPlot[j]) ) 
    plsd <- c( plsd, rep(j, (nyrPlot[j]*ntrapPlot[j]) ) )
  }
  
  tree <- data.frame(plot = plotNames[plot], year = yearNames[year], tree = tree )
  tree$plot <- as.character(tree$plot)
  tree$tree <- as.character(tree$tree)
  tree <- tree[order(tree$plot, tree$tree, tree$year),]
  
  tree$treeID  <- columnPaste( tree$plot, tree$tree )
  
  id           <- unique(as.character(tree$treeID))
  species      <- sample(specNames, length(id), replace  = TRUE) 
  tree$species <- factor(species[match(tree$treeID,id)], levels = specNames)
  
  tree$dcol <- match(as.character(tree$treeID),id)
  
  tree$plotYr <- columnPaste(tree$plot, tree$year,'_' )
  plotyr      <- unique(tree$plotYr)
  tree$plotyr <- match(tree$plotYr, plotyr )
  
  years <- sort(unique(tree$year))
  
  trap <- data.frame(plot = plotNames[plsd], year = yearNames[yrsd], trap )
  trap$plot <- as.character(trap$plot)
  trap$trap <- as.character(trap$trap)
  trap <- trap[order(trap$plot, trap$trap, trap$year),]
  
  trap$trapID  <- columnPaste( trap$plot, trap$trap )
  trapid       <- as.character( trap$trapID )
  drow         <- unique(trapid)
  trap$drow    <- match(trapid,drow)

  trap$plotYr <- columnPaste(trap$plot, trap$year,'_' )
  plotyr      <- unique(trap$plotYr)
  trap$plotyr <- match(trap$plotYr, plotyr )
  
  n <- nrow(tree)
  
  xfec <- round( matrix( .tnorm(n*(Q-1), 5, 50, 35, 5), n, (Q-1) ), 3)
  xnames <- paste('x',1:(Q-1),sep='')
  xnames[1] <- 'diam'
  colnames(xfec) <-  xnames
 
  xdata   <- data.frame( species = tree$species, xfec)
  
  formulaFec <- formula( ~ diam )
  formulaRep <- formula( ~ diam )
  
  if(SPECS){
    formulaFec <- formula( ~ species*diam  ) 
    formulaRep <- formula( ~ species*diam )
  }
  xfec    <- model.matrix(formulaFec, xdata)
  Qf      <- ncol(xfec)
  xrep    <- model.matrix(formulaRep, xdata)
  Qr      <- ncol(xrep)
  
  if(!SPECS){
    ss     <- paste('species',specNames,sep='')
    xnames <- paste(ss,colnames(xfec),sep=':')
    xnames[1] <- .replaceString(xnames[1],':(Intercept)','')
    colnames(xfec) <- xnames
  }
    
  xytree <- xytrap <- distall <- numeric(0)
  
  for(j in 1:nplot){
    
    xy1 <- matrix( runif(2*ntreePlot[j],0,plotWide), ntreePlot[j],2)  
    xy2 <- matrix( runif(2*ntrapPlot[j],0,plotWide), ntrapPlot[j],2)
    xy1 <- round(xy1,1)
    xy2 <- round(xy2,1)
    
    rownames(xy1) <- columnPaste(rep(plotNames[j],ntreePlot[j]), 
                                 c(1:ntreePlot[j]), '-')
    rownames(xy2) <- columnPaste(rep(plotNames[j],ntrapPlot[j]), 
                                 c(1:ntrapPlot[j]), '-')
    xytree  <- rbind(xytree, xy1)
    xytrap  <- rbind(xytrap, xy2)
  }
  
  xdata <- xdata[,!colnames(xdata) %in% colnames(tree),drop = FALSE]
  
  treeData <- cbind(tree,xdata)
  count    <- matrix(0, nrow(trap), length(seedNames))
  colnames(count) <- seedNames
  seedData <- data.frame(trap, count)
  
  colnames(xytree) <- colnames(xytrap) <- c('x','y')
  xy <- columnSplit(rownames(xytree),'-')
  colnames(xy) <- c('plot','tree')
  xytree <- data.frame( xy, xytree )
  xytree$tree <- as.character(xytree$tree)
  wws <- match(rownames(xytree), 
               as.character(treeData$treeID))
  xytree$species <- as.character( treeData$species[ wws ] )
  
  xy <- columnSplit(rownames(xytrap),'-')
  colnames(xy) <- c('plot','trap')
  xytrap <- data.frame( xy, xytrap )
  xytrap$trap <- as.character(xytrap$trap)
  
  xytree$treeID <- rownames(xytree)
  xytrap$trapID <- rownames(xytrap)
  
  seedData$active <- 1
  seedData$area   <- 1
  
  ntree    <- nrow(xytree)
  nyr      <- max(nyrPlot)
  dmat     <- matrix(runif(ntree*nyr, .2, .5), ntree, nyr )
  
  dmat[,1] <- .tnorm(ntree, 1, 70, 35, 40) 
  small    <- sample(ntree,round(ntree/3))
  dmat[small,1] <- .tnorm(length(small), 1, 40, 7, 20)
  if(nyr > 1)dmat <- round(t( apply(dmat,1,cumsum) ), 2)
  
  
  tyindex  <- cbind(treeData$dcol, match(treeData$year,years))
  treeData$diam <- dmat[ tyindex ]
  
  treeData$obs <- 1
  seedData$obs <- 1
  
  plotTreeYr <- columnPaste(as.character(treeData$treeID), 
                            as.character(tree$year), '_')
  rownames(treeData) <- plotTreeYr
  rownames(seedData) <- columnPaste(seedData$trapID, 
                                    as.character(seedData$year), '_')
  
  tmp <- .setupData(formulaFec, formulaRep, tdata = treeData, sdata = seedData, 
                    xytree, xytrap, specNames, seedNames, 
                    AR = FALSE, YR = FALSE, yearEffect, minDiam, maxDiam, TREESONLY = FALSE, 
                    notFit = NULL, maxF=maxFec, CONES = FALSE)
  
  treeData  <- tmp$tdata
  seedData  <- tmp$sdata
  distall   <- tmp$distall
  zmat      <- tmp$zmat
  zknown    <- tmp$zknown
  xytree    <- tmp$xytree
  xytrap    <- tmp$xytrap
  plotNames <- tmp$plotNames
  plots     <- tmp$plots
  years     <- tmp$years
  xfec      <- tmp$xfec
  xrep      <- tmp$xrep
  nseed     <- nrow(seedData)
  scode     <- tmp$scode
  nplot     <- length(plotNames)
  n         <- nrow(xfec)
  ntobs     <- table(treeData$plot) 
  nsobs     <- table(seedData$plot)
  ttab      <- table(treeData$plot, treeData$year)
  wtab      <- which(ttab > 0, arr.ind  = TRUE) 
  xfecMiss <- tmp$xfecMiss
  xrepMiss <- tmp$xrepMiss
  xfecCols <- tmp$xfecCols
  xrepCols <- tmp$xrepCols
  ntree    <- nrow(xytree)
  ntrap    <- nrow(xytrap)
  xfecU    <- tmp$xfecU; xfecT <- tmp$xfecT
  xrepU    <- tmp$xrepU; xrepT <- tmp$xrepT
  xfecs2u  <- tmp$xfecs2u
  xfecu2s  <- tmp$xfecu2s
  xreps2u  <- tmp$xreps2u
  xrepu2s  <- tmp$xrepu2s
  nspec    <- length(specNames)
  matYr    <- tmp$matYr
  last0first1 <- tmp$last0first1
  
  tyindex  <- cbind(as.character(treeData$treeID), as.character(treeData$year))
  
  zknownVec <- zknown[ tyindex ]
  z <- zmat[ tyindex ]
  
  
  dmin <- dmax <- numeric(0)
  for(k in 1:nspec){
    wk <- which(treeData$species == specNames[k])
    wm <- which.min( (treeData$diam[wk] - minDiam)^2 )
    wf <- xfec[wk[wm],]
    dmin <- c(dmin, wf[!wf %in% c(0,1)][1])
    
    wm <- which.min( (treeData$diam[wk] - maxDiam)^2 )
    wf <- xfec[wk[wm],]
    dmax <- c(dmax, wf[!wf %in% c(0,1)][1])
  }
  
  
  
   dcols <- grep('diam',colnames(xfec))
   dd    <- rowSums(xfec[,dcols, drop = FALSE])
   
   wcol  <- grep('diam',colnames(xfec)) 
   slope <- log(maxFec)/(max(xfec[,wcol]) - min(xfec[,wcol]))
   slope <- .tnorm(nspec, 1, 2, 1.5, 1) + slope
   int   <- -slope*dmin + .tnorm(nspec, 3, 7, 5, 1)
   betaFec <- matrix( c(int, slope), ncol=1)  # 2, 4
   rownames(betaFec) <- colnames(xfec)
   
   fec <- xfec%*%betaFec
   
   ######################## propose z 
   tmp <- .propZ(zmat, last0first1, matYr)
   zmat  <- tmp$zmat
   matYr <- tmp$matYr
   z <- zmat[tyindex]

   
   species <- as.factor(treeData$species)
   
   form <- formula( z ~ treeData$diam) 
   if(nspec > 1) form <- formula( z ~ species*treeData$diam)
   
   br <- suppressWarnings(
     glm(form, family=binomial("probit") )$coefficients
   )
   
  names(br) <- .replaceString(names(br),'treeData$','')
  names(br) <- .replaceString(names(br),'species','')
  slopes <- grep('diam',names(br))
  ints   <- which( !c(1:length(br)) %in% slopes)
  
  rspec <- specNames[ which(!specNames %in% names(br)) ]
  
  ints <- br[ints]
  ints[-1] <- ints[-1] + br[ '(Intercept)' ]
  names(ints)[1] <- rspec
  slopes <- br[slopes]
  slopes[-1] <- slopes[-1] + br['diam']
  names(slopes)[1] <- paste(rspec,':diam',sep='')
  ints <- ints[sort(names(ints))]
  slopes <- slopes[sort(names(slopes))]
  btmp <- matrix(c(ints,slopes),ncol=1)                # unstandardized diam
  rownames(btmp) <- c(names(ints),names(slopes))
  
  betaRep <- xrepu2s%*%btmp                            # standardized diam
  
  ztrue <- z
  q <- which(ztrue == 1)
  
  hi  <- 0*fec + log(maxFec)
  lo <- -hi/3
  hi[z == 0] <- 0
  lo[z ==1] <- 0
  
  for(j in 1:10){
    fec <- .tnorm(nrow(xfec), lo, hi, xfec%*%betaFec, .01)
    betaFec <- solve( crossprod(xfec) )%*%crossprod(xfec,fec)
  }
  
 
  zmat[ sample(n,n/20) ] <- NA
  
  treeData$repr <- zmat[ tyindex ]

  seedData$active <- seedData$area <- 1
  
  fec <- exp(fec)
  
  #unstandardized diam
  bfecSave <- xfecs2u%*%betaFec      
  brepSave <- xreps2u%*%betaRep

  # R matrix
  
  fill <- 0
  if(length(seedNames) == 1)fill <- 1
  
  R <- matrix(fill, length(plots)*nspec, length(seedNames))
  colnames(R) <- seedNames
  rr <- as.vector( outer(specNames, plots, paste, sep='-') )
  rownames(R) <- rr
    
  wun <- grep('UNKN', seedNames)
  if(length(wun) > 0){
    kk <- c(1:length(seedNames))[-wun]
    for(k in kk){
      wsp <- grep(seedNames[k],rownames(R))
      R[wsp,seedNames[k]] <- 1
    }
    R[,wun] <- 2
    R <- sweep(R, 1, rowSums(R), '/')
  }else{
    sr <- columnSplit(rr, '-')[,1]
    ir <- match(sr, seedNames)
    R[ cbind(1:length(rr), ir) ] <- 1
  }
  tmp <- columnSplit(rownames(R),'-')
  attr(R, 'species') <- tmp[,1]
  attr(R, 'plot')    <- tmp[,2]
  
  treeData$specPlot <- columnPaste(treeData$species,treeData$plot)
  
  obsRows <- which(treeData$fit == 1)
  
  
  lambda <- .getLambda(tdat1 = treeData[obsRows,c('specPlot','year','plotyr','dcol')], 
                       sdat1 = seedData[,c('year','plotyr','drow')], 
                       AA = seedData$area, ug = upar, 
                       fz = fec[obsRows]*z[obsRows], R,
                       SAMPR, USPEC, 
                       distance = distall, yrs = years, PERAREA = FALSE) 
  rownames(lambda) <- rownames(seedData)
  lambda <- lambda + 1e-12
  ss     <- matrix(rpois(length(lambda),lambda),nrow(lambda), ncol(lambda))
  seedData[,seedNames] <-   ss
  seedData$active <- 1
  seedData$area   <- 1
  
  seedData$active <- 1
  seedData$area   <- 1
  
  stab <- with( seedData, table(plot, year) )
  ttab <- with( treeData, table(plot, year) )
  sc   <- colSums(stab)
  stab <- stab[,sc > 0, drop = FALSE]
  ttab <- ttab[,sc > 0, drop = FALSE]
  
  form <- as.character( formulaFec )
  form <- .replaceString( form, 'species *', '')
  formulaFec <- as.formula( paste(form, collapse = ' ') )
  
  form <- as.character( formulaRep )
  form <- .replaceString( form, 'species *', '')
  formulaRep <- as.formula( paste(form, collapse = ' ') )
  
  names(fec) <- names(ztrue) <- rownames(xfec)
  
  trueValues <- list(fec = fec, repr = ztrue, betaFecStnd = betaFec,
                     betaRepStnd = betaRep, betaFec = bfecSave, 
                     betaRep = brepSave,upar = upar, R = R)
  
  rownames(treeData) <- columnPaste(treeData$treeID,treeData$year,'_')
  rownames(seedData) <- columnPaste(seedData$trapID,seedData$year,'_')
  
  treeData <- treeData[,c('plot','tree','year','species','diam','repr','repMu')]
  seedData <- seedData[,c('plot','trap','year','area','active',seedNames)]
  xytree   <- xytree[,c('plot','tree','x','y')]
  xytrap   <- xytrap[,c('plot','trap','x','y')]
  
  
  out <- list(trueValues = trueValues, treeData = treeData, seedData = seedData, 
       distall = distall, xytree = xytree, xytrap = xytrap, formulaFec = formulaFec,
       formulaRep = formulaRep, plots = plots, years = years,
       sim = sim,seedNames = seedNames, specNames = specNames, R = R)
  orr <- order(names(out))
  out <- out[orr]
  out
}
      
.seedFormat <- function(sfile, lfile, trapFile = NULL, seedNames = NULL, 
                        specNames, genusName = NULL, omitNames = NULL, plot,
                        newplot = plot, trapID = 'trap', monthYr = 7, active = 1,
                        area = .5, verbose = FALSE ) {
  
  # always include specNames due to ambiguous substr(genusName, 1, 4)
  
  if(is.null(newplot))newplot <- plot
  
  scols      <- c('site','plot','trap','trapID','trapnum','X','Y',
                  'month','day','year',
                  'UTMx','UTMy')
  
  hcols <- c('plot','trap','basket','month','day','year',
             'site','trapID','trapNum','trapName','X','Y','UTMx','UTMy')
  
  midCD <- c( 273890.6, 3938623.3 )  #plot center for (x,y) at GSNP_CD
  
  loc  <- read.csv(lfile, stringsAsFactors = F)
  
  loc  <- loc[loc$plot == plot,]
    xy <- loc[,c('UTMx','UTMy')]
    xy <- round(sweep(xy,2,colMeans(xy),'-'),1)
    loc$x <- xy[,1]
    loc$y <- xy[,2]
 
  ww <- which( is.finite(loc[,'UTMx']) & is.finite(loc[,'UTMy']) )
  
  if(length(ww) == 0){
    if(verbose){
      print('plot without seed trap locations:')
      print(lfile)
    }
    return( numeric(0) )
  }
  
  loc  <- loc[ww,]
  pcol <- rep(plot, nrow(loc))
  id   <- apply( cbind(plot, loc[,trapID]), 1, paste0, collapse='-')
  loc  <- data.frame(trapID = id, trap = loc[,trapID], 
                     loc[,!colnames(loc) == trapID])
  loc$plot <- pcol

  if(plot == "GSNP_CD" | plot == "GRSM_CD"){
    loc$x <- loc$UTMx - midCD[1]
    loc$y <- loc$UTMy - midCD[2]
  }
  
  counts <- read.csv(sfile, stringsAsFactors = F)
  counts <- counts[counts$plot == plot,]
  
  if('area' %in% colnames(counts))area <- counts$area[1]
    
  
  counts[is.na(counts$month),'month'] <- 3
  counts[counts[,'month'] < monthYr,'year'] <- 
    counts[counts[,'month'] < monthYr,'year'] - 1
  
  #all NA
  
  mcols <- c(hcols, 'Notes','notes')
  dmat <- as.matrix( counts[,!colnames(counts) %in% mcols] )
  dmat[is.finite(dmat)] <- 1
  miss <- which(rowSums(dmat, na.rm  = TRUE) == 0)
  
  
  specNames <- c( specNames, paste( substr(gcode, 1, 4), 'UNKN', sep='') )  
  seedNames <- sort(unique(c(seedNames, specNames)))
  
  
  sj <- numeric(0)
  for(k in 1:length(specNames)){
    sj <- cbind( sj, as.matrix( counts[,grep(specNames[k], 
                                             colnames(counts)), drop = FALSE]) )
  }
  seedNames <- colnames(sj)
  # }
  
  if(!is.null(omitNames)){
    sj <- sj[,!colnames(sj) %in% omitNames, drop = FALSE]
    seedNames <- colnames(sj)
  }
  
  if(length(sj) == 0){
    sj <- matrix(0,nrow(counts),1)
    seedNames <- colnames(sj) <- paste( substr(genusName, 1, 4) ,'UNKN',sep='')
  }
  
  
  yr <- sort(unique(counts[,'year']))
  tn <- sort(unique(counts[,trapID]))
  jj <- match(counts[,'year'],yr)
  ii <- match(counts[,trapID],tn)
  smat <- matrix(0, length(tn), length(yr))
  rownames(smat) <- tn
  colnames(smat) <- yr
  
  seedj <- numeric(0)
  
  for(k in 1:ncol(sj)){
    
    smat <- tmat <- smat*0
    
    ck <- sj[,k]
    
    # seed counts
    
    wk <- which(is.finite(ck))
    ck <- ck[wk]
    ik <- tn[ii[wk]]
    jk <- jj[wk]
    ky <- tapply(ck, list(trap = ik, year = jk), sum, na.rm  = TRUE)
    colnames(ky) <- yr[ as.numeric(colnames(ky)) ]
    smat[rownames(ky), colnames(ky)] <- ky
    ky <- smat
    
    # missing values
    ck <- nk <- sj[,k]*0+1
    wk <- which(is.na(ck))
    ck[wk] <- 0
    nk[wk] <- 1
    ik <- tn[ii]
    jk <- jj
    ny <- tapply(ck, list(ik, jk), sum, na.rm  = TRUE) #active intervals
    my <- tapply(nk, list(ik, jk), sum, na.rm  = TRUE) #total intervals
    colnames(ny) <-  colnames(my) <- yr[ as.numeric(colnames(ny)) ]
    smat[rownames(ny), colnames(ny)] <- ny
    tmat[rownames(my), colnames(my)] <- my
    
    ny <- smat
    my <- tmat
    
    colnames(ny) <- yr
    rownames(ny) <- tn
    
    seedj <- cbind(seedj, as.vector(ky))
    
    if(k == 1){
      active <- round(ny/my, 2)
      active[!is.finite(active)] <- 0
    }
  }
  
  colnames(seedj) <- colnames(sj)
  seed <- matrix(0, nrow(seedj), length(seedNames))
  colnames(seed) <- seedNames
  
  active <- as.vector(active)
  
  
  if(!is.null(trapFile)){
    taa <- read.table(trapFile,header  = TRUE)[,c('plot','seedArea')]
    area <- taa[taa$plot == plot,'seedArea']
    if(length(area) == 0){
      site <- columnSplit(plot, '_')[1]
      area <- taa[grep(site, taa$plot),'seedArea']
    }
  }
 
  area <- area[1]
  
  seed[,colnames(seedj)] <- seedj
  
  year <- rep(yr, each = length(tn))
  trap <- rep(tn, length(yr) )
  tr <- apply(cbind(plot, trap), 1, paste0, collapse='-')
  
  sd   <- data.frame( plot, trapID=tr, trap=trap, 
                      year=year,active=active, area = area, stringsAsFactors = F) 
  seed <- cbind(sd, seed)
  rownames(seed) <- 
    apply( cbind(plot, trap, year), 1, paste0, collapse='-')
  
  seed <- seed[seed$trapID %in% loc$trapID,]
  
  if('x' %in% colnames(loc) & 'UTMx' %in% colnames(loc)){
    loc$x <- loc$UTMx
    loc$y <- loc$UTMy
  }
  
  seed$trap <- as.character(seed$trap)
  loc$trap  <- as.character(loc$trap)
  
  list(counts = seed, xy = loc[,c('plot','trap','x','y')], active = active, 
       seedNames = seedNames )
}

multiStemBA <- function(diam){
  
  # effective diameter for BA matching diam for multiple stems
  # can be a single vector for one tree or a matrix, with each row being a different tree
  
  if(is.data.frame(diam)){
    diam <- as.matrix(diam)
  }
  if(is.matrix(diam)){
    return( round( sqrt(apply( diam^2, 1, sum, na.rm=T )), 1) )
  }
  
  round( sqrt( sum(diam^2, na.rm=T) ), 1 )
}

.treeFormat <- function( tfile, specNames = NULL, genusName = NULL, 
                         changeNames = NULL, plot,  
                         years = NULL, yrCols, MUSTHAVE = FALSE, keepCols = NULL, 
                         verbose = FALSE){
  
  # if MUSTHAVE, tree rejected if any yrCols are missing
  
  newplot <- plot
    
  midCD <- c( 273890.6, 3938623.3 )  #plot center for (x,y) at GSNP_CD
  
  site <- unlist(strsplit(plot,'_'))[1]
  
  CSV <- endsWith(tfile,'.csv')
  if(CSV){
    trees <- read.csv(tfile, stringsAsFactors = FALSE)
  }else{
    trees  <- read.table( tfile, header  = TRUE, stringsAsFactors = FALSE)
  }
  
  wd <- grep('diam', colnames(trees))
  
  kk <- which( sapply(trees[,wd], is.character) )
  
  if(length(kk) > 0){
    stop( paste('character in column',names(kk)) )
  }
  
  
  if('plotID' %in% colnames(trees))trees$plot <- trees$plotID
  if('individualID' %in% colnames(trees))trees$tree <- trees$individualID
  
  trees <- trees[which(trees$plot == plot),]
  if(nrow(trees) == 0)return( list(yrDat = numeric(0)) )
  
  if(is.null(specNames) & is.null(genusName)){
    specNames <- sort(unique(as.character(trees$species)))
  }
  if(is.null(specNames) & !is.null(genusName)){
    ss <- which( substr(trees$species, 1, 4) == substr(genusName, 1, 4) )
    specNames <- sort(unique(as.character(trees$species[ss])))
  }

  trees$ID   <- as.character(trees$tree)
  
  if('UTMX' %in% colnames(trees))colnames(trees)[colnames(trees) == 'UTMX'] <- 'UTMx'
  if('UTMY' %in% colnames(trees))colnames(trees)[colnames(trees) == 'UTMY'] <- 'UTMy'
  
  if(!'UTMx' %in% colnames(trees))trees$UTMx <- trees$x
  if(!'UTMy' %in% colnames(trees))trees$UTMy <- trees$y
  
  
 # if(!'x' %in% colnames(trees))trees$x <- trees$UTMx
 # if(!'y' %in% colnames(trees))trees$y <- trees$UTMy
  
  if(plot == "GSNP_CD"){
    trees$x <- trees$UTMx - midCD[1]
    trees$y <- trees$UTMy - midCD[2]
  }
  
  trees$x <- trees$UTMx
  trees$y <- trees$UTMy
  
  trees$plot <- newplot
  trees$species <- as.character(trees$species)
  
  if( is.null(genusName) ){
    wj     <- which(trees$species %in% specNames & is.finite(trees$x) &
                      is.finite(trees$y))
  }else{
    ww <- grep(genusName,as.character(trees$species) )
    wj <- ww[which(is.finite(trees$x[ww]) & is.finite(trees$y[ww]))]
    specNames <- sort(unique(trees$species[wj]))
  }
  if(length(wj) == 0)return( list(yrDat = numeric(0)) )
  
  trees  <- trees[drop = FALSE,wj,]
  
  mcol <- grep('diam', colnames(trees))
  dmm  <- suppressWarnings( max(trees[,mcol], na.rm=T) )
  if(dmm < 0)return( list(yrDat = numeric(0)) )
  
  wchange <- which( as.character(trees$species) %in% changeNames[,1])
  if(length(wchange) > 0){
    ts <- as.character(trees$species)
    wm <- match(ts[wchange],changeNames[,1])
    ts[wchange] <- changeNames[wm,2]
    trees$species[wchange] <- ts[wchange]
  }
  
  yrc <- yrCols
  ww  <- grep('crop',yrc)
  if(length(ww) > 0)yrc <- yrc[-ww]
  
  if(MUSTHAVE){
    
    for(k in 1:length(yrc)){                      
      
      tk <- as.matrix( trees[,grep(yrc[k],colnames(trees))] )
      wmin <- suppressWarnings( apply(tk,1,min, na.rm  = TRUE) )
      wna <- which( !is.finite( wmin ) )
      if(length(wna) > 0)trees <- trees[-wna,]
    }
    if(nrow(trees) == 0)return( list(yrDat = numeric(0)) )
  }
  
  # trees with multiple stems
  ID      <- as.character(trees$tree)
  dup     <- which(duplicated(ID))
  diamCol <- grep('diam',colnames(trees))
  
  if(length(dup) > 0){
    
    omit <- numeric(0)
    
    for(j in dup){
      
      idup <- which(ID == ID[j])
      tj   <- apply(trees[idup,diamCol,drop=F], 1, max, na.rm  = TRUE)
      tj   <- which.max(tj)
      omit <- c(omit, idup[-tj])
    }
    omit <- unique(omit)
    trees <- trees[-omit,]
  }
  
  ser <- grep( 'serotinous', colnames(trees) )
  if( length(ser) > 0 ){
    trees[,ser] <- lapply( trees[,ser], as.character )
    
    ww <- which(!trees$serotinous %in% c('B','NS','S',NA))
    if(length(ww) > 0 & verbose){
      print('serotinous includes:')
      print(trees$serotinous[ww])
    }
    serotinous <- rep(0, nrow(trees))
    serotinous[trees$serotinous == 'S'] <- 1
    serotinous[trees$serotinous == 'B'] <- .5
    trees$serotinous <- serotinous
  }
  
  dataYr <- as.numeric(columnSplit(colnames(trees)[diamCol],'diam')[,2] )
  
  wg <- grep('survey',colnames(trees))   #NEON plots
  if(length(wg) > 0){
    if(!'censinyr' %in% colnames(trees)){
      trees$censinyr <- dataYr[1]
    }
    if(!'deathyr' %in% colnames(trees)){
      trees$deathyr <- NA
    }
    if(!'censoryr' %in% colnames(trees)){
      trees$censoryr <- max(dataYr) + 1
    }
  }
  
  
  colnames(trees) <- .replaceString(colnames(trees), 'cropFractionSD', 'cropFractionSd')
  
  wk <- grep('cropCount',colnames(trees))
  if(length(wk) > 0)yrCols <- c(yrCols, 'cropCount')
  wk <- grep('cropFraction',colnames(trees))
  if(length(wk) > 0)yrCols <- c(yrCols, 'cropFraction')
  wk <- grep('cropFractionSd',colnames(trees))
  if(length(wk) > 0)yrCols <- c(yrCols, 'cropFractionSd')
  
  yrCols <- yrCols[!duplicated(yrCols)]
  
  yrange <- numeric(0)
  
  for(k in 1:length(yrCols)){
    
    wk <- which( startsWith(colnames(trees), yrCols[k]) )
    
    if(length(wk) == 0)next
    
    yk <- columnSplit(colnames(trees)[wk],yrCols[k])[,2]
    ww <- grep('Sd',yk)
    if(length(ww) > 0)next
    
    yk <- .replaceString(yk, 'Sd', '')
    yrange <- c(yrange, as.numeric(yk))
  }
  
  wk <- which( startsWith(colnames(trees), 'sex') )
  
  if(length(wk) > 0){
    yk <- columnSplit(colnames(trees)[wk],'sex')[,2]
    yrange <- c(yrange, as.numeric(yk))
  }
  years <- min(yrange):max(yrange)
  
  # reproduction
  
  snew <- matrix(NA, nrow(trees), length(years))
  colnames(snew) <- years
  
  scols <- as.matrix(trees[,grep('sex',colnames(trees)), drop = FALSE])
  
  if(ncol(scols) == 0)scols <- NULL
  
  if( !is.null(scols) ){
    
    ccol <- matrix( unlist( strsplit(colnames(scols),'sex') ), ncol=2,byrow=2)[,2]
    colnames(scols) <- ccol
    
    snot <- which(scols == 'N', arr.ind  = TRUE)
    srep <- which(scols == 'R' | scols == 'F', arr.ind  = TRUE)
    sfem <- which(scols == 'F', arr.ind  = TRUE)
    smal <- which(scols == 'M', arr.ind  = TRUE)
    
    smal <- smal[!smal[,1] %in% sfem[,1],]  # female overrides male
    if(!is.matrix(smal))smal <- matrix(smal,1)
    
    matr <- repr <- matrix(NA, nrow(trees), ncol(scols)) 
    colnames(matr) <- colnames(scols)
    
    colnames(matr) <- as.numeric(ccol)
    if(length(smal) > 0)repr[smal[,1],] <- 0
    if(length(sfem) > 0)repr[sfem[,1],] <- 1
    
    matr[snot] <- 0
    matr[srep] <- 1
    
    snew[,colnames(matr)] <- matr
  }
  
  
  yrDat  <- as.data.frame(matrix(NA, length(years)*nrow(trees), 
                                 length(yrCols)),
                          stringsAsFactors = F)
  id     <- sort( unique(as.character(trees$ID)) )
  tindex <- rep(id, each=length(years) )
  yindex <- rep(years, length(id) )
  
  index <- apply(cbind(tindex,yindex),1,paste0,collapse='-')
  rownames(yrDat) <- index
  colnames(yrDat) <- yrCols
  
  
  if( !'censinyr' %in% colnames(trees) ){
    trees$censinyr <- min(years,na.rm  = TRUE)
  }
  if( !'censoryr' %in% colnames(trees) ){
    trees$censoryr <- NA
  }
  if( !'deathyr' %in% colnames(trees) ){
    trees$deathyr <- NA
  }
  
  
  elev <- NA
  if('elev' %in% colnames(trees))elev <- trees[,'elev']
  
  xytree <- trees[,c('UTMx','UTMy')]
  colnames(xytree) <- c('x','y')
  id     <- apply( cbind(newplot, as.character(trees$ID)), 1, paste0, 
                   collapse='-')
  xytree <- data.frame(treeID = id, plot = newplot, tree = trees$ID, xytree)
  xytree$elev <- elev
  
  yrCols <- yrCols[yrCols != 'sex']

  fcols <- grep('crop',yrCols)
  
  for(k in 1:length(yrCols)){                      # yrCol2017
    
    tk <- as.matrix( trees[,grep(yrCols[k],colnames(trees)),drop = FALSE] )
    rownames(tk) <- trees$tree
    if(ncol(tk) == 0){
      yrDat[,k] <- NA
      next
    }
    
    if( yrCols[k] == 'cropFraction' & length( grep('Sd', colnames(tk)) ) > 0 ){
      sdcol <- grep('Sd', colnames(tk))
      tk    <- tk[,-sdcol,drop = FALSE]
    }
    
    dmat <- matrix(NA, nrow(tk), length(years))
    rownames(dmat) <- rownames(tk)
    colnames(dmat) <- years
    
    kk <- matrix( unlist(strsplit(colnames(tk), yrCols[k])), ncol=2, byrow  = TRUE)[,2]
                       
    yk <- as.numeric(kk)
    if(min(yk) < min(years) | max(yk) > max(years))
      stop('\nyears do not span diam years\n')
    yseq <- min(yk):max(yk)
    
      
    syr  <- trees[,'censinyr']                             # left censored
    imat <- matrix(years, nrow(trees), length(years), byrow  = TRUE) 
    wmin <- which(is.finite(tk),arr.ind  = TRUE)                 # measurements exist
   
    tmp <- tapply(yk[wmin[,2]], list(tree=wmin[,1], yr=wmin[,1]*0+1), FUN=min)
    t2  <- matrix(NA, nrow(trees), ncol(tmp))
    t2[as.numeric(rownames(tmp)),] <- tmp
    tmp <- t2
    
    syr  <- pmin(syr, tmp, na.rm  = TRUE)
    
    syr  <- matrix(syr, nrow(trees), length(years))           # start
    eyr  <- suppressWarnings( apply(trees[,c('deathyr','censoryr')],
                                    1, min, na.rm  = TRUE) )
    eyr  <- matrix(eyr, nrow(trees), length(years))               # end
    eyr[eyr == Inf] <- max(years) 
    imat[imat < syr | imat > eyr] <- NA
    tmat <- matrix(trees[,'ID'], nrow(tk), length(years))
    
    icol <- match(yk, years)
    irow <- rep(1:nrow(tk),each=length(icol))
    icol <- rep(icol, nrow(tk))
    jcol <- rep(1:length(yk),nrow(tk))
    dmat[ cbind(irow,icol) ] <- tk[ cbind(irow, jcol) ]   # interpolate here
    
    start <- match( suppressWarnings( apply(imat,1,min,na.rm  = TRUE)), years)
    end   <- match( suppressWarnings( apply(imat,1,max,na.rm  = TRUE)), years)
    
    wf <- which(is.finite(start) & is.finite(end) & end > start)
    
    if( !k %in% fcols ){
      
      if(length(wf) > 0){
        minVal <- 0
        maxVal <- 400
        rr     <- 1
        tinySlope <- .0001
        if(yrCols[k] == 'canopy'){
          minVal <- suppressWarnings( apply( dmat[drop = FALSE,wf,], 1, min, na.rm = T) )
          maxVal <- suppressWarnings( apply( dmat[drop = FALSE,wf,], 1, max, na.rm = T) )
          minVal[ minVal < 1 | !is.finite(minVal) ] <- 1
          maxVal[ maxVal > 5 | !is.finite(maxVal) ] <- 5
          rr     <- 0
          tinySlope <- NULL
        } 
        
        tmp <- .interpRows(dmat[drop = FALSE,wf,],startIndex=start[wf],endIndex=end[wf],
                           INCREASING = FALSE,minVal=minVal,maxVal=maxVal,
                           defaultValue=NULL,tinySlope = tinySlope)
        dmat[wf,] <- round(tmp, rr)
      }
      
      if( yrCols[k] == 'canopy' ){
        maxVal <- suppressWarnings( apply(dmat[drop = FALSE,wf,], 1, max,na.rm  = TRUE) )
        wmm    <- which(is.finite(maxVal))
        mmat   <- matrix(maxVal[wmm],length(wmm),ncol(tmp))
        ww     <- which(tmp[wmm,] > mmat)
        tmp[wmm,][ww] <- mmat[ww]
        dmat[wf,] <- round(tmp)
      }
    }
    
    dvec <- dmat[is.finite(imat), drop=F]
    
    it   <- tmat[is.finite(imat)]
    iy   <- imat[is.finite(imat)]
    inn  <- apply(cbind(it,iy),1,paste0,collapse='-')
    mm   <- match(inn, index)
    yrDat[mm,k] <- dvec
    
  }
  
  mm <- match(tindex,trees[,'tree'])

  spec <- trees[mm,'species']
  repr <- rep(NA, nrow(yrDat))
  repr[match(inn, index)] <- snew[is.finite(imat)]
  if('cropCount' %in% colnames(yrDat)){
    repr[yrDat$cropCount > 0] <- 1
  }
  
  treeID <- columnPaste(newplot,tindex)
  yrDat <- data.frame(plot = newplot, treeID = treeID, tree = tindex, species = spec, 
                      year = yindex, 
                      repr = repr, yrDat, stringsAsFactors = F)
  if(length(keepCols) > 0)
      yrDat <- cbind(yrDat,trees[,keepCols])
  
  if('taxonID' %in% colnames(trees))yrDat$taxonID <- trees[mm,'taxonID']
  if(!'survey' %in% colnames(trees))trees$survey <- 'DUKE'
  if('survey' %in% colnames(trees)){
    yrDat$survey <- trees[mm,'survey']
    yrDat$x <- trees[mm,'x']
    yrDat$y <- trees[mm,'y']
    yrDat$UTMx <- trees[mm,'UTMx']
    yrDat$UTMy <- trees[mm,'UTMy']
  }
  if('serotinous' %in% colnames(trees))yrDat$serotinous[mm] <- trees$serotinous[mm]
  
  if(MUSTHAVE)yrDat <- yrDat[is.finite(yrDat$diam) & yrDat$diam > 1,]
  xytree <- xytree[xytree$treeID %in% yrDat$treeID,]
  
  rownames(yrDat) <- NULL
  
 # yrDat <- yrDat[drop = FALSE, is.finite(yrDat$diam),]
  xytree <- xytree[xytree$treeID %in% yrDat$treeID,]
  
  rownames(yrDat) <- columnPaste(yrDat$treeID, yrDat$year, '_')
  
  if('censoryr' %in% colnames(trees)){
    
    wc <- which(is.finite(trees$censoryr))
    tid <- columnPaste(trees$plot[wc], trees$tree[wc])
    tid <- columnPaste(tid, trees$censoryr[wc], '_')
    yrDat <- yrDat[!rownames(yrDat) %in% tid,]
    xytree <- xytree[xytree$treeID %in% yrDat$treeID,]
  }
  
  list(yrDat = yrDat, xytree = xytree[,c('plot','tree','treeID','x','y','elev')])
}

.fac2num <- function(xx){ 
  
  dims <- dn <- NULL
  
  if(!is.null(ncol(xx))){
    dims <- dim(xx)
    dn   <- dimnames(xx)
  }
  xx <- if(is.list(xx))unlist(xx)
  xx <- as.numeric(as.character(xx)) 
  if(!is.null(dims))xx <- matrix(xx, dims[1], dims[2], 
                                 dimnames = dn)
  xx
}

.replaceString <- function(xx, now='_', new=' '){  #replace now string in vector with new
  
  if(!is.character(xx))xx <- as.character(xx)
  
  ww <- grep(now,xx,fixed  = TRUE)
  if(length(ww) == 0)return(xx)
  
  if(length(new) == 1){
    for(k in ww){
      s  <- unlist( strsplit(xx[k],now,fixed  = TRUE) )
      ss <- s[1]
      if(length(s) == 1)ss <- paste( ss,new,sep='')
      if(length(s) > 1)for(kk in 2:length(s)) ss <- paste( ss,s[kk],sep=new)
      xx[k] <- ss
    }
  }else{
    # new is a vector
    s  <- unlist( strsplit(xx,now,fixed  = TRUE) )
    nx <- length(xx)
    nc <- length(s)/length(xx)

    ss <- matrix(s, ncol=nc, byrow  = TRUE)
    nl <- nchar(ss[1,])
    
    if(nl[1] == 0)ss <- paste(new, ss[,2], sep='')
    if(nl[2] == 0)ss <- paste(ss[,1], new, sep='')
      
    xx <- ss
  }
  xx
}

.Iformat2Var <- function(iname){
  
  tt <- .replaceString(iname, 'I(','')
  tt <- .replaceString(tt, 'log(','')
  tt <- .replaceString(tt, 'sqrt(','')
  tt <- .replaceString(tt, '^2','')
  tt <- .replaceString(tt, ')','')
  tt <- .replaceString(tt, ' ','')
  tt
}

.get.model.frame <- function(formula, data){
  
  tmp <- model.frame(formula, data, na.action=NULL)
  
  wchar <- which( sapply(tmp, is.character) )
  
  if(length(wchar) == 0)return( tmp )
  
  for(k in 1:length(wchar)){
    data[,wchar[k]] <- as.character( data[,wchar[k]] )
  }
  model.frame(formula, data, na.action=NULL)
}
  
.getDesign <- function(formula, data, verbose = FALSE){
  
  # one set of columns for each tree species, retain NAs
  
  specNames <- sort(unique(as.character(data$species)))
  nspec     <- length(specNames)
  
  f1 <- paste0( as.character(formula), collapse='')
  
  if(f1  == '~1' & nspec == 1){
    x  <-  matrix(1, nrow(data), 1)
    colnames(x) <- "(Intercept)"
    return( list(x = x, missx = integer(0), specCols = numeric(0) ) )
  }
  
  data$species <- factor(data$species)
  
  attr(data$species,'contrasts') <- contrasts(data$species, contrasts = FALSE)
  
  tmp1 <- .get.model.frame(formula, data )
  tn1  <- attr( terms(tmp1), 'dataClasses' )
  sp1  <- names(tn1)[tn1 == 'numeric' | tn1 == 'nmatrix.1']
  sp1  <- .Iformat2Var(sp1)
  sp1  <- unique(sp1)
  miss <- which(is.na(data[,sp1,drop=F]),arr.ind  = TRUE)
  
  if(length(miss) > 0){
    xmean <- colMeans(data[,sp1,drop=F], na.rm  = TRUE)
    data[,sp1][miss] <- 1e+10
    if(verbose)cat('\nNote: missing values in xfec filled with mean values\n')
  }
  
  x  <- model.matrix(formula, data)
  if(nspec > 1){
    ws <- grep('species',colnames(x))
    x  <- x[,ws]
  }
  missx <- which(x > 1e+9,arr.ind  = TRUE)
  
  if(length(missx) > 0){
    data[,sp1][miss] <- xmean[miss[,2]]
    
    mux <- apply(x, 2, mean, na.rm  = TRUE)
    x[missx] <- mux[missx[,2]]
  }
  x  <- model.matrix(formula, data)
  if(nspec > 1){
    ws <- grep('species',colnames(x))
    x  <- x[,ws]
  }
  
  specCols  <- numeric(0)
  if(nspec > 1){
    for(j in 1:length(specNames)){
      specCols <- rbind(specCols, grep( paste('species',specNames[j],sep=''),
                                        colnames(x)))
    }
    rownames(specCols) <- specNames
  }
  
  list(x = x, missx = missx, specCols = specCols)
}

columnSplit <- function(vec, sep='_', ASFACTOR = F, ASNUMERIC = FALSE,
                        LASTONLY = FALSE){
  
  vec <- as.character(vec)
  nc  <- length( strsplit(vec[1], sep, fixed  = TRUE)[[1]] )
  
  mat <- matrix( unlist( strsplit(vec, sep, fixed  = TRUE) ), ncol=nc, byrow  = TRUE )
  if(LASTONLY & ncol(mat) > 2){
    rnn <- mat[,1]
    for(k in 2:(ncol(mat)-1)){
      rnn <- columnPaste(rnn,mat[,k])
    }
    mat <- cbind(rnn,mat[,ncol(mat)])
  }
  if(ASNUMERIC){
    mat <- matrix( as.numeric(mat), ncol=nc )
  }
  if(ASFACTOR){
    mat <- data.frame(mat)
  }
  if(LASTONLY)mat <- mat[,2]
  mat
}

columnPaste <- function(c1, c2, sep='-', NOSPACE = FALSE){
  
  c1    <- as.character(c1)
  c2    <- as.character(c2)
  if(NOSPACE){
    c1   <- .replaceString(c1, ' ', '')
    c2   <- .replaceString(c2, ' ', '')
  }
  c12   <- apply( cbind(c1, c2) , 1, paste0, collapse=sep)

  c12
}

specPriorVector <- function(pVar, tdata){
  
  if(length(pVar) > 1){
    pVar <- pVar[tdata$species]
  }else{
    pVar <- rep(pVar, nrow(tdata))
  }
  pVar
}



setupDistMat <- function(tdata1, sdata1, xytree1, xytrap1, verbose){
  
  # set up distall
  # sdata$drow are rows in distall
  # 
  
  pord <- sort(tdata1$plot)
  pord <- pord[!duplicated(pord)]
  
  treeIDs <- unique(as.character(tdata1$treeID))  # do not sort

  plotRm <- plotKp <- character(0)
  distTreeID <- numeric(0)
  tdata1$dcol <- sdata1$drow <- NA
  
  for(j in 1:length(pord)){
    
    tj <- which(xytree1$plot == pord[j] & xytree1$fit == 1)
    
    if(length(tj) == 0){
      plotRm <- c(plotRm, pord[j])
      next
    }
  }
  
  if(length(plotRm) > 0 & verbose){
    
    kk <- which(!sdata1$plot %in% plotRm)
    sdata1 <- sdata1[kk,]
    kk <- which(!xytrap1$plot %in% plotRm)
    xytrap1 <- xytrap1[kk,]
    
    kp <- paste0( plotRm, collapse=', ')
    kwords <- paste(' Trees absent from xytree in',kp)
    cat( paste('\nNote: ',kwords, '\n') )
  }
  
  trapIDs <- unique(as.character(sdata1$trapID))
  pord <- sort(sdata1$plot)
  pord <- pord[!duplicated(pord)]
  
  for(j in 1:length(pord)){
    
    tj <- which(xytree1$plot == pord[j] & xytree1$fit == 1)
    sj <- which(xytrap1$plot == pord[j])
    
    if(length(tj) == 0){
      plotRm <- c(plotRm, pord[j])
      next
    }
    if(length(sj) == 0)stop( paste('plot', pord[j] ,'has no traps in xytrap') )
    
    xy1     <- xytree1[tj,]
    xy2     <- xytrap1[sj,]
    
    species <- as.character(xytree1$species[tj])
    da      <- .distmat(xy1[,'x'],xy1[,'y'],xy2[,'x'],xy2[,'y']) 
    da      <- round(da, 1)
    sa      <- matrix(species, nrow(da), ncol(da), byrow=T)
    rownames(da) <- xy2$trapID
    
    distTreeID <- append( distTreeID, list(xytree1$treeID[tj]))
    plotKp     <- c(plotKp, pord[j])
    
    jj <- match(tdata1$treeID, xytree1$treeID[tj])
    wf <- which(is.finite(jj))
    tdata1$dcol[wf] <- jj[wf]
    
    if(j == 1){
      distall <- da
      specall <- sa
      
    }else{
      if( ncol(da) > ncol(distall) ){
        nnc     <- ncol(da) - ncol(distall)
        newCols <- matrix(NA, nrow(distall), nnc )
        distall <- cbind(distall,newCols)
        newCols <- matrix(species[1], nrow(distall), nnc )
        specall <- cbind(specall, newCols)
      }
      if( ncol(distall) > ncol(da) ){
        nnc     <- ncol(distall) - ncol(da)
        newCols <- matrix(NA, nrow(da), nnc)
        da <- cbind(da,newCols)
        newCols <- matrix(species[1], nrow(da), nnc)
        sa <- cbind(sa,newCols)
      }
      distall <- rbind(distall, da)
      specall <- rbind(specall, sa)
    }
  }
  rownames(specall) <- rownames(distall)
  names(distTreeID) <- plotKp
  
  sdata1$drow <- match(sdata1$trapID, rownames(distall))
  
  trapid  <- rownames(distall)
  attr(distall,'species') <- specall
  
  list(tdata = tdata1, sdata = sdata1, distall = distall, distTreeID = distTreeID)
}


kernYrR <- function(dmat, fec, seedrow, treecol, plotyrs,
                    treeplotYr, seedplotYr){
  
  # kernYrRcpp in R
  
  ny <- length(plotyrs)
  nf <- ncol(fec)
  lambda <- matrix(0, nr, nf)
  
  for(j in 1:ny){
    
    ws <- which(seedplotYr == plotyrs[j])
    wt <- which(treeplotYr == plotyrs[j])
    ds <- seedrow[ ws ]
    dt <- treecol[ wt ]
    
    lambda[ws,] <- dmat[ds,dt]%*%fec[wt,]
  }
  lambda
}

trimData <- function(treeData, xytree, seedData = NULL, xytrap = NULL,
                     formulaFec = NULL, formulaRep = NULL, 
                     specNames = NULL, seedNames = NULL){
  
  # when treeData is trimmed, clean other data.frames
  
  if(is.null(specNames))specNames <- sort(unique(treeData$species))
  
  nspec <- length(specNames)
  
  xid      <- columnPaste(xytree$plot,xytree$tree)
  tid      <- columnPaste(treeData$plot,treeData$tree)
  xytree   <- xytree[xid %in% tid, ]
  
  if(!is.null(seedData)){
    seedData <- seedData[seedData$plot %in% treeData$plot, ] 
    xytrap   <- xytrap[xytrap$plot %in% treeData$plot, ]
    if(nrow(seedData) == 0)seedData <- xytrap <- NULL
    
    if(nrow(xytree) == 0){
      xytree <- NULL
    }else{
      mm <- match(tid, xid)        # add locations to treeData
      treeData$x <- xytree$x[mm]
      treeData$y <- xytree$y[mm]
      xid      <- columnPaste(xytree$plot,xytree$tree)
      tid      <- columnPaste(treeData$plot,treeData$tree)
      xytree   <- xytree[xid %in% tid, ]
    }
  }
  if(!is.null(seedNames))seedNames <- seedNames[ seedNames %in% colnames(seedData) ]
  nold  <- nspec
  nspec <- length(specNames)
  
  if(nspec == 1){
    ff <- as.character(formulaFec)
  
    if( length( grep('species', ff[2]) ) > 0 ){
      ff <- .replaceString( ff, 'species * ', '')
      formulaFec <- as.formula( paste0( ff[1], ff[2], collapse='') )
      ff <- as.character(formulaRep)
      ff <- .replaceString( ff, 'species * ', '')
      formulaRep <- as.formula( paste0( ff[1], ff[2], collapse='') )
    }
  }
  specNames <- sort(unique(treeData$species))
  
  list(treeData = treeData, xytree = xytree, seedData = seedData, xytrap = xytrap,
       formulaFec = formulaFec, formulaRep = formulaRep, 
       specNames = specNames, seedNames = seedNames)
}


checkPlotName <- function(pdata){
  
  nt <- stringr::str_count(pdata[,'plot'], "_")
  w0 <- which(nt != 1)
  if(length(w0) > 0){
    print("plot names without one and only one '_'")
    stop( pdata[w0,'plot'] )
  }
}


treeSeedPlots <- function( tdata, sdata, xytree, xytrap, 
                           seedNames = NULL, minTrees = 5,
                           minYears = 2,
                           formulaFec = NULL, formulaRep = NULL, 
                           skipSpec = character(0), omitNames = character(0), 
                           combineSeeds = NULL, combineSpecs = NULL, 
                           cropCountCols = c('fecMin','cropCount'), 
                           verbose = FALSE){
  SEEDDATA <- FALSE
  
  if( !is.null(sdata) ){
    if(nrow(sdata) > 1)SEEDDATA <- T
  }
  
 # checkPlotName(tdata)
 # checkPlotName(xytree)
  
  # sufficient years
  pvec  <- tdata$plot
  yvec  <- tdata$year
  if(SEEDDATA){
    pvec <- c(pvec, sdata$plot)
    yvec <- c(yvec, sdata$year)
  }
  
  tplus <- names( which( sapply( tapply( yvec, pvec, range ), diff ) >= (minYears-1) ) )
  ww    <- which( tdata$plot %in% tplus )
  
  tdata <- tdata[ tdata$species %in% tdata$species[ww], ]
  
  if(nrow(tdata) == 0){
    print( 'no trees with enough years' )
    return( list( treeData = NULL ) )
  }
  
  tid   <- columnPaste(tdata$plot, tdata$tree)
  maxd  <- tapply(tdata$diam, tid, max, na.rm=T)
  maxd  <- maxd[ tid ]
  maxy  <- max(table(tid))
  
  mm    <- match(tdata$species, priorTable[,code])  
  mind  <- priorTable$minDiam[ mm ]
  sites <- columnSplit(tdata$plot, '_')[,1]
  
  wi    <- which(!duplicated(tid) & maxd > mind/2 )
  sbys  <- table(tdata$species[wi], sites[wi])
  stot  <- rowSums(sbys)                           # trees > mind
  
  if(verbose){
    print( 'number trees above min diam:' )
    print(stot)
  }
  
  specNames <- names(stot)[ stot > minTrees ]      # minimum large trees
  if(length(specNames) == 0){
    print( 'no trees big enough' )
    return( list( treeData = NULL ) )
  }
  
  tdata <- tdata[ tdata$species %in% specNames, ]
  tid   <- columnPaste(tdata$plot, tdata$tree)
  maxd  <- tapply(tdata$diam, tid, max, na.rm=T)
  maxd  <- maxd[ tid ]
  
  mm    <- match(tdata$species, priorTable[,code])  
  mind  <- priorTable$minDiam[ mm ]
  sites <- columnSplit(tdata$plot, '_')[,1]
  
  
  if( SEEDDATA ){
    
    if( is.null(seedNames) ){
      ss <- colnames(sdata)[startsWith(colnames(sdata),genCode[t])]
      ss <- .replaceString(ss,'_min','')
      ss <- .replaceString(ss,'_max','')
      seedNames <- sort(unique(ss))
    }
    sdata <- sdata[, !colnames(sdata) %in% skipSpec]
    gen4 <- substr(tdata$species[1], 1, gee)
    ws   <- which( startsWith(colnames(sdata), gen4)  )
    cc   <- colnames(sdata)[ws]
    wk   <- grep('_min', colnames(sdata))
    wm   <- grep('_max', colnames(sdata))
    ws   <- ws[ !ws %in% c(wk, wm) ]
    
    seedNames <- colnames(sdata)[ws]
    totSeed   <- sum(as.vector(sdata[,cc]), na.rm=T)
  #  checkPlotName(sdata)
  #  checkPlotName(xytrap)
    
    onames <- c( omitNames, paste(seedNames, 'cones', sep='_'),
                 paste(seedNames, 'emptyseeds', sep='_'))
    
    womit <- which(seedNames %in% onames)
    if(length(womit) > 0)seedNames <- seedNames[-womit]
    
    specNames <- sort(unique(as.character(tdata$species)))
    
    wc <- grep('_caps', colnames(sdata))
    if(length(wc) > 0){
      from <- colnames(sdata)[wc]
      to   <- columnSplit(from, '_')[,1]
      wu    <- which(nchar(to) == 4)
      if(length(wu) > 0){
        if(length(wu) > 1){
          print( 'multiple short seedNames' )
          return( list(treeData = NULL) )
        }
        to[wu] <- paste(to[wu], 'UNKN',sep='')
      }
      cmore <- cbind(from,to)
      combineSeeds <- rbind(combineSeeds, cmore)
    }
    
    tmp <- combineSeedNames(sdata, seedNames, 
                            rbind(combineSeeds,combineSpecs))
    sdata     <- tmp$seedData
    seedNames <- tmp$seedNames
    
    wg <- unlist( gregexpr("[A-Z]", specNames) ) # position where species starts
    
    als <- c(specNames, paste(specNames, '_min', sep=''),
             paste(specNames,'_max',sep=''),
             paste( substr(specNames, 1, wg-1), 'UNKN',sep='') )
    ttt <- seedNames[seedNames %in% als]
    
    if(length(ttt) == 0){
      print(specNames)
      print(seedNames)
      print('seedNames not in specNames')
      return( list(treeData = NULL) )
    }
    seedNames <- ttt
  }else{
    sdata <- NULL
  }
  
  tmp <- combineSpecies(tdata$species, specNames, combineSpecs)
  tdata$species <- tmp$species
  specNames     <- tmp$specNames
  nspec <- length(specNames)
  
  # remove a species if no cropCount, no seeds on any plot, less than 2 yr
  CROP <- FALSE
  
  wcc <- which( cropCountCols %in% colnames(tdata) )
  if( length(wcc) > 0 ){
    treeSeeds <- rowSums(tdata[,cropCountCols[wcc],drop=F], na.rm=T)
    CROP <- TRUE
  }

  allp <- sort(unique( c(tdata$plot, sdata$plot) ))
  allt <- unique( c(specNames, seedNames) )
  tpy  <- matrix( 0, length(allp), length(allt) )   #tree years by plot/species
  colnames(tpy) <- allt
  rownames(tpy) <- allp
  spy <- tpy                                        #seeds counted
  
  for(j in 1:length(allp)){
    
    ssum <- tpy[1,]*0
    
    wt <- which(tdata$plot == allp[j] &  (mind < maxd | treeSeeds > 0) )
    ws <- which(sdata$plot == allp[j])
    ttab <- table(tdata$species[wt], tdata$year[wt])
    ttab[ ttab > 0 ] <- 1
    tyr  <- rowSums(ttab)
    if(CROP){
      crop <- tdata[wt, cropCountCols[wcc] ]
      if( length(wcc) > 1 )crop <- rowSums(crop, na.rm=T)
      crop <- tapply( crop, tdata$species[wt], sum, na.rm=T)
      crop[ crop > 1 ] <- 1
      ssum[ names(crop) ] <- crop
    }

    if(length(ws) > 0){
      
      stab <- suppressWarnings( 
        tapply( as.matrix(sdata[ws,seedNames, drop=FALSE]), 
                list( year = rep(sdata$year[ws], length(seedNames)), 
                      type = rep(seedNames, each = length(ws))), max, na.rm=T) )
      stab[ stab > 1 ] <- 1
      ssum <- colSums(stab)  # any seeds?
      if(CROP)ssum[ names(crop) ] <- ssum[ names(crop) ] + crop
      
      syr  <- ssum*0 + nrow(stab)
      syr  <- syr[ names(syr) %in% names(tyr) ]
      tyr[ names(syr) ] <- apply( rbind(tyr[ names(syr) ],syr), 2, max )
    }
    spy[j, names(ssum)] <- ssum
    tpy[ j, names(tyr) ] <- tyr
  }
  # plots/species with no tree years | (tree years but never any seeds)
  
  w0 <- which( colSums(spy) == 0 & colSums(tpy) == 0 )  # never any seeds
  wu <- which( endsWith( colnames(spy), "UNKN" ) )      # could be seeds in UNKN class
  w0 <- w0[ !w0 %in% wu ]
  
  if( length(w0) > 0 ){
    tpy[, -w0, drop=F] <- 0    # never seeds
  }

  keep <- outer( rownames(tpy), colnames(tpy), paste, sep='-' )
  pltr <- columnPaste( tdata$plot, tdata$species )
  keep <- which( pltr %in% keep )
  
  tdata <- tdata[keep,]
  
  if(nrow(tdata) == 0){
    print( 'no trees produced seeds' )
    return( list( treeData = NULL ) )
  }
  
  ttt   <- trimData(tdata, xytree, sdata, xytrap,
                    formulaFec, formulaRep, specNames, seedNames)
  tdata <- ttt$treeData
  sdata <- ttt$seedData
  xytree <- ttt$xytree
  xytrap <- ttt$xytrap
  specNames <- ttt$specNames
  seedNames <- ttt$seedNames
  formulaFec <- ttt$formulaFec
  formulaRep <- ttt$formulaRep
  
  # species by site
  site <- columnSplit(tdata$plot,'_')[,1]
  wg    <- grep('.', site, fixed=T)
  if(length(wg) > 0)site[wg] <- columnSplit(site[wg],'.')[,1]
  if(verbose)  print( table(site,tdata$species) )
  
  
  list( treeData = tdata, seedData = sdata, xytree = xytree, xytrap = xytrap, 
        specNames = specNames, seedNames = seedNames, formulaFec = formulaFec,
        formulaRep = formulaRep)
}


.setupData <- function(formulaFec, formulaRep, mainEffects, tdata, sdata,
                       xytree, xytrap, specNames, seedNames, AR, YR, 
                       yearEffect, minDiam, maxDiam, TREESONLY, maxFec, CONES,
                       notFit, seedTraits = NULL, verbose = FALSE){
  
  SEEDDATA <- TRUE
  if(is.null(sdata))SEEDDATA <- FALSE
  
  # formulas have 'species *' already
  arList <- numeric(0)
  plag <- p <- 0
  
  if(!is.null(seedTraits)){
    ww <- which(!specNames %in% rownames(seedTraits))
    if(length(ww) > 0)
      stop( paste('\nspecNames not in rownames(seedTraits): ',specNames[ww],
                  sep='') )
  }
  
  if(!is.null(yearEffect)){
    if('p' %in% names(yearEffect))p <- plag <- yearEffect$p
  }
  
  notCols <- designTable <- NULL
  words <- character(0)
  
  plotNames <- sort(unique( as.character(tdata$plot) ))
  years <- sort(unique( tdata$year[tdata$obs == 1]) )
  if(SEEDDATA)years <- sort(unique(c( sdata$year[sdata$obs == 1], years )))
  years <- c(min(years):max(years))
  
  nplot <- length(plotNames)
  nspec <- length(specNames)
  
  # if no crop count and seed never observed, remove
  
  ttt <- treeSeedPlots( tdata, sdata, xytree, xytrap, 
                        seedNames = seedNames, 
                        formulaFec = formulaFec, formulaRep = formulaRep )
  tdata  <- ttt$treeData
  
  if(length(tdata) == 0)stop( 'insufficient data' )
  
  sdata  <- ttt$seedData
  xytree <- ttt$xytree
  xytrap <- ttt$xytrap
  specNames  <- ttt$specNames
  seedNames  <- ttt$seedNames
  formulaFec <- ttt$formulaFec
  formulaRep <- ttt$formulaRep
  
  nspec <- length(specNames)
  nseed <- length(seedNames)
  
    
  # note: reorganizes tdata to put known immature and serotinous at end
  # trees in plots with seeds, but only small trees: decrease minDiam
  
  minD  <- specPriorVector(minDiam, tdata)
  maxD  <- specPriorVector(maxDiam, tdata)
  wd    <- which(tdata$diam > minD)
  
  ttab <- table(tdata$plot, tdata$species)
  dtab <- ttab*0
  tmp  <- table(tdata$plot[wd],tdata$species[wd])
  dtab[rownames(tmp),colnames(tmp)] <- tmp
  ww <- which(ttab > 0 & dtab == 0, arr.ind  = TRUE)
  
  if( SEEDDATA & length(ww) > 0 ){  # if there are seeds and nothing mature
                                    # then largest tree on plot initialized as mature
    for(k in 1:nrow(ww)){
      
      wk <- which( tdata$plot == rownames(ww)[k] &
                     tdata$species == colnames(ttab)[ww[k,2]] )
      qk <- quantile(tdata$diam[wk], .75)
      
      if( min(minDiam) > qk )minD[wk] <- qk
      
      tdata$repr[wk][tdata$diam[wk] >= qk] <- 1
      if( 'lastRepr' %in% colnames(tdata) ){
        tdata$lastRepr[wk][tdata$diam[wk] >= qk] <- 1
        tdata$lastFec[wk][tdata$diam[wk] >= qk] <- 1.1
      }
      if( 'serotinous' %in% colnames(tdata) ){
        tdata$serotinous[wk][tdata$diam[wk] > qk] <- 0
      }
    }
  }
  
  if(verbose) print(seedTraits)
  
  tmp <- setupZ(tdata, xytree, specNames, years, minD, maxD, maxFec, CONES,  
                seedTraits, verbose)
  z           <- tmp$z
  zmat        <- tmp$zmat
  zknown      <- tmp$zknown
  matYr       <- tmp$matYr 
  last0first1 <- tmp$last0first1
  tdata       <- tmp$tdata       # fecMin, fecMax included
  fecMinCurrent <- tmp$fecMinCurrent
  fecMaxCurrent <- tmp$fecMaxCurrent
  fstart        <- tmp$fstart
  seedTraits    <- tmp$seedTraits
  tdata         <- tmp$tdata
  
  if(verbose){
    cat('\nMaximum diameter:\n')
    ww <- which.max(tdata$diam)
    print(tdata[ww,])
  }
  
  xytree$fit <- last0first1[xytree$treeID, 'fit']
  
  if(SEEDDATA){
    tdata$obsTrap <- addObsTrap(tdata, sdata)
  }else{
    tdata$obsTrap <- 0
  }
  
  mm <- last0first1[tdata$treeID,'fit']
  tdata$obsTrap <- tdata$obsTrap*mm
  tdata$fit      <- mm
  
  # no more reordering
  
  if(is.character(tdata$species))tdata$species <- as.factor(tdata$species)
  
  yeGr <- specNames[1]
  
  tdata$group <- 1
  
  if( YR | AR ){
    
    gnames <- character(0)
    if('groups' %in% names(yearEffect)){
      gnames <- yearEffect$groups
    }else{
      tdata$group <- 1
      gnames <- 'group'
    }
    
    wg <- which(!gnames %in% colnames(tdata))
    if(length(wg) > 0){
      words <- paste('columns specified in yearEffect not found in treeData:\n"',
                     gnames[wg],'"',sep='' )
      stop( words )
    }
    
    group <- sapply(tdata[,gnames],as.character)
    
    if(length(gnames) > 1){
      ws <- which(gnames == 'species')     # 'species' must be last
      if(length(ws) > 0){
        gnames <- c('species',gnames)
        gnames <- gnames[!duplicated(gnames)]
        gnames <- rev(gnames)
      }
      group <- apply( tdata[,gnames], 1, paste0, collapse='-')
    }
    group <- .replaceString(group, ' ', '')
    tdata$groupName <- group
    yeGr <- sort(unique(as.character(group)))
    tdata$group <- match(as.character(group), yeGr)
  }
  
  if( !'group' %in% colnames(tdata) )tdata$group <- 1
  
  allYears <- min(tdata$year):max(tdata$year)
  
  if( AR ){        
    
    plag <- yearEffect$p
    tmp  <- msarSetup(tdata, plag, icol='treeID', jcol = 'year',
                       gcol = 'groupName', yeGr, verbose = verbose)
    groupByInd <- tmp$groupByInd
    betaYr     <- tmp$betaYr
    ngroup     <- length(yeGr)
    
    ttt <- tmp$xdata
    rownames(ttt) <- columnPaste(ttt$treeID, ttt$year, '_')
    ttt <- ttt[rownames(ttt) %in% rownames(tdata),]
    ttt$obs <- tdata$obs
    ttt$obsTrap <- tdata$obsTrap
    
    if('cropCount' %in% colnames(tdata))ttt$cropCount <- tdata$cropCount
    
    tdata <- ttt
    
    rm(ttt)
    
    tdata$plotYr <- columnPaste(tdata$plot,tdata$year,'_')
    allYears <- sort(unique(tdata$year))
    tdata$times  <- match(tdata$year, allYears)
    plotYrs      <- sort(unique(tdata$plotYr))
    tdata$plotyr <- match(tdata$plotYr, plotYrs)
    tdata$group <- match(as.character(group), yeGr)
    
    lagMat <- msarLagTemplate(plag, tdata, icol='treeID', jcol='year', 
                              gcol='group', ocol='obs', yeGr,
                              verbose = verbose)
    arList <- list(times = tdata$times, groupByInd = groupByInd, betaYr = betaYr,
                   yeGr = rownames(betaYr), ngroup = length(yeGr),
                   lagMatrix = lagMat$matrix, lagGroup = lagMat$group)
  }
  
  times   <- match(tdata$year, allYears)
  yrIndex <- cbind(tdata$group, tdata$tnum, times)
  colnames(yrIndex) <- c('group','tnum','year')
  
  
  seedSummary <- NULL
  if(SEEDDATA)seedSummary <- with(sdata, table(plot, year) )
  treeSummary <- with(tdata, table(plot, year) )
  plotNames   <- rownames(treeSummary)
  
  if(nspec == 1){
    fc <- .replaceString( as.character(formulaFec), 'species *','')
    fr <- .replaceString( as.character(formulaRep), 'species *','')
    formulaFec <- as.formula( paste(fc, collapse = ' ') )
    formulaRep <- as.formula( paste(fr, collapse = ' ') )
  }
  
  tunstand   <- .get.model.frame( formulaFec, tdata ) # check only, not yet standardized
 
  if(ncol(tunstand) == 1){
    checkNA <- range(tunstand)
    if(is.na(checkNA[1])){
      pmiss <- colnames(tunstand)
      if(verbose){
        cat('\nFix missing values in:\n')
        print(pmiss)
      }
    }
  }else{
    inn     <- which( !sapply(tunstand, is.factor) )
    checkNA <- sapply(tunstand[inn], range)
    wna <- which(is.na(checkNA[drop = FALSE,1,]))
    if(length(wna) > 0){
      pmiss <- paste0(colnames(checkNA)[wna],collapse=', ')
      if(verbose){
        cat('\nNote: fix missing values in these variables:\n')
        print(pmiss)
      }
    }
  }
  
  if(length(tunstand) == 0)tunstand <- numeric(0)
  tmp1  <- .get.model.frame(formulaRep, tdata)
  
  if(length(tmp1) > 0){
    
    wnew <- which(!colnames(tmp1) %in% colnames(tunstand))
    
    if(length(wnew) > 0){                               # all unique columns
      tunstand <- cbind(tunstand,tmp1[,wnew])
    }
  }
  xallNames <- colnames(tunstand)
  
  scode <- names(tunstand[ which(sapply( tunstand, is.factor )) ])
  if(length(scode) > 0){
    for(j in 1:length(scode)) tdata[,scode[j]] <- droplevels(tdata[,scode[j]])
  }
  ccode <- names(tunstand[ which(sapply( tunstand, is.character )) ])
  scode <- c(ccode, scode)
  
  specNames <- sort(unique(as.character(tdata$species)))
  
  standX <- character(0)
  xmean <- xsd <- numeric(0)
  
  wstand <- which(!colnames(tunstand) %in% scode)
  
  # standardize columns in tdata
  if(length(wstand) > 0){
    
    standX <- colnames(tunstand)[wstand]
    
    wlog <- grep( "log(", standX, fixed  = TRUE)
    if(length(wlog) > 0)standX <- standX[ -wlog ]
    wlog <- grep( "sqrt(", standX, fixed  = TRUE)
    if(length(wlog) > 0)standX <- standX[ -wlog ]
    wlog <- grep( "^2", standX, fixed  = TRUE)
    if(length(wlog) > 0)standX <- standX[ -wlog ]
    
    if(length(standX) > 0){
      treeUnstand <- tdata[,standX, drop = FALSE] # original scale
      
      xmean <- colMeans(tdata[,standX, drop = FALSE],na.rm  = TRUE)
      xsd   <- apply(tdata[,standX, drop = FALSE],2, sd, na.rm  = TRUE)
      xss   <- t( (t(tdata[,standX, drop = FALSE]) - xmean)/xsd )
      tdata[,colnames(treeUnstand)] <- xss
    }
  }
  
  tmp  <- .getDesign(formulaFec, tdata, verbose = verbose)
  xfec <- tmp$x
  if(nspec > 1)xfec <- xfec[,grep('species',colnames(xfec)),drop = FALSE]
  xfecMiss <- tmp$missx
  xfecCols <- tmp$specCols
  
  tmp  <- .getDesign(formulaRep, tdata, verbose = verbose)
  xrep <- tmp$x
  
  if(nspec > 1)xrep <- xrep[,grep('species',colnames(xrep)),drop = FALSE]
  xrepMiss <- tmp$missx
  xrepCols <- tmp$specCols
  
  rank <- qr(xfec)$rank
  if(rank < ncol(xfec)){
    for(m in 1:nspec){
      sname <- "(Intercept)"
      if(nspec > 1){
        sname <- paste('species',specNames[m],sep='')
        wk <- grep( sname, colnames(xfec) )
        wn <- which(colnames(xfec)[wk] != sname & !colnames(xfec)[wk] %in% notFit)
        wk <- wk[wn]
        snew <- colnames(xfec)[wk]
        snew <- .replaceString( snew, paste(sname,':',sep=''), '')
      }else{
        wk <- colnames(xfec)[colnames(xfec) != sname]
        snew <- colnames(xfec)[-1]
      }
      wr <- which( xfec[,sname] == 1 & z == 1 )
      
      xfecm <- xfec[wr,wk]
      rspec <- qr(xfecm)$rank
      
      if(rspec < ncol(xfecm))stop( paste(specNames[m],'not full rank') )
    }
  }
  
  if( is.null(notFit) ){
    notFit <- character(0)
  }else{
    if(nspec == 1)notFit <- .replaceString(notFit, paste('species',specNames,':',sep=''), '')
  }
  
  notFull <- character(0)
  
  diamVec <- tdata$diam
  
 # notFit <- notFit[ notFit %in% colnames(xfec) ]
  
  nff <- character(0)
  ncc <- numeric(0)
  
  
  if( ncol(xfec)/nspec > 2 ){   # more than intercept and slope
    
    for(m in 1:nspec){
      
      notw <- character(0)
      
      sname <- "(Intercept)"
      pname <- paste('species',specNames[m],sep='')
      
      if(nspec > 1){
        
        nff <- c(nff, notFit[ notFit %in% colnames(xfec) ])
        ncc <- c(ncc, match(nff, colnames(xfec)))
        
        sname <- pname
        wk <- grep( sname, colnames(xfec) )
        wn <- which(colnames(xfec)[wk] != sname & !colnames(xfec)[wk] %in% notFit)
        wk <- wk[wn]
        snew <- colnames(xfec)[wk]
        snew <- .replaceString( snew, paste(sname,':',sep=''), '')
      }else{
        wk <- colnames(xfec)[colnames(xfec) != sname]
        snew <- colnames(xfec)[-1]
        
        pname <- paste(pname, ':',colnames(xfec),sep='')
        nff   <- c(nff, notFit[ notFit %in% pname ] )
        ncc   <- c(ncc, match(nff, pname))
      }
      
      wr <- which( xfec[,sname] == 1 )
      
      xfecm <- xfec[wr,wk]
      
      # quadratic or interactions that lack main effects
      
      gw <- grep("^2)", colnames(xfecm), fixed  = TRUE)
      if(length(gw) > 0){
        
        for(i in gw){
          fi <- colnames(xfecm)[i]
          f2 <- .replaceString(fi,'I(','')
          f2 <- .replaceString(f2,'^2)','')
          if( !f2 %in% colnames(xfecm) ){
            notw <- c(notw, fi)
          }
        }
        xfecm <- xfecm[,-gw,drop=F]
        snew <- snew[ -gw ]
      }
      
      gw <- grep(":", snew, fixed  = TRUE)
      if(length(gw) > 0){
        for(i in gw){
          fi <- colnames(xfecm)[i]
          si <- columnSplit(snew[i], ':' )
          if( !si[1] %in% snew | !si[2] %in% snew )notw <- c(notw, fi)
        }
        xfecm <- xfecm[,-gw,drop=F]
        snew <- snew[ -gw ]
      }
      
      colnames(xfecm) <- snew
      cc <- suppressWarnings( cor(xfecm) )
      diag(cc) <- 0
      
      ww <- which(abs(cc) > .95 | is.na(cc), arr.ind  = TRUE)
      if(length(ww) > 0){                                  
        mvars <- unique( c(rownames(cc)[ww[,1]],colnames(cc)[ww[,2]]) )
        kvars <- rep(0, length(mvars))
        #unique values
        for(k in 1:length(mvars)){
          
          xcc <- xfecm[,!colnames(xfecm) == mvars[k]]
          xcc <- cbind( 1, xcc )
          colnames(xcc)[1] <- 'intercept'
          xcheck <- .checkDesign(xcc)
          VIF <- xcheck$VIF
          kvars[k] <- sum(VIF)
        }
        m0 <- mvars[ which.min(kvars) ]
        m1   <- paste( 'species',specNames[m], ':', m0, sep='' )
        m2   <- paste( 'species',specNames[m], ':I(', m0, '^2)', sep='' )
        notw <- c(m0, m1, m2)
        notw <- notw[ notw %in% colnames(xfec) ]
      }
      
      if( length(notw) > 0 ){
        nff <- c(nff, notw)
      }
    }
    
  }
  
 # notFit <- notFit[ notFit %in% colnames(xfec) ]
  notFit <- unique( c(notFit, nff) )
  
  
  if( length(notFit) > 0 ){
    kwords <- paste0(notFit, collapse='\n')
    kwords <- paste(' Fecundity is not full rank, omitted columns:\n',kwords)
    if(verbose)cat( paste('\nNote: ',kwords, '\n') )
    
    nff <- character(0)
    for(k in 1:length(notFit)){  # formula contains variables without species
      kk <- columnSplit(notFit[k], ':')
      vk <- kk[2]
      sk <- kk[1]
      sk <- .replaceString(sk, 'species', '')
      sr <- specNames[ !specNames == sk ]   #species not in notFit
      if(length(sr) == 0){                  #remove from formula
        ff <- as.character(formulaFec)[2]
        ff <- .replaceString(ff, paste('+',vk), '' )
        formulaFec <- as.formula( paste('~',ff, collapse = ' ') )
        nff <- c(nff, notFit[k])
      }
    }
    if(length(nff) > 0)notFit <- notFit[ notFit %in% nff ]
      
    if(length(notFit) > 0){
      notCols <- match(notFit, colnames(xfec))
      words <- paste(words, kwords)
    }
  }
  
  rank <- qr(xrep)$rank
  if(rank < ncol(xrep))stop('\nmaturation design not full rank\n')
  
  xfecU <- xfec
  xrepU <- xrep
  xfecT <- xrepT <- NULL
  
  xfecs2u <- diag(ncol(xfec))
  xreps2u <- diag(ncol(xrep))
  colnames(xfecs2u) <- rownames(xfecs2u) <- colnames(xfec)
  colnames(xreps2u) <- rownames(xreps2u) <- colnames(xrep)
  
  xfecu2s <- xfecs2u
  xrepu2s <- xreps2u
  
  if( length(xmean) > 0 ){   # unstandardized
    
    tdata[,standX] <- treeUnstand
    
    tmp <- .unstandBeta(formula = formulaFec, xdata = tdata, xnow = xfec, 
                        xmean = xmean, notCols = notCols, notFit = notFit,
                        specNames = specNames)
    xfecU   <- tmp$x          
    xfecs2u <- tmp$s2u      # multiply by beta_s to get beta_u
    xfecu2s <- tmp$u2s      # multiply by beta_u to get beta_s
    
    tmp <- .unstandBeta(formula = formulaRep, xdata = tdata, xnow = xrep, 
                        xmean = xmean, notCols = NULL, specNames = specNames)
    xrepU   <- tmp$x
    xreps2u <- tmp$s2u
    xrepu2s <- tmp$u2s
    
    xmean[abs(xmean) < 1e-10] <- 0
  }
  
  
  tdata <- cleanFactors(tdata)
  
  distTreeID <- NULL
  
  if(SEEDDATA){
    
    sdata <- cleanFactors(sdata)
    tmp   <- setupDistMat(tdata, sdata, xytree, xytrap, verbose)
    tdata <- tmp$tdata
    sdata <- tmp$sdata
    distall <- tmp$distall
    distTreeID <- tmp$distTreeID
    
    # CHECK THAT TREESONLY ARE AT END
    
    fitTrees <- unique(tdata$treeID[tdata$fit == 1])
    
    if( !'obsTrap' %in% colnames(tdata) ){           # observation period for traps
      tdata$obsTrap <- addObsTrap(tdata, sdata)
    }
    tdata$obsTrap <- tdata$obsTrap*tdata$fit         # only include fit
    
    trapRows  <- which(tdata$fit == 1)
    seedNames <- seedNames[seedNames %in% colnames(sdata)]
    
    keepCol <- c('plot','trap','trapID','year','plotYr','plotyr','drow',
                 'area','active','obs',seedNames)
    sdata   <- sdata[,keepCol]
  }
  
  tdata$species <- as.character(tdata$species)
  
  ynow    <- colnames(yrIndex)
  gy      <- columnPaste(tdata$group, yrIndex[,'year'])
  gyall   <- unique(gy)
  groupYr <- match(gy,gyall)
  if( !'dcol' %in% colnames(tdata) )tdata$dcol <- 0
  
  yrIndex <- cbind(yrIndex, tdata$dcol, groupYr)
  colnames(yrIndex) <- c(ynow, 'dcol','groupYr')
  
  
  
  specPlot  <- columnPaste( as.character( tdata$species), 
                            as.character(tdata$plot) )
  tdata$specPlot <- specPlot
  specPlots <- unique(as.character(specPlot))
  specPlot  <- match(specPlot,specPlots)
  yrIndex   <- cbind(yrIndex,specPlot)
  
  tdata$species <- as.factor(tdata$species)
  
  vif <- numeric(0)
  
  if(ncol(xfec)/nspec > 2 & verbose){
    cat('\n\nVariance Inflation Factors, range, and correlation matrix:\n')
  }
  
  for(k in 1:nspec){
    
    if(nspec == 1){
      gg <- 1:ncol(xfec)
      xt <- xfec
      xu <- xfecU
      if(ncol(xt) <= 2)next
    }else{
      sk <- paste('species',specNames[k],sep='')
      st <- paste(sk,':',sep='')
      gg <- grep( sk, colnames(xfec))
      rr <- which(xfec[,sk] == 1)
      xt <- xfec[rr,gg]
      xu <- xfecU[rr,gg]
      if(length(notFit) > 0){
        xt <- xt[,!colnames(xt) %in% notFit]
        xu <- xu[,!colnames(xu) %in% notFit]
      }
      
      colnames(xt)[1] <- 'intercept'
      colnames(xt)    <- .replaceString(colnames(xt), st, '')
    }
    
    s2 <- grep('^2',colnames(xt), fixed = TRUE)
    if(length(s2) > 0)xt <- xt[,-s2]
    
    if(ncol(xt) < 3)next
    
    xcheck <- .checkDesign(xt)
    VIF <- xcheck$VIF
    designTable <- xcheck$design
    xrange <- xcheck$range
    
    cfx <- round( xcheck$correlation,2 )
    
    rownames(cfx) <- colnames(cfx) <- names(VIF) <- .coeffNames(rownames(cfx))
    
    colnames(cfx) <- substr(colnames(cfx), 1, 6)
    cfx <- cfx[,1:(ncol(cfx)-1),drop=0]
    
    
    vif <- append(vif, list( cbind(VIF, xrange, cfx) ) )
    names(vif)[length(vif)] <- specNames[k]
    
    if(verbose){
      cat( paste('\n',specNames[k],':\n', sep='') )
      print( cbind(VIF, xrange, cfx) )
      
      if(max(xcheck$VIF, na.rm  = TRUE) > 10 | is.na(max(xcheck$VIF)))
        cat('VIF too high or undefined--too many predictors\n')
    }
  }
  
  tdata$species <- as.character(tdata$species)
  tdata$treeID  <- as.factor(tdata$treeID)
  if(!'repMu' %in% colnames(tdata))tdata$repMu <- .5
  
  tdata$repMu[!is.finite(tdata$repMu)] <- .5
  
  
  plotYears <- sort(unique( as.character(tdata$plotYr)) )
  if(SEEDDATA)plotYears <- sort(unique( c(plotYears, 
                                          as.character(sdata$plotYr))) )
  
  out <- list(tdata = tdata, z = z, zmat = zmat, zknown = zknown,
       specNames = specNames, arList = arList, plotYears = plotYears,
       plotNames = plotNames, plots = plotNames, years = years, 
       xfec = xfec, xrep = xrep, scode = scode,
       xfecMiss = xfecMiss, yeGr = yeGr,
       xrepMiss = xrepMiss, xfecCols = xfecCols, xrepCols = xrepCols,
       xmean = xmean, xsd = xsd, xfecU = xfecU, 
       xfecs2u = xfecs2u, xfecu2s = xfecu2s, xreps2u = xreps2u, xrepu2s = xrepu2s,
       xrepU = xrepU, xrepT = xrepT, specPlots = specPlots, yrIndex = yrIndex, 
       notFit = notFit, notCols = notCols,
       VIF = vif, fstart = fstart,
       fecMinCurrent = fecMinCurrent, fecMaxCurrent = fecMaxCurrent,
       matYr = matYr, last0first1 = last0first1, minDiam = minDiam,
       distTreeID = distTreeID)
  if(SEEDDATA){
    out$sdata <- sdata
    out$seedNames <- seedNames
    out$distall <- distall
    out$xytree <- xytree
    out$xytrap <- xytrap
    out$trapRows <- trapRows
  }
  out
}
 

.unstandBeta <- function(formula, xdata, xnow, xmean=NULL, notCols=NULL, notFit = NULL,
                         specNames){
  
  # xnow  - current standardized matrix
  # xdata - data.frame, unstandardized variables
  
  nspec <- length(specNames)
  
  xdata$species <- factor(xdata$species)
  
  tmp <- .get.model.frame(formula, xdata)  #unstandardized
  
  xterm <- names( tmp )
  st    <- grep('species',xterm)
  if(length(st) > 0)xterm <- xterm[-st]
  
  xfu  <- .getDesign(formula, xdata)$x    # unstandardized
  
  if(length(st) > 0 & length(xterm) > 1)
    xfu  <- xfu[,grep('species',colnames(xfu))]
  
  xuu <- xfu
  xss <- xnow  # standardized
  
  
  if( is.null(notFit) & !is.null(notCols))notFit <- colnames(xfu)[notCols]
  
  if(length(notFit) > 0){
    
    xuu <- xfu[, !colnames(xfu) %in% notFit, drop  = TRUE]
    xss <- xnow[, !colnames(xfu) %in% notFit, drop  = TRUE]
  }
  
  u2s <- solveRcpp(crossprod(xss))%*%crossprod(xss,xuu)
  s2u <- solveRcpp(crossprod(xuu))%*%crossprod(xuu,xss)
  
  list(x = xfu, s2u = s2u, u2s = u2s)
}

.blockDiag <- function(mat1,mat2){
  
  #creates block diagional
  
  if(length(mat1) == 0)return(mat2)
  
  if(!is.matrix(mat1))mat1 <- as.matrix(mat1)
  if(!is.matrix(mat2))mat2 <- as.matrix(mat2)
  
  namesc <- c(colnames(mat1),colnames(mat2))
  namesr <- c(rownames(mat1),rownames(mat2))
  
  nr1 <- nrow(mat1)
  nr2 <- nrow(mat2)
  nc1 <- ncol(mat1)
  nc2 <- ncol(mat2)
  nr  <- nr1 + nr2
  nc  <- nc1 + nc2
  
  new <- matrix(0,nr,nc)
  new[ 1:nr1, 1:nc1 ] <- mat1
  new[ (nr1+1):nr, (nc1+1):nc ] <- mat2
  colnames(new) <- namesc
  rownames(new) <- namesr
  new
}

.getLambda <- function(tdat1, sdat1, AA, ug, fz, R, SAMPR, USPEC, 
                       distance, yrs, PERAREA = FALSE, SPECPRED = F){
  
  # tdat needs species, year, dcol
  # sdat needs year, drow
  
  # PERAREA - from per-trap to per-area
  # if length(AA === 1) then it must equal 1
  # SPECPRED - predict species rather than seed types
  
  nf <- length(fz)
  
  if(SAMPR | length(R) > 1){
    if(SPECPRED){
      ff <- matrix(0, nf, nrow(R))
      jj <- match(as.character(tdat1$specPlot),rownames(R)) # predict species
      ff[ cbind(1:nf,jj) ] <- fz
    }else{
      ff <- matrix(fz,length(fz),ncol=ncol(R))*
            R[drop = FALSE,tdat1$specPlot,] # predict types
    }
  }else{
    ff <- matrix(fz,ncol=1)
  }
  
  uvec <- ug[1]
#  if(USPEC) uvec <- ug[ attr(distance,'species') ]
  
  if(USPEC){
    uvec <- matrix( ug[attr(distance, 'species')], nrow(distance), ncol(distance) )
  }
 
  dmat <- uvec/pi/(uvec + distance^2)^2
  dmat[dmat < 1e-8] <- 0
#  dmat[is.na(dmat)] <- 0
  
  plotyrs <- unique(sdat1$plotyr)

  lambda <- kernYrRcpp(dmat, ff, seedrow = sdat1[,'drow'],
                    treecol = tdat1[,'dcol'], plotyrs, 
                    treeplotYr = tdat1[,'plotyr'], seedplotYr = sdat1[,'plotyr'])
  
  if(SPECPRED){
    colnames(lambda) <- rownames(R)
    sname  <- sort(unique(attr(R,'species')))
    ii     <- rep( c(1:nrow(lambda)), ncol(lambda) )
    jj     <- match(attr(R,'species'),sname)
    jj     <- rep(jj, each=nrow(lambda))
    lambda <- .myBy(as.vector(lambda), ii, jj, fun='sum')
    colnames(lambda) <- sname
    
  }else{
    colnames(lambda) <- colnames(R)
  }
  
  if(PERAREA | length(AA) == 1) return(lambda)    # per area
  
  lambda*matrix( AA, nrow(lambda), ncol(lambda) )  # per trap
}

.tnorm <- function(n,lo,hi,mu,sig, tiny=0){   
  
  #normal truncated lo and hi
  
  if(length(lo) == 1 & length(mu) > 1)lo <- rep(lo,length(mu))
  if(length(hi) == 1 & length(mu) > 1)hi <- rep(hi,length(mu))
  
  if(length(lo) != length(mu)){
    print(length(lo))
    print(length(mu))
    stop()
  }
  
  q1 <- pnorm(lo,mu,sig)
  q2 <- pnorm(hi,mu,sig) 
  
  z <- runif(n,q1,q2)
  z <- qnorm(z,mu,sig)
  
  z[z == Inf]  <- lo[z == Inf] + tiny
  z[z == -Inf] <- hi[z == -Inf] + tiny
  z
}

rtpois <- function(lo, hi, mu){
  
  #Poisson truncated lo and hi
  
  # lo, hi, and mu are matrices
  
  ww <- which(hi > lo, arr.ind  = TRUE)  #only update where there is an interval
  xx <- lo
  
  p1 <- ppois(lo[ww],mu[ww])
  p2 <- ppois(hi[ww],mu[ww]) 
  
  xx[ ww[p2 == 0,] ] <- hi[ ww[p2 == 0,] ]
  
  vv <- which(p1 < 1)
  qq <- runif(length(vv),p1[vv],p2[vv])
  zz <- qpois(qq,mu[ww[vv,]]) - 1
  zz[zz < 0] <- 0
  
  wx <- ww[drop = FALSE,vv[zz == Inf],]
  if(length(wx) > 0)zz[zz == Inf]  <- lo[wx]
  
  wx <- ww[drop = FALSE,vv[zz == -Inf],]
  if(length(wx) > 0)zz[zz == -Inf] <- hi[wx]
  
  xx[ ww[vv,] ] <- zz
  xx
}

dtpois <- function(lo, hi, mu, index = NULL, tiny = 1e-10){

  # Poisson probability interval censusored lo and hi
  # index used for matrices
  
  xx <- lo*0 + 1
  
  if(!is.null(index)){
    hi <- hi[index]
    lo <- lo[index]
    mu <- mu[index]
    pr <-  ppois(hi, mu) - ppois(lo, mu)
    pr[pr < tiny] <- tiny
    xx[index] <- pr
    return(xx)
  }
  
  ppois(hi, mu) - ppois(lo, mu)
}

.getPlotLayout <- function(np, WIDE  = TRUE){
  
  # np - no. plots
  
  if(np == 1)return( list( mfrow=c(1,1), left=1, bottom=c(1,2) ) )
  if(np == 2){
    if(WIDE)return( list( mfrow=c(1,2), left = 1, bottom = c(1,2) ) )
    return( list( mfrow=c(2,1), left = c(1,2), bottom = 2 ) )
  }
  
  if(np == 3){
    if(WIDE)return( list( mfrow=c(1,3), left = 1, bottom = c(1:3) ) )
    return( list( mfrow=c(3,1), left = 1:3, bottom = 3 ) )
  }
  if(np <= 4)return( list( mfrow=c(2,2), left = c(1, 3), bottom = c(3:4) ) )
  if(np <= 6){
    if(WIDE)return( list( mfrow=c(2,3), left = c(1, 4), bottom = c(4:6) ) )
    return( list( mfrow=c(3,2), left = c(1,3,5), bottom = 5:6 ) )
  }
  if(np <= 9)return( list( mfrow=c(3,3), left = c(1, 4, 7), bottom = c(7:9) ) )
  if(np <= 12){
    if(WIDE)return( list( mfrow=c(3,4), left = c(1, 5, 9), bottom = c(9:12) ) )
    return( list( mfrow=c(4,3), left = c(1,4,7,10), bottom = 10:12 ) )
  }
  if(np <= 16)return( list( mfrow=c(4,4), left = c(1, 5, 9, 13), 
                            bottom = c(13:16) ) )
  if(np <= 20){
    if(WIDE)return( list( mfrow=c(4,5), left = c(1, 6, 11, 15), 
                            bottom = c(15:20) ) )
    return( list( mfrow=c(5,4), left = c(1,5,9,13), bottom = 17:20 ) )
  }
  if(np <= 25)return( list( mfrow=c(5,5), left = c(1, 6, 11, 15, 20), 
                            bottom = c(20:25) ) )
  if(np <= 25){
    if(WIDE)return( list( mfrow=c(5,6), left = c(1, 6, 11, 15, 20, 25), 
                            bottom = c(25:30) ) )
    return( list( mfrow=c(6,5), left = c(1,6,11,16,21,26), bottom = 26:30 ) )
  }
  if(np <= 36){
    return( list( mfrow=c(6,6), left = c(1, 7, 13, 19, 25, 31), bottom = c(31:36) ) )
  }
  return( list( mfrow=c(7,6), left = c(1, 7, 13, 19, 25, 31, 37), bottom = c(37:42) ) )
}

.seedProb <- function(tdat1, ug, fz, distall, sdat1, 
                      seedNames, R, SAMPR, USPEC, year1){
  
  lambda <- .getLambda(tdat1, sdat1, sdat1$area, ug, fz, R, SAMPR, USPEC, 
                       distall, year1, PERAREA = FALSE)
  lambda <- lambda + 1e-10
  ss     <- as.matrix(sdat1[,seedNames])
  dpois(ss, lambda, log  = TRUE)
}
  
.myBy <- function(x, i, j, summat=matrix(0,max(i),max(j)), 
                    totmat=summat, fun='mean'){  
  
  nn <- length(x)
  if( nn != length(i) | nn != length(j) )
    stop('\nvectors unequal in byFunctionRcpp\n')
  if( nrow(summat) < max(i) | ncol(summat) < max(j) )
    stop('\nmatrix too small\n')
  
  ww <- which(is.na(x))
  if(length(ww) > 0){
    x <- x[-ww]
    i <- i[-ww]
    j <- j[-ww]
  }
  
  frommat <- cbind(i,j,x)
  
  nr  <- nrow(frommat)
  
  maxmat <- summat*0 - Inf
  minmat <- summat*0 + Inf
  
  tmp <- byRcpp(nr, frommat, totmat, summat, minmat, maxmat)
  
  if(fun == 'sum')return(tmp$sum)
  if(fun == 'mean'){
    mu <- tmp$sum/tmp$total
    mu[is.na(mu)] <- 0
    return(mu)
  }
  if(fun == 'min'){
    return( tmp$min )
  }
  tmp$max
}

.tnormMVNmatrix <- function(avec, muvec, smat, 
                            lo=matrix(-1000,nrow(muvec),ncol(muvec)), 
                            hi=matrix(1000,nrow(muvec),ncol(muvec)),
                            whichSample = c(1:nrow(smat))){
  
  # lo, hi must be same dimensions as muvec,avec
  # each sample is a row
  
  lo[lo < -1000] <- -1000
  hi[hi > 1000]  <- 1000
  
  if(max(whichSample) > length(muvec))
    stop('\nwhichSample outside length(muvec)\n')
  
  whichSample <- sample(whichSample) # randomize order
  
  nd <- dim(avec)
  
  r <- avec
  a <- trMVNmatrixRcpp(avec, muvec, smat, 
                       lo, hi, whichSample, 
                       idxALL = c(0:(nrow(smat)-1)) ) 
  r[,whichSample] <- a[,whichSample]
  r
}

.initEM <- function(last0first1, yeGr, distall, priorU,
                    tdata, sdata, specNames, seedNames, R, 
                    SAMPR, USPEC, years, trapRows,
                    plotYears, z, xfec, fstart, verbose,
                    nsim=100){
  FSTART <- FALSE
  tiny <- 1e-4
  id <- unique(as.character(tdata$treeID))
  tt <- tdata[ match(id,as.character(tdata$treeID)), ]  #unique trees
  zobs <- tdata$repr
  
  tdata$fecMin[zobs == 0 & tdata$fecMin > tiny] <- tiny
  tdata$fecMin[zobs == 1 & tdata$fecMin < 1] <- 1 + tiny
  tdata$fecMax[zobs == 0 & tdata$fecMax > 1] <- 1
  tdata$fecMax[zobs == 1 & tdata$fecMax < (1 + tiny)] <- 1 + tiny
  
  nspec  <- length(specNames)
  tnum   <- numeric(0)
#  wdcol  <- unique(tdata$tnum[trapRows])
  
  plotyrs <- unique(tdata$plotyr[trapRows])
  ny      <- length(plotyrs)
  
  for(j in 1:ny){
    
    ws <- which(sdata$plotyr == plotyrs[j])
    wt <- which(tdata$plotyr == plotyrs[j])
    wt <- wt[wt %in% trapRows]
    
    ds <- sdata$drow[ ws ]
    dt <- tdata$dcol[ wt ]
    
    ys <- suppressWarnings( apply( distall[drop=F, ds, dt], 2, min, na.rm=T ) )
    
    close <- which(ys < 25)
    
    dj   <- wt[ close  ]
    
    if(length(dj) > 100){
      ww <- order(ys, decreasing = F)[1:100]
      dj <- dt[ww]
    }
    if(length(dj) == 0)dj <- wt
    
    if( is.na(range(dj)[1]))stop()
    
    tnum <- c(tnum, dj )
  }
  
  wtree <- which(tdata$tnum %in% tnum & zobs != 0 & tdata$fit == 1)
  
  fg <- rep(.5,nrow(tdata))
  
  ss <- sdata[,seedNames,drop = FALSE]
  ff <- sapply(ss,is.factor)
  if(!all(!ff))ss <- .fac2num( ss )
  ss <- rowSums(ss, na.rm  = TRUE)
  
  wfix <- which(is.finite(fstart))
  if(length(wfix) > 0)FSTART <- TRUE
  
  for(j in 1:length(plotyrs)){
    
    i  <- which(tdata$plotyr == plotyrs[j] )
    m  <- i[ which(i %in% wtree) ]
    
    k  <- which(sdata$plotyr == plotyrs[j])
    sj <- sum(ss[k])
    
    if(length(k) == 0)next
    if(length(i) == 0)next
    if(length(m) == 0 & sj > 0){
      m <- which(tdata$plotyr == plotyrs[j] & tdata$fit == 1)
      if(length(m) > 50)m <- sample(m, 50)
    }
    if(length(m) == 0)next
    
 #   d <- unique(tdata$tnum[m])
 #   dj <- d[d %in% tnum]
 #   ij <- m[d %in% tnum]
 #   dcol <- tdata$dcol[dj]
    
    dcol <- tdata$dcol[m]
    if(length(dj) < 1){
      ij <- which(tdata$plotyr == plotyrs[j] & tdata$tnum %in% d)
   #   dj <- tdata$tnum[ij]
      dcol <- tdata$dcol[ij]
    }
 
    dk     <- distall[ sdata[k,'drow'], dcol, drop = FALSE ]
    
    kern   <- priorU/pi/(priorU + dk^2)^2
    fg[m] <- .getF( kern, gg = ss[k]/(.1 + sdata$area[k]) )
  }
  
  fg[!is.finite(fg)] <- .1
  
  nn <- length(fg)
  
  lo <- tdata$fecMin
  hi <- tdata$fecMax
  lo[lo < tiny] <- tiny
  
  fg[fg < lo] <- lo[fg < lo]
  fg[fg > hi] <- hi[fg > hi]
  
  fg <- .tnorm(nn, lo, hi, fg, .1)
  if(FSTART)fg[wfix] <- fstart[wfix]
  fg[fg < tiny] <- tiny
  
  propF <- fg/20
  propU <- .1
  
  # by plot-yr
  sm <- matrix(0, max(c(tdata$plotyr, sdata$plotyr)), 1)
  
  pcheck <- seq(1,nsim,by=20)  
  
  if(verbose)cat("\ninitializing\n")
  
  pbar <- txtProgressBar(min=1,max=nsim,style=1)
  
  
  ug <- rep(priorU, nspec )
  names(ug) <- specNames
  nn <- length(trapRows)
  
  pall <- -1e+10
  count <- 0
  
  trapFix <- which( trapRows %in% wfix )
  
  fk <- fg[trapRows]
  zk <- z[trapRows]
  pf <- propF[trapRows]
  lk <- lo[trapRows]
  hk <- hi[trapRows]
  
  fss <- fstart[trapRows]
  fss[fss < tiny] <- tiny
  
  wwf <- which(is.finite(fss))
  
  
  for(g in 1:nsim){
    
    fnew <- .tnorm(nn, lk, hk, fk, rexp(nn,1/pf))
    fnew[trapFix] <- fk[trapFix]
    
    if(FSTART)fnew[wwf] <- fss[wwf]

    pnow <- .seedProb(tdat1 = tdata[trapRows,c('specPlot','year','plotyr','dcol')],
                      ug, fz = fk*zk, distall, sdat1 = sdata, 
                      seedNames, R, SAMPR, USPEC, years)
    pnew <- .seedProb(tdata[trapRows,c('specPlot','year','plotyr','dcol')],
                      ug, fnew*zk, distall, sdata, 
                      seedNames, R, SAMPR, USPEC, years)
    
    pnow[pnow < -1e+8] <- -1e+8   # intensity parameter is zero
    pnew[pnew < -1e+8] <- -1e+8
    
    # by plot-yr
    ii <- sdata[,'plotyr']
    ii <- rep(ii, length(seedNames))
    
    pdif <- .myBy(as.vector(pnew - pnow), ii, ii*0 + 1, summat = sm*0, fun='sum')

    if(g == 1)accept <- pdif*0
    
    a  <- exp(pdif)        #wt on seed data
    az  <- runif(length(a),0,1)
    aw  <- which(az < a)
    
    if(length(aw) > 0){
      wa <- which(tdata[trapRows,'plotyr'] %in% aw)
      fk[ wa ] <- fnew[ wa ]
      accept[aw] <- accept[aw] + 1
    }
    
    if(g %in% pcheck){
      whi <- which(accept > g/2)
      if(length(whi) > 0)pf[whi] <- pf[whi]*2
      wlo <- which(accept < g/5)
      if(length(wlo) > 0)pf[wlo] <- pf[wlo]/2
     
      pq   <- sum(pnow)
      dl   <- pq - pall
      pall <- pq
      
      
      
      
      if(dl < 0){
        count <- count + 1
        if(count > 4)break
      }
      
    }
    setTxtProgressBar(pbar,g)
  }
  fg[fg < tiny] <- tiny
  fg[trapRows] <- fk

  XX <- crossprod(xfec[wtree,])
  diag(XX) <- diag(XX) + .000001
  
  bf <- solve(XX)%*%crossprod(xfec[wtree,],log(fg[wtree]))
  mu <- exp(xfec%*%bf)
  mu[mu < lo] <- lo[mu < lo]
  mu[mu > hi] <- hi[mu > hi]
  mu[wtree] <- fg[wtree]
  if(FSTART)mu[wfix] <- fstart[wfix]
  
  
  fstart <- .tnorm(length(mu),lo,hi,mu,.1)
  fstart[wtree] <- fg[wtree]
  fstart[fstart >= .95 & z == 0] <- .95
  fstart[fstart < 1 & z == 1] <- 1.01
  fstart[fstart < 1e-4] <- 1e-4
  fstart
}

.checkDesign <- function( x, intName='intercept', xflag=':', 
                          isFactor = character(0) ){  # 
  
  # xflag - indicates that variable is an interaction
  # isFactor - character vector of factor names returned if not supplied
  
  p <- ncol(x)
  
  if(ncol(x) < 3){
    return( list(VIF = 0, correlation = 1, rank = 2, p = 2, isFactor=isFactor) )
  }
  
  if(is.null(colnames(x))){
    colnames(x) <- paste('x',c(1:p),sep='_')
  }
  xrange      <- apply(x,2,range,na.rm  = TRUE)
  wi          <- which(xrange[1,] == 1 & xrange[2,] == 1)
  if(length(wi) > 0)colnames(x)[wi] <- 'intercept'
  
  wx <- grep(xflag,colnames(x))
  wi <- which(colnames(x) == intName )
  wi <- unique(c(wi,wx))
  
  xname <- colnames(x)
  
  wmiss <- which(is.na(x),arr.ind  = TRUE)
  
  if(length(wmiss) > 0){
    rowTab <- table( table(wmiss[,1]) )
    colTab <- table(wmiss[,2])
  }
  
  VIF <- rep(NA,p)
  names(VIF) <- xname
  
  GETF <- FALSE
  if(length(isFactor) > 0)GETF <- TRUE
  
  rr <- matrix(NA, p, 2)
  
  for(k in 1:p){
    
    if(xname[k] %in% wi)next
    
    notk <- xname[xname != xname[k] & !xname %in% xname[wi]]
    ykk  <- x[,xname[k]]
    xkk  <- x[,notk,drop = FALSE]
    
    if(ncol(xkk) == 0)next
    
    wna <- which(is.na(ykk) | is.na(rowSums(xkk)))
    if(length(wna) > 0){
      ykk <- ykk[-wna]
      xkk <- xkk[-wna,]
    }
    
    ttt <- suppressWarnings( lm(ykk ~ xkk) )
    
    tkk <- suppressWarnings( summary(ttt)$adj.r.squared )
    VIF[k] <- 1/(1 - tkk)
    
    xu <- sort( unique(x[,k]) )
    tmp <- identical(c(0,1),xu)
    if(GETF)if(tmp)isFactor <- c(isFactor,xname[k])
    rr[k,] <- range(ykk, na.rm  = TRUE)
  }
  
  VIF <- VIF[-wi] 
  rr  <- rr[drop = FALSE, -wi, ]
  
  rownames(rr) <- names(VIF)
  colnames(rr) <- c('min','max')
  
  corx <- suppressWarnings( cor(x[,-wi], use="complete.obs") )
  if(length(wna) == 0){
    rankx <- qr(x)$rank
  } else {
    rankx <- qr(x[-wna,])$rank
  }
  corx[upper.tri(corx,diag  = TRUE)] <- NA
  
  findex <- rep(0,p)
  
  findex[xname %in% isFactor] <- 1
  
  designTable <- list('table' = rbind( round(VIF,2),findex[-wi],round(corx,2)) )
  rownames(designTable$table) <- c('VIF','factor',xname[-wi])
  
  designTable$table <- designTable$table[-3,]
  
  if(p == rankx)designTable$rank <- paste('full rank:',rankx,'= ncol(x)')
  if(p < rankx) designTable$rank <- paste('not full rank:',rankx,'< ncol(x)')
  
  list(VIF = round(VIF,2), correlation = round(corx,2), rank = rankx, p = p,
       isFactor = isFactor, designTable = designTable, range = rr)
}

.wrapperBeta <- function( rvp, rVPI, priorB, priorIVB, SAMPR, obsRows,
                          tdata, xfecCols, xrepCols, last0first1, ntree, nyr, 
                          betaPrior, years, YR, AR, yrIndex,
                          RANDOM, reIndex, xrandCols, RANDYR, fitCols, 
                          specNames, FECWT){
  
  function(pars, xfec, xrep, w, z, zmat, matYr, muyr){
    
    fg        <- pars$fg
    sg        <- pars$sg
    nspec     <- length(specNames)
    bgFec     <- pars$bgFec
    bgRep     <- pars$bgRep
    betaYrF   <- pars$betaYrF
    betaYrR   <- pars$betaYrR
    alphaRand <- pars$alphaRand
    Arand     <- pars$Arand
    ngroup    <- length(pars$ug)
    qr <- ncol(xrep)
    
    accept <- 0

    ONEF <- ONER <- ONEA <- FALSE
    if(ncol(xfec) == 1)ONEF <- TRUE
    if(ncol(xrep) == 1)ONER <- TRUE
    if(length(Arand) == 1)ONEA <- TRUE
    
    nxx <- length(obsRows)
    yg  <- log(fg)[obsRows]
    
    yeffect <- reffect <- 0
    
    
    xfz <- xfec[drop = FALSE,obsRows,]
    bgf <- bgFec
    
    if(FECWT)weight <- tdata$fecWt[obsRows]
    
    w0  <- which( colSums(xfz) == 0 )  
    
    apply(xfz[,fitCols], 2, sd)
    
    fitCols <- fitCols[!fitCols %in% w0]
    qf <- length(fitCols)
    xfz <- xfz[,fitCols, drop = FALSE]
    bgf <- bgf[drop = FALSE,fitCols,]
    
    if(YR){                            # yeffect in mean
      yg <- yg - betaYrF[yrIndex[obsRows,'year']] 
      if(RANDYR) yg <- yg - betaYrR[yrIndex[obsRows,c('group','year')]]
    }
    
    if(AR){
      yg <- yg - muyr[obsRows]
    }
    
    if(RANDOM){                        
      reffect <- xfec[,xrandCols]*alphaRand[reIndex,]
      if(!ONEA)reffect <- rowSums( reffect )
      yg <- yg - reffect[obsRows]
    }
    
    
    if(!is.null(betaPrior)){
      lims <- betaPrior$fec[fitCols,]
    }
    
    zrow <- z[obsRows]
    
    if(ONEF){
      
      if(FECWT){
        sw <- weight[zrow == 1]
        zz <- sum(zrow*sw)
        yy <- sum(yg[zrow == 1]*sw)
      }else{
        zz <- sum(zrow)
        yy <- sum(yg[zrow == 1])
      }
      
      V <- 1/(zz/sg + .1)
      v <- yy/sg
      bgf <- matrix( rnorm(1,V*v,sqrt(V)), 1)
      
    }else{
      
      if(FECWT){
        xw <- xfz[zrow == 1,]*weight[zrow == 1]
        yw <- yg[zrow == 1]*weight[zrow == 1]
        XX <- 1/sg*crossprod(xw) + diag(.1, qf)
        v  <- 1/sg*crossprod(xw,yw) 
      }else{
        XX <- 1/sg*crossprod(xfz[zrow == 1,]) + diag(.1, qf)
        v  <- 1/sg*crossprod(xfz[zrow == 1,], yg[zrow == 1]) 
      }
      
      testv <- try( chol(XX) ,T)
      if( inherits(testv,'try-error') ){
        diag(XX)  <- diag(XX)*1.00001
        testv <- try(chol(XX, pivot  = TRUE),T)
      }
      V  <- chol2inv(testv)
      
      if(is.null(betaPrior)){
        bgf  <- t( rmvnormRcpp(1, V%*%v, V) )
        ww   <- which( abs(bgf) > 10 )
        if(length(ww) > 0){
          bgf <- t(.tnormMVNmatrix( avec=t(bgf), muvec=t(bgf), smat=V,
                             lo=matrix(-10,1, length(bgf)), 
                             hi=matrix(10, 1, length(bgf)),
                             whichSample = ww) )
        }
      }else{
        diag(V) <- diag(V)*1.000000001
        ma   <- t(V%*%v)
        bgf  <- t(.tnormMVNmatrix( avec=ma, muvec=ma, smat=V,
                                   lo=matrix(lims[,1],1), 
                                   hi=matrix(lims[,2],1)))
      }
    }
    
    bgFec <- bgFec*0
    bgFec[fitCols,] <- bgf
    
    ww <- which(bgFec < betaPrior$fec[,1])
    vv <- which(bgFec > betaPrior$fec[,2])
    
    if(length(ww) > 0 | length(vv) > 0){
      print(bgFec)
      stop( 'bgfec error' )
    }
    
    # maturation
    if(ONER){
      
      V <- 1/(nxx + .1)
      v <- sum(w[obsRows] - .1*1)
      bgRep <- rnorm(1,V*v,sqrt(V))
      
    }else{
      
      V  <- solve(crossprod(xrep[obsRows,]) + rVPI)
      v  <- crossprod(xrep[obsRows,], w[obsRows]) + rVPI%*%rvp
      if(is.null(betaPrior)){
        bgRep <- t( rmvnormRcpp(1, V%*%v, V) )
      }else{
        lims <- betaPrior$rep
        bgRep <- t(.tnormMVNmatrix( avec=matrix(bgRep,1), muvec=t(V%*%v), 
                                    smat=V,
                                    lo=matrix(lims[,1],1), 
                                    hi=matrix(lims[,2],1)))
      }
    }
    rownames(bgFec) <- colnames(xfec) 
    rownames(bgRep) <- colnames(xrep) 
    
    list(bgFec = bgFec, bgRep = bgRep)
  }
}
    
.wrapperU <- function(distall, tdata, minU, maxU, priorU, priorVU,
                      seedNames, nspec, trapRows, obsRowSeed, obsYr, 
                      tau1, tau2, SAMPR, RANDYR, USPEC){
  
  tdat <- tdata[trapRows,c('specPlot','year','plotyr','dcol')]
  
  function(pars, z, propU, sdata){
                    
    fg <- pars$fg
    ug <- pars$ug
    umean <- pars$umean
    uvar  <- pars$uvar
    R     <- pars$R
    
    if(USPEC){
      unew <- .tnorm(nspec, minU, maxU, ug, rexp(nspec, 1/propU) )
      names(unew) <- names(ug)
    }else{
      unew <- ug
      unew[1] <- .tnorm(1, minU[1], maxU[1], ug[1], rexp(1, 1/propU) )
    }
    
    if(!RANDYR){
      umean <- priorU
      uvar  <- priorVU
    }
    
    pnow <- .seedProb(tdat1 = tdat, ug, fz = fg[trapRows]*z[trapRows], distall, 
                      sdat1 = sdata[obsRowSeed,], seedNames, 
                      R, SAMPR, USPEC, year1 = obsYr)
    pnew <- .seedProb(tdat1 = tdat, unew, fz = fg[trapRows]*z[trapRows], distall, 
                      sdat1 = sdata[obsRowSeed,], seedNames, 
                      R, SAMPR, USPEC, year1 = obsYr) 
    
    pnow <- sum(pnow) + sum( dnorm(ug, umean, sqrt(uvar),log  = TRUE) )
    pnew <- sum(pnew) + sum( dnorm(unew, umean, sqrt(uvar),log  = TRUE) )
    pdif <- pnew - pnow
    
    a <- exp(pdif)
    if(is.finite(a)){
      if( runif(1,0,1) < a){
        ug <- unew
        propU <- min(c( ug/4, propU*2) )
      }else{
        propU <- propU*.9 + .02
      }
    }
    
    if(USPEC){  # prior
      V <- 1/(nspec/uvar + 1/priorVU)
      v <- 1/uvar*sum(ug) + priorU/priorVU
      umean <- .tnorm(1, min(minU), max(maxU), V*v, sqrt(V)) 
      uvar  <- 1/rgamma(1, tau1 + nspec/2, tau2 + .5*sum( (ug - umean)^2 ) )
    }
    
    list(ug = ug, umean = umean, uvar = uvar, propU = propU)
  }
}

############ not in use:
.lambda <- function(ug, ff, zz, R, tdata, sdata, obsRows, obsRowSeed, obsYr, 
                    distall, AA, SAMPR, PERAREA = FALSE, SPECPRED = F){
  
  # PERAREA - from per-trap to per-area
  # if length(AA === 1) then it must equal 1
  # SPECPRED - predict species rather than seed types
  #AA   - area, vector or one number
  
  nf <- length(ff)
  fz <- ff*zz
  
  if( SAMPR | length(R) > 1 ){
    if(SPECPRED){
      fk <- matrix(0, nf, nrow(R))
      jj <- match(as.character(tdata$specPlot[obsRows]),rownames(R)) 
      fk[ cbind(1:nf,jj) ] <- fz
      fz <- fk
      colnames(fz) <- rownames(R)
    }else{
      fz <- matrix(fz,length(ff),ncol=ncol(R))*
        R[drop = FALSE,as.character(tdata$specPlot)[obsRows],] 
    }
  }else{
    fz <- matrix(fz,ncol=1)
  }
  
  uvec <- ug[ attr(distall,'species') ]
  
  dmat <- t(uvec/pi/(uvec + t(distall)^2)^2)
  dmat[dmat < 1e-8] <- 0
  
  plotyrs <- unique(sdata$plotyr[obsRows])
  
  lambda <- kernYrRcpp(dmat, fz, seedrow = sdata$drow[obsRowSeed],
                       treecol = tdata$dcol[obsRows], plotyrs,
                       treeplotYr = tdata$plotyr[obsRows], 
                       seedPlotYr = sdata$plotyr[obsRowSeed])
  if(SPECPRED){
    colnames(lambda) <- rownames(R)
    
    sname <- sort(unique(attr(R,'species')))
    ii <- rep( c(1:nrow(lambda)), ncol(lambda) )
    jj <- match(attr(R,'species'),sname)
    jj <- rep(jj, each=nrow(lambda))
    
    lambda <- .myBy(as.vector(lambda), ii, jj, fun='sum')
    colnames(lambda) <- sname
    
  }else{
    colnames(lambda) <- colnames(R)
  }
  
  if( PERAREA | length(AA) == 1 ) return(as.matrix(lambda))   # per area
  
  if( length(AA) > 1 )AA <- AA[obsRowSeed]
  
  as.matrix( lambda*matrix( AA, nrow(lambda), ncol(lambda) ) )  # per trap
}

.wrapperStates <- function( SAMPR, USPEC, RANDOM, SEEDDATA, obsTimes, plotYears,
                            sdata, tdat, seedNames, last0first1, distall, 
                            YR, AR, trapRows, obsRows, obsTrapRows, obsYr, predYr, obsRowSeed, 
                            ntree, years, nyr, xrandCols, reIndex, yrIndex, 
                            plag, groupByInd, RANDYR, 
                            updateProp, seedTraits, pHMC){
  
 # maxF <- specPriorVector(maxFec, tdat)
  
  function(g, pars, xfec, xrep, propF, z, zmat, matYr, muyr, 
           muran, epsilon = .00001){
    
    # tdat    - species, dcol, year, plotyr
    # pHMC    - fraction of steps that are Hamiltonian
    
    fg        <- pars$fg
    ug        <- pars$ug
    sg        <- pars$sg
    bgFec     <- pars$bgFec
    bgRep     <- pars$bgRep
    betaYrF   <- pars$betaYrF
    betaYrR   <- pars$betaYrR
    alphaRand <- pars$alphaRand
    Arand     <- pars$Arand
    R         <- pars$R
    fecMinCurrent <- tdat$fecMin
    fecMaxCurrent <- tdat$fecMax
    
    fecMinCurrent[z == 1 & fecMinCurrent < 1] <- 1
    fecMinCurrent[z == 0 & fecMinCurrent > 1e-4] <- 1e-4
    fecMaxCurrent[z == 0 & fecMaxCurrent > .999] <- .999

    
    nxx       <- length(fg)
    ngroup    <- nrow(betaYrR)
    bottom    <- -15
    
    accept  <- 0
      
    nspec  <- nrow(R)
    
    ONEF <- ONER <- ONEA    <- FALSE
    if(ncol(xfec) == 1)ONEF <- TRUE
    if(ncol(xrep) == 1)ONER <- TRUE
    if(length(Arand) == 1)ONEA <- TRUE
    
    if(AR)lindex <- 1:plag
    
    ww <- which(fg > fecMaxCurrent)
    vv <- which(fg < fecMinCurrent)
    
    if(length(ww) > 0)fg[ww] <- fecMaxCurrent[ww]
    if(length(vv) > 0)fg[vv] <- fecMinCurrent[vv]
    propF[propF > .1*fg] <- .1*fg[propF > .1*fg]
    
    yg <- log(fg)
    yg[yg < bottom] <- bottom
    
    yeffect <- reffect <- 0
    
    xfz <- xfec
    bgf <- bgFec 
    
    w0  <- which( colSums(xfz) == 0 )  
    if(length(w0) > 0){
      xfz <- xfz[,-w0]
      bgf <- bgf[drop = FALSE,-w0,]
    }
    
    if(YR & !AR){                            # yeffect in mean
      yeffect <- betaYrF[yrIndex[,'year']]
      if(RANDYR)yeffect <- yeffect + betaYrR[yrIndex[,c('group','year')]]
    }
    
    if(RANDOM){                        
      reffect <- xfec[,xrandCols]*alphaRand[reIndex,]
      if(!ONEA)reffect <- rowSums( reffect )
    }
    
    lmu           <- xfec%*%bgFec
    if(YR)lmu     <- lmu + yeffect
    if(RANDOM)lmu <- lmu + reffect
    nall <- length(fg)
    
    tt <- rbinom(1,1,pHMC)
    
    if(tt == 1 & SEEDDATA){
      
      fg[fg < 1e-6] <- 1e-6
      
      tmp <- HMC(ff = fg[obsTrapRows], fMin = fecMinCurrent[obsTrapRows], 
                 fMax = fecMaxCurrent[obsTrapRows], 
                 ep = epsilon[obsTrapRows],   
                 L = 4, tree = tdat[obsTrapRows,], sdat = sdata[obsRowSeed,], ug, 
                 mu = lmu[obsTrapRows], sg, zz = z[obsTrapRows], R, SAMPR, 
                 distance = distall, obsTrapRows, obsYr, seedNames, USPEC)
      fg[obsTrapRows] <- tmp$fg
      epsilon[obsTrapRows] <- tmp$epsilon
      
      
      ww <- which(fg > fecMaxCurrent)
      vv <- which(fg < fecMinCurrent)
      
      if(length(ww) > 0)fg[ww] <- fecMaxCurrent[ww]
      if(length(vv) > 0)fg[vv] <- fecMinCurrent[vv]
      
      fg[fg < 1e-6] <- 1e-6
      
      return( list(fg = fg, fecMinCurrent = fecMinCurrent, fecMaxCurrent = fecMaxCurrent,
                   z = z, zmat = zmat, matYr = matYr, propF = propF,
                   accept = tmp$accept, epsilon = epsilon) )
    }
    
    # maturation
    pr  <- pnorm( xrep%*%bgRep )
    iss <- sdata$plotyr
    ii  <- rep(iss, length(seedNames))
    
    if( !AR ){               
      
      tmp      <- .propZ(zmat, last0first1, matYr)
      zmatNew  <- tmp$zmat
      znew     <- zmatNew[ yrIndex[,c('tnum','year')] ] 
      
    #  znew[ tdat$repr == 1 ] <- 1
    #  znew[ tdat$repr == 0 ] <- 0
      
      zmatNew[ yrIndex[,c('tnum','year')] ] <- znew
      
      matYrNew <- tmp$matYr 
      mnow <- z*log(pr) + (1 - z)*log(1 - pr)
      mnew <- znew*log(pr) + (1 - znew)*log(1 - pr)
      
  #    dz <- znew - z
      lo <- fecMinCurrent
      hi <- fecMaxCurrent
      lo[znew == 0 & lo > 1] <- 1e-4
      lo[znew == 1 & lo < 1] <- 1
      hi[znew == 0 & hi > 1] <- 1
      hi[znew == 1 & hi < 1] <- 1
      
      
  #    lo[dz == 1] <- 1
  #    hi[dz == 1] <- tdat[,'fecMax'][dz == 1]
  #    lo[dz == -1] <- 1e-4
  #    hi[dz == -1] <- 1
      
      fnew <- .tnorm(nall, lo, hi, fg, rexp(nxx,1/propF), .001)
      ynew <- log(fnew)
      
      # fecundity model
      bnow <- bnew <- fg*0
      
      ss   <- sqrt(sg)
      
      bnow[z == 1]    <- dnorm(yg[z == 1], lmu[z == 1], ss, log  = TRUE) 
      bnew[znew == 1] <- dnorm(ynew[znew == 1], lmu[znew == 1], ss, log  = TRUE) 
      
      if( !is.null(seedTraits) ){
        
        wc <- which(z == 1 & is.finite(tdat$cropCount))
        
        if(length(wc) > 0){
          
          fc <- round(fg[wc])
          
          oss <- tdat$cropCount[wc]*seedTraits[tdat$species[wc],'seedsPerFruit']
          
          cnow <- cnew <- fc*0
          tf <- tdat$cropFraction[wc]
          ts <- tdat$cropFractionSd[wc]
          cnow <- dbetaBinom(oss, fc, tf, ts, log  = TRUE)
          
          bnow[wc] <- bnow[wc] + cnow
          
          fn <- round(fnew[wc])
          
          cnew <- dbetaBinom(oss, fn, tf, ts, log  = TRUE)
          bnew[wc] <- bnew[wc] + cnew 
        }
      }
      
      w0 <- which(z == 0 | znew == 0)
      bnew[w0] <- bnow[w0] <- 0
      
      pdif <- 0
      sm   <- matrix(0, max( tdat$plotyr ), 1)
      
      ii   <- tdat$plotyr     #consider bnow[znow == 0] = 0
      bdif <- .myBy(bnew - bnow, i = ii, j = ii*0 + 1, summat = sm*0, fun='sum')
      mdif <- .myBy(mnew - mnow, ii, ii*0+1,summat = sm, fun='sum')
      
      if(SEEDDATA){
        pnow <- pnew <- matrix(0,nrow(sdata),length(seedNames))
        pnow[obsRowSeed,] <- .seedProb(tdat[obsTrapRows,c('specPlot','year','plotyr','dcol')],
                                       ug, fg[obsTrapRows]*z[obsTrapRows], 
                                       distall, sdata[obsRowSeed,], seedNames, R, 
                                       SAMPR, USPEC, obsYr)
        pnew[obsRowSeed,] <- .seedProb(tdat[obsTrapRows,c('specPlot','year','plotyr','dcol')], 
                                       ug, fnew[obsTrapRows]*znew[obsTrapRows], 
                                       distall, sdata[obsRowSeed,], seedNames, R, 
                                       SAMPR, USPEC, obsYr)
        pnow[pnow < -1e+10] <- -1e+10   # intensity parameter is zero
        pnew[pnew < -1e+10] <- -1e+10
        
        # by plot-yr
        
        ii <- sdata$plotyr
        ii <- rep(ii, length(seedNames))
        
        sm   <- matrix(0, max( c(tdat$plotyr, sdata$plotyr) ), 1)
        pdif <- .myBy(as.vector(pnew - pnow), ii, ii*0 + 1, summat = sm*0, fun='sum')
      }
      

      a  <- exp( pdif + bdif + mdif )        
      az  <- runif(length(a),0,1)
      aw  <- which(az < a)
      
      accept <- accept + length(aw)
      
      propF <- propF/2
      
      if(length(aw) > 0){
        
        wa <- which(tdat$plotyr %in% aw)
        
        yg[ wa ] <- ynew[ wa ]
        z[ wa ]  <- znew[ wa ]
        fecMinCurrent[wa] <- lo[wa]
        fecMaxCurrent[wa] <- hi[wa]
        zmat[ yrIndex[,c('tnum','year')] ] <- z  
        
        tmp <- apply(zmat,1,which.max)
        tmp[rowSums(zmat) == 0] <- ncol(zmat)
        matYr <- tmp
        
        if(g %in% updateProp){
          propF[wa] <- propF[wa]*5
        }
        
      }else{
        propF <- propF*.5
      }
      
    }else{         # AR
       
      #independence sampler
      
      p2s <- rep(0,length(plotYears))
      
      yy <- mu <- matrix(0, ntree, nyr)
      yy[ yrIndex[,c('tnum','year')] ] <- yg
      mu[ yrIndex[,c('tnum','year')] ] <- lmu
      
      #prior for backcast based on obsRows
      
      yprior <- .myBy(yg[obsRows], yrIndex[obsRows,'tnum'], 
                      obsRows*0+1, fun='mean')
      
      for(t in 1:length(predYr)){        
        
        tii  <- which(yrIndex[,'year'] == t)  # row in tdat
        oii  <- tii[ tdat$obs[tii] == 1]      # with observations
        yii  <- yrIndex[tii,'tnum']           # location in ytmp
        nii  <- length(yii)
        fmin <- log(tdat$fecMin[tii])       # from data
        fmax <- log(tdat$fecMax[tii])
        
        zt <- zmat[yii,t]
        zp <- rbinom(length(zt),1,.5)
        
        #new and current
        ln <- lo <- fmin
        hn <- hi <- fmax
        
        wp <- last0first1[yii,'first1'] <= t | last0first1[yii,'all1'] == 1
        wn <- last0first1[yii,'last0'] >= t | last0first1[yii,'all0'] == 1
        if(t > 1)  wp <- wp | zmat[yii,t-1] == 1
        if(t < nyr)wn <- wn | zmat[yii,t+1] == 0
        zp[wp] <- 1
        zp[wn] <- 0
        
        dz <- zp - zt
        ln[dz == 1]  <- 0                #from zero to one
        hn[dz == 1]  <- fmax[dz == 1]
        ln[dz == -1] <- -fmax[dz == -1]  #from one to zero
        hn[dz == -1] <- 0
        
        mnow <- zt*log(pr[tii]) + (1 - zt)*log(1 - pr[tii])
        mnew <- zp*log(pr[tii]) + (1 - zp)*log(1 - pr[tii])
        
        a  <- exp( mnew - mnow )        
        az  <- runif(nii,0,1)
        aw  <- which(az < a)
        
        accept <- accept + length(aw)
        
        pindex <- t - (1:plag)
        w0     <- which(pindex > 0)
        mt     <- mu[,t]
        VI     <- rep(1,ntree)       # prior
        if(t > 1){                   # m_t and VI
          pindex <- pindex[w0]
          byr <- matrix(betaYrF[w0],ntree,length(w0),byrow  = TRUE)
          if(RANDYR)byr <- byr + betaYrR[groupByInd,w0,drop = FALSE]
          
          bv <- rowSums( byr*yy[,pindex] )
          mt <- mt + bv 
          
          bv <- rowSums( byr^2 )
          VI <- VI + bv
        }
        V <- sg/VI
        
        if(t > max(obsTimes)){    ###### predict forward
          
          if(length(aw) > 0){
            z[ tii[aw] ] <- zp[ aw ]
            zmat[ yii[aw],t]   <- zp[ aw ]
            fecMinCurrent[ tii[aw] ]  <- exp(ln[aw])
            fecMaxCurrent[ tii[aw] ]  <- exp(hn[aw])
          }
          yy[yii,t] <- .tnorm(nii,lo,hi,mt[yii],sqrt(sg))
          next
        }
        
        vt <- mt            # v_t
        
        for(k in 1:plag){

          tindex <- t + k - lindex
          wt  <- which(lindex != k & tindex > 0)
          
          if(length(wt) == 0)next
          
          byr1 <- matrix(betaYrF[wt],ntree,length(wt),byrow  = TRUE) 
          byr2 <- betaYrF[k]
          if(RANDYR){
            byr1 <- byr1 + betaYrR[groupByInd,wt,drop = FALSE]
            byr2 <- byr2 + betaYrR[groupByInd,k]
          }
          
          ntl <- yy[,t+k] - mu[,t+k] - rowSums( byr1*yy[,tindex[wt]] )
          vt  <- vt + ntl*byr2
        }
        vt <- vt/sg 
        
        if(t <= plag){       #### impute backward
          
          if(length(aw) > 0){
            z[ tii[aw] ] <- zp[ aw ]
            zmat[ yii[aw],t] <- zp[ aw ]
            fecMinCurrent[ tii[aw] ]  <- exp(ln[aw])
            fecMaxCurrent[ tii[aw] ]  <- exp(hn[aw])
            
          }
          V[yii]  <- 1/(1/V[yii] + 1/.5)
          vt[yii] <- vt[yii] + yprior[yii]/.5
          yy[yii,t] <- .tnorm(nii,lo,hi,(V*vt)[yii],sqrt(V[yii]))
          
          next
        }
        
        # seed data
        sii  <- which(sdata$year == years[t]) # year row in seedData
        
        ynew <- .tnorm(nii,ln,hn,(V*vt)[yii],sqrt(V[yii])) # from conditional
        
        tree <- tdat[tii,]
        seed <- sdata[sii,]
        spy  <- seed$plotyr
        
        wtt <- which(!tree$plotyr %in% spy) 
        
        if(length(wtt) > 0){         #year before seed data, draw from conditional
          if(length(aw) > 0){
            aww <- aw[aw %in% wtt]   #year before and update
            z[ tii[aww] ] <- zp[ aww ]
            zmat[ yii[aww],t] <- zp[ aww ]
            fecMinCurrent[ tii[aww] ]  <- exp(ln[aww])
            fecMaxCurrent[ tii[aww] ]  <- exp(hn[aww])
          }
          yy[yii[wtt],t] <- .tnorm(length(wtt),lo[wtt],hi[wtt],
                                   (V*vt)[yii[wtt]],sqrt(V[yii[wtt]]))
          wtk <- which(tree$plotyr %in% spy)
          if(length(wtk) > 0){
            tii <- tii[ wtk ]
            yii <- yii[ wtk ]
            
            tree <- tdat[tii,]
            mnow  <- mnow[wtk]
            mnew  <- mnew[wtk]
            ynew  <- ynew[wtk]
            zp    <- zp[wtk]
          }
        }
        
        if(length(spy) > 0){
          
          cnow <- cnew <- 0
          
          if(!is.null(seedTraits)){
            cc <- is.finite(tree$cropCount)
            wc <- which(z[tii] == 1 & cc)
            
            if(length(wc) > 0){
              
              fc <- round(fg[tii[wc]])
              
              oss <- tree$cropCount[wc]*seedTraits[tree$species[wc],'seedsPerFruit']
              
              cnow <- cnew <- fc*0
              tf   <- tree$cropFraction[wc]
              ts   <- tree$cropFractionSd[wc]
              cnow <- dbetaBinom(oss, round(fc), tf, ts, log  = TRUE)

              fn   <- round(exp(ynew[wc]))
              cnew <- dbetaBinom(oss, round(fn), tf, ts, log  = TRUE)
              ip   <- match(tree$plotYr[wc],plotYears)
              sm   <- matrix(0,max(ip),1)
              cdif <- .myBy(cnew - cnow, ip, ip*0+1,summat = sm, fun='sum')
            }
          }
          
          #########################
          
          if(SEEDDATA){
            wii <- which(tii %in% obsTrapRows)
            qii <- tii[wii]
            
            pnow <- .seedProb(tree[wii,c('specPlot','year','plotyr','dcol')],
                              ug, fg[qii]*z[qii], distall, seed,
                              seedNames, R, SAMPR, USPEC, years[t])
            pnew <- .seedProb(tree[wii,c('specPlot','year','plotyr','dcol')],
                              ug, exp(ynew[wii])*zp[wii], distall, seed,
                              seedNames, R, SAMPR, USPEC, years[t]) ###############z[tii]
            pnow[pnow < -1e+10] <- -1e+10   # intensity parameter is zero
            pnew[pnew < -1e+10] <- -1e+10
            
            iy   <- seed$plotyr
            iy   <- rep(iy,length(seedNames))
            pdif <- .myBy(as.vector(pnew - pnow), iy, iy*0 + 1, fun='sum')
            
            pyID <- unique(iy)       # plot yr for pnew/pnow
            p2s[pyID] <- p2s[pyID] + pdif[pyID]
            
          }
          
          ip <- match(tree$plotYr,plotYears)
          sm <- matrix(0,max(ip),1)
          
          mdif <- .myBy(mnew - mnow, ip, ip*0+1,summat = sm, fun='sum')
          
          myID <- unique(ip)
          p2s[myID] <- p2s[myID] +  mdif[myID] 
          if(!is.null(seedTraits) & length(cdif) > 1)
            p2s[myID] <- p2s[myID] + cdif[myID] 
          
          a  <- exp( p2s )        
          az  <- runif(length(a),0,1)
          aw  <- which(az < a)
          
          accept <- length(aw)  # no. plot-years
          
          if(length(aw) > 0){
            wa <- which(tree$plotyr %in% aw)  #rows in tdat[tii,]
            yy[ yii[wa],t] <- ynew[wa]
            z[ tii[wa] ]  <- zp[ wa ]
            fecMinCurrent[ tii[wa] ]  <- exp(ln[wa])
            fecMaxCurrent[ tii[wa] ]  <- exp(hn[wa])
            #      zmat[wa,t] <- z[ tii[wa] ]
          }
        }
      }
      
      yg <- yy[ cbind(tdat$dcol, tdat$times) ]
      tmp <- apply(zmat,1,which.max)
      tmp[rowSums(zmat) == 0] <- ncol(zmat)
      matYr <- tmp
    }
    
    fg <- exp(yg)
    
    wf <- which(propF < fg/100)
    
    if(length(wf) > 0)propF[wf] <- fg[wf]/100
    
    fecMinCurrent[fecMinCurrent < 1e-6] <- 1e-6
    fecMaxCurrent[fecMaxCurrent < 1] <- 1
    
    fg[fg > fecMaxCurrent] <- fecMaxCurrent[fg > fecMaxCurrent]
    fg[fg < fecMinCurrent] <- fecMinCurrent[fg < fecMinCurrent]
    fg[fg < 1e-6] <- 1e-6
            
    list(fg = fg, fecMinCurrent = fecMinCurrent, fecMaxCurrent = fecMaxCurrent, 
         z = z, zmat = zmat, matYr = matYr, propF = propF, 
         epsilon = epsilon, accept = accept) 
  } 
}
 
.getF <- function(kern, gg ){
  
  tiny <- .0001
  
  fec <- rep(0, ncol(kern))
  
  kk <- kern
  K   <- crossprod(kk)
  K   <- K + diag(tiny*diag(K), nrow(K), nrow(K))
  fec <- solve(K)%*%crossprod(kk, gg)
  
  fec[fec < tiny] <- tiny
  fec
}

.specFormula <- function(formula, NOINTERCEPT = FALSE){
  
  form <- paste0( as.character(formula), collapse=' ')
  form <- .replaceString(form, '~', '~ species*')
  form <- .replaceString(form, ' + ', '+ species*')
  form <- .replaceString(form, '* 1','')
  form <- .replaceString(form, '*1','')
  if(NOINTERCEPT) form <- paste(form, '-1')
  as.formula( paste(form, collapse = ' ') )
}

.getBetaPrior <- function(betaPrior, bgFec, bgRep, specNames){
  
  fecHi <- bgFec*0 + 5
  fecLo <- bgFec*0 - 5
  repHi <- bgRep*0 + 5
  repLo <- bgRep*0 - 5
  nspec <- length(specNames)
  
  if('pos' %in% names(betaPrior)){
    for(j in 1:length(betaPrior$pos)){
      jn <- character(0)
      if(nspec > 1)jn <- paste('species',specNames,':',sep='')
      if(betaPrior$pos[j] != 'intercept')jn <- paste(jn,betaPrior$pos[j],sep='')
      fecLo[rownames(fecLo) %in% jn] <- 0
      repLo[rownames(repLo) %in% jn] <- 0
    }
  }
  if('neg' %in% names(betaPrior)){
    for(j in 1:length(betaPrior$neg)){
      jn <- character(0)
      if(nspec > 1)jn <- paste('species',specNames,':',sep='')
      if(betaPrior$neg[j] != 'intercept')jn <- paste(jn,betaPrior$neg[j],sep='')
      fecHi[rownames(fecHi) %in% jn] <- 0
      repHi[rownames(repHi) %in% jn] <- 0
    }
  }
  
  # maturation intercept
 # j <- 1:nspec
 # repHi[j] <- 0

  list(fec = cbind(fecLo, fecHi), rep = cbind(repLo, repHi) )
}

.updateBetaYr <- function(yg, z, sg, sgYr, betaYrF, betaYrR, yrIndex, yeGr,
                          RANDYR, obs){
  
  #fixed effects
  
  wz <- which(z == 1 & obs == 1)
  nk <- max(yrIndex[,'year'])              # no. years
  yk <- yrIndex[wz,c('group','year')]      # year groups, years
  G  <- length(yeGr)
  bf <- betaYrF*0
  
  yfix <- yg[wz] 
 # if(RANDYR)yfix <- yfix - betaYrR[yrIndex[wz,c('group','year')]]
  
  if(!RANDYR){
    
    ygroup <- .myBy(yfix, yk[,2]*0+1, yk[,2], 
                    summat=matrix(0, 1, nk), fun='sum')
    ngr  <- .myBy(yfix*0+1, yk[,2]*0+1, yk[,2], 
                  summat=matrix(0, 1, nk), fun='sum')
    v <- ygroup/sg
    V <- 1/(ngr/sg + .1)
    bf <- matrix( .tnorm(length(v), -2, 2, V*v, sqrt(V)), 1, nk)
    bf <- bf - mean(bf)                                   # sum to zero
    
    if(!RANDYR)return( list(betaYrF = bf, betaYrR = bf*0, sgYr = sgYr, 
                            wfinite = 1:nk) )
  }
  
  # random effects
  
  ygroup <- .myBy(yfix, yk[,1], yk[,2], 
                  summat=matrix(0, G, nk), fun='sum')
  ngr  <- .myBy(wz*0+1, yk[,1], yk[,2], 
                summat=matrix(0, G, nk), fun='sum')
  
  v  <- ygroup/sg 
  V  <- 1/(ngr /sg + matrix(1/sgYr, G, nk, byrow  = TRUE))
  nc <- ngr 
  nc[nc > 1] <- 1
  ns <- colSums(nc)
  nc[,ns == 1] <- 0     # no group effect if only one group
  
  br <- matrix( .tnorm(length(v), -3, 3, V*v, sqrt(V)), G, nk )
  br <- br*nc
  
  rs <- rowSums(nc)
  ws <- which(rs > 0)
  
  br[ws,] <- sweep(br[drop = FALSE,ws,], 1, rowSums(br[drop = FALSE,ws,])/rs[ws], '-')*nc[drop = FALSE,ws,]
  
  sgYr <- 1/rgamma(nk, 1 + ns/2, 1 + .5* colSums(br^2))
  
  list(betaYrF = bf, betaYrR = br, sgYr = sgYr, wfinite = which(nc > 0))
}

.multivarChainNames <- function(rowNames,colNames){
  as.vector( t(outer(colNames,rowNames,paste,sep='_')) )
}

.updateR <- function(ug, fz, SAMPR, USPEC, distall, sdata, seedNames, 
                     tdat, R, priorR, priorRwt, years, posR, plots){
  
  mnew <- R
  mnew[posR] <- .tnorm(length(posR), 0, 1, R[posR], rexp(length(posR),50))
  
  mnew <- sweep(mnew, 1, rowSums(mnew,na.rm  = TRUE), '/')
 # mnew[-posR] <- 0

  qnow <- 2*priorRwt*log(R)
  qnow[-posR] <- 0
  qnew <- 2*priorRwt*log(mnew)
  qnew[-posR] <- 0
  
  jj <- rep(attr(R,'plot'), ncol(R))
  qnow <- tapply(as.vector(qnow), jj, sum, na.rm  = TRUE)
  qnew <- tapply(as.vector(qnew), jj, sum, na.rm  = TRUE)

  tnow <- .seedProb(tdat[,c('specPlot','year','plotyr','dcol')],
                    ug, fz, distall, sdata, seedNames,
                        R, SAMPR, USPEC, years)
  tnew <- .seedProb(tdat[,c('specPlot','year','plotyr','dcol')],
                    ug, fz, distall, sdata, seedNames,
                    mnew, SAMPR, USPEC, years)
  tnow[!is.finite(tnow)] <- -10   # intensity parameter is zero
  tnew[!is.finite(tnew)] <- -10
  
  ps <- plots[plots %in% sdata$plot]
  ii <- match(sdata$plot, ps)
  pnow <- .myBy( rowSums(tnow), ii, ii*0+1, fun='sum')[,1]
  pnew <- .myBy( rowSums(tnew), ii, ii*0+1, fun='sum')[,1]
  names(pnow) <- names(pnew) <- ps
  
  pdif <- pnew - pnow
  qdif <- qnew - qnow
  
  qdif[ps] <- qdif[ps] <- pdif
  
  a <- exp(qdif) 
  wa <- which(a > runif(length(a),0,1))
  if(length(wa) > 0){
    R[attr(R,'plot') %in% plots[wa],] <- mnew[attr(R,'plot') %in% plots[wa],]
  }
  R
}

.distmat <- function(x1,y1,x2,y2){
    xd <- outer(x1,x2,function(x1,x2) (x1 - x2)^2)
    yd <- outer(y1,y2,function(y1,y2) (y1 - y2)^2)
    t(sqrt(xd + yd)) 
}

.updateVariance <- function(yy,mu,s1=1,s2=1){
  
  ss <- (yy - mu)^2
  u1 <- s1 + length(yy)/2
  u2 <- s2 + .5*sum( ss )
  1/rgamma(1,u1,u2) 
  
}

sqrtSeq <- function(maxval){ #labels for sqrt scale
  
  # maxval on sqrt scale
  
  by   <- signif(maxval^1.7, 1)
  labs <- seq(0, maxval^2, by = by)
  at   <- sqrt(labs)

  list(at = at, labs = labs)
}

.plotObsPred <- function(obs, yMean, ySE=NULL, opt = NULL){
  
  nbin <- nPerBin <- xlimit <- ylimit <- NULL
  add <- log <- SQRT <- FALSE
  xlabel <- 'Observed'
  ylabel <- 'Predicted'
  trans <- .4
  col <- 'black'
  bins <- NULL
  atx <- aty <- labx <- laby <- NULL
  ptcol <- 'black'
  
  for(k in 1:length(opt))assign( names(opt)[k], opt[[k]] )
  
  if(!is.null(bins))nbin <- length(bins)
  
  if(log & SQRT)stop('\ncannot have both log and SQRT scale\n')
  
  yMean <- as.matrix(yMean)
  obs   <- as.matrix(obs)
  
  if(SQRT){
    if(is.null(xlimit))xlimit <- c(0, max(opt$labx))
    if(is.null(ylimit))ylimit <- c(0, max(opt$laby))
    
    xlim <- sqrt(xlimit)
    ylim <- sqrt(ylimit)
    obs   <- as.vector(sqrt(obs))
    yMean <- as.vector(sqrt(yMean))
    if(!is.null(bins))bins <- sqrt(bins)
    xlimit <- sqrt(range(obs,na.rm  = TRUE))
    xlimit[2] <- xlimit[2]*2
    ylimit <- sqrt(range(yMean,na.rm  = TRUE))
    ylimit[2] <- 1.2*ylimit[2]
 
    maxy <- max(yMean,na.rm  = TRUE)
    maxx   <- max(obs,na.rm  = TRUE)
    maxval <- max( c(maxx, maxy) )
    
    tt   <- sqrtSeq(1.2*maxx)
    if(is.null(atx))atx   <- tt$at
    if(is.null(labx))labx <- tt$labs
    
    if(ylimit[2] < xlimit[2]) ylimit[2] <- xlimit[2]
    if(xlimit[2] < xlim[2])   xlimit[2] <- xlim[2]
    if(ylimit[2] < ylim[2])   ylimit[2] <- ylim[2]
    
    tt   <- sqrtSeq(1.5*ylimit[2])
    if(is.null(aty))aty   <- tt$at
    if(is.null(laby))laby <- tt$labs

  }
    
  if(is.null(xlimit))xlimit <- range(obs)
  if(is.null(ylimit) & !add){                      # can only happen if !SQRT
    if(!log){
      plot(obs,yMean,col=.getColor(ptcol,.2),cex=.2, xlim=xlimit,
           xlab=xlabel,ylab=ylabel)
      if(log) suppressWarnings( plot(obs,yMean,col=.getColor('black',.2),cex=.3,
                                     xlim=xlimit,xlab=xlabel,ylab=ylabel,log='xy') )
    }
  }
    
  if( !is.null(ylimit) ){
    if(!log & !add){
      if(!SQRT){
        plot(obs,yMean,col=.getColor(ptcol,trans),cex=.2,
                 xlim=xlimit,xlab=xlabel,ylab=ylabel,ylim=ylimit)
      }else{
        plot(obs,yMean,col=.getColor(ptcol,trans),cex=.2,
             xlim=xlimit,xlab=xlabel,ylab=ylabel,ylim=ylimit,
             xaxt='n',yaxt='n')
        
        axis(1, at = atx, labels = labx)
        axis(2, at = aty, labels = laby, las=2)
      }
    }
    if(log & !add) plot(obs,yMean,col=.getColor(ptcol,trans),cex=.2,
                 xlim=xlimit,xlab=xlabel,ylab=ylabel,log='xy',ylim=ylimit)
  }
  if(!is.null(ySE)){
    ylo <- yMean - 1.96*ySE
    yhi <- yMean + 1.96*ySE
    for(i in 1:length(obs))lines(c(obs[i],obs[i]),c(ylo[i],yhi[i]),
                                 col='grey',lwd=2)
  }
  
  if( !is.null(nbin) | !is.null(nPerBin) ){
    
    if(is.null(bins)){
      nbin <- 20
      bins <- seq(min(obs,na.rm  = TRUE),max(obs,na.rm  = TRUE),length=nbin)
    }else{
      nbin <- length(bins)
    }
    
    if(!is.null(nPerBin)){
      nbb <- nPerBin/length(obs)
      nbb <- seq(0,1,by=nbb)
      if(max(nbb) < 1)nbb <- c(nbb,1)
      bins <- quantile(obs,nbb,na.rm  = TRUE)
      bins <- bins[!duplicated(bins)]
      nbin <- length(bins)
    }
    
    xxk <- findInterval(obs,bins)
    
    if(SQRT & is.null(bins)){
      opos <- obs[obs > 0]
      qq <- seq(0, 1, length=15)
      bins <- quantile(opos, qq)
      dbb  <- diff(bins)
      bins <- c(bins[1], bins[-1][dbb > .01])
      
      nbin <- length(bins)
      
      xxk <- findInterval(obs,bins)
      xxk[xxk == max(xxk)] <- max(xxk) - 1
    }
    xxk[xxk == nbin] <- nbin - 1
    
    wide <- diff(bins)/2
    db   <- 1
    for(k in 2:(nbin-1)){
      
      qk <- which(is.finite(yMean) & xxk == k)
      q  <- quantile(yMean[qk],c(.5,.025,.158,.841,.975),na.rm  = TRUE)
      
      if(!is.finite(q[1]))next
      if(q[1] == q[2])next
      
      ym <- mean(yMean[qk])
      xx <- mean(bins[k:(k+1)])
      rwide <- wide[k]
      
      if(k > 1)db <- bins[k] - bins[k-1]
      
      if( xx > (bins[k] + db) ){
        xx <- bins[k] + db
        rwide <- wide[ max(c(1,k-1)) ]
      }
      
      suppressWarnings(
        arrows(xx, q[2], xx, q[5], lwd=2, angle=90, code=3, col=.getColor(col,.8),
               length=.05)
      )
      lines(c(xx-.5*rwide,xx+.5*rwide),q[c(1,1)],lwd=2, 
            col=.getColor(col,.8))
      rect(xx-.4*rwide,q[3],xx+.4*rwide,q[4], col=.getColor(col,.5), border=col)
    }
  }
  invisible( list(atx = atx, labx = labx, aty = aty, laby = laby) )
}

.getKern <- function(u,dij){
  
  uvec <- u[ attr(dij,'group') ]
  kk <- t(uvec/pi/(uvec + t(dij)^2)^2)
  
 # kk <- u/pi/(u + dij^2)^2
  kk[is.na(kk)] <- 0
  kk
}

.mapSpec <- function(x, y, z, mapx=range(x), mapy=range(y), scale=0,
                     add = FALSE, sym='circles',
                     colVec=rep(1,length(x)), fill = FALSE){
  
  fillCol <- NA
  if(is.logical(fill))fillCol <- colVec
  if(is.character(fill))fillCol <- fill
  
  opin <- par()$pin
  
  if(scale > 0).mapSetup(mapx,mapy,scale)
  if(!add){
    plot(NA, xlim=mapx, ylim=mapy, axes = F, xlab='', ylab='')
    Axis(side=1, labels = FALSE)
    Axis(side=2, labels = FALSE)
    add <- TRUE
  }
  
  if(sym == 'circles'){

    symbols(x,y,circles=z/10,inches = FALSE,
                              xlim=mapx,ylim=mapy,fg=colVec,bg=fillCol,
                              lwd=2,add=add)
  }
  if(sym == 'squares'){
    symbols(x,y,squares=z/10,inches = FALSE,
                              xlim=mapx,ylim=mapy,fg=colVec,bg=fillCol,
                              lwd=2,add=add)
  }
  par(pin = opin)
}

scaleBar = function(label, value = 1, fromLeft = .5, yadj = .1, 
                    lwd = 3, cex = 1) {
  
  xl <- par("usr")[1:2]
  yl <- par("usr")[3:4]
  
  xm <- xl[1] + fromLeft*diff(xl)
  x1 <- xm - value/2
  x2 <- xm + value/2
  
  y  <- yl[1] + .05*diff(yl)
  ym <- y + yadj*diff(yl)
    
  lines(c(x1,x2),c(y,y), lwd=lwd + 2, col='white')
  lines(c(x1,x2),c(y,y), lwd=lwd)
  
  lab <- paste(value, label)
  text(xm, ym, lab, cex = cex)
}

.mapSetup<- function(xlim,ylim,scale){  #scale is m per inch

  px   <- diff(xlim)/scale
  py   <- diff(ylim)/scale
  pin  <- c(px, py)
  par(pin = pin)
  invisible( pin )
}

.getColor <- function(col,trans){
  
  # trans - transparency fraction [0, 1]
  
  tmp <- col2rgb(col)
  rgb(tmp[1,], tmp[2,], tmp[3,], maxColorValue = 255, 
      alpha = 255*trans, names = paste(col,trans,sep='_'))
}

.interp <- function(y,INCREASING = FALSE,minVal=-Inf,maxVal=Inf,defaultValue=NULL,
                   tinySlope=NULL){  #interpolate vector x
  
  if(is.null(defaultValue))defaultValue <- NA
  
  tiny <- .00001
  if(!is.null(tinySlope))tiny <- tinySlope
  
  y[y < minVal] <- minVal
  y[y > maxVal] <- maxVal
  
  n  <- length(y)
  wi <- which(is.finite(y))
  
  if(length(wi) == 0)return(rep(defaultValue,n))
  if(length(wi) == 1)ss <- tiny
  
  xx  <- c(1:n)
  z  <- y
  
  if(wi[1] != 1) wi <- c(1,wi)
  if(max(wi) < n)wi <- c(wi,n)
  
  ss <- diff(z[wi])/diff(xx[wi])
  
  ss[is.na(ss)] <- 0
  
  if(length(ss) > 1){
    if(length(ss) > 2)ss[1] <- ss[2]
    ss[length(ss)] <- ss[length(ss)-1]
  }
  if(INCREASING)ss[ss < tiny] <- tiny
  
  if(is.na(y[1]))  z[1] <- z[wi[2]] - xx[wi[2]]*ss[1]
  if(z[1] < minVal)z[1] <- minVal
  if(z[1] > maxVal)z[1] <- maxVal
  
  for(k in 2:length(wi)){
    
    ki <- c(wi[k-1]:wi[k])
    yk <- z[wi[k-1]] + (xx[ki] - xx[wi[k-1]])*ss[k-1]
    yk[yk < minVal] <- minVal
    yk[yk > maxVal] <- maxVal
    z[ki] <- yk
  }
  z
}

.interpRows <- function(x, startIndex=rep(1,nrow(x)), endIndex=rep(ncol(x),nrow(x)),
                       INCREASING = FALSE, minVal=-Inf, maxVal=Inf,
                       defaultValue=NULL,tinySlope=.001){  
  #interpolate rows of x subject to increasing
  
  nn  <- nrow(x)
  p  <- ncol(x)
  xx <- c(1:p)
  
  if(length(minVal) == 1)minVal <- rep(minVal,nn)
  if(length(maxVal) == 1)maxVal <- rep(maxVal,nn)
  
  ni   <- rep(NA,nn)
  flag <- numeric(0)
  
  z <- x
  
  for(i in 1:nn){
    if(startIndex[i] == endIndex[i]){
      z[i,-startIndex[i]] <- NA
      next
    }
    z[i,startIndex[i]:endIndex[i]] <- .interp(x[i,startIndex[i]:endIndex[i]],
                                             INCREASING,minVal[i],maxVal[i],
                                             defaultValue,tinySlope)
  }
  
  z
}

.shadeInterval <- function(xvalues,loHi,col='grey',PLOT  = TRUE,add  = TRUE,
                           xlab=' ',ylab=' ', xlim = NULL, ylim = NULL, 
                           LOG = FALSE, trans = .5){
  
  #draw shaded interval
  
  tmp <- smooth.na(xvalues,loHi)

  xvalues <- tmp[,1]
  loHi <- tmp[,-1]
  
  xbound <- c(xvalues,rev(xvalues))
  ybound <- c(loHi[,1],rev(loHi[,2]))
  if(is.null(ylim))ylim <- range(as.numeric(loHi))
  if(is.null(xlim))xlim <- range(xvalues)
  
  if(!add){
    if(!LOG)plot(NULL, xlim = xlim, ylim=ylim, 
                 xlab=xlab, ylab=ylab)
    if(LOG)suppressWarnings( plot(NULL,  xlim = xlim, ylim=ylim, 
                xlab=xlab, ylab=ylab, log='y') )
  }
 
  
  if(PLOT)polygon(xbound,ybound, border=NA,col=.getColor(col, trans))
  
  invisible(cbind(xbound,ybound))
  
}

smooth.na <- function(x,y){   
  
  #remove missing values
  #x is the index
  #y is a matrix with rows indexed by x
  
  if(!is.matrix(y))y <- matrix(y,ncol=1)
  
  wy <- which(!is.finite(y),arr.ind   = TRUE)
  if(length(wy) == 0)return(cbind(x,y))
  wy <- unique(wy[,1])
  ynew <- y[-wy,]
  xnew <- x[-wy]
  
  return(cbind(xnew,ynew))
}

.boxplotQuant <- function( xx, ..., boxfill=NULL, omit.na = TRUE ){
  
  pars <- list(...)
  
  #print(pars)
  
  q    <- pnorm(c(-1.96, -1, 0, 1, 1.96))
  qfec <- apply( xx, 2, quantile, q, na.rm=T )
  
  wf   <- which( is.finite(qfec[1,]) )
  if(omit.na){
    xx <- xx[,wf, drop=F]
    qfec <- qfec[,wf, drop=F]
    if('border' %in% names(list))border = border[wf]
    if('whiskcol' %in% names(list))whiskcol = whiskcol[wf]
    if('boxfill' %in% names(list))boxfill = boxfill[wf]
  }else{
    qfec[,wf, drop=F] <- 0
  }
  
 # print('i2')
  
  tmp <- boxplot( xx, ..., na.rm=T, plot = FALSE)
 # ss  <- apply( xx, 2, quantile, q), na.rm  = TRUE ) 
  tmp$stats <- qfec
  

  if( 'col' %in% names(pars) )boxfill <- pars$col
  
  bxp( tmp, ..., boxfill = boxfill )
  
  invisible( tmp )
}

.fitText2Fig <- function(xx, width  = TRUE, fraction=1, cex.max=1){
  
  # returns cex to fit xx within fraction of the current plotting device
  # width - horizontal labels stacked vertically
  #!width - vertical labels plotted horizontally
  
  px <- par('pin')[1]
  py <- par('pin')[2]
  cl <- max( strwidth(xx, units='inches') )
  ch <- strheight(xx, units='inches')[1]*length(xx)  # ht of stacked vector
  
  if(width){              #horizontal labels stacked vertically
    xf <- fraction*px/cl
    yf <- fraction*py/ch
  } else {                #vertical labels plotted horizontally
    xf <- fraction*px/ch
    yf <- fraction*py/cl
  }
  
  cexx <- min(c(xf,yf))
  if(cexx > cex.max)cexx <- cex.max
  cexx
}



