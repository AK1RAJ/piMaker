piMaker_theme <- theme(
  plot.title = element_text(size=15),
  axis.text.x = element_text(size=10),
  axis.text.y = element_text(size=10),
  axis.title.x = element_text(size=10),
  axis.title.y = element_text(size=10),
  panel.background = element_rect(fill = NA),
  axis.line = element_line(linetype = "solid", linewidth = 1),
  legend.background = element_rect(fill = NA),
  legend.key = element_rect(fill = NA),
  legend.position = "top"
)

#count matrix
countMatrix <- function(x, Count = parent.frame()$Count){
  pat <- x
  count <- Count
  #make a blank matrix to store the values
  p <- matrix( nrow = length(count) )
  dimnames(p)<- list(count)
  #set values to 0
  var <- c(1:length(count))
  p[,1] <- c(0)
  
  #count 
  for(i in 1:length(count)){
    
    var[i] <- sum(pat == count[i], na.rm = T)
    
  }
  
  p <- cbind(p, var)
  p <- p[,-1]
  
}

#this function builds a matrix for each input position, with each nucleotide covered by Length-1
#the minus 1 is essential as the first covered position is always given by the 5' start position.
coverMatrix <- function(x, Length){
  dat <- x
  l <- Length
  l <- l-1
  #get the data from the list
  #remove the unwanted bits
  dat <- dat[,c(4,5)]
  #add n to each position to make a matrix of covered nucleotides
  for (i in 1:(l)) {
    n = i
    x = dat$pos
    dat[paste0(n)] <- (x + i)
  }
  return(dat)
}
#this function builds a matrix for each input position of piRNAs, with each nucleotide covered by Length-1
#the minus 1 is essential as the first covered position is always given by the 5' start position.
coverPiMatrix <- function(x){
  In. <- dat
  #get the data from the list
  #remove the unwanted bits
  In <- In.[,c(4,6,5)]
  
  res <- c()
  
  for(I in 1:nrow(In)){
    p <- c()
    Start <- In$pos[I]
    Length <- In$qwidth[I]
    
    for (i in 1:(Length)) {
      n = i
      
      p[paste0(n)] <- (Start + i)
    }
    
    res <- bind_rows(res,p)
  }
  
  In <- bind_cols(In,res)
  
  return(In)
}

#read in files if multiple
filesInMultiple <- function(files = NULL, tidyname = TRUE, what = parent.frame()$what){
  #reads in the files according to the parameters 'what'
  filename <- names(files[i])
  ifelse(tidyname, namefile <- tidyName(filename), namefile <- filename)
  data <- as.data.frame(scanBam(paste(filename)), param = Param)
  return( namefile <- data )
}      
#read in files if single
filesInSingle <- function(files = NULL, tidyname = TRUE, what = parent.frame()$what){
  
  filename <- names(files)
  ifelse(tidyname, namefile <- tidyName(filename), namefile <- filename)
  data <- as.data.frame(scanBam(paste(filename)), param = Param)
  return( namefile <- data )
}
#read in files
filesIn <- function(files = NULL, tidyname = TRUE, what = c("qname", "rname", "strand", "pos", "qwidth", "seq")){
  if(length(files)>1){
    files <- filesInMultiple(files)
  }
  if(length(files)==1){
    files <- filesInSingle(files)
  }
}

#calculates the coverage of each nucleotide in the genome
getCoverage <- function(x, GenSize){
  res <- x
  Len = GenSize
  Count <- GenSize
  #filter the positives
  datSP <- split(res,res$strand)
  datSP <- datSP[-3]
  #make the size matrix
  for(sp in 1:length(datSP)){
    dSP <- datSP[[sp]]
    if(names(datSP[sp]) == "+"){
      tabP <- countMatrix(dSP, Count = Len)
    }else{
      tabN <- -countMatrix(dSP, Count = Len) 
    }
  }
  tab <- bind_cols(tabP,tabN)
  colnames(tab) <- c("Pos", "Neg")
  tab$x <- as.numeric(rownames(tab))
  return(data.frame(tab))
}

getGenLength <- function(x){
  dat = x
  gm <- aggregate(pos~rname, dat, FUN = max)
  ifelse(exists("readLength"), (readLength <- rbind(readLength, gm)), (readLength <- gm))
  readLength <- aggregate(.~rname, readLength, FUN = max)
  return(readLength)
}

#get the size distribution of mapped reads
getSize <- function(x, Count){
  dat <- x
  count <- Count
  #filter the positives
  datP <- dplyr::filter(dat, dat$strand == "+")
  #make the size matrix
  tabP <- countMatrix(datP$qwidth)
  #filter the negatives
  datN <- dplyr::filter(dat, dat$strand == "-")
  #get the size matrix
  tabN <- -countMatrix(datN$qwidth)
  #combine the positive and negative
  tab <- bind_rows(tabP,tabN)
  tab <- data.frame(tab)
  row.names(tab) <- c("Pos", "Neg")
  colnames(tab) <- c(count)
  return(tab)
}

#function to get the genome length
makeRsq <- function(x){
  nam <- x
  #retrieve the length from either input or the BAM file
  rSq <- 
    if(exists("GenLength")){
      rSq <-   GenLength[grepl(nam, GenLength$filename),]
    }else{
      rSq <-   readLength[grepl(nam, readLength$filename),]
    }
  rsq <- rSq$dat
  return(rsq)
}

#get max read count for plots
maxCount <- function(Min, Max){
  cPos <- Max
  cNeg <- Min
  c <- cbind(max(cPos),max(abs(cNeg)))
  colnames(c) <- c("Pos", "Neg")
  return(c)
}

#Mean and SD calculations
meanCalc <-function(x){
  Sumry <- x
  pos <- Sumry[grepl("Pos", colnames(Sumry))]
  Sumry$Pos_Mean <-  rowMeans(pos)
  Sumry$Pos_SD <- apply(pos[,-(ncol(pos))],1, sd)
  Sumry$Pos_E_min <- (Sumry$Pos_Mean - Sumry$Pos_SD)
  Sumry$Pos_E_max <- (Sumry$Pos_Mean + Sumry$Pos_SD)
  neg <- Sumry[grepl("Neg", colnames(Sumry))]
  Sumry$Neg_Mean <-  rowMeans(neg)
  Sumry$Neg_SD <- apply(neg[,-(ncol(neg))],1, sd)
  Sumry$Neg_E_min <- (Sumry$Neg_Mean - Sumry$Neg_SD)
  Sumry$Neg_E_max <- (Sumry$Neg_Mean + Sumry$Neg_SD)
  Sumry$x <- 1:nrow(Sumry)
  return(Sumry)
}

#normalising data
Normalise <- function(x){
  data <-x
  ifelse(sign(data)>0, ((data - min(data))/(max(data) - min(data))), 
         -((abs(data) - (min(abs(data))))/(max(abs(data)) - min(abs(data)))))
}
#normalising data to a given max value
NormaliseTo <- function(x, minvalue, maxvalue){
  data <-x
  ifelse(sign(data)>0, ((data - minvalue)/(maxvalue - minvalue)), 
         -((abs(data) - (abs(minvalue)))/(abs(maxvalue) - abs(minvalue))))
}

overlapMatrixMake<- function(x, GenLength, Count){
  dat <- x
  GenPosn <- GenLength
  Count <- Count
  for(t in 1:length(GenPosn)){
    targ <- data.frame(GenPosn[t])#sets a target for a single nucleotide in the genome
    pat <- (dplyr::filter(dat, dat$pos == targ))#returns all reads at that target position
    if(is.na(pat[1,1])){
      pat[1,] <-  c(paste0(rSeq), paste0(spt), (paste0(GenPosn[t])),(0L),(0L))
    }#sets a grid of 0 if there are no reads at the target position
    pat <- pat[,3:4]#removes excess input
    res <- (countMatrix(pat))#returns sum of the count of read length at the target
    res <- data.frame(t(res))
    row.names(res) <- targ
    res$pos <- as.numeric(row.names(res))
    res <- res[, c( ncol(res), (1:(ncol(res)-1))  )]
    ifelse(exists("result"), result <- dplyr::bind_rows(result, res), result <- res)
  }
  return(result)
  rm(result)
}

overLaps <- function(x, Overlap){
  #get input
  Lap <-  x
  lap <- data.frame(Lap[, grepl("By", colnames(Lap))])
  oLap <- c()
  #count the overlap totals from the corresponding column
  for (o in 1:length(Overlap)){
    
    v <- sum(lap[,o], na.rm = TRUE)
    
    oLap[paste0("by_",o)] <- as.numeric(v)
    
  }
  #arrange the data
  oLaps <- as.data.frame(as.numeric(oLap))
  oLaps$x <- as.numeric(row.names(oLaps))
  oLaps <- oLaps[,c(2,1)]
  colnames(oLaps) <- c("x", "count")
  #return the data
  return(oLaps)
}

piCatcher <- function(x, Length_In , Target_Length , Overlap, Split){
  
  input <- x
  lap <-c()
  #split the data into the opposing strands
  datIn <- data.frame(input[, grepl(paste0(Split), colnames(input))])
  datIn <- data.frame(datIn[, names(dplyr::select(datIn, matches(c( as.character(Length_In)))))])#makes sure we only keep the lengths we want
  
  #make two frames from each strand with the total reads starting at each position
  datIn$Total <- rowSums(datIn)
  datIn$pos <- row.names(datIn)
  lap <- datIn[, c( (ncol(datIn)), (ncol(datIn)-1) ), drop = F]
  row.names(lap) <- lap$pos
  
  #split the data into the opposing strands
  datO <- data.frame(input[, !grepl(paste0(Split), colnames(input))])
  datO <- data.frame(datO[,names(dplyr::select(datO, matches(c( as.character(Target_Length)))))])#only keep the targets
  
  #make two frames from each strand with the total reads starting at each position
  datO$Total <- rowSums(datO)
  datO$pos <- row.names(datO)
  datO <- datO[, c( (ncol(datO)), (ncol(datO)-1) ), drop = F]
  
  result <- c()
  
  for (i in 1:nrow(lap)){
    
    max <- lap[i,2]
    
    #find the range of rows which are above 0 
    range <- ( ( paste(i-length(Overlap)) ) : ( paste(i) ) )
    ind <- which( ( paste(i-length(Overlap)) ) : ( paste(i) ) > 0)
    #set the target indexes     
    targ <- range[ind]
    #get the target    
    target <- datO[c( targ ) , "Total"]
    var <- c()
    #count the overlaps    
    
    for(l in 1:length(Overlap)){
      
      var[l] <- target[(l)]
      #this now limits the overlaps counted to the number of input sequences 
      
      if(is.na(var[l])){var[l] <- NA
      }else{
        if((var[l]) == 0 ){var[l] <- 0
        }else{
          if(var[l] > max){var[l] <- max}
        }
      }
    }
    #return the result    
    res <- t(data.frame(var))
    colnames(res) <- paste0("By_", 1:length(Overlap))
    rownames(res) <- i
    
    ifelse( ( exists("result") ), result <- rbind( result, res ), result <- res) 
    
  }
  
  
  Lap <- cbind(lap, result)
  rm(result, Target_Length, Length_In, Overlap)
  return(Lap)
} 

#gets the frequency of piRNAs by position
piFreq <- function(x){
  datr <- x
  datr <- datr
  datrPos <- dplyr::filter(datr, datr$strand == "+")
  datrNeg <- dplyr::filter(datr, datr$strand == "-")
  datrNeg$FivPr <- datrNeg$pos + datrNeg$qwidth
  sumPos <- as.data.frame(xtabs(~pos, data = datrPos))
  row.names(sumPos) <- sumPos$pos
  sumPos$x <- as.numeric(row.names(sumPos))
  sumNeg <- as.data.frame(xtabs(~FivPr, data = datrNeg))
  row.names(sumNeg) <- sumNeg$FivPr
  sumNeg$x <- as.numeric(row.names(sumNeg))
  sumNeg$Freq <- -sumNeg$Freq
  plot <- as.data.frame(c(-1 : rsq+1))
  row.names(plot) <- plot[,1]
  plot <- merge(plot, sumPos, by.x = 1, by.y =0, all.x = TRUE, all.y = TRUE)
  row.names(plot) <- plot[,1]
  plot <- plot[,c(1,3)]
  plot[is.na(plot)] <- 0 
  colnames(plot) <- c("Pos", "Freq")
  plotNeg <- as.data.frame(c(-1 : rsq+1))
  row.names(plotNeg) <- plotNeg[,1]
  plotNeg <- merge(plotNeg, sumNeg, by.x = 1, by.y =0, all.x = TRUE, all.y = TRUE)
  plotNeg <- plotNeg[,c(1,3)]
  row.names(plotNeg) <- plotNeg[,1]
  plotNeg[is.na(plotNeg)] <- 0
  colnames(plotNeg) <- c("Neg", "Freq")
  res <- merge(plot, plotNeg, by.x = 1, by.y = 1, all.x = TRUE, all.y = TRUE)
  colnames(res) <- c("pos", "Pos", "Neg")
  return(res)
}

#summarise the data for plots
SizDistSum <- function(x){
  res <- x
  #get the positives
  resPos <- res[grepl(("_Pos"), colnames(res))]
  #make means
  resPos$Pos_Mean <- rowMeans(resPos)
  #make SD from all columns but the last(means)
  resPos$Pos_SD <- apply(resPos[,1:(ncol(resPos)-1)], 1, sd, na.rm = T)
  #error bar minimum
  resPos$Pos_E_Min <- resPos$Pos_Mean - resPos$Pos_SD
  #error bar maximum
  resPos$Pos_E_Max <- resPos$Pos_Mean + resPos$Pos_SD
  #remove the input data
  resPos <- resPos[, c("Pos_Mean","Pos_SD","Pos_E_Min","Pos_E_Max")]
  #repeat for the negative values
  resNeg <- res[grepl(("_Neg"), colnames(res))]
  resNeg$Neg_Mean <- rowMeans(resNeg)
  resNeg$Neg_SD <- apply(resNeg[,1:(ncol(resNeg)-1)], 1, sd, na.rm = T)
  resNeg$Neg_E_Min <- resNeg$Neg_Mean - resNeg$Neg_SD
  resNeg$Neg_E_Max <- resNeg$Neg_Mean + resNeg$Neg_SD
  resNeg <- resNeg[, c("Neg_Mean","Neg_SD","Neg_E_Min","Neg_E_Max")]
  #bind the positive and negative together
  resB <- bind_cols(resPos, resNeg)
  #add x values for plot
  resB$x <- rownames(resB)
  return(resB)
}

#rename odd filenames from RNA Seq data
tidyName <- function(x){
  x2 <- ifelse(is.na(str_extract(x,
                                 '_[A-Za-z].*[0-9]+[A-Za-z]')),
               sub('\\..[^\\.]*$', '', x),
               str_extract(x,
                           '_[A-Za-z].*[0-9]+[A-Za-z]') %>% 
                 str_extract('[A-Za-z].*[0-9]'))
  return(x2)
}

#calculate the z-score
Z_Score <- function(x){
  
  dat <- x 
  dat$Z <- ((dat[,2]- 
               mean(dat[,2]))/
              sd(dat[,2]))
  return(dat)
}

#compact code to save images from plots
saveImage <- function(x,  SVG = T, PNG = T, scale = 1, width = 14, height = 14, units = "cm"){
  
  scale. = scale
  width. = width
  height. = height
  units. = units
  
  
  
  if (savefiles == T){
    
    filename = x
    directory = OUT
    
    if(SVG == T){
      ggsave(paste0(directory, "/",filename, ".svg"), plot = last_plot(), device = svg, scale = scale., width = width., height = height., units = units., 
             dpi = 600, limitsize = TRUE, bg = NULL)
    }
    
    if(PNG == T){
      ggsave(paste0(directory, "/",filename, ".png"), plot = last_plot(), device = png, scale = scale., width = width., height = height., units = units., 
             dpi = 600, limitsize = TRUE, bg = NULL)
    }
    
  }
}

#calculate overlap probability per position
piProb <- function(x, .Overlap = parent.frame()$Overlap) {
  dat <- x
  max <- ncol(dat)
  min <- max - length(.Overlap)
  range <- ((min+1):(max))
  res <- dat[,c(1,2)]
  res$Lap_total <- rowSums(dat[,c(3:length(.Overlap))], na.rm=TRUE)
  
  for (i in 1:length(range)){
    
    pos <- range[i]
    nam <- paste0("By_",.Overlap[i])
    res[[nam]] <- ifelse(res$Lap_total > 0,
                         ifelse( is.na(dat[,(pos)]/res$Lap_total),0, dat[,(pos)]/res$Lap_total ),
                         0
    )
  }
  
  return(res)
}

piProbSum <- function(x, .Overlap = parent.frame()$Overlap){ #clone of overLaps function but 
  #calculates probability by dividing input by total
  Lap <-  x
  lap <- data.frame(Lap[, grepl("By", colnames(Lap))])
  oLap <- c()
  #count the overlap totals from the corresponding column
  for (o in 1:length(.Overlap)){
    
    v <- sum(lap[,o], na.rm = TRUE)
    
    oLap[paste0("by_",o)] <- as.numeric(v)
    
  }
  #arrange the data
  oLaps <- as.data.frame(as.numeric(oLap))
  oLaps$x <- as.numeric(row.names(oLaps))
  oLaps <- oLaps[,c(2,1)]
  colnames(oLaps) <- c("x", "Probability")
  oLaps$Probability <- oLaps$Probability/sum(oLaps$Probability)
  #return the data
  return(oLaps)
}

#calculate weighted overlap probability per position
piProbWeighted <- function(x, .Overlap = parent.frame()$Overlap) {
  dat <- x
  max <- ncol(dat)
  min <- max - length(.Overlap)
  range <- ((min+1):(max))
  res <- dat[,c(1,2)]
  res$Lap_total <- rowSums(dat[,c(3:length(.Overlap))], na.rm=TRUE)
  
  res$weighting <- res$Total/sum(res$Total)
  
  for (i in 1:length(range)){
    
    pos <- range[i]
    nam <- paste0("By_",.Overlap[i])
    res[[nam]] <- ifelse(res$Lap_total > 0,
                         ifelse( is.na(dat[,(pos)]/res$Lap_total),0, (dat[,(pos)]/res$Lap_total)*res$weighting ),
                         0
    )
  }
  
  return(res)
}

piCatcherPos <- function(x, .Length_In = parent.frame()$Length_In, 
                         .Target_Length = parent.frame()$Target_Length, 
                         .Overlap = parent.frame()$Overlap, 
                         Split = "Pos"){
  
  input <- x
  lap <-c()
  #split the data into the opposing strands
  datIn <- data.frame(input[, grepl(paste0(Split), colnames(input))])
  datIn <- data.frame(datIn[, names(dplyr::select(datIn, matches(c( as.character(.Length_In)))))])#makes sure we only keep the lengths we want
  
  #make two frames from each strand with the total reads starting at each position
  datIn$Total <- rowSums(datIn)
  datIn$pos <- row.names(datIn)
  lap <- datIn[, c( (ncol(datIn)), (ncol(datIn)-1) ), drop = F]
  row.names(lap) <- lap$pos
  
  #split the data into the opposing strands
  datO <- data.frame(input[, !grepl(paste0(Split), colnames(input))])
  datO <- data.frame(datO[,names(dplyr::select(datO, matches(c( as.character(.Target_Length)))))])#only keep the targets
  
  #make two frames from each strand with the total reads starting at each position
  datO$Total <- rowSums(datO)
  datO$pos <- row.names(datO)
  datO <- datO[, c( (ncol(datO)), (ncol(datO)-1) ), drop = F]
  
  result <- c()
  
  for (i in 1:nrow(lap)){
    
    max <- lap[i,2]
    
    #find the range of rows which are above 0 
    range <- ( ( paste(i-length(.Overlap)) ) : ( paste(i) ) )
    ind <- which( ( paste(i-length(.Overlap)) ) : ( paste(i) ) > 0)
    #set the target indexes     
    targ <- range[ind]
    #get the target    
    target <- datO[c( targ ) , "Total"]
    var <- c()
    #count the overlaps    
    
    for(l in 1:length(.Overlap)){
      
      var[l] <- target[(l)]
      #this now limits the overlaps counted to the number of input sequences 
      
      if(is.na(var[l])){var[l] <- NA
      }else{
        if((var[l]) == 0 ){var[l] <- 0
        }else{
          if(var[l] > max){var[l] <- max}
        }
      }
    }
    #return the result    
    res <- t(data.frame(var))
    colnames(res) <- paste0("By_", 1:length(.Overlap))
    rownames(res) <- i
    
    ifelse( ( exists("result") ), result <- rbind( result, res ), result <- res) 
    
  }
  
  
  Lap <- cbind(lap, result)
  
  return(Lap)
} 

piCatcherNeg <- function(x, .Length_In = parent.frame()$Length_In, 
                         .Target_Length = parent.frame()$Target_Length, 
                         .Overlap = parent.frame()$Overlap, 
                         Split = "Neg"){
  
  input <- x
  lap <-c()
  #split the data into the opposing strands
  datIn <- data.frame(input[, grepl(paste0(Split), colnames(input))])
  datIn <- data.frame(datIn[, names(dplyr::select(datIn, matches(c( as.character(.Length_In)))))])#makes sure we only keep the lengths we want
  
  #make two frames from each strand with the total reads starting at each position
  datIn$Total <- rowSums(datIn)
  datIn$pos <- row.names(datIn)
  lap <- datIn[, c( (ncol(datIn)), (ncol(datIn)-1) ), drop = F]
  row.names(lap) <- lap$pos
  
  #split the data into the opposing strands
  datO <- data.frame(input[, !grepl(paste0(Split), colnames(input))])
  datO <- data.frame(datO[,names(dplyr::select(datO, matches(c( as.character(.Target_Length)))))])#only keep the targets
  
  #make two frames from each strand with the total reads starting at each position
  datO$Total <- rowSums(datO)
  datO$pos <- row.names(datO)
  datO <- datO[, c( (ncol(datO)), (ncol(datO)-1) ), drop = F]
  
  result <- c()
  
  for (i in 1:nrow(lap)){
    
    max <- lap[i,2]
    
    #find the range of rows which are above 0 
    range <- ( ( paste(i-length(.Overlap)) ) : ( paste(i) ) )
    ind <- which( ( paste(i-length(.Overlap)) ) : ( paste(i) ) > 0)
    #set the target indexes     
    targ <- range[ind]
    #get the target    
    target <- datO[c( targ ) , "Total"]
    var <- c()
    #count the overlaps    
    
    for(l in 1:length(.Overlap)){
      
      var[l] <- target[(l)]
      #this now limits the overlaps counted to the number of input sequences 
      
      if(is.na(var[l])){var[l] <- NA
      }else{
        if((var[l]) == 0 ){var[l] <- 0
        }else{
          if(var[l] > max){var[l] <- max}
        }
      }
    }
    #return the result    
    res <- t(data.frame(var))
    colnames(res) <- paste0("By_", 1:length(.Overlap))
    rownames(res) <- i
    
    ifelse( ( exists("result") ), result <- rbind( result, res ), result <- res) 
    
  }
  
  
  Lap <- cbind(lap, result)
  
  return(Lap)
} 

piCatcherDual <- function(x, Length_In , Target_Length , Overlap = c(1:21)){
  #define and set inputs
  input <- x
  datPos <- c()
  datNeg <- c()
  
  datPos <- piCatcherPos(input)
  datNeg <- piCatcherNeg(input)
  oLaps <- overLaps(datPos, Overlap = parent.frame()$Overlap)
  oLapsO <- overLaps(datNeg, Overlap = parent.frame()$Overlap)
  
  oLapsFinal <- cbind(oLaps,oLapsO)
  oLapsFinal <- oLapsFinal[,c(1,2,4)]
  colnames(oLapsFinal) <- c("x", "Pos_count", "Neg_count")
  
  #get the Z-scores
  
  oLapsFinal$Pos_Z <- ((oLapsFinal$Pos_count - 
                          mean(oLapsFinal$Pos_count))/
                         sd(oLapsFinal$Pos_count))
  
  oLapsFinal$Neg_Z <- ((oLapsFinal$Neg_count - 
                          mean(oLapsFinal$Neg_count))/
                         sd(oLapsFinal$Neg_count))
  
  #get overlap probability for each position
  piLapProbPos <- piProb(datPos)
  piLapProbNeg <- piProb(datNeg)
  #get weighted probability for each position
  piLapProbWeightedPos <- piProbWeighted(datPos)
  piLapProbWeightedNeg <- piProbWeighted(datNeg)
  #calculate probability for the overlap length
  piProbTotPos <- piProbSum(piLapProbPos)
  piProbTotNeg <- piProbSum(piLapProbNeg)
  
  piLapProb <- cbind(piProbTotPos,piProbTotNeg)
  piLapProb <- piLapProb[,c(2,4)]
  colnames(piLapProb) <- c("Pos_probability", "Neg_probability")
  #calculate weighted probability for the overlap length
  piProbWeightedTotPos <- piProbSum(piLapProbWeightedPos)
  piProbWeightedTotNeg <- piProbSum(piLapProbWeightedNeg)
  
  piProbWeight <- cbind(piProbWeightedTotPos, piProbWeightedTotNeg)
  piProbWeight <- piProbWeight[,c(2,4)]
  colnames(piProbWeight) <- c("Pos_weighted_probability", "Neg_weighted_probability")
  
  #get results
  ResFin <- cbind(oLapsFinal, piLapProb, piProbWeight)
  
  #get totals
  ResFin$Overlaps <- ResFin$Pos_count + ResFin$Neg_count
  ResFin$Z_Score <- ((ResFin$Overlaps - 
                        mean(ResFin$Overlaps))/
                       sd(ResFin$Overlaps))
  ResFin$Probability <- (ResFin$Pos_probability+ResFin$Neg_probability)/sum(ResFin$Pos_probability,ResFin$Neg_probability)
  
  ResFin$Weighted_Probability <- (ResFin$Pos_weighted_probability+ResFin$Neg_weighted_probability)/sum(ResFin$Pos_weighted_probability,ResFin$Neg_weighted_probability)
  
  return(ResFin)
  
}

makeTally <- function(x, Genome_Length, Target_Length){
  #define the parameters
  data <- x
  Count <- Target_Length
  GenPosn <- Genome_Length
  
  for(t in 1:length(GenPosn)){
    targ <- as.numeric(GenPosn[t])#sets a target for a single nucleotide in the genome t=25
    pat <- (dplyr::filter(datS, datS$pos == targ))#returns all reads at that target position
    if(is.na(pat[1,1])){
      pat[1,] <-  c(paste0(rSeq), paste0(spt), (paste0(GenPosn[t])),(0L),(0L))
    }#sets a grid of 0 if there are no reads at the target position
    pat <- pat[,3:4]#removes excess input
    res <- countMatrix(pat, Count)#returns sum of the count of read length at the target
    res <- data.frame(t(res))
    row.names(res) <- targ
    res$pos <- as.numeric(row.names(res))
    res <- res[, c( ncol(res), (1:(ncol(res)-1))  )]
    ifelse(exists("resultMakeTally"), resultMakeTally <- rbind(resultMakeTally, res), resultMakeTally <- res)
  }
  
  return(resultMakeTally)
}

CoverageSummaryPlot <- function(GenSummary_Plot,  siRNA_Summary, piRNA_Summary){
  
  filename <- "Final_Coverage_Summary_Plot"
  
  namTopLeft <- names(GenSummary_Plot[1])
  datTopLeft <- GenSummary_Plot[[1]]
  
  namTopRight <- names(GenSummary_Plot[2])
  datTopRight <- GenSummary_Plot[[2]]
  
  namMidLeft <- names(siRNA_Summary[1])
  namMidLeft <- str_extract(namMidLeft, '(?<=_)[A-Z0-9]*')
  datMidLeft <- siRNA_Summary[[1]]
  
  namMidRight <- names(siRNA_Summary[2])
  namMidRight <- str_extract(namMidRight, '(?<=_)[A-Z0-9]*')
  datMidRight <- siRNA_Summary[[2]]
  
  namBotLeft <- names(piRNA_Summary[1])
  namBotLeft <- str_extract(namBotLeft, '(?<=_)[A-Z0-9]*')
  datBotLeft <- piRNA_Summary[[1]]
  
  namBotRight <- names(piRNA_Summary[2])
  namBotRight <- str_extract(namBotRight, '(?<=_)[A-Z0-9]*')
  datBotRight <- piRNA_Summary[[2]]
  
  TopLeft <- ggplot(data = datTopLeft)+
    geom_col( aes(x = x, y = Pos_Mean), fill = group.colours["Pos"], linewidth=0.5, colour="black", alpha=0.9)+
    geom_col( aes(x = x, y = Neg_Mean), fill = group.colours["Neg"],linewidth=0.5, colour="black", alpha=0.9)+
    geom_errorbar( aes(x = x, ymin = Pos_E_Min, ymax = Pos_E_Max), linewidth=1, colour="black", alpha=0.9 )+
    geom_errorbar( aes(x = x, ymin = Neg_E_Min, ymax = Neg_E_Max), linewidth=1, colour="black", alpha=0.9 )+
    xlab ("Size")+
    ylab ("Count")+
    ylim (-(max(GenMaxCount)*1.1),(max(GenMaxCount)*1.1))+
    ggtitle(namTopLeft)+
    piMaker_theme+
    theme(legend.position = "none")
  
  TopRight <- ggplot(data = datTopLeft)+
    geom_col( aes(x = x, y = Pos_Mean), fill = group.colours["Pos"], linewidth=0.5, colour="black", alpha=0.9)+
    geom_col( aes(x = x, y = Neg_Mean), fill = group.colours["Neg"],linewidth=0.5, colour="black", alpha=0.9)+
    geom_errorbar( aes(x = x, ymin = Pos_E_Min, ymax = Pos_E_Max), linewidth=1, colour="black", alpha=0.9 )+
    geom_errorbar( aes(x = x, ymin = Neg_E_Min, ymax = Neg_E_Max), linewidth=1, colour="black", alpha=0.9 )+
    xlab ("Size")+
    ylab ("Count")+
    ylim (-(max(GenMaxCount)*1.1),(max(GenMaxCount)*1.1))+
    ggtitle(namTopRight)+
    piMaker_theme+
    theme(legend.position = "none")
  
  MidLeft <-  ggplot(data = datMidLeft)+
    geom_ribbon(aes(x = x, ymin = Pos_E_min, ymax = Pos_E_max), colour = group.colours["Pos"],
                fill = group.colours["Pos"], alpha = 0.2, linewidth = 0.1)+
    geom_line(aes(x = x, y = Pos_Mean),colour = group.colours["Pos"], linewidth = 0.5)+
    geom_ribbon(aes(x = x, ymin = Neg_E_min, ymax = Neg_E_max), colour = group.colours["Neg"],
                fill = group.colours["Neg"], alpha = 0.2, linewidth = 0.1)+
    geom_line(aes(x = x, y = Neg_Mean), colour = group.colours["Neg"], linewidth = 0.5)+
    geom_hline(yintercept = 0, linetype = "solid", colour = "grey")+
    ggtitle(paste0(namMidLeft, "_siRNA_Coverage"))+
    ylim(-max(NormScale), max(NormScale))+
    xlab ("nt position")+
    guides(guide_legend, fill = NULL)+
    piMaker_theme+
    theme(legend.position = "none")
  
  MidRight <- ggplot(data = datMidRight)+
    geom_ribbon(aes(x = x, ymin = Pos_E_min, ymax = Pos_E_max), colour = group.colours["Pos"],
                fill = group.colours["Pos"], alpha = 0.2, linewidth = 0.1)+
    geom_line(aes(x = x, y = Pos_Mean),colour = group.colours["Pos"], linewidth = 0.5)+
    geom_ribbon(aes(x = x, ymin = Neg_E_min, ymax = Neg_E_max), colour = group.colours["Neg"],
                fill = group.colours["Neg"], alpha = 0.2, linewidth = 0.1)+
    geom_line(aes(x = x, y = Neg_Mean), colour = group.colours["Neg"], linewidth = 0.5)+
    geom_hline(yintercept = 0, linetype = "solid", colour = "grey")+
    ggtitle(paste0(namMidRight, "_siRNA_Coverage"))+
    ylim(-max(NormScale), max(NormScale))+
    xlab ("nt position")+
    guides(guide_legend, fill = NULL)+
    piMaker_theme+
    theme(legend.position = "none")
  
  BotLeft <- ggplot(data = datBotLeft)+
    geom_ribbon(aes(x = x, ymin = Pos_E_min, ymax = Pos_E_max), colour = group.colours["Pos"],
                fill = group.colours["Pos"], alpha = 0.2, linewidth = 0.1)+
    geom_line(aes(x = x, y = Pos_Mean),colour = group.colours["Pos"], linewidth = 0.5)+
    geom_ribbon(aes(x = x, ymin = Neg_E_min, ymax = Neg_E_max), colour = group.colours["Neg"],
                fill = group.colours["Neg"], alpha = 0.2, linewidth = 0.1)+
    geom_line(aes(x = x, y = Neg_Mean), colour = group.colours["Neg"], linewidth = 0.5)+
    geom_hline(yintercept = 0, linetype = "solid", colour = "grey")+
    ggtitle(paste0(namBotLeft, "_piRNA_Coverage"))+
    ylim(-max(NormScale), max(NormScale))+
    xlab ("nt position")+
    guides(guide_legend, fill = NULL)+
    piMaker_theme+
    theme(legend.position = "none")
  
  BotRight <- ggplot(data = datBotRight)+
    geom_ribbon(aes(x = x, ymin = Pos_E_min, ymax = Pos_E_max), colour = group.colours["Pos"],
                fill = group.colours["Pos"], alpha = 0.2, linewidth = 0.1)+
    geom_line(aes(x = x, y = Pos_Mean),colour = group.colours["Pos"], linewidth = 0.5)+
    geom_ribbon(aes(x = x, ymin = Neg_E_min, ymax = Neg_E_max), colour = group.colours["Neg"],
                fill = group.colours["Neg"], alpha = 0.2, linewidth = 0.1)+
    geom_line(aes(x = x, y = Neg_Mean), colour = group.colours["Neg"], linewidth = 0.5)+
    geom_hline(yintercept = 0, linetype = "solid", colour = "grey")+
    ggtitle(paste0(namBotRight, "_piRNA_Coverage"))+
    ylim(-max(NormScale), max(NormScale))+
    xlab ("nt position")+
    guides(guide_legend, fill = NULL)+
    piMaker_theme+
    theme(legend.position = "none")
  
  CoverageSummary <- ggarrange(TopLeft ,
                               TopRight +rremove("ylab") +rremove("y.text"),
                               MidLeft  +rremove("xlab") +rremove("x.text"),
                               MidRight +rremove("xlab") +rremove("x.text") +rremove("ylab") +rremove("y.text"),
                               BotLeft,
                               BotRight +rremove("ylab") +rremove("y.text"),
                               nrow = 3, ncol = 2, widths = 1, heights = 1, align = "hv" )
  plot(CoverageSummary)
  saveImage(paste(filename))
  
}

piMapper <- function(x, Scale){
  
  dat_2429. <- x
  
  Scale. = Scale
  
  for(n in 24:29){
    dat_2429.[paste0("Pos_",n, "_x")] <- as.numeric(ifelse( (dat_2429.[(paste0("Pos_",n))] >0), 
                                                            dat_2429.$Pos_pos, NA ))
    dat_2429.[paste0("Pos_",n, "_xend")] <- dat_2429.[paste0("Pos_",n, "_x")] + n
    dat_2429.[paste0("Neg_",n, "_x")] <- as.numeric(ifelse((dat_2429.[(paste0("Neg_",n))] >0), 
                                                           dat_2429.$Neg_pos, NA))
    dat_2429.[paste0("Neg_",n, "_xend")] <- dat_2429.[paste0("Neg_",n, "_x")] + n
  }
  
  piMap <- ggplot(data = dat_2429.)+
    
    geom_hline(yintercept = 0, linetype = "solid", colour = "grey")+
    
    geom_rect(data = dat_2429., aes(ymin = 0, ymax = Pos_24, xmin = Pos_24_x,
                                    xmax = Pos_24_xend, fill = "24"), colour = piRNA.colours["24"], alpha = 0.4)+
    geom_rect(data = dat_2429., aes(ymin = 0, ymax = - Neg_24, xmin = Neg_24_x,
                                    xmax = Neg_24_xend, fill = "24"), colour = piRNA.colours["24"], alpha = 0.4)+ 
    
    geom_rect(data = dat_2429., aes(ymin = 0, ymax = Pos_25, xmin = Pos_25_x,
                                    xmax = Pos_25_xend, fill = "25"), colour = piRNA.colours["25"],  alpha = 0.4)+
    geom_rect(data = dat_2429., aes(ymin = 0, ymax = - Neg_25, xmin = Neg_25_x,
                                    xmax = Neg_25_xend, fill = "25"), colour = piRNA.colours["25"],  alpha = 0.4)+
    
    geom_rect(data = dat_2429., aes(ymin = 0, ymax = Pos_26, xmin = Pos_26_x,
                                    xmax = Pos_26_xend, fill = "26"), colour = piRNA.colours["26"],  alpha = 0.4)+
    geom_rect(data = dat_2429., aes(ymin = 0, ymax = - Neg_26, xmin = Neg_26_x,
                                    xmax = Neg_26_xend, fill = "26"), colour = piRNA.colours["26"],  alpha = 0.4)+
    
    geom_rect(data = dat_2429., aes(ymin = 0, ymax = Pos_27, xmin = Pos_27_x,
                                    xmax = Pos_27_xend, fill = "27"), colour = piRNA.colours["27"],  alpha = 0.4)+
    geom_rect(data = dat_2429., aes(ymin = 0, ymax = - Neg_27, xmin = Neg_27_x,
                                    xmax = Neg_27_xend, fill = "27"), colour = piRNA.colours["27"],  alpha = 0.4)+
    
    geom_rect(data = dat_2429., aes(ymin = 0, ymax = Pos_28, xmin = Pos_28_x,
                                    xmax = Pos_28_xend, fill = "28"), colour = piRNA.colours["28"],  alpha = 0.4)+
    geom_rect(data = dat_2429., aes(ymin = 0, ymax = - Neg_28, xmin = Neg_28_x,
                                    xmax = Neg_28_xend, fill = "28"), colour = piRNA.colours["28"],  alpha = 0.4)+
    
    geom_rect(data = dat_2429., aes(ymin = 0, ymax = Pos_29, xmin = Pos_29_x,
                                    xmax = Pos_25_xend, fill = "29"), colour = piRNA.colours["29"],  alpha = 0.4)+
    geom_rect(data = dat_2429., aes(ymin = 0, ymax = - Neg_29, xmin = Neg_29_x,
                                    xmax = Neg_25_xend, fill = "29"), colour = piRNA.colours["29"],  alpha = 0.4)+
    ylim(-max(Scale),max(Scale))+
    ggtitle(paste(nam))+
    piMaker_theme+
    scale_fill_manual("Size", values = c(piRNA.colours))
  
}
