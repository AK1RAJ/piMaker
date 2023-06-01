#Themes and Libraries####
# Package names
packages <- c("BiocManager", "Biostrings", "tidyverse", "Rsamtools", "ggseqlogo", "ggpubr" )

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
# Packages loading
if (!require("BiocManager", quietly = TRUE)){
  invisible(lapply(packages, library, character.only = TRUE))
}

source("https://raw.githubusercontent.com/AK1RAJ/piMaker/main/piMaker_functions.R")

#set working directory and get the files####
#make a project folder with three subfolders for 1- the BAM files (BAM), 2- the reference sequences (refSeq)
#3- the output (Output)
DIR <- "F:/TidyCode"
BAM <- paste0(DIR,"/BAM")
REF <- paste0(DIR,"/refSeq")
OUT <- paste0(DIR,"/Output")
setwd(paste(REF))
#get the genome length, if this is not done, the length will be calculated from 
#the BAM file, but may miss uncovered regions
setwd(paste(REF))
refSeq <- sub('\\..[^\\.]*$', '',(c(list.files(full.names = FALSE, include.dirs = FALSE, no.. = FALSE ))))
if(exists("GenLength")){rm(GenLength)}
for (i in 1:length(refSeq)) {
  filename <- refSeq[i]
  dat <- nchar(str_flatten(readDNAStringSet(file = paste0(filename,".fa"))))
  gm <- cbind(filename,dat)
  ifelse(exists("GenLength"), (GenLength <- rbind(GenLength, gm)), (GenLength <- gm))
  GenLength <- data.frame(GenLength)
  if(i == length(refSeq)){rm(dat, gm)}
}
#make the file list
setwd(BAM)
files <- BamFileList((c(list.files(full.names = FALSE, #list of the BAM files
                                   include.dirs = TRUE, no.. = FALSE ))))
#assign the samples here
samples <- c("MOCK", "rTOSV")
#read in the files
if(exists("BAMList")){rm(BAMList)}
for (i in 1:length(files)){
  namefile <- names(files[i])
  namefile <- tidyName(namefile)
  data <- filesIn(files[i])
  MaxCount <- maxCount(data, as.data.frame(sum(data$strand == "+"), row.names = namefile),
                     as.data.frame(sum(data$strand == "-"), row.names = namefile))
  readLength <- getGenLength(data)
  if (exists("BAMList")){ 
      BAMList[[namefile]] <- data
    }else{  
      BAMList <-list()
      BAMList[[namefile]] <- data
    }
  rm(data)
}
#make size distributions rm(Count)
if(exists("SizeDistribution")){rm(SizeDistribution)}
for (i in 1:length(BAMList)){
  filename <- names(BAMList[i])
  data <- (BAMList[[i]])
  nam <- paste0(filename, "_Size_Distribution")
  dsz <- getSize(data, Count = c(18:40))
      if(exists("SizeDistribution")){
        SizeDistribution[[nam]] <- dsz
      }else{
        SizeDistribution <- list()
        SizeDistribution[[nam]] <- dsz
        }
  rm(data, dsz) 
}
#plot size distributions
for(i in 1:length(SizeDistribution)){
  nam <- names(SizeDistribution[i])
  data <- SizeDistribution[[i]]
  dat <- data.frame(t(data))
  dat$x <- row.names(dat)
  gg <- ggplot(dat, aes(x = x) )+
    geom_col(aes( y = Pos))+
    geom_col(aes( y = Neg))+
    xlab ("Size")+
    ylab ("Count")+
    ylim((0-max(MaxCount)),(0+max(MaxCount)))+
    ggtitle(paste0(nam))
  plot(gg)
rm(nam,dat)
}
#make mean size distributions
if(exists("Summary_Plot")){rm(Summary_Plot)}
if(length(samples>1)){
for(i in 1:length(samples)){
  namefile <- samples[i]
  filename <- paste0(namefile, "_SizeDistSum")
  data <- (SizeDistribution[ grepl((paste0(namefile)), names(SizeDistribution))])
  for (i in 1:length(data)) {
    nam <- names(data[i])
    dat <- data.frame(t(data[[i]]))
    colnames(dat) <- c(paste0(nam, "_Pos"), paste0(nam, "_Neg"))
    ifelse(exists("res"),
           res <- bind_cols(res, dat),
           res <- dat)
  }
  res <- SizDistSum(res)
  if(exists("Summary_Plot")){
    Summary_Plot[[filename]] <- res
    rm(res)
  }else{
    Summary_Plot <- list()
    Summary_Plot[[filename]] <- res
    rm(res)
  }
}
}
#plot the means
if(exists("Summary_Plot")){
for(i in 1:length(Summary_Plot)){
  data <- Summary_Plot[[i]]
  namefile <- names(Summary_Plot[i])
gg <- ggplot(data = data)+
  geom_col( aes(x = data$x, y = data$Pos_Mean, fill = "black"), linewidth=0.5, colour="black", alpha=0.9)+
  geom_col( aes(x = data$x, y = data$Neg_Mean, fill = "red"), linewidth=0.5, colour="black", alpha=0.9)+
  geom_errorbar( aes(x = data$x, ymin = data$Pos_E_Min, ymax = data$Pos_E_Max), linewidth=1, colour="black", alpha=0.9 )+
  geom_errorbar( aes(x = data$x, ymin = data$Neg_E_Min, ymax = data$Neg_E_Max), linewidth=1, colour="black", alpha=0.9 )+
  xlab ("Size")+
  ylab ("Count")+
  ylim (-(max(MaxCount)*1.1),(max(MaxCount)*1.1))+
  theme(legend.position = "none")+
  ggtitle(paste0(namefile))
plot(gg)
}
}
#split the data in the 21nt siRNAs 
if(exists("siRNA_Coverage_Plot")){rm(siRNA_Coverage_Plot, CovCount)}
for (i in 1:length(BAMList)){
  filename <- names(BAMList[i])
  data <- (BAMList[[i]])
  d21 <- filter(data, data$qwidth == 21)
  d21 <- split(d21, d21$rname)
  for (n in 1:length(d21)){
    rSeq <- names(d21[n])
    rsq <- makeRsq(rSeq)
    namb <- paste0(names(BAMList[i]), "_", names(d21[n]))
    dat <- d21[[n]]
    datc <- coverMatrix(dat, Length = 21)
    cov <- getCoverage(datc, Count = 1:(paste0(rsq)))
    CovCount <- maxCount(cov, max(cov$Pos), abs(min(cov$Neg)))
    if(exists("siRNA_Coverage_Plot")){
      siRNA_Coverage_Plot[[namb]] <- cov
      rm(cov)
    }else{
      siRNA_Coverage_Plot <- list()
      siRNA_Coverage_Plot[[namb]] <- cov
      rm(cov)
    }
  }
  rm(data,d21,dat,namb,datc)
}
#plot the coverage
for (i in 1:(length(siRNA_Coverage_Plot))) {
  filename <- names(siRNA_Coverage_Plot[i])
  cvr <- data.frame(siRNA_Coverage_Plot[[i]])
  gg<- ggplot()+
    geom_line(data = cvr, aes(x= x, y = Pos ), colour = "black", linewidth = 0.5)+
    geom_line(data = cvr, aes(x= x, y = Neg), colour = "red", linewidth = 0.5)+
    ggtitle(paste0(filename))+
    ylim(-(max(CovCount)*1.1),(max(CovCount)*1.1))+
    xlab ("nt position")
  plot(gg)
  rm(cvr)
}
#normalise the data to take unequal read counts in account
for (i in 1:(length(siRNA_Coverage_Plot))){
  filename = names(siRNA_Coverage_Plot[i])
  data <- (siRNA_Coverage_Plot[[i]])
  data$Pos_Normalised <- Normalise(data$Pos)
  data$Neg_Normalised <- Normalise(data$Neg)
  if (exists("siRNA_Coverage_Plot_Normalised")){
    siRNA_Coverage_Plot_Normalised[[filename]] <- data
  } else {
    siRNA_Coverage_Plot_Normalised <- list()
    siRNA_Coverage_Plot_Normalised[[filename]] <- data
  }
  if(i == length(siRNA_Coverage_Plot)){rm(data)}
}
#plot the normalised data
for (i in 1:(length(siRNA_Coverage_Plot_Normalised))) { 
  filename <- names(siRNA_Coverage_Plot_Normalised[i])
  namefile <- (paste0(filename, "_Normalised"))
  data <- siRNA_Coverage_Plot_Normalised[[i]]
  gg<- ggplot(data)+
    geom_line(aes(x= x, y = Pos_Normalised), colour = "black", linewidth = 0.5)+
    geom_line(aes(x= x, y = Neg_Normalised), colour = "red", linewidth = 0.5)+
    ggtitle(paste0(namefile))+
    ylim(-1,1)+
    xlab ("nt position")
  plot(gg)
}
#make summary data
if(exists("siRNA_Summary")){rm(siRNA_Summary)}
for (i in 1:length(samples)) {
  namefile = samples[i]
  data <- siRNA_Coverage_Plot_Normalised[grepl(namefile, names(siRNA_Coverage_Plot_Normalised))]
  for (j in 1:length(refSeq)) {
    rSeq <- refSeq[j]
    datRseq <- data[grepl(rSeq, names(data))] %>% keep( ~ !is.null(.) ) 
      for (n in 1:length(datRseq)) {
        namen <- names(datRseq[n])
        datn <- as.data.frame(datRseq [[n]] [,c("Pos_Normalised","Neg_Normalised")] ) 
        colnames(datn) <- paste0(namen, "_", colnames(datn))
        ifelse(exists("Sumry"), 
               Sumry <- cbind(Sumry, datn),
               Sumry <- (datn))
      }
    filename <- paste0(namefile, "_", rSeq)
    Sumry <- meanCalc(Sumry)
    Sumry <- Sumry[,grepl(("Mean|SD|min|max"), colnames(Sumry))]
    Sumry$x <- as.numeric(row.names(Sumry))
    if(exists("siRNA_Summary")){
      siRNA_Summary[[filename]] <- Sumry
      rm(Sumry)
    }else{
      siRNA_Summary <- list()
      siRNA_Summary[[filename]] <- Sumry
      rm(Sumry)
    }
  }
  if(i == length(samples)){rm(data,datRseq,datn)}
}  
#plot summary data        
for(i in 1:length(siRNA_Summary)){
  filename <- names(siRNA_Summary[i])
  data <- siRNA_Summary[[i]]
  gg<- ggplot(data)+
    geom_ribbon(aes(x = x, ymin = Pos_E_min, ymax = Pos_E_max, colour = "SDpos"), alpha = 0.5, linewidth = NULL)+
    geom_line(aes(x = x, y = Pos_Mean, colour = "Meanpos"), linewidth = 1)+
    geom_ribbon(aes(x = x, ymin = Neg_E_min, ymax = Neg_E_max, colour = "SDneg"), alpha = 0.5, linewidth = NULL)+
    geom_line(aes(x = x, y = Neg_Mean, colour = "Meanneg"), linewidth = 1)+
    scale_colour_manual("", 
                        breaks = c("SDpos", "Meanpos", "SDneg", "Meanneg"),
                        values = c("black", "black", "red", "red"))+
    ggtitle(paste0(filename))+
    ylim(-1,1)+
    xlab ("nt position")+
    guides(guide_legend, fill = NULL)+
    theme(legend.position = "none")
  plot(gg)
}        

#piRNAs####
#load functions and colours
#function to find overlaps in piRNAs piFlinger and piCatcher are required
#these colours match standard peak viewers
DNA_col_scheme =  make_col_scheme(chars=c('A', 'T', 'G', 'C'), groups= NULL, 
                                  cols=c('green3', 'red2', 'black', 'blue'), name='DNA_col_scheme')

Sizes <- c(24:29) #put in the piRNA sizes we are looking for here
#21nt siRNAs are dealt with separately for overlap and z score calculations
setwd(paste(OUT))
#make long frames of the sequences, first filter the data sets by size
if(exists("piList")){rm(piList)}
if(exists("siList")){rm(siList)}
for (i in 1:length(samples)) {
  namefile = samples[i]
  data <- BAMList[grepl(namefile, names(BAMList))]
  for (l in 1:(length(data))) {
    filename <- names(data[l])
    datl <- data[[l]]
    datl <- datl[,c(3,4:6,12)]
    datp <- filter(datl, between( datl$qwidth,  as.numeric(min(Sizes)),  as.numeric(max(Sizes)) ) )
    dats <- filter(datl,  datl$qwidth == 21)
    nam <- (paste0(filename, "_piSequences"))
    if(exists("piList")){
      piList[[nam]] <- datp
    }else{
      piList <- list()
      piList[[nam]] <- dats
    }
    if(exists("siList")){
      siList[[nam]] <- dats
    }else{
      siList <- list()
      siList[[nam]] <- dats
    }
  }   
  rm(data, datl, datp, dats)
}
#now combine into single sets for each sample
if(exists("piSeqList")){rm(piSeqList)}
for (i in 1:length(samples)) {
  namefile = samples[i]
  filename <- paste0(namefile, "_piRNAs")
  data <- piList[grepl(namefile, names(piList))]
  for (l in 1:length(data)){
    dat <- data[[l]]
    ifelse(exists("pi_res"), pi_res <- rbind(pi_res,dat), pi_res <- dat)
  }
    if(exists("piSeqList")){
      piSeqList[[filename]] <- pi_res
      rm(pi_res)
    }else{
      piSeqList <- list()
      piSeqList[[filename]] <- pi_res
      rm(pi_res)
    }
  }
#make piRNA position matrix
if(exists("piMatrix")){rm(piMatrix)}
for (i in 1:(length(piSeqList))) { 
  namfile <- names(piSeqList[i])
  data <- piSeqList[[i]]
  for (r in 1:(length(refSeq))) {
    rSeq <- refSeq[r]
    nam <- paste0(namfile, "_", rSeq)
    rsq <- makeRsq(rSeq)
    datr <- filter(data, data$rname == rSeq)
    freq <- piFreq(datr)
    if(exists("piMatrix")){
      piMatrix[[nam]] <- freq
    }else{
      piMatrix <- list()
      piMatrix[[nam]] <- freq
    }    
  }
  if(i == length(samples)){rm(data,datr,freq)}
}
#plot piRNA position matrix
for (i in 1: length(piMatrix)){
  data <- piMatrix[[i]]
  filename <- names(piMatrix[i])
  rsq <- str_extract(filename,
                     '([A-Za-z].)([0-9]+)')
  rsq <- makeRsq(rsq)
gg <- ggplot(data)+
  geom_hline(yintercept = 0, linetype = "solid", colour = "grey")+
  geom_col(aes(x = pos, y = Pos, colour = "black"))+
  geom_col(aes(x = pos, y = Neg, colour = "red"))+
  ggtitle(paste0(filename))+
  theme(legend.position = "none")+
  ylim(-100,100)+
  xlim(0,rsq)
plot(gg)
}
#now split the lists into respective piRNAs
if(exists("piRNA_List")){rm(piRNA_List)}
for (i in 1:(length(piSeqList))) {
  namfile <- names(piSeqList[i])
  data <- piSeqList[[i]]
  for (s in 1:length(Sizes)) {
    sz = Sizes[s]
    datsz <- dplyr::filter(data, data$qwidth == sz)
    for (j in 1:length(refSeq)) {
      rSeq = refSeq[j]
      namj = paste0(namfile, "_", rSeq, "_", sz)
      datj <- dplyr::filter(datsz, datsz$rname == rSeq)
      namP <- paste0(namj,"_Pos")
      namN <- paste0(namj, "_Neg")
      datjPos <- dplyr::filter(datj, datj$strand == "+")
      datjNeg <- dplyr::filter(datj, datj$strand == "-")
      datjNeg$RC <- reverseComplement(DNAStringSet(datjNeg$seq))
      piPos <- as.data.frame(substr(datjPos$seq, 1, 20))
      piNeg <- as.data.frame(substr(datjNeg$RC, 1, 20))
      if(exists("piRNA_List")){
        piRNA_List[[namP]] <- piPos
        piRNA_List[[namN]] <- piNeg
      } else{
        piRNA_List <- list()
        piRNA_List[[namP]] <- piPos
        piRNA_List[[namN]] <- piNeg
      }
      #and make the sequence plots
      gp <- ggplot()+
        geom_logo(piPos, method = 'prob', seq_type = "dna", font = "helvetica_bold", col_scheme = DNA_col_scheme)+
        theme_logo()+
        ggtitle(namP)
      gn <- ggplot()+
        geom_logo(piNeg, method = 'prob', seq_type = "dna", font = "helvetica_bold", col_scheme = DNA_col_scheme)+
        theme_logo()+
        ggtitle(namN)
      fig <- ggarrange(gp,gn, nrow = 2)
      plot(fig)
    }
  }
  if(i == length(samples)){rm(data, datj, datjNeg, datjPos, datsz, piPos, piNeg)}
} 
#make the position matrix for the reads, counts how many of each read length at each position
#this can take a while! rm(Count) s=1
if(exists("Tally_List")){rm(Tally_List)}
for (i in 1:(length(piSeqList))) {
  namfile <- names(piSeqList[i])
  data <- piSeqList[[i]]
  for (r in 1:(length(refSeq))) {
    rSeq <- refSeq[r]
    datr <- dplyr::filter(data, data$rname == rSeq)
    datr <- split(datr, datr$strand)
    datr <- datr[-3]
    GenPosn <- c(1:makeRsq(rSeq))
    for (s in 1:length(datr)) {
      spt <- names(datr[s])
      ifelse(names(datr[s]) == "+", namC <- "Pos", namC <- "Neg")
      namB <- paste0(namfile,"_", rSeq,"_", namC, "_Tally")
      datS <- datr[[s]]
      if(names(datr[s]) == "-"){
        datS$pos <- datS$pos + datS$qwidth#this makes sure the 5' of the negative strand is in the correct position
      }
        for(t in 1:length(GenPosn)){
          targ <- as.numeric(GenPosn[t])#sets a target for a single nucleotide in the genome t=25
          pat <- (dplyr::filter(datS, datS$pos == targ))#returns all reads at that target position
          if(is.na(pat[1,1])){
            pat[1,] <-  c(paste0(rSeq), paste0(spt), (paste0(GenPosn[t])),(0L),(0L))
          }#sets a grid of 0 if there are no reads at the target position
          pat <- pat[,3:4]#removes excess input
          res <- countMatrix(pat, Count = c(24:29))#returns sum of the count of read length at the target
          res <- data.frame(t(res))
          row.names(res) <- targ
          res$pos <- as.numeric(row.names(res))
          res <- res[, c( ncol(res), (1:(ncol(res)-1))  )]
          ifelse(exists("result"), result <- rbind(result, res), result <- res)
        }
      if(exists("Tally_List")){
        colnames(result) <- c("pos", c(24:29))
        Tally_List[[namB]] <- result
        rm(result,res)
      }else{
        Tally_List <- list()
        colnames(result) <- c("pos", c(24:29))
        Tally_List[[namB]] <- result
        rm(result,res)
      }
    }
  }
  rm(datS)
}
#join the positive and negative reads into a single matrix
if(exists("Tally_List_Matrix")){rm(Tally_List_Matrix)}
for (i in 1:length(samples)){
  smp <- samples[i]
  namS <- paste0(smp)
  for (r in 1:length(refSeq)){
    rSeq <- refSeq[r]
    namRP <- paste0(smp, "_piRNAs_", rSeq,  "_Pos_Tally")
    namRN <- paste0(smp, "_piRNAs_", rSeq, "_Neg_Tally")
    namRS <- paste0(smp, "_", rSeq, "_Tally_Matrix")
    datP <- Tally_List[[(paste0(namRP))]]
    colnames(datP) <- paste0("Pos_", colnames(datP))
    datN <- Tally_List[[(paste0(namRN))]]
    colnames(datN) <- paste0("Neg_", colnames(datN))
    datC <- cbind(datN, datP)
    if(exists("Tally_List_Matrix")){
      Tally_List_Matrix[[namRS]] <- datC
    } else{
      Tally_List_Matrix <- list()
      Tally_List_Matrix[[namRS]] <- datC
    }
  }
  rm(smp,rSeq,namRP,namRN,namRS,datP,datN,datC)
}
#map the piRNA overlaps
for(i in 1:length(Tally_List_Matrix)){
  namfile <- names(Tally_List_Matrix[i])
  dat <- Tally_List_Matrix[[i]]
  nam <- paste0(namfile, "_2429")
  dat_2429 <- dat
  
  for(n in 24:29){
    dat_2429[paste0("Pos_",n, "_x")] <- as.numeric(ifelse( (dat_2429[(paste0("Pos_",n))] >0), 
                                               dat_2429$Pos_pos, NA ))
    dat_2429[paste0("Pos_",n, "_xend")] <- dat_2429[paste0("Pos_",n, "_x")] + n
    dat_2429[paste0("Neg_",n, "_x")] <- as.numeric(ifelse((dat_2429[(paste0("Neg_",n))] >0), 
                                               dat_2429$Neg_pos, NA))
    dat_2429[paste0("Neg_",n, "_xend")] <- dat_2429[paste0("Neg_",n, "_x")] + n
  }
  
  gg <- ggplot()+
    geom_hline(yintercept = 0, linetype = "solid", colour = "grey")+
    geom_rect(data = dat_2429, aes(ymin = 0, ymax = Pos_24, xmin = Pos_24_x,
                                   xmax = Pos_24_xend, colour = "red", fill = "red", group = "24",  alpha = 0.4))+
    geom_rect(data = dat_2429, aes(ymin = 0, ymax = - Neg_24, xmin = Neg_24_x,
                                   xmax = Neg_24_xend, colour = "red", fill = "red",  group = "24", alpha = 0.4))+ 
    
    
    geom_rect(data = dat_2429, aes(ymin = 0, ymax = Pos_25, xmin = Pos_25_x,
                                   xmax = Pos_25_xend, colour = "black", fill = "black",  alpha = 0.4))+
    geom_rect(data = dat_2429, aes(ymin = 0, ymax = - Neg_25, xmin = Neg_25_x,
                                   xmax = Neg_25_xend, colour = "black", fill = "black",  alpha = 0.4))+
    
    
    geom_rect(data = dat_2429, aes(ymin = 0, ymax = Pos_26, xmin = Pos_26_x,
                                   xmax = Pos_26_xend, colour = "green", fill = "green",  alpha = 0.4))+
    geom_rect(data = dat_2429, aes(ymin = 0, ymax = - Neg_26, xmin = Neg_26_x,
                                   xmax = Neg_26_xend, colour = "green", fill = "green",  alpha = 0.4))+
    
    
    geom_rect(data = dat_2429, aes(ymin = 0, ymax = Pos_27, xmin = Pos_27_x,
                                   xmax = Pos_27_xend, colour = "orange", fill = "orange",  alpha = 0.4))+
    geom_rect(data = dat_2429, aes(ymin = 0, ymax = - Neg_27, xmin = Neg_27_x,
                                   xmax = Neg_27_xend, colour = "orange", fill = "orange",  alpha = 0.4))+
    
    
    geom_rect(data = dat_2429, aes(ymin = 0, ymax = Pos_28, xmin = Pos_28_x,
                                   xmax = Pos_28_xend, colour = "blue", fill = "blue",  alpha = 0.4))+
    geom_rect(data = dat_2429, aes(ymin = 0, ymax = - Neg_28, xmin = Neg_28_x,
                                   xmax = Neg_28_xend, colour = "blue", fill = "blue",  alpha = 0.4))+
    
    geom_rect(data = dat_2429, aes(ymin = 0, ymax = Pos_29, xmin = Pos_29_x,
                                   xmax = Pos_25_xend, colour = "purple", fill = "purple",  alpha = 0.4))+
    geom_rect(data = dat_2429, aes(ymin = 0, ymax = - Neg_29, xmin = Neg_29_x,
                                   xmax = Neg_25_xend, colour = "purple", fill = "purple",  alpha = 0.4))+
    ylim(-50,50)+
    guides( fill = FALSE, alpha = FALSE)+
    scale_colour_discrete(name = "Size",labels=c('24', '25', "26", "27", "28", "29"))+
    guides(color=guide_legend(override.aes=list(fill=NA)))+
    ggtitle(paste(nam))
  
  plot(gg)
}    
#make the overlap matrix
overlap <- c(1:21) #will count overlaps of this length
Split <- c("Pos", "Neg")#this is the two dimensions of the sequencing reads - doesn't need defining really!
if(exists("siLap_List")){rm(siLap_List)}+
if(exists("piLap_List")){rm(piLap_List)}+
if(exists("Z_List")){rm(Z_List)}
for(i in 1:length(Tally_List_Matrix)){
  namfile <- names(Tally_List_Matrix[i])
  namfile <- gsub("_Tally_Matrix","",namfile)
  dat <- Tally_List_Matrix[[i]]
  for(s in 1:length(c(Split))){
    spt <- Split[s]
    namC <- paste0(namfile, spt, "_siRNA_Overlap")
    namD <- paste0(namfile, spt, "_piRNA_Overlap")
    #get the overlaps
    siLap <- piCatcher(dat, Length_In = c(21), Target_Length = c(21), Overlap = c(1:21), Split = spt)
    piLap <- piCatcher(dat, Length_In = c(24:29), Target_Length = c(24:29), Overlap = c(1:21), Split = spt)
    #count the overlaps
    overSiLap <- overLaps(siLap, Overlap = c(1:21))
    overPiLap <- overLaps(piLap, Overlap = c(1:21))
      if(exists("siLap_List")){
        siLap_List[[namC]] <- overSiLap
        rm(siLap)
      }else{
        siLap_List <- list()
        siLap_List[[namC]] <- overSiLap
        rm(siLap)
      }
    if(exists("piLap_List")){
      piLap_List[[namD]] <- overPiLap
      rm(piLap)
    }else{
      piLap_List <- list()
      piLap_List[[namD]] <- overPiLap
      rm(piLap)
    }
#now calculate the z-score from the overlaps
siZscore <- Z_Score(overSiLap)
piZscore <- Z_Score(overPiLap)
    if(exists("Z_List")){
      Z_List[[namC]] <- siZscore
      Z_List[[namD]] <- piZscore
      rm(siZscore, piZscore)
    }else{
      Z_List <- list()
      Z_List[[namC]] <- siZscore
      Z_List[[namD]] <- piZscore
      rm(siZscore, piZscore)
    }
  }
} 
#plot the overlap and z-score
i=1
for (i in 1:length(samples)){
  namfile <- samples[i]
  for(r in 1:length(refSeq))
    Rnam <- refSeq[r]
  for(s in 1:length(Split))
    spt <- Split[s]
  dat <- piLap_List [[grep(paste0(namfile, "_", Rnam, spt, "_piRNA_Overlap"), names(piLap_List))]]
  namd <- paste0(namfile, "_", Rnam, "_", spt)
  datz <- Z_List [[grep(paste0(namfile, "_", Rnam, spt, "_piRNA_Overlap"), names(Z_List))]]

gg <- ggplot(data = datz, aes(x = x, y = Z))+
  geom_line(aes(colour = "red"), linewidth = 1 )+
    ggtitle(paste0(namd))+
  ylim(-5,5)+
  ylab("Z-Score")+
  xlab("Overlap (nt)")+
  theme(legend.position="none")+
  guides(colour = FALSE)+
  theme(aspect.ratio = 0.25:1)
plot(gg)

gg <- ggplot(data = dat, aes(x = x, y = count))+
  geom_col(aes(colour = "black"), linewidth = 1)+
  ggtitle(paste0(namd))+
  ylim(0,150)+
  xlab("Overlap (nt)")+
  ylab("No. of pairs")+
  guides(colour = FALSE)
plot(gg)

}
    
    
    







