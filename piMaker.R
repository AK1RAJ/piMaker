#Themes and Libraries####
# Package names
#this code uses the stringr, ggplot2, purrr and dplyr packages 
#from the tidyverse so it might be easier to install the whole tidyverse if wanted....
packages <- c( "BiocManager", "stringr" , "dplyr", "ggplot2", "purrr", "Biostrings",  "Rsamtools", "ggseqlogo", "ggpubr", "devtools" )

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
if (any(installed_packages ==FALSE)) {
  BiocManager::install(packages[!installed_packages])
}
# Packages loading
invisible(lapply(packages, library, character.only = TRUE))
# Packages might require forced installation on newer versions of R : BiocManager::install("ggseqlogo", force = T)
#devtools::install_github("omarwagih/ggseqlogo")
#load the functions from here
source("https://raw.githubusercontent.com/AK1RAJ/piMaker/main/piMaker_functions.R")

#test Bam files are at https://github.com/AK1RAJ/piMaker/tree/main/BAM
#test reference sequences are at: https://github.com/AK1RAJ/piMaker/tree/main/refSeq

savefiles = T #T/F option on whether to save the output plots or not

#set working directory and get the files####
#make a project folder with three subfolders for 1- the BAM files (BAM), 2- the reference sequences (refSeq)
#3- the output (Output)

DIR <- "C:/Users/Documents/TidyCode"                  #"C:/" put your directory here
BAM <- paste0(DIR,"/Mock_BAM") #location of the BAM alignment files
REF <- paste0(DIR,"/refSeq") #location of the reference sequences if you have them, if not do not run the first bit of code
OUT <- paste0(DIR,"/Output") #this directory is only used if you set the savefiles to true

#assign colours for the plots to differentiate positive reads, negative reads and combined:
group.colours <- c(Pos = "red3", Neg = "blue2", Overall = "skyblue4")
#assign colours to the piRNA sizes:
piRNA.colours <- c("24" = "magenta3", "25" = "dodgerblue2", "26" = "purple3",
                   "27" = "darkgoldenrod3", "28" = "aquamarine4", "29" = "orangered3")

#assigning names to the reference sequences
refNames <- as.matrix(c(c("RNA1_pp9ad"="RNA1"),
              c("RNA2_pp9ad"="RNA2")))
#putting in known protein coding regions here to plot on final graphs

KnownCoding <- list()
KnownCoding[["RNA1"]] <- data.frame(number = c(1,2,3),
                                     name = c("protein A",  "protein B1", "protein B2"),
                                     start = c(40, 2728, 2738),
                                     end = c(3036, 3036, 3058))
KnownCoding[["RNA2"]] <- data.frame(number = c(1),
                                     name = c("coat protein precursor"),
                                     start = c(22),
                                     end = c(1245))


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
  rm(dat, gm)
}

#make the file list
setwd(BAM)
files <- BamFileList((c(list.files(full.names = FALSE, #list of the BAM files
                                   include.dirs = TRUE, no.. = FALSE ))))
#assign the samples here
samples <- c("Mock")

#read in the files warning - do not have any variables in the environment called max, they will cause an error!
#in case you need it: rm(max)
if(exists("BAMList")){rm(BAMList, MaxCount)}
for (i in 1:length(files)){
  namefile <- names(files[i])
  namefile <- tidyName(namefile)#tidies the name if needed
  data <- filesIn(files[i])
  maxC <- maxCount(as.data.frame(sum(data$strand == "+"), row.names = namefile),
                   as.data.frame(sum(data$strand == "-"), row.names = namefile))
  ifelse(exists("MaxCount"), MaxCount <- rbind(MaxCount, maxC), MaxCount <- maxC)
  readLength <- getGenLength(data)
  if (exists("BAMList")){ 
    BAMList[[namefile]] <- data
  }else{  
    BAMList <-list()
    BAMList[[namefile]] <- data
  }
  rm(data)
}
#make size distributions 
if(exists("SizeDistribution")){rm(SizeDistribution, NormaliseBy)}
for (i in 1:length(BAMList)){
  filename <- names(BAMList[i])
  data <- (BAMList[[i]])
  nam <- paste0(filename, "_Size_Distribution")
  dsz <- getSize(data, Count = c(18:40))
  NormFactor <- data.frame(sum(abs(dsz)))
  row.names(NormFactor) <- filename
  ifelse(exists("NormaliseBy"), 
         NormaliseBy <- bind_rows(NormaliseBy, NormFactor),
         NormaliseBy <- NormFactor)
  if(exists("SizeDistribution")){
    SizeDistribution[[nam]] <- dsz
  }else{
    SizeDistribution <- list()
    SizeDistribution[[nam]] <- dsz
  }
  rm(data, dsz, NormFactor) 
}
#plot size distributions
for(i in 1:length(SizeDistribution)){
  nam <- names(SizeDistribution[i])
  data <- SizeDistribution[[i]]
  dat <- data.frame(t(data))
  dat$x <- row.names(dat)
  gg <- ggplot(dat, aes(x = x) )+
    geom_col(aes( y = Pos),colour = "black", fill = group.colours["Pos"])+
    geom_col(aes( y = Neg),colour = "black", fill = group.colours["Neg"])+
    xlab ("Size")+
    ylab ("Count")+
    ylim((0-max(MaxCount)*1.1),(0+max(MaxCount)*1.1))+
    ggtitle(paste0(nam))+
    piMaker_theme
  plot(gg)
  #saveImage(paste0(nam))
  rm(nam,dat)
}
#make mean size distributions grouped by sample - average across the replicates per sample
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
#plot the means for the sample
if(exists("Summary_Plot")){
  for(i in 1:length(Summary_Plot)){
    data <- Summary_Plot[[i]]
    namefile <- names(Summary_Plot[i])
    gg <- ggplot(data = data)+
      geom_col( aes(x = data$x, y = Pos_Mean), fill = group.colours["Pos"], linewidth=0.5, colour="black", alpha=0.9)+
      geom_col( aes(x = data$x, y = Neg_Mean), fill = group.colours["Neg"],linewidth=0.5, colour="black", alpha=0.9)+
      geom_errorbar( aes(x = data$x, ymin = Pos_E_Min, ymax = Pos_E_Max), linewidth=1, colour="black", alpha=0.9 )+
      geom_errorbar( aes(x = data$x, ymin = Neg_E_Min, ymax = Neg_E_Max), linewidth=1, colour="black", alpha=0.9 )+
      xlab ("Size")+
      ylab ("Count")+
      ylim (-(max(MaxCount)*1.1),(max(MaxCount)*1.1))+
      ggtitle(paste0(namefile))+
      piMaker_theme+
      theme(legend.position = "none")
    plot(gg)
    #saveImage(paste(namefile))
  }
}
#make size distributions grouped by genome segment - splits data by genome segment 
if(exists("GenSizeDistribution")){rm(GenSizeDistribution, NormaliseByGen)}
for (i in 1:length(BAMList)){
  filename <- names(BAMList[i])
  data <- (BAMList[[i]])
  dSP <- split(data,data$rname)
  for(r in 1:length(dSP)){
    dat <- dSP[[r]]
    namfile <- names(dSP[r])
    nam <- paste0(filename, "_", namfile, "_Size_Distribution")
    dsz <- getSize(dat, Count = c(18:40))
    NormFactor <- data.frame(sum(abs(dsz)))
    row.names(NormFactor) <- namfile
    ifelse(exists("NormaliseByGen"), 
           NormaliseByGen <- bind_rows(NormaliseByGen, NormFactor),
           NormaliseByGen <- NormFactor)
    if(exists("GenSizeDistribution")){
      GenSizeDistribution[[nam]] <- dsz
    }else{
      GenSizeDistribution <- list()
      GenSizeDistribution[[nam]] <- dsz
    }
  }
  rm(data, dsz, NormFactor) 
}
#make mean size distributions grouped by genome segment - calculates means for each genome segment
if(exists("GenSummary_Plot")){rm(GenSummary_Plot, GenMaxCount)}
if(length(samples>1)){
  for(i in 1:length(refSeq)){
    namefile <- refSeq[i]
    filename <- paste0(namefile, "_Size_Distribution")
    data <- (GenSizeDistribution[ grepl((paste0(namefile)), names(GenSizeDistribution))])
    res <- c()
    for (i in 1:length(data)) {
      nam <- names(data[i])
      dat <- data.frame(t(data[[i]]))
      colnames(dat) <- c(paste0(nam, "_Pos"), paste0(nam, "_Neg"))
      ifelse(exists("res"),
             res <- bind_cols(res, dat),
             res <- dat)
      maxC <- max(res)
      ifelse(exists("GenMaxCount"), GenMaxCount <- rbind(GenMaxCount, maxC), GenMaxCount <- maxC)
    }
    res <- SizDistSum(res)
    if(exists("GenSummary_Plot")){
      GenSummary_Plot[[filename]] <- res
      rm(res)
    }else{
      GenSummary_Plot <- list()
      GenSummary_Plot[[filename]] <- res
      rm(res)
    }
  }
}
#plot the means by genome segment 
if(exists("GenSummary_Plot")){
  for(i in 1:length(GenSummary_Plot)){
    data <- GenSummary_Plot[[i]]
    filename <- names(GenSummary_Plot[i])
    gg <- ggplot(data = data)+
      geom_col( aes(x = data$x, y = Pos_Mean), fill = group.colours["Pos"], linewidth=0.5, colour="black", alpha=0.9)+
      geom_col( aes(x = data$x, y = Neg_Mean), fill = group.colours["Neg"],linewidth=0.5, colour="black", alpha=0.9)+
      geom_errorbar( aes(x = data$x, ymin = Pos_E_Min, ymax = Pos_E_Max), linewidth=1, colour="black", alpha=0.9 )+
      geom_errorbar( aes(x = data$x, ymin = Neg_E_Min, ymax = Neg_E_Max), linewidth=1, colour="black", alpha=0.9 )+
      xlab ("Size")+
      ylab ("Count")+
      ylim (-(max(GenMaxCount)*1.1),(max(GenMaxCount)*1.1))+
      ggtitle(paste0(filename))+
      piMaker_theme+
      theme(legend.position = "none")
    plot(gg)
    #saveImage(paste(filename))
  }
}
#split the data by the 21nt siRNAs and gets the coverage 
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
    cov <- getCoverage(datc, GenSize = 1:(paste0(rsq))) 
    covC <- c(maxCount(cov$Neg, cov$Pos))
    ifelse(exists("siCovCount"), siCovCount <- rbind(siCovCount, covC), siCovCount <- covC)
    if(exists("siRNA_Coverage_Plot")){
      siRNA_Coverage_Plot[[namb]] <- cov
      rm(cov)
    }else{
      siRNA_Coverage_Plot <- list()
      siRNA_Coverage_Plot[[namb]] <- cov
      rm(cov)
    }
  }
  rm(data,d21,dat,namb,datc,covC)
}
#plots the coverage 
for (i in 1:(length(siRNA_Coverage_Plot))) {
  filename <- names(siRNA_Coverage_Plot[i])
  #namefile <- 
  cvr <- data.frame(siRNA_Coverage_Plot[[i]])
  gg<- ggplot()+
    geom_line(data = cvr, aes(x= x, y = Pos ), colour = group.colours["Pos"], linewidth = 0.5)+
    geom_line(data = cvr, aes(x= x, y = Neg), colour = group.colours["Neg"], linewidth = 0.5)+
    geom_hline(yintercept = 0, linetype = "solid", colour = "grey")+
    ggtitle(paste0(filename))+
    ylim(-(max(siCovCount)*1.1),(max(siCovCount)*1.1))+
    xlab ("nt position")+
    piMaker_theme
  plot(gg)
  #saveImage(paste(filename))
  rm(cvr)
}
#normalises the data to take unequal read counts in account
if(exists("siRNA_Coverage_Plot_Normalised")){rm(siRNA_Coverage_Plot_Normalised, NormScale)}
for (i in 1:(length(siRNA_Coverage_Plot))){
  filename <- names(siRNA_Coverage_Plot[i])
  nam <- str_extract(filename, '[A-Za-z]+_[0-9]')
  NormTo <- NormaliseBy[grepl(nam, row.names(NormaliseBy)),]
  data <- siRNA_Coverage_Plot[[i]]
  data$Pos_Normalised <- NormaliseTo(data$Pos, 0, NormTo)
  data$Neg_Normalised <- NormaliseTo(data$Neg, 0, NormTo)
  NorScal <- max(data$Pos_Normalised, abs(data$Neg_Normalised))
  ifelse(exists("NormScale"), NormScale <- rbind(NormScale, NorScal), NormScale <- NorScal)
  if (exists("siRNA_Coverage_Plot_Normalised")){
    siRNA_Coverage_Plot_Normalised[[filename]] <- data
  } else {
    siRNA_Coverage_Plot_Normalised <- list()
    siRNA_Coverage_Plot_Normalised[[filename]] <- data
  }
  rm(data, NormTo, NorScal)
}
#plots the normalised data
for (i in 1:(length(siRNA_Coverage_Plot_Normalised))) { 
  filename <- names(siRNA_Coverage_Plot_Normalised[i])
  namefile <- (paste0(filename, "_Normalised"))
  data <- siRNA_Coverage_Plot_Normalised[[i]]
  gg<- ggplot(data)+
    geom_line(aes(x= x, y = Pos_Normalised), colour = group.colours["Pos"], linewidth = 0.5)+
    geom_line(aes(x= x, y = Neg_Normalised), colour = group.colours["Neg"], linewidth = 0.5)+
    geom_hline(yintercept = 0, linetype = "solid", colour = "grey")+
    ggtitle(paste0(namefile))+
    ylim(-max(NormScale)*1.1, max(NormScale)*1.1)+
    xlab ("nt position")+
    ylab(NULL)+
    piMaker_theme
  plot(gg)
  #saveImage(paste(namefile))
}
#make summary data - the mean and SD of the coverage for each genome segment by sample
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
  rm(data,datRseq,datn)
}  
#plots summary data        
if(exists("Summary_Plot")){
  for(i in 1:length(siRNA_Summary)){
    name <- names(siRNA_Summary[i])
    name <- str_extract(name, '(?<=_)[A-Z0-9]*')
    filename <- paste0(name, "_siSummary")
    data <- siRNA_Summary[[i]]
    gg<- ggplot(data)+
      geom_ribbon(aes(x = x, ymin = Pos_E_min, ymax = Pos_E_max), colour = group.colours["Pos"],
                  fill = group.colours["Pos"], alpha = 0.2, linewidth = 0.1)+
      geom_line(aes(x = x, y = Pos_Mean),colour = group.colours["Pos"], linewidth = 0.5)+
      geom_ribbon(aes(x = x, ymin = Neg_E_min, ymax = Neg_E_max), colour = group.colours["Neg"],
                  fill = group.colours["Neg"], alpha = 0.2, linewidth = 0.1)+
      geom_line(aes(x = x, y = Neg_Mean), colour = group.colours["Neg"], linewidth = 0.5)+
      geom_hline(yintercept = 0, linetype = "solid", colour = "grey")+
      ggtitle(paste0(filename))+
      ylim(-max(NormScale), max(NormScale))+
      xlab ("nt position")+
      guides(guide_legend, fill = NULL)+
      piMaker_theme+
      theme(legend.position = "none")
    plot(gg)
    #saveImage(paste(filename))
  }        
}
#add genome information to summary plot if wanted 
for(i in 1:length(siRNA_Summary)){
  name <- names(siRNA_Summary[i])
  name <- str_extract(name, '(?<=_)[A-Z0-9]*')
  filename <- paste0(name, "_siSummary")
  data <- siRNA_Summary[[i]]
  gendata <- KnownCoding[[name]]
  gg<- ggplot(data)+
    geom_ribbon(aes(x = x, ymin = Pos_E_min, ymax = Pos_E_max), colour = group.colours["Pos"],
                fill = group.colours["Pos"], alpha = 0.2, linewidth = 0.1)+
    geom_line(aes(x = x, y = Pos_Mean),colour = group.colours["Pos"], linewidth = 0.5)+
    geom_ribbon(aes(x = x, ymin = Neg_E_min, ymax = Neg_E_max), colour = group.colours["Neg"],
                fill = group.colours["Neg"], alpha = 0.2, linewidth = 0.1)+
    geom_line(aes(x = x, y = Neg_Mean), colour = group.colours["Neg"], linewidth = 0.5)+
    geom_hline(yintercept = 0, linetype = "solid", colour = "grey")+
    geom_rect(data = gendata, aes(xmin = start, xmax = end, ymin = ((0.5*max(NormScale)*(0.45*number))),
                                  ymax = ((0.5*max(NormScale)*(0.45*number))+(0.1*max(NormScale))), colour = "black"),
              fill = "lightblue", alpha = 0.25)+
    geom_text(data = gendata, label = gendata$name, 
              x = gendata$start, 
              y = ((0.5*(0.25*gendata$number)*max(NormScale))+0.005),
              nudge_x = 1, nudge_y = 1,
              check_overlap = F)+
    ggtitle(paste0(filename))+
    ylim(-max(NormScale), max(NormScale))+
    xlab ("nt position")+
    guides(guide_legend, fill = NULL)+
    piMaker_theme+
    theme(legend.position = "none")
  
  plot(gg)
  saveImage(paste(filename))
}        

#piRNAs####
#load functions and colours
#function to find overlaps in piRNAs 
#these colours match standard peak viewers
DNA_col_scheme =  make_col_scheme(chars=c('A', 'U', 'G', 'C'), groups= NULL, 
                                  cols=c('green3', 'red2', 'black', 'blue'), name='DNA_col_scheme')

Sizes <- c(24:29) #put in the piRNA sizes we are looking for here
#21nt siRNAs are dealt with separately for overlap and z score calculations
#make long frames of the sequences, first filter the data sets by read length to reduce computing requirements
if(exists("piList")){rm(piList)}+
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
    namp <- (paste0(filename, "_piSequences"))
    nams <- (paste0(filename, "_siSequences"))
    if(exists("piList")){
      piList[[namp]] <- datp
    }else{
      piList <- list()
      piList[[namp]] <- datp
    }
    if(exists("siList")){
      siList[[nams]] <- dats
    }else{
      siList <- list()
      siList[[nams]] <- dats
    }
  }   
  rm(data,datl,datp,dats,namp,nams)
}
#This section deals with all the reads combined - 'overall data', if you want individual plots and mean/SD data, go to the next section!####
#now combine into single sets for each sample for the piRNAs
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
  rm(dat, data)
}
#now combine into single sets for each sample for the siRNAs
if(exists("siSeqList")){rm(siSeqList)}
for (i in 1:length(samples)) {
  namefile = samples[i]
  filename <- paste0(namefile, "_siRNAs")
  data <- siList[grepl(namefile, names(siList))]
  for (l in 1:length(data)){
    dat <- data[[l]]
    ifelse(exists("si_res"), si_res <- rbind(si_res,dat), si_res <- dat)
  }
  if(exists("siSeqList")){
    siSeqList[[filename]] <- si_res
    rm(si_res)
  }else{
    siSeqList <- list()
    siSeqList[[filename]] <- si_res
    rm(si_res)
  }
  rm(dat,data)
}
#make piRNA position matrix 
if(exists("piMatrix")){rm(piMatrix, piCount)}
for (i in 1:(length(piSeqList))) { 
  namfile <- names(piSeqList[i])
  data <- piSeqList[[i]]
  for (r in 1:(length(refSeq))) {
    rSeq <- refSeq[r]
    nam <- paste0(namfile, "_", rSeq)
    rsq <- makeRsq(rSeq)
    datr <- filter(data, data$rname == rSeq)
    freq <- piFreq(datr)
    piC <- maxCount(freq$Neg, freq$Pos)
    ifelse(exists("piCount"), piCount <- rbind(piCount, piC), piCount <- piC)
    if(exists("piMatrix")){
      piMatrix[[nam]] <- freq
    }else{
      piMatrix <- list()
      piMatrix[[nam]] <- freq
    }    
  }
  rm(data,datr,freq,piC)
}
#make siRNA position matrix
if(exists("siMatrix")){rm(siMatrix, siCount)}
for (i in 1:(length(siSeqList))) { 
  namfile <- names(siSeqList[i])
  data <- siSeqList[[i]]
  for (r in 1:(length(refSeq))) {
    rSeq <- refSeq[r]
    nam <- paste0(namfile, "_", rSeq)
    rsq <- makeRsq(rSeq)
    datr <- filter(data, data$rname == rSeq)
    freq <- piFreq(datr)
    siC <- maxCount(freq$Neg, freq$Pos)
    ifelse(exists("siCount"), siCount <- rbind(siCount, siC), siCount <- siC)
    if(exists("siMatrix")){
      siMatrix[[nam]] <- freq
    }else{
      siMatrix <- list()
      siMatrix[[nam]] <- freq
    }    
  }
  rm(data,datr,freq)
}
# get the coverage of the piRNA of size Sizes 
if(exists("piRNA_Coverage_Plot")){rm(piRNA_Coverage_Plot, piCovCount)}
for (i in 1:length(BAMList)){
  filename <- names(BAMList[i])
  data <- (BAMList[[i]])
  dpiRNA <- dplyr::filter(data, data$qwidth %in% c(Sizes))
  dpiRNA <- split(dpiRNA, dpiRNA$rname)
  for (n in 1:length(dpiRNA)){
    rSeq <- names(dpiRNA[n])
    rsq <- makeRsq(rSeq)
    namb <- paste0(names(BAMList[i]), "_", names(dpiRNA[n]))
    dat <- dpiRNA[[n]]
    datc <- coverPiMatrix(dat)
    datc <- datc[,-2]
    cov <- getCoverage(datc, GenSize = 1:paste(rsq)) 
    covP <- c(maxCount(cov$Neg, cov$Pos))
    ifelse(exists("piCovCount"), piCovCount <- rbind(piCovCount, covP), piCovCount <- covP)
    if(exists("piRNA_Coverage_Plot")){
      piRNA_Coverage_Plot[[namb]] <- cov
      rm(cov)
    }else{
      piRNA_Coverage_Plot <- list()
      piRNA_Coverage_Plot[[namb]] <- cov
      rm(cov)
    }
  }
  rm(data,dpiRNA,dat,namb,datc,covP)
}
#plot the coverage 
for (i in 1:(length(piRNA_Coverage_Plot))) {
  name <- names(piRNA_Coverage_Plot[i])
  filename <- paste0(name, "_piCoverage")
  cvr <- data.frame(piRNA_Coverage_Plot[[i]])
  gg<- ggplot()+
    geom_line(data = cvr, aes(x= x, y = Pos ), colour = group.colours["Pos"], linewidth = 0.5)+
    geom_line(data = cvr, aes(x= x, y = Neg), colour = group.colours["Neg"], linewidth = 0.5)+
    geom_hline(yintercept = 0, linetype = "solid", colour = "grey")+
    ggtitle(paste0(filename))+
    ylim(-(max(piCovCount)*1.1),(max(piCovCount)*1.1))+
    xlab ("nt position")+
    piMaker_theme
  plot(gg)
  #saveImage(paste(filename))
  rm(cvr)
}
#normalise the data to take unequal read counts in account 
if (exists("piRNA_Coverage_Plot_Normalised")){rm(piRNA_Coverage_Plot_Normalised)}
for (i in 1:(length(piRNA_Coverage_Plot))){
  filename = names(piRNA_Coverage_Plot[i])
  nam <- str_extract(filename, '[A-Za-z]+_[0-9]')
  NormTo <- NormaliseBy[grepl(nam, row.names(NormaliseBy)),]
  data <- (piRNA_Coverage_Plot[[i]])
  data$Pos_Normalised <- NormaliseTo(data$Pos, 0, NormTo)
  data$Neg_Normalised <- NormaliseTo(data$Neg, 0, NormTo)
  if (exists("piRNA_Coverage_Plot_Normalised")){
    piRNA_Coverage_Plot_Normalised[[filename]] <- data
  } else {
    piRNA_Coverage_Plot_Normalised <- list()
    piRNA_Coverage_Plot_Normalised[[filename]] <- data
  }
  rm(data)
}
#plot the normalised data
for (i in 1:(length(piRNA_Coverage_Plot_Normalised))) { 
  filename <- names(piRNA_Coverage_Plot_Normalised[i])
  namefile <- (paste0(filename, "_pi_Normalised"))
  data <- siRNA_Coverage_Plot_Normalised[[i]]
  gg<- ggplot(data)+
    geom_line(aes(x= x, y = Pos_Normalised), colour = group.colours["Pos"], linewidth = 0.5)+
    geom_line(aes(x= x, y = Neg_Normalised), colour = group.colours["Neg"], linewidth = 0.5)+
    geom_hline(yintercept = 0, linetype = "solid", colour = "grey")+
    ggtitle(paste0(namefile))+
    ylim(-max(NormScale)*1.1, max(NormScale)*1.1)+
    xlab ("nt position")+
    ylab(NULL)+
    piMaker_theme
  plot(gg)
  #saveImage(paste(namefile))
}
#make summary data
if(exists("piRNA_Summary")){rm(piRNA_Summary)}
for (i in 1:length(samples)) {
  namefile = samples[i]
  data <- piRNA_Coverage_Plot_Normalised[grepl(namefile, names(piRNA_Coverage_Plot_Normalised))]
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
    if(exists("piRNA_Summary")){
      piRNA_Summary[[filename]] <- Sumry
      rm(Sumry)
    }else{
      piRNA_Summary <- list()
      piRNA_Summary[[filename]] <- Sumry
      rm(Sumry)
    }
  }
  rm(data,datRseq,datn)
}  
#plot summary data       
if(exists("Summary_Plot")){
  for(i in 1:length(piRNA_Summary)){
    name <- names(piRNA_Summary[i])
    name <- str_extract(name, '(?<=_)[A-Z0-9]*')
    gendata <- KnownCoding[[name]]
    filename <- paste0(name, "_piSummary")
    data <- piRNA_Summary[[i]]
    gg<- ggplot(data)+
      geom_ribbon(aes(x = x, ymin = Pos_E_min, ymax = Pos_E_max), colour = group.colours["Pos"],
                  fill = group.colours["Pos"], alpha = 0.2, linewidth = 0.1)+
      geom_line(aes(x = x, y = Pos_Mean),colour = group.colours["Pos"], linewidth = 0.5)+
      geom_ribbon(aes(x = x, ymin = Neg_E_min, ymax = Neg_E_max), colour = group.colours["Neg"],
                  fill = group.colours["Neg"], alpha = 0.2, linewidth = 0.1)+
      geom_line(aes(x = x, y = Neg_Mean), colour = group.colours["Neg"], linewidth = 0.5)+
      geom_hline(yintercept = 0, linetype = "solid", colour = "grey")+
      geom_rect(data = gendata, aes(xmin = start, xmax = end, ymin = ((0.0005*number)+0.0045),
                                    ymax = ((0.0005*number)+0.005), colour = "black"),
                fill = "lightblue", alpha = 0.25)+
      geom_text(data = gendata, label = gendata$name, 
                x = gendata$start, 
                y = ((0.0005*gendata$number)+0.005),
                nudge_x = 1, nudge_y = 1,
                check_overlap = F)+
      ggtitle(paste0(filename))+
      ylim(-max(NormScale), max(NormScale))+
      xlab ("nt position")+
      guides(guide_legend, fill = NULL)+
      piMaker_theme+
      theme(legend.position = "none")
    plot(gg)
    #saveImage(paste(filename))
  }
}

#this is the first analysis now complete - size distribution and coverage of siRNA and piRNA
#adjust this as necessary. This plot is currently designed for 2 genome segments but the code
#is easily modified to include more by copy/pasting and changing the layout
#make a plot of this data:

finalPlot <- CoverageSummaryPlot(GenSummary_Plot,  siRNA_Summary, piRNA_Summary)

#plot piRNA position matrix - we don't really need to do this for the siRNAs as this is covered previously!
for (i in 1: length(piMatrix)){
  data <- piMatrix[[i]]
  filename <- names(piMatrix[i])
  rsq <- str_extract(filename,
                     '([A-Za-z].)([0-9]+)')
  rsq <- makeRsq(rsq)
  gg <- ggplot(data)+
    geom_hline(yintercept = 0, linetype = "solid", colour = "grey")+
    geom_col(aes(x = pos, y = Pos), colour = group.colours["Pos"])+
    geom_col(aes(x = pos, y = Neg), colour = group.colours["Neg"])+
    ggtitle(paste0(filename))+
    ylim(-max(piCount,siCount),max(piCount,siCount))+
    xlim(0,rsq)+
    piMaker_theme+
    theme(legend.position = "none")
  plot(gg)
  #saveImage(paste0(filename))
}
#but the code is here anyway
for (i in 1: length(siMatrix)){
  data <- siMatrix[[i]]
  filename <- names(siMatrix[i])
  rsq <- str_extract(filename,
                     '([A-Za-z].)([0-9]+)')
  rsq <- makeRsq(rsq)
  gg <- ggplot(data)+
    geom_hline(yintercept = 0, linetype = "solid", colour = "grey")+
    geom_col(aes(x = pos, y = Pos), colour = group.colours["Pos"])+
    geom_col(aes(x = pos, y = Neg), colour = group.colours["Neg"])+
    ggtitle(paste0(filename))+
    ylim(-max(piCount,siCount),max(piCount,siCount))+
    xlim(0,rsq)+
    piMaker_theme+
    theme(legend.position = "none")
  plot(gg)
  #saveImage(paste0(filename))
}
#now split the lists into respective piRNAs and plot, this will show if any sizes have the characteristic A/U overlap
if(exists("piRNA_List")){rm(piRNA_List)}
for (i in 1:(length(piSeqList))) {
  namfile <- names(piSeqList[i])
  data <- piSeqList[[i]]
  data$seq <- gsub("T", "U", data$seq)
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
      datjNeg$RC <- reverseComplement(RNAStringSet(datjNeg$seq))
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
    }
  }
  rm(data, datj, datjNeg, datjPos, datsz, piPos, piNeg)
} 
#make the sequence logos 
for (i in 1:(length(samples))){
  snam <- samples[[i]]
  for (r in 1:length(refSeq)){
    rname <- refSeq[[r]]
   # rname <- refNames[rnam]
    for (sz in 1:length(Sizes)){ 
      sznam <- paste0("_",Sizes[[sz]])
      namj <- paste0(rname,"_",sznam)
      namP <- paste0(snam, "_", rname, sznam, "_Pos")
      namN <- paste0(snam, "_", rname, sznam, "_Neg")
      dat <- piRNA_List[grepl(snam, names(piRNA_List)) & 
                          grepl(rname, names(piRNA_List)) &
                          grepl(sznam, names(piRNA_List))]
      piPos <- as.data.frame(dat[grepl("Pos", names(dat))])
      piNeg <- as.data.frame(dat[grepl("Neg", names(dat))])
      rm(dat)
      
      gp <- ggplot()+
        geom_logo(piPos, method = 'prob', seq_type = "rna", font = "helvetica_bold", col_scheme = DNA_col_scheme)+
        theme_logo()+
        ggtitle(namP)+
        piMaker_theme
      gn <- ggplot()+
        geom_logo(piNeg, method = 'prob', seq_type = "rna", font = "helvetica_bold", col_scheme = DNA_col_scheme)+
        theme_logo()+
        ggtitle(namN)+
        piMaker_theme
      fig <- ggarrange(gp,gn, nrow = 2)
      plot(fig)
      saveImage(paste0(namj))
    }
  }
}
#make the position matrix for the reads, counts how many of each read length at each position
#this can take a while! 
if(exists("piTally_List")){rm(piTally_List)}
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
        datS$pos <- (datS$pos + datS$qwidth)-1#this makes sure the 5' of the negative strand is in the correct position
      }
      result <- makeTally(datS, GenPosn, c(24:29))
      colnames(result) <- c("pos", c(24:29))
      if(exists("piTally_List")){
        piTally_List[[namB]] <- result
        rm(result)
      }else{
        piTally_List <- list()
        piTally_List[[namB]] <- result
        rm(result)
      }
    }
  }
  rm(datS)
}
#and do the same for the siRNAs -  this is a lot quicker!
if(exists("siTally_List")){rm(siTally_List)}
for (i in 1:(length(siSeqList))) {
  namfile <- names(siSeqList[i])
  data <- siSeqList[[i]]
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
        datS$pos <- (datS$pos + datS$qwidth)-1#this makes sure the 5' of the negative strand is in the correct position
      }
      result <- makeTally(datS, GenPosn, c(21))
      colnames(result) <- c("pos", c(21))
      if(exists("siTally_List")){
        siTally_List[[namB]] <- result
        rm(result)
      }else{
        siTally_List <- list()
        siTally_List[[namB]] <- result
        rm(result)
      }
    }
  }
  rm(datS)
}
#join the positive and negative reads into a single matrix for the piRNAs 
if(exists("piTally_List_Matrix")){rm(piTally_List_Matrix, piMax)}
for (i in 1:length(samples)){
  smp <- samples[i]
  namS <- paste0(smp)
  for (r in 1:length(refSeq)){
    rSeq <- refSeq[r]
    namRP <- paste0(smp, "_piRNAs_", rSeq,  "_Pos_Tally")
    namRN <- paste0(smp, "_piRNAs_", rSeq, "_Neg_Tally")
    namRS <- paste0(smp, "_", rSeq, "_Tally_Matrix")
    datP <- piTally_List[[(paste0(namRP))]]
    colnames(datP) <- paste0("Pos_", colnames(datP))
    datN <- piTally_List[[(paste0(namRN))]]
    colnames(datN) <- paste0("Neg_", colnames(datN))
    datC <- cbind(datN, datP)
    piTalMax <- maxCount( max(datC[,c(2:7)]), max(datC[,c(9:14)]))
    ifelse(exists("piMax"), piMax <- rbind(piMax,piTalMax), piMax <- piTalMax)
    if(exists("piTally_List_Matrix")){
      piTally_List_Matrix[[namRS]] <- datC
    } else{
      piTally_List_Matrix <- list()
      piTally_List_Matrix[[namRS]] <- datC
    }
  }
  rm(smp,rSeq,namRP,namRN,namRS,datP,datN,datC,piTalMax)
}
#join the positive and negative reads into a single matrix for the siRNAs
if(exists("siTally_List_Matrix")){rm(siTally_List_Matrix, siMax)}
for (i in 1:length(samples)){
  smp <- samples[i]
  namS <- paste0(smp)
  for (r in 1:length(refSeq)){
    rSeq <- refSeq[r]
    namRP <- paste0(smp, "_siRNAs_", rSeq,  "_Pos_Tally")
    namRN <- paste0(smp, "_siRNAs_", rSeq, "_Neg_Tally")
    namRS <- paste0(smp, "_", rSeq, "_Tally_Matrix")
    datP <- siTally_List[[(paste0(namRP))]]
    colnames(datP) <- paste0("Pos_", colnames(datP))
    datN <- siTally_List[[(paste0(namRN))]]
    colnames(datN) <- paste0("Neg_", colnames(datN))
    datC <- cbind(datN, datP)
    siTalMax <- maxCount( max(datC[,2]), max(datC[,4]))
    ifelse(exists("siMax"), siMax <- rbind(siMax, siTalMax), siMax <- siTalMax)
    if(exists("siTally_List_Matrix")){
      siTally_List_Matrix[[namRS]] <- datC
    } else{
      siTally_List_Matrix <- list()
      siTally_List_Matrix[[namRS]] <- datC
    }
  }
  rm(smp,rSeq,namRP,namRN,namRS,datP,datN,datC,siTalMax)
}
#map the piRNA overlaps
for(i in 1:length(piTally_List_Matrix)){
  namfile <- names(piTally_List_Matrix[i])
  dat <- piTally_List_Matrix[[i]]
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
                                   xmax = Pos_24_xend, fill = "24"), colour = piRNA.colours["24"], alpha = 0.4)+
    geom_rect(data = dat_2429, aes(ymin = 0, ymax = - Neg_24, xmin = Neg_24_x,
                                   xmax = Neg_24_xend, fill = "24"), colour = piRNA.colours["24"], alpha = 0.4)+ 
    
    
    geom_rect(data = dat_2429, aes(ymin = 0, ymax = Pos_25, xmin = Pos_25_x,
                                   xmax = Pos_25_xend, fill = "25"), colour = piRNA.colours["25"],  alpha = 0.4)+
    geom_rect(data = dat_2429, aes(ymin = 0, ymax = - Neg_25, xmin = Neg_25_x,
                                   xmax = Neg_25_xend, fill = "25"), colour = piRNA.colours["25"],  alpha = 0.4)+
    
    
    geom_rect(data = dat_2429, aes(ymin = 0, ymax = Pos_26, xmin = Pos_26_x,
                                   xmax = Pos_26_xend, fill = "26"), colour = piRNA.colours["26"],  alpha = 0.4)+
    geom_rect(data = dat_2429, aes(ymin = 0, ymax = - Neg_26, xmin = Neg_26_x,
                                   xmax = Neg_26_xend, fill = "26"), colour = piRNA.colours["26"],  alpha = 0.4)+
    
    
    geom_rect(data = dat_2429, aes(ymin = 0, ymax = Pos_27, xmin = Pos_27_x,
                                   xmax = Pos_27_xend, fill = "27"), colour = piRNA.colours["27"],  alpha = 0.4)+
    geom_rect(data = dat_2429, aes(ymin = 0, ymax = - Neg_27, xmin = Neg_27_x,
                                   xmax = Neg_27_xend, fill = "27"), colour = piRNA.colours["27"],  alpha = 0.4)+
    
    
    geom_rect(data = dat_2429, aes(ymin = 0, ymax = Pos_28, xmin = Pos_28_x,
                                   xmax = Pos_28_xend, fill = "28"), colour = piRNA.colours["28"],  alpha = 0.4)+
    geom_rect(data = dat_2429, aes(ymin = 0, ymax = - Neg_28, xmin = Neg_28_x,
                                   xmax = Neg_28_xend, fill = "28"), colour = piRNA.colours["28"],  alpha = 0.4)+
    
    geom_rect(data = dat_2429, aes(ymin = 0, ymax = Pos_29, xmin = Pos_29_x,
                                   xmax = Pos_25_xend, fill = "29"), colour = piRNA.colours["29"],  alpha = 0.4)+
    geom_rect(data = dat_2429, aes(ymin = 0, ymax = - Neg_29, xmin = Neg_29_x,
                                   xmax = Neg_25_xend, fill = "29"), colour = piRNA.colours["29"],  alpha = 0.4)+
    ylim(-max(piMax),max(piMax))+
    ggtitle(paste(nam))+
    piMaker_theme+
    scale_fill_manual("Size", values = c(piRNA.colours))
  
  plot(gg)
  #saveImage(paste0(nam))
}    

#Next section uses the data to make summary figures combining all results
#making individual data sets for each read to give final mean and SD plots####

#define the overlaps to be counted:
Overlap <- c(1:21) #will count overlaps of this length 
#make piRNA position matrix for each sample 
if(exists("piMatrixIndividual")){rm(piMatrixIndividual, piCountInd)}
for (i in 1:(length(piList))) { 
  namfile <- names(piList[i])
  data <- piList[[i]]
  for (r in 1:(length(refSeq))) {
    rSeq <- refSeq[r]
    nam <- paste0(namfile, "_", rSeq)
    rsq <- makeRsq(rSeq)
    datr <- dplyr::filter(data, data$rname == rSeq)
    freq <- piFreq(datr)
    piC <- maxCount(freq$Neg, freq$Pos)
    ifelse(exists("piCountInd"), piCountInd <- rbind(piCountInd, piC), piCountInd <- piC)
    if(exists("piMatrixIndividual")){
      piMatrixIndividual[[nam]] <- freq
    }else{
      piMatrixIndividual <- list()
      piMatrixIndividual[[nam]] <- freq
    }    
  }
  rm(data,datr,freq,piC)
}  
#make siRNA position matrix for each sample 
if(exists("siMatrixIndividual")){rm(siMatrixIndividual, siCountInd)}
for (i in 1:(length(siList))) { 
  namfile <- names(siList[i])
  data <- siList[[i]]
  for (r in 1:(length(refSeq))) {
    rSeq <- refSeq[r]
    nam <- paste0(namfile, "_", rSeq)
    rsq <- makeRsq(rSeq)
    datr <- dplyr::filter(data, data$rname == rSeq)
    freq <- piFreq(datr)
    piC <- maxCount(freq$Neg, freq$Pos)
    ifelse(exists("siCountInd"), piCountInd <- rbind(piCountInd, piC), piCountInd <- piC)
    if(exists("siMatrixIndividual")){
      siMatrixIndividual[[nam]] <- freq
    }else{
      siMatrixIndividual <- list()
      siMatrixIndividual[[nam]] <- freq
    }    
  }
  rm(data,datr,freq,piC)
}   
#this can take a while! 
if(exists("piTally_List_Individual")){rm(piTally_List_Individual)}
for (i in 1:(length(piList))) {
  namfile <- names(piList[i])
  data <- piList[[i]]
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
        datS$pos <- (datS$pos + datS$qwidth)-1#this makes sure the 5' of the negative strand is in the correct position
      }
      result <- makeTally(datS, GenPosn, c(24:29))
      colnames(result) <- c("pos", c(24:29))
      if(exists("piTally_List_Individual")){
        piTally_List_Individual[[namB]] <- result
        rm(result)
      }else{
        piTally_List_Individual <- list()
        piTally_List_Individual[[namB]] <- result
        rm(result)
      }
    }
  }
  rm(datS)
}
if(exists("siTally_List_Individual")){rm(siTally_List_Individual)}
for (i in 1:(length(siList))) {
  namfile <- names(siList[i])
  data <- siList[[i]]
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
        datS$pos <- (datS$pos + datS$qwidth)-1#this makes sure the 5' of the negative strand is in the correct position
      }
      result <- makeTally(datS, GenPosn, c(21))
      colnames(result) <- c("pos", c(21))
      if(exists("siTally_List_Individual")){
        siTally_List_Individual[[namB]] <- result
        rm(result)
      }else{
        siTally_List_Individual <- list()
        siTally_List_Individual[[namB]] <- result
        rm(result)
      }
    }
  }
  rm(datS)
}
#join the positive and negative reads into a single matrix for the piRNAs 
if(exists("piTally_List_Matrix_Individual")){rm(piTally_List_Matrix_Individual, piMaxInd)}
for (i in 1:length(samples)){
  smp <- samples[i]
  namS <- paste0(smp)
  for (n in 1:3){
    num <- n
    for (r in 1:length(refSeq)){
      rSeq <- refSeq[r]
      namRP <- paste0(smp, "_", n, "_piSequences_", rSeq,  "_Pos_Tally")
      namRN <- paste0(smp, "_", n, "_piSequences_", rSeq, "_Neg_Tally")
      namRS <- paste0(smp, "_", n, "_", rSeq, "_Tally_Matrix")
      datP <- piTally_List_Individual[[(paste0(namRP))]]
      colnames(datP) <- paste0("Pos_", colnames(datP))
      datN <- piTally_List_Individual[[(paste0(namRN))]]
      colnames(datN) <- paste0("Neg_", colnames(datN))
      datC <- cbind(datN, datP)
      piTalMax <- maxCount( max(datC[,c(2:7)]), max(datC[,c(9:14)]))
      ifelse(exists("piMaxInd"), piMaxInd <- rbind(piMaxInd,piTalMax), piMaxInd <- piTalMax)
      
      if(exists("piTally_List_Matrix_Individual")){
        piTally_List_Matrix_Individual[[namRS]] <- datC
      } else{
        piTally_List_Matrix_Individual <- list()
        piTally_List_Matrix_Individual[[namRS]] <- datC
      }
    }
  }
  rm(smp,rSeq,namRP,namRN,namRS,datP,datN,datC,piTalMax)
}
#join the positive and negative reads into a single matrix for the siRNAs 
if(exists("siTally_List_Matrix_Individual")){rm(siTally_List_Matrix_Individual, siMaxInd)}
for (i in 1:length(samples)){
  smp <- samples[i]
  namS <- paste0(smp)
  for (n in 1:3){
    num <- n
    for (r in 1:length(refSeq)){
      rSeq <- refSeq[r]
      namRP <- paste0(smp, "_", n, "_siSequences_", rSeq,  "_Pos_Tally")
      namRN <- paste0(smp, "_", n, "_siSequences_", rSeq, "_Neg_Tally")
      namRS <- paste0(smp, "_", n, "_", rSeq, "_Tally_Matrix")
      datP <- siTally_List_Individual[[(paste0(namRP))]]
      colnames(datP) <- paste0("Pos_", colnames(datP))
      datN <- siTally_List_Individual[[(paste0(namRN))]]
      colnames(datN) <- paste0("Neg_", colnames(datN))
      datC <- cbind(datN, datP)
      siTalMax <- maxCount( max(datC[,c(2)]), max(datC[,c(4)]))
      ifelse(exists("siMaxInd"), siMaxInd <- rbind(siMaxInd,siTalMax), siMaxInd <- siTalMax)
      
      if(exists("siTally_List_Matrix_Individual")){
        siTally_List_Matrix_Individual[[namRS]] <- datC
      } else{
        siTally_List_Matrix_Individual <- list()
        siTally_List_Matrix_Individual[[namRS]] <- datC
      }
    }
  }
  rm(smp,rSeq,namRP,namRN,namRS,datP,datN,datC,siTalMax)
}
#map the piRNA overlaps for the individual samples -  including this as I haven't put it before!
for(i in 1:length(piTally_List_Matrix_Individual)){
  namfile <- names(piTally_List_Matrix_Individual[i])
  dat <- piTally_List_Matrix_Individual[[i]]
  nam <- paste0(namfile, "_2429")
  piPlot <- piMapper(dat, piMaxInd)
  plot(piPlot)
  #saveImage(paste0(nam))
}    
#run piCatcherDual for the pi
for(i in 1:length(piTally_List_Matrix_Individual)){
  namp <- names(piTally_List_Matrix_Individual[i])
  namp <- gsub("_Tally_Matrix","",namp)
  datp <- piTally_List_Matrix_Individual[[i]]
  results <- piCatcherDual(datp, Length_In = c(24:29), Target_Length = c(24:29), Overlap = c(1:21))
  ifelse( exists("piProMaxI"), 
          piProMaxI <- rbind( piProMaxI, max(results[,c(6,7)]) ),
          piProMaxI <- c( max(results[,c(6,7)]) ) )
  ifelse( exists("piProWeiMaxI"),
          piProWeiMaxI <- rbind( piProWeiMaxI, max(results[,c(8,9)]) ),
          piProWeiMaxI <- c( max(results[,c(8,9)]) )  )
  if(exists("pi_Signatures_Individual")){
    pi_Signatures_Individual[[namp]] <- results
    rm(results)
  }else{
    pi_Signatures_Individual <- list()
    pi_Signatures_Individual[[namp]] <- results
    rm(results)
  }
}
#run piCatcherDual for the si
for(i in 1:length(siTally_List_Matrix_Individual)){
  namp <- names(siTally_List_Matrix_Individual[i])
  namp <- gsub("_Tally_Matrix","",namp)
  datp <- siTally_List_Matrix_Individual[[i]]
  results <- piCatcherDual(datp, Length_In = c(21), Target_Length = c(21), Overlap = c(1:21))
  ifelse( exists("siProMaxI"), 
          siProMaxI <- rbind( siProMaxI, max(results[,c(6,7)]) ),
          siProMaxI <- c( max(results[,c(6,7)]) ) )
  ifelse( exists("siProWeiMaxI"),
          siProWeiMaxI <- rbind( siProWeiMaxI, max(results[,c(8,9)]) ),
          siProWeiMaxI <- c( max(results[,c(8,9)]) )  )
  if(exists("si_Signatures_Individual")){
    si_Signatures_Individual[[namp]] <- results
    rm(results)
  }else{
    si_Signatures_Individual <- list()
    si_Signatures_Individual[[namp]] <- results
    rm(results)
  }
}
#and make the barplot dataset
#make bar chart data
for (i in 1:length(pi_Signatures_Individual)){
  namp <- names(pi_Signatures_Individual[i])
  datp <- pi_Signatures_Individual[[i]]
  barMaxCountInd <- max(datp$Overlaps)
  ifelse(exists("barMaxInd"), barMaxInd <- rbind(barMaxInd, barMaxCountInd), barMaxInd <- barMaxCountInd)
  datp <- datp[,c(1:3)]
  dat1 <- datp[,c(1,2)]
  colnames(dat1) <- c("x", "Count")
  dat1$Group <- "Pos"
  dat2 <- datp[,c(1,3)]
  colnames(dat2) <- c("x", "Count")
  dat2$Group <- "Neg"
  res <- bind_rows(dat1, dat2)
  if(exists("pi_barplot_Ind")){
    pi_barplot_Ind[[namp]] <- res
    rm(res)
  }else{
    pi_barplot_Ind <- list()
    pi_barplot_Ind[[namp]] <- res
    rm(res)
  }
  rm(datp,dat1,dat2, barMaxCountInd)
}
for (i in 1:length(si_Signatures_Individual)){
  namp <- names(si_Signatures_Individual[i])
  datp <- si_Signatures_Individual[[i]]
  barMaxCountInd <- max(datp$Overlaps)
  ifelse(exists("sarMaxInd"), sarMaxInd <- rbind(sarMaxInd, barMaxCountInd), sarMaxInd <- barMaxCountInd)
  datp <- datp[,c(1:3)]
  dat1 <- datp[,c(1,2)]
  colnames(dat1) <- c("x", "Count")
  dat1$Group <- "Pos"
  dat2 <- datp[,c(1,3)]
  colnames(dat2) <- c("x", "Count")
  dat2$Group <- "Neg"
  res <- bind_rows(dat1, dat2)
  if(exists("si_barplot_Ind")){
    si_barplot_Ind[[namp]] <- res
    rm(res)
  }else{
    si_barplot_Ind <- list()
    si_barplot_Ind[[namp]] <- res
    rm(res)
  }
  rm(datp,dat1,dat2, barMaxCountInd)
}
#and plot for pi
for (i in 1:length(BAMList)){
  namfile <- names(BAMList[i])
  for(r in 1:length(refSeq)){
    Rnam <- refSeq[r]
    Line <- pi_Signatures_Individual [[grep(paste0(namfile, "_", Rnam ), names(pi_Signatures_Individual))]]
    Bar <- pi_barplot_Ind [[grep(paste0(namfile, "_", Rnam ), names(pi_barplot_Ind))]]
    namd <- paste0(namfile, "_", Rnam, "_pi")
    fig <-  MatrixPlot(Line, Bar, piProMaxI, piProWeiMaxI, barMaxInd)
    plot(fig)
    #saveImage(paste0(namd))
  }
}
#and plot for si
for (i in 1:length(BAMList)){
  namfile <- names(BAMList[i])
  for(r in 1:length(refSeq)){
    Rnam <- refSeq[r]
    Line <- si_Signatures_Individual [[grep(paste0(namfile, "_", Rnam ), names(pi_Signatures_Individual))]]
    Bar <- si_barplot_Ind [[grep(paste0(namfile, "_", Rnam ), names(pi_barplot_Ind))]]
    namd <- paste0(namfile, "_", Rnam, "_pi")
    fig <- MatrixPlot(Line, Bar, siProMaxI, siProWeiMaxI, sarMaxInd)
    plot(fig)
    #saveImage(paste0(namd))
  }
}
#Summary data with mean and standard deviation####
#making the summary and standard deviations for the pi
if(exists("pi_Final_Plot")){rm(pi_Final_Plot, piProMaxF, piProWeiMaxF)}
for (i in 1:length(refSeq)){
  ref <- refSeq[i]
  dat <- (pi_Signatures_Individual[grepl(ref, names(pi_Signatures_Individual))])
  nama <- paste0(ref)
  res <- data.frame(c(1:21))
  sz <- length(dat)
  len <- ncol(dat[[1]])
  for (l in 1:len){
    nam <- colnames(dat[[1]][l])
    for(s in 1:sz){
      resa <- dat[[s]][l]
      ifelse(exists("resb"), resb <- cbind(resb, resa), resb <- resa)
      if(s == sz){
        res[paste0(nam, "_Mean")] <- rowMeans(resb)
        res[paste0(nam, "_SD")] <- apply(resb,1,sd)
        res[paste0(nam, "_SD_Plus")] <-  res[paste0(nam, "_Mean")] + res[paste0(nam, "_SD")]
        res[paste0(nam, "_SD_Minus")] <-  res[paste0(nam, "_Mean")] - res[paste0(nam, "_SD")]
        rm(resb)
      }
    }
  }
  res <- res[,-c(1,3,4,5)]
  ifelse( exists("piProMaxF"), 
          piProMaxF <- rbind( piProMaxF, max(res[,c(18:25)]) ),
          piProMaxF <- c( max(res[,c(18:25)]) ) )
  ifelse( exists("piProWeiMaxF"),
          piProWeiMaxF <- rbind( piProWeiMaxF, max(res[,c(26:33)]) ),
          piProWeiMaxF <- c( max(res[,c(26:33)]) )  )
  FinProbScalePi <- cbind(piProMaxF,piProWeiMaxF)
  if(exists("pi_Final_Plot")){
    pi_Final_Plot[[nama]] <- res
    write.csv(res, paste0(nama, "_Pi_Final_Statistics.csv"), row.names = T)
    rm(res)
  }else{
    pi_Final_Plot <- list()
    pi_Final_Plot[[nama]] <- res
    write.csv(res, paste0(nama, "_Pi_Final_Statistics.csv"), row.names = T)
    rm(res)
  }
}
#making the summary and standard deviations for the si
if(exists("si_Final_Plot")){rm(si_Final_Plot, siProMaxF, siProWeiMaxF)}
for (i in 1:length(refSeq)){
  ref <- refSeq[i]
  dat <- (si_Signatures_Individual[grepl(ref, names(si_Signatures_Individual))])
  nama <- paste0(ref)
  res <- data.frame(c(1:21))
  sz <- length(dat)
  len <- ncol(dat[[1]])
  for (l in 1:len){
    nam <- colnames(dat[[1]][l])
    for(s in 1:sz){
      resa <- dat[[s]][l]
      ifelse(exists("resb"), resb <- cbind(resb, resa), resb <- resa)
      if(s == sz){
        res[paste0(nam, "_Mean")] <- rowMeans(resb)
        res[paste0(nam, "_SD")] <- apply(resb,1,sd)
        res[paste0(nam, "_SD_Plus")] <-  res[paste0(nam, "_Mean")] + res[paste0(nam, "_SD")]
        res[paste0(nam, "_SD_Minus")] <-  res[paste0(nam, "_Mean")] - res[paste0(nam, "_SD")]
        rm(resb)
      }
    }
  }
  res <- res[,-c(1,3,4,5)]
  ifelse( exists("siProMaxF"), 
          siProMaxF <- rbind( siProMaxF, max(res[,c(18:25)]) ),
          siProMaxF <- c( max(res[,c(18:25)]) ) )
  ifelse( exists("siProWeiMaxF"),
          siProWeiMaxF <- rbind( siProWeiMaxF, max(res[,c(26:33)]) ),
          siProWeiMaxF <- c( max(res[,c(26:33)]) )  )
  FinProbScaleSi <- cbind(siProMaxF,siProWeiMaxF)
  if(exists("si_Final_Plot")){
    si_Final_Plot[[nama]] <- res
    write.csv(res, paste0(nama, "_Si_Final_Statistics.csv"), row.names = T)
    rm(res)
  }else{
    si_Final_Plot <- list()
    si_Final_Plot[[nama]] <- res
    write.csv(res, paste0(nama, "_Si_Final_Statistics.csv"), row.names = T)
    rm(res)
  }
}
#making summary mappings of the piRNA to better match the summary mapping of the siRNAs if needed 
if(exists("pi_Final_Map_Plot")){rm(pi_Final_Map_Plot, piMapMaxF)}
for (i in 1:length(refSeq)){
  ref <- refSeq[i]
  dat <- (piTally_List_Matrix_Individual[grepl(ref, names(piTally_List_Matrix_Individual))])
  nama <- paste0(ref)
  sz <- length(dat)
  len <- ncol(dat[[1]])
  res <- data.frame(c(1:nrow(dat[[1]])))
  for (l in 1:len){
    nam <- colnames(dat[[1]][l])
    for(s in 1:sz){
      resa <- dat[[s]][l]
      ifelse(exists("resb"), resb <- cbind(resb, resa), resb <- resa)
      if(s == sz){
        res[paste0(nam, "_Mean")] <- rowMeans(resb)
        res[paste0(nam, "_SD")] <- apply(resb,1,sd)
        res[paste0(nam, "_SD_Plus")] <-  res[paste0(nam, "_Mean")] + res[paste0(nam, "_SD")]
        res[paste0(nam, "_SD_Minus")] <-  res[paste0(nam, "_Mean")] - res[paste0(nam, "_SD")]
        rm(resb)
      }
    }
  }
  res <- res[,-c(1,3:5,31:33)]
  ifelse( exists("piMapMaxF"), 
          piMapMaxF <- rbind( piMapMaxF, max(res[,c(2:25,27:50)]) ),
          piMapMaxF <- c( max(res[,c(2:25,27:50)]) ) )
  if(exists("pi_Final_Map_Plot")){
    pi_Final_Map_Plot[[nama]] <- res
    rm(res)
  }else{
    pi_Final_Map_Plot <- list()
    pi_Final_Map_Plot[[nama]] <- res
    rm(res)
  }
}
#and map
for(i in 1:length(pi_Final_Map_Plot)){
  namfile <- names(pi_Final_Map_Plot[i])
  dat <- pi_Final_Map_Plot[[i]]
  nam <- paste0(namfile, "2429_Summary")
  dat_2429 <- dat
  
  gg <- ggplot(data = dat_2429,)+
    geom_hline(yintercept = 0, linetype = "solid", colour = "grey")+
    
    geom_ribbon(aes(x = Pos_pos_Mean, ymin = Pos_24_SD_Minus, ymax = Pos_24_SD_Plus,
                    fill = "24"), colour = piRNA.colours["24"], alpha = 0.4)+
    geom_line(aes(x = Pos_pos_Mean, y = Pos_24_Mean, colour = "24"), colour = piRNA.colours["24"], alpha = 0.4)+ 
    geom_ribbon(aes(x = Neg_pos_Mean, ymin = Neg_24_SD_Minus, ymax = Neg_24_SD_Plus,
                    fill = "24"), colour = piRNA.colours["24"], alpha = 0.4)+
    geom_line(aes(x = Neg_pos_Mean, y = Neg_24_Mean, colour = "24"), colour = piRNA.colours["24"], alpha = 0.4)+
    
    
    geom_ribbon(aes(x = Pos_pos_Mean, ymin = Pos_25_SD_Minus, ymax = Pos_25_SD_Plus,
                    fill = "25"), colour = piRNA.colours["25"], alpha = 0.4)+
    geom_line(aes(x = Pos_pos_Mean, y = Pos_25_Mean, colour = "25"), colour = piRNA.colours["25"], alpha = 0.4)+ 
    geom_ribbon(aes(x = Neg_pos_Mean, ymin = Neg_25_SD_Minus, ymax = Neg_25_SD_Plus,
                    fill = "25"), colour = piRNA.colours["25"], alpha = 0.4)+
    geom_line(aes(x = Neg_pos_Mean, y = Neg_25_Mean, colour = "25"), colour = piRNA.colours["25"], alpha = 0.4)+
    
    
    geom_ribbon(aes(x = Pos_pos_Mean, ymin = Pos_26_SD_Minus, ymax = Pos_26_SD_Plus,
                    fill = "26"), colour = piRNA.colours["26"], alpha = 0.4)+
    geom_line(aes(x = Pos_pos_Mean, y = Pos_26_Mean, colour = "26"), colour = piRNA.colours["26"], alpha = 0.4)+ 
    geom_ribbon(aes(x = Neg_pos_Mean, ymin = Neg_26_SD_Minus, ymax = Neg_26_SD_Plus,
                    fill = "26"), colour = piRNA.colours["26"], alpha = 0.4)+
    geom_line(aes(x = Neg_pos_Mean, y = Neg_26_Mean, colour = "26"), colour = piRNA.colours["26"], alpha = 0.4)+
    
    
    geom_ribbon(aes(x = Pos_pos_Mean, ymin = Pos_27_SD_Minus, ymax = Pos_27_SD_Plus,
                    fill = "27"), colour = piRNA.colours["27"], alpha = 0.4)+
    geom_line(aes(x = Pos_pos_Mean, y = Pos_27_Mean, colour = "27"), colour = piRNA.colours["27"], alpha = 0.4)+ 
    geom_ribbon(aes(x = Neg_pos_Mean, ymin = Neg_27_SD_Minus, ymax = Neg_27_SD_Plus,
                    fill = "27"), colour = piRNA.colours["27"], alpha = 0.4)+
    geom_line(aes(x = Neg_pos_Mean, y = Neg_27_Mean, colour = "27"), colour = piRNA.colours["27"], alpha = 0.4)+
    
    
    geom_ribbon(aes(x = Pos_pos_Mean, ymin = Pos_28_SD_Minus, ymax = Pos_28_SD_Plus,
                    fill = "28"), colour = piRNA.colours["28"], alpha = 0.4)+
    geom_line(aes(x = Pos_pos_Mean, y = Pos_28_Mean, colour = "28"), colour = piRNA.colours["28"], alpha = 0.4)+ 
    geom_ribbon(aes(x = Neg_pos_Mean, ymin = Neg_28_SD_Minus, ymax = Neg_28_SD_Plus,
                    fill = "28"), colour = piRNA.colours["28"], alpha = 0.4)+
    geom_line(aes(x = Neg_pos_Mean, y = Neg_28_Mean, colour = "28"), colour = piRNA.colours["28"], alpha = 0.4)+
    
    geom_ribbon(aes(x = Pos_pos_Mean, ymin = Pos_29_SD_Minus, ymax = Pos_29_SD_Plus,
                    fill = "29"), colour = piRNA.colours["29"], alpha = 0.4)+
    geom_line(aes(x = Pos_pos_Mean, y = Pos_29_Mean, colour = "29"), colour = piRNA.colours["29"], alpha = 0.4)+ 
    geom_ribbon(aes(x = Neg_pos_Mean, ymin = Neg_29_SD_Minus, ymax = Neg_29_SD_Plus,
                    fill = "29"), colour = piRNA.colours["29"], alpha = 0.4)+
    geom_line(aes(x = Neg_pos_Mean, y = Neg_29_Mean, colour = "29"), colour = piRNA.colours["29"], alpha = 0.4)+
    
    ylim(-max(piMapMaxF),max(piMapMaxF))+
    ggtitle(paste(nam))+
    piMaker_theme+
    scale_fill_manual("Size", values = c(piRNA.colours))
  
  plot(gg)
  # saveImage(paste0(nam))
} 
#get the data for the final bar plots, error bars must be done separate unfortunately
if(exists("pi_barplot_Fin")){rm(pi_barplot_Fin, barMaxFin)}
for (i in 1:length(pi_Final_Plot)){
  namp <- names(pi_Final_Plot[i])
  datp <- pi_Final_Plot[[i]]
  barMaxCount <- max(datp$Overlaps_Mean)
  ifelse(exists("barMaxFin"), barMaxFin <- rbind(barMaxFin, barMaxCount), barMaxFin <- barMaxCount)
  dat1 <- datp[,c(1,2,3)]
  dat1$Cumulative <- dat1$Pos_count_Mean
  dat2 <- datp[,c(1,6,7)]
  dat2$Cumulative <- dat1$Pos_count_Mean + dat2$Neg_count_Mean
  colnames(dat1) <- c("x", "Count", "SD", "Cumulative")
  dat1$Group <- "Pos"
  colnames(dat2) <- c("x", "Count", "SD", "Cumulative")
  dat2$Group <- "Neg"
  res <- bind_rows(dat1, dat2)
  if(exists("pi_barplot_Fin")){
    pi_barplot_Fin[[namp]] <- res
    write.csv(res, paste0(namp, "_Pi_Final_Overlaps.csv"), row.names = T)
    rm(res)
  }else{
    pi_barplot_Fin <- list()
    pi_barplot_Fin[[namp]] <- res
    write.csv(res, paste0(namp, "_Pi_Final_Overlaps.csv"), row.names = T)
    rm(res)
  }
  rm(datp,dat1,dat2, barMaxCount)
}
#and for the siRNA
if(exists("si_barplot_Fin")){rm(si_barplot_Fin, sarMaxFin)}
for (i in 1:length(si_Final_Plot)){
  namp <- names(si_Final_Plot[i])
  datp <- si_Final_Plot[[i]]
  sarMaxCount <- max(datp$Overlaps_Mean) + max(datp$Overlaps_SD)
  ifelse(exists("sarMaxFin"), sarMaxFin <- rbind(sarMaxFin, sarMaxCount), sarMaxFin <- sarMaxCount)
  dat1 <- datp[,c(1,2,3)]
  dat1$Cumulative <- dat1$Pos_count_Mean
  dat2 <- datp[,c(1,6,7)]
  dat2$Cumulative <- dat1$Pos_count_Mean + dat2$Neg_count_Mean
  colnames(dat1) <- c("x", "Count", "SD", "Cumulative")
  dat1$Group <- "Pos"
  colnames(dat2) <- c("x", "Count", "SD", "Cumulative")
  dat2$Group <- "Neg"
  res <- bind_rows(dat1, dat2)
  if(exists("si_barplot_Fin")){
    si_barplot_Fin[[namp]] <- res
    write.csv(res, paste0(namp, "_Si_Final_Overlaps.csv"), row.names = T)
    rm(res)
  }else{
    si_barplot_Fin <- list()
    si_barplot_Fin[[namp]] <- res
    write.csv(res, paste0(namp, "_Si_Final_Overlaps.csv"), row.names = T)
    rm(res)
  }
  rm(datp,dat1,dat2, sarMaxCount)
}
#plot the final pi data
for(r in 1:length(refSeq)){
  Rnam <- refSeq[r]
  Line <- pi_Final_Plot [[grep(paste0(Rnam ), names(pi_Final_Plot))]]
  Bar <- pi_barplot_Fin [[grep(paste0(Rnam ), names(pi_barplot_Fin))]]
  namd <- paste0(Rnam, "_piRNAs")
  fig <- MatrixPlotSD(Line, Bar, FinProbScalePi, FinProbScalePi, barMaxFin)
  plot(fig)
  saveImage(paste0(namd))
}
#plot the final si data
for(r in 1:length(refSeq)){
  Rnam <- refSeq[r]
  Line <- si_Final_Plot [[grep(paste0(Rnam ), names(si_Final_Plot))]]
  Bar <- si_barplot_Fin [[grep(paste0(Rnam ), names(si_barplot_Fin))]]
  namd <- paste0(Rnam, "_siRNAs")
  fig <- MatrixPlotSD(Line, Bar, FinProbScaleSi, FinProbScaleSi, sarMaxFin)
  plot(fig)
  saveImage(paste0(namd))
}
#calculate the ratio of positive to negative sense siRNAs 
if(exists("siRatios")){rm(siRatios)}
for(i in 1:length(siTally_List_Matrix_Individual)){
  Rnam <- names(siTally_List_Matrix_Individual[i])
  dat <- siTally_List_Matrix_Individual [[i]]
  pos <- sum(dat$Pos_21)
  neg <- sum(dat$Neg_21)
  ratP <- pos/neg
  ratN <- neg/pos
  res <- cbind(Rnam, pos, neg, ratP, ratN)
 if(exists("siRatios")){
   siRatios <- rbind(siRatios, res)
 }else{
   siRatios <- res
 }
}
#calculate the ratio of positive to negative sense piRNAs 
if(exists("piRatios")){rm(piRatios)}
for(i in 1:length(piTally_List_Matrix_Individual)){
  Rnam <- names(piTally_List_Matrix_Individual[i])
  dat <- piTally_List_Matrix_Individual [[i]]
  pos <- sum(dat[,9:14])
  neg <- sum(dat[,2:7])
  ratP <- pos/neg
  ratN <- neg/pos
  res <- cbind(Rnam, pos, neg, ratP, ratN)
  if(exists("piRatios")){
    piRatios <- rbind(piRatios, res)
  }else{
    piRatios <- res
  }
}
#make the means and SD 
for(i in 1:length(refSeq)){
  refNam <- refSeq[i]
  
  siDat <- data.frame(siRatios[grep(refNam, siRatios),])
  piDat <- data.frame(piRatios[grep(refNam, piRatios),])
  
  siMePos <- mean(as.numeric(siDat$ratP))
  siMeNeg <- mean(as.numeric(siDat$ratN))
  piMePos <- mean(as.numeric(piDat$ratP))
  piMeNeg <- mean(as.numeric(piDat$ratN))
  
  siSDPos <- sd(as.numeric(siDat$ratP))
  siSDNeg <- sd(as.numeric(siDat$ratN))
  piSDPos <- sd(as.numeric(piDat$ratP))
  piSDNeg <- sd(as.numeric(piDat$ratN))
  
  res <- cbind(refNam, siMePos, siSDPos, siMeNeg, siSDNeg, piMePos, piSDPos, piMeNeg, piSDNeg)
  
  if(exists("Final_Ratios")){
    Final_Ratios <- rbind(Final_Ratios, res)
  }else{
    Final_Ratios <- res
  }
}
#FIN
