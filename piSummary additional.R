
if(exists("piRNA_Coverage_Plot")){rm(piRNA_Coverage_Plot, piCovCount, NormCount)}
for (i in 1:length(BAMList)){
  filename <- names(BAMList[i])
  data <- (BAMList[[i]])
  Posdata <- dplyr::filter(data, data$strand == "+")
  MinMaxPos <- data.frame(countMatrix(Posdata$qwidth, Count = c(24:29)))
  colnames(MinMaxPos) <- paste0(filename, "_Pos")
  Negdata <- dplyr::filter(data, data$strand == "-")
  MinMaxNeg <- data.frame(countMatrix(Negdata$qwidth, Count = c(24:29)))
  colnames(MinMaxNeg) <- paste0(filename, "_Neg")
  MinMax <- cbind(MinMaxPos,MinMaxNeg)
  ifelse(exists("NormCount"), NormCount <- cbind(NormCount, MinMax), NormCount <- MinMax)
  for (l in 1:length(Sizes)){
    sz <- Sizes[[l]]
    datSz <- filter(data, data$qwidth == sz)
    datrSeq <- split(datSz, datSz$rname)
    for (n in 1:length(datrSeq)){
      rSeq <- names(datrSeq[n])
      rsq <- makeRsq(rSeq)
      namb <- paste0(names(BAMList[i]), "_", names(datrSeq[n]), "_", sz)
      dat <- datrSeq[[n]]
      datc <- coverMatrix(dat, Length = sz)
      cov <- getCoverage(datc, Count = 1:(paste0(rsq))) 
      covC <- c(maxCount(cov$Neg, cov$Pos))
      ifelse(exists("piCovCount"), piCovCount <- rbind(piCovCount, covC), piCovCount <- covC)
      
      if(exists("piRNA_Coverage_Plot")){
        piRNA_Coverage_Plot[[namb]] <- cov
        rm(cov)
      }else{
        piRNA_Coverage_Plot <- list()
        piRNA_Coverage_Plot[[namb]] <- cov
        rm(cov)
      }
    }
  }
  rm(data,dat,namb,datc,covC,MinMaxNeg,MinMaxPos,MinMax,Negdata,Posdata)
}
#normalise the data to take unequal read counts in account
if(exists("piRNA_Coverage_Plot_Normalised")){rm(piRNA_Coverage_Plot_Normalised, piNormPlotCount)}
for (i in 1:(length(piRNA_Coverage_Plot))){
  filename = names(piRNA_Coverage_Plot[i])
  namefile <- str_extract(filename, '[A-Za-z]+_[0-9]')
  sz <- str_extract(filename,'_[0-9]{2}')
  sz <- str_extract(sz,'[0-9]{2}')
  normTo <- NormCount[grepl(namefile, colnames(NormCount))] 
  normTo <- normTo[sz,] 
  data <- (piRNA_Coverage_Plot[[i]])
  data$Pos_Normalised <- NormaliseTo(data$Pos, 0, max(normTo))
  data$Neg_Normalised <- NormaliseTo(data$Neg, 0, max(normTo))
  ifelse(exists("piNormPlotCount"), piNormPlotCount <- rbind(piNormPlotCount, max(abs(data[,4:5]))),
         piNormPlotCount <- max(abs(data[,4:5])))
  if (exists("piRNA_Coverage_Plot_Normalised")){
    piRNA_Coverage_Plot_Normalised[[filename]] <- data
  } else {
    piRNA_Coverage_Plot_Normalised <- list()
    piRNA_Coverage_Plot_Normalised[[filename]] <- data
  }
  rm(data,sz,normTo)
}
if(exists("piRNA_Summary")){rm(piRNA_Summary)}
for (i in 1:length(samples)) {
  namefile = samples[i]
  data <- piRNA_Coverage_Plot_Normalised[grepl(namefile, names(piRNA_Coverage_Plot_Normalised))] 
     
     for (j in 1:length(refSeq)) {
        rSeq <- refSeq[j]
        datRseq <- data[grepl(rSeq, names(data))] %>% keep( ~ !is.null(.) ) 
        
        for(s in 1:length(Sizes)){
          sz <- paste0("_", Sizes[[s]])
          datSz <- datRseq[grepl(sz, names(datRseq))] 
          
          for (n in 1:length(datSz)) {
            namen <- names(datSz[n])
            namen <- str_extract(namen, '[A-Za-z0-9]+_[0-9]')
            datn <- as.data.frame(datSz [[n]] [,c("Pos_Normalised","Neg_Normalised")] ) 
            colnames(datn) <- paste0(namen, "_", colnames(datn))
            ifelse(exists("Sumry"), 
                   Sumry <- cbind(Sumry, datn),
                   Sumry <- (datn))
          }
        filename <- paste0(namefile, "_", rSeq, "_", sz)
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
  }
  rm(data,datRseq,datn)
}  
#plot summary data        
  for(i in 1:length(refSeq)){
    rSeq <- refSeq[[i]]
    filename <- paste0(rSeq, "_piRNAs")
    data <- piRNA_Summary[grepl(rSeq, names(piRNA_Summary))]
    
    dat24 <- data[[grep("_24", names(data))]]
    dat25 <- data[[grep("_25", names(data))]]
    dat26 <- data[[grep("_26", names(data))]]
    dat27 <- data[[grep("_27", names(data))]]
    dat28 <- data[[grep("_28", names(data))]]
    dat29 <- data[[grep("_29", names(data))]]

    for(n in 24:29){
      dat_2429[paste0("Pos_",n, "_x")] <- as.numeric(ifelse( (dat_2429[(paste0("Pos_",n))] >0), 
                                                             dat_2429$Pos_pos, NA ))
      dat_2429[paste0("Pos_",n, "_xend")] <- dat_2429[paste0("Pos_",n, "_x")] + n
      dat_2429[paste0("Neg_",n, "_x")] <- as.numeric(ifelse((dat_2429[(paste0("Neg_",n))] >0), 
                                                            dat_2429$Neg_pos, NA))
      dat_2429[paste0("Neg_",n, "_xend")] <- dat_2429[paste0("Neg_",n, "_x")] + n
    }
    
    
    
    
    
    gg<- ggplot()+
      
      geom_hline(yintercept = 0, linetype = "solid", colour = "grey")+
      
      geom_ribbon(data = dat24, aes(x = x, ymin = Pos_E_min, ymax = Pos_E_max), colour = piRNA.colours["24"],
                  fill = piRNA.colours["24"],  alpha = 0.5, linewidth = NULL)+
      geom_line(data = dat24, aes(x = x, y = Pos_Mean), colour = piRNA.colours["24"], linewidth = 1)+
      geom_ribbon(data = dat24, aes(x = x, ymin = Neg_E_min, ymax = Neg_E_max), colour = piRNA.colours["24"],
                  fill = piRNA.colours["24"], alpha = 0.5, linewidth = NULL)+
      geom_line(data = dat24, aes(x = x, y = Neg_Mean), colour = piRNA.colours["24"], linewidth = 1)+
      
      geom_ribbon(data = dat25, aes(x = x, ymin = Pos_E_min, ymax = Pos_E_max), colour = piRNA.colours["25"],
                  fill = piRNA.colours["25"],  alpha = 0.5, linewidth = NULL)+
      geom_line(data = dat25, aes(x = x, y = Pos_Mean), colour = piRNA.colours["25"], linewidth = 1)+
      geom_ribbon(data = dat25, aes(x = x, ymin = Neg_E_min, ymax = Neg_E_max), colour = piRNA.colours["25"],
                  fill = piRNA.colours["25"], alpha = 0.5, linewidth = NULL)+
      geom_line(data = dat25, aes(x = x, y = Neg_Mean), colour = piRNA.colours["25"], linewidth = 1)+
      
      geom_ribbon(data = dat26, aes(x = x, ymin = Pos_E_min, ymax = Pos_E_max), colour = piRNA.colours["26"],
                  fill = piRNA.colours["26"],  alpha = 0.5, linewidth = NULL)+
      geom_line(data = dat26, aes(x = x, y = Pos_Mean), colour = piRNA.colours["26"], linewidth = 1)+
      geom_ribbon(data = dat26, aes(x = x, ymin = Neg_E_min, ymax = Neg_E_max), colour = piRNA.colours["26"],
                  fill = piRNA.colours["26"], alpha = 0.5, linewidth = NULL)+
      geom_line(data = dat26, aes(x = x, y = Neg_Mean), colour = piRNA.colours["26"], linewidth = 1)+
      
      geom_ribbon(data = dat27, aes(x = x, ymin = Pos_E_min, ymax = Pos_E_max), colour = piRNA.colours["27"],
                  fill = piRNA.colours["27"],  alpha = 0.5, linewidth = NULL)+
      geom_line(data = dat27, aes(x = x, y = Pos_Mean), colour = piRNA.colours["27"], linewidth = 1)+
      geom_ribbon(data = dat27, aes(x = x, ymin = Neg_E_min, ymax = Neg_E_max), colour = piRNA.colours["27"],
                  fill = piRNA.colours["27"], alpha = 0.5, linewidth = NULL)+
      geom_line(data = dat27, aes(x = x, y = Neg_Mean), colour = piRNA.colours["27"], linewidth = 1)+
      
      geom_ribbon(data = dat28, aes(x = x, ymin = Pos_E_min, ymax = Pos_E_max), colour = piRNA.colours["28"],
                  fill = piRNA.colours["28"],  alpha = 0.5, linewidth = NULL)+
      geom_line(data = dat28, aes(x = x, y = Pos_Mean), colour = piRNA.colours["28"], linewidth = 1)+
      geom_ribbon(data = dat28, aes(x = x, ymin = Neg_E_min, ymax = Neg_E_max), colour = piRNA.colours["28"],
                  fill = piRNA.colours["28"], alpha = 0.5, linewidth = NULL)+
      geom_line(data = dat28, aes(x = x, y = Neg_Mean), colour = piRNA.colours["28"], linewidth = 1)+
      
      geom_ribbon(data = dat29, aes(x = x, ymin = Pos_E_min, ymax = Pos_E_max), colour = piRNA.colours["29"],
                  fill = piRNA.colours["29"],  alpha = 0.5, linewidth = NULL)+
      geom_line(data = dat29, aes(x = x, y = Pos_Mean), colour = piRNA.colours["29"], linewidth = 1)+
      geom_ribbon(data = dat29, aes(x = x, ymin = Neg_E_min, ymax = Neg_E_max), colour = piRNA.colours["29"],
                  fill = piRNA.colours["29"], alpha = 0.5, linewidth = NULL)+
      geom_line(data = dat29, aes(x = x, y = Neg_Mean), colour = piRNA.colours["29"], linewidth = 1)+
      
      
      geom_hline(yintercept = 0, linetype = "solid", colour = "grey")+
      ggtitle(paste0(filename))+
      ylim(-max(piNormPlotCount)*1.1,max(piNormPlotCount)*1.1)+
      xlab ("nt position")+
      ylab ("Normalised Coverage")+
      piMaker_theme+
      guides(guide_legend)+
      theme(legend.position="top")+
      scale_colour_manual("Size", values = c(piRNA.colours))
    plot(gg)
    saveImage(paste(filename))
  }        
}



for(i in 1:length(piTally_List_Matrix)){
  namfile <- names(piTally_List_Matrix[i])
  dat <- piTally_List_Matrix[[i]]
  nam <- paste0(namfile, "_piSummary")
  dat_2429 <- dat
  
  #split data into the component parts for the plots
  
  dat_Mean <- (dat_2429[grepl("Mean", colnames(dat_2429))])
  dat_SD <- (dat_2429[grepl("SD", colnames(dat_2429))])
  dat_SD$x <- rownames(dat_SD)
  
  for(n in 24:29){
    dat_Mean[paste0("Pos_",n, "_x")] <- as.numeric(ifelse( (dat_Mean[(paste0("Pos_",n, "_Mean"))] >0), 
                                                           dat_Mean$Pos_pos_Mean, NA ))
    dat_Mean[paste0("Pos_",n, "_xend")] <- dat_Mean[paste0("Pos_",n, "_x")] + n
    dat_Mean[paste0("Neg_",n, "_x")] <- as.numeric(ifelse((dat_Mean[(paste0("Neg_",n, "_Mean"))] >0), 
                                                          dat_Mean$Neg_pos_Mean, NA))
    dat_Mean[paste0("Neg_",n, "_xend")] <- dat_Mean[paste0("Neg_",n, "_x")] + n
  }
  
  gg <- ggplot()+
    
    geom_hline(yintercept = 0, linetype = "solid", colour = "grey", linewidth = 1)+
    
    geom_ribbon(data = dat_SD, aes(x = 1:nrow(dat_SD), ymin = Pos_24_SD_Minus, ymax = Pos_24_SD_Plus), colour = piRNA.colours["24"],
                fill = piRNA.colours["24"], alpha = 0.5, linewidth = 1)+
    

    #geom_rect(data = dat_Mean, aes(ymin = 0, ymax = - Neg_24_Mean, xmin = Neg_24_x,
    #                               xmax = Neg_24_xend, fill = "24"), colour = piRNA.colours["24"], alpha = 0.2)+ 
    #geom_rect( aes(ymin = dat_SD$Neg_24_SD_Minus, ymax = dat_SD$Neg_24_SD_Plus, xmin = dat_Mean$Neg_24_x,
    #               xmax = dat_Mean$Neg_24_xend, fill = "24"), colour = piRNA.colours["24"], alpha = 0.4)+
  
    

    ggtitle(paste(nam))+
    piMaker_theme+
    scale_fill_manual("Size", values = c(piRNA.colours))
  
  plot(gg)  
  
  geom_rect( aes(ymin = 0, ymax = dat_Mean$Pos_24_Mean, 
                 xmin = dat_Mean$Pos_24_x,
                 xmax = dat_Mean$Pos_24_xend, fill = "24"), 
             colour = piRNA.colours["24"], alpha = 0.4)+
    
    geom_rect( aes(ymin = dat_SD$Pos_24_SD_Minus, ymax = dat_SD$Pos_24_SD_Plus, 
                   xmin = dat_Mean$Pos_24_x, 
                   xmax = dat_Mean$Pos_24_xend), 
               fill = NA,linewidth = 0.5, colour = piRNA.colours["24"], alpha = 0.2)+
    ylim(-max(piMax),max(piMax))+
  
    
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
  saveImage(paste0(nam))
}    
