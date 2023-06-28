
MatrixPlot <- function(Linedata, Bardata, ProbMax, ProbWieghtedMax, BarMax, z_Max = 5){

  linedata <- Linedata
  bardata <- Bardata
  probmax <- ProbMax
  probmax2 <- ProbWieghtedMax
  barmax <- BarMax
  zmax <- z_Max

gz <- ggplot(data = linedata, aes(x = x))+
  geom_line(aes(y = Pos_Z, colour = "Pos"), linewidth = 0.75 )+
  geom_line(aes(y = Neg_Z, colour = "Neg"), linewidth = 0.75 )+
  geom_line(aes(y = Z_Score, colour = "Overall"),  linewidth = 1 )+
  ggtitle("Z-score")+
  ylim(-zmax,zmax)+
  ylab("Z-Score")+
  xlab("Overlap (nt)")+
  theme(legend.position="none")+
  guides(colour = "none")+
  scale_colour_manual(values=group.colours)+
  #theme(aspect.ratio = 0.25:1)+
  piMaker_theme
#plot(gz)

go <- ggplot(data = bardata, aes(x = x))+
  geom_bar(stat = "identity", aes(y = Count, fill = Group, alpha = 0.5), colour = "darkslategrey",  linewidth = 0.75)+
  ggtitle("Overlapping pairs")+
  ylim(0, max(barmax)*1.1)+
  xlab("Overlap (nt)")+
  ylab("No. of pairs")+
  guides(colour = FALSE)+
  piMaker_theme+
  theme(legend.position="none")+
  #scale_colour_manual(values=group.colours)
  scale_fill_manual(values = group.colours)
#plot(go)

gp <- ggplot(data = linedata, aes(x = x))+
  geom_line(aes(y = Pos_probability, colour = "Pos"), linewidth = 0.75 )+
  geom_line(aes(y = Neg_probability, colour = "Neg"), linewidth = 0.75 )+
  geom_line(aes(y = Probability, colour = "Overall"), linewidth = 1 )+
  ggtitle("Probability")+
  ylim(0, max(probmax)*1.1)+
  ylab("Probability")+
  xlab("Overlap (nt)")+
  theme(legend.position="none")+
  guides(colour = FALSE)+
  #theme(aspect.ratio = 0.25:1)+
  piMaker_theme+
  scale_colour_manual( values = group.colours) 
#plot(gp)

gpw <- ggplot(data = linedata, aes(x = x))+
  geom_line(aes(y = Pos_weighted_probability, colour = "Pos"), linewidth = 0.75 )+
  geom_line(aes(y = Neg_weighted_probability, colour = "Neg"), linewidth = 0.75 )+
  geom_line(aes(y = Weighted_Probability, colour = "Overall"), linewidth = 1 )+
  ggtitle("Weighted probability")+
  ylim(0, max(probmax2)*1.1)+
  ylab("Weighted probability")+
  xlab("Overlap (nt)")+
  theme(legend.position="none")+
  guides(colour = FALSE)+
  #theme(aspect.ratio = 0.25:1)+
  piMaker_theme+
  scale_colour_manual( values = group.colours) 
#plot(gpw)

fig <- ggarrange(gz +rremove("xlab") +rremove("x.text"),  gp +rremove("xlab") +rremove("x.text"), go  ,gpw, nrow = 2, ncol = 2 )
annotate_figure(fig, top = text_grob(paste0(namd), 
                                     color = "darkslategrey", face = "bold", size = 14))
}

MatrixPlotSD <- function(Linedata, Bardata, ProbMax, ProbWieghtedMax, BarMax, z_Max = 5){
  
  linedata <- Linedata
  bardata <- Bardata
  probmax <- ProbMax
  probmax2 <- ProbWieghtedMax
  barmax <- BarMax
  zmax <- z_Max


gz <- ggplot(data = linedata, aes(x = x_Mean))+
  geom_line(aes(y = Pos_Z_Mean, colour = "Pos" ), linewidth = 0.75 )+
  geom_ribbon(aes(x = x_Mean, ymin = Pos_Z_SD_Minus, 
                  ymax = Pos_Z_SD_Plus, colour = "Pos" ), 
              fill = group.colours["Pos"], 
              alpha = 0.25, linewidth = 0)+
  geom_line(aes(y = Neg_Z_Mean, colour = "Neg"), linewidth = 0.75 )+
  geom_ribbon(aes(x = x_Mean, ymin = Neg_Z_SD_Minus, 
                  ymax = Neg_Z_SD_Plus, colour = "Neg" ), 
              fill = group.colours["Neg"], 
              alpha = 0.25, linewidth = 0)+
  geom_line(aes(y = Z_Score_Mean, colour = "Overall"),  linewidth = 1 )+
  geom_ribbon(aes(x = x_Mean, ymin = Z_Score_SD_Minus, 
                  ymax = Z_Score_SD_Plus, colour = "Overall" ), 
              fill = group.colours["Overall"], 
              alpha = 0.25, linewidth = 0)+
  ggtitle("Z-score")+
  ylim(-zmax,zmax)+
  ylab("Z-Score")+
  xlab("Overlap (nt)")+
  theme(legend.position="none")+
  guides(colour = "none")+
  scale_colour_manual(values=group.colours)+
  scale_fill_manual(values=group.colours)+
  #theme(aspect.ratio = 0.25:1)+
  piMaker_theme
#plot(gz)

go <- ggplot(data = Bardata, aes(x = x))+
  geom_bar(stat = "identity", aes(y = Count, fill = Group, alpha = 0.5), colour = "darkslategrey",  linewidth = 0.75)+
  geom_errorbar( aes(ymin = Cumulative - SD, ymax = Cumulative + SD), position = "identity", linewidth=0.5, colour="black", alpha=0.7 )+
  #geom_errorbar( aes(x = data$x, ymin = data$Neg_E_Min, ymax = data$Neg_E_Max), linewidth=1, colour="black", alpha=0.9 )+
  ggtitle("Overlapping pairs")+
  ylim(0, (max(barmax)*1.1))+
  xlab("Overlap (nt)")+
  ylab("No. of pairs")+
  guides(colour = FALSE)+
  piMaker_theme+
  theme(legend.position="none")+
  #scale_colour_manual(values=group.colours)
  scale_fill_manual(values = group.colours)
#plot(go)

gp <- ggplot(data = linedata, aes(x = x_Mean))+
  geom_line(aes(y = Pos_probability_Mean, colour = "Pos"), linewidth = 0.75 )+
  geom_ribbon(aes(x = x_Mean, ymin = Pos_probability_SD_Minus, 
                  ymax = Pos_probability_SD_Plus, colour = "Pos" ), 
              fill = group.colours["Pos"], 
              alpha = 0.25, linewidth = 0)+
  geom_line(aes(y = Neg_probability_Mean, colour = "Neg"), linewidth = 0.75 )+
  geom_ribbon(aes(x = x_Mean, ymin = Neg_probability_SD_Minus, 
                  ymax = Neg_probability_SD_Plus, colour = "Neg" ), 
              fill = group.colours["Neg"], 
              alpha = 0.25, linewidth = 0)+
  geom_line(aes(y = Probability_Mean, colour = "Overall"), linewidth = 1 )+
  geom_ribbon(aes(x = x_Mean, ymin = Probability_SD_Minus, 
                  ymax = Probability_SD_Plus, colour = "Overall" ), 
              fill = group.colours["Overall"], 
              alpha = 0.25, linewidth = 0)+
  ggtitle("Probability")+
  ylim(0, (max(probmax)*1.1))+
  ylab("Probability")+
  xlab("Overlap (nt)")+
  theme(legend.position="none")+
  guides(colour = FALSE)+
  #theme(aspect.ratio = 0.25:1)+
  piMaker_theme+
  scale_colour_manual( values = group.colours) 
#plot(gp)

gpw <- ggplot(data = linedata, aes(x = x_Mean))+
  geom_line(aes(y = Pos_weighted_probability_Mean, colour = "Pos"), linewidth = 0.75 )+
  geom_ribbon(aes(x = x_Mean, ymin = Pos_weighted_probability_SD_Minus, 
                  ymax = Pos_weighted_probability_SD_Plus, colour = "Pos" ), 
              fill = group.colours["Pos"], 
              alpha = 0.25, linewidth = 0)+
  geom_line(aes(y = Neg_weighted_probability_Mean, colour = "Neg"), linewidth = 0.75 )+
  geom_ribbon(aes(x = x_Mean, ymin = Neg_weighted_probability_SD_Minus, 
                  ymax = Neg_weighted_probability_SD_Plus, colour = "Pos" ), 
              fill = group.colours["Neg"], 
              alpha = 0.25, linewidth = 0)+
  geom_line(aes(y = Weighted_Probability_Mean, colour = "Overall"), linewidth = 1 )+
  geom_ribbon(aes(x = x_Mean, ymin = Weighted_Probability_SD_Minus, 
                  ymax = Weighted_Probability_SD_Plus, colour = "Overall" ), 
              fill = group.colours["Overall"], 
              alpha = 0.25, linewidth = 0)+
  ggtitle("Weighted probability")+
  ylim(0, (max(probmax2)*1.1))+
  ylab("Weighted probability")+
  xlab("Overlap (nt)")+
  theme(legend.position="none")+
  guides(colour = FALSE)+
  #theme(aspect.ratio = 0.25:1)+
  piMaker_theme+
  scale_colour_manual( values = group.colours) 
#plot(gpw)

fig <- ggarrange(gz +rremove("xlab") +rremove("x.text"),  gp +rremove("xlab") +rremove("x.text"), go  ,gpw, nrow = 2, ncol = 2, align = "hv" )
annotate_figure(fig, top = text_grob(paste0(namd), 
                                     color = "darkslategrey", face = "bold", size = 14))
}

piMapper <- function(x, Sizes = c(24:29) ){
  
  dat_2429. <- x
  
  for(n in 24:29){
    dat_2429.[paste0("Pos_",n, "_x")] <- as.numeric(ifelse( (dat_2429.[(paste0("Pos_",n))] >0), 
                                                           dat_2429.$Pos_pos, NA ))
    dat_2429.[paste0("Pos_",n, "_xend")] <- dat_2429.[paste0("Pos_",n, "_x")] + n
    dat_2429.[paste0("Neg_",n, "_x")] <- as.numeric(ifelse((dat_2429.[(paste0("Neg_",n))] >0), 
                                                          dat_2429.$Neg_pos, NA))
    dat_2429.[paste0("Neg_",n, "_xend")] <- dat_2429.[paste0("Neg_",n, "_x")] + n
  }
  
  gg <- ggplot(data = dat_2429.)+
  
    geom_hline(yintercept = 0, linetype = "solid", colour = "grey")+
    
    ylim(-max(piMaxInd),max(piMaxInd))+
    ggtitle(paste(nam))+
    piMaker_theme+
    scale_fill_manual("Size", values = c(piRNA.colours))+
    
    for (i in 1:length(Sizes)){
      n <- Sizes[[i]]
      
      PyMax <- paste0("Pos_", n)
      PyMin <- paste0("Pos_", n, "_x")
      PxMax <- paste0("Pos_", n, "_end")
      NyMax <- paste0("Neg_", n)
      NyMin <- paste0("Neg_", n, "_x")
      NxMax <- paste0("Neg_", n, "_end")
      
  gg <- gg + geom_rect(data = dat_2429., aes(ymin = 0, ymax = paste(PyMax), xmin = paste(PyMin),
                                      xmax = paste(PxMax), fill = paste(n)), colour = piRNA.colours[paste(n)], alpha = 0.4)
  gg <- gg + geom_rect(data = dat_2429., aes(ymin = 0, ymax = - paste(NyMax), xmin = paste(NyMin),
                                        xmax = paste(NxMax), fill = paste(n)), colour = piRNA.colours[paste(n)], alpha = 0.4)
      
    }
  plot(gg)
  
  ylim(-max(piMaxInd),max(piMaxInd))+
    ggtitle(paste(nam))+
    piMaker_theme+
    scale_fill_manual("Size", values = c(piRNA.colours))
    
    
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
    ylim(-max(piMaxInd),max(piMaxInd))+
    ggtitle(paste(nam))+
    piMaker_theme+
    scale_fill_manual("Size", values = c(piRNA.colours))
  
}