# plot data from NET-Seq experiments

# get date
today <- Sys.Date()

# plotting function
locus_plot <-
  function(sample_input,
           sample_IP = c(" ", 1, 0, 0, 0),
           libname,
           date=Sys.Date()) {
    # construct title and filename
    main_title <- paste(libname, sample_input$name[1])
    filename <- paste(libname,sample_input$name[1], date,".pdf", sep="")
    
    # axes_limits
    max1 <- max(sample_input$sense)
    max2 <- max(sample_IP$sense)
    ymax <- max(max1, max2)
    
    min1 <- min(sample_input$antisense)
    min2 <- min(sample_IP$antisense)
    ymin <- min(min1, min2)
    
    # plot values
    plot(
      sample_input$position,
      # use input in case the ip sample is undefined
      sample_IP$sense,
      type = "l",
      col = "black",
      ylim = c(ymin, ymax),
      main = main_title,
      xlab = "",
      ylab = ""
    )
    # add axis labels using mtext
    mtext(text="position (bp)", side=1, line=2)
    mtext(text="reads (ppm)", side=2, line=2)
    lines(sample_input$position, sample_IP$antisense, col = "red")
    lines(sample_input$position, sample_input$sense, col = "gray")
    lines(sample_input$position, sample_input$antisense, col = "orange")
    # add legend
    legend(
      "topright",
      inset=0.02,
      legend = c("PRO-Seq sense", "PRO-Seq as", "PRO-cap sense", "PRO-cap as"),
      col = c("gray", "orange", "black", "red"),
      lty = 1,
      cex = 0.7,
      box.lty = 0
    )
    dev.copy(pdf,filename, width=8, height=4)
    dev.off()
  }


# load in data
# Proseq
sevenSK_Proseq_N <- read.delim("binned_loci_ppm_proseq/pro-seq_N_mapping_loci_200920_7SK_RNA.binppm", header = FALSE, skip=1)
colnames(sevenSK_Proseq_N) <- c("name", "position", "sense", "antisense", "total")
buffer <- -(sevenSK_Proseq_N$antisense)
for (i in 1:(length(sevenSK_Proseq_N$position)-43)){ 
  sevenSK_Proseq_N$antisense[i+43] <- -(sevenSK_Proseq_N$sense[i])
}
sevenSK_Proseq_N$antisense[1:43] <- 0

for (i in 1:(length(sevenSK_Proseq_N$position)-43)){ 
  sevenSK_Proseq_N$sense[i+43] <- buffer[i]
}
sevenSK_Proseq_N$sense[1:43] <- 0


sevenSK_Procap_2 <- read.delim("binned_loci_ppm_proseq/pro-cap_2_mapping_loci_200920_7SK_RNA.binppm", header = FALSE, skip=1)
colnames(sevenSK_Procap_2) <- c("name", "position", "sense", "antisense", "total")

for (i in 1:(length(sevenSK_Procap_2$position)-43)){ 
  sevenSK_Procap_2$antisense[i+43] <- sevenSK_Procap_2$antisense[i]
}
sevenSK_Procap_2$antisense[0:43] <- 0 
for (i in 43:(length(sevenSK_Procap_2$position))){ 
  sevenSK_Procap_2$sense[i-43] <- sevenSK_Procap_2$sense[i]
}
sevenSK_Procap_2$sense[(length(sevenSK_Procap_2$sense)-43):length(sevenSK_Procap_2$sense)] <- 0 

# plot values
locus_plot(sevenSK_Proseq_N,sevenSK_Proseq_N, "PRO-Seq_data_only",today)



