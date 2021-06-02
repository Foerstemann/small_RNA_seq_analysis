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
      legend = c("input sense", "input as", "IP sense", "IP as"),
      col = c("gray", "orange", "black", "red"),
      lty = 1,
      cex = 0.7,
      box.lty = 0
    )
    dev.copy(pdf,filename, width=8, height=4)
    dev.off()
  }


# load in data
# RB10
act5C_RB10_TCTP_input <- read.delim("binned_loci_ppm/RB10_TCTP_input_trim3_20_28_mapping_loci_200920_Act5C_genomic_region.binppm", header = FALSE, skip=1)
colnames(act5C_RB10_TCTP_input) <- c("name", "position", "sense", "antisense", "total")
act5C_RB10_TCTP_IP_pol2 <- read.delim("binned_loci_ppm/RB10_TCTP_IP-pol2_trim3_20_28_mapping_loci_200920_Act5C_genomic_region.binppm", header = FALSE, skip=1)
colnames(act5C_RB10_TCTP_IP_pol2) <- c("name", "position", "sense", "antisense", "total")

act5C_RB10_CG15098_input <- read.delim("binned_loci_ppm/RB10_CG15098_input_trim3_20_28_mapping_loci_200920_Act5C_genomic_region.binppm", header = FALSE, skip=1)
colnames(act5C_RB10_CG15098_input) <- c("name", "position", "sense", "antisense", "total")
act5C_RB10_CG15098_IP_pol2 <- read.delim("binned_loci_ppm/RB10_CG15098_IP-pol2_trim3_20_28_mapping_loci_200920_Act5C_genomic_region.binppm", header = FALSE, skip=1)
colnames(act5C_RB10_CG15098_IP_pol2) <- c("name", "position", "sense", "antisense", "total")

# RB12-13
act5C_RB12_13_TCTP_input_pol2 <- read.delim("binned_loci_ppm/RB12-13_TCTP_input-pol2_trim3_20_28_mapping_loci_200920_Act5C_genomic_region.binppm", header = FALSE, skip=1)
colnames(act5C_RB12_13_TCTP_input_pol2) <- c("name", "position", "sense", "antisense", "total")
act5C_RB12_13_TCTP_IP_pol2 <- read.delim("binned_loci_ppm/RB12-13_TCTP_IP-pol2_trim3_20_28_mapping_loci_200920_Act5C_genomic_region.binppm", header = FALSE, skip=1)
colnames(act5C_RB12_13_TCTP_IP_pol2) <- c("name", "position", "sense", "antisense", "total")

act5C_RB12_13_TCTP_input_pol3 <- read.delim("binned_loci_ppm/RB12-13_TCTP_input-pol3_trim3_20_28_mapping_loci_200920_Act5C_genomic_region.binppm", header = FALSE, skip=1)
colnames(act5C_RB12_13_TCTP_input_pol3) <- c("name", "position", "sense", "antisense", "total")
act5C_RB12_13_TCTP_IP_pol3 <- read.delim("binned_loci_ppm/RB12-13_TCTP_IP-pol3_trim3_20_28_mapping_loci_200920_Act5C_genomic_region.binppm", header = FALSE, skip=1)
colnames(act5C_RB12_13_TCTP_IP_pol3) <- c("name", "position", "sense", "antisense", "total")

act5C_RB12_13_CG15098_input_pol2 <- read.delim("binned_loci_ppm/RB12-13_CG15098_input-pol2_trim3_20_28_mapping_loci_200920_Act5C_genomic_region.binppm", header = FALSE, skip=1)
colnames(act5C_RB12_13_CG15098_input_pol2)<- c("name", "position", "sense", "antisense", "total")
act5C_RB12_13_CG15098_IP_pol2 <- read.delim("binned_loci_ppm/RB12-13_CG15098_IP-pol2_trim3_20_28_mapping_loci_200920_Act5C_genomic_region.binppm", header = FALSE, skip=1)
colnames(act5C_RB12_13_CG15098_IP_pol2) <- c("name", "position", "sense", "antisense", "total")

act5C_RB12_13_CG15098_input_pol3 <- read.delim("binned_loci_ppm/RB12-13_CG15098_input-pol3_trim3_20_28_mapping_loci_200920_Act5C_genomic_region.binppm", header = FALSE, skip=1)
colnames(act5C_RB12_13_CG15098_input_pol3) <- c("name", "position", "sense", "antisense", "total")
act5C_RB12_13_CG15098_IP_pol3 <- read.delim("binned_loci_ppm/RB12-13_CG15098_IP-pol3_trim3_20_28_mapping_loci_200920_Act5C_genomic_region.binppm", header = FALSE, skip=1)
colnames(act5C_RB12_13_CG15098_IP_pol3) <- c("name", "position", "sense", "antisense", "total")

# RB15-16
act5C_RB15_16_TCTP_input_pol2 <- read.delim("binned_loci_ppm/RB15-16_TCTP_input-pol2_trim3_20_28_mapping_loci_200920_Act5C_genomic_region.binppm", header = FALSE, skip=1)
colnames(act5C_RB15_16_TCTP_input_pol2) <- c("name", "position", "sense", "antisense", "total")
act5C_RB15_16_TCTP_IP_pol2 <- read.delim("binned_loci_ppm/RB15-16_TCTP_IP-pol2_trim3_20_28_mapping_loci_200920_Act5C_genomic_region.binppm", header = FALSE, skip=1)
colnames(act5C_RB15_16_TCTP_IP_pol2) <- c("name", "position", "sense", "antisense", "total")

act5C_RB15_16_TCTP_input_pol3 <- read.delim("binned_loci_ppm/RB15-16_TCTP_input-pol3_trim3_20_28_mapping_loci_200920_Act5C_genomic_region.binppm", header = FALSE, skip=1)
colnames(act5C_RB15_16_TCTP_input_pol3) <- c("name", "position", "sense", "antisense", "total")
act5C_RB15_16_TCTP_IP_pol3 <- read.delim("binned_loci_ppm/RB15-16_TCTP_IP-pol3_trim3_20_28_mapping_loci_200920_Act5C_genomic_region.binppm", header = FALSE, skip=1)
colnames(act5C_RB15_16_TCTP_IP_pol3) <- c("name", "position", "sense", "antisense", "total")

act5C_RB15_16_CG15098_input_pol2 <- read.delim("binned_loci_ppm/RB15-16_CG15098_input-pol2_trim3_20_28_mapping_loci_200920_Act5C_genomic_region.binppm", header = FALSE, skip=1)
colnames(act5C_RB15_16_CG15098_input_pol2)<- c("name", "position", "sense", "antisense", "total")
act5C_RB15_16_CG15098_IP_pol2 <- read.delim("binned_loci_ppm/RB15-16_CG15098_IP-pol2_trim3_20_28_mapping_loci_200920_Act5C_genomic_region.binppm", header = FALSE, skip=1)
colnames(act5C_RB15_16_CG15098_IP_pol2) <- c("name", "position", "sense", "antisense", "total")

act5C_RB15_16_CG15098_input_pol3 <- read.delim("binned_loci_ppm/RB15-16_CG15098_input-pol3_trim3_20_28_mapping_loci_200920_Act5C_genomic_region.binppm", header = FALSE, skip=1)
colnames(act5C_RB15_16_CG15098_input_pol3) <- c("name", "position", "sense", "antisense", "total")
act5C_RB15_16_CG15098_IP_pol3 <- read.delim("binned_loci_ppm/RB15-16_CG15098_IP-pol3_trim3_20_28_mapping_loci_200920_Act5C_genomic_region.binppm", header = FALSE, skip=1)
colnames(act5C_RB15_16_CG15098_IP_pol3) <- c("name", "position", "sense", "antisense", "total")


# plot values
locus_plot(act5C_RB10_TCTP_input,act5C_RB10_TCTP_IP_pol2,"RB10_TCTP-cut",today)
locus_plot(act5C_RB10_CG15098_input,act5C_RB10_CG15098_IP_pol2,"RB10_CG15098-cut",today)

locus_plot(act5C_RB12_13_TCTP_input_pol2,act5C_RB12_13_TCTP_IP_pol2,"RB12-13_TCTP-cut_pol2",today)
locus_plot(act5C_RB12_13_CG15098_input_pol2,act5C_RB12_13_CG15098_IP_pol2,"RB12-13_CG15098-cut_pol2",today)
locus_plot(act5C_RB12_13_TCTP_input_pol3,act5C_RB12_13_TCTP_IP_pol3,"RB12-13_TCTP-cut_pol3",today)
locus_plot(act5C_RB12_13_CG15098_input_pol3,act5C_RB12_13_CG15098_IP_pol3,"RB12-13_CG15098-cut_pol3",today)

locus_plot(act5C_RB15_16_TCTP_input_pol2,act5C_RB15_16_TCTP_IP_pol2,"RB15-16_TCTP-cut_pol2",today)
locus_plot(act5C_RB15_16_CG15098_input_pol2,act5C_RB15_16_CG15098_IP_pol2,"RB15-16_CG15098-cut_pol2",today)
locus_plot(act5C_RB15_16_TCTP_input_pol3,act5C_RB15_16_TCTP_IP_pol3,"RB15-16_TCTP-cut_pol3",today)
locus_plot(act5C_RB15_16_CG15098_input_pol3,act5C_RB15_16_CG15098_IP_pol3,"RB15-16_CG15098-cut_pol3",today)




