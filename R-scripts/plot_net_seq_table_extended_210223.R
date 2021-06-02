# read and plot Netseq data
net_seq_table <- read.csv2("Net_Seq_table_190221.csv", stringsAsFactors=T , dec = ",")
# net_seq_table <- relevel(net_seq_table$gene,"Act5C")
# net_seq_table <- relevel(net_seq_table$gene,"CG15099")
# net_seq_table <- relevel(net_seq_table$gene,"TCTP")
# net_seq_table <- relevel(net_seq_table$gene,"CG15098")
net_seq_table$gene <- factor(net_seq_table$gene, levels(net_seq_table$gene)[c(2,4,3,1)])


input_samples <- subset.data.frame(net_seq_table, net_seq_table$sample_type=="input")
input_samples_CGcut <- subset.data.frame(input_samples, input_samples$cut=="CG15098")
input_samples_CGcut_sense <- subset.data.frame(input_samples_CGcut, input_samples_CGcut$orientation=="sense")
input_samples_CGcut_antisense <- subset.data.frame(input_samples_CGcut, input_samples_CGcut$orientation=="antisense")

input_samples_TCTPcut <- subset.data.frame(input_samples, input_samples$cut=="TCTP")
input_samples_TCTPcut_sense <- subset.data.frame(input_samples_TCTPcut, input_samples_TCTPcut$orientation=="sense")
input_samples_TCTPcut_antisense <- subset.data.frame(input_samples_TCTPcut, input_samples_TCTPcut$orientation=="antisense")




pol2_samples <- subset.data.frame(net_seq_table, net_seq_table$sample_type=="NET_pol_2")
pol2_samples_CGcut <- subset(pol2_samples, pol2_samples$cut=="CG15098")
pol2_samples_CGcut_sense <- subset.data.frame(pol2_samples_CGcut, pol2_samples_CGcut$orientation=="sense")
pol2_samples_CGcut_antisense <- subset.data.frame(pol2_samples_CGcut, pol2_samples_CGcut$orientation=="antisense")

pol2_samples_TCTPcut <- subset(pol2_samples, pol2_samples$cut=="TCTP")
pol2_samples_TCTPcut_sense <- subset.data.frame(pol2_samples_TCTPcut, pol2_samples_TCTPcut$orientation=="sense")
pol2_samples_TCTPcut_antisense <- subset.data.frame(pol2_samples_TCTPcut, pol2_samples_TCTPcut$orientation=="antisense")

pol3_samples <- subset.data.frame(net_seq_table, net_seq_table$sample_type=="Net_pol_3")
pol3_samples_CGcut <- subset(pol3_samples, pol3_samples$cut=="CG15098")
pol3_samples_CGcut_sense <- subset.data.frame(pol3_samples_CGcut, pol3_samples_CGcut$orientation=="sense")
pol3_samples_CGcut_antisense <- subset.data.frame(pol3_samples_CGcut, pol3_samples_CGcut$orientation=="antisense")

pol3_samples_TCTPcut <- subset(pol3_samples, pol3_samples$cut=="TCTP")
pol3_samples_TCTPcut_sense <- subset.data.frame(pol3_samples_TCTPcut, pol3_samples_TCTPcut$orientation=="sense")
pol3_samples_TCTPcut_antisense <- subset.data.frame(pol3_samples_TCTPcut, pol3_samples_TCTPcut$orientation=="antisense")

# par(mfrow=c(3,4))
# plot(input_samples_CGcut_antisense$reads_ppm~input_samples_CGcut_antisense$gene, xlab="", ylab="reads (ppm)", main="input antisense", ylim=c(0,2), las=2 )
# points(input_samples_CGcut_antisense$reads_ppm~input_samples_CGcut_antisense$gene)
# mtext("CG15098 cut", cex=0.7)
# plot(input_samples_TCTPcut_antisense$reads_ppm~input_samples_TCTPcut_antisense$gene, xlab="", ylab="",main="input antisense", ylim=c(0,2), las=2)
# points(input_samples_TCTPcut_antisense$reads_ppm~input_samples_TCTPcut_antisense$gene)
# mtext("TCTP cut", cex=0.7)
# 
# plot(input_samples_CGcut_sense$reads_ppm~input_samples_CGcut_sense$gene, xlab="", ylab="", main="input sense", ylim=c(0,150), las=2)
# points(input_samples_CGcut_sense$reads_ppm~input_samples_CGcut_sense$gene)
# mtext("CG15098 cut", cex=0.7)
# plot(input_samples_TCTPcut_sense$reads_ppm~input_samples_TCTPcut_sense$gene, xlab="", ylab="", main="input sense", ylim=c(0,150), las=2)
# points(input_samples_TCTPcut_sense$reads_ppm~input_samples_TCTPcut_sense$gene)
# mtext("TCTP cut", cex=0.7)
# 
# plot(pol2_samples_CGcut_antisense$reads_ppm~pol2_samples_CGcut_antisense$gene, xlab="", ylab="reads (ppm)", main="NET-Seq pol2 antisense", ylim=c(0,12), las=2)
# points(pol2_samples_CGcut_antisense$reads_ppm~pol2_samples_CGcut_antisense$gene)
# mtext("CG15098 cut", cex=0.7)
# plot(pol2_samples_TCTPcut_antisense$reads_ppm~pol2_samples_TCTPcut_antisense$gene, xlab="", ylab="", main="NET-Seq pol2 antisense",  ylim=c(0,12), las=2)
# points(pol2_samples_TCTPcut_antisense$reads_ppm~pol2_samples_TCTPcut_antisense$gene)
# mtext("TCTP cut", cex=0.7)
# 
# plot(pol2_samples_CGcut_sense$reads_ppm~pol2_samples_CGcut_sense$gene, xlab="", ylab="",main="NET-Seq pol2 sense", ylim=c(0,400), las=2)
# points(pol2_samples_CGcut_sense$reads_ppm~pol2_samples_CGcut_sense$gene)
# mtext("CG15098 cut", cex=0.7)
# plot(pol2_samples_TCTPcut_sense$reads_ppm~pol2_samples_TCTPcut_sense$gene,xlab="", ylab="",main="NET-Seq pol2 sense", ylim=c(0,400), las=2)
# points(pol2_samples_TCTPcut_sense$reads_ppm~pol2_samples_TCTPcut_sense$gene)
# mtext("TCTP cut", cex=0.7)
# 
# plot(pol3_samples_CGcut_antisense$reads_ppm~pol3_samples_CGcut_antisense$gene, xlab="", ylab="reads (ppm)", main="NET-Seq pol3 antisense", ylim=c(0,12), las=2)
# points(pol3_samples_CGcut_antisense$reads_ppm~pol3_samples_CGcut_antisense$gene)
# mtext("CG15098 cut", cex=0.7)
# plot(pol3_samples_TCTPcut_antisense$reads_ppm~pol3_samples_TCTPcut_antisense$gene,xlab="", ylab="", main="NET-Seq pol3 antisense", ylim=c(0,12), las=2)
# points(pol3_samples_TCTPcut_antisense$reads_ppm~pol3_samples_TCTPcut_antisense$gene)
# mtext("TCTP cut", cex=0.7)
# 
# plot(pol3_samples_CGcut_sense$reads_ppm~pol3_samples_CGcut_sense$gene, xlab="", ylab="", main="NET-Seq pol3 sense", ylim=c(0,400), las=2)
# points(pol3_samples_CGcut_sense$reads_ppm~pol3_samples_CGcut_sense$gene)
# mtext("CG15098 cut", cex=0.7)
# plot(pol3_samples_TCTPcut_sense$reads_ppm~pol3_samples_TCTPcut_sense$gene,xlab="", ylab="", main="NET-Seq pol3 sense", ylim=c(0,400), las=2)
# points(pol3_samples_TCTPcut_sense$reads_ppm~pol3_samples_TCTPcut_sense$gene)
# mtext("TCTP cut", cex=0.7)

# dev.copy(pdf,"Plot_NET_Seq_table_210219.pdf", width=9, height=6)
# dev.off()
pol2_samples_CGcut_ratios <- pol2_samples_CGcut_sense
pol2_samples_CGcut_ratios$reads_ppm <- pol2_samples_CGcut_antisense$reads_ppm/(pol2_samples_CGcut_antisense$reads_ppm+pol2_samples_CGcut_sense$reads_ppm)
pol2_samples_TCTPcut_ratios <- pol2_samples_TCTPcut_sense
pol2_samples_TCTPcut_ratios$reads_ppm <- pol2_samples_TCTPcut_antisense$reads_ppm/(pol2_samples_TCTPcut_antisense$reads_ppm+pol2_samples_TCTPcut_sense$reads_ppm)

#par(mfrow=c(1,2))
plot(pol2_samples_CGcut_ratios$reads_ppm~pol2_samples_CGcut_ratios$gene, xlab="", ylab="",main="NET-Seq pol2 ratios", ylim=c(0,0.2), las=2)
points(pol2_samples_CGcut_ratios$reads_ppm~pol2_samples_CGcut_ratios$gene)
mtext("CG15098 cut")
plot(pol2_samples_TCTPcut_ratios$reads_ppm~pol2_samples_TCTPcut_ratios$gene,xlab="", ylab="",main="NET-Seq pol2 ratios", ylim=c(0,0.2), las=2)
points(pol2_samples_TCTPcut_ratios$reads_ppm~pol2_samples_TCTPcut_ratios$gene)
mtext("TCTP cut")