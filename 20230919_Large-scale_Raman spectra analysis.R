# Set working directory to the folder containing spectra in .txt
setwd("/Users/putingdong/Documents/2023_Fn tsRNA manuscripts/20230830_Fn galk and HMKO_Raman spectra")
filelist <- list.files(pattern = ".*.txt")
datalist <- lapply(filelist, function(x)read.table(x, header=T, sep = "\t")) 
# the name of certain columns
filelist <- as.data.frame(filelist)

# Create a data frame containing all spectra
combinedata <- do.call("cbind", datalist)
Ramanshift <- combinedata[,1]

# exact the Raman intensity data from the combined data by removing the odd columns
row_odd <- seq_len(ncol(combinedata)) %% 2
spectra <- combinedata[, row_odd == 0]

# annotate the filename to individual column by using the names function
names(spectra) <- t(filelist)
# plot the original plots without baseline correction
matplot(Ramanshift, spectra, type = "l")

# remove baseline
spectra <- as.matrix(t(spectra))
library(baseline)
spectra.bl<- baseline.modpolyfit(spectra, degree = 5, tol = 0.01, rep=100)

# extract the spectral information after baseline correction
spectra.corr <- t(spectra.bl[["corrected"]])

# plot the spectra after baseline correction
matplot(Ramanshift, spectra.corr, type = "l")

# annotate the filename to individual column by using the names function
spectra.corr <- as.data.frame(spectra.corr)
names(spectra.corr) <- t(filelist)
spectra.corr <- t(spectra.corr)

# remove cosmic ray
library(hyperSpec.utils)
library(hyperSpec)
spectra.corr <- as.hyperSpec(spectra.corr)
spectra.corr.crr <- crr(spectra.corr, threshold = 20)

# plot the spectra after cosmic ray removal
spectra.corr.crr.data <- spectra.corr.crr@data$spc
matplot(Ramanshift, t(spectra.corr.crr.data), type = "l")

# normalization (through vector normalization)
X <- t(spectra.corr.crr.data)
final.dataset <- scale(X, center = T, scale = T)
matplot(Ramanshift, final.dataset, type = "l")


# calculate the mean value from each group
log.phase.Fn.galk <- final.dataset[,1:52]
mean.lp.Fn.galk <- rowMeans(log.phase.Fn.galk)

log.phase.Fn.hmko <- final.dataset[,53:103]
mean.lp.Fn.hmko <- rowMeans(log.phase.Fn.hmko)

stationary.phase.Fn.galk <- final.dataset[,104:154] #stationary-phase Fn galk ko has apparently higher amount of RNA compared to stationary-phase Fn hmko
mean.sp.Fn.galk <- rowMeans(stationary.phase.Fn.galk)

stationary.phase.Fn.hmko <- final.dataset[,155:213] 
mean.sp.Fn.hmko <- rowMeans(stationary.phase.Fn.hmko)

mean.matrix <-cbind(mean.lp.Fn.galk, mean.lp.Fn.hmko, mean.sp.Fn.galk, mean.sp.Fn.hmko)
matplot(Ramanshift, mean.matrix, type = "l")

#quantification of the amount of DNA/RNA from both log-phase and stationary-phase F. nucleatum
sp.Fn.galk.rna <- colSums(stationary.phase.Fn.galk[212:223,]) #from wavenumber 770 cm-1 to 789 cm-1
sp.Fn.hmko.rna <- colSums(stationary.phase.Fn.hmko[212:223,])
lp.Fn.galk.rna <- colSums(log.phase.Fn.galk[212:223,])
lp.Fn.hmko.rna <- colSums(log.phase.Fn.hmko[212:223,])

# calculate the p-value from the above groups
t.test(sp.Fn.galk.rna, sp.Fn.hmko.rna) # a p value (5.414e-08) smaller than 0.0001 was obtained.
t.test(lp.Fn.galk.rna, lp.Fn.hmko.rna) # a p value (0.05269) was obtained, not significant difference.

# plot the averaged spectra
matplot(Ramanshift, mean.matrix[,3:4], type = "l",  
        xlab = "Raman shift", ylab = "Averaged Raman Intensity", xlim= c(400, 1800))
legend('topleft',inset=0.02,c("Stationary-phase Fn galk ko","Stationary-phase Fn hmko")
       ,lty=c(1,2),col=c("black","red"))
grid(nx = NULL, ny = 5, col = "lightgray", lty = "dotted")
abline(v=780, lwd = 1, lty = "dashed", col = "blue")

# Quantification of nucleic acid amount based on the integrated Raman intensity at 780 cm-1
row1 <- rep(c(1), 52)
row2 <- rep(c(2), 51)
data <- data.frame (name = c(row1, row2), value = c(lp.Fn.galk.rna, lp.Fn.hmko.rna)) #has to be data.frame format
data$name<- as.factor(data$name) 
library(ggplot2)
library(tidyverse)
library(hrbrthemes)
library(viridis)
ggplot(data, aes(x = name, y = value, fill = name)) +
  geom_boxplot(alpha = 0.80, size = 0.5) +
  geom_point(aes(colour = name), size = 5, shape = 21, position = position_jitterdodge()) +
  ylim(-5, 6)+ xlim(-2,3)+ylab("Integrated Raman Intensity")+
  theme(text = element_text(size = 24),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(), 
        legend.position = "none") + 
  scale_x_discrete(labels=c("Log-phase Fn galk ko", "Log-phase Fn hmko"))+
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1))+
  theme(panel.grid.major = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))
# ggsave("Quantification of RNA amount from log-phase Fn galk and hmko.png", dpi = 1200)

# principal component analysis
final.dataset <- as.data.frame(final.dataset)
# regroup all the dataset
group1 <- c(rep("log-phase Fn galk", 52))
group2 <- c(rep("log-phase Fn hmko", 51))
group3 <- c(rep("stationary-phase Fn galk", 51))
group4 <- c(rep("stationary-phase Fn hmko", 59))
data.pca <- t(final.dataset)
data.pca <- as.data.frame(data.pca)
group <- c(group1, group2, group3, group4)
data.pca <- cbind(data.pca, group)

library("factoextra")
# stationary-phase and log-phase matrix
data.pca.stationary <- data.pca[104:213,]
data.pca.log <- data.pca[1:103,]

# log-phase pca plot
final.dataset.pca.log <- prcomp(data.pca.log[, -856], scale = T)
fviz_pca_ind(final.dataset.pca.log, pointsize=2, geom="point", habillage = data.pca.log$group, 
             addEllipses = TRUE, ellipse.level = 0.68) + theme_minimal()

# stationary-phase pca plot
final.dataset.pca.stationary <- prcomp(data.pca.stationary[, -856], scale = T)
fviz_pca_ind(final.dataset.pca.stationary, pointsize=2, geom="point", habillage = data.pca.stationary$group, 
             addEllipses = TRUE, ellipse.level = 0.68) + theme_minimal()

