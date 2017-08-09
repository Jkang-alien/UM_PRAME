library(ssGSEA4TCGA)
library(ConsensusClusterPlus)
library(NMF)
library(Cairo)
library(survival)
library(rms)

rdate <- getFirehoseRunningDates(last = NULL)

dset <- getFirehoseDatasets()
dset <- 'UVM'
gs <- gs_gmt('h.all.v6.0.symbols.gmt')

data <- pancancer_ssGSEA(dset, rdate[1], gs)
data_t <- scale(t(data[[1]]), center = TRUE, scale = TRUE)

colnames(data_t)
colnames(data_t) <- gsub('http://www.broadinstitute.org/gsea/msigdb/cards/HALLMARK_', '', colnames(data_t))

results_col = ConsensusClusterPlus(data_t,maxK=6,reps=5000,pItem=0.8,pFeature=1,
                                   title='consensus_col',
                                   clusterAlg="hc",
                                   innerLinkage = "ward.D2",
                                   finalLinkage = "ward.D2",
                                   distance="euclidean",
                                   plot="pdf")

results_row = ConsensusClusterPlus(data[[1]],maxK=6,reps=5000,pItem=0.8,pFeature=1,
                                   title='consensus_row',
                                   clusterAlg="hc",
                                   innerLinkage = "ward.D2",
                                   finalLinkage = "ward.D2",
                                   distance="euclidean",
                                   plot="pdf")

ann_col <- data.frame(immune_class = factor(results_col[[3]]$consensusClass))
ann <- data.frame(class = as.factor( results_row[[3]]$consensusClass))
rownames(ann) <- gsub('-', '\\.', gsub('-...-...-....-..','',rownames(ann)))

data <- read.delim('./gdac.broadinstitute.org_UVM-TP.Aggregate_AnalysisFeatures.Level_4.2016012800.0.0/UVM-TP.transferedsamplefeatures.txt')

colnames(data)
data[1:5,1:5]
data_gdac <- t(data[,c(-1,-82, -83)])
dim(data_gdac)
colnames(data_gdac) <- data[,1]
rownames(data_gdac) <- colnames(data)[c(-1,-82, -83)]
data <- data.frame(data_gdac)
colnames(data)
data[, c(1:4, 15:251)] <- sapply(data[, c(1:4, 15:251)], as.character)

data[, c(1:4, 15:251)] <- sapply(data[, c(1:4, 15:251)], as.numeric)
data$ID <- rownames(data)

summary(data)
colnames(data)

ID_3pDel <- data$ID[data$CN_3p_Del <= -0.2]
ID_3qDel <- data$ID[data$CN_3q_Del <= -0.2]
ID_3pNor <- data$ID[data$CN_3p_Del > -0.2]
ID_3qNor <- data$ID[data$CN_3q_Del > -0.2]

sum(ID_3qNor == ID_3pNor)
sum(ID_3qDel == ID_3pDel)

Chr3 <- rownames(ann) %in% ID_3qNor
Chr3 <- factor(Chr3, 
               levels = c(TRUE, FALSE),
               labels = c('Disomy', 'Monosomy'))

ann$Chr3 <- Chr3

CairoPDF(file = 'cluster.pdf',
         width =7.5, height = 7.5, pointsize = 16)
aheatmap(data_t,
         hclustfun=function(d) hclust(dist(d, method = 'euclidean'), method = "ward.D2"),
         annRow = ann,
         annCol = ann_col,
         Colv = results_col[[3]]$consensusTree,
         Rowv = results_row[[3]]$consensusTree,
         labRow = rep('',dim(data_t)[1]))
dev.off()

sum(is.na(data_t))
data <- data.frame(data, ann[match(rownames(ann),data$ID),])
surv_months <- pmax(data$CLI_days_to_death, 
                    data$CLI_days_to_last_followup,
                    na.rm = TRUE)/30.4
data$surv_months <- surv_months

diff = survdiff(Surv(surv_months, CLI_vital_status == 1)~ class, 
                data = data)
diff

diff = survdiff(Surv(surv_months, CLI_vital_status == 1)~ class, 
                data = subset(data, Chr3 == 'Disomy'))
diff

diff = survdiff(Surv(surv_months, CLI_vital_status == 1)~ class, 
                data = subset(data, Chr3 == 'Monosomy'))
diff


fit = npsurv(Surv(surv_months, CLI_vital_status == 1)~ class, 
             data = subset(data, Chr3 == 'Monosomy'))

fit = npsurv(Surv(surv_months, CLI_vital_status == 1)~ class, 
             data = data)


strata = levels(data$class)

CairoPDF(file = 'survival_monosomy.pdf',
         width =7.5, height = 7.5, pointsize = 16)
survplot(fit,
         time.inc = 12,
         xlab = 'Months',
         lty = c(1:5),
         conf="none", add=FALSE, 
         label.curves=FALSE, abbrev.label=FALSE,
         levels.only=TRUE, lwd=par('lwd'),
         col=data$class,
         col.fill=gray(seq(.95, .75, length=5)),
         loglog=FALSE,n.risk=TRUE,logt=FALSE,
         dots=FALSE,
         grid=FALSE,
         srt.n.risk=0, sep.n.risk=0.04, adj.n.risk=0.5, 
         y.n.risk=-0.3, cex.n.risk=0.6, pr=FALSE       
)

legend(60, 1.0, strata, lty = c(1:3), cex = 0.8,
       xjust = 0, yjust = 1, x.intersp = 1, y.intersp = 1,
       trace = TRUE,
       bty = 'n')

legend(5.5*10, 0.8, 'P-value: 0.62', cex = 0.8,
       xjust = 0, yjust = 1, x.intersp = 1, y.intersp = 1,
       trace = TRUE,
       bty = 'n')
dev.off()

fit = npsurv(Surv(surv_months, CLI_vital_status == 1)~ Chr3, 
             data = data)

strata = levels(data$Chr3)

survplot(fit,
         time.inc = 12,
         xlab = 'Months',
         lty = c(1:2),
         conf="none", add=FALSE, 
         label.curves=FALSE, abbrev.label=FALSE,
         levels.only=TRUE, lwd=par('lwd'),
         col=data$class,
         col.fill=gray(seq(.95, .75, length=5)),
         loglog=FALSE,n.risk=TRUE,logt=FALSE,
         dots=FALSE,
         grid=FALSE,
         srt.n.risk=0, sep.n.risk=0.04, adj.n.risk=0.5, 
         y.n.risk=-0.3, cex.n.risk=0.6, pr=FALSE       
)


boxplot(data)

colnames(data)

data <- merge (data, 
               data.frame(ID = gsub('-', '\\.', gsub('-...-...-....-..','',rownames(data_t))), data_t), 
               by = 'ID')


pca <- prcomp(data[,256:283],
              center = FALSE,
              scale. = FALSE) 
predict(pca, 
        newdata=tail(data[,256:283], 2))

library(devtools)
install_github("ggbiplot", "vqv")

library(ggbiplot)
g <- ggbiplot(pca, obs.scale = 1, var.scale = 1, 
              groups = data$class, ellipse = TRUE, 
              circle = TRUE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')


print(g)

g <- ggbiplot(pca, obs.scale = 1, var.scale = 1, 
              groups = data$Chr3, ellipse = TRUE, 
              circle = TRUE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')


print(g)

rna <- read.table('20160128-UVM-RNAseq2GeneNorm.txt', sep="\t", dec=".", skip = 2)
con <- file("20160128-UVM-RNAseq2GeneNorm.txt","r")
first_line <- readLines(con,n=1)
close(con)
ID <- unlist(strsplit(first_line, '\t'))[-1]
PRAME <- t(rna[grep('PRAME\\|', rna[,1]),][-1])

PRAME <- c(t(as.numeric(as.character(rna[grep('PRAME\\|', rna[,1]),-1]))))
hist(log(PRAME))
