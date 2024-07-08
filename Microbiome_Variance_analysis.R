# =============================================================
# =============================================================
## The source coda for CEMC Microbiome Project in Southwest China ####
## Phenotypes associated with interindividual variation of gut microbiome ####
# =============================================================
# =============================================================

# >> load libraries ==================
install.packages("")

library(ggplot2)
library(vegan)
library(psych)
library(plyr)
library(dplyr)
library(reshape2)
library(corrplot)
library(stringr)
library(viridis)
library(RColorBrewer) 
library(colourpicker) 
library(ggsci)
library(scales)
library(pheatmap) 
library(patchwork) 

setwd()



# >> envfit analysis with 10,000 permutations =======
## An envfit analysis with 10,000 permutations was performed to fit phenotypes and 
## α-diversity index onto PCoA ordination
## envfit function from the R package vegan

load("~/phe_data.Rdata")
load("~/MetaPhlAn.Species.Rdata") # from file <Microbiome.analysis.R>
load("~/Alpha.result.Rdata") # from file <Microbiome.analysis.R>


phe <- cbind(Alpha_index,phe_data)
phe <- data.frame(apply(phe, 2,function(x) as.numeric(as.character(x))))
rownames(phe) <- rownames(phe_data)

envfitVarsTouse <- colnames(phe)
phe$ID <- rownames(phe)
MetaPhlAn.S$ID <- rownames(MetaPhlAn.S)

#Enfitresult <- NULL
for (i in envfitVarsTouse) {
   print (paste('  >> collecting complete cases for',i))
   inPhe <- phe[, c(i,"ID")] 
   allDF <- merge(inPhe,MetaPhlAn.S,by="ID") 
   rownames(allDF) <- allDF$ID
   allDF$ID <- NULL 
   allDF <- allDF[complete.cases(allDF),] 
   av <- allDF[[i]] 
   allDF[[i]] <- NULL 
   nrRows <- length(av) 

   print ('  >> calculating B/C distance')
   inBC <- vegdist(allDF,method = "bray")
   pCoa <- cmdscale(inBC, eig = T,k = 2 )
   
   #envfit
   print(paste0('>> preparing envfit for ',colnames(phe[i])))
   env <- envfit(pCoa,av,permutations=10000)
   
   #summary
   Phenotype <- colnames(phe[i])
   NR_nonNA = nrRows 
   NAs = 921-NR_nonNA
   Dim <- env$vectors$arrows
   r2 <- env$vectors$r
   pval <- env$vectors$pvals
   df <- cbind(Phenotype,NR_nonNA,NAs,Dim,r2,pval)
   Enfitresult <- data.frame(rbind(Enfitresult,df))
}

Enfitresult$NR_nonNA <- as.numeric(Enfitresult$NR_nonNA)
Enfitresult$NAs <- as.numeric(Enfitresult$NAs)
Enfitresult$Dim1 <- as.numeric(Enfitresult$Dim1)
Enfitresult$Dim2 <- as.numeric(Enfitresult$Dim2)
Enfitresult$r2 <- as.numeric(Enfitresult$r2) #r2 
Enfitresult$pval <- as.numeric(Enfitresult$pval)

# P-value correction FDR.BH
Enfitresult$FDR <- p.adjust(Enfitresult$pval,method = "BH")
Enfitresult$Significant <- ifelse(Enfitresult$FDR<0.05,"Yes","No")


enfit.taxa.result <- Enfitresult[order(Enfitresult$r2,decreasing = T),]
rownames(enfit.taxa.result) <- NULL

#write
write.csv(enfit.taxa.result,row.names = F,
          file = "~\\Phenotype_enfit_taxa.csv")




# >> permutational multivariate analysis of variance (PERMANOVA) =======
## The proportion of variance of Bray-Curtis distance that can be explained by each phenotype
## adonis function in R package vegan

load("~/phe_data.Rdata") 
load("~/MetaPhlAn.Species.Rdata") 


adonisVarsTouse <- colnames(phe_data)
phe_data$ID <- rownames(phe_data)
MetaPhlAn.S$ID <- rownames(MetaPhlAn.S)

#adonisResults <- NULL
for (i in adonisVarsTouse) {
   print (paste('  >> collecting complete cases for',i))
   inPhe <- phe_data[, c(i,"ID")]
   allDF <- merge(inPhe,MetaPhlAn.S,by="ID") 
   rownames(allDF) <- allDF$ID
   allDF$ID <- NULL 
   allDF <- allDF[complete.cases(allDF),] 
   av <- allDF[[i]] 
   allDF[[i]] <- NULL 
   nrRows <- length(av) 
   
   print ('  >> calculating B/C distance')
   inBC <- vegdist(allDF,method = "bray")
   
   #univariate adonis
   print ('  >> doing adonis')
   ad <- adonis2(inBC ~ av,permutations=10000,parallel=4) 
   
   #summary
   Phenotype=i
   NR_nonNA = nrRows 
   NAs = 921-NR_nonNA
   DF = ad[1,1]
   SumsOfSqs = ad[1,2]
   R2 = ad[1,3] 
   FModel = ad[1,4]
   pval = ad[1,5] 
   oneRow <- cbind.data.frame(Phenotype,NR_nonNA,NAs,DF,SumsOfSqs,FModel,R2,pval)
   adonisResults <- rbind.data.frame(adonisResults,oneRow)
}


adonisResults$NR_nonNA <- as.numeric(adonisResults$NR_nonNA)
adonisResults$NAs <- as.numeric(adonisResults$NAs)
adonisResults$SumsOfSqs <- as.numeric(adonisResults$SumsOfSqs)
adonisResults$FModel <- as.numeric(adonisResults$FModel)
adonisResults$R2 <- as.numeric(adonisResults$R2)
adonisResults$pval <- as.numeric(adonisResults$pval)

## P-value correction FDR.BH
adonisResults$FDR <- p.adjust(adonisResults$pval, method = "BH")
adonisResults$Significant <- ifelse(adonisResults$FDR<0.05,"Yes","No")

adonis.taxa.result <- adonisResults[order(adonisResults$R2,decreasing = T),]
rownames(adonis.taxa.result) <- NULL


## write
write.csv(adonis.taxa.result,row.names = F,
          file = "~\\Phenotype_adonis_taxa.csv")
save(adonis.taxa.result,file = "~\\adonis.taxa.result.Rdata")



# >> multivariate PERMANOVA ======
## total proportion of variance in microbiome composition explained by each phenotype group 
## function adonis2 in R package vegan


load("~/MetaPhlAn.Species.Rdata")
load("~/phe_data.Rdata")
load("~/adonis.taxa.result.Rdata")


### Collinearity phenotypes(Spearman |r| > 0.8) that showed significant association (FDR<0.05) 
### with microbiome composition in the univariate adonis analyses
adonis.FDR.05 <- adonis.taxa.result[adonis.taxa.result$FDR<0.05,]
corr_data <- phe_data[,colnames(phe_data) %in% adonis.FDR.05$Phenotype]
corr_data <- apply(corr_data, 2,function(x) as.numeric(as.character(x)))

# spearman
corr <- corr.test(corr_data,method = "spearman")
spearman <- corr$ci
spearman$p.adj <- corr$ci2$p.adj
spearman <- cbind(data.frame(t(combn(rownames(corr$r), 2))),spearman)
spearman <- plyr::rename(spearman,c("X1"="Phenotype1","X2"="Phenotype2"))

# spearman|r|>0.8
collinear <- dplyr::filter(spearman,abs(r)>0.8)
df <- collinear[,c("Phenotype1","Phenotype2","r")]
df1 <- adonis.FDR.05[,c("Phenotype","Group","R2")]
coll <- merge(df,df1,by.x = "Phenotype1",by.y = "Phenotype")
colnames(coll)[4:5] <- c("Phenotype1.Group","Phenotype1.R2")

coll <- merge(coll,df1,by.x = "Phenotype2",by.y = "Phenotype")  
colnames(coll)[6:7] <- c("Phenotype2.Group","Phenotype2.R2")
coll$Variable_excluded_for_subsequent_analyses <- ifelse(coll$Phenotype1.R2>coll$Phenotype2.R2,
                                                         coll$Phenotype2,coll$Phenotype1)

## write
write.csv(coll,file ="~\\Collinearity_taxa.csv",row.names = F)



### Multivariate ADONIS 
adonis.FDR.sig <- adonis.taxa.result[adonis.taxa.result$FDR<0.05,]
adonis.FDR.sig <- adonis.FDR.sig[!adonis.FDR.sig$Phenotype %in% coll$Variable_excluded_for_subsequent_analyses,]

sigPhenos <- adonis.FDR.sig$Phenotype
phe_data$ID <- rownames(phe_data)
MetaPhlAn.S$ID <- rownames(MetaPhlAn.S)
group <- unique(adonis.FDR.sig$Group)

#resAdonisMVT <- NULL
for (gg in group) {
   df <- dplyr::filter(adonis.FDR.sig,Category==gg)
   inPhenos <- phe_data[,c(df$Phenotype,"ID")]
   allDF <- merge(inPhenos,MetaPhlAn.S,by="ID") 
   rownames(allDF) <- allDF$ID 
   allDF$ID <- NULL

   print (paste('  >> collecting complete cases for group',gg))
   allDF <- allDF[complete.cases(allDF),]
   av <- as.data.frame(allDF[,colnames(allDF) %in% sigPhenos])
   if (ncol(av) == 1) {colnames(av) <- colnames(inPhenos[1])}
   allDF <- allDF[,!(colnames(allDF) %in% sigPhenos)]

   print ('  >> calculating B/C distance')
   inBC <- vegdist(allDF,method = "bray")
   
   #multivariate adonis
   print ('  >> doing multivariate adonis')
   print(paste0('length',length(av[[1]]),'unique',length(unique(av[[1]]))))
   frm <- reformulate(termlabels=colnames(av),response='inBC')
   print ('  >> running summary adonis')
   ad <- adonis2(formula=frm,data=av,permutations=10000,by=NULL)
   print(ad)
   inR <- ad
   R2 <- ad$R2[1]
   pV <- ad$`Pr(>F)`[1]
   resAdonisMVT <- rbind.data.frame(resAdonisMVT,data.frame(Group=gg,R2=R2,pV=pV))
   print (paste0('--- DONE with group ',gg))
}

sum(resAdonisMVT$R2)



### make plot
resAdonisMVT$Group <- factor(as.character(resAdonisMVT$Group),
                             levels = resAdonisMVT$Group[order(resAdonisMVT$R2,decreasing = F)])
resAdonisMVT <- resAdonisMVT[order(resAdonisMVT$R2, decreasing = T),]
resAdonisMVT$csumR2 <- resAdonisMVT$R2[1]
for (rn in c(2:nrow(resAdonisMVT))) {resAdonisMVT$csumR2[rn] <- resAdonisMVT$csumR2[rn-1] + resAdonisMVT$R2[rn]}

resAdonisMVT$ypos <- resAdonisMVT$csumR2
for (rn in c(1:nrow(resAdonisMVT))) {resAdonisMVT$ypos[rn] <- resAdonisMVT$ypos[rn] - resAdonisMVT$R2[rn]/2}
resAdonisMVT$R2lbl <- paste0(format(resAdonisMVT$R2*100,digits = 1),'%')
resAdonisMVT$Data <- "Taxa"
rownames(resAdonisMVT) <- resAdonisMVT$Group

#ggplot2
gg1 <- ggplot(resAdonisMVT,aes(x = Data, y=R2,col=Group,fill=Group))+ 
   geom_col() + 
   theme_classic()+
   geom_text(aes(y = ypos, label = R2lbl), color = "black") + 
   labs(x="",y="Explained variation in Bray-Curtis distance (R2)")+
   ylim(0,0.55)+
   theme(legend.title = element_text(size=12,face ="bold",family="serif"),
         legend.text =element_text(size=10,face ="bold",family="serif"),
         axis.text.y=element_text(size=12,face ="bold",family="serif"),
         axis.text.x=element_text(size=12,face ="bold",family="serif"),
         axis.title.x=element_text(size = 14,face ="bold",family="serif"),
         axis.title.y=element_text(size = 14,face ="bold",family="serif"),
   )

gg1 #400*600/600*800

