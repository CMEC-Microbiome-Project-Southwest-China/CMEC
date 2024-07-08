# =============================================================
# =============================================================
## The source coda for CEMC Microbiome Project in Southwest China ####
# =============================================================
# =============================================================

# >> load libraries ==================
install.packages("")

library(tidyverse)
library(ggplot2)
library(ggrepel)
library(microbiome)
library(stringr)
library(dplyr)
library(vegan)
library(iNEXT) 
library(GMCM) 
library(matrixStats) 
library(psych) 
library(clusterSim)
library(ade4)
library(oddsratio) 
library(igraph)
library(foreach)
library(ggridges) 
library(RColorBrewer) 
library(colourpicker) 
library(viridis) 
library(gghalves)
library(pheatmap) 
library(patchwork)
library(colourpicker)
library(igraph)
library(scales)
library(reshape2)
library(gridExtra)
source('~/R_Microbiome_scripts.R')



# >> 100x sampling based core microbiome selection =====================
source('~/R_Microbiome_scripts.R')
load("~/MetaPhlAn.result.Rdata")

# Species
MetaPhlAn.S <- filterMetaGenomeDF(inDF=t(MetaPhlAn.result),presPerc=-1,minMRelAb=-1,minMedRelAb=-1,
                                  keepLevels=c("S"),rescaleTaxa = F,keepDomains = "All")
MetaPhlAn.S <- MetaPhlAn.S/100

# Genus
MetaPhlAn.G <- filterMetaGenomeDF(inDF=t(MetaPhlAn.result),presPerc=-1,minMRelAb=-1,minMedRelAb=-1,
                                  keepLevels=c("G"),rescaleTaxa = F,keepDomains = "All")
MetaPhlAn.G <- MetaPhlAn.G/100

# filter
MetaPhlAn.S.core <- MetaPhlAn.S[rowSums(MetaPhlAn.S)>0,colSums(MetaPhlAn.S>0)>=nrow(MetaPhlAn.S)*0.05]
#core_members(t(MetaPhlAn.S),detection = 0,prevalence = 5/100,include.lowest=FALSE)


## 100x sampling
sp <- MetaPhlAn.S.core
result=NULL
sd=NULL
percent=seq(0.01,1,by=0.01)

for (i in percent) {
   print(i)
   tmp_percent=NULL
   for(j in 1:100){
      permutation=sample(row.names(sp),round(nrow(sp)*i,0))
      per_sp=sp[permutation,]
      tmp_percent=rbind(tmp_percent,colSums(per_sp>0)/round(nrow(sp)*i,0))
   }
   result=cbind(result,colMeans(tmp_percent))
   sd=cbind(sd,colSds(tmp_percent))
}
colnames(result)=1:100
colnames(sd)=1:100

res2 <- result
row.names(res2) <- str_split_fixed(row.names(result),".s__",2)[,2]
row.names(res2)[row.names(res2)==""] <- row.names(result)[row.names(res2)==""]
res2 <- cbind(res2,mean_rate=rowMeans(res2))
colnames(res2)[101]="mean_rate"


## plot prevalence VS resampling for each microbiome feature
pdf("~/core_microbiome_sampling.pdf",width = 10,height = 10,useDingbats = F)
par(mfrow=c(5,2),mgp=c(1,2,0.2))
for(i in 1:nrow(res2)){
   plot(as.numeric(colnames(res2)[1:100]),res2[i,1:100],type = "b",cex=0.4,frame = FALSE, pch = 19, col = "lightpink", xlab = "Sampling percentage (%)", ylab = "Presence rate",main = row.names(res2)[i],cex.main=0.8,xlim = c(0,100),ylim = c(min(res2[i,1:100])-max(sd[i,1:100]),max(res2[i,1:100])+max(sd[i,1:100])))
   segments(as.numeric(colnames(res2)[1:100])-0.02,res2[i,1:100]-sd[i,1:100],as.numeric(colnames(res2)[1:100])+0.02,res2[i,1:100]+sd[i,1:100],col = "gray60",lty = 1,cex=0.02)
}
dev.off()



data=data.frame(t(res2[,1:100]))
data$id=row.names(data)
data=melt(data,id="id")
data$value=data$value*100
data$id=as.numeric(as.character(data$id))
p1=ggplot(data, aes(id, value, group=variable)) +
   geom_line(color="black", linewidth=0.05,alpha=0.8)+
   xlab(label = "Sampling percentage (%)") + 
   ylab(label = "Presence rate (%)")+
   guides(linetype="none",shape="none")+
   theme_bw()+
   theme(panel.grid=element_blank())+ 
   theme(legend.position = "none")+
   scale_x_continuous(limits = c(0,100),breaks = c(0,20,40,60,80,100)
                      ) 


data=data.frame(res2)
data=data[order(data$mean_rate,decreasing = T),]
data$order=1:nrow(data)
data$mean_rate=data$mean_rate*100
data_tmp=data
data_tmp$name=row.names(data_tmp)
data_tmp$name[which(data_tmp$mean_rate<90)]=NA
p2=ggplot(data_tmp,aes(order,mean_rate))+
   geom_point(shape=20,alpha=1,cex=0.2)+
   geom_hline(yintercept = 80,linetype=3,colour="mediumpurple1",alpha=0.5)+
   labs(x="Rank of core microbiome",y="Presence rate (%)")+
   theme_bw()+
   theme(panel.grid=element_blank(),axis.title=element_text(size=8))+
   theme(legend.position="none")+
   geom_text_repel(aes(order,mean_rate,label=factor(data_tmp$name)),
             size=2.5,alpha=1,colour="black",direction = "both",max.overlaps=50,
             segment.size=0.1,segment.alpha = 0.7,segment.color="red")


res2=data.frame(res2)
res2$mean_rate=res2$mean_rate*100
data=data.frame(percent=1:100,number=NA)
for(i in 1:100){
   data$number[i]=length(which(res2$mean_rate>=i))
}


p3=ggplot(data, aes(percent, number)) +
   geom_line(color="black", linewidth=0.5,alpha=0.8)+
   xlab(label = "Presence rate (%)") + 
   ylab(label = "Number of core microbiome")+
   guides(linetype=FALSE,shape=FALSE)+
   theme_bw()+
   theme(panel.grid=element_blank())+ 
   theme(legend.position = "none")+
   scale_x_continuous(limits = c(0,100),breaks = c(0,20,40,60,80,100)) 

svg(file="~/core_microbiome_path.svg")
grid.arrange(p3,p2,ncol=1,nrow = 2)
dev.off()



# >> taxa abundance description =====================
source('~/R_Microbiome_scripts.R')
load("~/MetaPhlAn.result.Rdata")

## function obtain descriptive features
do_MetaPhlAn_statistical <- function(microdata,Taxonomic_level="Species"){
   #'@microdata
   
   Taxon = rownames(microdata)
   Taxonomic_level = rep(Taxonomic_level,length(Taxon))
   microdata = microdata/100
   Summary = t(apply(microdata, 1, FUN=summary))
   colnames(Summary)=c("Min.Ab","Q1.Ab","Median.Ab","Mean.Ab","Q3.Ab","Max.Ab")
   SD.Ab = apply(microdata, 1, FUN=sd)
   Prevalence = apply(microdata, 1, function(x) length(which(x>0)))/ncol(microdata)
   Summary = cbind(Taxon,Taxonomic_level,Summary,SD.Ab,Prevalence)
   Summary = data.frame(Summary[,c("Taxon","Taxonomic_level","Mean.Ab","SD.Ab","Min.Ab","Q1.Ab","Median.Ab",
                                   "Q3.Ab","Max.Ab","Prevalence")])
   Summary[,3:10] = apply(Summary[,3:10],2,function(x) as.numeric(x))
   Summary = Summary[!Summary$Prevalence==0,]
   Summary
}

## Kingdom
MetaPhlAn.K <- filterMetaGenomeDF(inDF=t(MetaPhlAn.result),presPerc=-1,minMRelAb=-1,minMedRelAb=-1,
                                  keepLevels=c("K"),rescaleTaxa = F,keepDomains = "All")
#colnames(MetaPhlAn.K) <- str_split_fixed(colnames(MetaPhlAn.K),"k__",2)[,2]

Statistical.K <- do_MetaPhlAn_statistical(microdata = t(MetaPhlAn.K),Taxonomic_level="Kingdom")
Statistical.K <- Statistical.K[order(Statistical.K$Prevalence,decreasing = T),]


## Phylum
MetaPhlAn.P <- filterMetaGenomeDF(inDF=t(MetaPhlAn.result),presPerc=-1,minMRelAb=-1,minMedRelAb=-1,
                                  keepLevels=c("P"),rescaleTaxa = F,keepDomains = "All")
#colnames(MetaPhlAn.P) <- str_split_fixed(colnames(MetaPhlAn.P),".p__",2)[,2]

Statistical.P <- do_MetaPhlAn_statistical(microdata = t(MetaPhlAn.P),Taxonomic_level="Phylum")
Statistical.P <- Statistical.P[order(Statistical.P$Prevalence,decreasing = T),]


## Class
MetaPhlAn.C <- filterMetaGenomeDF(inDF=t(MetaPhlAn.result),presPerc=-1,minMRelAb=-1,minMedRelAb=-1,
                                  keepLevels=c("C"),rescaleTaxa = F,keepDomains = "All")
#colnames(MetaPhlAn.C) <- str_split_fixed(colnames(MetaPhlAn.C),".c__",2)[,2]

Statistical.C <- do_MetaPhlAn_statistical(microdata = t(MetaPhlAn.C),Taxonomic_level="Class")
Statistical.C <- Statistical.C[order(Statistical.C$Prevalence,decreasing = T),]


## Order
MetaPhlAn.O <- filterMetaGenomeDF(inDF=t(MetaPhlAn.result),presPerc=-1,minMRelAb=-1,minMedRelAb=-1,
                                  keepLevels=c("O"),rescaleTaxa = F,keepDomains = "All")
#colnames(MetaPhlAn.O) <- str_split_fixed(colnames(MetaPhlAn.O),".o__",2)[,2]
MetaPhlAn.O <- MetaPhlAn.O[,-which(colSums(MetaPhlAn.O)==0)]

Statistical.O <- do_MetaPhlAn_statistical(microdata = t(MetaPhlAn.O),Taxonomic_level="Order")
Statistical.O <- Statistical.O[order(Statistical.O$Prevalence,decreasing = T),]

## Family
MetaPhlAn.F <- filterMetaGenomeDF(inDF=t(MetaPhlAn.result),presPerc=-1,minMRelAb=-1,minMedRelAb=-1,
                                  keepLevels=c("F"),rescaleTaxa = F,keepDomains = "All")
#colnames(MetaPhlAn.F) <- str_split_fixed(colnames(MetaPhlAn.F),".f__",2)[,2]
MetaPhlAn.F <- MetaPhlAn.F[,-which(colSums(MetaPhlAn.F)==0)]

Statistical.F <- do_MetaPhlAn_statistical(microdata = t(MetaPhlAn.F),Taxonomic_level="Family")
Statistical.F <- Statistical.F[order(Statistical.F$Prevalence,decreasing = T),]

## Genus
MetaPhlAn.G <- filterMetaGenomeDF(inDF=t(MetaPhlAn.result),presPerc=-1,minMRelAb=-1,minMedRelAb=-1,
                                  keepLevels=c("G"),rescaleTaxa = F,keepDomains = "All")
#colnames(MetaPhlAn.G) <- str_split_fixed(colnames(MetaPhlAn.G),".g__",2)[,2]
MetaPhlAn.G <- MetaPhlAn.G[,-which(colSums(MetaPhlAn.G)==0)]

Statistical.G <- do_MetaPhlAn_statistical(microdata = t(MetaPhlAn.G),Taxonomic_level="Genus")
Statistical.G <- Statistical.G[order(Statistical.G$Prevalence,decreasing = T),]


## Speices
MetaPhlAn.S <- filterMetaGenomeDF(inDF=t(MetaPhlAn.result),presPerc=-1,minMRelAb=-1,minMedRelAb=-1,
                                  keepLevels=c("S"),rescaleTaxa = F,keepDomains = "All")
#colnames(MetaPhlAn.S) <- str_split_fixed(colnames(MetaPhlAn.S),".s__",2)[,2]
MetaPhlAn.S <- MetaPhlAn.S[,!colSums(MetaPhlAn.S)==0]

Statistical.S <- do_MetaPhlAn_statistical(microdata = t(MetaPhlAn.S),Taxonomic_level="Species")
Statistical.S <- Statistical.S[order(Statistical.S$Prevalence,decreasing = T),]


## write
Statistical.results <- rbind(Statistical.K,Statistical.P,Statistical.C,Statistical.O,Statistical.F,
                             Statistical.G,Statistical.S)
write.csv(Statistical.results,row.names = F,
          file = "~/Statistical.results.csv")
save(MetaPhlAn.S,file = "~/MetaPhlAn.Species.Rdata")



# >> Percentage Stacking Chart =====================
load("~/MetaPhlAn.result.Rdata")
source('~/R_Microbiome_scripts.R')

MetaPhlAn.P <- filterMetaGenomeDF(inDF=t(MetaPhlAn.result),presPerc=-1,minMRelAb=-1,minMedRelAb=-1,
                                  keepLevels=c("P"),rescaleTaxa = F,keepDomains = "All")
colnames(MetaPhlAn.P) <- str_split_fixed(colnames(MetaPhlAn.P),".p__",2)[,2]

Phylum.ave <- apply(MetaPhlAn.P, 2, FUN=mean) 
Phyla_data <- rbind(MetaPhlAn.P,Phylum.ave)[,order(Phylum.ave, decreasing=TRUE)]
P_Topten <- cbind(Phyla_data[-nrow(Phyla_data),1:10],Other=(100-apply(Phyla_data[-nrow(Phyla_data),1:10], 1, sum)))
P_Topten[P_Topten<0] <- 0
P_Topten <- P_Topten[order(P_Topten[,1],decreasing = T),]

P_Topten <- data.frame(t(P_Topten))
P_Topten$Taxonomy <- rownames(P_Topten)
P_Topten$Taxonomy <- factor(P_Topten$Taxonomy,levels=P_Topten$Taxonomy)
P_Topten_melt <- melt(P_Topten,value.name="Abundance",id.vars="Taxonomy", variable.name="ID")


# ggplot2
#my_col <- brewer.pal(11,"Set3") 
my_col <- c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5",
            "#D9D9D9","#BC80BD","#CCEBC5")

gg <- ggplot(P_Topten_melt,aes(x=ID, y=Abundance, fill = Taxonomy)) +
   geom_col(position = 'stack', width = 0.6) + 
   scale_y_continuous(expand=c(0, 0))+
   geom_bar(stat="identity",width=1)+
   scale_fill_manual(values=my_col,name="Phylum" ) + 
   labs(y="Relative abundance (%)",family="serif") + 
   theme_classic()+
   theme(legend.title = element_text(size=14,family ="serif",face="bold"), 
         legend.text = element_text(size = 12,family ="serif",face="bold"),
         axis.text.y = element_text(size = 12,family ="serif",face="bold"), 
         axis.title.y  = element_text(size = 14,family ="serif",face="bold"), 
         axis.title.x=element_blank(),
         axis.text.x=element_blank(),
         axis.ticks.x=element_blank(),
   )

gg #1200*400



# >> Alpha diversity==================
load("~/MetaPhlAn.Species.Rdata")

Shannon <- vegan::diversity(MetaPhlAn.S,index = "shannon")
Observed_species <- specnumber(MetaPhlAn.S,MARGIN = 1)
Alpha_index <- cbind.data.frame(Shannon,Observed_species)

# save
save(Alpha_index,file ="~/Alpha.result.Rdata")

# ggplot2
gg1 <- ggplot(Alpha_index,aes(x=Shannon)) + 
   geom_histogram(binwidth=0.10,fill="#69b3a2",color="#e9ecef",alpha=0.9)+
   scale_x_continuous(expand=c(0,0))+
   scale_y_continuous(expand=c(0,0))+
   theme_bw()+
   labs(x="Shannon diversity index",y="Frequency")+
   theme(panel.grid = element_blank(),
         axis.title.x = element_text(size=14,family ="serif",face="bold"),
         axis.title.y = element_text(size=14,family ="serif",face="bold"),
         axis.text = element_text(size=12,family ="serif",face="bold"),
         )

gg1 #400*350



# >> PCoA analysis=====================
MetaPhlAn.S <- filterMetaGenomeDF(inDF=t(MetaPhlAn.result),presPerc=-1,minMRelAb=-1,minMedRelAb=-1,
                                  keepLevels=c("S"),rescaleTaxa = F,keepDomains = "All")
colnames(MetaPhlAn.S) <- str_split_fixed(colnames(MetaPhlAn.S),".s__",2)[,2]
MetaPhlAn.S <- MetaPhlAn.S[,!colSums(MetaPhlAn.S)==0]


Bray_curtis_Species <- vegdist(MetaPhlAn.S, method = 'bray')
pCoa <- cmdscale(Bray_curtis_Species, eig = T,k = 2 )
varExp <- (eigenvals(pCoa)/sum(eigenvals(pCoa)))[1:2]
xVar <- as.numeric(varExp[1]*100) 
yVar <- as.numeric(varExp[2]*100) 

pcoa_Axis <- as.data.frame(pCoa$points)
colnames(pcoa_Axis) <- paste0("PCo",c(1:2))
pcoa_Axis$ID <- row.names(pcoa_Axis)


## Correlation analysis between PCoA1, PCoA2 and microbiome features
pco_species <- NULL
for (i in c("PCo1","PCo2")) {
   for (j in colnames(MetaPhlAn.S)) {
      corr <- corr.test(pcoa_Axis[,i],MetaPhlAn.S[,j],method = "spearman",adjust = "fdr")
      Axis = i
      Species = j
      Type = "Taxon"
      result  = cbind(Axis,Species,Type,corr$ci2[,c(2,4,5)])
      pco_species = rbind(pco_species,result)
   }
}

rownames(pco_species) <- NULL
colnames(pco_species) <- c("Principal_Coordinate","Microbiome_feature","Type",
                           "Spearman_Correlation","Pvalue","FDR")
# write
write.csv(pco_species,row.names = F,
          file = "~/MetaPhlAn/PCo_Species_Correlation.csv")


## PCoA plot 
pCoa_Taxa <- MetaPhlAn.S[,c("Bacteroides_vulgatus","Prevotella_copri",
                            "Ruminococcus_gnavus","Faecalibacterium_prausnitzii")]
pCoa_Taxa <- cbind(pcoa_Axis,pCoa_Taxa)

g1 <- ggplot(pCoa_Taxa,aes(x=PCo1,y=PCo2,color=Bacteroides_vulgatus)) + 
   geom_point(size=1.25) + 
   theme_classic() + 
   scale_colour_viridis(direction = -1, name="Relative \nabundance %")+ 
   labs(x=paste0("PCo1 ",xVar),
        y=paste0("PCo2 ",yVar),
        title="Bacteroides (Phocaeicola) vulgatus") +  
   theme(plot.title = element_text(face = "bold.italic",family="serif"),
         text = element_text(size = 12),
         legend.text = element_text(size=12,family="serif"),
         axis.title = element_text(size=14,family="serif"),
         legend.title = element_text(size=12,face = "bold",family="serif"),
          )


g2 <- ggplot(pCoa_Taxa,aes(x=PCo1,y=PCo2,color=Prevotella_copri)) + 
   geom_point(size=1.25) + 
   theme_classic() + 
   scale_colour_viridis(direction = -1, name="Relative \nabundance %")+ 
   labs(x=paste0("PCo1 ",xVar),
        y=paste0("PCo2 ",yVar),
        title="Prevotella copri (Segatella copri)") +  
   theme(plot.title = element_text(face = "bold.italic",family="serif"),
         text = element_text(size = 12),
         legend.text = element_text(size=12,family="serif"),
         axis.title = element_text(size=14,family="serif"),
         legend.title = element_text(size=12,face = "bold",family="serif"),
         )


g3 <- ggplot(pCoa_Taxa,aes(x=PCo1,y=PCo2,color=Ruminococcus_gnavus)) + 
   geom_point(size=1.25) + 
   theme_classic() + 
   scale_colour_viridis(direction = -1, name="Relative \nabundance %")+ 
   labs(x=paste0("PCo1 ",xVar),
        y=paste0("PCo2 ",yVar),
        title="Ruminococcus gnavus") +  
   theme(plot.title = element_text(face = "bold.italic",family="serif"),
         text = element_text(size = 12),
         legend.text = element_text(size=12,family="serif"),
         axis.title = element_text(size=14,family="serif"),
         legend.title = element_text(size=12,face = "bold",family="serif"),
   )


g4 <- ggplot(pCoa_Taxa,aes(x=PCo1,y=PCo2,color=Faecalibacterium_prausnitzii)) + 
   geom_point(size=1.25) + 
   theme_classic() + 
   scale_colour_viridis(direction = -1, name="Relative \nabundance %")+ 
   labs(x=paste0("PCo1 ",xVar),
        y=paste0("PCo2 ",yVar),
        title="Faecalibacterium prausnitzii") +  
   theme(plot.title = element_text(face = "bold.italic",family="serif"),
         text = element_text(size = 12),
         legend.text = element_text(size=12,family="serif"),
         axis.title = element_text(size=14,family="serif"),
         legend.title = element_text(size=12,face = "bold",family="serif"),
   )


## Summary
g1+g2+g3+g4



## enfit analysis
load("~/Alpha.result.Rdata")
load("~/phe_data.Rdata")

vars <- c("City","Age","Gender","BMI","Defecation_frequency")
phedata <- phe_data[,vars]
phedata$Shannon <- Alpha_index$Shannon

phedata <- as.data.frame(apply(phedata,2,function(x){as.numeric(as.character(x))}))
rownames(phedata) <- rownames(phe_data)


env <- envfit(pCoa~.,phedata,na.rm=T,permutations=10000)
en_coord_cont <- as.data.frame(scores(env, "vectors")) * ordiArrowMul(env)
coord2 <- en_coord_cont/2
pCoa_Taxa$Shannon <- Alpha_index$Shannon

# ggplot2
g2 <- ggplot(pCoa_Taxa,aes(x=PCo1,y=PCo2,color=Shannon)) + 
   geom_point(size=1.25) + 
   theme_classic() + 
   scale_color_viridis(direction = -1,name="Shannon \ndiversity") +
   labs(x=paste0("PCo1 ",xVar),y=paste0("PCo2 ",yVar)) +
   geom_segment(data = coord2,aes(x = 0, y = 0, xend = Dim1, yend = Dim2),linewidth=1, colour = "black",
                arrow = arrow(length = unit(0.02, "npc")) ) +
   geom_label_repel(data = coord2, aes(x = Dim1, y = Dim2), 
                    colour = "black", family="serif",
                    label = row.names(coord2),size=4,segment.colour = NA) + 
   theme(legend.title = element_text(size=12,face ="bold",family="serif"),
         legend.text =element_text(size=12,family="serif"),
         axis.text.y=element_text(size=14,family="serif"),
         axis.text.x=element_text(size=14,family="serif"),
         axis.title.x=element_text(size = 16,family="serif"),
         axis.title.y=element_text(size = 16,family="serif"),
   )

print(g2) #450*350



# >> Top species =====================
load("~/MetaPhlAn.Species.Rdata")

Speices.ave <- apply(MetaPhlAn.S, 2, FUN=mean) 
Speices_data <- rbind(MetaPhlAn.S,Speices.ave)[,order(Speices.ave, decreasing=TRUE)]
S_Topten <- Speices_data[-nrow(Speices_data),1:10]


S_Topten <- data.frame(t(S_Topten))
S_Topten$Taxonomy <- rownames(S_Topten)
S_Topten$Taxonomy <- factor(S_Topten$Taxonomy,levels=S_Topten$Taxonomy)
S_Topten_melt <- melt(S_Topten,value.name="Abundance",id.vars="Taxonomy", variable.name="ID")

S_Topten_melt <- S_Topten_melt[!S_Topten_melt$Abundance==0,]
S_Topten_melt$Relative <- -log2(S_Topten_melt$Abundance)

## Ridge map
palette <- c("#D5DCB9","#80CED5","#CFBAD1","#DC7B6C","#DB73C8","#D7C56A","#A744DE","#A3E85B",
             "#878EDA","#78DFA4")

gg <- ggplot(S_Topten_melt,aes(x=Relative,y=Taxonomy,fill=Taxonomy,color=Taxonomy))+ 
   theme_classic() +
   geom_density_ridges(alpha=0.8,scale = 2)+ 
   scale_fill_manual(values = palette) +
   scale_color_manual(values = palette)+
   labs(x="Relative abundance (-log2)",y="")+
   theme(axis.text.y=element_text(size=12,face ="bold.italic",family="serif"),
         axis.text.x=element_text(size=12,face ="bold",family="serif"),
         axis.title.x=element_text(size =14,face ="bold",family="serif"),
         legend.position  ="None",)

gg #600*300



# >> Enterotype analysis =====================
source('~/R_Microbiome_scripts.R')
load("~/MetaPhlAn.result.Rdata")
load("~/phe_data.Rdata")

## function calculator for JSD distance
dist.JSD <- function(inMatrix, pseudocount=0.000001, ...) {
   KLD <- function(x,y) sum(x *log(x/y))
   JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
   matrixColSize <- length(colnames(inMatrix))
   matrixRowSize <- length(rownames(inMatrix))
   colnames <- colnames(inMatrix)
   resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
   
   inMatrix = apply(inMatrix,1:2,function(x) ifelse (x==0,pseudocount,x))
   
   for(i in 1:matrixColSize) {
      for(j in 1:matrixColSize) { 
         resultsMatrix[i,j]=JSD(as.vector(inMatrix[,i]),
                                as.vector(inMatrix[,j]))
      }
   }
   colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
   as.dist(resultsMatrix)->resultsMatrix
   attr(resultsMatrix, "method") <- "dist"
   return(resultsMatrix) 
}


## function partition around medoid (PAM) clustering
pam.clustering=function(x,k) { # x is a distance matrix and k the number of clusters
   cluster = as.vector(pam(as.dist(x), k, diss=TRUE)$clustering)
   return(cluster)
}


## function CLR-transformed
do_MetaPhlAn_CLR <- function(microbedata){
   #'@microbedata
   if(any(microbedata==0)) microbedata = microbedata + min(microbedata[microbedata>0])/2 
   gm_mean = function(x, na.rm=TRUE){exp(sum(log(x), na.rm=na.rm)/length(x))}
   Gmean = apply(microbedata, 1, gm_mean)
   data_prepared = cbind(Gmean,microbedata)
   data_transformed = t(apply(data_prepared,1,function(x){log(x/x[1])[-1]}))
   colnames(data_transformed) = colnames(data_transformed)
   rownames(data_transformed) = rownames(data_transformed)
   
   return(data_transformed)
}


### genus clustering using PAM
MetaPhlAn.G <- filterMetaGenomeDF(inDF=t(MetaPhlAn.result),presPerc=-1,minMRelAb=-1,minMedRelAb=-1,
                                  keepLevels=c("G"),rescaleTaxa = F,keepDomains = "All")
colnames(MetaPhlAn.G) <- str_split_fixed(colnames(MetaPhlAn.G),".g__",2)[,2]

MetaPhlAn.data <- MetaPhlAn.G
MetaPhlAn.data <- MetaPhlAn.data[,!colSums(MetaPhlAn.data)==0]
MetaPhlAn.data <- MetaPhlAn.data/100


PAM.dist <- dist.JSD(t(MetaPhlAn.data))
nclusters=NULL
for (k in 1:20) {
   if (k==1) {
      nclusters[k]=NA
   } else {
      cluster_temp=pam.clustering(PAM.dist, k)
      nclusters[k]=index.G1(MetaPhlAn.data,cluster_temp,d=PAM.dist,
                            centrotypes = "medoids")
   }
}
plot(nclusters, type="h", xlab="k clusters", ylab="CH index")

cluster <- pam.clustering(PAM.dist, k=2)
table(cluster)


# save
cluster_data <- data.frame(row.names = rownames(MetaPhlAn.data),Cluster=cluster)
write.csv(cluster_data,"~/Cluster.result.csv")
save(MetaPhlAn.data,PAM.dist,cluster_data,file = "~/Enterotype.Rdata")


###  BCA between-class analysis
load("~/Enterotype.Rdata")

pca <- dudi.pca(MetaPhlAn.data, scannf=F, nf=10)
bca <- bca(pca,fac=as.factor(cluster_data$Cluster),scannf=F,nf=2) 
pcoa <- dudi.pco(PAM.dist, scannf=F, nf=2)
driven <- t(bca$tab)


## 600*600 "#66CD00"
s.class(pcoa$li,fac=as.factor(cluster_data$Cluster),grid=F,col=c("#87CEFA", "#FFC125"))


## Violin figure
plot.data <- data.frame(do_MetaPhlAn_CLR(MetaPhlAn.data[,c("Bacteroides","Prevotella")]))
plot.data$ID <- rownames(plot.data)
plot.data <- melt(plot.data,value.name = "Abundance",id.vars = "ID",variable.name="Taxon")
plot.data <- merge(plot.data,cluster_data,by.x = "ID",by.y = "row.names",all.x = T)
plot.data$Cluster <- as.factor(plot.data$Cluster)

gg <- ggplot(data=plot.data, aes(x=Taxon,y=Abundance,color=Cluster)) + 
   geom_half_violin(aes(x=Taxon,y=Abundance,split=Cluster,fill=Cluster),
                    position = "identity",fill="white",linewidth=1)+
   geom_jitter(alpha=0.25,position=position_jitterdodge(jitter.width=0.15,
                                                        jitter.height=0,dodge.width=0.15))+ 
   geom_half_boxplot(data = dplyr::filter(plot.data,Cluster=="1"),aes(Taxon,Abundance),
                     width = 0.1,side = "l",outlier.colour = NA) +
   geom_half_boxplot(data = dplyr::filter(plot.data,Cluster=="2"),aes(Taxon,Abundance),
                     width = 0.1,side = "r",outlier.colour = NA) +
   scale_color_manual(values = c("1"="#87CEFA", 
                                 "2"="#FFC125"))+ 
   labs(y="Relative abundance (CLR transformation)",x="")+
   #stat_compare_means(method="wilcox.test",label="p.format",vjust=0.01,family="serif",size=5)+ 
   theme_classic()+
   theme(axis.title.x=element_text(size=14,face ="bold",family="serif"),
         axis.title.y=element_text(size =14,face ="bold",family="serif"),
         axis.text.x=element_text(size =11,face ="bold.italic",family="serif"),
         plot.subtitle = element_text(size =14,face ="bold",family="serif"),
         legend.position  ="right",)
gg # 500*400




# >> Species composition distribution of different enterotypes ====
load("~/Enterotype.Rdata")
load("~/MetaPhlAn.Species.Rdata")
load("~/MetaPhlAn.tree.Rdata")

cluster_data$rownames <- rownames(cluster_data)
MetaPhlAn.S <- MetaPhlAn.S/100
total.tree$Genus <- gsub("g__","",total.tree$Genus)
total.tree$Species <- gsub("s__","",total.tree$Species)


# >> genus

###  cluster 1 
cluster1 <- rownames(cluster_data[cluster_data$Cluster=="1",])

Genus_cluster1 <- t(MetaPhlAn.data[cluster1,])
Genus_ave <- apply(Genus_cluster1, 1, FUN=mean)
Genus_cluster1 <- cbind(Genus_cluster1, Genus_ave)[order(Genus_ave, decreasing=T),]
Genus_cluster1 <- data.frame(Genus_cluster1[1:10,-ncol(Genus_cluster1)])

Genus_cluster1$Taxonomy <- rownames(Genus_cluster1)
Genus_cluster1 <- melt(Genus_cluster1, id.vars="Taxonomy", variable.name="ID",
                       value.name="Abundance")
Genus_cluster1 <- merge(Genus_cluster1,unique(total.tree[,c(4,8)]),
                        by.x="Taxonomy",by.y = "Genus",all.x = T)
table(Genus_cluster1$Phylum,exclude = NULL)


# ggplot2
g1 <- ggplot(Genus_cluster1, aes(x=reorder(Taxonomy,-Abundance),
                                 y=Abundance,color=Phylum)) + 
   geom_boxplot()+
   scale_fill_manual(values ="white") +
   scale_color_manual(values = c("p__Actinobacteria"="#1B9E77",
                                 "p__Bacteroidetes"="#D95F02",
                                 "p__Firmicutes"="#7570B3",
                                 "p__Proteobacteria"="#E7298A"),
                      name="Phylum")+
   theme_classic()+
   labs(y="Relative abundance",x="",subtitle = "Bacteroides")+
   theme(axis.text.x=element_text(size=12,face ="bold.italic",family="serif",
                                  vjust=1,hjust =1.0,angle = 90),
         axis.text.y=element_text(size=14,face ="bold",family="serif"),
         axis.title.y=element_text(size =16,face ="bold",family="serif"),
         legend.title = element_text(size=14,face ="bold",family="serif"),
         legend.text = element_text(size=12,face ="bold.italic",family="serif"),
         legend.position  ="right",
         )

print(g1) # 600*600


### cluster 2 genus
cluster2 <- rownames(cluster_data[cluster_data$Cluster=="2",])

Genus_cluster2 <- t(MetaPhlAn.data[cluster2,])
Genus_ave <- apply(Genus_cluster2, 1, FUN=mean)
Genus_cluster2 <- cbind(Genus_cluster2, Genus_ave)[order(Genus_ave, decreasing=T),]
Genus_cluster2 <- data.frame(Genus_cluster2[1:10,-ncol(Genus_cluster2)])

Genus_cluster2$Taxonomy <- rownames(Genus_cluster2)
Genus_cluster2 <- melt(Genus_cluster2, id.vars="Taxonomy", variable.name="ID",
                       value.name="Abundance")
Genus_cluster2 <- merge(Genus_cluster2,unique(total.tree[,c(4,8)]),
                        by.x="Taxonomy",by.y = "Genus",all.x = T)
table(Genus_cluster2$Phylum,exclude = NULL)


# ggplot2
g2 <- ggplot(Genus_cluster2, aes(x=reorder(Taxonomy,-Abundance),
                                 y=Abundance,color=Phylum)) + 
   geom_boxplot()+
   scale_fill_manual(values ="white") +
   scale_color_manual(values = c("p__Actinobacteria"="#1B9E77",
                                 "p__Bacteroidetes"="#D95F02",
                                 "p__Firmicutes"="#7570B3",
                                 "p__Proteobacteria"="#E7298A"),
                      name="Phylum")+
   theme_classic()+
   labs(y="Relative abundance",x="",subtitle = "Prevotella")+
   theme(axis.text.x=element_text(size=12,face ="bold.italic",family="serif",vjust=1,
                                  hjust =1.0,angle = 90),
         axis.text.y=element_text(size=14,face ="bold",family="serif"),
         axis.title.y=element_text(size =16,face ="bold",family="serif"),
         legend.title = element_text(size=14,face ="bold",family="serif"),
         legend.text = element_text(size=12,face ="bold.italic",family="serif"),
         legend.position  ="right",
         )

print(g2) # 600*600



# >> speices
### cluster 1
cluster1 <- rownames(cluster_data[cluster_data$Cluster=="1",])

Species_cluster1 <- t(MetaPhlAn.S[cluster1,])
Species_ave <- apply(Species_cluster1, 1, FUN=mean)
Species_cluster1 <- cbind(Species_cluster1, Species_ave)[order(Species_ave, decreasing=T),]
Species_cluster1 <- data.frame(Species_cluster1[1:10,-ncol(Species_cluster1)])


Species_cluster1$Taxonomy <- rownames(Species_cluster1)
Species_cluster1 <- melt(Species_cluster1, id.vars="Taxonomy", variable.name="ID",
                         value.name="Abundance")
Species_cluster1 <- merge(Species_cluster1,unique(total.tree[,c(8,9)]),
                          by.x="Taxonomy",by.y = "Species",all.x = T)
table(Species_cluster1$Genus,exclude = NULL)


palette <- c(brewer.pal(7,"Dark2"),brewer.pal(6,"Set3"))
g3 <- ggplot(Species_cluster1, aes(x=reorder(Taxonomy,-Abundance),
                                   y=Abundance,color=Genus)) + 
   geom_boxplot()+
   scale_fill_manual(values ="white") +
   scale_color_manual(values = palette,name="Genus")+
   theme_classic()+
   labs(y="Relative abundance",x="",subtitle = "Bacteroides")+
   theme(axis.text.x=element_text(size=12,face ="bold.italic",family="serif",vjust=1,
                                  hjust =1.0,angle = 90),
         axis.title.y=element_text(size =16,face ="bold",family="serif"),
         legend.title = element_text(size=14,face ="bold",family="serif"),
         legend.text = element_text(size=12,face ="bold.italic",family="serif"),
         legend.position  ="right",
         )

print(g3) # 600*600


### cluster 2
cluster2 <- rownames(cluster_data[cluster_data$Cluster=="2",])

Species_cluster2 <- t(MetaPhlAn.S[cluster2,])
Species_ave <- apply(Species_cluster2, 1, FUN=mean)
Species_cluster2 <- cbind(Species_cluster2, Species_ave)[order(Species_ave, decreasing=T),]
Species_cluster2 <- data.frame(Species_cluster2[1:10,-ncol(Species_cluster2)])


Species_cluster2$Taxonomy <- rownames(Species_cluster2)
Species_cluster2 <- melt(Species_cluster2,id.vars="Taxonomy",variable.name="ID",
                         value.name="Abundance")
Species_cluster2 <- merge(Species_cluster2,unique(total.tree[,c(8,9)]),
                          by.x="Taxonomy",by.y = "Species",all.x = T)
table(Species_cluster2$Genus,exclude = NULL)


palette <- c(brewer.pal(7,"Dark2"),brewer.pal(6,"Set3"),brewer.pal(4,"Set1"))
g4 <- ggplot(Species_cluster2, aes(x=reorder(Taxonomy,-Abundance),
                                   y=Abundance,color=Genus)) + 
   geom_boxplot()+
   scale_fill_manual(values ="white") +
   scale_color_manual(values = palette,name="Genus")+
   theme_classic()+
   labs(y="Relative abundance",x="",subtitle = "Prevotella")+
   theme(axis.text.x=element_text(size=12,face ="bold.italic",family="serif",vjust=1,
                                  hjust =1.0,angle = 90),
         axis.title.y=element_text(size =16,face ="bold",family="serif"),
         legend.title = element_text(size=14,face ="bold",family="serif"),
         legend.text = element_text(size=12,face ="bold.italic",family="serif"),
         legend.position  ="right",
         )

print(g4) # 1200*600


# >> Summary
(g1+g3)/(g2+g4) #1400*1200




# >>  The association between enterotypes and diseases ========
load("~/phe_data.Rdata")
load("~/Variable.description.Rdata")
load("~/Enterotype.Rdata")

## Spearman correlation between enterotypes and phenotype
phe <- phe_data
phe$batch <- factor(phe$batch,labels=c(1:14))
phe <- data.frame(apply(phe,2,function(x) as.character(x)))
phe[phe=="Q1"] <- "1"
phe[phe=="Q2"] <- "2"
phe[phe=="Q3"] <- "3"
phe[phe=="Q4"] <- "4"
phe <- data.frame(apply(phe, 2,function(x) as.numeric(x)))
rownames(phe) <- rownames(phe_data)

enterotype_covar <- NULL
for (j in colnames(phe)) {
   corr <- corr.test(phe[,j],cluster_data[,1],method = "spearman")
   Phenotype = j
   result  = cbind(Phenotype,corr$ci2[,c(2,4)])
   enterotype_covar = rbind(enterotype_covar,result)
}

enterotype_covar$FDR <- p.adjust(enterotype_covar$p,method = "fdr")
enterotype_covar <- merge(enterotype_covar,Variable.description[,c(1,2,5,7)],
                          by = "Phenotype",all.x = T)
enterotype_covar <- enterotype_covar[,c(5:7,2:4,1)]
enterotype_covar <- enterotype_covar[order(enterotype_covar$p,decreasing = F),]

# save
write.csv(enterotype_covar[,!colnames(enterotype_covar)%in%c("Phenotype")],row.names = F,
          file = "~/Enterotype.phenotype.csv")


## enterotypes and diseases
vars0 <- c("Hypertension","Prehypertension","Diabetes","Prediabetes","Obesity","Overweight",
           "MetS","NAFLD","Dyslipidemia","Dyslipidemia","Health")
Disease.phe <- phe_data[,vars0]

vars1 <- c("age","A1","BMI","city","E4","batch","sample_month")
covar.phe <- phe_data[,vars1]

Disease <- data.frame(Phenotype=colnames(Disease.phe))

# No covariates
Disease_cluster <- NULL
Disease_cluster <-
   foreach(i=1:nrow(Disease),.combine = rbind) %do%  {
      Disease.sub=Disease.phe[,colnames(Disease.phe)==Disease$Phenotype[i],drop=F]
      Disease.sub=merge(Disease.sub,cluster_data,
                        by="row.names",all.x = T)
      if(nrow(Disease.sub)>100){
         fit_glm=glm(Disease.sub[,2]~factor(Disease.sub[,3]), family=binomial)
         pvalue=summary(fit_glm)$coefficients[2,4] 
         or=or_glm(data = Disease.sub, model = fit_glm, incr = 0.95)
         return.string = data.frame(Disease=Disease$Phenotype[i],Pvalue=pvalue,OR=or[1,2],CI_up=or[1,4],CI_low=or[1,3])
      }
   }

Disease_cluster$FDR <- p.adjust(Disease_cluster$Pvalue,method = "BH")
Disease_cluster <- Disease_cluster %>%
   mutate(Sig=ifelse(FDR<0.05,"Significant","Non-Significant"))
Disease_cluster <- Disease_cluster[order(Disease_cluster$OR),]

# save
write.csv(Disease_cluster,row.names = F,
          file = "~/Enterotype.N.disease.csv")


# with covariates
Disease_cluster <- NULL
for (i in 1:ncol(Disease.phe)) {
   Disease.sub=Disease.phe[,colnames(Disease.phe)==Disease$Phenotype[i],drop=F]
   Disease.sub=cbind(Disease.sub,covar.phe)
   Disease.sub=merge(Disease.sub,cluster_data,
                     by="row.names",all.x = T)
   colnames(Disease.sub)=c("ID","Disease",vars1,"Cluster")
   Disease.sub = Disease.sub[complete.cases(Disease.sub),]
   if(nrow(Disease.sub)>100){
      fit_glm=glm(Disease~factor(Cluster)+age+A1+city,family=binomial,data=Disease.sub)
      pvalue=summary(fit_glm)$coefficients[2,4]
      OR=exp(coef(fit_glm)[2])
      CI = exp(confint(fit_glm)[2, ])
      CI.LOW = CI[1]
      CI.HIGH = CI[2]
      return.string = data.frame(Disease=Disease$Phenotype[i],Pvalue=pvalue,OR=OR,
                                 CI_low= CI.LOW,CI_up=CI.HIGH)
      Disease_cluster <- rbind(Disease_cluster,return.string)
   }
}

Disease_cluster$FDR <- p.adjust(Disease_cluster$Pvalue,method = "BH")
Disease_cluster <- Disease_cluster %>%
   mutate(Sig=ifelse(FDR<0.05,"Significant","Non-Significant"))
Disease_cluster <- Disease_cluster[order(Disease_cluster$OR),]

# save
write.csv(Disease_cluster,row.names = F,
          file = "~/Enterotype.disease.csv")



## ggplot2
library(ggsci)
Disease_cluster$Disease <- factor(Disease_cluster$Disease,
                                  levels = Disease_cluster$Disease)
g4 <- ggplot(Disease_cluster, aes(x=OR, y=Disease, color=Sig)) + 
   geom_vline(aes(xintercept = 1), linewidth = .25, linetype = "dashed") + 
   geom_errorbarh(aes(xmax = `CI_up`, xmin = `CI_low`),linewidth = .3, height = 0.2) +
   geom_point(size = 3.5) +
   labs(x="Odds ratio",y="") + 
   scale_color_lancet(name="Significance")+
   theme_bw()+ 
   theme(panel.grid.minor = element_blank(),
         axis.title.x=element_text(size=14,face ="bold",family="serif"),
         axis.text.y=element_text(size=12,face ="bold",family="serif"),
         axis.text.x=element_text(size=12,face ="bold",family="serif"),
         legend.title = element_text(size=14,family="serif"),
         legend.text = element_text(size=12,family="serif"),
   )



# >> Species accumulation curve =====================
source('~/R_Microbiome_scripts.R')
load("~/MetaPhlAn.result.Rdata")

## function 
do_accumulation = function(MetaPhlAn,type="O"){
   #'@microdata 
   #'@type
   
   # debug
   if(!type%in%c("Order","Family","Genus","Species")){stop("type only in Order/Family/Genus/Species")}
   
   MetaPhlAn.data <- filterMetaGenomeDF(inDF=t(MetaPhlAn),presPerc=-1,minMRelAb=-1,minMedRelAb=-1,
                                     keepLevels=type,rescaleTaxa = F,keepDomains = "All")
   MetaPhlAn.data <- MetaPhlAn.data[,!colSums(MetaPhlAn.data)==0]
   MetaPhlAn.data[MetaPhlAn.data>0] <- 1
   # collector", "random", "exact", "rarefaction", "coleman"
   Accumulation <- specaccum(MetaPhlAn.data,permutations = 100,method = 'rarefaction')
   Curves <- data.frame(richness=Accumulation$richness,sites=round(Accumulation$sites),
                        sd=Accumulation$sd)
   Curves$type <- type
   
   return(Curves)
}

# Species
Curves_S <- do_accumulation(MetaPhlAn.result,type ="Species")
# Genus
Curves_G <- do_accumulation(MetaPhlAn.result,type ="Genus")
# Family
Curves_F <- do_accumulation(MetaPhlAn.result,type ="Family")
# Combined
MetaPhlAn_Curves_data <- rbind(Curves_S,Curves_G,Curves_F)
MetaPhlAn_Curves_filter <- MetaPhlAn_Curves_data[MetaPhlAn_Curves_data$sites%in%c(1,2,4,8,10,20,30,40,50,
                                                                                  75,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,
                                                                                  921),]
MetaPhlAn_Curves_filter$type[MetaPhlAn_Curves_filter$type=="S"] <- "Species"
MetaPhlAn_Curves_filter$type[MetaPhlAn_Curves_filter$type=="G"] <- "Genera"
MetaPhlAn_Curves_filter$type[MetaPhlAn_Curves_filter$type=="F"] <- "Families"
MetaPhlAn_Curves_filter$linetype <- "solid"


plot.data <- rbind(MetaPhlAn_Curves_filter,
                    Curves_data_filter[Curves_data_filter$type%in%c("Pathways"),])
plot.data2 <- Curves_data_filter[Curves_data_filter$type%in%c("KOs"),]

cbPalette <- c("KOs"="#CC79A7", "Genera"="#009E73","Pathways"="#D55E00",
               "Species"="#56B4E9","Families"="#E69F00")
gg1 <- ggplot(plot.data[!plot.data$type%in%c("KOs"),], aes(x=sites, y=richness,color=type))+
   geom_point(shape=21,size=2.0,fill="white")+
   theme_classic()+
   scale_x_continuous(breaks = c(0,250,500,750,921))+
   scale_color_manual(values = cbPalette,name="Taxon")+
   geom_line(linewidth=1.15,linetype=plot.data$linetype) +
   geom_errorbar(mapping=aes(ymin=richness-sd, ymax=richness+sd,colour=type),
                 linewidth=0.75,stat="identity")+
   labs(x="Number of samples",y="Number of features")+
   theme(legend.title = element_text(size=14,face ="bold",family="serif"),
         legend.text =element_text(size=12,family="serif"),
         axis.text.y=element_text(size=14,face ="bold",family="serif"),
         axis.text.x=element_text(size=14,face ="bold",family="serif"),
         axis.title.y=element_text(size = 16,face ="bold",family="serif"),
         axis.title.x=element_text(size = 16,face ="bold",family="serif"),
   ) 


gg2 <- ggplot(plot.data2, aes(x=sites, y=richness,color=type))+
   geom_point(shape=21,size=2.0,fill="white")+
   theme_classic()+
   scale_x_continuous(breaks = c(0,250,500,750,921))+
   scale_y_continuous(limits=c(3500,10100))+
   scale_color_manual(values = cbPalette,name="Taxon")+
   geom_line(linewidth=1.15,linetype=plot.data2$linetype) +
   geom_errorbar(mapping=aes(ymin=richness-sd, ymax=richness+sd,colour=type),
                 linewidth=0.75,stat="identity")+
   labs(x="Number of samples",y="Number of features")+
   theme(legend.title = element_text(size=14,face ="bold",family="serif"),
         legend.text =element_text(size=12,family="serif"),
         axis.text.y=element_text(size=14,face ="bold",family="serif"),
         axis.text.x=element_text(size=14,face ="bold",family="serif"),
         axis.title.y=element_text(size = 16,face ="bold",family="serif"),
         axis.title.x=element_text(size = 16,face ="bold",family="serif"),
   ) 

gg2/gg1 #800*600



# >> Rarefaction and extrapolation (R/E) sampling curves  ===============
source('~/R_Microbiome_scripts.R')
load("~/MetaPhlAn.result.Rdata")
load("~/MetaPhlAn.Species.Rdata")

## species
MetaPhlAn.S.total <- data.frame(t(MetaPhlAn.S/100))
MetaPhlAn.S.total[MetaPhlAn.S.total > 0] <- 1
Species.sum <- rowSums(MetaPhlAn.S.total)
Species_T_RE <- iNEXT(Species.sum, q=0,datatype ='abundance',knots = 200,
                      endpoint = sum(Species.sum)*3)


Species_T_RE$iNextEst$size_based$m <- Species_T_RE$iNextEst$size_based$m/sum(Species.sum)*ncol(MetaPhlAn.S.total)
View(Species_T_RE$iNextEst$size_based)

gg1 <- ggiNEXT(Species_T_RE, type=1, se=T,facet.var="None",color.var="Order.q",grey=FALSE)+ 
   theme_classic() + 
   scale_x_continuous(breaks = c(0,921,1842,2763))+
   scale_color_manual(values = "#56B4E9")+
   scale_fill_manual(values = "#56B4E9")+
   labs(y="Number of Species",x="Sample size",
        subtitle = "mean = 1,051\nsd = 21.1\nmean = 1,280\nsd = 57.0",) + 
   theme(legend.position = "None",
         plot.subtitle = element_text(size=14,family = "serif",hjust = 0.5,vjust = -40.5),
         text = element_text(size = 16,family = "serif"),
         axis.text=element_text(size=14,family="serif"),
         axis.title=element_text(size = 16,family="serif"),
   )

print(gg1)


## Genera 
MetaPhlAn.G <- filterMetaGenomeDF(inDF=t(MetaPhlAn.result),presPerc=-1,minMRelAb=-1,minMedRelAb=-1,
                                  keepLevels=c("G"),rescaleTaxa = F,keepDomains = "All")
colnames(MetaPhlAn.G) <- str_split_fixed(colnames(MetaPhlAn.G),".g__",2)[,2]
MetaPhlAn.G <- MetaPhlAn.G[,-which(colSums(MetaPhlAn.G)==0)]

MetaPhlAn.G[MetaPhlAn.G>0] <- 1
MetaPhlAn.G <- data.frame(t(MetaPhlAn.G))

Genus.sum <- rowSums(MetaPhlAn.G)
Genus_T_RE <- iNEXT(Genus.sum, q=0,datatype ='abundance',knots = 200,
                      endpoint = sum(Genus.sum)*3)

# ggplot2
Genus_T_RE$iNextEst$size_based$m <- Genus_T_RE$iNextEst$size_based$m/sum(Genus.sum)*ncol(MetaPhlAn.G)
View(Genus_T_RE$iNextEst$size_based)

gg2 <- ggiNEXT(Genus_T_RE, type=1, se=T,facet.var="None",color.var="Order.q",grey=FALSE)+ 
   theme_classic() + 
   scale_x_continuous(breaks = c(0,921,1842,2763))+
   scale_color_manual(values = "#009E73")+
   scale_fill_manual(values = "#009E73")+
   labs(y="Number of Genera",x="Sample size",
        subtitle = "mean = 305\nsd = 11.4\nmean = 363\nsd = 26.8",) + 
   theme(legend.position = "None",
         plot.subtitle = element_text(size=14,family = "serif",hjust = 0.5,vjust = -40.5),
         text = element_text(size = 16,family = "serif"),
         axis.text=element_text(size=14,family="serif"),
         axis.title=element_text(size = 16,family="serif"),
   )

print(gg2)

## summary
gg2+gg1




# >> SparCC ===============
load("~/MetaPhlAn.Species.Rdata")
load("~/kneaddata.summary.Rdata")

## Based on the processing results of kneaddata, obtain the number of paired reads for each sample, and 
## convert the annotation results of metaphlan into absolute abundance values
MetaPhlAn.S.core <- MetaPhlAn.S[rowSums(MetaPhlAn.S)>0,colSums(MetaPhlAn.S>0)>nrow(MetaPhlAn.S)*0.05]
MetaPhlAn.S.core <- MetaPhlAn.S.core/100

result_table <- result_table[rownames(result_table)%in%rownames(MetaPhlAn.S.core),]
result.reads <- data.frame(reads=result_table[,11])
rownames(result.reads) <- rownames(result_table)

MetaPhlAn.count <- MetaPhlAn.S.core * result.reads$reads
MetaPhlAn.count <- data.frame(t(MetaPhlAn.count))

## python analysis
cor_spracc_species <- read.table("~\\SparCC\\cor_sparcc.txt",
                                 header=TRUE,sep="\t",row.names = 1)
Species_pvalues <- read.table("~\\SparCC\\pvals_two_sided.txt",
                               header=TRUE,sep="\t",row.names = 1)
Species_pvalues_fdr <- data.frame(rowname=rownames(Species_pvalues))

for (i in 1:ncol(Species_pvalues)) {
   d <- p.adjust(Species_pvalues[,i],method = "fdr")
   Species_pvalues_fdr <- cbind(Species_pvalues_fdr,d)
}
rownames(Species_pvalues_fdr) <- Species_pvalues_fdr[,1]
Species_pvalues_fdr <- Species_pvalues_fdr[,-1]
colnames(Species_pvalues_fdr) <- rownames(Species_pvalues_fdr)


## keep result FDR < 0.05
Species_pvalues_fdr[Species_pvalues_fdr>=0.05] <- -1
Species_pvalues_fdr[Species_pvalues_fdr>=0 & Species_pvalues_fdr<0.05 ] <- 1
Species_pvalues_fdr[Species_pvalues_fdr == -1] <- 0


cor_spracc_species <- round(cor_spracc_species,4)
Species.adj <- as.matrix(cor_spracc_species) * as.matrix(Species_pvalues_fdr)
diag(Species.adj) <- 0


Species.sparcc.graph <- graph.adjacency(Species.adj,mode="undirected",
                                        weighted = TRUE, diag = FALSE)
range(E(Species.sparcc.graph)$sparcc)
E(Species.sparcc.graph)$sparcc <- E(Species.sparcc.graph)$weight
E(Species.sparcc.graph)$weight <- abs(E(Species.sparcc.graph)$weight)

# graphml
write.graph(Species.sparcc.graph, 
            file ="~\\SparCC\\MetaPhlAn.network.graphml", 
            format = 'graphml')
