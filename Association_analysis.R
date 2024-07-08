# =============================================================
# =============================================================
## The source coda for CEMC Microbiome Project in Southwest China ####
## Calculation of microbial features associated with phenotypes ####
# =============================================================
# =============================================================


# >> load libraries ==================
install.packages("")

library(ggplot2)
library(reshape2) 
library(dplyr)
library(plyr)
library(lme4) 
library(lmerTest)
library(partR2) 
library(compositions)
library(pheatmap)
library(data.table)
library(pheatmap)
library(foreach) 
library(doParallel)
library(parallel)

setwd()


# >> Microbiome-phenotype associations (Multivariable linear regression) =====


load("~/MetaPhlAn.Species.Rdata")
load("~/phe_data.Rdata")


### CLR transformation
do_MetaPhlAn_CLR <- function(microbedata){
   #'@microbedata Relative abundance data
   
   if(any(microbedata==0)) microbedata = microbedata + min(microbedata[microbedata>0])/2 
   gm_mean = function(x, na.rm=TRUE){exp(sum(log(x), na.rm=na.rm)/length(x))}
   Gmean = apply(microbedata, 1, gm_mean)
   data_prepared = cbind(Gmean,microbedata)
   data_transformed = t(apply(data_prepared,1,function(x){log(x/x[1])[-1]}))
   colnames(data_transformed) = colnames(data_transformed)
   rownames(data_transformed) = rownames(data_transformed)
   
   return(data_transformed)
}

MetaPhlAn_trans <- do_clr_externalWeighting(MetaPhlAn.S)


### phenotype data
## adjust for city, defecation frequency, batch, and sample month
covar <- c("city","Defecation_frequency","batch","sample_month")

## adjust for age, gender, BMI, city, defecation frequency, batch, and sample month
#covar <- c("Age","gender","BMI","city","Defecation_frequency","batch","sample_month")

## adjust for annual household income, city, defecation frequency, batch, and sample month
#covar <- c("annual_household_income","city","Defecation_frequency","batch","sample_month")

covardata <- phe_data[,covar]
pheno <- phe_data[,!colnames(phe_data)%in%covar]


### association analysis
detectCores()
registerDoSEQ() 
registerDoParallel(makeCluster(8))

ftrs_transformed <- MetaPhlAn_trans
result_ftrs <- NULL
result_ftrs = foreach(i = 1:ncol(pheno),.combine = rbind) %:% 
   foreach(j = 1:ncol(ftrs_transformed),.combine = rbind) %dopar% { 
      predictors = data.frame(covardata[!is.na(pheno[,i]),],
                              model.matrix(as.formula(paste0("~ ",colnames(pheno)[i])),data = pheno)[,-1,drop = F]) 

      cleaned_data = predictors[complete.cases(predictors),]
      rn <- rownames(cleaned_data)
      rn <- rn[rn %in% rownames(ftrs_transformed)]
      ftrs.cleaned = ftrs_transformed[rn,]
      cleaned_data = cleaned_data[rn,]
      if (nrow(cleaned_data) > 3) {
         
         # make model
         s1 = lm(
            as.formula(paste("ftrs.cleaned[,j] ~ ",paste(collapse=" + ",colnames(cleaned_data)))),
            data = cleaned_data)
         
         # make model with extra covariates 
         s0 = lm(
            as.formula(paste("ftrs.cleaned[,j] ~ ",paste(collapse=" + ",colnames(cleaned_data)[1:ncol(covardata)]))),
            data = cleaned_data)
         
         # compare models 
         an1 = anova(s1,s0)
         output = data.frame(
            Phenotype = colnames(pheno)[i],
            taxon = colnames(ftrs.cleaned)[j],
            Nsamples = nrow(cleaned_data),
            levels = if(class(pheno[,i]) == "factor") paste(collapse=":",levels(pheno[,i])) else "Not Applicable",
            levels_SampleSize = 
               if(class(pheno[,i]) == "factor" | length(table(pheno[,i]))==2) paste(collapse= ":",table(pheno[,i])) else "Not Applicable",
            
            effect.size = if(class(pheno[,i]) == "factor") {paste(collapse = ":",c(0,round(digits = 5, s1$coef[grep(colnames(pheno)[i],names(s1$coef))])))
            } else round(digits = 5, s1$coef[colnames(pheno)[i]]),
            
            SEs = if(class(pheno[,i]) == "factor") {paste(collapse = ":",c(0,round(digits = 5, summary(s1)$coef[grep(colnames(pheno)[i],names(s1$coef)),"Std. Error"])))
            } else round(digits = 5, summary(s1)$coef[colnames(pheno)[i],"Std. Error"]),
            
            R2 = summary(s1)$r.squared- summary(s0)$r.squared,
            F.stat = an1[2,5],
            Pvalue = an1[2,6]
         )
      }
   }

#close parallelization
registerDoSEQ()
stopCluster(makeCluster(8))

MetaPhlAn_result <- result_ftrs
rownames(MetaPhlAn_result) <- NULL
MetaPhlAn_result$FDR <- p.adjust(MetaPhlAn_result$Pvalue,method = "BH")
MetaPhlAn_result_out <- MetaPhlAn_result[MetaPhlAn_result$Pvalue<0.05,]
   
# write
write.table(MetaPhlAn_result_out,row.names = F,sep = "\t",
            file ="~/MetaPhlAn_association_overview.txt")




# >> Microbiome-phenotype associations (Linear mixed-effects models, LMMs) =====
 
load("~/MetaPhlAn.Species.Rdata")
load("~/phe_data.Rdata")

### CLR transformation
do_MetaPhlAn_CLR <- function(microbedata){
   #'@microbedata Relative abundance data
   
   if(any(microbedata==0)) microbedata = microbedata + min(microbedata[microbedata>0])/2 
   gm_mean = function(x, na.rm=TRUE){exp(sum(log(x), na.rm=na.rm)/length(x))}
   Gmean = apply(microbedata, 1, gm_mean)
   data_prepared = cbind(Gmean,microbedata)
   data_transformed = t(apply(data_prepared,1,function(x){log(x/x[1])[-1]}))
   colnames(data_transformed) = colnames(data_transformed)
   rownames(data_transformed) = rownames(data_transformed)
   
   return(data_transformed)
}

MetaPhlAn_trans <- do_clr_externalWeighting(MetaPhlAn.S)


### phenotype data
## adjust for age, gender, BMI, defecation frequency, batch, and sample month with city as a random effect
covar <- c("Age","gender","BMI","city","Defecation_frequency","batch","sample_month")
covardata <- phe_data[,covar]
pheno <- phe_data[,!colnames(phe_data)%in%covar]



### association analysis
registerDoSEQ() 
registerDoParallel(makeCluster(8))


ftrs_transformed <- MetaPhlAn_trans 
result_ftrs <- NULL

for (i in 1:ncol(pheno)) {
   foreach(j = 1:ncol(ftrs_transformed),.combine = rbind) %do% { 
      
      predictors = data.frame(covardata[!is.na(pheno[,i]),],
                              model.matrix(as.formula(paste0("~ ",colnames(pheno)[i])),data = pheno)[,-1,drop = F]) 
     
      cleaned_data = predictors[complete.cases(predictors),]
      rn <- rownames(cleaned_data)
      rn <- rn[rn %in% rownames(ftrs_transformed)]
      ftrs.cleaned = ftrs_transformed[rn,]
      cleaned_data = cleaned_data[rn,]
      
      # make model
      s1 = lmer(
         as.formula(paste("ftrs.cleaned[,j] ~ ",paste(collapse=" + ",colnames(cleaned_data)[2:ncol(cleaned_data)]),"+ (1|city)")),
         data = cleaned_data,control=lmerControl(check.conv.singular="ignore"))
      
      # make model with extra covariates
      s0 = lmer(
         as.formula(paste("ftrs.cleaned[,j] ~ ",paste(collapse=" + ",colnames(cleaned_data)[2:ncol(covardata)]),"+ (1|city)")),
         data = cleaned_data,control=lmerControl(check.conv.singular="ignore"))
      
      # partR2()
      partR2.R2 = data.frame(partR2(s1,partvars=colnames(cleaned_data)[grep(colnames(pheno)[i],colnames(cleaned_data))],
                                    R2_type = "marginal")$R2)
      
      # compare models
      an1 = anova(s1,s0,refit=FALSE)
      
      output = data.frame(
         Phenotype = colnames(pheno)[i],
         taxon = colnames(ftrs.cleaned)[j],
         Nsamples = nrow(cleaned_data),
         levels = if(class(pheno[,i]) == "factor") paste(collapse=":",levels(pheno[,i])) else "Not Applicable",
         levels_SampleSize = 
            if(class(pheno[,i]) == "factor" | length(table(pheno[,i]))==2) paste(collapse= ":",table(pheno[,i])) else "Not Applicable",
         
         effect.size = if(class(pheno[,i]) == "factor") {
            paste(collapse = ":",c(0,round(digits = 5, fixef(s1)[grep(colnames(pheno)[i],names(fixef(s1)))])))
         }else round(digits = 5, fixef(s1)[colnames(pheno)[i]]),
         
         SEs = if(class(pheno[,i]) == "factor") {paste(collapse = ":",c(0,round(digits = 5, summary(s1)$coef[grep(colnames(pheno)[i],rownames(summary(s1)$coef)),"Std. Error"])))
         } else round(digits = 5, summary(s1)$coef[colnames(pheno)[i],"Std. Error"]),
         
         Marginal.R2 = partR2.R2[nrow(partR2.R2),2],
         Statistic = an1[2,6],
         Pvalue = an1[2,8]
      )  
      
      result_ftrs <- rbind(result_ftrs,output)
   }
}

#close parallelization
registerDoSEQ()
stopCluster(makeCluster(8))

MetaPhlAn_result <- result_ftrs
rownames(MetaPhlAn_result) <- NULL
MetaPhlAn_result$FDR <- p.adjust(MetaPhlAn_result$Pvalue,method = "BH")
MetaPhlAn_result_out <- MetaPhlAn_result[MetaPhlAn_result$Pvalue<0.05,]

#write 
write.table(MetaPhlAn_result_out,row.names = F,sep="\t",
            file = "~/MetaPhlAn_Mixture_association_disease.txt")



