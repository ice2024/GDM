#### MR ####
library(VariantAnnotation)
library(gwasglue)
library(TwoSampleMR)
load("expo_list.rda")
dat_list_res=do.call(rbind,expo_list)
dat_list_res1=dat_list_res[which(dat_list_res$pval.exposure<1e-5),]
dat=dat_list_res1
N=dat[1,"samplesize.exposure"]    
dat=transform(dat,R2=((beta.exposure)^2/((beta.exposure)^2+((se.exposure)^2)*N))) 
dat=transform(dat,F=(N-2)*R2/(1-R2))
dat_list_res1=dat
save(dat_list_res1,file = "expo_list_final.F.rda")
dat_list_res1=dat_list_res1[which(dat_list_res1$F>10),]
table(dat_list_res1$id.exposure)
length(unique(dat_list_res1$id.exposure))
exp_dat=dat_list_res1

exposure_idall1 <- unique(exp_dat$id.exposure)
exposure_idall=exposure_idall1
outcomeID <- c('GEST_DIABETES')
outcomeName <- c('Gestational diabetes')
dat_list <- list()
res_single_list=list()
res_single_OR_list=list()
for (exposureID in exposure_idall) {
  
  exposure_dat=exp_dat[which(exp_dat$id.exposure==exposureID),]
  outcomeTab<-merge(exposure_dat, outcomeData1, by.x="SNP", by.y="SNP")
  if(dim(outcomeTab)[1]>2){
    write.csv(outcomeTab[,-(2:ncol(exposure_dat))], file="outcome.csv")
    outcome_dat<-read_outcome_data(snps=exposure_dat$SNP,
                                   filename="outcome.csv", sep = ",",
                                   snp_col = "SNP",
                                   beta_col = "beta.outcome",
                                   se_col = "se.outcome",
                                   effect_allele_col = "effect_allele.outcome",
                                   other_allele_col = "other_allele.outcome",
                                   pval_col = "pval.outcome",
                                   eaf_col = "eaf.outcome")
    
    outcome_dat$outcome=outcomeName
    outcome_dat$id.outcome=outcomeID
    dat <- TwoSampleMR::harmonise_data(
      exposure_dat = exposure_dat,
      outcome_dat = outcome_dat
    )
    
    res <- mr(dat)
    if(dim(res)[1]>2){
      res_single_list[[exposureID]] <- res
      OR <-generate_odds_ratios(res)
      res_single_OR_list[[exposureID]] <- OR
      dat_list[[exposureID]] <- dat
    }
  }
}
save(dat_list,file = "exposure_outcome_harmonise_all.rda")
save(res_single_list,file = "exposure_outcome_res_list.rda")
save(res_single_OR_list,file = "exposure_outcome_res_OR_list.rda")
res_single_list1=do.call(rbind,res_single_list)
res_single_OR_list1=do.call(rbind,res_single_OR_list)
dat_list1=do.call(rbind,dat_list)
write.csv(dat_list1,"exposure_outcome_harmonise.csv",row.names = FALSE)
write.csv(res_single_list1,"result_MR.csv",row.names = FALSE)
write.csv(res_single_OR_list1,"result_MR_OR.csv",row.names = FALSE)

id=unique(res_single_list1$id.exposure)
exposure_idall=id
library(MRPRESSO)
pleio_list <- list()
hete_list <- list()
presso_list=list()
presso_list_main=list()
for (exposureID in exposure_idall) {
  dat <- dat_list[[exposureID]]
  het <- mr_heterogeneity(dat)
  hete_list[[exposureID]] <- het
  pleio <- mr_pleiotropy_test(dat)
  if(dim(dat)[1]>3){
    presso <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                        OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = dat, NbDistribution = 1000,  
                        SignifThreshold = 0.05)
    presso$`Main MR results`$result=exposureID
    presso_list[[exposureID]]=presso
    presso_list_main[[exposureID]]=presso$`Main MR results`
    
    pleio$PRESSO_P <- presso$`MR-PRESSO results`$`Global Test`$Pvalue
    pleio$RSSobs <- presso$`MR-PRESSO results`$`Global Test`$RSSobs
  }else{
    pleio$PRESSO_P <-""
    pleio$RSSobs <-""
  }
  pleio_list[[exposureID]] <- pleio
}

pleio_merge <- do.call(rbind,pleio_list)
hete_merge <- do.call(rbind,hete_list)
presso_list_main_merge=do.call(rbind,presso_list_main)
write.csv(pleio_merge,"pleiotropy_test.csv",row.names = FALSE)
write.csv(hete_merge,"heterogeneity_test.csv",row.names = FALSE)
write.csv(presso_list_main_merge,"mr_presso_main.csv")
save(presso_list,file = "mr_presso_result.rda")

for (exposureID in exposure_idall){
  library(patchwork)
  library(ggplot2)
  single <- mr_leaveoneout(dat=dat_list[[exposureID]])
  write.csv(single,paste0("leaveone_",exposureID,".csv"),row.names = FALSE)
  pa <- mr_leaveoneout_plot(single)[[1]]
  # pa$theme$line$linewidth <- 3
  pa <- pa+scale_color_brewer(palette = 'Set1')+
    theme_classic(base_size = 18)+
    theme(panel.border = element_rect(size = 1.7,fill = 'transparent'),
          axis.ticks = element_line(size = 1),
          legend.position = 'none')+
    update_geom_defaults("line", list(size = 5))+
    labs(x=paste(dat_list[[exposureID]]$exposure,paste0(' on ',dat_list[[exposureID]]$outcome),sep = "\n"))
  pdf(paste0("leaveone_",exposureID,".pdf"),width = 15,height = 30)
  print(pa)
  dev.off()
  
  pa <- mr_scatter_plot(dat=dat_list[[exposureID]],mr_results = res_single_list[[exposureID]])[[1]]
  pa$labels$x <- unlist(strsplit(pa$labels$x,split = ' \\|\\| '))[1]
  pa$labels$y <- paste0("SNP effect on ",dat_list[[exposureID]]$outcome)
  pa <- pa+scale_color_brewer(palette = 'Set1')+
    theme_classic(base_size = 18)+
    theme(panel.border = element_rect(size = 1.7,fill = 'transparent'),
          axis.ticks = element_line(size = 1),
          legend.position = 'top')+
    update_geom_defaults("line", list(size = 5))+
    labs(x=paste0("SNP effect on ",dat_list[[exposureID]]$exposure))
  pdf(paste0("scatter_",exposureID,".pdf"),height = 8,width = 8)
  print(pa)
  dev.off()
  
  res_single_2 <- mr_singlesnp(dat=dat_list[[exposureID]])
  write.csv(res_single_2,paste0("forest_",exposureID,".csv"),row.names = FALSE)
  pa <- mr_forest_plot(res_single_2)[[1]]
  pa <- pa+scale_color_brewer(palette = 'Set1')+
    theme_classic(base_size = 18)+
    theme(panel.border = element_rect(size = 1.7,fill = 'transparent'),
          axis.ticks = element_line(size = 1),
          legend.position = 'none')+
    update_geom_defaults("line", list(size = 5))+
    labs(x=paste(dat_list[[exposureID]]$exposure,paste0(' on ',dat_list[[exposureID]]$outcome),sep = "\n"))
  pdf(paste0("forest_",exposureID,".pdf"),width = 12,height = 30)
  print(pa)
  dev.off()
  
  ## 漏斗图
  res_single_2 <- mr_singlesnp(dat=dat_list[[exposureID]])
  nrow(res_single_2)
  res_single_2=res_single_2[-nrow(res_single_2),]
  pa <- mr_funnel_plot(res_single_2)[[1]]
  pa <- pa+scale_color_brewer(palette = 'Set1')+
    theme_classic(base_size = 18)+
    theme(panel.border = element_rect(size = 1.7,fill = 'transparent'),
          axis.ticks = element_line(size = 1),
          legend.position = 'top')+
    update_geom_defaults("line", list(size = 5))
  pdf(paste0("funnel_",exposureID,".pdf"),height = 8,width = 8)
  print(pa)
  dev.off()
  
}

library(VariantAnnotation)
library(gwasglue)
library(TwoSampleMR)
exp_dat=read.csv("IEU_Gut.microbiota.abundance_snp_p1e-5_clump.csv",header = T,stringsAsFactors = FALSE,check.names = FALSE,row.names = 1)
colnames(exp_dat)
exp_dat$exposure=NULL
expo=read.csv("Gut microbiota abundance_GWAS_database.csv",header = T,stringsAsFactors = FALSE,check.names = FALSE,row.names = 1)
colnames(expo)[1]="id.exposure"
colnames(expo)[2]="exposure"
expo=expo[,c(1,2)]

exp_dat=merge(expo,exp_dat,by="id.exposure")
dat=exp_dat
dat$samplesize.exposure=as.numeric(dat$samplesize.exposure)
N=dat[1,"samplesize.exposure"]     
dat=transform(dat,R2=((beta.exposure)^2/((beta.exposure)^2+((se.exposure)^2)*N))) #
dat=transform(dat,F=(N-2)*R2/(1-R2))
dat_list_res1=dat
dat_list_res1=dat_list_res1[which(dat_list_res1$F>10),]

length(unique(dat_list_res1$id.exposure))
exp_dat=dat_list_res1
save(exp_dat,file = "IEU_Gut.microbiota.abundance_snp_p1e-5_clump_F.rda")
exposure_idall1 <- unique(exp_dat$id.exposure)
exposure_idall=exposure_idall1

outcomeID <- c('GEST_DIABETES')
outcomeName <- c('Gestational diabetes')
dat_list <- list()
res_single_list=list()
res_single_OR_list=list()
for (exposureID in exposure_idall) {
  
  exposure_dat=exp_dat[which(exp_dat$id.exposure==exposureID),]
  outcomeTab<-merge(exposure_dat, outcomeData1, by.x="SNP", by.y="SNP")
  if(dim(outcomeTab)[1]>2){
    write.csv(outcomeTab[,-(2:ncol(exposure_dat))], file="outcome.csv")
    
    outcome_dat<-read_outcome_data(snps=exposure_dat$SNP,
                                   filename="outcome.csv", sep = ",",
                                   snp_col = "SNP",
                                   beta_col = "beta.outcome",
                                   se_col = "se.outcome",
                                   effect_allele_col = "effect_allele.outcome",
                                   other_allele_col = "other_allele.outcome",
                                   pval_col = "pval.outcome",
                                   eaf_col = "eaf.outcome")
    
    outcome_dat$outcome=outcomeName
    outcome_dat$id.outcome=outcomeID

    dat <- TwoSampleMR::harmonise_data(
      exposure_dat = exposure_dat,
      outcome_dat = outcome_dat
    )
    
    res <- mr(dat)
    if(dim(res)[1]>2){
      res_single_list[[exposureID]] <- res
      OR <-generate_odds_ratios(res)
      res_single_OR_list[[exposureID]] <- OR
      dat_list[[exposureID]] <- dat
    }
  }
}

save(dat_list,file = "exposure_outcome_harmonise_all.rda")
save(res_single_list,file = "exposure_outcome_res_list.rda")
save(res_single_OR_list,file = "exposure_outcome_res_OR_list.rda")
res_single_list1=do.call(rbind,res_single_list)
res_single_OR_list1=do.call(rbind,res_single_OR_list)
dat_list1=do.call(rbind,dat_list)
write.csv(dat_list1,"exposure_outcome_harmonise.csv",row.names = FALSE)
write.csv(res_single_list1,"result_MR.csv",row.names = FALSE)
write.csv(res_single_OR_list1,"result_MR_OR.csv",row.names = FALSE)

id=unique(res_single_list1$id.exposure)
exposure_idall=id
library(MRPRESSO)
pleio_list <- list()
hete_list <- list()
presso_list=list()
presso_list_main=list()
for (exposureID in exposure_idall) {
  dat <- dat_list[[exposureID]]
  het <- mr_heterogeneity(dat)
  hete_list[[exposureID]] <- het
  pleio <- mr_pleiotropy_test(dat)
  if(dim(dat)[1]>3){
    presso <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                        OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = dat, NbDistribution = 1000,  
                        SignifThreshold = 0.05)
    presso$`Main MR results`$result=exposureID
    presso_list[[exposureID]]=presso
    presso_list_main[[exposureID]]=presso$`Main MR results`
    
    pleio$PRESSO_P <- presso$`MR-PRESSO results`$`Global Test`$Pvalue
    pleio$RSSobs <- presso$`MR-PRESSO results`$`Global Test`$RSSobs
  }else{
    pleio$PRESSO_P <-""
    pleio$RSSobs <-""
  }
  pleio_list[[exposureID]] <- pleio
}

pleio_merge <- do.call(rbind,pleio_list)
hete_merge <- do.call(rbind,hete_list)
presso_list_main_merge=do.call(rbind,presso_list_main)
write.csv(pleio_merge,"pleiotropy_test.csv",row.names = FALSE)
write.csv(hete_merge,"heterogeneity_test.csv",row.names = FALSE)
write.csv(presso_list_main_merge,"mr_presso_main.csv")
save(presso_list,file = "mr_presso_result.rda")

for (exposureID in exposure_idall){
 
  library(patchwork)
  library(ggplot2)
  single <- mr_leaveoneout(dat=dat_list[[exposureID]])
  write.csv(single,paste0("leaveone_",exposureID,".csv"),row.names = FALSE)
  pa <- mr_leaveoneout_plot(single)[[1]]
  pa <- pa+scale_color_brewer(palette = 'Set1')+
    theme_classic(base_size = 18)+
    theme(panel.border = element_rect(size = 1.7,fill = 'transparent'),
          axis.ticks = element_line(size = 1),
          legend.position = 'none')+
    update_geom_defaults("line", list(size = 5))+
    labs(x=paste(dat_list[[exposureID]]$exposure,' on ',outcomeName,sep = "\n"))
  pdf(paste0("leaveone_",exposureID,".pdf"))
  print(pa)
  dev.off()
  ## 散点图
 
  pa <- mr_scatter_plot(dat=dat_list[[exposureID]],mr_results = res_single_list[[exposureID]])[[1]]
  pa$labels$x <- unlist(strsplit(pa$labels$x,split = ' \\|\\| '))[1]
  pa$labels$y <- paste('SNP effect on ',outcomeName)##
  pa <- pa+scale_color_brewer(palette = 'Set1')+
    theme_classic(base_size = 18)+
    theme(panel.border = element_rect(size = 1.7,fill = 'transparent'),
          axis.ticks = element_line(size = 1),
          legend.position = 'top')+
    update_geom_defaults("line", list(size = 5))+
    labs(x=paste0("SNP effect on ",dat_list[[exposureID]]$exposure))
  pdf(paste0("scatter_",exposureID,".pdf"),height = 8,width = 8)
  print(pa)
  dev.off()
  ## 森林图
  
  res_single_2 <- mr_singlesnp(dat=dat_list[[exposureID]])
  write.csv(res_single_2,paste0("forest_",exposureID,".csv"),row.names = FALSE)
  pa <- mr_forest_plot(res_single_2)[[1]]
  pa <- pa+scale_color_brewer(palette = 'Set1')+
    theme_classic(base_size = 18)+
    theme(panel.border = element_rect(size = 1.7,fill = 'transparent'),
          axis.ticks = element_line(size = 1),
          legend.position = 'none')+
    update_geom_defaults("line", list(size = 5))+
    labs(x=paste(dat_list[[exposureID]]$exposure,' on ',outcomeName,sep = "\n"))
  pdf(paste0("forest_",exposureID,".pdf"),width = 8)
  print(pa)
  dev.off()
  
  ## 漏斗图
  
  res_single_2 <- mr_singlesnp(dat=dat_list[[exposureID]])
  nrow(res_single_2)
  res_single_2=res_single_2[-nrow(res_single_2),]
  pa <- mr_funnel_plot(res_single_2)[[1]]
  pa <- pa+scale_color_brewer(palette = 'Set1')+
    theme_classic(base_size = 18)+
    theme(panel.border = element_rect(size = 1.7,fill = 'transparent'),
          axis.ticks = element_line(size = 1),
          legend.position = 'top')+
    update_geom_defaults("line", list(size = 5))
  pdf(paste0("funnel_",exposureID,".pdf"),height = 8,width = 8)
  print(pa)
  dev.off()
  
}


library(TwoSampleMR)

jieju=read.csv("result_MR_OR_ivw_p0.05.csv",header = T,stringsAsFactors = FALSE,check.names = FALSE)
jieju=jieju[which(jieju$method=="Inverse variance weighted"),]
jieju=jieju[which(jieju$pval<0.05),]
jieju=jieju[,c(1,4)]
jieju$id.exposure=substr(jieju$id.exposure,7,18)
exposure_idIDall <- c('ebi-a-GCST006098')
outcomeIDall <- jieju$id.exposure

exp_dat=dat_list_res1
dat_list <- list()
res_single_list=list()
res_single_OR_list=list()

for (exposureID in exposure_idIDall) {

  exposure_dat=exp_dat[which(exp_dat$id.exposure==exposureID),]

for (outcomeID in outcomeIDall) {
  
  outcomeFile=paste0(outcomeID,"_buildGRCh37.tsv")
  biofsci.1 <- vroom(outcomeFile, col_names = TRUE)
  colnames(biofsci.1)
  biofsci.1 <- biofsci.1 %>%
    rename(
      SNP = variant_id,
      CHR = chromosome,
      BP = base_pair_location,
      effect_allele = effect_allele,
      other_allele = other_allele,
      P = p_value,
      BETA = beta,
      SE = standard_error,
      N=n
    ) %>%
    select(SNP, CHR, BP, effect_allele, other_allele, P, BETA, SE,N) %>%
    mutate(P = as.numeric(P))
  outcomeData1=biofsci.1
  
  outcomeTab<-merge(exposure_dat, outcomeData1, by.x="SNP", by.y="SNP")
  if(dim(outcomeTab)[1]>2){
    write.csv(outcomeTab[,-(2:ncol(exposure_dat))], file="outcome.csv")
    
    outcome_dat<-read_outcome_data(snps=exposure_dat$SNP,
                                   filename="outcome.csv", sep = ",",
                                   snp_col = "SNP",
                                   beta_col = "BETA",
                                   se_col = "SE",
                                   effect_allele_col = "effect_allele",
                                   other_allele_col = "other_allele",
                                   pval_col = "P",
                                   #eaf_col = "EAF",
                                   samplesize_col = "N",
                                   chr_col = "CHR",
                                   pos_col = "BP")
    outcome_dat$id.outcome=outcomeID
  outcomeName=jieju[which(jieju$id.exposure==outcomeID),2]
  outcome_dat$outcome=outcomeName
  dat <- TwoSampleMR::harmonise_data(
    exposure_dat = exposure_dat,
    outcome_dat = outcome_dat
  )
  res <- mr(dat)
  if(dim(res)[1]>2){
  res_single_list[[paste(exposureID,outcomeID,sep = "_")]] <- res
  OR <-generate_odds_ratios(res)
  res_single_OR_list[[paste(exposureID,outcomeID,sep = "_")]] <- OR
  dat_list[[paste(exposureID,outcomeID,sep = "_")]] <- dat
  }
}
}
}

save(dat_list,file = "exposure_outcome_harmonise_all.rda")
save(res_single_list,file = "exposure_outcome_res_list.rda")
save(res_single_OR_list,file = "exposure_outcome_res_OR_list.rda")
res_single_list1=do.call(rbind,res_single_list)
res_single_OR_list1=do.call(rbind,res_single_OR_list)
dat_list1=do.call(rbind,dat_list)
library(MRPRESSO)
pleio_list <- list()
hete_list <- list()
presso_list=list()
presso_list_main=list()
for (outcomeID in outcomeIDall){
  for (exposureID in exposure_idIDall) {
    dat <- dat_list[[paste(exposureID,outcomeID,sep = "_")]]

    het <- mr_heterogeneity(dat)
    hete_list[[paste(exposureID,outcomeID,sep = "_")]] <- het
    pleio <- mr_pleiotropy_test(dat)
    if(dim(dat)[1]>3){
      presso <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                          OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = dat, NbDistribution = 1000,  
                          SignifThreshold = 0.05)
      presso$`Main MR results`$result=paste(exposureID,outcomeID,sep = "_")
      presso_list[[paste(exposureID,outcomeID,sep = "_")]]=presso
      presso_list_main[[paste(exposureID,outcomeID,sep = "_")]]=presso$`Main MR results`
      
      pleio$PRESSO_P <- presso$`MR-PRESSO results`$`Global Test`$Pvalue
      pleio$RSSobs <- presso$`MR-PRESSO results`$`Global Test`$RSSobs
    }else{
      pleio$PRESSO_P <-""
      pleio$RSSobs <-""
    }
    pleio_list[[paste(exposureID,outcomeID,sep = "_")]] <- pleio
  }
}

pleio_merge <- do.call(rbind,pleio_list)
hete_merge <- do.call(rbind,hete_list)
presso_list_main_merge=do.call(rbind,presso_list_main)
write.csv(pleio_merge,"pleiotropy_test.csv",row.names = FALSE)
write.csv(hete_merge,"heterogeneity_test.csv",row.names = FALSE)
write.csv(presso_list_main_merge,"mr_presso_main.csv")
save(presso_list,file = "mr_presso_result.rda")


for (outcomeID in outcomeIDall) {
  for (exposureID in exposure_idIDall){

    library(patchwork)
    library(ggplot2)
    single <- mr_leaveoneout(dat=dat_list[[paste(exposureID,outcomeID,sep = "_")]])
    single$samplesize="NA"
    write.csv(single,paste0("leaveone_",paste(exposureID,outcomeID,sep = "_"),".csv"),row.names = FALSE)
    pa <- mr_leaveoneout_plot(single)[[1]]
    pa <- pa+scale_color_brewer(palette = 'Set1')+
      theme_classic(base_size = 18)+
      theme(panel.border = element_rect(size = 1.7,fill = 'transparent'),
            axis.ticks = element_line(size = 1),
            legend.position = 'none')+
      update_geom_defaults("line", list(size = 5))+
      labs(x=paste(dat_list[[paste(exposureID,outcomeID,sep = "_")]]$exposure,paste0(' on ',dat_list[[paste(exposureID,outcomeID,sep = "_")]]$outcome),sep = "\n"))
    pdf(paste0("leaveone_",paste(exposureID,outcomeID,sep = "_"),".pdf"),width = 8,height = 15)
    print(pa)
    dev.off()
    
    
    pa <- mr_scatter_plot(dat=dat_list[[paste(exposureID,outcomeID,sep = "_")]],mr_results = res_single_list[[paste(exposureID,outcomeID,sep = "_")]])[[1]]
    pa$labels$x <- unlist(strsplit(pa$labels$x,split = ' \\|\\| '))[1]
    pa$labels$y <- paste0("SNP effect on ",dat_list[[paste(exposureID,outcomeID,sep = "_")]]$outcome)
    pa <- pa+scale_color_brewer(palette = 'Set1')+
      theme_classic(base_size = 18)+
      theme(panel.border = element_rect(size = 1.7,fill = 'transparent'),
            axis.ticks = element_line(size = 1),
            legend.position = 'top')+
      update_geom_defaults("line", list(size = 5))+
      labs(x=paste0("SNP effect on ",dat_list[[paste(exposureID,outcomeID,sep = "_")]]$exposure))
    pdf(paste0("scatter_",paste(exposureID,outcomeID,sep = "_"),".pdf"),height = 8,width = 8)
    print(pa)
    dev.off()
    
    res_single_2 <- mr_singlesnp(dat=dat_list[[paste(exposureID,outcomeID,sep = "_")]])
    res_single_2$samplesize="NA"
    write.csv(res_single_2,paste0("forest_",paste(exposureID,outcomeID,sep = "_"),".csv"),row.names = FALSE)
    pa <- mr_forest_plot(res_single_2)[[1]]
    pa <- pa+scale_color_brewer(palette = 'Set1')+
      theme_classic(base_size = 18)+
      theme(panel.border = element_rect(size = 1.7,fill = 'transparent'),
            axis.ticks = element_line(size = 1),
            legend.position = 'none')+
      update_geom_defaults("line", list(size = 5))+
      labs(x=paste(dat_list[[paste(exposureID,outcomeID,sep = "_")]]$exposure,paste0(' on ',dat_list[[paste(exposureID,outcomeID,sep = "_")]]$outcome),sep = "\n"))
    pdf(paste0("forest_",paste(exposureID,outcomeID,sep = "_"),".pdf"),width = 10,height = 15)
    print(pa)
    dev.off()
    
    res_single_2 <- mr_singlesnp(dat=dat_list[[paste(exposureID,outcomeID,sep = "_")]])
    nrow(res_single_2)
    res_single_2=res_single_2[-nrow(res_single_2),]
    pa <- mr_funnel_plot(res_single_2)[[1]]
    pa <- pa+scale_color_brewer(palette = 'Set1')+
      theme_classic(base_size = 18)+
      theme(panel.border = element_rect(size = 1.7,fill = 'transparent'),
            axis.ticks = element_line(size = 1),
            legend.position = 'top')+
      update_geom_defaults("line", list(size = 5))
    pdf(paste0("funnel_",paste(exposureID,outcomeID,sep = "_"),".pdf"),height = 8,width = 8)
    print(pa)
    dev.off()
  }
}

#### edgeR ####
a=read.table("GSE249311_count_matrix.tsv", header = T, sep = "\t")
a[1:2,1:4]
length(unique(a$gene))
exprSet=a
table(duplicated(exprSet$gene))

index <- order(rowMeans(exprSet[,-1]),decreasing = T)
exprSet_ordered <- exprSet[index,]
keep <- !duplicated(exprSet_ordered$gene)
exprSet_max <- exprSet_ordered[keep,]
table(duplicated(exprSet_max$gene))
row.names(exprSet_max)=exprSet_max$gene
exprSet_max$gene=NULL
max(exprSet_max)
min(exprSet_max)
write.table(exprSet_max, "GSE249311_counts.txt", quote = T,
            sep = "\t", col.names = NA)

juzhen1=read.delim("GSE249311_counts.txt",header = T,stringsAsFactors = F,check.names = F,row.names = 1)
fenzu=read.delim("group.GSE249311.txt",header = T,stringsAsFactors = FALSE,check.names = FALSE,row.names = 1)
juzhen1=juzhen1[,row.names(fenzu)]
all(colnames(juzhen1)==row.names(fenzu))
library(edgeR)
data=juzhen1
group=read.delim("group.GSE249311.txt",header = T,stringsAsFactors = FALSE,check.names = FALSE)
all(group$sample==colnames(data))
data=data[,group$sample]


dge <- DGEList(counts = data)
dge <- calcNormFactors(dge)
logCPM <- cpm(dge, log=TRUE, prior.count=3)

all(group$sample==colnames(data))
group=factor(group$group)
design<-model.matrix(~-1+group)
colnames(design)=levels(group)
row.names(design)=colnames(data)

fit <- lmFit(logCPM, design)
fit <- eBayes(fit)
contrasts<-makeContrasts(GDM-Control,levels=design)

fit2<-contrasts.fit(fit,contrasts)
fit2<-eBayes(fit2)
x1<-topTable(fit2,coef=1,number=nrow(data),adjust.method="BH",sort.by="B",resort.by="M")
write.csv(x1,"GSE249311_diff.matrix.csv")
write.csv(logCPM,"GSE249311_exp.cpm.csv")    
s=which(x1$adj.P.Val<0.05&abs(x1$logFC)>0.5)
x2=x1[s,]
write.csv(x2, file = "GSE249311.diff.adjp0.05.lfc0.5.csv")
a=which(row.names(x2)%in%junqun$GeneSymbol)
x3=x2[a,]
b=which(row.names(x2)%in%yundong$GeneSymbol)
x4=x2[b,]

library(ggVolcano)
deg_data=read.csv("GSE249311_diff.matrix.csv",header = T,stringsAsFactors = FALSE,check.names = FALSE,row.names = 1)
colnames(deg_data)
deg_data$label=""
gradual_volcano(deg_data, x = "logFC", y = "adj.P.Val",
                x_lab = "log2(FoldChange)",y_lab = "-log10(adj.P.Value)",legend_title = "-log10(adj.P.Value)",
                log2FC_cut = 0.5,FDR_cut = 0.05,
                add_label = TRUE,label = "label",label_number = nrow(deg_data),
                output = TRUE,filename = "volcano_plot"
)+labs(title="Tumor vs Normal")

library(pheatmap)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)
library(dendextend)
library(dendsort)
counts <- read.csv("GSE249311_exp.cpm.csv",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
deg <- read.table("GSE249311.diff.adjp0.05.lfc0.5.csv",sep = ",",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
logFC_cutoff=0.5
type1 = (deg$adj.P.Val<0.05)&(deg$logFC < -logFC_cutoff)
type2 = (deg$adj.P.Val<0.05)&(deg$logFC > logFC_cutoff)
deg$change = ifelse(type1,"Down",ifelse(type2,"Up","Not"))
table(deg$change)
deg$Label = ""  
deg <- deg[order(deg$adj.P.Val), ]   #
deg$Gene <- rownames(deg)
up.genes <- head(deg$Gene[which(deg$change == "Up")], 50)
down.genes <- head(deg$Gene[which(deg$change == "Down")], 50)
deg.top5.genes <- c(as.character(up.genes), as.character(down.genes))
deg$Label[match(deg.top5.genes, deg$Gene)] <- deg.top5.genes


deg <- deg[deg.top5.genes,]
deg$id <- rownames(deg)
ann_row=as.data.frame(deg[,c(7:8)])
ann_row <- ann_row[which(ann_row$change!="Not"),]
ann_row <- ann_row[order(ann_row$change,decreasing = T),]
ann_row1 <- as.data.frame(ann_row[,1])
rownames(ann_row1) <- rownames(ann_row)
colnames(ann_row1) <-"Change"
group<- read.delim("group.GSE249311.txt",check.names = F,stringsAsFactors = F,header = T)
ann_col <- group[order(group$group,decreasing = T),]
row.names(ann_col)=ann_col$sample
ann_col$sample=NULL
ann_col$id <- rownames(ann_col)
ann_col1 <- as.data.frame(ann_col[,1])
rownames(ann_col1) <- ann_col$id
colnames(ann_col1) <- "Group"
ha1 = HeatmapAnnotation(df = ann_col1,col = list(Group=c(Control="#0072b5",GDM="#bc3c29")))


data_pheat <- t(scale(t(as.matrix(counts[rownames(ann_row1),rownames(ann_col1)]))))
pdf("heatmap_plot.pdf", width = 6, height = 6)
densityHeatmap(data_pheat,top_annotation = ha1,title = "Distribution as heatmap", ylab = "Values") %v%
  Heatmap(data_pheat, col = colorRampPalette(c("navy", "white", "firebrick3"))(50),name = "Exp", 
          height = unit(8, "cm"),show_row_names = F,show_column_names=F,cluster_rows =T)
dev.off()


library(WGCNA)  
dataExpr=read.delim("GSE249311_exp.cpm.txt",header=T,row.names=1,stringsAsFactors = FALSE)
group=read.delim("group.wgcna.txt",header=T,row.names=1,stringsAsFactors = FALSE)
dataExpr=dataExpr[,row.names(group)]
all(colnames(dataExpr)==row.names(group))
dataExprVar=dataExpr
diff_exp=dataExprVar
mydata <- as.data.frame(t(dataExprVar))

nGenes = ncol(mydata)
nSamples = nrow(mydata)
mytree=hclust(dist(mydata), method="complete") 
pdf(file="1.sampleClustering.pdf",width = 20, height = 9)
par(mar = c(5,5,5,5))
plot(mytree, main = "Sample clustering", sub="", xlab="")
dev.off()

powers = seq(1,20)
sft = pickSoftThreshold(mydata, powerVector=powers, verbose=5)  
pdf(file = "2.softPower.pdf", width = 9, height = 5)
par(mfrow = c(1,2)) 
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylim = c(-1,1), ylab="Scale Free Topology Model Fit,signed R^2",type="n", main = "Scale independence")
abline(h=0.85,col="blue")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=0.6,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="n",main = "Mean connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=0.6,col="red")
dev.off()
softPower=sft$powerEstimate
adjacency=adjacency(mydata,power=softPower)
TOM=TOMsimilarity(adjacency)
save(TOM,file = "TOM.rda")
dissTOM=1-TOM
geneTree = hclust(as.dist(dissTOM), method = "average");
minModuleSize =70

dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize) #

dynamicColors = labels2colors(dynamicMods)

MEList = moduleEigengenes(mydata, colors = dynamicColors)
MEs = MEList$eigengenes 
MEDiss = 1-cor(MEs)
METree = hclust(as.dist(MEDiss), method = "average")
pdf(file = "4.eigengeneClustering.pdf", width = 12, height = 9)
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "") #
dev.off()

pdf("5.Eigengene_adjacency_heatmap.pdf")
plotEigengeneNetworks(MEs, 
                      "Eigengene adjacency heatmap", 
                      marHeatmap = c(3,4,2,2), 
                      plotDendrograms = FALSE, 
                      xLabelsAngle = 90)
dev.off()
MEDissThres=0.3
merge=mergeCloseModules(mydata,dynamicColors,cutHeight=MEDissThres,verbose=3)
mergedColors=merge$colors
mergedMEs=merge$newMEs
pdf("6.ModuleTree.pdf", width = 12, height = 9)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()

pdf("7.mergeModuleTree.pdf", width = 12, height = 9)
plotDendroAndColors(geneTree,mergedColors, "Mergeddynamic",dendroLabels=FALSE,hang=0.03,addGuide=TRUE,guideHang=0.05)
dev.off()
moduleColors = mergedColors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder) - 1
MEs = mergedMEs
write.table(paste(colnames(mydata), moduleColors, sep = "\t"), file = "4.netcolor2gene.xls", row.names=FALSE, col.names=FALSE, quote=FALSE)
gene_color=read.table("4.netcolor2gene.xls",row.names=1)
all(row.names(gene_color)==colnames(mydata))
color=unique(gene_color[,1])
color=color[color!="grey"]

traitColors = numbers2colors(group, signed = TRUE,centered=TRUE);
pdf("9.Sample_dendrogram_and_trait_heatmap.pdf",width = 20, height = 9)
plotDendroAndColors(mytree, 
                    traitColors, 
                    groupLabels = names(group), 
                    rowTextAlignment = "right-justified",
                    addTextGuide = TRUE ,
                    hang = 0.03,
                    dendroLabels = NULL, #
                    addGuide = FALSE,  #
                    guideHang = 0.05,
                    main = "Sample dendrogram and trait heatmap")
dev.off()
all(colnames(diff_exp)==row.names(group))
moduleTraitCor_noFP <- cor(mergedMEs, group, use = "p");
moduleTraitPvalue_noFP = corPvalueStudent(moduleTraitCor_noFP, nSamples); 
textMatrix_noFP <-paste(signif(moduleTraitCor_noFP, 2), "\n(", signif(moduleTraitPvalue_noFP, 1), ")", sep = ""); 
dim(textMatrix_noFP) = dim(moduleTraitCor_noFP)
pdf("10.Module-trait_relationships.pdf",width = 6, height = 8)
par(mar = c(4, 6, 5, 3)); 
labeledHeatmap(Matrix = moduleTraitCor_noFP, 
               xLabels = colnames(group), 
               yLabels = names(mergedMEs), 
               ySymbols = names(mergedMEs), 
               colorLabels = FALSE, 
               colors = blueWhiteRed(50), 
               textMatrix = textMatrix_noFP,
               setStdMargins = FALSE, 
               cex.text = 0.65, 
               zlim = c(-1,1), 
               main = paste("Module-trait relationships"))
dev.off()

#### GSVA ####
library(genefilter)
library(GSVA)
library(Biobase)
library(stringr)
uni_matrix=read.csv("GSE249311_exp.cpm.csv",header = T,stringsAsFactors = FALSE,check.names = FALSE,row.names = 1)
gene_immune=read.delim("gene.set.txt",header=T,stringsAsFactors = FALSE)##
list<- split(as.matrix(gene_immune)[,1], gene_immune[,2])
gsva_matrix<- gsva(as.matrix(uni_matrix), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
immune_score=as.data.frame(t(gsva_matrix))
phenotype_file=read.delim("group.GSE249311.txt",header = T,stringsAsFactors = FALSE,check.names = FALSE)
immune_score$sample=row.names(immune_score)
imm=merge(phenotype_file,immune_score,by="sample")
row.names(imm)=imm$sample
imm$sample=NULL
write.csv(imm,"GSE249311.score.csv")
data=imm
test=as.list(data[,-1])
height<-stack(test)
Group=rep(data$group,ncol(data)-1)     
df=as.data.frame(cbind(Group,height))
colnames(df)=c("group","Infiltration_level","cells")
library(ggplot2)
library(ggpubr)
library(ggsignif)
df$group=factor(df$group,levels=c("Lowrisk","Highrisk"))
pdf("immune.cells_boxplot.pdf",height=6,width=14)
ggplot(data=df,aes(x=cells,y=Infiltration_level,fill=group))+
  geom_boxplot(width=0.6)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_fill_manual(values = c("#4dbbd5", "#e64b35"))+
  stat_compare_means(aes(group=group),
                     label="p.signif",
                     method = "wilcox.test",
                     hide.ns = TRUE)+
  theme(axis.text.x = element_text(size=12,colour="black",angle=45,hjust = 1),axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=12),axis.title.y = element_text(size=12))

dev.off() 
#### venn ####
diff=read.csv("GSE249311_diff.matrix.csv",header = T,stringsAsFactors = F,check.names = F,row.names = 1)
diff1=diff[which(diff$adj.P.Val<0.05&abs(diff$logFC)>0.5),]
wgcna=read.delim("wgcna/result/hub.genes.txt",header = T,stringsAsFactors = F,check.names = F)
s=which(row.names(diff1)%in%wgcna$gene)
diff2=diff1[s,]
diff3=diff2[which(row.names(diff2)%in%junqun$GeneSymbol),]
diff4=diff2[which(row.names(diff2)%in%yundong$GeneSymbol),]
diff5=rbind(diff3,diff4)
#### venn ####
library(ggvenn)
gene<-list(LRGs=junqun$GeneSymbol,
  VPARGs=yundong$GeneSymbol,
           DEGs=row.names(diff1),
           WGCNA=wgcna$gene)
ggvenn(gene,fill_color=c("#0073C2FF","#EFC000FF","#CD534CFF","skyblue"),
       stroke_size=1, text_size = 4.5,
       set_name_size=6)
ggsave("venn_gene29.pdf")


library(clusterProfiler)
library(enrichplot)
library(tidyverse)
library(ggstatsplot)
library(ggnewscale)
library(org.Hs.eg.db)
library(DOSE)

diff=read.csv("diff29.csv",header = T,stringsAsFactors = FALSE,check.names = FALSE,row.names = 1)
gene_up=unique(row.names(diff))
gene_diff_entrez <- as.character(na.omit(bitr(gene_up, 
                                              fromType="SYMBOL", #
                                              toType="ENTREZID", # 
                                              OrgDb="org.Hs.eg.db")[,2])) 
kegg_enrich_results <- enrichKEGG(gene  = gene_diff_entrez,
                                  organism  = "hsa", #
                                  pvalueCutoff = 0.05,
                                  qvalueCutoff = 1,pAdjustMethod = "none")

kk_read <- DOSE::setReadable(kegg_enrich_results, 
                             OrgDb="org.Hs.eg.db", 
                             keyType='ENTREZID')#

a=kk_read@result
a=a[which(a$pvalue<0.05),]
write.csv(a,'KEGG_gene_enrichresults.csv')
save(kk_read, file ='KEGG_gene_enrichresults.Rdata')

####GO###
go_enrich_results <- enrichGO(gene = gene_diff_entrez,
                              OrgDb = "org.Hs.eg.db",
                              ont   = "MF"  ,     
                              pvalueCutoff  = 0.05,
                              qvalueCutoff  = 0.2,pAdjustMethod = "none",
                              readable      = TRUE)
b=go_enrich_results@result
b=b[which(b$pvalue<0.05),]
write.csv(b, 'GO_gene_MF_enrichresults.csv') 
save(go_enrich_results, file ='GO_gene_MF_enrichresults.Rdata')

pt <- pairwise_termsim(go_enrich_results)
treep <- treeplot(pt,color = "pvalue",
                  showCategory = 30)
ggsave(treep, filename = 'treeplot_MF.pdf', width=14, height=8)
diff <- read.csv("diff29.csv",header = T,stringsAsFactors = FALSE,check.names = FALSE,row.names = 1)
genelist=as.numeric(diff$logFC)
names(genelist) <- row.names(diff)
pdf("kegg_cnetplot.pdf",width = 12,height = 8)
cnetplot(kk_read,  
         foldChange = genelist,
         showCategory = 10,
         node_label = 'all',
         circular = T, 
         colorEdge = T)
dev.off()


library(tidyverse)
library(glmnet)
source('msvmRFE.R')   
library(e1071)
library(caret)
library(randomForest)

train <- read.csv("lasso_input.csv",row.names = 1, 
                  as.is = F) 
dim(train)
x <- as.matrix(train[,-1])  
y <- ifelse(train$group == "Control", 0,1) 
library(glmnet)
set.seed(123456)
fit = glmnet(x, y, family = "binomial", alpha = 1, lambda = NULL)
pdf("A_lasso_model.pdf", width = 5, height = 5)
plot(fit, xvar = "lambda", label = TRUE)
dev.off()
cvfit = cv.glmnet(x, y, 
                  nfold=10, 
                  family = "binomial", type.measure = "class")

pdf("A_lasso_cvfit1.pdf", width = 6, height = 4)
plot(cvfit)
dev.off()
cvfit$lambda.min
myCoefs <- coef(cvfit, s="lambda.min");
lasso_fea <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
(lasso_fea <- lasso_fea[-1])
write.csv(lasso_fea,"feature_lasso_0.2400862.csv")

library(randomForest)
df_train <- read.csv("lasso_input.csv",sep=",",header=T,check.names=F,row.names = 1)
set.seed(123456)
df_train$group <- factor(df_train$group)

train.forest <- randomForest(group ~ ., 
                             data = df_train, 
                             ntree = 500, # 
                             importance = TRUE,# 
                             proximity = TRUE) # 
train.forest
pdf("RF_trees.pdf",width = 5,height = 5)
plot(train.forest,main="Random Forest")   #
dev.off()

summary(train.forest)
importance_otu <- train.forest$importance
head(importance_otu)
importance_otu <- data.frame(importance(train.forest))

importance_otu <- importance_otu[order(importance_otu$MeanDecreaseAccuracy,decreasing = T),]
head(importance_otu)
write.csv(importance_otu,"RF_importance_otu.csv")

pdf("RF_top10.pdf",width = 6.5,height = 6.5)
varImpPlot(train.forest,n.var = min(10, nrow(train.forest$importance)), main = 'Top 10 - variable importance')
dev.off()
importance_otu=read.csv("RF_importance_otu.csv",header = T,stringsAsFactors = FALSE,check.names = FALSE,row.names = 1)
importance_otu$Feature=row.names(importance_otu)
importance_otu=importance_otu[1:10,]
colnames(importance_otu)
ggplot(importance_otu, aes(x= reorder( Feature,MeanDecreaseAccuracy), y=MeanDecreaseAccuracy,fill=Feature)) +
  geom_bar(stat="identity") +
  theme_classic() +
  guides(fill=FALSE)+
  scale_fill_manual(values=c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#D9D9D9","#BC80BD","#CCEBC5","#FFED6F","#377EB8","#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#D9D9D9","#BC80BD","#CCEBC5","#FFED6F","#377EB8"))+
  coord_flip()+
  theme_bw()+
  ggtitle('Top 10 - variable importance')+
  theme(plot.title = element_text(size=18,color='black'),##
        axis.title.x =element_text(size=16,color='black'),
        axis.text.x =element_text(size=14, color='black'),
        axis.title.y =element_blank(),
        axis.text.y=element_text(size=14,   color='black'),
        legend.title=element_text(size=16, color='black' ),
        legend.text=element_text(size=14, color='black'),
        title=element_text(size=20, color='black'),
        strip.text = element_text(size = 14))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  labs(x="gene",y="MeanDecreaseAccuracy",fill="")
ggsave('RF_importance_TOP10.pdf',w=10,h=10)

library(VennDiagram)
lasso=read.csv("feature_lasso_0.2400862.csv",header = T,stringsAsFactors = FALSE,check.names = FALSE,row.names = 1)
lasso_fea=lasso$x
rf=read.csv("RF_importance_otu.csv",header = T,stringsAsFactors = FALSE,check.names = FALSE,row.names = 1)
rf=rf[1:10,]
rf_fea=row.names(rf)
pdf("Venn_lasso_RF.pdf", width = 5, height = 3)
venn.plot <- venn.diagram(list(LASSO = lasso_fea, 
                               RF=rf_fea), NULL, 
                          fill = c("#E7B800",
                                   "#E31A1C"), 
                          cex = 2, cat.fontface=4, 
                          category.names = c("LASSO", 
                                             "RandomForest"), 
                          main = "Overlap")
grid.draw(venn.plot)
dev.off()
lasso_fea1=lasso_fea
s=which(lasso_fea1%in%rf_fea)
lasso_fea2=lasso_fea1[s]
write.csv(lasso_fea2,"venn_lasso_RF.csv")


library(foreign)
library(survival)
library(rms)
bc=read.csv("../lasso_svm_rf/lasso_input.csv",header = T,stringsAsFactors = FALSE,row.names = 1)
gene=read.csv("../lasso_svm_rf/venn_lasso_RF.csv",header = T,stringsAsFactors = F,check.names = F)
bc=bc[,c(1,which(colnames(bc)%in%gene$x))]
write.csv(bc,"train_gene2.csv")
bc$group=ifelse(bc$group == "Control", 0,1)
bc <- na.omit(bc)
dd <- datadist(bc)
options(datadist="dd")
colnames(bc)
formula1<-as.formula(group~ SRSF10+RPL11)
fit1<-lrm(formula1,data = bc,x=T,y=T)
summary(fit1)
nom1<-nomogram(fit1,
               fun=function(x)1/(1+exp(-x)),
               lp=F,
               fun.at = c(0.1,0.5,0.9),
              funlabel = "Risk of GDM")
pdf("nomogram.pdf",height = 4,width = 6)
plot(nom1)
dev.off() 

cal1<-calibrate(fit1,method = "boot",B=1000)
pdf("nomogram_predict.pdf")
plot(cal1,xlim=c(0,1.0),ylim=c(0,1.0),
     xlab = "Nomogram Predicted Probability", ylab = "Actual Probability")
dev.off()

library(DecisionCurve)
nom <- nomogram(fit1)
library(nomogramFormula)
results<-formula_rd(nomogram=nom)
bc$points<-points_cal(formula = results$formula,rd=bc)#
pre<-bc$points
library(pROC)
plot.roc(bc$group, pre,
         main="Nomogram ROC Curve", percent=F,
         print.auc=T,
         ci=F, of="thresholds",
         thresholds="best",print.thres=FALSE,col="blue")

library(corrplot)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(psych)
uniq.symbol_sva=read.csv("../GSE249311_exp.cpm.csv",header = T,stringsAsFactors = FALSE,check.names = FALSE,row.names = 1)
group=read.delim("../group.GSE249311.txt",header = T,stringsAsFactors = FALSE,check.names = FALSE)
group1=group[which(group$group=="GDM"),]
data_FPKM=uniq.symbol_sva[,group1$sample]
rt <- data_FPKM

for (ii in c( "SRSF10",
"RPL11")){  
  tar.exp <- rt[ii,]
  y <- as.numeric(tar.exp)
  data1 <- data.frame()
  for (i in rownames(rt)) {
    dd  <- corr.test(as.numeric(rt[i,]), y , method="pearson",adjust = "fdr")
    data1 = rbind(data1,data.frame(gene=i,cor=dd$r,p.value=dd$p))
  }
  data1 <- data1[order(data1$cor,decreasing = T),]
  write.csv(data1,file = paste0(ii,".pearson.cor.csv"))
}
  ii="RPL11"
  data1=read.csv(file = paste0(ii,".pearson.cor.csv"),header = T,stringsAsFactors = FALSE,check.names = FALSE)
  gene_cor=data1
  geneList <- gene_cor$cor
  names(geneList) = gene_cor$gene
  geneList = sort(geneList, decreasing = TRUE)
  head(geneList)
  tail(geneList)
  
  library(clusterProfiler)
  marks <- read.gmt("c5.go.bp.v7.4.symbols.gmt")
  y <- GSEA(geneList,TERM2GENE =marks)
  save(y,file=paste0(ii,"-GSEA-GOBP.rda"))
  write.table(y,file=paste0(ii,"-GSEA-GOBP.txt"),sep="\t",quote=F,row.names = F)   

  library(GseaVis)
  library(ggplot2)
  pdf(file = paste0(ii, "-GSEA-GOBP-TOP10.pdf"),width = 10, height = 8)
  gseaNb(y, geneSetID = y@result$Description[1:10],
         curveCol=c("#0000FF","#00FFFF","#04FF8A","#7CFC00","#AEFFB3",
                    "#F5DEB3","#F4A460","#FF6803","#FF6969","#FF0000"),
         pvalX = 0.9,pvalY = 0.9,subPlot = 3,rmHt = T)
  dev.off()
  
  
  marks <- read.gmt("c2.cp.kegg.v7.4.symbols.gmt")
  y1 <- GSEA(geneList,TERM2GENE =marks)
  save(y1,file=paste0(ii,"-GSEA-KEGG.rda"))
  write.table(y1,file=paste0(ii,"-GSEA-KEGG.txt"),sep="\t",quote=F,row.names = F)
  
  
  pdf(file = paste0(ii, "-GSEA-KEGG-TOP10.pdf"),width = 10, height = 8)
  gseaNb(y1, geneSetID = y1@result$Description[1:10],
         curveCol=c("#0000FF","#00FFFF","#04FF8A","#7CFC00","#AEFFB3",
                    "#F5DEB3","#F4A460","#FF6803","#FF6969","#FF0000"),
         pvalX = 0.9,pvalY = 0.9,subPlot = 3,rmHt = T)
  dev.off()


library(dplyr)
library(ggplot2)
library(magrittr)
library(gtools)
library(stringr)
library(Matrix)
library(tidyverse)
library(patchwork)
library(data.table)
library(RColorBrewer)
library(ggpubr)
suppressMessages(library(Seurat))
suppressMessages(library(SingleR))
suppressMessages(library(Matrix))
suppressMessages(library(dplyr))
suppressMessages(library(sctransform))
suppressMessages(library(stringr))
suppressMessages(library(progress))
suppressMessages(library(edgeR))
suppressMessages(library(BiocParallel))
suppressMessages(library(scater))
ii <- c('GSM5261695',
'GSM5261696',
'GSM5261697',
'GSM5261698')

seurat_obj <- list()

for (i in ii) {
  

  data.dir = i
  barcode.path <- paste0(data.dir, "/barcodes.tsv.gz")
  genes.path <- paste0(data.dir, "/features.tsv.gz")
  matrix.path <- paste0(data.dir, "/matrix.mtx.gz")


  zebrafish.data <- readMM(file = matrix.path)  

  gene.names = read.delim(genes.path, header = FALSE, stringsAsFactors = FALSE)
  barcode.names = read.delim(barcode.path, header = FALSE, stringsAsFactors = FALSE)
  colnames(zebrafish.data) = barcode.names$V1
  rownames(zebrafish.data) = gene.names$V2  
  dim(zebrafish.data)  

  write.table(zebrafish.data, paste0(i, "_data.txt"), quote = T,
             sep = "\t", col.names = NA)
  
}

if (!file.exists(paste0("1.CreatSeuratObject/"))) {
  dir.create("1.CreatSeuratObject")
}

samplename <- c('GSM5261695',
                'GSM5261696',
                'GSM5261697',
                'GSM5261698')

path <- "1.CreatSeuratObject"

pb <- progress_bar$new(format = "  完成百分比 [:bar] :percent 执行时间 :elapsed",
                       total = length(filelist), clear = FALSE, width = 60)
dimlist <- list()
object.list <- list()

for (i in 1:length(samplename)) {
  #
  file = paste0(samplename[i],"_data.txt")
  #
  if (!file.exists(paste0(path, "/1.QC/"))) {
    dir.create(paste0(path, "/1.QC/"))
  }
  if (!file.exists(paste0(path, "/1.QC/", samplename[i]))) {
    dir.create(paste0(path, "/1.QC/", samplename[i]))
  }
  message("")
  message("read file: ", samplename[i])
  # read data-read expression matrix-BD expression matrix csv file
  data <- read.delim(file, header = T, row.names = 1, check.names = F,stringsAsFactors = FALSE)
  message("creat seurat object: ", samplename[i])
  object <- CreateSeuratObject(counts = data, project = samplename[i], min.cells = 3,
                               min.features = 200)
  #object$orig.ident=samplename[i]
  #object@active.ident=factor(samplename[i])
  dimlist[[samplename[i]]] <- c(dim(data)[2], dim(object)[2])
  message("dim of object:", dim(object)[1], " ", dim(object)[2])
  message("QC: ", samplename[i])

  object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^MT-")
  print(summary(object[["percent.mt"]]))
  mt.percnet_summary <- summary(object[["percent.mt"]])
  capture.output(mt.percnet_summary, file = paste0(path, "/1.QC/", samplename[i],
                                                   "/mt.percent_statitics.txt"))
  
  HB.genes_total <- c("HBA1", "HBA2", "HBB", "HBD", "HBG1", "HBQ1")
  HB_m <- match(HB.genes_total, rownames(object@assays$RNA))
  HB.genes <- rownames(object@assays$RNA)[HB_m]
  HB.genes <- HB.genes[!is.na(HB.genes)]
  object[["percent.HB"]] <- PercentageFeatureSet(object, features = HB.genes)
  # QC-statistics of mt.percent per cell
  print(summary(object[["percent.HB"]]))
  HB.percnet_summary <- summary(object[["percent.HB"]])
  capture.output(HB.percnet_summary, file = paste0(path, "/1.QC/",samplename[i], "/HB.percent_statitics.txt"))
  
  summary(colSums(object))
  UMI_counts_summary <- summary(colSums(object))
  capture.output(UMI_counts_summary, file = paste0(path, "/1.QC/", samplename[i],
                                                   "/UMI_counts_statitics.txt"))
  UMI_counts_summary.show <- data.frame(Quantile=names(UMI_counts_summary),summary=as.numeric(UMI_counts_summary))
  write.csv(UMI_counts_summary.show,file=paste0(path,'/1.QC/',samplename[i],'_细胞文库大小统计.csv'))
  
  
  p <- hist(colSums(object), breaks = 200)
  pdf(file = paste0(path, "/1.QC/", samplename[i], "/Distribution of UMI counts.pdf"))
  # QC-histogram of UMI counts distribution
  hist(colSums(object), breaks = 200, main = "Distribution of UMI counts", xlab = "UMI counts",
       ylab = "Cell number", col = "#2b8cbe", xlim = c(0, max(p$breaks) * 1.1), ylim = c(0,
                                                                                         max(p$counts) * 1.1))
  dev.off()
  
  gene_counts_per_cell <- as.numeric(object$nFeature_RNA)
  message("QC-calculate the gene number of each cell")
  print(summary(gene_counts_per_cell))
  gene_counts_summary <- summary(gene_counts_per_cell)
  capture.output(gene_counts_summary, file = paste0(path, "/1.QC/", samplename[i],
                                                    "/gene_counts_statitics.txt"))
  gene_counts_summary.show <- data.frame(Quantile=names(gene_counts_summary),summary=as.numeric(gene_counts_summary))
  write.csv(gene_counts_summary.show,file=paste0(path,'/1.QC/',samplename[i],'_每个细胞基因数统计.csv'))											
  
  
  
  p <- hist(gene_counts_per_cell, breaks = 200)
  pdf(file = paste0(path, "/1.QC/", samplename[i], "/Distribution of gene counts.pdf"))
  # QC-histogram of gene counts distribution
  hist(gene_counts_per_cell, breaks = 200, main = "Distribution of gene counts",
       xlab = "Gene counts", ylab = "Cell number", col = "#2b8cbe", xlim = c(0, max(p$breaks) *
                                                                               1.1), ylim = c(0, max(p$counts) * 1.1))
  dev.off()

  p <- VlnPlot(object, col = "#2b8cbe", features = c("nFeature_RNA", "nCount_RNA",
                                                     "percent.mt"
                                                     
  ), ncol = 3, pt.size = 0.01)
  #print(p)
  pdf(file = paste0(path, "/1.QC/", samplename[i], "/Distribution of UMI counts, gene counts, mitochondria pct.pdf"))
  print(p)
  dev.off()
  
  p.hb <- VlnPlot(object, col = "#2b8cbe", features = c("nFeature_RNA", "nCount_RNA",
                                                        "percent.HB"
                                                        
  ), ncol = 3, pt.size = 0.01)
  #print(p)
  pdf(file = paste0(path, "/1.QC/", samplename[i], "/Distribution of UMI counts, gene counts, HB pct.pdf"))
  print(p.hb)
  dev.off() 
  plot1 <- FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
                          cols = "black", pt.size = 0.01) + NoLegend()
  plot2 <- FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "percent.mt",
                          cols = "black", pt.size = 0.01) + NoLegend()
  plot4 <- FeatureScatter(object, feature1 = "nFeature_RNA", feature2 = "percent.mt",
                          cols = "black", pt.size = 0.01) + NoLegend()
  pdf(file = paste0(path, "/1.QC/", samplename[i], "/Distribution of UMI counts, gene counts, mitochondria pct, in feature scatter.pdf"))
  p <- CombinePlots(plots = list(plot1,plot2, plot4), ncol = 3, legend = "none")
  print(p)
  dev.off()
  
  
  plot3 <- FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
                          cols = "black", pt.size = 0.01) + NoLegend()
  plot5 <- FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "percent.HB",
                          cols = "black", pt.size = 0.01) + NoLegend()
  plot6 <- FeatureScatter(object, feature1 = "nFeature_RNA", feature2 = "percent.HB",
                          cols = "black", pt.size = 0.01) + NoLegend()
  pdf(file = paste0(path, "/1.QC/", samplename[i], "/Distribution of UMI counts, gene counts, HB pct in feature scatter.pdf"))
  p <- CombinePlots(plots = list(plot3,plot5, plot6), ncol = 3, legend = "none")
  print(p)
  dev.off()

  libsize.drop <- isOutlier(object$nCount_RNA, nmads = 2, type = "both", log = TRUE)
  feature.drop <- isOutlier(object$nFeature_RNA, nmads = 2, type = "both", log = TRUE)
  mt.drop <- isOutlier(object$percent.mt, nmads = 2, type = "higher")
  HB.drop <- isOutlier(object$percent.HB, nmads = 2, type = "higher")
  print(attr(libsize.drop, "thresholds"))
  print(attr(feature.drop, "thresholds"))
  print(attr(mt.drop, "thresholds"))
  print(attr(HB.drop, "thresholds"))
  keep <- !(libsize.drop | feature.drop | mt.drop|HB.drop)
  filter=data.frame(ByLibSize = sum(libsize.drop), ByFeature = sum(feature.drop), ByMt = sum(mt.drop),ByHB = sum(HB.drop), 
                    Remaining = sum(keep))
  write.csv(filter, file = paste0("1.CreatSeuratObject", "/1.QC/", samplename[i],
                                  "/细胞过滤情况.csv"))
  
  message("Subset seurat object: ", samplename[i])
  object.list[[samplename[i]]] <- object[, keep]
  
  message("dim of subset object:", dim(object.list[[samplename[i]]])[1], dim(object.list[[samplename[i]]])[2])
 
  lost <- calculateAverage(object@assays$RNA@counts[, !keep])
  kept <- calculateAverage(object@assays$RNA@counts[, keep])
  capped.lost <- pmax(lost, min(lost[lost > 0]))
  capped.kept <- pmax(kept, min(kept[kept > 0]))
  is.mito <- grepl(pattern = "^MT-", rownames(object))
  
  pdf(file = paste0(path, "/1.QC/", samplename[i], "/Average count of discarded and retained cells.pdf"))
  plot(capped.lost, capped.kept, xlab = "Average count (discarded)", ylab = "Average count (retained)", 
       log = "xy", pch = 16) + points(capped.lost[is.mito], capped.kept[is.mito], 
                                      col = "#2b8cbe", pch = 16)
  dev.off()
  y <- cbind(lost, kept)
  y <- DGEList(y, group = c(2, 1))
  design <- model.matrix(~group, data = y$samples)
  predlfc <- predFC(y, design)[, 2]
  info <- data.frame(logFC = predlfc, Lost = lost, Kept = kept, diff=lost-kept,row.names = rownames(object))
  message("top 20 genes in discarded cells :")
  print(head(info[order(info$diff, decreasing = TRUE), ], 20))
  
  pdf(file = paste0(path, "/1.QC/", samplename[i], "/Distribution of UMI counts, gene counts, mitochondria pct, HB pct for newseurat.pdf"))
  p <- VlnPlot(object.list[[samplename[i]]], col = "#2b8cbe", features = c("nFeature_RNA", 
                                                                           "nCount_RNA"), ncol = 2, pt.size = 0.01)
  print(p)
  dev.off()
  #print(p)
  
  message("creat seurat object ", samplename[i], " finished!")
  Sys.sleep(1/1000)
}


dimlist <- do.call(rbind, dimlist)
save(object.list, dimlist, file = "object.list.rda")
sample.cell=data.frame(dimlist)
sample.cell$sample=row.names(sample.cell)
group=read.delim("单细胞fenzu.txt",header = T,stringsAsFactors = FALSE,check.names = FALSE)
sample.cell1=merge(group,sample.cell,by="sample")
write.csv(sample.cell1,"sample.cell.size.csv")


if (!file.exists(paste0("2.integrate_and_cluster/"))) {
  dir.create("2.integrate_and_cluster")
}

object.list <- lapply(X = object.list, FUN = function(x) {
  samplename <- as.character(x@meta.data$orig.ident[1])
  if (grepl("-1", colnames(x@assays$RNA@counts)[1])) {
    x <- RenameCells(x, new.names = str_replace_all(
      colnames(x@assays$RNA@counts),
      "-1", paste0("-", samplename)
    ))
  } else {
    x <- RenameCells(x, new.names = paste0(colnames(x@assays$RNA@counts), paste0(
      "-",
      samplename
    )))
  }
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000, mean.cutoff = c(
    0.1,
    10
  ), dispersion.cutoff = c(1.25, 20))
})

features <- SelectIntegrationFeatures(object.list = object.list)
ifnb.anchors <- FindIntegrationAnchors(object.list = object.list, anchor.features = features)
save(ifnb.anchors,file=paste0('Step3.ifnb.anchors.rda'))
for (i in 1:length(object.list)) {
  object.subset <- object.list[[i]]
  samplename <- as.character(object.subset@meta.data$orig.ident[1])
  message("visualize high varible genes: ", samplename)
  if (!file.exists(paste0("2.integrate_and_cluster/", "1.visualize_high_varible_genes"))) {
    dir.create(paste0("2.integrate_and_cluster/", "1.visualize_high_varible_genes"))
  }
  top20 <- head(VariableFeatures(object.subset), 20)
  plot6 <- VariableFeaturePlot(object.subset, pt.size = 2, cols = c("black", "#2b8cbe"))
  plot7 <- LabelPoints(plot = plot6, points = top20, repel = TRUE) + NoLegend()
  pdf(paste0(
    "2.integrate_and_cluster/", "1.visualize_high_varible_genes/", samplename,
    " visualize high varible genes.pdf"
  ), width = 16, height = 8)
  p <- CombinePlots(plots = list(plot6, plot7, legend = "bottom"))#, legend = "bottom"
  print(p)
  dev.off()

  Sys.sleep(1 / 1000)
}

sum(sample.cell$X2)

PBMC.combined <- IntegrateData(anchorset = ifnb.anchors)
DefaultAssay(PBMC.combined) <- "integrated"
save(PBMC.combined, file = paste0("Step3.PBMC.combined.rda"))
dim(PBMC.combined)
PBMC.combined <- ScaleData(PBMC.combined, verbose = FALSE, vars.to.regress = c(
  "nCount_RNA",
  "nFeature_RNA"
))

if (!file.exists(paste0("2.integrate_and_cluster/", "2.PCA_analysis"))) {
  dir.create(paste0("2.integrate_and_cluster/", "2.PCA_analysis"))
}

PBMC.combined <- RunPCA(PBMC.combined,
                        npcs = 50, verbose = FALSE, do.print = TRUE,
                        pcs.print = 1:5, genes.print = 5
)

pdf(file="2.integrate_and_cluster/2.PCA_analysis/RunPCA.pdf",width=7,height=7)
DimPlot(object = PBMC.combined, reduction = "pca")
dev.off()
# Examine and visualize PCA results with printing the major contributors
major.PC.related.gene <- print(PBMC.combined[["pca"]], dims = 1:50, nfeatures = 5)

pdf(file = paste0("2.integrate_and_cluster/", "2.PCA_analysis", "/top9 PCs related genes.pdf"),
    width = 14, height = 14)
VizDimLoadings(PBMC.combined, dims = 1:9, reduction = "pca")
dev.off()

pdf(file = paste0("2.integrate_and_cluster/", "2.PCA_analysis",  "/top9 PCs related genes pheatmap.pdf"), 
    width = 14, height = 14)
DimHeatmap(PBMC.combined, dims = 1:9, cells = 100, reduction = "pca", ncol = 3, balanced = TRUE)
dev.off()



pdf(
  file = paste0("2.integrate_and_cluster/", "2.PCA_analysis", "/Examine and visualize PCA results with ElbowPlot.pdf"),
  width = 9, height = 9
)
ElbowPlot(PBMC.combined, ndims = 50, reduction = "pca")
dev.off()
save(PBMC.combined, file = paste0("Step4.PBMC.combined.rda"))

suppressMessages(library(Seurat))
suppressMessages(library(SingleR))
suppressMessages(library(Matrix))
suppressMessages(library(dplyr))
suppressMessages(library(sctransform))
suppressMessages(library(scater))
suppressMessages(library(stringr))
suppressMessages(library(edgeR))
suppressMessages(library(ggplot2))

suppressMessages(library(progress))
suppressMessages(library(cluster))
suppressMessages(library(dynamicTreeCut))
suppressMessages(library(tidytree))
suppressMessages(library(ggtree))
suppressMessages(library(clustree))

if (!file.exists(paste0("2.integrate_and_cluster/", "3.cluster_analysis"))) {
  dir.create(paste0("2.integrate_and_cluster/", "3.cluster_analysis"))
}

PBMC.combined <- FindNeighbors(PBMC.combined, reduction = "pca", dims = 1:30)
PBMC.combined <- FindClusters(object = PBMC.combined, reduction = "pca", resolution = c(1:10) / 10,
                              dims = 1:30, save.SNN = TRUE)

pdf(file = paste0("2.integrate_and_cluster/", "2.PCA_analysis", "/03_Cells_cluster.pdf"), width = 40, height = 6)
p <- DimPlot(PBMC.combined, reduction = "pca", split.by = "orig.ident", label = TRUE, repel = FALSE)
print(p)
dev.off()
save(PBMC.combined, file = "Step5.FindNeighbors.resolution.PBMC.combined.Rdata")

n <- ncol(PBMC.combined)
set.seed(2023)
n <- ifelse(n > ncol(PBMC.combined), ncol(PBMC.combined), n)
message("number of cell sample: ", n)

sample <- sample(1:ncol(PBMC.combined), n)
pcs <- PBMC.combined@reductions$pca@cell.embeddings[sample, 1:30]

resolution <- 0.5
message("resolution :", resolution)
integrated_snn_res <- paste0("integrated_snn_res.", resolution)

Idents(PBMC.combined) <- PBMC.combined@meta.data[[integrated_snn_res]]
PBMC.combined$seurat_clusters <- PBMC.combined@meta.data[[integrated_snn_res]]

library(RColorBrewer)
mypalette <- colorRampPalette(brewer.pal(8,"Set1"))
library("ggsci")

# run non linear dimensionality reduction-UMAP
PBMC.combined <- RunUMAP(PBMC.combined, reduction = "pca", dims = 1:30)
#DimPlot(PBMC.combined, reduction = "umap", label = TRUE)

pdf(file = paste0("2.integrate_and_cluster/", "3.cluster_analysis", "/3-cluster with  umap.pdf"),
    width = 7, height = 6)
DimPlot(PBMC.combined, reduction = "umap", label = TRUE,cols = mypalette(28),pt.size = 0.5)
dev.off()
PBMC.combined <- RunTSNE(PBMC.combined, dims = 1:30, perplexity = 20)


pdf(file = paste0("2.integrate_and_cluster/", "3.cluster_analysis", "/3-cluster with  tsne.pdf"),
    width = 7, height = 6)
DimPlot(PBMC.combined, reduction = "tsne", label = TRUE,cols = mypalette(28),pt.size = 0.5)
dev.off()
table(PBMC.combined$seurat_clusters)

p1 <- DimPlot(PBMC.combined, reduction = "umap", label = F, pt.size = 0.5,group.by = "orig.ident",cols = mypalette(28))
p2 <- DimPlot(PBMC.combined, reduction = "umap", pt.size = 0.5,label = TRUE,cols = mypalette(28))

pdf(file = paste0("2.integrate_and_cluster/", "3.cluster_analysis", "/integrate and cluster with  umap.pdf"),
    width = 15, height = 6)
print(CombinePlots(plots = list(p1, p2)))
dev.off()

# Visualization
p1 <- DimPlot(PBMC.combined, reduction = "tsne", label = F, pt.size = 0.5, group.by = "orig.ident",cols = mypalette(28))
p2 <- DimPlot(PBMC.combined, reduction = "tsne", pt.size = 0.5, label = TRUE,cols = mypalette(28))
pdf(file = paste0("2.integrate_and_cluster/", "3.cluster_analysis", "/integrate and cluster with  tsne.pdf"),
    width = 15, height = 6)
print(CombinePlots(plots = list(p1, p2)))
dev.off()

pdf(file = paste0("2.integrate_and_cluster/", "3.cluster_analysis", "/integrate and cluster with  tsne individually.pdf"),
    width = 10, height = 4)
print(DimPlot(PBMC.combined, reduction = "tsne", split.by = "orig.ident", ncol = 4,pt.size = 1,cols = mypalette(28)))
dev.off()

pdf(file = paste0("2.integrate_and_cluster/", "3.cluster_analysis", "/integrate and cluster with  umap individually.pdf"),
    width = 10, height = 4)
print(DimPlot(PBMC.combined, reduction = "umap", split.by = "orig.ident", ncol = 4,pt.size = 1,cols = mypalette(28)))
dev.off()

save(PBMC.combined,file = '04.PBMC.combined.Rdata')

if (!file.exists(paste0("2.integrate_and_cluster/", "4.cluster_markers"))) {
  dir.create(paste0("2.integrate_and_cluster/", "4.cluster_markers"))
}

PBMC.combined.all.markers <- FindAllMarkers(PBMC.combined, only.pos = TRUE, min.pct = 0.2,
                                            logfc.threshold = 0.5, test.use = "wilcox")
PBMC.combined.all.markers %>%
  group_by(cluster) %>%
  top_n(n = 2, wt = avg_log2FC)

# write the marker genes to the txt file
write.table(PBMC.combined.all.markers, file = paste0("2.integrate_and_cluster/",
                                                     "4.cluster_markers", "/markers_allcluster.txt"), row.names = T, quote = FALSE,
            sep = "\t")
library(RColorBrewer)
mypalette <- colorRampPalette(brewer.pal(8,"Set1"))
genes=c()
var=unique(PBMC.combined.all.markers$cluster)
for (i in 1:length(var)) {
  
  PBMC.combined.all.markers1=PBMC.combined.all.markers[which(PBMC.combined.all.markers$cluster==var[i]),]
  fea=PBMC.combined.all.markers1$gene[1:5]
  genes=c(genes,fea)
}

pdf(file="2.integrate_and_cluster/4.cluster_markers/markers-_heatmap_new.pdf",width=18,height=12)
DoHeatmap(object = PBMC.combined, features =genes,group.by = "seurat_clusters",group.colors = mypalette(28))+scale_fill_gradientn(colors = c("navy","white","firebrick3"))
dev.off()

library(scRNAtoolVis)
library(ggsci)
library(RColorBrewer)
mypalette <- colorRampPalette(brewer.pal(8,"Set1"))
scRNA=PBMC.combined

pdf(file="2.integrate_and_cluster/4.cluster_markers/markers_AverageHeatmap.lfc.pdf",width=15,height=18)
averageHeatmap(object = scRNA,
               markerGene = genes,group.by = 'seurat_clusters')
dev.off()
load("04.PBMC.combined.Rdata")
DefaultAssay(PBMC.combined) <- "RNA"
dim(PBMC.combined)

VCT=c('CDH1','MET','PARP1')
SCT=c('CYP19A1', 'CGA', 'LGALS13')
EVT=c( 'PAPPA2', 'MMP2', 'MMP11')
Granulocyte=c('FCGR3B', 'CXCL8', 'SELL')
Myelocyte=c('TCN1', 'CEACAM8', 'CAMP')
NK_cell=c('GZMA', 'GZMK')
T_cell=c('CD3G', 'CD3D')
B_cell=c('CD79A', 'CD79B')
Monocytes=c('CD14', 'CD300E')
Macrophages=c('CD14', 'CD163', 'CD209')
Fibroblasts=c('DCN','LUM')
Mast_cell=c('CPA3', 'HDC', 'KIT')
features <- unique(c(VCT,SCT,EVT,Granulocyte,Myelocyte,T_cell,NK_cell,B_cell,
                     Monocytes,Macrophages,Fibroblasts,Mast_cell))

dim(PBMC.combined)
pdf(file = "2.integrate_and_cluster/4.cluster_markers/05_markers-DotPlot_final.pdf",width = 12,height = 6)
print(DotPlot(PBMC.combined, features=features,
              cols = c("lightgrey", "red")) + RotatedAxis())

dev.off()
PBMC.combined@meta.data$seurat_clusters

if (!file.exists(paste0("2.integrate_and_cluster/", "5.cluster_named"))) {
  dir.create(paste0("2.integrate_and_cluster/", "5.cluster_named"))
}
table(PBMC.combined@meta.data$orig.ident)
PBMC.combined@meta.data<-PBMC.combined@meta.data %>%
  mutate(Group= if_else(orig.ident == "GSM5261695"|orig.ident == "GSM5261696",
                        true = 'Control',
                        false ='GDM'))

table(PBMC.combined$Group)
table(PBMC.combined$seurat_clusters)

new.cluster.ids <- c('Macrophages',
'Granulocyte',
'VCT',
'VCT',
'EVT',
'T cell',
'Myelocyte',
'EVT',
'Monocytes',
'VCT',
'NK cell',
'Myelocyte',
'VCT',
'VCT',
'SCT',
'VCT',
'Fibroblasts',
'Macrophages',
'Mast cell',
'Macrophages',
'EVT',
'Fibroblasts',
'B cell',
'Mast cell',
'Fibroblasts')
names(new.cluster.ids) <- levels(PBMC.combined)
PBMC.combined <- RenameIdents(PBMC.combined,new.cluster.ids)
PBMC.combined$celltype <- Idents(PBMC.combined) 
table(PBMC.combined$celltype)

pdf('2.integrate_and_cluster/5.cluster_named/05-细胞类型注释-umap.pdf',width = 8,height = 6)
DimPlot(PBMC.combined, reduction = "umap", label = TRUE, pt.size = 0.5,cols = mypalette(28)) 
dev.off()

pdf("2.integrate_and_cluster/5.cluster_named/05-细胞类型注释-sample-umap.pdf",width = 9,height = 5)
DimPlot(PBMC.combined, reduction = "umap",split.by = 'Group', label = TRUE, pt.size = 0.5,cols = mypalette(28))
dev.off()

pdf("2.integrate_and_cluster/5.cluster_named/05-细胞类型注释-tsne.pdf",width = 8,height = 6)
DimPlot(PBMC.combined, reduction = "tsne", label = TRUE, pt.size = 0.5,cols = mypalette(28))
dev.off()

pdf("2.integrate_and_cluster/5.cluster_named/05-细胞类型注释-sample-tsne.pdf",width = 9,height = 5)
DimPlot(PBMC.combined, reduction = "tsne",split.by = 'Group', label = TRUE, pt.size = 0.5,cols = mypalette(28))
dev.off()

dim(PBMC.combined@meta.data)
dim(PBMC.combined)
DefaultAssay(PBMC.combined) <- "integrated"
save(PBMC.combined, file = "05_PBMC.combined.Rdata")#

features <- unique(c(B_cell,Mast_cell,Fibroblasts,SCT,NK_cell,Monocytes,Myelocyte,
                     T_cell,EVT,VCT,Granulocyte,
                     Macrophages))


dim(PBMC.combined@meta.data)
DefaultAssay(PBMC.combined) <- "RNA"
pdf(file = "2.integrate_and_cluster/5.cluster_named/05-细胞类型注释_markers-DotPlot.pdf",width = 9,height = 4.5)
print(DotPlot(PBMC.combined, features=features,
              cols = c("lightgrey", "red")) + RotatedAxis())

dev.off()

pdf(file="2.integrate_and_cluster/5.cluster_named/05-细胞类型注释_markers-pheatmap.pdf",width=10,height=8)
mypalette <- colorRampPalette(brewer.pal(8,"Set1"))
DoHeatmap(object = PBMC.combined, features =features,group.by = "celltype",slot = "counts",group.colors = mypalette(28))+scale_fill_gradientn(colors = c("navy","white","firebrick3"))
dev.off()

pdf(file="2.integrate_and_cluster/5.cluster_named/05-markers_AverageHeatmap.lfc.pdf",width=6,height=6)
averageHeatmap(object = PBMC.combined,
               markerGene = features,group.by = 'celltype')
dev.off()

library(RColorBrewer)
get_palette=colorRampPalette(brewer.pal(6,"Set1"))
Cell_type_stat_table <- table(PBMC.combined@meta.data[, c("orig.ident", "celltype")])
write.csv(Cell_type_stat_table, file = "2.integrate_and_cluster/5.cluster_named/05_细胞数量统计.csv")
sum(Cell_type_stat_table)
Cell_type_stat_table = Cell_type_stat_table/sum(Cell_type_stat_table)
write.csv(Cell_type_stat_table, file = "2.integrate_and_cluster/5.cluster_named/05_细胞比例统计.csv")
library(ggplot2)

library(RColorBrewer)
library(ggsci)
get_palette=colorRampPalette(brewer.pal(6,"Set1"))
cols=get_palette(30)
bili=as.data.frame(table(PBMC.combined@meta.data$celltype,PBMC.combined@meta.data$Group))
colnames(bili)=c("celltype","sample","number")
pdf(file = "2.integrate_and_cluster/5.cluster_named/05-celltype_sample.pdf",width = 6,height = 5)
ggplot(data=bili,aes(y=celltype,x=number,fill=sample))+
  geom_bar(stat = "identity",position = "fill",width = 0.6)+theme_classic()+
  theme(axis.text.x = element_text(size = 15),axis.title.y = element_text(size = 15),axis.title = element_text(size = 15))+scale_fill_jco()
dev.off()
##
library(RColorBrewer)
get_palette=colorRampPalette(brewer.pal(6,"Set1"))
cols=get_palette(30)
bili=as.data.frame(table(PBMC.combined@meta.data$celltype,PBMC.combined@meta.data$Group))
colnames(bili)=c("celltype","sample","number")
pdf(file = "2.integrate_and_cluster/5.cluster_named/05-sample_celltype.pdf",width = 6,height = 5)
ggplot(data=bili,aes(y=sample,x=number,fill=celltype))+
  geom_bar(stat = "identity",position = "fill",width = 0.4)+theme_classic()+
  theme(axis.text.x = element_text(size = 15),axis.title.y = element_text(size = 15),axis.title = element_text(size = 15))#+scale_fill_jco()
dev.off()



library(AUCell)
library(Seurat)
gene=c("SRSF10",
"RPL11")
load("05_PBMC.combined.Rdata")
sce2=PBMC.combined
DefaultAssay(sce2) <- "RNA"
dim(sce2)
s=which(row.names(sce2)%in%gene)
HYPOXIA_AUC=gene
cells_rankings <- AUCell_buildRankings(sce2@assays$RNA@data,nCores = 10) 
save(cells_rankings,file = "cells_rankings.rda")

cells_AUC <- AUCell_calcAUC(HYPOXIA_AUC, cells_rankings, nCores = 10,
                            aucMaxRank=nrow(cells_rankings)*0.05)


geneSet <- "LLPS_AUC"
AUCell_auc <- as.numeric(getAUC(cells_AUC))
sce2$AUCell <- AUCell_auc
head(sce2@meta.data)
dim(sce2@meta.data)
colnames(sce2@meta.data)[19] <- 'AUCell'
save(sce2,file = "06_PBMC.combined_AUCell.rda")
pdf("AUCell_viovin_all.pdf",width = 12,height = 6)
VlnPlot(sce2,features = 'AUCell', 
        pt.size = 0,group.by = "celltype")
dev.off() 
library(ggpubr)
pdf("AUCell_boxplot_all.pdf",width = 13,height = 5)
colour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#FFFF00",  
         "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
         "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
         "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
ggboxplot(sce2@meta.data, x="celltype", y="AUCell", width = 0.6, 
          color = "black",#
          fill="celltype",#
          palette = colour,
          xlab = F, #
          bxp.errorbar=T,
          bxp.errorbar.width=0.5, 
          size=0.5, 
          legend = "right") #
dev.off() 

library(igraph)
library(CellChat)
library(dplyr)
load('05_PBMC.combined.Rdata')
seurat.obj_PD <- PBMC.combined
dim(seurat.obj_PD)
data.input = seurat.obj_PD@assays$RNA@data
meta.data =  seurat.obj_PD@meta.data
meta.data = meta.data[!is.na(meta.data$celltype),]
data.input = data.input[,row.names(meta.data)]

table(meta.data$celltype)
meta.data$celltype = factor(meta.data$celltype,
                            levels = c("Macrophages", "Granulocyte", "VCT", "EVT","T cell","Myelocyte",
                                       "Monocytes" ,"NK cell","SCT","Fibroblasts","Mast cell","B cell"))
cellchat <- createCellChat(object = data.input, 
                           meta = meta.data, 
                           group.by = "celltype")
cellchat <- addMeta(cellchat, meta = meta.data)

levels(cellchat@idents) # show factor levels of the cell labels

groupSize <- as.numeric(table(cellchat@idents)) # 
CellChatDB <- CellChatDB.human #
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)

CellChatDB.use <- CellChatDB # simply use the default CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat,population.size = F)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
write.csv(df.net,"cellchat_result.csv")
head(df.net)

levels(cellchat@idents)


cellchat <- computeCommunProbPathway(cellchat)
head(cellchat@net)
head(cellchat@netP)
cellchat <- aggregateNet(cellchat)
save(cellchat,file = "cellchat_result.rda")

groupSize <- as.numeric(table(cellchat@idents))
pdf("cellchat.pdf",width = 10, height = 8)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize,
                 weight.scale = T, label.edge= F,
                 title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize,
                 weight.scale = T, label.edge= F,
                 title.name = "Interaction weights/strength")
dev.off()
mat <- cellchat@net$weight
pdf("cellchat_singlecell_weight.pdf",width = 14, height = 12)
par(mfrow=c(3,4),xpd=T)
for (i in 1:nrow(mat)) {
  mat1 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat1[i, ] <- mat[i, ]
  netVisual_circle(mat1, vertex.weight = groupSize, weight.scale = T, arrow.width = 0.2,arrow.size = 0.1,
                   edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()


mat <- cellchat@net$count
pdf("cellchat_singlecell_count.pdf",width = 14, height = 12)
par(mfrow=c(3,4),xpd=T)
for (i in 1:nrow(mat)) {
  mat1 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat1[i, ] <- mat[i, ]
  netVisual_circle(mat1, vertex.weight = groupSize, weight.scale = T, arrow.width = 0.2,arrow.size = 0.1,
                   edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

levels(cellchat@idents)
p=netVisual_bubble(cellchat,sources.use = c(3),targets.use = c(1,2,4,5,6,7,8,9,10,11,12),remove.isolate = F)
ggsave("cellchat_source_target.pdf",p,height = 8,width = 6)
p=netVisual_bubble(cellchat,sources.use = c(1,2,4,5,6,7,8,9,10,11,12),targets.use = c(3),remove.isolate = F)
ggsave("cellchat_target_source.pdf",p,height = 10,width = 6)

p=netVisual_bubble(cellchat,sources.use = c(3),targets.use = c(1,2,4,5,6,7,8,9,10),remove.isolate = F)
ggsave("cellchat_source_target1.pdf",p,height = 12,width = 8)
p=netVisual_bubble(cellchat,sources.use = c(1,2,4,5,6,7,8,9,10),targets.use = c(3),remove.isolate = F)
ggsave("cellchat_target_source1.pdf",p,height = 15,width = 8)


p=netVisual_bubble(cellchat,sources.use = c(1),targets.use = c(2,3,4,5,6,7,8,9,10),remove.isolate = F)
ggsave("cellchat_source_target2.pdf",p,height = 15,width = 8)
p=netVisual_bubble(cellchat,sources.use = c(2,3,4,5,6,7,8,9,10),targets.use = c(1),remove.isolate = F)
ggsave("cellchat_target_source2.pdf",p,height = 15,width = 8)

library(igraph)
packageVersion("igraph")
library(monocle)
load('05_PBMC.combined.Rdata')
DefaultAssay(PBMC.combined) <- "RNA"
table(PBMC.combined$celltype)
data1 = PBMC.combined[,which(PBMC.combined@active.ident == 'VCT')]
table(data1$Group)
seurat.obj_monocle<-subset(x =data1)
dim(seurat.obj_monocle)
seurat.obj_monocle1 = seurat.obj_monocle[rowSums(seurat.obj_monocle)>0,]
dim(seurat.obj_monocle1)

data <- as(as.matrix(seurat.obj_monocle1@assays$RNA@counts), 'sparseMatrix')

pd <- new('AnnotatedDataFrame', data = seurat.obj_monocle@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

cds <- newCellDataSet(data,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
diff=read.csv("VCT-差异分析.csv",header = T,stringsAsFactors = F,check.names = F)

save(diff,file = "monocle.diff.gene.rda")


load("monocle.diff.gene.rda")
deg <- subset(diff, p_val< 0.05) 
deg <- deg[order(deg$p_val,decreasing=F),]

ordergene <- unique(deg$gene)
cds <- setOrderingFilter(cds, ordergene)  
save(cds,file = "monocle.cds.rda")

pdf("monocle__DSP.pdf")
plot_ordering_genes(cds)
dev.off()

cds <- reduceDimension(cds, max_components = 2,method = "DDRTree")

dim(cds)
cds <- orderCells(cds)   
colnames(pData(cds))
save(cds,file = "pseudotime.rda")
library(ggsci)

pdf("monocle.pseudotime.pdf",width = 7,height = 7,onefile=FALSE)
plot_cell_trajectory(cds,color_by="Pseudotime", size=1,show_backbone=TRUE)
dev.off()

pdf("monocle.group.pdf",width = 7,height = 7)
plot_cell_trajectory(cds,color_by="Group", size=1,show_backbone=TRUE)
dev.off()

pdf("monocle.celltype.pdf",width = 7,height = 7)
plot_cell_trajectory(cds,color_by="celltype", size=1,show_backbone=TRUE)
dev.off()

pdf("monocle.State.pdf",width = 7,height = 7)
plot_cell_trajectory(cds, color_by = "State")
dev.off()

pdf("monocle.sample.pdf",width = 8,height = 8)
plot_cell_trajectory(cds, color_by = "State") + facet_wrap(~orig.ident,nrow = 4)
dev.off()
pdf("monocle.sample.group.pdf",width = 6,height = 7)
plot_cell_trajectory(cds, color_by = "State") + facet_wrap(~Group,nrow = 4)
dev.off()

blast_genes <- row.names(subset(fData(cds), gene_short_name %in% 
                                  c('SRSF10',
'RPL11')))
pdf("monocle.pseudotime.gene.pdf",width = 7,height = 7)
plot_genes_jitter(cds[blast_genes, ], grouping = "State", min_expr = 0.1)
dev.off()


ordering_genes= c('SRSF10',
                  'RPL11')

pdf("monocle.pseudotime.gene.pheatmap.pdf",width = 7,height = 7)
plot_pseudotime_heatmap(cds[ordering_genes, ],num_clusters = 5, show_rownames = T)
dev.off()

my_pseudotime_gene=ordering_genes
pdf("monocle.pseudotime.gene.by_state.pdf",width = 7,height = 7)
plot_genes_in_pseudotime(cds[my_pseudotime_gene,])+ scale_color_npg()
dev.off()
