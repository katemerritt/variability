#compression  = 'study',
run_variance<- function(mydata_selection, compression){
 # have variables in this order:  <- %>% <- <- <- <- <- <- <- <- <- <- <- <- %>% <- <- <- <- <- )
  #whole group SD mean correlation
  sdmean_corr <-  cor.test(log(c(mydata_selection[[1]], mydata_selection[[2]])), log(c(mydata_selection[[3]], mydata_selection[[4]])))$estimate
  
  # Get a summary estimate for each study (rather than each drug/dose) to make forest less cluttered
  VRestimates <-  VRvariances <- CVRestimates <- CVRvariances <-c()
  i <- 1
  if(compression =='study'){
      dataframe_list <- vector('list', length(unique(mydata_selection[[7]])))
      for (study in unique(mydata_selection[[7]])){
          dataframe_list[[i]] <- mydata_selection %>% filter(studyid==study)
          i=i+1
      }}
  i <- 1
  if(compression =='drug'){
      dataframe_list <- vector('list', length(unique(mydata_selection[[8]])))
      for (drugname in unique(mydata_selection[[8]])){
          dataframe_list[[i]] <- mydata_selection %>% filter(drug==drugname)
          i=i+1
      }}
  
  i<-1
  for(i in 1:length(dataframe_list)){
      dataframe <- dataframe_list[[i]]
    VR_tot <- Calc.lVR(dataframe[[3]], dataframe[[5]], dataframe[[4]],  dataframe[[6]])  # note functions take control parameters first
    var_VR_tot = Calc.var.lVR(dataframe[[5]],  dataframe[[6]])
    tot_VR_meta = rma.uni(yi=VR_tot,vi=var_VR_tot,method='REML',control=list(maxiter=1000,stepadj=0.5))
    VRestimates[i] <- as.numeric(coef(tot_VR_meta))
    VRvariances[i] <- as.numeric(vcov(tot_VR_meta))
    
    
    CVR_tot <- Calc.lnCVR(dataframe[[1]], dataframe[[3]],  dataframe[[5]], dataframe[[2]], dataframe[[4]], dataframe[[6]])  # note functions take control parameters first
    # need to use the correlation from whole group (which Calc.var.lnCVR normally calculates automatically)
    var_CVR_tot = Calc.var.lnCVRrob(dataframe[[1]], dataframe[[3]],  dataframe[[5]], dataframe[[2]], dataframe[[4]], dataframe[[6]], sdmean_corr)
    tot_CVR_meta = rma.uni(yi=CVR_tot,vi=var_CVR_tot,method='REML',control=list(maxiter=1000,stepadj=0.5))
    CVRestimates[i] <- as.numeric(coef(tot_CVR_meta))
    CVRvariances[i] <- as.numeric(vcov(tot_CVR_meta))
    i <- i+1
  }
  
  # VR
  VR_tot <- Calc.lVR(mydata_selection[[3]], mydata_selection[[5]], mydata_selection[[4]],  mydata_selection[[6]])  # note functions take control parameters first
  var_VR_tot = Calc.var.lVR(mydata_selection[[5]], mydata_selection[[6]])
  tot_VR_meta = rma.uni(yi=VR_tot,vi=var_VR_tot,method='REML',control=list(maxiter=1000,stepadj=0.5))
  
  # CVR
  CVR_tot <- Calc.lnCVR(mydata_selection[[1]], mydata_selection[[3]],  mydata_selection[[5]], mydata_selection[[2]], mydata_selection[[4]], mydata_selection[6])  # note functions take control parameters first
  var_CVR_tot = Calc.var.lnCVR(mydata_selection[[1]], mydata_selection[[3]],  mydata_selection[[5]], mydata_selection[[2]], mydata_selection[[4]], mydata_selection[6]) 
  tot_CVR_meta = rma.uni(yi=CVR_tot[,1],vi=var_CVR_tot[,1],method='REML',control=list(maxiter=1000,stepadj=0.5))
  
  #return(VRestimates, VRvariances, VRmeta, CVRestimates, CVRvariances, CVR_meta)
  return(list(VRestimates, VRvariances, CVRestimates, CVRvariances, tot_VR_meta, tot_CVR_meta))
}

run_smd<- function(mydata_selection, compression){
  # have variables in this order: 'pla_mn_tot_ch_adj', 'ap_mn_tot_ch_adj', 'pla_sd_tot_ch','ap_sd_tot_ch', 'pla_n_adjus', 'ap_n', 'studyid', 'drug')
  #whole group SD mean correlation
  # Get a summary estimate for each study (rather than each drug/dose) to make forest less cluttered
  SMDestimates <-  SMDvariances <-c()
  i <- 1
  if(compression =='study'){
    dataframe_list <- vector('list', length(unique(mydata_selection[[7]])))
    for (study in unique(mydata_selection[[7]])){
      dataframe_list[[i]] <- mydata_selection %>% filter(studyid==study)
      i=i+1
    }}
  i <- 1
  if(compression =='drug'){
    dataframe_list <- vector('list', length(unique(mydata_selection[[8]])))
    for (drugname in unique(mydata_selection[[8]])){
      dataframe_list[[i]] <- mydata_selection %>% filter(drug==drugname)
      i=i+1
    }}
  
  i<-1
  for(i in 1:length(dataframe_list)){
    dataframe <- dataframe_list[[i]]
    dataframe <- escalc(measure="SMD", m1i=dataframe[[1]], sd1i=dataframe[[3]], n1i=dataframe[[5]], m2i=dataframe[[2]], sd2i=dataframe[[4]], n2i=dataframe[[6]])
    tot_smd_meta = rma.uni(yi=yi,vi=vi,method='REML',control=list(maxiter=1000,stepadj=0.5), data=dataframe)
    SMDestimates[i] <- as.numeric(coef(tot_smd_meta))
    SMDvariances[i] <- as.numeric(vcov(tot_smd_meta))
    i <- i+1
  }
  mydata_selection_es <- escalc(measure="SMD", m1i=mydata_selection[[1]], sd1i=mydata_selection[[3]], n1i=mydata_selection[[5]], m2i=mydata_selection[[2]], sd2i=mydata_selection[[4]], n2i=mydata_selection[[6]])
  tot_SMD_meta = rma.uni(yi=yi,vi=vi,method='REML',control=list(maxiter=1000,stepadj=0.5), data=mydata_selection_es)
  
  return(list(SMDestimates, SMDvariances, SMDestimates, SMDvariances, tot_SMD_meta, tot_SMD_meta))
}



myforest_smd<- function(total_results, mydata_tot, compression){
  
  label = 'SMD'
  sorting_order <- order(total_results[[1]])
  yi <-total_results[[1]][sorting_order] ;vi <-total_results[[2]][sorting_order] ;rma <- total_results[[5]]
  
  #Calculate N's'
  df <- data.frame(study=character(), apn=integer(), plan=integer(), totn=integer(),stringsAsFactors = FALSE)
  for (study in unique(mydata_tot[[compression]])){
    #look for the study
    if (compression=='studyid'){
      df2<- mydata_tot %>% filter(studyid==study)}
    if (compression=='drug'){
      df2<- mydata_tot %>% filter(drug==study)}
    apn=sum(df2$ap_n)
    plan=as.integer(sum(df2$pla_n_adjus))
    totn = apn+plan
    df[nrow(df)+1,] = list(study, apn, plan, totn)
  }
  
  forest(yi, vi, slab = (unique(mydata_tot[[compression]]))[sorting_order], xlab= label, alim=c(-0.25, 1.25),
         ylim = c(0, length(total_results[[1]])+3),xlim=c(-1, 1.7),at=seq(-0.25,1.25, by=0.25), cex=0.75, refline=0,
         ilab = df$totn[sorting_order], ilab.xpos = -0.4, lty=c('solid', 'blank'))
  text(c(-0.8,-0.4, 1.5), length(total_results[[1]])+2, c("Study","N", "SMD (95% CI)"))
  text(c(0.0,1.0), length(total_results[[1]])+3, c("Greater efficacy \n in placebo", "Greater efficacy \n in antipsychotic"))
  if(compression=='drug'){
    addpoly(rma, row = 0, cex=1,  mlab="")}
  if(compression=='studyid'){
    addpoly(rma, row = -1, cex=1,  mlab="")}
  segments(-0.95,length(total_results[[1]])+1,-0.3,length(total_results[[1]])+1)
  segments(1.35,length(total_results[[1]])+1,1.75,length(total_results[[1]])+1)
}


#function for drawing forest plot
myforest<- function(total_results, mydata_tot, compression, vcvr){

    label = vcvr#paste("        Greater  in placebo               ", vcvr, "              Greater in antipsychotic")
    #compression = 'drug or 'studyid', vcvr = 'VR' or 'CVR'
    if (vcvr == 'VR'){
        sorting_order <- order(total_results[[1]])
        yi <-total_results[[1]][sorting_order] ;vi <-total_results[[2]][sorting_order] ;rma <- total_results[[5]]}
    if (vcvr == 'CVR'){
      sorting_order <- order(total_results[[3]])
        yi <-total_results[[3]][sorting_order] ;vi <-total_results[[4]][sorting_order] ;rma <- total_results[[6]]}

    #Calculate N's'
    df <- data.frame(study=character(), apn=integer(), plan=integer(), totn=integer(),stringsAsFactors = FALSE)
    for (study in unique(mydata_tot[[compression]])){
      #look for the study
      if (compression=='studyid'){
        df2<- mydata_tot %>% filter(studyid==study)}
      if (compression=='drug'){
        df2<- mydata_tot %>% filter(drug==study)}
      apn=sum(df2$ap_n)
      plan=as.integer(sum(df2$pla_n_adjus))
      totn = apn+plan
      df[nrow(df)+1,] = list(study, apn, plan, totn)
    }

    df3 <- data.frame(study=character(), apn=integer(), plan=integer(), totn=integer(),stringsAsFactors = FALSE)
    for (study in unique(mydata_tot[[compression]])){
      #look for the study
      if (compression=='studyid'){
        df4<- mydata_tot %>% filter(studyid==study)}
      if (compression=='drug'){
        df4<- mydata_tot %>% filter(drug==study)}
      apn=sum(df4$ap_n)
      plan=as.integer(sum(df4$pla_n_adjus))
      totn = apn+plan
      df3[nrow(df)+1,] = list(study, apn, plan, totn)
    }
    
    forest(yi, vi, slab = (unique(mydata_tot[[compression]]))[sorting_order], xlab= label, alim=c(0.5, 1.25),
           ylim = c(0, length(total_results[[1]])+3),xlim=c(0.15, 1.5),at=seq(0.5,1.25, by=0.25), transf=exp, cex=0.75, refline=1,
           ilab = df$totn[sorting_order], ilab.xpos = 0.43, lty=c('solid', 'blank'))
    text(c(0.25,0.43, 1.425), length(total_results[[1]])+2, c("Study","N", "CVR (95% CI)"))
    text(c(0.85,1.15), length(total_results[[1]])+3, c("Greater variability \n in placebo", "Greater variability \n in antipsychotic"))
    if(compression=='drug'){
    addpoly(rma, row = 0, transf=exp, cex=1,  mlab="")}
    if(compression=='studyid'){
      addpoly(rma, row = -1, transf=exp, cex=1,  mlab="")}
    segments(0.15,length(total_results[[1]])+1,0.45,length(total_results[[1]])+1)
    segments(1.35,length(total_results[[1]])+1,1.5,length(total_results[[1]])+1)
    }


myforest2<- function(total_results, mydata_tot, compression, vcvr){
  
  label = vcvr#paste("        Greater  in placebo               ", vcvr, "              Greater in antipsychotic")
  #compression = 'drug or 'studyid', vcvr = 'VR' or 'CVR'
  if (vcvr == 'VR'){
    sorting_order <- order(total_results[[1]])
    yi <-total_results[[1]][sorting_order] ;vi <-total_results[[2]][sorting_order] ;rma <- total_results[[5]]}
  if (vcvr == 'CVR'){
    sorting_order <- order(total_results[[3]])
    yi <-total_results[[3]][sorting_order] ;vi <-total_results[[4]][sorting_order] ;rma <- total_results[[6]]}
  
  #Calculate N's'
  df <- data.frame(study=character(), apn=integer(), plan=integer(), totn=integer(),stringsAsFactors = FALSE)
  for (study in unique(mydata_tot[[compression]])){
    #look for the study
    if (compression=='studyid'){
      df2<- mydata_tot %>% filter(studyid==study)}
    if (compression=='drug'){
      df2<- mydata_tot %>% filter(drug==study)}
    apn=sum(df2$ap_n)
    plan=as.integer(sum(df2$pla_n_adjus))
    totn = apn+plan
    df[nrow(df)+1,] = list(study, apn, plan, totn)
  }
  
  for (study in unique(mydata_tot[[compression]])){
    #look for the study
    if (compression=='studyid'){
      df4<- mydata_tot %>% filter(studyid==study)}
    if (compression=='drug'){
      df4<- mydata_tot %>% filter(drug==study)}
    arm_n[i]=length(unique(df4$studyid))
    i=i+1
  }
  
  forest(yi, vi, slab = (unique(mydata_tot[[compression]]))[sorting_order], xlab= label, alim=c(0.5, 1.25),
         ylim = c(0, length(total_results[[1]])+3),xlim=c(0.15, 1.5),at=seq(0.5,1.25, by=0.25), transf=exp, cex=0.75, refline=1,
         ilab = cbind(df$totn[sorting_order], arm_n[sorting_order]), ilab.xpos = c(0.37, 0.46),  lty=c('solid', 'blank'))
  text(c(0.21,0.37,0.44, 1.405), length(total_results[[1]])+2, c("Study","N", "Trials", "CVR (95% CI)"))
  text(c(0.85,1.15), length(total_results[[1]])+3, c("Greater variability \n in placebo", "Greater variability \n in antipsychotic"))
  if(compression=='drug'){
    addpoly(rma, row = 0, transf=exp, cex=1,  mlab="")}
  if(compression=='studyid'){
    addpoly(rma, row = -1, transf=exp, cex=1,  mlab="")}
  segments(0.15,length(total_results[[1]])+1,0.45,length(total_results[[1]])+1)
  segments(1.35,length(total_results[[1]])+1,1.5,length(total_results[[1]])+1)
}

#function for metaregression
metaregress<- function(year_results, mydata_year){
    res3 <- rma(year_results[[6]]$yi, year_results[[6]]$vi, data=mydata_year, method="REML", mods=year)
    wi <- 0.5/sqrt(year_results[[6]]$vi)
    preds <- predict(res3, addx=TRUE)
    plot(mydata_year$year, year_results[[6]]$yi, cex=wi, xlab="year", ylab="Effect Size", main = "Meta-regression - Publication year VR")
    lines(mydata_year$year, preds$pred)
    
    graph_data <- data.frame(exp(year_results[[6]]$yi), mydata_year$year,wi)
    year_graph <- graph_data %>% ggplot(aes_string(x=names(graph_data)[2], y=names(graph_data)[1], size='wi'))+
      geom_point(color='dodgerblue4', alpha=0.4)+geom_smooth(method='lm', color='dodgerblue4', fill='dodgerblue4')+
      xlab('Publication Year')+ylab('CVR')+theme_minimal()+theme(legend.position="none")
    return(list(res3,year_graph))
    }

#function for SMD-variance comparison
smd_regress <- function(mydata_tot){
  results_tot <- run_variance(mydata_tot, 'study')
  res3 <- rma(results_tot[[6]]$yi, results_tot[[6]]$vi,  method="REML", mods=as.vector(mydata_tot$yi))
  wi <- 0.5/sqrt(results_tot[[6]]$vi)
  graph_data <- data.frame(exp(results_tot[[6]]$yi),mydata_tot$yi,  wi)
  smd_graph <- graph_data %>% ggplot(aes_string(x=names(graph_data)[2], y=names(graph_data)[1], size='wi'))+
    geom_point(color='dodgerblue4', alpha=0.4)+geom_smooth(method='lm', color='dodgerblue4', fill='dodgerblue4')+
    xlab('SMD (Hedges g)')+ylab('CVR')+theme_minimal()+theme(legend.position="none")
  return(list(res3,smd_graph))
  }

#function for SMD-ch_Adj comparison
smdchadj_regress <- function(mydata_tot, pos_tot_neg){
  results_tot <- run_variance(mydata_tot, 'study')
  res3 <- rma(mydata_tot$yi, mydata_tot$vi,  method="REML", mods=as.vector(mydata_tot[pos_tot_neg]))
  wi <- 0.5/sqrt(mydata_tot$vi)
  graph_data <- data.frame(mydata_tot$yi, mydata_tot[pos_tot_neg], wi)
  smd_graph <- graph_data %>% ggplot(aes_string(x=names(graph_data)[2], y=names(graph_data)[1], size='wi'))+
    geom_point(color='dodgerblue4', alpha=0.4)+geom_smooth(method='lm', color='dodgerblue4', fill='dodgerblue4')+
    xlab('Symptom Change (Adjusted)')+ylab('SMD')+theme_minimal()+theme(legend.position="none")
  return(list(res3,smd_graph))
}

#function for CVR-ch_Adj comparison
cvrchadj_regress <- function(mydata_tot, pos_tot_neg){
  results_tot <- run_variance(mydata_tot, 'study')
  res3 <- rma(results_tot[[5]]$yi, results_tot[[5]]$vi,  method="REML", mods=as.vector(mydata_tot[pos_tot_neg]))
  wi <- 0.5/sqrt(results_tot[[5]]$vi)
  graph_data <- data.frame(mydata_tot[pos_tot_neg], results_tot[[5]]$yi, wi)
  smd_graph <- graph_data %>% ggplot(aes_string(x=names(graph_data)[2], y=names(graph_data)[1], size='wi'))+
    geom_point(color='dodgerblue4', alpha=0.4)+geom_smooth(method='lm', color='dodgerblue4', fill='dodgerblue4')+
    xlab('CVR')+ylab('Symptom Change (Adjusted)')+theme_minimal()+theme(legend.position="none")
  return(list(res3,smd_graph))
}

#function for CVR-age comparison
age_regress <- function(mydata_tot){
  results_tot <- run_variance(mydata_tot, 'study')
  res3 <- rma(results_tot[[6]]$yi, results_tot[[6]]$vi,  method="REML", mods=as.vector(mydata_tot$total_age))
  wi <- 0.5/sqrt(results_tot[[6]]$vi)
  graph_data <- data.frame(exp(results_tot[[6]]$yi),mydata_tot$total_age,  wi)
  smd_graph <- graph_data %>% ggplot(aes_string(x=names(graph_data)[2], y=names(graph_data)[1], size='wi'))+
    geom_point(color='dodgerblue4', alpha=0.4)+geom_smooth(method='lm', color='dodgerblue4', fill='dodgerblue4')+
    xlab('Age (Years)')+ylab('CVR')+theme_minimal()+theme(legend.position="none")
  return(list(res3,smd_graph))
}

#function for CVR-gender comparison
gender_regress <- function(mydata_tot){
  results_tot <- run_variance(mydata_tot, 'study')
  res3 <- rma(results_tot[[6]]$yi, results_tot[[6]]$vi,  method="REML", mods=as.vector(mydata_tot$total_male))
  wi <- 0.5/sqrt(results_tot[[6]]$vi)
  graph_data <- data.frame(exp(results_tot[[6]]$yi),mydata_tot$total_male,  wi)
  smd_graph <- graph_data %>% ggplot(aes_string(x=names(graph_data)[2], y=names(graph_data)[1], size='wi'))+
    geom_point(color='dodgerblue4', alpha=0.4)+geom_smooth(method='lm', color='dodgerblue4', fill='dodgerblue4')+
    xlab('% Male')+ylab('CVR')+theme_minimal()+theme(legend.position="none")
  return(list(res3,smd_graph))
}

#function for CVR-gender comparison
dose_regress <- function(mydata_tot){
  results_tot <- run_variance(mydata_tot, 'study')
  res3 <- rma(results_tot[[6]]$yi, results_tot[[6]]$vi,  method="REML", mods=as.vector(mydata_tot$olanz))
  wi <- 0.5/sqrt(results_tot[[6]]$vi)
  graph_data <- data.frame(exp(results_tot[[6]]$yi),mydata_tot$olanz,  wi)
  smd_graph <- graph_data %>% ggplot(aes_string(x=names(graph_data)[2], y=names(graph_data)[1], size='wi'))+
    geom_point(color='dodgerblue4', alpha=0.4)+geom_smooth(method='lm', color='dodgerblue4', fill='dodgerblue4')+
    xlab('Dose (Olanzapine Equivalents mg/day)')+ylab('CVR')+theme_minimal()+theme(legend.position="none")
  return(list(res3,smd_graph))
}

severity_regress <- function(mydata_tot){
  results_tot <- run_variance(mydata_tot, 'study')
  res3 <- rma(results_tot[[6]]$yi, results_tot[[6]]$vi,  method="REML", mods=as.vector(mydata_tot$combined_tot_bl))
  wi <- 0.5/sqrt(results_tot[[6]]$vi)
  graph_data <- data.frame(exp(results_tot[[6]]$yi),mydata_tot$combined_tot_bl,  wi)
  smd_graph <- graph_data %>% ggplot(aes_string(x=names(graph_data)[2], y=names(graph_data)[1], size='wi'))+
    geom_point(color='dodgerblue4', alpha=0.4)+geom_smooth(method='lm', color='dodgerblue4', fill='dodgerblue4')+
    xlab('Baseline PANSS')+ylab('CVR')+theme_minimal()+theme(legend.position="none")
  return(list(res3,smd_graph))
}

length_regress <- function(mydata_tot){
  results_tot <- run_variance(mydata_tot, 'study')
  res3 <- rma(results_tot[[6]]$yi, results_tot[[6]]$vi,  method="REML", mods=as.vector(mydata_tot$week_of_assessment))
  wi <- 0.5/sqrt(results_tot[[6]]$vi)
  graph_data <- data.frame(exp(results_tot[[6]]$yi),mydata_tot$week_of_assessment,  wi)
  smd_graph <- graph_data %>% ggplot(aes_string(x=names(graph_data)[2], y=names(graph_data)[1], size='wi'))+
    geom_point(color='dodgerblue4', alpha=0.4)+geom_smooth(method='lm', color='dodgerblue4', fill='dodgerblue4')+
    xlab('Week of Assessment')+ylab('CVR')+theme_minimal()+theme(legend.position="none")
  return(list(res3,smd_graph))
}
