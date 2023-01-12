data_clean<- function(mydata){
# Remove ineligible rows
mydata <- mydata %>% filter(!is.na(sd_se))

# Make numeric
mydata <- mydata %>% mutate_at(vars(ap_n:pla_sd_tot_ch), funs(as.numeric(as.character(.))))
mydata$year <- as.numeric(mydata$year)

# Convert SE to SD
mydata <- mydata %>%
  mutate_at(vars(ap_sd_pos_ch, ap_sd_neg_ch,ap_sd_tot_ch), funs(ifelse(sd_se=="se",(.*ap_n^0.5), .))) 
mydata <- mydata %>%
  mutate_at(vars(pla_sd_pos_ch, pla_sd_neg_ch, pla_sd_tot_ch), funs(ifelse(sd_se=="se",(.*pla_n_orig^0.5), .))) 

#Calculate age & gender
mydata<- mutate(mydata, total_age = ifelse(is.na(total_age), ((ap_age*ap_n+pla_age*pla_n_orig)/(ap_n+pla_n_orig)), total_age))
mydata<- mutate(mydata, total_male = ifelse(is.na(total_male), ((ap_male*ap_n+pla_male*pla_n_orig)/(ap_n+pla_n_orig)), total_male))

#calculate olanzapine equivalents
mydata$olanz <- NA
drug_list <- c("amisulpride","aripiprazole","asenapine","blonanserin","chlorpromazine","haloperidol","lurasidone","olanzapine","paliperidone","quetiapine","risperidone","sertindole","ziprasidone","zotepine")
#"brexpiprazole","cariprazine","iloperidone",
olanz_list <- c(40, 1.5, 2, 1.6, 30,0.8, 6, 1, 0.6, 40, 0.5, 1.6, 8, 20)
i=1
for(drugname in drug_list){
  mydata<- mutate(mydata, olanz = ifelse((drug==drugname), (dose/olanz_list[i]), olanz))
  i=i+1
}

# Make zero based floor
#make sure we are not subtracting from the SANS studies
mydata<- mutate(mydata, ap_mn_tot_bl = ifelse(((non_zero_floor=='non zero floor' | is.na(non_zero_floor)) & total_scale=="PANSS"), (ap_mn_tot_bl-30), ap_mn_tot_bl))
mydata<- mutate(mydata, pla_mn_tot_bl = ifelse(((non_zero_floor=='non zero floor' | is.na(non_zero_floor)) & total_scale=="PANSS"), (pla_mn_tot_bl-30), pla_mn_tot_bl))
mydata<- mutate(mydata, ap_mn_tot_bl = ifelse(((non_zero_floor=='non zero floor' | is.na(non_zero_floor)) & total_scale=="bprs"), (ap_mn_tot_bl-24), ap_mn_tot_bl))
mydata<- mutate(mydata, pla_mn_tot_bl = ifelse(((non_zero_floor=='non zero floor' | is.na(non_zero_floor)) & total_scale=="bprs"), (pla_mn_tot_bl-24), pla_mn_tot_bl))
mydata<- mutate(mydata, ap_mn_tot_bl = ifelse(((non_zero_floor=='non zero floor' | is.na(non_zero_floor)) & total_scale=="bprs18"), (ap_mn_tot_bl-18), ap_mn_tot_bl))
mydata<- mutate(mydata, pla_mn_tot_bl = ifelse(((non_zero_floor=='non zero floor' | is.na(non_zero_floor)) & total_scale=="bprs18"), (pla_mn_tot_bl-18), pla_mn_tot_bl))
mydata<- mutate(mydata, ap_mn_tot_bl = ifelse(((non_zero_floor=='non zero floor' | is.na(non_zero_floor)) & total_scale=="bprs16"), (ap_mn_tot_bl-16), ap_mn_tot_bl))
mydata<- mutate(mydata, pla_mn_tot_bl = ifelse(((non_zero_floor=='non zero floor' | is.na(non_zero_floor)) & total_scale=="bprs16"), (pla_mn_tot_bl-16), pla_mn_tot_bl))

mydata<- mutate(mydata, ap_mn_pos_bl = ifelse(((non_zero_floor=='non zero floor')&(positive_scale=="PANSS_p")), (ap_mn_pos_bl-7), ap_mn_pos_bl))
mydata<- mutate(mydata, ap_mn_pos_bl = ifelse(((non_zero_floor=='non zero floor')&(positive_scale=="bprs4")), (ap_mn_pos_bl-4), ap_mn_pos_bl))
mydata<- mutate(mydata, ap_mn_pos_bl = ifelse(((non_zero_floor=='non zero floor')&(positive_scale=="bprs")), (ap_mn_pos_bl-4), ap_mn_pos_bl))
mydata<- mutate(mydata, ap_mn_pos_bl = ifelse(((non_zero_floor=='non zero floor')&(positive_scale=="PANSS_marder")), (ap_mn_pos_bl-8), ap_mn_pos_bl))
mydata<- mutate(mydata, pla_mn_pos_bl = ifelse(((non_zero_floor=='non zero floor')&(positive_scale=="PANSS_p")), (pla_mn_pos_bl-7), pla_mn_pos_bl))
mydata<- mutate(mydata, pla_mn_pos_bl = ifelse(((non_zero_floor=='non zero floor')&(positive_scale=="bprs4")), (pla_mn_pos_bl-4), pla_mn_pos_bl))
mydata<- mutate(mydata, pla_mn_pos_bl = ifelse(((non_zero_floor=='non zero floor')&(positive_scale=="bprs")), (pla_mn_pos_bl-4), pla_mn_pos_bl))
mydata<- mutate(mydata, pla_mn_pos_bl = ifelse(((non_zero_floor=='non zero floor')&(positive_scale=="PANSS_marder")), (pla_mn_pos_bl-8), pla_mn_pos_bl))

mydata<- mutate(mydata, ap_mn_neg_bl = ifelse(((non_zero_floor=='non zero floor')&(negative_scale=="PANSS_n")), (ap_mn_neg_bl-7), ap_mn_neg_bl))
mydata<- mutate(mydata, ap_mn_neg_bl = ifelse(((non_zero_floor=='non zero floor')&(negative_scale=="PANSS_marder")), (ap_mn_neg_bl-7), ap_mn_neg_bl))
mydata<- mutate(mydata, ap_mn_neg_bl = ifelse(((non_zero_floor=='non zero floor')&(negative_scale=="bprs3")), (ap_mn_neg_bl-3), ap_mn_neg_bl))
mydata<- mutate(mydata, pla_mn_neg_bl = ifelse(((non_zero_floor=='non zero floor')&(negative_scale=="PANSS_n")), (pla_mn_neg_bl-7), pla_mn_neg_bl))
mydata<- mutate(mydata, pla_mn_neg_bl = ifelse(((non_zero_floor=='non zero floor')&(negative_scale=="PANSS_marder")), (pla_mn_neg_bl-7), pla_mn_neg_bl))
mydata<- mutate(mydata, pla_mn_neg_bl = ifelse(((non_zero_floor=='non zero floor')&(negative_scale=="bprs3")), (pla_mn_neg_bl-3), pla_mn_neg_bl))


#Impute missing baselines - total_symptoms
#For all studies if they dont have baseline they dont have FU either: summary(is.na(mydata$ap_mn_pos_bl)  & !is.na(mydata$ap_mn_pos_fu))
#Although some have change scores but no baseline: summary(is.na(mydata$ap_mn_tot_bl)  & !is.na(mydata$ap_mn_tot_ch))

#For those without total baseline take mean of the entire sample
bprs18_ap_mean <- mydata %>% filter(total_scale=="bprs18") %>% .$ap_mn_tot_bl %>% mean(na.rm=TRUE)
bprs18_pla_mean <- mydata %>% filter(total_scale=="bprs18") %>% .$pla_mn_tot_bl %>% mean(na.rm=TRUE)
bprs18_mean <- (bprs18_ap_mean+bprs18_pla_mean)/2

PANSS_ap_mean <- mydata %>% filter(total_scale=="PANSS") %>% .$ap_mn_tot_bl %>% mean(na.rm=TRUE)
PANSS_pla_mean <- mydata %>% filter(total_scale=="PANSS") %>% .$pla_mn_tot_bl %>% mean(na.rm=TRUE)
PANSS_mean <- (PANSS_ap_mean+PANSS_pla_mean)/2

#Replace empty cells with mean
mydata$missing_baseline=0
mydata$missing_baseline[is.na(mydata$ap_mn_tot_bl)|is.na(mydata$pla_mn_tot_bl)] <- 1
mydata$ap_mn_tot_bl[is.na(mydata$ap_mn_tot_bl)&mydata$total_scale=="PANSS"] <- PANSS_mean
mydata$ap_mn_tot_bl[is.na(mydata$ap_mn_tot_bl)&mydata$total_scale=="bprs18"] <- bprs18_mean
mydata$pla_mn_tot_bl[is.na(mydata$pla_mn_tot_bl)&mydata$total_scale=="PANSS"] <- PANSS_mean
mydata$pla_mn_tot_bl[is.na(mydata$pla_mn_tot_bl)&mydata$total_scale=="bprs18"] <- bprs18_mean

mydata$combined_tot_bl <- (mydata$ap_n*mydata$ap_mn_tot_bl + mydata$pla_n_orig*mydata$pla_mn_tot_bl)/(mydata$pla_n_orig+mydata$ap_n)

#Positive
PANSS_ap_pos_mean <- mydata %>% filter(positive_scale=="PANSS_p") %>% .$ap_mn_pos_bl %>% mean(na.rm=TRUE)
PANSS_pla_pos_mean <- mydata %>% filter(positive_scale=="PANSS_p") %>% .$pla_mn_pos_bl %>% mean(na.rm=TRUE)
PANSS_pos_mean <- (PANSS_ap_pos_mean+PANSS_pla_pos_mean)/2

#Replace empty cells with mean
mydata$ap_mn_pos_bl[is.na(mydata$ap_mn_pos_bl)&mydata$positive_scale=="PANSS_p"] <- PANSS_pos_mean
mydata$pla_mn_pos_bl[is.na(mydata$pla_mn_pos_bl)&mydata$positive_scale=="PANSS_p"] <- PANSS_pos_mean


#Negative
bprs_ap_neg_mean <- mydata %>% filter(total_scale=="bprs") %>% .$ap_mn_neg_bl %>% mean(na.rm=TRUE)
bprs_pla_neg_mean <- mydata %>% filter(total_scale=="bprs") %>% .$pla_mn_neg_bl %>% mean(na.rm=TRUE)
bprs_neg_mean <- (bprs_ap_neg_mean+bprs_pla_neg_mean)/2

PANSS_ap_neg_mean <- mydata %>% filter(total_scale=="PANSS_n") %>% .$ap_mn_neg_bl %>% mean(na.rm=TRUE)
PANSS_pla_neg_mean <- mydata %>% filter(total_scale=="PANSS_n") %>% .$pla_mn_neg_bl %>% mean(na.rm=TRUE)
PANSS_neg_mean <- (PANSS_ap_neg_mean+PANSS_pla_neg_mean)/2

#Replace empty cells with mean
mydata$ap_mn_neg_bl[is.na(mydata$ap_mn_neg_bl)&mydata$total_scale=="PANSS_n"] <- PANSS_neg_mean
mydata$ap_mn_neg_bl[is.na(mydata$ap_mn_neg_bl)&mydata$total_scale=="bprs"] <- bprs_neg_mean
mydata$pla_mn_neg_bl[is.na(mydata$pla_mn_neg_bl)&mydata$total_scale=="PANSS"] <- PANSS_neg_mean
mydata$pla_mn_neg_bl[is.na(mydata$pla_mn_neg_bl)&mydata$total_scale=="bprs_n"] <- bprs_neg_mean


#NB need to sort out SD situation for Corp, sliwa and umbricht
mydata$ap_mn_tot_ch_adj <- abs(mydata$ap_mn_tot_ch)+mydata$ap_mn_tot_bl
mydata$pla_mn_tot_ch_adj <- abs(mydata$pla_mn_tot_ch)+mydata$pla_mn_tot_bl
mydata$ap_mn_pos_ch_adj <- abs(mydata$ap_mn_pos_ch)+mydata$ap_mn_pos_bl
mydata$pla_mn_pos_ch_adj <- abs(mydata$pla_mn_pos_ch)+mydata$pla_mn_pos_bl
mydata$ap_mn_neg_ch_adj <- abs(mydata$ap_mn_neg_ch)+mydata$ap_mn_neg_bl
mydata$pla_mn_neg_ch_adj <- abs(mydata$pla_mn_neg_ch)+mydata$pla_mn_neg_bl

# Add study identifier
mydata$studyid <- paste(mydata$author, mydata$year)

return(mydata)
}