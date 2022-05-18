##########################
# allocation probabilities
##########################

alloca <- function (K,rule,Tmax, Nsim, success, failure){
  alp <- array(0,c(1,K)) 
  prior <- array(1,c(1,K));
  success <- success+prior 
  failure <- failure+ prior
  n_patients <- (success+failure) #patients in each arm
  t <- sum(n_patients)-2*K
  #print(n_patients)
  #print(t)
  be <- array(0,c(Nsim,K)) 
  for (j in 1:K) {
    be[,j] <- rbeta(Nsim,success[j],failure[j])
  }##################################
  
  #Thompson sampling
  if (rule == 'Thompson') {
    bmax <- apply(be,1,which.max)
    c <- ((t+1)/(2*Tmax))
    for (j in 1:K){
      alp[j] <- (sum(bmax == j)/Nsim)^c                    
    }
    
    
    #print(c)
    #print(t)
    constant <- sum(alp)
    alp <- alp/constant
    print(alp)
    print(constant)
    allocation_prob <- t(alp)
  }#################################### 
  # Trippa et al. Procedure
  else if (rule == 'Trippa') {
    factor <- array(0,c(1,K)) 
    elp <- array(0,c(1,K)) 
    # alp <- array(0,c(1,K)) 
    gamma <- 13*((t+1)/Tmax)^0.75
    nu <- ((t+1)/Tmax)*0.25
    for (j in 2:K) {
      factor[j] <- mean(be[,j]>be[,1])^gamma
    }
    total <- sum(factor)
    for (j in 2:K){
      elp[j] <- factor[j]/total
    }
    constant <- max(n_patients[2:K])
    elp[1] <- 1/K*exp(constant-n_patients[1])^nu
    for (j in 1:K){
      alp[j] <- elp[j]/sum(elp)
    }
    allocation_prob <- t(alp)
  }
}


##########################
# power calculations
##########################
                          
do_inference = function(ss = 20, delta = 0.3, Nsim_alloc = 500, replicas = 5000, b = 4, K = 3, alpha = 0.025, A){
  # details of the simulation studies
  # ss=20; #max sample size
  # delta=0.3; #30% increase
  # replica=5000;
  # b= 4; #block size
  # K=3; # arms = number of treatment arms in the trial (including the control)
  # Nsim=500; #simulations for RAR algorithm
  # A = real data of a biomarker to evaluate
  
  interims= ss/b
  Tmax =ss
  
  pf0 <-array(0,c(1,replicas))  
  pf1 <-array(0,c(1,replicas)) 
  pf2 <-array(0,c(1,replicas)) 
  mpower <- array(0, c(1, replicas))
  error <- array(0, c(1, replicas))

  success <- array(0,c(replicas,K)) 
  failure <- array(0,c(replicas,K))
  outcome <- array(0,c(replicas,ss))
  action <- array(0,c(replicas,ss))
  
  ### generate sample and trial replicas with our proposed RAR algorithm
  for (i in 1:replicas){
    nb=sample(A, ss, replace = TRUE);
    nfc=sample(A, ss, replace = TRUE);
    nft=sample((1+delta)*A, ss, replace = TRUE);

    ##############################################
    for (block in 1:interims){
      #get probabilities for block j. Start with equal randomisation. 
      allocation_prob <- alloca(K, 'Trippa',Tmax, Nsim_alloc, success[i,], failure[i,])
      coin<- runif(b)

      for (pts in 1:b){
        if (coin[pts]< allocation_prob[1]){
          action[i,(block-1)*b+pts] <- 1
          if ((nfc[(block-1)*b+pts]-nb[(block-1)*b+pts])/nb[(block-1)*b+pts] >= delta){
            outcome[i,(block-1)*b+pts] <- 1}
          success[i,1] <-  success[i,1]+outcome[i,(block-1)*b+pts]
          failure[i,1]<-  failure[i,1]+(1-outcome[i,(block-1)*b+pts])
        } else if (allocation_prob[1] <= coin[pts]  & coin[pts]< (allocation_prob[1]+allocation_prob[2])) {
          #outcome[i,(block-1)*b+pts] <- floor((nfc[(block-1)*b+pts]-nb[(block-1)*b+pts])/nb[(block-1)*b+pts]-.3)+1
          action[i,(block-1)*b+pts] <- 2
          if ((nfc[(block-1)*b+pts]-nb[(block-1)*b+pts])/nb[(block-1)*b+pts]>= delta){
            outcome[i,(block-1)*b+pts] <- 1}
          success[i,2] <-  success[i,2]+outcome[i,(block-1)*b+pts]
          failure[i,2]<-  failure[i,2]+ (1-outcome[i,(block-1)*b+pts])
        }
        else if ((allocation_prob[1]+allocation_prob[2]) <= coin[pts]){
          #outcome[i,(block-1)*b+pts] <- floor((nft[(block-1)*b+pts]-nb[(block-1)*b+pts])/nb[(block-1)*b+pts]-.3)+1
          action[i,(block-1)*b+pts] <- 3
          if ((nft[(block-1)*b+pts]-nb[(block-1)*b+pts])/nb[(block-1)*b+pts]>= delta){
            outcome[i,(block-1)*b+pts] <- 1}
          success[i,3] <-  success[i,3]+ outcome[i,(block-1)*b+pts]
          failure[i,3]<-  failure[i,3]+(1-outcome[i,(block-1)*b+pts])
        }
      }
    }
  
    pf0[i]=sum(action[i,]==1)/ss
    pf1[i]=sum(action[i,]==2)/ss
    pf2[i]=sum(action[i,]==3)/ss
    wilcox.test(nfc, jitter(nft, amount=0.001), paired = TRUE, alternative = "less")
    ## power
    mpower[i]=wilcox.test(nb, jitter(nft, amount=0.000001), paired = TRUE, alternative = "less")$p.value
    error[i]=wilcox.test(nb, jitter(nfc, amount=0.000001), paired = TRUE, alternative = "less")$p.value
    #test[i]=wilcox.test(nfc, nft, paired = TRUE, alternative = "two.sided", exact = FALSE)
  }
  
  return(list(alloc_arm0 = mean(pf0), alloc_arm1 = mean(pf1), alloc_arm2 = mean(pf2),
              power = mean(mpower<alpha), error = mean(error<alpha)))
}


##########################
# data
##########################

# current data: Stratosph part1 & COHORT
library(readxl)
library(dplyr)
#DB_2T <- read_excel("/Volumes/GoogleDrive/.shortcut-targets-by-id/1TnabQdwBZI4EkLOUMvh3kIQi0eV_8HG_/Sofia Villar - Nina Deliu/STRATOSPHERE/Prior Data/Stratosphere_qPCR/Stratosphere_qPCR_baseline 4mth Only.xlsx")

#RJ_expr <- read_excel("/Volumes/GoogleDrive/.shortcut-targets-by-id/1TnabQdwBZI4EkLOUMvh3kIQi0eV_8HG_/Sofia Villar - Nina Deliu/STRATOSPHERE/Prior Data/COHORT_RNA/RJ_expression data Allan-Mark.xlsx")



# All_biomarkers_qPCR_RU = list(ID3 = subset(DB_2T$ID3.RU, (DB_2T$Case=="Case")),
#                               SMAD1 = subset(DB_2T$SMAD1.RU, (DB_2T$Case=="Case")),
#                               SMAD5 = subset(DB_2T$SMAD5.RU, (DB_2T$Case=="Case")),
#                               NOTCH1 = subset(DB_2T$NOTCH1.RU, (DB_2T$Case=="Case")),
#                               NOTCH2 = subset(DB_2T$NOTCH2.RU, (DB_2T$Case=="Case")),
#                               ID2 = subset(DB_2T$ID2.RU, (DB_2T$Case=="Case")),
#                               ARL4C = subset(DB_2T$ARL4C.RU, (DB_2T$Case=="Case")),
#                               PTGS2 = subset(DB_2T$PTGS2.RU, (DB_2T$Case=="Case")))
# 
# All_biomarkers_qPCR_dCT = list(ID3 = subset(DB_2T$ID3.dCT, (DB_2T$Case=="Case")),
#                               SMAD1 = subset(DB_2T$SMAD1.dCT, (DB_2T$Case=="Case")),
#                               SMAD5 = subset(DB_2T$SMAD5.dCT, (DB_2T$Case=="Case")),
#                               NOTCH1 = subset(DB_2T$NOTCH1.dCT, (DB_2T$Case=="Case")),
#                               NOTCH2 = subset(DB_2T$NOTCH2.dCT, (DB_2T$Case=="Case")),
#                               ID2 = subset(DB_2T$ID2.dCT, (DB_2T$Case=="Case")),
#                               ARL4C = subset(DB_2T$ARL4C.dCT, (DB_2T$Case=="Case")),
#                               PTGS2 = subset(DB_2T$PTGS2.dCT, (DB_2T$Case=="Case")))
# 
# 
# All_biomarkers_RNA = list(ID3 = subset(RJ_expr$ID3, (RJ_expr$Group=="Case")),
#                           SMAD1 = subset(RJ_expr$SMAD1, (RJ_expr$Group=="Case")),
#                           SMAD5 = subset(RJ_expr$SMAD5, (RJ_expr$Group=="Case")),
#                           NOTCH1 = subset(RJ_expr$NOTCH1, (RJ_expr$Group=="Case")),
#                           NOTCH2 = subset(RJ_expr$NOTCH2, (RJ_expr$Group=="Case")),
#                           ARL4C = subset(RJ_expr$ARL4C, (RJ_expr$Group=="Case")),
#                           PTGS2 = subset(RJ_expr$PTGS2, (RJ_expr$Group=="Case")),
#                           ID2 = subset(RJ_expr$ID2, (RJ_expr$Group=="Case")),
#                           MYB = subset(RJ_expr$MYB, (RJ_expr$Group=="Case")),
#                           BMPR2 = subset(RJ_expr$BMPR2, (RJ_expr$Group=="Case")))


# qPCR_res_RU = matrix(NA, nrow = 5, ncol = length(All_biomarkers_qPCR_RU))
# for (i in 1:length(All_biomarkers_qPCR_RU)){
#   qPCR_res_RU[,i] = unlist(do_inference(ss = 20, delta = 0.3, Nsim_alloc = 500, replicas = 5000, b = 4, K = 3, alpha = 0.025, A = All_biomarkers_qPCR_RU[[i]]))
# }
# 
# round(qPCR_res_RU, 2)
# 
# 
# 
# qPCR_res_dCT = matrix(NA, nrow = 5, ncol = length(All_biomarkers_qPCR_dCT))
# for (i in 1:length(All_biomarkers_qPCR_dCT)){
#   qPCR_res_dCT[,i] = unlist(do_inference(ss = 20, delta = 0.3, Nsim_alloc = 500, replicas = 5000, b = 4, K = 3, alpha = 0.025, A = All_biomarkers_qPCR_dCT[[i]]))
# }
# 
# round(qPCR_res_dCT, 2)
# 
# colnames(qPCR_res_RU) = colnames(qPCR_res_dCT) =c("ID3", "SMAD1", "SMAD5", "NOTCH1", "NOTCH2", "ID2", "ARL4C", "PTGS2")
# 
# 
# RNA_res = matrix(NA, nrow = 5, ncol = length(All_biomarkers_RNA))
# for (i in 1:length(All_biomarkers_RNA)){
#   RNA_res[,i] = unlist(do_inference(ss = 20, delta = 0.3, Nsim_alloc = 500, replicas = 5000, b = 4, K = 3, alpha = 0.025, A = All_biomarkers_RNA[[i]]))
# }
# 
# colnames(RNA_res) =c("ID3", "SMAD1", "SMAD5", "NOTCH1", "NOTCH2", "ID2", "ARL4C", "PTGS2", "MYB", "BMPR2")
# 
# round(RNA_res, 2)


#A (real) data example is shared to document the power analysis

biom_A = c(3.651169, 4.681626, 3.942112, 4.715114, 5.549290, 6.255736, 4.651201, 4.491330, 5.660505, 5.075578, 9.009672, 
           7.855070, 6.037752, 5.571764, 7.875398, 5.948716, 5.911751, 6.204240, 6.813782, 5.647305, 5.238439, 8.421222,
           7.890520, 6.771511, 7.474281, 9.466536, 6.425828, 5.956931, 5.880286, 6.229845, 7.550006, 6.571509, 8.228874,
           6.697354, 5.567237, 7.345511, 5.052818, 4.654433, 5.335524, 6.844448)


do_inference(ss = 20, delta = 0.3, Nsim_alloc = 500, replicas = 5000, b = 4, K = 3, alpha = 0.025, A = biom_A)
