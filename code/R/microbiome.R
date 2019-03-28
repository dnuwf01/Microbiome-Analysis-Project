


## read the relative abundance file for taxa
## check correlation
## centered log ratio transform
## apply lasso and predict msm/msw 
msm_regressor = read.csv(file = "../../data/msm_tt_rel_abundance_r.csv",header = TRUE, sep = ",")
msw_regressor = read.csv(file = "../../data/msw_tt_rel_abundance_r.csv",header = TRUE, sep = ",")


corb_m = round(cor(msm_regressor, method = 'spearman'),2)
corb_w = round(cor(msw_regressor, method = 'spearman'),2)


## correlation plot for the bacteria
ggcorrplot(corb_m, outline.color = "white", ggtheme = ggplot2::theme_gray,
           , insig="blank")
ggcorrplot(corb_w, outline.color = "white", ggtheme = ggplot2::theme_gray,
           , insig="blank")

## centered log-ratio transform for regressors along each row
library(compositions)

for(i in 1:43){
  msm_regressor[i,] = as.list(clr(msm_regressor[i,]))
  msw_regressor[i,] = as.list(clr(msw_regressor[i,]))
}

## add a column for msm and msw status and then 
## add a status bar
msm_regressor$status = rep(1,43)
msw_regressor$status = rep(0,43)


combined_reg = rbind(msm_regressor,msw_regressor) 

## write it into a csv file
write.csv(combined_reg,file = "C:/Users/Debarghya Nandi/Desktop/research/mehta/deesha_dulal_mehta/bacteria_lt.csv")


## ----------------------------cytokines ----------------------------------
#read the csv file from the given location
out_msm = read.csv(file = "../Documents/GitHub/Microbiome-Analysis-Project/data/MSM_urinary_cytokines.csv",header = TRUE, sep = ",")
out_msw = read.csv(file = "../Documents/GitHub/Microbiome-Analysis-Project/data/MSW_urinary_cytokines.csv",header = TRUE, sep = ",")
cytokines = c('tnfa','il1b','il8','il10','ip10')


## select the important cytokines
cytokine_frame_m = out_msm[cytokines]
cytokine_frame_w= out_msw[cytokines]





#cormatt = rcorr(as.matrix(cytokine_frame),type = "spearman")

cor_m = round(cor(cytokine_frame_m,method = 'spearman'),2)
cor_w = round(cor(cytokine_frame_w, method = 'spearman'),2)


ggcorrplot(cor_m, hc.order = TRUE, outline.color = "white", ggtheme = ggplot2::theme_gray,
          , insig="blank")
ggcorrplot(cor_w, hc.order = TRUE, outline.color = "white", ggtheme = ggplot2::theme_gray,
           , insig="blank")




## carrying out a t-test to determine difference of meanstyep
df = rbind(data.frame(group = 'msm_tnfa',outcome = log_cytokines_m$tnfa),data.frame(group = 'msw_tnfa',outcome = log_cytokines_w$tnfa))
summary(df)


#read in the data for the linear regression


msm_regressor = read.csv(file = "../../data/msm_tt_rel_abundance_r.csv",header = TRUE, sep = ",")
msw_regressor = read.csv(file = "../../data/msw_tt_rel_abundance_r.csv",header = TRUE, sep = ",")


## centered log-ratio transform for regressors along each row
library(compositions)

for(i in 1:43){
  msm_regressor[i,] = as.list(clr(msm_regressor[i,]))
  msw_regressor[i,] = as.list(clr(msw_regressor[i,]))
}

## log transform the cytokine data
log_cytokine_frame_m = log(cytokine_frame_m)
log_cytokine_frame_w = log(cytokine_frame_w)



## add a status bar
msm_regressor$status = rep(1,43)
msw_regressor$status = rep(0,43)

### set index for msm
msm_regressor$index = c(1:43)
log_cytokine_frame_m$index = c(1:43)


## set index for msw
msw_regressor$index = c(1:43)
log_cytokine_frame_w$index = c(1:43)


## merge male reg and out
merge_msm = merge(msm_regressor,log_cytokine_frame_m,by="index")
merge_msw = merge(msw_regressor,log_cytokine_frame_w,by="index")


## bind the msm/msw regressors and msm/msw outcomes
combined_reg = rbind(msm_regressor,msw_regressor) 
combined_out = rbind(log_cytokine_frame_m,log_cytokine_frame_w)


## add index for merge
combined_reg$index = c(1:86)
combined_out$index = c(1:86)






## create a merged dataframe for msm and msw
merged_data = merge(combined_reg,combined_out,by = "index")



write.csv(merged_data, file = "final_reg_out.csv")


## replace all inf values by 0
merged_data$tnfa[which(!is.finite(merged_data$tnfa))] = 0
merged_data$ip10[which(!is.finite(merged_data$ip10))] = 0


## linear regression for msm

il10_lr = lm(il10 ~ d5__lactobacillus + d5__corynebacterium1+ d5__gardnerella +d5__sneathia +d5__finegoldia+d5__anaerococcus+d5__peptoniphilus+d5__staphylococcus+d5__streptococcus+d5__veillonella+Other , data = merge_msw)
summary(il10_lr)


il8_lr = lm(il8 ~ d5__lactobacillus + d5__corynebacterium1+ d5__gardnerella +d5__sneathia +d5__finegoldia+d5__anaerococcus+d5__peptoniphilus+d5__staphylococcus+d5__streptococcus+d5__veillonella+Other , data = merge_msw)
summary(il8_lr)

il1b_lr = lm(il1b ~ d5__lactobacillus + d5__corynebacterium1+ d5__gardnerella +d5__sneathia +d5__finegoldia+d5__anaerococcus+d5__peptoniphilus+d5__staphylococcus+d5__streptococcus+d5__veillonella+Other , data = merge_msw)
summary(il1b_lr)


tnfa_lr = lm(tnfa ~ d5__lactobacillus + d5__corynebacterium1+ d5__gardnerella +d5__sneathia +d5__finegoldia+d5__anaerococcus+d5__peptoniphilus+d5__staphylococcus+d5__streptococcus+d5__veillonella+Other , data = merge_msw)
summary(tnfa_lr)

ip10_lr = lm(ip10 ~ d5__lactobacillus + d5__corynebacterium1+ d5__gardnerella +d5__sneathia +d5__finegoldia+d5__anaerococcus+d5__peptoniphilus+d5__staphylococcus+d5__streptococcus+d5__veillonella+Other , data = merge_msw)
summary(ip10_lr)


## linear regression for merged_total
feat = colnames(msm_regressor)

il10_lr = lm(il10 ~ d5__lactobacillus + d5__corynebacterium1+ d5__gardnerella +d5__sneathia +d5__finegoldia+d5__anaerococcus+d5__peptoniphilus+d5__staphylococcus+d5__streptococcus+d5__veillonella+Other , data = merged_data)
summary(il10_lr)


il8_lr = lm(il8 ~ d5__lactobacillus + d5__corynebacterium1+ d5__gardnerella +d5__sneathia +d5__finegoldia+d5__anaerococcus+d5__peptoniphilus+d5__staphylococcus+d5__streptococcus+d5__veillonella+Other , data = merged_data)
summary(il8_lr)

il1b_lr = lm(il1b ~ d5__lactobacillus + d5__corynebacterium1+ d5__gardnerella +d5__sneathia +d5__finegoldia+d5__anaerococcus+d5__peptoniphilus+d5__staphylococcus+d5__streptococcus+d5__veillonella+Other , data = merged_data)
summary(il1b_lr)


tnfa_lr = lm(tnfa ~ d5__lactobacillus + d5__corynebacterium1+ d5__gardnerella +d5__sneathia +d5__finegoldia+d5__anaerococcus+d5__peptoniphilus+d5__staphylococcus+d5__streptococcus+d5__veillonella+Other , data = merged_data)
summary(tnfa_lr)

ip10_lr = lm(ip10 ~ d5__lactobacillus + d5__corynebacterium1+ d5__gardnerella +d5__sneathia +d5__finegoldia+d5__anaerococcus+d5__peptoniphilus+d5__staphylococcus+d5__streptococcus+d5__veillonella+Other , data = merged_data)
summary(ip10_lr)


## add indicator variable  for msm_msw
a = c(rep(1,43),rep(0,43))
merged_data$status = a
merged_data$status


## logistic regression for msm and msw using regressors
status_bac_lr = glm(status ~ d5__lactobacillus + d5__corynebacterium1+ d5__gardnerella +d5__sneathia +d5__finegoldia+d5__anaerococcus+d5__peptoniphilus+d5__staphylococcus+d5__streptococcus+d5__veillonella+Other , data = merged_data,family = "binomial")
summary(status_bac_lr)

status_cyt_lr = glm(status ~ tnfa + il8 + il1b + il10 + ip10,data = merged_data,family = "binomial" )
summary(status_cyt_lr)



## carry out linear regression using bacteria and status as regressors and cytokines as response

il10_lr_s = lm(il10 ~ status + d5__lactobacillus + d5__corynebacterium1+ d5__gardnerella +d5__sneathia +d5__finegoldia+d5__anaerococcus+d5__peptoniphilus+d5__staphylococcus+d5__streptococcus+d5__veillonella+Other , data = merged_data)
summary(il10_lr_s)
## forward selection of 


il8_lr_s = lm(il8 ~ status + d5__lactobacillus + d5__corynebacterium1+ d5__gardnerella +d5__sneathia +d5__finegoldia+d5__anaerococcus+d5__peptoniphilus+d5__staphylococcus+d5__streptococcus+d5__veillonella+Other , data = merged_data)
summary(il8_lr_s)

il1b_lr_s = lm(il1b ~ status + d5__lactobacillus + d5__corynebacterium1+ d5__gardnerella +d5__sneathia +d5__finegoldia+d5__anaerococcus+d5__peptoniphilus+d5__staphylococcus+d5__streptococcus+d5__veillonella+Other , data = merged_data)
summary(il1b_lr_s)


tnfa_lr_s = lm(tnfa ~ status + d5__lactobacillus + d5__corynebacterium1+ d5__gardnerella +d5__sneathia +d5__finegoldia+d5__anaerococcus+d5__peptoniphilus+d5__staphylococcus+d5__streptococcus+d5__veillonella+Other , data = merged_data)
summary(tnfa_lr_s)

ip10_lr_s = lm(ip10 ~  status + d5__lactobacillus + d5__corynebacterium1+ d5__gardnerella +d5__sneathia +d5__finegoldia+d5__anaerococcus+d5__peptoniphilus+d5__staphylococcus+d5__streptococcus+d5__veillonella+Other , data = merged_data)
summary(ip10_lr_s)






### question 1 :  difference in cytokines MSM vs MSW
qs1x = rbind(cytokine_frame_m,cytokine_frame_w)
qs1x$status = c(rep(0,43),rep(1,43))

qs1x_tnfa = lm(tnfa ~ status, data = qs1x)
summary(qs1x_tnfa)

qs1x_il10 = lm(il10 ~ status, data = qs1x)
summary(qs1x_il10)

qs1x_il1b = lm(il1b ~ status, data = qs1x)
summary(qs1x_il1b)

qs1x_il8 = lm(il8 ~ status, data = qs1x)
summary(qs1x_il8)

qs1x_ip10 = lm(ip10 ~ status, data = qs1x)
summary(qs1x_ip10)


### question 2 does bacteria account for the cytokines
## read all 42 features and create a regression model

msm_bact_all  = read.csv(file = '../Documents/GitHub/Microbiome-Analysis-Project/data/msm_tt_rel_abundance_all.csv', head = TRUE, sep = ",")
msw_bact_all = read.csv(file = '../Documents/GitHub/Microbiome-Analysis-Project/data/msw_tt_rel_abundance_all.csv', head = TRUE, sep = ",")

## centered log ratio transformation

library(compositions)

for(i in 1:43){
  msm_bact_all[i,] = as.list(clr(msm_bact_all[i,]))
  msw_bact_all[i,] = as.list(clr(msw_bact_all[i,]))
}

# join msm/msw bacteria
bact_all = rbind(msm_bact_all,msw_bact_all)


# index the two types of bacteria
bact_all$index = c(1:86)
qs1x$index = c(1:86)


## merge bacteria and outcome
qs2 = merge(bact_all,qs1x, by = "index")

features = paste("feature",1:42,sep = "")

## regression


qs2_tnfa = lm(
  as.formula(paste(paste(colnames(qs2)[43]), "~",
                   paste(colnames(qs2)[c(1:42)], collapse = "+"),
                   sep = ""
  )),
  data=qs2
)
summary(qs2_tnfa)


qs2_il1b = lm(
  as.formula(paste(paste(colnames(qs2)[44]), "~",
                   paste(colnames(qs2)[c(1:42)], collapse = "+"),
                   sep = ""
  )),
  data=qs2
)
summary(qs2_il1b)


qs2_il8 = lm(
  as.formula(paste(paste(colnames(qs2)[45]), "~",
                   paste(colnames(qs2)[c(1:42)], collapse = "+"),
                   sep = ""
  )),
  data=qs2
)
summary(qs2_il8)


qs2_il10 = lm(
  as.formula(paste(paste(colnames(qs2)[46]), "~",
                   paste(colnames(qs2)[c(1:42)], collapse = "+"),
                   sep = ""
  )),
  data=qs2
)
summary(qs2_il10)

qs2_ip10 = lm(
  as.formula(paste(paste(colnames(qs2)[47]), "~",
                   paste(colnames(qs2)[c(1:42)], collapse = "+"),
                   sep = ""
  )),
  data=qs2
)
summary(qs2_ip10)



### question 3 does bacteria account for cytokines when controlled for status
qs3_tnfa = lm(as.formula(paste(paste(colnames(qs2)[43]), "~", 
          paste(paste(colnames(qs2)[c(2:42)],collapse ="+"),
                paste(colnames(qs2[48]),colnames(qs2)[c(2:42)],sep="*",collapse="+"),
                sep = "+"),sep = "")),data = qs2)
summary(qs3_tnfa)


qs3_il1b = lm(as.formula(paste(paste(colnames(qs2)[44]), "~", 
                               paste(paste(colnames(qs2)[c(2:42)],collapse ="+"),
                                     paste(colnames(qs2[48]),colnames(qs2)[c(2:42)],sep="*",collapse="+"),
                                     sep = "+"),sep = "")),data = qs2)
summary(qs3_il1b)


qs3_il8 = lm(as.formula(paste(paste(colnames(qs2)[45]), "~", 
                               paste(paste(colnames(qs2)[c(2:42)],collapse ="+"),
                                     paste(colnames(qs2[48]),colnames(qs2)[c(2:42)],sep="*",collapse="+"),
                                     sep = "+"),sep = "")),data = qs2)
summary(qs3_il8)



qs3_il10 = lm(as.formula(paste(paste(colnames(qs2)[46]), "~", 
                              paste(paste(colnames(qs2)[c(2:42)],collapse ="+"),
                                    paste(colnames(qs2[48]),colnames(qs2)[c(2:42)],sep="*",collapse="+"),
                                    sep = "+"),sep = "")),data = qs2)
summary(qs3_il10)


qs3_ip10 = lm(as.formula(paste(paste(colnames(qs2)[47]), "~", 
                               paste(paste(colnames(qs2)[c(2:42)],collapse ="+"),
                                     paste(colnames(qs2[48]),colnames(qs2)[c(2:42)],sep="*",collapse="+"),
                                     sep = "+"),sep = "")),data = qs2)
summary(qs3_ip10)


### question 4





