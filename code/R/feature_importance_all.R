## Random Forest and Lasso applied to all the selected taxa after filtering
## relative importance of each of the feature is obtained


### @ Deb

df_msm  = read.csv(file = '../Documents/GitHub/Microbiome-Analysis-Project/data/msm_tt_rel_abundance_all.csv', head = TRUE, sep = ",")
df_msw = read.csv(file = '../Documents/GitHub/Microbiome-Analysis-Project/data/msw_tt_rel_abundance_all.csv', head = TRUE, sep = ",")





### Select taxa that accounts for appx 90% of the total relative abundance.
##MSM 
taxa_90p = c('d5__lactobacillus','d5__corynebacterium1','d5__gardnerella','d5__sneathia',
             'd5__finegoldia','d5__anaerococcus','d5__peptoniphilus','d5__staphylococcus',
             'd5__streptococcus','d5__veillonella','d5__prevotella','d5__ezakiella',
             'd5__shuttleworthia','d4__prevotellaceaeother','d5__atopobium')



df_msm = df_msm[taxa_90p]
df_msw = df_msw[taxa_90p]


## add a status to the dataframes and then  merge it

df_msm$status = rep(1,43)
df_msw$status = rep(0,43)


df_combined = rbind(df_msm,df_msw)

## replace all missing values by zero
df_combined[is.na(df_combined)] = 0







## random forest

rf = randomForest(status ~ ., data = df_combined, importance = TRUE,proximity = TRUE)
print(rf)

round(importance(rf),2)

library(glmnet)
## lasso 
drops = c("status")
x = as.matrix(df_combined[, !(names(df_combined) %in% drops)])
y = as.double(as.matrix(df_combined[,16]))


lasso = cv.glmnet(x,y, alpha = 1,family = "binomial",type.measure = "auc")
coef(lasso, s=lasso$lambda.min)


elnet = cv.glmnet(x,y, alpha = 0.5,family = "binomial",type.measure = "auc")
coef(elnet, s=elnet$lambda.min)


## regression of bacteria vs cytokines
## ----------------------------cytokines ----------------------------------
#read the csv file from the given location
## change the file name according to what need (continuous or upper quantile)
out_msm = read.csv(file = "../Documents/GitHub/Microbiome-Analysis-Project/data/MSM_urinary_cytokines.csv",header = TRUE, sep = ",")
out_msw = read.csv(file = "../Documents/GitHub/Microbiome-Analysis-Project/data/MSW_urinary_cytokines.csv",header = TRUE, sep = ",")
cytokines = c('tnfa','il1b','il8','il10','ip10')


## select the important cytokines
cytokine_frame_m = out_msm[cytokines]
cytokine_frame_w= out_msw[cytokines]




## centered log-ratio transform for regressors along each row
library(compositions)

for(i in 1:43){
  df_msm[i,] = as.list(clr(df_msm[i,]))
  df_msw[i,] = as.list(clr(df_msw[i,]))
}

##
## add a status bar
df_msm$status = rep(1,43)
df_msw$status = rep(0,43)




## log transform the cytokine data
## NOT REQUIRED for upper quantile
log_cytokine_frame_m = log(cytokine_frame_m+0.5)
log_cytokine_frame_w = log(cytokine_frame_w+0.5)


## NOT NEEDED
### set index for msm
df_msm$index = c(1:43)
log_cytokine_frame_m$index = c(1:43)


## set index for msw
df_msw$index = c(1:43)
log_cytokine_frame_w$index = c(1:43)

#####




## bind the msm/msw regressors and msm/msw outcomes
combined_reg = rbind(df_msm,df_msw) 
## for continuous
combined_out = rbind(log_cytokine_frame_m,log_cytokine_frame_w)
## for upper quantile
#combined_out = rbind(cytokine_frame_m, cytokine_frame_w) 


## add index for merge
combined_reg$index = c(1:86)
combined_out$index = c(1:86)

## create a merged dataframe for msm and msw
merged_data = merge(combined_reg,combined_out,by = "index")


## remove all nan values
merged_data[is.na(merged_data)] = 0


## linear regression for msm
## for stepwise selection
library(MASS)

il10_lr = lm(il10 ~ d5__lactobacillus+d5__peptoniphilus+d5__staphylococcus + d5__ezakiella + d4__prevotellaceaeother + d5__shuttleworthia + d5__prevotella, data = merged_data)
summary(il10_lr)
stepil10 = stepAIC(il10_lr,direction = "both")
stepil10$anova


il8_lr = lm(il8 ~ d5__lactobacillus+d5__peptoniphilus+d5__staphylococcus + d5__ezakiella + d4__prevotellaceaeother + d5__shuttleworthia + d5__prevotella, data = merged_data)
summary(il8_lr)
stepil8 = stepAIC(il8_lr,direction = "both")
stepil8$anova

il1b_lr = lm(il1b ~ d5__lactobacillus+d5__peptoniphilus+d5__staphylococcus + d5__ezakiella + d4__prevotellaceaeother + d5__shuttleworthia + d5__prevotella, data = merged_data)
summary(il1b_lr)
stepil1b = stepAIC(il1b_lr,direction = "both")
stepil1b$anova




tnfa_lr = lm(tnfa ~ d5__lactobacillus+d5__peptoniphilus+d5__staphylococcus + d5__ezakiella + d4__prevotellaceaeother + d5__shuttleworthia + d5__prevotella, data = merged_data)
summary(tnfa_lr)
steptnfa = stepAIC(tnfa_lr,direction = "both")
steptnfa$anova


ip10_lr = lm(ip10 ~ d5__lactobacillus+d5__peptoniphilus+d5__staphylococcus + d5__ezakiella + d4__prevotellaceaeother + d5__shuttleworthia + d5__prevotella, data = merged_data)
summary(ip10_lr)
stepip10 = stepAIC(ip10_lr,direction = "both")
stepip10$anova


#### including status in the model
il10_lr = lm(il10 ~ status + d5__lactobacillus+d5__peptoniphilus+d5__staphylococcus + d5__ezakiella + d4__prevotellaceaeother + d5__shuttleworthia + d5__prevotella, data = merged_data)
summary(il10_lr)
stepil10 = stepAIC(il10_lr,direction = "both")
stepil10$anova


il8_lr = lm(il8 ~  status + d5__lactobacillus+d5__peptoniphilus+d5__staphylococcus + d5__ezakiella + d4__prevotellaceaeother + d5__shuttleworthia + d5__prevotella, data = merged_data)
summary(il8_lr)
stepil8 = stepAIC(il8_lr,direction = "both")
stepil8$anova


il1b_lr = lm(il1b ~ status +  d5__lactobacillus+d5__peptoniphilus+d5__staphylococcus + d5__ezakiella + d4__prevotellaceaeother + d5__shuttleworthia + d5__prevotella, data = merged_data)
summary(il1b_lr)
stepil1b = stepAIC(il1b_lr,direction = "both")
stepil1b$anova


tnfa_lr = lm(tnfa ~ status + d5__lactobacillus+d5__peptoniphilus+d5__staphylococcus + d5__ezakiella + d4__prevotellaceaeother + d5__shuttleworthia + d5__prevotella, data = merged_data)
summary(tnfa_lr)
steptnfa = stepAIC(tnfa_lr,direction = "both")
steptnfa$anova


ip10_lr = lm(ip10 ~ status +  d5__lactobacillus+d5__peptoniphilus+d5__staphylococcus + d5__ezakiella + d4__prevotellaceaeother + d5__shuttleworthia + d5__prevotella, data = merged_data)
summary(ip10_lr)
stepip10 = stepAIC(ip10_lr,direction = "both")
stepip10$anova


### including INTERACTION terms in the model
il10_lr = lm(il10 ~ status + d5__lactobacillus+d5__peptoniphilus+d5__staphylococcus + d5__ezakiella + d4__prevotellaceaeother + d5__shuttleworthia + d5__prevotella + status*d5__lactobacillus+status*d5__peptoniphilus+status*d5__staphylococcus + status*d5__ezakiella + status*d4__prevotellaceaeother + status*d5__shuttleworthia + status*d5__prevotella, data = merged_data)
summary(il10_lr)
stepil10 = stepAIC(il10_lr,direction = "both")
stepil10$anova

il8_lr = lm(il8 ~ status + d5__lactobacillus+d5__peptoniphilus+d5__staphylococcus + d5__ezakiella + d4__prevotellaceaeother + d5__shuttleworthia + d5__prevotella + status*d5__lactobacillus+status*d5__peptoniphilus+status*d5__staphylococcus + status*d5__ezakiella + status*d4__prevotellaceaeother + status*d5__shuttleworthia + status*d5__prevotella, data = merged_data)
summary(il8_lr)
stepil8 = stepAIC(il8_lr,direction = "both")
stepil8$anova



il1b_lr = lm(il1b ~ status + d5__lactobacillus+d5__peptoniphilus+d5__staphylococcus + d5__ezakiella + d4__prevotellaceaeother + d5__shuttleworthia + d5__prevotella + status*d5__lactobacillus+status*d5__peptoniphilus+status*d5__staphylococcus + status*d5__ezakiella + status*d4__prevotellaceaeother + status*d5__shuttleworthia + status*d5__prevotella, data = merged_data)
summary(il1b_lr)
stepil1b = stepAIC(il1b_lr,direction = "both")
stepil1b$anova


tnfa_lr = lm(tnfa ~ status + d5__lactobacillus+d5__peptoniphilus+d5__staphylococcus + d5__ezakiella + d4__prevotellaceaeother + d5__shuttleworthia + d5__prevotella + status*d5__lactobacillus+status*d5__peptoniphilus+status*d5__staphylococcus + status*d5__ezakiella + status*d4__prevotellaceaeother + status*d5__shuttleworthia + status*d5__prevotella, data = merged_data)
summary(tnfa_lr)
steptnfa = stepAIC(tnfa_lr,direction = "both")
steptnfa$anova

ip10_lr = lm(ip10 ~ status + d5__lactobacillus+d5__peptoniphilus+d5__staphylococcus + d5__ezakiella + d4__prevotellaceaeother + d5__shuttleworthia + d5__prevotella + status*d5__lactobacillus+status*d5__peptoniphilus+status*d5__staphylococcus + status*d5__ezakiella + status*d4__prevotellaceaeother + status*d5__shuttleworthia + status*d5__prevotella, data = merged_data)
summary(ip10_lr)
stepip10 = stepAIC(ip10_lr,direction = "both")
stepip10$anova