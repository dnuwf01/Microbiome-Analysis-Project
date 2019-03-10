# this is a sample r code to read a csv file and then create a spearman heatmap out of it.

msw <- read.csv(file='C:/Users/Debarghya Nandi/Desktop/research/mehta/deesha_dulal_mehta/msw43.csv', header=TRUE,sep=",")
msw <- data.frame(msw)
head(msw)

# select the specific columns from the data
feature_list = c('d5__lactobacillus','d5__corynebacterium1','d5__gardnerella','d5__sneathia','d5__finegoldia','d5__anaerococcus','d5__peptoniphilus','d5__staphylococcus','d5__streptococcus',
                 'd5__veillonella','Other')

msw <- msw[feature_list]
head(msw)


# centered log transformation
library(compositions)
(tmp <- clr(msw))
(tmp)



cormat <- round(cor(tmp, use="complete.obs", method="spearman"),2)
head(cormat)





library(reshape2)
melted_cormat <- melt(cormat, na.rm = TRUE)


library(ggplot2)
ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + geom_tile(color="white")+  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1,1), space = "Lab", name="Spearman\nCorrelation")+   theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()  



