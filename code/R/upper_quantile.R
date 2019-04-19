

out_msm = read.csv(file = "../Documents/GitHub/Microbiome-Analysis-Project/data/MSM_urinary_cytokines.csv",header = TRUE, sep = ",")
out_msw = read.csv(file = "../Documents/GitHub/Microbiome-Analysis-Project/data/MSW_urinary_cytokines.csv",header = TRUE, sep = ",")
cytokines = c('tnfa','il1b','il8','il10','ip10')


msm_cyt = out_msm[cytokines]
msw_cyt = out_msw[cytokines]


## tnfa msw
for (i in 1:43){
  if (msw_cyt$tnfa[i] > quantile(msw_cyt$tnfa)[4]){
    msw_cyt$tnfa[i] = 1
  }
  else msw_cyt$tnfa[i] = 0
}

## il1b msw
for (i in 1:43){
  if (msw_cyt$il1b[i] > quantile(msw_cyt$il1b)[4]){
    msw_cyt$il1b[i] = 1
  }
  else msw_cyt$il1b[i] = 0
}

## il8 msw
for (i in 1:43){
  if (msw_cyt$il8[i] > quantile(msw_cyt$il8)[4]){
    msw_cyt$il8[i] = 1
  }
  else msw_cyt$il8[i] = 0
}

## il10 msw
for (i in 1:43){
  if (msw_cyt$il10[i] > quantile(msw_cyt$il10)[4]){
    msw_cyt$il10[i] = 1
  }
  else msw_cyt$il10[i] = 0
}

## ip10 msw
for (i in 1:43){
  if (msw_cyt$ip10[i] > quantile(msw_cyt$ip10)[4]){
    msw_cyt$ip10[i] = 1
  }
  else msw_cyt$ip10[i] = 0
}


## tnfa msm
for (i in 1:43){
  if (msm_cyt$tnfa[i] > quantile(msm_cyt$tnfa)[4]){
    msm_cyt$tnfa[i] = 1
  }
  else msm_cyt$tnfa[i] = 0
}

## il1b msm
for (i in 1:43){
  if (msm_cyt$il1b[i] > quantile(msm_cyt$il1b)[4]){
    msm_cyt$il1b[i] = 1
  }
  else msm_cyt$il1b[i] = 0
}

## il8 msm
for (i in 1:43){
  if (msm_cyt$il8[i] > quantile(msm_cyt$il8)[4]){
    msm_cyt$il8[i] = 1
  }
  else msm_cyt$il8[i] = 0
}

## il10 msm
for (i in 1:43){
  if (msm_cyt$il10[i] > quantile(msm_cyt$il10)[4]){
    msm_cyt$il10[i] = 1
  }
  else msm_cyt$il10[i] = 0
}


## ip10 msm
for (i in 1:43){
  if (msm_cyt$ip10[i] > quantile(msm_cyt$ip10)[4]){
    msm_cyt$ip10[i] = 1
  }
  else msm_cyt$ip10[i] = 0
}

write.csv(msm_cyt,"../Documents/GitHub/Microbiome-Analysis-Project/data/MSM_Cyt_Upper_Quantile.csv")
write.csv(msw_cyt,"../Documents/GitHub/Microbiome-Analysis-Project/data/MSW_Cyt_Upper_Quantile.csv")

