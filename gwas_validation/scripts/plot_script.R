library(data.table)
library(ggplot2)
library(tidyverse)
#define function for plotting
model_and_plot <- function(data,plot_title,pheno,x_label,y_label) {

#if b1 < 0, b1*-1, b2*-1
sig = sign(data$b1)
data$b1 = sig*data$b1
data$b2 = sig*data$b2
#models
if (nrow(data)>=3){ 
betamodel = lm(b2~0+b1,data=data)
summary(betamodel)
mlogpmodel = lm(mlogp2~0+mlogp1,data=data)
summary(mlogpmodel)
#extract params
br2 = summary(betamodel)$r.squared
bslope = summary(betamodel)$coefficients[1]

mr2 = summary(mlogpmodel)$r.squared
mslope = summary(mlogpmodel)$coefficients[1]
}
else {
bslope=0
br2=0
mslope=0
mr2=0
}
#plot limits
b_min = min(c( min(data$b1),min(data$b2)))
b_max = max(c( max(data$b1),max(data$b2)))
mlogp_min = min(c( min(data$mlogp1),min(data$mlogp2)))
mlogp_max = max(c( max(data$mlogp1),max(data$mlogp2)))

bplot <- ggplot(data,aes(x=b1,y=b2))+geom_point()+
geom_abline(intercept=0,slope=bslope,linetype=5)+
geom_abline(intercept=0,slope=1,linetype=2,color="gray")+
annotate("text",x=b_min+0.2*(b_max-b_min),y=b_min+0.8*(b_max-b_min),label=paste("r2:",format(round(br2,2),nsmall=2)),size=8)+
xlim(b_min,b_max)+ylim(b_min,b_max)+ggtitle(label=paste(plot_title,"Betas",pheno))+
xlab(x_label)+ylab(y_label)

mplot <- ggplot(data,aes(x=mlogp1,y=mlogp2))+geom_point()+
geom_abline(intercept=0,slope=mslope,linetype=5)+
geom_abline(intercept=0,slope=1,linetype=2,color="gray")+
annotate("text",x=mlogp_min +0.2*(mlogp_max-mlogp_min),y=mlogp_min + 0.8*(mlogp_max-mlogp_min),label=paste("r2:",format(round(mr2,2),nsmall=2)),size=8)+
xlim(mlogp_min,mlogp_max)+ylim(mlogp_min,mlogp_max)+ggtitle(label=paste(plot_title,"mlogp",pheno))+
xlab(x_label)+ylab(y_label)

print(bplot)
print(mplot)
}
pheno <- "~{phenotype}"
prev_label <- "~{prev_label}"
current_label <- "~{current_label}"

data_1e5 <- fread("~{phenotype}_1e5_current.tsv")
data_old_gws <- fread("~{phenotype}_5e8_previous.tsv")
data_new_hits <- fread("~{phenotype}_new_hits.tsv")
data_new_gws_not_in_old <- data_1e5 %>%filter(p2<5e-8 & p1>5e-8)

pdf("~{phenotype}_plot.pdf")
model_and_plot(data_1e5,"new 1e-5 vars",pheno,prev_label,current_label)
model_and_plot(data_old_gws,"Old GWS vars",pheno,prev_label,current_label)
model_and_plot(data_new_gws_not_in_old,"Vars GWS in new but not old",pheno,prev_label,current_label)
model_and_plot(data_new_hits,"New hits",pheno,prev_label,current_label)

dev.off()