
##########
###########  maternal height vs gestational age study  (2017)
##########   => PTD risk

# estimate the most extreme risk generating threshold values (for height and BW zscore) and validate it.
# human-friendly results and illustrations (risks, means, p-values)

library(binom)

##### 1) PREPARE DATA
###########################################################################

load("~/Biostuff/SEMFR_HEIGHT-GESTAGE/working_dataset_TWISIN_20170613_n1159478.RData")
dat = dat[which(dat$PARITET_F==1),]
dat = dat[which((dat$excl==0)|(is.na(dat$excl))),]
dat = dat[which(dat$MALDER %in% 18:45),]

twi = dat[dat$status=="twins",]
sin = dat[dat$status=="singl",]
sin = sin[which((!is.na(sin$BVIKTBS))&(!is.na(sin$GRDBS))),]


marsal = function(ga,birthweight,sex){
        if(is.na(sex)) sex=0
        if(sex==1) {  ## boys
                mw=-1.907345e-6*ga^4 + 1.140644e-3*ga^3 - 1.336265e-1*ga^2 + 1.976961*ga + 2.410053e+2
        } else if(sex==2) { ## girls
                mw=-2.761948e-6*ga^4 + 1.744841e-3*ga^3 - 2.893626e-1*ga^2 + 18.91197*ga - 4.135122e+2
        } else { ## ?
                mw=-2.278843e-6*ga^4 + 1.402168e-3*ga^3 - 2.008726e-1*ga^2 + 9.284121*ga - 4.125956e+1
        }
        prb = pnorm(birthweight, mean = mw, sd = abs(0.12*mw)) # probability
        qnorm(prb)  # z-score
}
sin$bw_zsc_Marsal=sapply(1:nrow(sin), function(x) marsal(sin$GRDBS[x],sin$BVIKTBS[x],sin$KON[x])) # sex "unknown"
twi$bw_zsc_Marsal=NA

# split data into testing and validation subsets, 
# cherry-pick the highest PTD risk in dataset 1, 
# validate it in dataset 2

### define datasets

twi = twi[order(twi$MLANGD),]
grps = rep(c(0,1),floor(nrow(twi)/2))
if(length(grps)!=nrow(twi)) grps = c(grps,0)
twi$grps = grps; rm(grps)

sin = sin[order(sin$MLANGD),]
grps = rep(c(0,1),floor(nrow(sin)/2))
if(length(grps)!=nrow(sin)) grps = c(grps,0)
sin$grps = grps; rm(grps)

dat = rbind(sin,twi)
set.seed(1)
dat = dat[sample(nrow(dat),nrow(dat),replace = F),]

tra = dat[which(dat$grps==0),]  # training dataset
val = dat[which(dat$grps==1),]  # validation dataset

### OPTIONAL (if no training/validation is being used): 
tra = dat
val = dat

######################################################################   FETAL GROWTH RATE  (model 2)

###  chosen maternal height threshold is 160cm (inclusive)
mh_thr = 161

### find the BWzscore threshold

### find most extreme mrisk ratio (using TEST DATASET)
sd_thrs = seq(0,3,0.05)  # birtheight z-score threshold (>)
rez = NULL
for (sd_thr in sd_thrs) {
                ix_rsk = which((tra$status=="singl")&(tra$bw_zsc_Marsal>sd_thr)&(tra$MLANGD<mh_thr))  # risk group
                ix_ref = which((tra$status=="singl")&((tra$bw_zsc_Marsal<=lga_thr)|(tra$MLANGD>=mh_thr)))  # risk group
                
                
                mnga_rsk = mean(tra$GRDBS[ix_rsk])
                rptd_rsk = mean(tra$GRDBS[ix_rsk]<259)
                Nmin_rsk = sum(tra$GRDBS[ix_rsk]<259)
                
                s1 = sum(tra$GRDBS[ix_rsk]<259)
                s2 = sum(tra$GRDBS[ix_rsk]>=259)
                s3 = sum(tra$GRDBS[ix_ref]<259)
                s4 = sum(tra$GRDBS[ix_ref]>=259)
                odds_ratio = (s1/s2)/(s3/s4)
                
                tmp = data.frame(mh=mh_thr,lga=sd_thr,Nmin_rsk,mnga_rsk,rptd_rsk,odds_ratio,stringsAsFactors = F)
                rez = rbind(rez,tmp)
                rm(ix_rsk,ix_ref,mnga_rsk,rptd_rsk,Nmin_rsk,s1,s2,s3,s4,odds_ratio,tmp)
        }
rez = rez[which(rez$Nmin_rsk>5),]
plot(rez$odds_ratio ~ rez$lga,xlab="SD threshold to define LGA (>)",ylab="oddsRatio PTD",type="l")
grid()
plot(rez$rptd_rsk ~ rez$lga,xlab="SD threshold to define LGA (>)",ylab="PTD risk",type="l")
grid()


### choose he most extreme value
rez[which(rez$odds_ratio==max(rez$odds_ratio)),]
#mh_thr = 161
lga_thr = 2.75 #2.9

# estimate interaction value using LOGISTIC regression
sub =    val[which(val$status=="singl"),]
sub$PTD = as.numeric(sub$GRDBS<259)
sub$GRWTH = as.numeric(sub$bw_zsc_Marsal>lga_thr)
sub$MHGH = as.numeric(sub$MLANGD<mh_thr)
sub = sub[,c("PTD","GRWTH","MHGH")]
sub$gr = sub$MHGH + sub$GRWTH*2  #  explicit categories
table(sub$gr)
head(sub)
df = group_by(sub,gr) %>% summarise(nPTD=sum(PTD),nTRM=sum(PTD==0),nTotal=n(),ptdRisk=mean(PTD)) %>% ungroup()
df = as.data.frame(df)
tmp = binom.confint(df$nPTD,df$nTotal,conf.level = 0.95, methods = "agresti-coull")
df$lower = tmp$lower
df$upper = tmp$upper
df$ptdRisk = df$ptdRisk*100
df$lower = df$lower*100
df$upper = df$upper*100


odds = NULL
for (i in 4:2) {
m = df[c(i,1),2:3]
ftst = fisher.test(m)
cint = as.numeric(ftst$conf.int)
tmp = data.frame(gr=i,OR=as.numeric(ftst$estimate),ciLow=cint[1],ciUpp=cint[2],stringsAsFactors = F)
odds = rbind(odds,tmp); rm(m,ftst,cint,tmp)
}
odds





mod0 = glm(PTD ~ MHGH*GRWTH,data=sub,family = "binomial")
cfs0 = coef(summary(mod0)); cfs0  # interaction P = 4.623 e-4
exp(cfs0[,1])

mod = glm(PTD ~ factor(gr),data=sub,family = "binomial")
cfs = coef(summary(mod)); cfs
exp(cfs[,1])
sum(cfs0[2:4,1])

#### how would the non-interaction in teh 4'th group look like?
fun = function() {
        gr1 = rnorm(1e4,mean=cfs[2,1],sd=cfs[2,2]) # maternal height
        gr2 = rnorm(1e4,mean=cfs[3,1],sd=cfs[3,2]) # uterine load
        gr12 = exp(gr1)*exp(gr2) # odds ratios are multiplied
        quantile(gr12,probs = c(0.025,0.5,0.975)) # 95% CI adn the midpoint OR
}
rez = replicate(100,fun())
noIntX = apply(rez,1,median) # expected OddsRatio without interaction


sq = 2:4
plot(NA,ylim=c(1.7,5),xlim=c(0,max(odds$ciUpp)*1.1),yaxt="n",bty="n",xlab="PTD odds ratio")
#### add shaded are to teh plot
xs = c(noIntX[1],noIntX[3],noIntX[3],noIntX[1],noIntX[1])
ys = c(2-0.2,2-0.2,2+0.2,2+0.2,2-0.2) 
polygon(xs,ys,border = NA,angle = 60,density = 15,col="grey")
#segments(y0 = 2-0.1,y1=2+.1,x0=noIntX[2],x1=noIntX[2],lwd=3,col="darkgrey")
segments(y0 = sq-0.05,y1=sq+0.05,x0=odds$OR,x1=odds$OR,lwd=2)
segments(y0 = sq,y1=sq,x0=odds$ciLow,x1=odds$ciUpp,col="black",lwd=2)
abline(v=1,col="grey")


##########################################
#### investigate the interaction in theory
# set.seed(1754)
# fun  = function(x) {
# n = 1e4
# a = data.frame(PTD=sample(c(0,1),n,replace=T),MHGH=0,GRWTH=0,stringsAsFactors = F)
# ix1 = sample(which(a$PTD==1),n*0.5*0.5)
# ix0 = sample(which(a$PTD==0),n*0.5*0.2)
# a$MHGH[c(ix1,ix0)] = 1
# ix1 = sample(which(a$PTD==1),n*0.5*0.5)
# ix0 = sample(which(a$PTD==0),n*0.5*0.2)
# a$GRWTH[c(ix1,ix0)] = 1
# a$gr = a$MHGH + a$GRWTH*2  #  explicit categories
# 
# df = group_by(a,gr) %>% summarise(nPTD=sum(PTD),nTRM=sum(PTD==0),nTotal=n(),ptdRisk=mean(PTD)) %>% ungroup()
# df = as.data.frame(df)
# tmp = binom.confint(df$nPTD,df$nTotal,conf.level = 0.95, methods = "agresti-coull")
# 
# mod = glm(PTD ~ MHGH*GRWTH,data=a,family = "binomial")
# cfs = coef(summary(mod)) # cfs
# #return(cfs[4,4])
# return(df[4,5])
# }
# 
# rsk = replicate(100,fun())
# hist(rsk)
# 
# mod = glm(PTD ~ factor(gr),data=a,family = "binomial")
# cfs = coef(summary(mod)); cfs
# exp(cfs[,1])
# sum(cfs[2:4,1])
# 
# 
# df = group_by(a,gr) %>% summarise(nPTD=sum(PTD),nTRM=sum(PTD==0),nTotal=n(),ptdRisk=mean(PTD)) %>% ungroup()
# df = as.data.frame(df)







######################################################################   TWIN PREGNANCIES  (model 1)

###  chosen maternal height threshold is 160cm (inclusive)
mh_thr = 161 #161

sub =    val
sub$PTD = as.numeric(sub$GRDBS<259) #259
sub$TWIN = as.numeric(sub$status=="twins") 
sub$MHGH = as.numeric(sub$MLANGD<mh_thr) 
sub = sub[,c("PTD","TWIN","MHGH")]
sub$gr = sub$MHGH + sub$TWIN*2  #  explicit categories


df = group_by(sub,gr) %>% summarise(nPTD=sum(PTD),nTRM=sum(PTD==0),nTotal=n(),ptdRisk=mean(PTD)) %>% ungroup()
df = as.data.frame(df)
tmp = binom.confint(df$nPTD,df$nTotal,conf.level = 0.95, methods = "agresti-coull")
df$lower = tmp$lower
df$upper = tmp$upper
df$ptdRisk = df$ptdRisk*100
df$lower = df$lower*100
df$upper = df$upper*100

odds = NULL
for (i in 4:2) {
        m = df[c(i,1),2:3]
        ftst = fisher.test(m)
        cint = as.numeric(ftst$conf.int)
        tmp = data.frame(gr=i,OR=as.numeric(ftst$estimate),ciLow=cint[1],ciUpp=cint[2],stringsAsFactors = F)
        odds = rbind(odds,tmp); rm(m,ftst,cint,tmp)
}
odds






head(sub)
mod0 = glm(PTD ~ MHGH*TWIN,data=sub,family = "binomial")
summary(mod0)
#mod1 = lm(PTD ~ MHGH*TWIN,data=sub)  # strangely, this one shows signif interaction!
#summary(mod1)


cfs0 = coef(summary(mod0)); cfs0  
exp(cfs0[,1])

mod = glm(PTD ~ factor(gr),data=sub,family = "binomial")
cfs = coef(summary(mod)); cfs

#### how would the non-interaction in teh 4'th group look like?
fun = function() {
        gr1 = rnorm(1e4,mean=cfs[2,1],sd=cfs[2,2]) # maternal height
        gr2 = rnorm(1e4,mean=cfs[3,1],sd=cfs[3,2]) # uterine load
        gr12 = exp(gr1)*exp(gr2) # odds ratios are multiplied
        quantile(gr12,probs = c(0.025,0.5,0.975)) # 95% CI adn the midpoint OR
}
rez = replicate(100,fun())
noIntX = apply(rez,1,median) # expected OddsRatio without interaction


sq = 2:4
plot(NA,ylim=c(1.7,5),xlim=c(0,max(odds$ciUpp)*1.05),yaxt="n",bty="n",xlab="PTD odds ratio")
#### add shaded are to teh plot
xs = c(noIntX[1],noIntX[3],noIntX[3],noIntX[1],noIntX[1])
ys = c(2-0.2,2-0.2,2+0.2,2+0.2,2-0.2) 
polygon(xs,ys,border = NA,angle = 60,density = 15,col="grey")
#segments(y0 = 2-0.1,y1=2+.1,x0=noIntX[2],x1=noIntX[2],lwd=3,col="darkgrey")
segments(y0 = sq-0.05,y1=sq+0.05,x0=odds$OR,x1=odds$OR,lwd=2)
segments(y0 = sq,y1=sq,x0=odds$ciLow,x1=odds$ciUpp,col="black",lwd=2)
abline(v=1,col="grey")

