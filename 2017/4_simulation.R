# project: "Maternal height - gestational age" study  (aka uterine distention hypothesis)
# description: a simulation of predicted interaction pattern (GA=MH*UL)
# method: simulates births occuring due to two competing birth-risks (baseline and interaction)
# baseline risks - any other (than interaction) genetic and environmental factors that trigger the birth
# interaction risk - the one due to the space constraint arising from fetal growth and uterine size
# created: 2017-10-02
# reviewed: 2018-09-01
# author: Jonas Bacelis

## load packages:
library(flexsurv)
library(ggplot2)

## load functions "fun_interactionPlots", "fun_gradHist", "fun_gradHist", "fun_gradHist"
source("simulation_functions.R") # dir: ~/Dropbox/GIT/SEMFR_HEIGHT-GESTAGE/2017/

## parameters:
n_ind = 1e4  # sample size

fun_run = function(base_shift,densit,scaPlot,graHistA,graHistB,graHistX) {
# base_shift = how many days earlier births would be due to baseline risk, as compared to interaction
# densit = to-plot-or-not the densities (logical TRUE/FALSE)
# scaPlot = to-plot-or-not main scatterplot (logical TRUE/FALSE) 
# graHistA = to-plot-or-not final simulated gestational age due to competing baseline and interaction risks (logical TRUE/FALSE) 
# graHistB = to-plot-or-not counterfactual gestational age due to baseline risks (logical TRUE/FALSE) 
# graHistX = to-plot-or-not counterfactual gestational age due to interaction (logical TRUE/FALSE) 
                
## generate variables
gas = rgompertz(n_ind*1.1,0.1212,1.026*1e-8) + 150 # approximation of natural gest.age-at-birth distribution
# truncate the tail (for not wasting precious gradient colors on invisible bars in histogram)
gas = gas[which(gas>220)]
gas = sample(gas,n_ind,replace = F)

fgr = rnorm(n_ind,mean=10,sd=1) # fetal growth rate (arbitrary units)
mus = rnorm(n_ind,mean=10,sd=1) # maternal uterine size (arbitrary units)
rat =  fgr/mus  # ratio, determining interaction-imposed birth time

# create dependency between maternal height and uterine size
mhs_src = rnorm(n_ind,mean=165,sd=10) # maternal height (source)
hlp = mus + rnorm(n_ind) # helper
rnk = match(rank(mhs_src),rank(hlp))
mhs = mhs_src[order(rnk)]; rm(rnk,hlp,mhs_src)

# interaction("rat")-imposed getational age data (gas_X)
tmp = (length(gas)+1) - rank(gas)
rnk = match(tmp,rank(rat))
gas_X = gas[order(rnk)]; rm(rnk)

# baseline-imposed getational age data (gas_B)
gas_B = sample(gas)  + base_shift 

# create a dataframe 
dat = data.frame(fgr,mus,mhs,rat,gas_B,gas_X,stringsAsFactors = F)

### determine the actual Gestational age: whichever risk triggered the birth first
dat$gas = apply(dat[,c("gas_B","gas_X")],1,min)
#dat$gas = dat$gas_B
#dat$gas = dat$gas_X

# two distict groups of fast and slow fetal growth
dat$sta = NA
dat$sta[which(dat$fgr<=quantile(dat$fgr,probs = 0.1))] = 0  # SGA
dat$sta[which(dat$fgr>=quantile(dat$fgr,probs = 0.9))] = 1  # LGA

mod = lm(gas ~ sta*mhs,data=dat)
ix = which(!is.na(dat$sta))

if (densit==TRUE) {
d1 = density(dat$gas_B)
d2 = density(dat$gas_X)
plot(d1,lwd=2,col="olivedrab",xlim=range(c(dat$gas_B,dat$gas_X)),
     xlab="gestational age at birth")
points(d2,lwd=2,type="l",col="orange")
}

if (scaPlot==TRUE) {
fun_interactionPlots(mod,mhs = dat$mhs[ix],sta = dat$sta[ix],gas = dat$gas[ix],1)
}
if (graHistA==TRUE) {
        fun_gradHist(dat$gas,dat$rat,20,"gestational age at birth","counts")
}
if (graHistB) {
fun_gradHist(dat$gas_B,dat$rat,20,"gestational age at birth","counts")
}
if (graHistX) {
fun_gradHist(dat$gas_X,dat$rat,20,"gestational age at birth","counts")
}

mod = lm(gas ~ sta*mhs,data=dat)
print(coef(summary(mod)))

mod = lm(gas ~ fgr*mhs,data=dat)
print(coef(summary(mod)))
} # end of function



# run the simulation:
par(mfrow=c(1,1))
fun_run(0,F,T,F,F,F) # base_shift,densit,scaPlot,graHistA,graHistB,graHistX

        





####################
### counterfactual scenario: 
### Maternal uterine size and fetal growth rate 
### manifest their effects independently
# i.e., what pattern should we expect if there is no physical interaction
####################


n_ind = 5e3
gas = round(rgompertz(n_ind,0.1212,1.026*1e-8) + 150,0)
fgr = rnorm(n_ind,mean=10,sd=1)
mus = rnorm(n_ind,mean=10,sd=1)

# rat = ratio  (in this case this is not a fgr/mus ratio at all)
rat = rep(NA,n_ind)
rnd = sample(n_ind,n_ind/2,replace = F)
rat[rnd] = fgr[rnd] # introduce correlation with fetal growth rate
rat[-rnd] = mus[-rnd]  # introduce correlation with maternal uterine size

# create dependence between maternal height and uterine size
mhs_src = rnorm(n_ind,mean=165,sd=10)
hlp = mus + rnorm(n_ind)
rnk = match(rank(mhs_src),rank(hlp))
mhs = mhs_src[order(rnk)]; rm(rnk)
hist(mhs)

dat = data.frame(fgr,mus,mhs,rat,GAbase=gas,stringsAsFactors = F)
dat = dat[order(dat$rat,decreasing = T),]
dat$GAintrx = sort(gas)

dat$GA = apply(dat[,c("GAbase","GAintrx")],1,min)
# what was the cause of birth?
dat$cause = NA
dat$cause[which(dat$GAintrx < dat$GAbase)] = "X" # births due to interaction
dat$cause[which(dat$GAbase < dat$GAintrx)] = "B" # births due to baseline risks
dat$cause[which(dat$GAintrx == dat$GAbase)] = "N" # rare event when no particular risk can be blamed
table(dat$cause)

qnts = quantile(dat$fgr,probs = c(0.1,0.9))
dat$sta = NA
dat$sta[which(dat$fgr<=qnts[1])] = 0  # SGA
dat$sta[which(dat$fgr>=qnts[2])] = 1  # LGA

summary(lm(GA ~ fgr*mhs,data=dat))
summary(lm(GA ~ sta*mhs,data=dat))

ggplot(dat, aes(GA, fill = cause)) + geom_histogram(bins=20)

source("simulation_functions.R") # ~/Dropbox/GIT/SEMFR_HEIGHT-GESTAGE/2017/
mod = lm(GA ~ sta*mhs,data=dat)
summary(mod)
fun_interactionPlots(mod,dat$mhs,dat$sta,dat$GA,1)
summary(mod)
table(dat$sta)

summary(lm(rat ~ fgr))
summary(lm(rat ~ mus))
summary(lm(rat ~ fgr + mus))
summary(lm(rat ~ fgr * mus))



