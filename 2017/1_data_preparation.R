### Study: Explore the causality between Maternal Height and Gestatinal age
###  mnGestAge vs matHgh for specific subgroups: twin vs singleton pregnancies
#  rationale: twins occupy more space than singletons, thus
#  interaction of matHeight, gestAge and type of pregnancy should be detected

# 2016 Nov 01. by Jonas Bacelis
# update: 2017 Apr 10-18
# update: 2017 Jun 10-13


###########################
########################### CLEANING
########################### AND SUBSEQUENT TWIN/SINGLETON IDENTIFICATION
###########################


setwd("~/Dropbox/GIT/SEMFR_HEIGHT-GESTAGE")
#library(dplyr)
#library(tidyr)
library(RcppRoll)

# load mfr data 
load("~/Biostuff/SEMFR_DATA/verena_files/SOS/Death Cancer Patient Birthregistry 150202/mfr150202/mfr_150202_maxi_20170414.RData") # years 1973-2012
load("~/Biostuff/SEMFR_DATA/verena_files/SOS/Death Cancer Patient Birthregistry 150202/mfr150202/postalcodes2013_verenasIDs_maxi_20170414.RData") # additional year 2013
dat = rbind(dat,postalcodes2013); rm(postalcodes2013)
colnames(dat); nrow(dat)

#save(list = "dat",file="~/Biostuff/SEMFR_HEIGHT-GESTAGE/starting_dataset_1973-2013.RData")
#load("~/Biostuff/SEMFR_HEIGHT-GESTAGE/starting_dataset_1973-2013.RData")
#rm(list = ls()[-grep("^dat$",ls())])

orig = dat # temporarily save the raw (original) version for twin QC and identification procedures


# load helper-functions
source("~/Dropbox/GIT/SEMFR_CLEAN_THE_DATA/1_cleaning_modules.R")
source("~/Dropbox/GIT/SEMFR_CLEAN_THE_DATA/2_renumber_parity_to_parityF.R")


###########################
########################### standard cleaning
###########################

### 0)  initiate the year matrix (for plotting exclusions etc)
year_matrix = NULL; generate_year_counts(dat, stage="initial", show=T)

### 1) MANDATORY CLANING BLOCK
dat = fun_momID(dat)
dat = fun_momkidID(dat)
dat = fun_origin(dat=dat, orig=orig, countryBlocks="nordic",
                 parentVariabl=c("MFODLAND","MNAT"),strict=T) # should be run before the BMI filter
#save(list = c("dat","year_matrix"),
#     file="~/Biostuff/SEMFR_HEIGHT-GESTAGE/deleteme_20170414_3rd-stage.RData")
#load("~/Biostuff/SEMFR_HEIGHT-GESTAGE/deleteme_20170414_3rd-stage.RData")
#fun_visualize_exclusions_by_year(year_matrix)



### 2) recalculate the parity to PARITET_F
new_parity = recalculateParity(dat,thr_d_low = 31,thr_d_upp = 200) # takes about 5 min
dat$PARITET_F = new_parity; rm(new_parity) # note, NA's get introduced, but no rows removed!
#save(list = c("dat","year_matrix"),
#     file="~/Biostuff/SEMFR_HEIGHT-GESTAGE/deleteme_20170414_4th-stage.RData")
#load("~/Biostuff/SEMFR_HEIGHT-GESTAGE/deleteme_20170414_4th-stage.RData")


### 3) HEIGHT CLEANING BLOCK
dat = fun_mHghQC(dat=dat,low=140,upp=200,discordantClean = T,transferHeight = T) # this must be run first
dat = fun_mWghQC(dat,duringPregTransf = TRUE) # recovering some weight values (does not reduce the file size)
dat = fun_mBmiQC(dat) # bivariate cleaning (weight vs height)
#save(list = c("dat","year_matrix"),
#     file="~/Biostuff/SEMFR_HEIGHT-GESTAGE/deleteme_20170414_5th-stage.RData")
#load("~/Biostuff/SEMFR_HEIGHT-GESTAGE/deleteme_20170414_5th-stage.RData")
fun_visualize_exclusions_by_year(year_matrix)

### 4) DEADBORNS
dat = fun_deadborn(dat) 
#save(list = c("dat","year_matrix"),
#     file="~/Biostuff/SEMFR_HEIGHT-GESTAGE/deleteme_20170414_6th-stage.RData")
#load("~/Biostuff/SEMFR_HEIGHT-GESTAGE/deleteme_20170414_6th-stage.RData")
fun_visualize_exclusions_by_year(year_matrix)

### 5) GESTATIONAL AGE THINGS
dat = fun_GAmiss(dat)
#save(list = c("dat","year_matrix"),
#     file="~/Biostuff/SEMFR_HEIGHT-GESTAGE/deleteme_20170414_7th-stage.RData")
#load("~/Biostuff/SEMFR_HEIGHT-GESTAGE/deleteme_20170414_7th-stage.RData")
fun_visualize_exclusions_by_year(year_matrix)

dat = fun_GAdating(dat,c(1)) #5:7
#save(list = c("dat","year_matrix"),
#     file="~/Biostuff/SEMFR_HEIGHT-GESTAGE/deleteme_20170414_8th-stage.RData")
#load("~/Biostuff/SEMFR_HEIGHT-GESTAGE/deleteme_20170414_8th-stage.RData")
fun_visualize_exclusions_by_year(year_matrix)

dat = fun_spont1990(dat)
dat = fun_noPROM(dat) # could be trimmed, esp icd9 658C and icd10 O756 (too many excluded)
dat = fun_currentCS(dat)
dat = fun_previousCS(dat)
#save(list = c("dat","year_matrix"),
#     file="~/Biostuff/SEMFR_HEIGHT-GESTAGE/deleteme_20170414_9th-stage.RData")
#load("~/Biostuff/SEMFR_HEIGHT-GESTAGE/deleteme_20170414_9th-stage.RData")
fun_visualize_exclusions_by_year(year_matrix)


dat = fun_matPrecond(dat,c("DIABETES","HYPERTON","SLE","ULCOLIT","NJURSJUK"))
#save(list = c("dat","year_matrix"),
#     file="~/Biostuff/SEMFR_HEIGHT-GESTAGE/deleteme_20170414_10th-stage.RData")
#load("~/Biostuff/SEMFR_HEIGHT-GESTAGE/deleteme_20170414_10th-stage.RData")
fun_visualize_exclusions_by_year(year_matrix)


dat = fun_ICDcodes(dat,c(T,T,T,T,T,T,T, T,T,T)) # uses icd9-10, thus eliminates years before 1987
#save(list = c("dat","year_matrix"),
#     file="~/Biostuff/SEMFR_HEIGHT-GESTAGE/deleteme_20170414_11th-stage.RData")
#load("~/Biostuff/SEMFR_HEIGHT-GESTAGE/deleteme_20170414_11th-stage.RData")
fun_visualize_exclusions_by_year(year_matrix)


### preview the PTD rate per year
#df = group_by(dat,AR) %>% summarise(n=n(),PTD = mean(GRDBS<257)) %>% ungroup()
#barplot(height = df$PTD,names.arg = df$AR)


###########################
########################### BIRTHWEIGH CLEANING (masking, not deleting)
###########################


fun_hexbplot = function(X,Y,LOG,XLAB,YLAB) {
        h=hexbin(Y~X); x=h@xcm; y=h@ycm; s=h@count
        plot(x,y,pch=1,cex=log(s,LOG)+0.2,xlab=XLAB,ylab=YLAB,
             xlim=range(X,na.rm=T),ylim=range(Y,na.rm=T))
}
fun_hexbplot(dat$GRDBS,dat$BVIKTBS,200,"gestAge","birthweight")
ix1 = which( dat$BVIKTBS>(dat$GRDBS*20-2200) & dat$BVIKTBS>(dat$GRDBS*70-12000) )
ix2 = which( dat$BVIKTBS<(dat$GRDBS*20-3900) | dat$BVIKTBS<(dat$GRDBS*10-1400) )
ix3 = which( dat$BVIKTBS>6200 )
bad_ix = unique(c(ix1,ix2,ix3)); rm(ix1,ix2,ix3)
points(dat$GRDBS[bad_ix],dat$BVIKTBS[bad_ix],pch=19,col="red",cex=0.5)
# abline(-12000,70,col="red")
# abline(-2200,20,col="red")
# abline(-3900,20,col="blue")
# abline(-1400,10,col="blue")
# abline(h=6200,col="green")
# mask suspicious values:
dat$BVIKTBS[bad_ix] = NA
rm(bad_ix)
#save(list = c("dat","year_matrix"),
#     file="~/Biostuff/SEMFR_HEIGHT-GESTAGE/deleteme_20170414_12th-stage.RData")
#load("~/Biostuff/SEMFR_HEIGHT-GESTAGE/deleteme_20170414_12th-stage.RData")


###########################
########################### reduce the file size
###########################

use_cols = c("lpnr_BARN","lpnr_mor","BFODDAT","AR","BORDF2","BORDNRF2",
             "MVIKT","MLANGD","PARITET","PARITET_F","MALDER","KON","GRDBS","GRMETOD","BVIKTBS")
dat = dat[,use_cols]; rm(use_cols)
# backup
dat_clean = dat # for data merge/recovery at later stages


###########################
###########################  detection of twins
###########################
library(dplyr)

dat = dat[which((dat$BORDF2==2)&(dat$BORDNRF2 %in% c("12","22"))),] # explicitly declarated twins
dat$midar = paste(dat$AR,dat$lpnr_mor,sep="_"); length(unique(dat$midar)) # year-motherID code

# explanation for "midar" necessity:
table(table(dat$lpnr_mor))  # -> bad choice
# since some mothers have two twin pairs, "midar" variable better describes unique twin pregnancy
table(table(dat$midar))     # -> good

# make sure that twins are really twins by crosref'ing with PARITY_F
ee = group_by(dat,midar) %>% summarise(n=n(),nP=length(unique(PARITET_F))) %>% ungroup()
table(ee$nP) # there should be no such case where one "midar" has two or more parity_F values

# ensure that variables are consistent for the same mother at the same pregnancy...
# also ensure that each pregnancy has non-missing values
df = data.frame(group_by(dat,midar) %>% 
                        summarise(n=n(),minH=min(MLANGD,na.rm=T),maxH=max(MLANGD,na.rm=T),meanH=mean(MLANGD,na.rm=T),
                                  minGA=min(GRDBS,na.rm=T),maxGA=max(GRDBS,na.rm=T),meanGA=mean(GRDBS,na.rm=T),
                                  minBW=min(BVIKTBS),maxBW=max(BVIKTBS)) %>% # important NOT to use "na.rm=T" here
                        ungroup())
df$difH = df$maxH - df$minH  # max dif between two height of same mom at same pregnancy
df$difGA = df$maxGA - df$minGA  # max dif between two gestAge values of same mom at same pregnancy
df$difBW = df$maxBW - df$minBW  # max dif between two birthweight values of same mom at same pregnancy

head(df)
table(df$n)
table(df$difH)
table(df$difGA)

# restrict the object to contain only FULL data for both twins (note: this is conservative!)
df = df[which(df$n==2),] # data from both twins must be present (e.g., if=1 -> stilbirth, if=4 -> two twin pairs)
df = df[which(!is.na(df$difBW)),] # no missing values
df = df[which(df$difH < 5),] # maternal height should not differ too much at the same pregnancy
df = df[which(df$difGA < 2),] # gestAge should not differ by more than 1d at the same pregnancy

# prune the dataset and assign new non-missing values 
df$meanH = as.integer(round(df$meanH,0))
df$meanGA = as.integer(round(df$meanGA,0))
df = df[,c("midar","meanH","meanGA")]
colnames(df) = c("midar","MLANGD","GRDBS")
head(df)
dat = dat[,-grep("^MLANGD$|^GRDBS$",colnames(dat))]
dat = merge(dat,df,by="midar",all=F) # it is ok to merge, because these are only two twins per mom
head(dat)


###### fix  erroneously assigned order-of-birth
df = data.frame(group_by(dat,midar) %>% 
                        summarise(n=n(),
                                  n12=sum(BORDNRF2=="12"),
                                  n22=sum(BORDNRF2=="22")) %>% 
                        ungroup())
# preview:
#df[which(df$n12==0),]
#df[which(df$n22==0),]
# at a random possition assign the oposite (remaining order)
to_fix_ids = df$midar[which( (df$n12==0)|(df$n22==0) )]
for (midar in to_fix_ids) {
dat[which(dat$midar == midar),"BORDNRF2"] = sample(c(12,22),2,replace = F) # order does not realy matter
}

#####################
#### convert it to wide format
#####################

library(tidyr)

####  convert to wide format the BVIKTBS values
tmp = dat[,c("midar","BORDNRF2","BVIKTBS")]
bvikt = spread(data=tmp,key=BORDNRF2,value=BVIKTBS)
colnames(bvikt) = c("midar","BVIKTBS_12","BVIKTBS_22")
head(bvikt); rm(tmp)

####  convert to wide format the PARITY values
tmp = dat[,c("midar","BORDNRF2","PARITET_F")]   # note! not PARITET, but PARITET_F !
parit = spread(data=tmp,key=BORDNRF2,value=PARITET_F)
parit = parit[,c("midar","12")]  # since both twins get the same value anyway
colnames(parit) = c("midar","PARITET_F") 
head(parit); rm(tmp)

####  convert to wide format the MLANGD values
tmp = dat[,c("midar","BORDNRF2","MLANGD")]
mlang = spread(data=tmp,key=BORDNRF2,value=MLANGD)
mlang = mlang[,c("midar","12")]  #  both values (for both twins) are the same
colnames(mlang) = c("midar","MLANGD")
head(mlang); rm(tmp)


####  convert to wide format the GRDBS values
tmp = dat[,c("midar","BORDNRF2","GRDBS")]
grdbs = spread(data=tmp,key=BORDNRF2,value=GRDBS)
grdbs = grdbs[,c("midar","12")]  #  both values (for both twins) are the same
colnames(grdbs) = c("midar","GRDBS")
head(grdbs); rm(tmp)

####  convert to wide format the lpnr_BARN values
tmp = dat[,c("midar","BORDNRF2","lpnr_BARN")]
barnid = spread(data=tmp,key=BORDNRF2,value=lpnr_BARN)
colnames(barnid) = c("midar","lpnr_BARN_12","lpnr_BARN_22")
head(barnid); rm(tmp)

mrg3 = merge(parit,mlang,by="midar",all=T)
mrg2 = merge(mrg3,bvikt,by="midar",all=T)
mrg1 = merge(mrg2,grdbs,by="midar",all=T)
mrg = merge(mrg1,barnid,by="midar",all=T)
head(mrg)
table(mrg$PARITET_F)



###################
################### identify birthweight-discordant twins
###################


set.seed(1)
ix = sample(c(0,1),nrow(mrg),replace=T) # randomize which twin is 1's
mrg$bw1 = mrg$bw2 = NA
mrg$bw1[which(ix==0)] = mrg$BVIKTBS_12[which(ix==0)]
mrg$bw1[which(ix==1)] = mrg$BVIKTBS_22[which(ix==1)]
mrg$bw2[which(ix==0)] = mrg$BVIKTBS_22[which(ix==0)]
mrg$bw2[which(ix==1)] = mrg$BVIKTBS_12[which(ix==1)]
#mrg = mrg[,c("midar","bw1","bw2","ga1","ga2")]
#mrg = mrg[,-grep("BVIKTBS",colnames(mrg))]
mrg$bwdif = mrg$bw1 - mrg$bw2
hist(mrg$bwdif,xlab="twin birthweight difference (g)",breaks=100,col="grey")
head(mrg)

# first, need to take into account gestAge of the twin-pair
mrg$ga_bin = mrg$GRDBS %/% 7 * 7 + 7/2 # devide into larger bins
table(mrg$ga_bin) # doublecheck
sds = NULL  # what is the SD of BWdif for each gestAge bin
nob = NULL # number of observations in each gestAge bin
gas = sort(unique(mrg$ga_bin))
for (ga in gas) {
        if (sum(mrg$ga_bin==ga)>=10) {
                difs = mrg$bwdif[which(mrg$ga_bin == ga)]
                sds = c(sds,sd(difs,na.rm=T))
                nob = c(nob,sum(!is.na(difs)))
        } else {
                sds = c(sds,NA)
                nob = c(nob,sum(mrg$ga_bin==ga))
        }
}
plot(sds ~ gas)
#mod = loess(sds ~ gas,span = 1,weights = log(nob)) # INTRODUCES NAs at extreme ends of GA
mod = lm(sds ~ gas) # option to use poly(x,2)
points(predict(mod,data.frame(gas=gas)) ~ gas,type="l",col="grey",lwd=0.5)
ga_vals = sort(unique(mrg$GRDBS))
devs = data.frame(ga=ga_vals,BWdev1sd=predict(mod,data.frame(gas=ga_vals))) # deviations
mrg = merge(mrg,devs,by.x="GRDBS",by.y="ga",all.x=T)
head(mrg)

SD_THR = 2
mrg$excl = as.numeric(abs(mrg$bwdif) > (mrg$BWdev1sd * SD_THR))
ix = which(mrg$excl==0)
plot(mrg$bwdif[ix] ~ mrg$GRDBS[ix],xlim=range(mrg$GRDBS),ylim=range(mrg$bwdif), # mrg$excl+1
     xlab="gestational age at birth (days)",ylab="birthweight difference (g)",col="darkgrey")
points(y=mrg$bwdif[-ix],x= mrg$GRDBS[-ix],col="indianred1",pch=19,cex=0.7) #brown1  lightcoral
legend(165,-700,legend = c("exclusion","analysis"),title="selected for:",pch=c(19,1),col=c("indianred1","darkgrey"),bty="n")
table(mrg$excl,useNA="a")
mean(mrg$excl)
head(mrg)



twins = mrg[,c("midar","GRDBS","PARITET_F","MLANGD","lpnr_BARN_12","lpnr_BARN_22",
               "BVIKTBS_12","BVIKTBS_22","bwdif","BWdev1sd","excl")]

#save(list = c("twins","mrg"),
#     file="~/Biostuff/SEMFR_HEIGHT-GESTAGE/deleteme_20170414_twins.RData")
#load("~/Biostuff/SEMFR_HEIGHT-GESTAGE/deleteme_20170414_twins.RData")




###########################
###########################  SINGLETONS
###########################

# broadest way of detecting twin(+) pregnacies - based on parity_f

load("~/Biostuff/SEMFR_HEIGHT-GESTAGE/deleteme_20170414_4th-stage.RData") # (right after PARITET_F inference)
dat = dat[,c("lpnr_mor","lpnr_BARN","PARITET_F")]
dat = dat[which(!is.na(dat$PARITET_F)),] # since duplicated(c(NA,NA)) will give TRUE
dat = as.data.frame(group_by(dat, lpnr_mor) %>% arrange(lpnr_mor,PARITET_F) %>%
        mutate(pfdif = PARITET_F - lag(PARITET_F)))
ix1 = which(dat$pfdif==0) # this is the "second" twin
ix = unique(c(ix1,ix1-1)) # these are all twins and triplets etc
twin_child_ids = dat$lpnr_BARN[ix] # convert indexes to IDs
rm(dat,year_matrix,ix1,ix)


# load cleaned data
load("~/Biostuff/SEMFR_HEIGHT-GESTAGE/deleteme_20170414_12th-stage.RData")
# reduce the file size
use_cols = c("lpnr_BARN","lpnr_mor","BFODDAT","AR","BORDF2","BORDNRF2",
             "MVIKT","MLANGD","PARITET","PARITET_F","MALDER","KON","GRDBS","GRMETOD","BVIKTBS")
dat = dat[,use_cols]; rm(use_cols)

# remove recently identified twins
dat = dat[which(!dat$lpnr_BARN %in% twin_child_ids),]
# remove twins based on the declared twin-status
dat = dat[which(dat$BORDF2==1),]

# rename
singl = dat; rm(dat)

# extra-extra: remove twins mentioned in the twin register
twi = read.table("~/Biostuff/SEMFR_DATA/verena_files/SOS/Twin Registry 140918/din_fil_/din_fil_.dat",h=T,sep=",",stringsAsFactors = F)
# note- a possible non-overlap conflict with the merged dataset where a part is used from postal codes file
sum(singl$lpnr_BARN %in% twi$lpnr) # 2-6 twins were not detected previously
singl = singl[which(!singl$lpnr_BARN %in% twi$lpnr),]

# save for future use
save(list = c("singl"),file="~/Biostuff/SEMFR_HEIGHT-GESTAGE/deleteme_20170414_singl.RData")
#load("~/Biostuff/SEMFR_HEIGHT-GESTAGE/deleteme_20170414_singl.RData")

#######  cleaning is finished


#######
#######  PREPARE DATA FOR ANALYSES
#######

load("~/Biostuff/SEMFR_HEIGHT-GESTAGE/deleteme_20170414_singl.RData")
load("~/Biostuff/SEMFR_HEIGHT-GESTAGE/deleteme_20170414_twins.RData")
load("~/Biostuff/SEMFR_HEIGHT-GESTAGE/deleteme_20170414_12th-stage.RData") # full data (to ensure the same variables are used)


twins$KON=NA
sup = dat[,c("lpnr_BARN","lpnr_mor","MFODDAT","AR","MFODLAND","MNAT","MALDER")] # KON
twins = merge(twins,sup,by.x="lpnr_BARN_12",by.y="lpnr_BARN",all.x=T)
colnames(twins)[which(colnames(twins)=="lpnr_BARN_12")] = "lpnr_BARN"
twins$BVIKTBS = round((twins$BVIKTBS_12 + twins$BVIKTBS_22) / 2,0)
twins$status = "twins"

singl$excl = NA
sup = dat[,c("lpnr_BARN","MFODDAT","MFODLAND","MNAT")]
singl = merge(singl,sup,by="lpnr_BARN",all.x=T)
singl$status = "singl"

col_names = c("lpnr_BARN","lpnr_mor","status","excl",
"AR","GRDBS","KON","BVIKTBS","MLANGD","MALDER","MFODDAT","MFODLAND","MNAT","PARITET_F")

twins = twins[,col_names]
singl = singl[,col_names]

dat = rbind(twins,singl)
dat = dat[sample(nrow(dat),nrow(dat),replace = F),]

## make sure that no pregnancies with misisng data remains
bad_rows = which((is.na(dat$BVIKTBS))&(dat$status=="singl"))
dat = dat[-bad_rows,]; rm(bad_rows)
nrow(dat)

save(list = c("dat"),file="~/Biostuff/SEMFR_HEIGHT-GESTAGE/working_dataset_TWISIN_20170613_n1159478.RData")
#load("~/Biostuff/SEMFR_HEIGHT-GESTAGE/working_dataset_TWISIN_20170613_n1159478.RData")

head(dat)

# doublecheck (should be 0)
sum(duplicated(dat$lpnr_mor[which(dat$PARITET_F==1)]))



