
# project: causality in maternal height <-> gest.Age association
# alternative title: uterine overdistension determines gest.Age
# prepared by Jonas B. 2017 July

# (this is a cleaned version for Ge to rerun this on USA data)

library(dplyr)

# load the data
load("~/Biostuff/SEMFR_HEIGHT-GESTAGE/working_dataset_TWISIN_20170613_n1159478.RData")

## STRUCTURE:   data frame with child-centered structure (one row per child);
##              only one twin from a pair is kept, the other is removed;
##              compulsory columns in the loaded dataframe (called "dat") are:
# PARITET_F     mother's parity (1 = this is her first pregnancy, 2= her second pregnancy, etc)
# MALDER        mother's age at pregnanc, in years
# MLANGD        motehr's height at pregnancy in cm
# GRDBS         gestational age in days (ultrasound-based)
# BVIKTBS       birtweight in grams
# KON           child's sex (1=boy, 2=girl)
# status        "twins" or "singl"  (one of these character strings as values)

# FILTERS and EXCLUSIONS (not reflected in this script, but applied previously) :
# .. only spontaneous deliveries (i.e., not induced)
# .. only ultrasound-evaluated gestational age
# .. no PROM deliveries (delivery starts with labor, not rupture of membranes)
# .. no maternal medical preconditions (diabetes etc)
# .. no pregnancy complications (placenta previa etc)
# .. fetal twins with very discordant weight are both excluded


# load the data
#dat = read.table("~/Biostuff/USA_BIRTH_REGISTER/Nat2014PublicUS.c20150514.r20151022.tab",h=T,stringsAsFactors = F,sep="\t")
#save(list=c("dat"),file = "~/Biostuff/USA_BIRTH_REGISTER/loaded_data.RData")
load("~/Biostuff/USA_BIRTH_REGISTER/loaded_data.RData")


# MRACE6 (1=white)
# MRACE15  (01=white)
# FRACE31 (01 = white) for father
# FBRACE  (1=white) for father
# PRIORLIVE
# PRIORDEAD
# PRIORTERM
# DPLURAL (1=single, 2=twin)
# MAGER (mother's age)
# OEGest_Comb  (gest age in weeks)
# DBWT  (birthweight in gramms)

####  Leaev only White Americans
#table(dat$MRACE6,dat$MRACE15,useNA = "a")
#table(dat$FRACE31,dat$FBRACE,useNA = "a")
dat = dat[which((dat$MRACE6==1)&(dat$MRACE15==1)),]
dat = dat[which((dat$FRACE31==1)&(dat$FBRACE==1)),]


#### LEAVE ONLY PARITY 0
#table(dat$PRIORTERM,useNA = "a")
#table(dat$PRIORLIVE,dat$PRIORDEAD,useNA = "a")
dat = dat[which((dat$PRIORLIVE==0)&(dat$PRIORDEAD==0)&(dat$PRIORTERM==0)),]
dim(dat)

### singletons and twins
table(dat$DPLURAL,useNA="a")
dat = dat[which(dat$DPLURAL %in% c(1,2)),]

### Leave only one child per pregnancy
#table(dat$SETORDER_R,dat$DPLURAL,useNA="a")
bad_rows = which( (dat$DPLURAL==2)&(dat$SETORDER_R==2))
dat = dat[-bad_rows,]; rm(bad_rows)

# Mother's age
table(dat$MAGER,useNA = "a")
dat = dat[which((dat$MAGER>=18)&(dat$MAGER<=45)),]
dim(dat)

## Gest Age and Birthweight
#hist(dat$OEGest_Comb,breaks=100,col="Grey")
#hist(dat$DBWT,breaks=100,col="Grey")
#plot(dat$OEGest_Comb,dat$DBWT)
dat$GRDBS = dat$OEGest_Comb * 7


# maternal height
hist(dat$M_Ht_In)
dat$MLANGD = dat$M_Ht_In * 2.54
hist(dat$MLANGD)


dat$status = "singl"
dat$status[which(dat$DPLURAL==2)] = "twins"
head(dat)

### no missing data
dat = dat[which(!is.na(dat$MLANGD)),]
dat = dat[which(!is.na(dat$GRDBS)),]

# data cleaning
hist(dat$MLANGD)
dat = dat[which(dat$MLANGD>=140),]

save(list=c("dat"),file = "~/Biostuff/USA_BIRTH_REGISTER/loaded_data_digested.RData")
load("~/Biostuff/USA_BIRTH_REGISTER/loaded_data_digested.RData")


###########
#dat = dat[,c("GRDBS","DBWT","DPLURAL","MLANGD")]




#########
#########  PLOT 1:  TWINS vs SINGLETONS
#########

######  polynomial linear regression model
mod = lm(GRDBS ~ poly(MLANGD,2) * status, data = dat)
mod = lm(GRDBS ~ MLANGD * status, data = dat)
summary(mod)


fun= function() {
        df = expand.grid(status=c("twins","singl"),MLANGD=140:200)
        df = df[order(df$status,df$MLANGD),]
        prd = predict(mod,df,interval = "confidence",level = 0.95) # predicted values and confidence intervals
        df$nw = prd[,1] # nw = new values
        
        tst = df$status=="singl"
        plot(df$MLANGD[tst],df$nw[tst],col="cornflowerblue",type="l",
             xlim=range(dat$MLANGD),
             ylim=c(min(df$nw)*0.98,max(df$nw)*1.02),
             ylab="child's gestational age at birth (days)",xlab="maternal height (cm)")
        points(df$MLANGD[!tst],df$nw[!tst],col="orange",type="l")
        grid()
        #clean the area for legends (white background)
        polygon(x=c(154,195,195,154,154),y=c(221,221,239,239,221),border = NA,col="white")
        
        ##### confidence intervals for the model curves
        # for singletons
        xs = c(df$MLANGD[tst],rev(df$MLANGD[tst]),df$MLANGD[tst][1])
        ys = c(prd[tst,2],rev(prd[tst,3]),prd[tst,2][1])
        clrcm = col2rgb("cornflowerblue") / 255 # color components
        ciclr_sin = rgb(clrcm[1,1],clrcm[2,1],clrcm[3,1],0.4) # color for conf. intervals
        polygon(xs,ys,col=ciclr_sin,border = NA)
        
        # for twins
        xs = c(df$MLANGD[!tst],rev(df$MLANGD[!tst]),df$MLANGD[!tst][1])
        ys = c(prd[!tst,2],rev(prd[!tst,3]),prd[!tst,2][1])
        clrcm = col2rgb("orange") / 255 # color components
        ciclr_twi = rgb(clrcm[1,1],clrcm[2,1],clrcm[3,1],0.4) # color for conf. intervals
        polygon(xs,ys,col=ciclr_twi,border = NA)
        
        
        ###### empirical plot (means in bins)
        dat$matHgh_cat = dat$MLANGD %/% 5 * 5 + 5/2
        df = group_by(dat,matHgh_cat,status) %>% 
                summarise(n=n(),mnGA=mean(GRDBS),seGA=sd(GRDBS)/sqrt(n)) %>% filter(n>10) %>% ungroup()
        df = as.data.frame(df)
        df$status[which(df$status=="twins")] = "orange"
        df$status[which(df$status=="singl")] = "cornflowerblue"
        # data points
        points(y = df$mnGA, x = df$matHgh_cat,pch=19,col=df$status, cex=log(df$n,50))
        points(y = df$mnGA, x = df$matHgh_cat,cex=log(df$n,50)) # "coat"
        # 95% Conf.intervals
        df$upp = df$mnGA + 1.96*df$seGA
        df$low = df$mnGA - 1.96*df$seGA
        segments(x0 = df$matHgh_cat,x1 = df$matHgh_cat,y0 = df$low,y1=df$upp,lwd=1.3,col="white")
        segments(x0 = df$matHgh_cat,x1 = df$matHgh_cat,y0 = df$low,y1=df$upp,lwd=1,col="darkgrey")
        
        
        # LEGEND 1
        x = 155; dx = 5.0; xs = c(x,x+dx,x+dx,x,x)
        y = 238; dy = 2; ys = c(y-dy,y-dy,y+dy,y+dy,y-dy)
        polygon(xs,ys,col=ciclr_sin,border = NA)
        segments(x,y,x+dx,y,col="cornflowerblue") # SGA
        polygon(xs,ys-5,col=ciclr_twi,border = NA)
        segments(x,y-5,x+dx,y-5,col="orange")
        text(x-1,y+4,"model, 95% CI",pos = 4,cex = 0.7)
        text(x+4,y,"singletons",pos = 4,cex = 0.7)
        text(x+4,y-5,"twins",pos = 4,cex = 0.7)
        
        # LEGEND 2
        x = 170; y = 239; dy = 0.5
        text(x-2,y+4,"mean, 95% CI",pos = 4,cex = 0.7)
        points(c(x,x-0.2),c(y,y-5),pch=19,cex=1.2,col=c("cornflowerblue","orange"))
        segments(x,y-2,x,y+2,col="darkgrey") # SINGLETONS
        segments(x-0.2,y-5-3,x-0.2,y-5+2,col="darkgrey") # TWINS
        text(x,y,"singletons",pos = 4,cex = 0.7)
        text(x,y-5,"twins",pos = 4,cex = 0.7)
        
        # LEGEND 3
        x = 183; y = 239
        text(x-2,y,"observations (n)",pos = 4,cex = 0.7)
        szs = c(1e2,1e3,1e4)
        legend(x,y,legend = szs,pch=1,bty = "n",pt.cex = log(szs,50))
}
fun()






#########
#########  PLOT 2:  LGA vs SGA
#########

sin = dat[which(dat$status=="singl"),]  # LEAVE ONLY SINGLETONS
sin = sin[which((!is.na(sin$BVIKTBS))&(!is.na(sin$GRDBS))),]  

# function estimating gest.Age-adjusted birtweight Z-score  (Marsal et al)
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

# estimate BWzsc (pretending that sex is unknown for all children)
sin$bw_zsc_Marsal=sapply(1:nrow(sin), function(x) marsal(sin$GRDBS[x],sin$BVIKTBS[x],sex=sin$KON[x]))
head(sin); dim(sin)
hist(sin$bw_zsc_Marsal)

# define two fetal growth strata
sin$cat = NA  
sin$cat[which(sin$bw_zsc_Marsal<(-1.5))] = "SGA"  
sin$cat[which(sin$bw_zsc_Marsal>1.5)] = "LGA"  
sin = sin[which(!is.na(sin$cat)),] # remove the "apropriate for gestational age"
table(sin$cat)

######  polynomial linear regression model
mod = lm(GRDBS ~ poly(MLANGD,2) * cat, data = sin)

fun = function() {
        
        df = expand.grid(cat=c("LGA","SGA"),MLANGD=140:196)
        df = df[order(df$cat,df$MLANGD),]
        prd = predict(mod,df,interval = "confidence",level = 0.95) # 
        df$nw = prd[,1] # new value (predicted one)
        tst = df$cat=="SGA"
        plot(df$MLANGD[tst],df$nw[tst],col="cornflowerblue",type="l",ylim=c(min(df$nw)*0.98,max(df$nw)*1.01),
             ylab="child's gestational age at birth (days)",xlab="maternal height (cm)")
        grid()
        points(df$MLANGD[!tst],df$nw[!tst],col="orange",type="l")
        #clean the area for legends (white background)
        polygon(x = c(151,191,191,151,151),y=c(265,265,272,272,265),col = "white",border = NA) 
        
        ##### confidence intervals for the model curves
        # for SGA
        xs = c(df$MLANGD[tst],rev(df$MLANGD[tst]),df$MLANGD[tst][1])
        ys = c(prd[tst,2],rev(prd[tst,3]),prd[tst,2][1])
        clrcm = col2rgb("cornflowerblue") / 255 # color components
        ciclr_sga = rgb(clrcm[1,1],clrcm[2,1],clrcm[3,1],0.4) # color for conf. intervals
        polygon(xs,ys,col=ciclr_sga,border = NA)
        # for LGA
        xs = c(df$MLANGD[!tst],rev(df$MLANGD[!tst]),df$MLANGD[!tst][1])
        ys = c(prd[!tst,2],rev(prd[!tst,3]),prd[!tst,2][1])
        clrcm = col2rgb("orange") / 255 # color components
        ciclr_lga = rgb(clrcm[1,1],clrcm[2,1],clrcm[3,1],0.4) # color for conf. intervals
        polygon(xs,ys,col=ciclr_lga,border = NA)
        
        # LEGEND 1
        x = 152; dx = 5.0; xs = c(x,x+dx,x+dx,x,x)
        y = 269; dy = 0.5; ys = c(y-dy,y-dy,y+dy,y+dy,y-dy)
        polygon(xs,ys,col=ciclr_sga,border = NA)
        segments(x,y,x+dx,y,col="cornflowerblue") # SGA
        polygon(xs,ys-1.3,col=ciclr_lga,border = NA)
        segments(x,y-1.3,x+dx,y-1.3,col="orange")
        text(x-2,y+1.5,"model, 95% CI",pos = 4,cex = 0.7)
        text(x+5,y,"SGA",pos = 4,cex = 0.7)
        text(x+5,y-1.3,"LGA",pos = 4,cex = 0.7)
        
        ###### data for empirical means
        sin$matHgh_cat = sin$MLANGD %/% 5 * 5 + 5/2 # more data = smaller bins
        df_sing = group_by(sin,cat,matHgh_cat) %>% 
                summarise(n=n(),mnGA=mean(GRDBS),seGA=sd(GRDBS)/sqrt(n)) %>% filter(n>10) %>% ungroup() 
        df_sing = as.data.frame(df_sing)
        df_sing$clr = NA
        df_sing$clr = ifelse(df_sing$cat=="LGA","orange","cornflowerblue")
        df_sing = df_sing[order(df_sing$cat,decreasing = T),]
        # small shift to sseparate LGA and SGA
        vls = df_sing$matHgh_cat[which(df_sing$cat=="LGA")]
        df_sing$matHgh_cat[which(df_sing$cat=="LGA")] = vls - 0.2
        
        # data points
        points(y = df_sing$mnGA, x = df_sing$matHgh_cat,pch=19,col=df_sing$clr, cex=log(df_sing$n,50))
        points(y = df_sing$mnGA, x = df_sing$matHgh_cat,cex=log(df_sing$n,50)) # "coat"
        # 95% Conf.intervals
        df_sing$upp = df_sing$mnGA + 1.96*df_sing$seGA
        df_sing$low = df_sing$mnGA - 1.96*df_sing$seGA
        segments(x0 = df_sing$matHgh_cat,x1 = df_sing$matHgh_cat,
                 y0 = df_sing$low,y1=df_sing$upp,lwd=1.5,col="white") # white background
        segments(x0 = df_sing$matHgh_cat,x1 = df_sing$matHgh_cat,
                 y0 = df_sing$low,y1=df_sing$upp,lwd=1,col="darkgrey") # line itself
        
        # LEGEND 2
        x = 168; y = 269; dy = 0.5
        text(x-5,y+1.5,"mean, 95% CI",pos = 4,cex = 0.7)
        points(c(x,x-0.2),c(y,y-1.3),pch=19,cex=1.2,col=c("cornflowerblue","orange"))
        segments(x,y-0.5,x,y+0.8,col="darkgrey") # SGA
        segments(x-0.2,y-1.3-0.8,x-0.2,y-1.3+0.5,col="darkgrey") # LGA
        text(x,y,"SGA",pos = 4,cex = 0.7)
        text(x,y-1.3,"LGA",pos = 4,cex = 0.7)
        
        # LEGEND 3
        x = 179; y = 271
        text(x-3,y,"observations (n)",pos = 4,cex = 0.7)
        szs = c(1e2,1e3,1e4)
        legend(x,y,legend = szs,pch=1,bty = "n",pt.cex = log(szs,50))
        
}
fun()




