
# project: causality in maternal height - gest.Age association
# alternative title: uterine overdistension determines gest.Age
# prepared by Jonas B. 2017 July

library(dplyr)

# load the data
load("~/Biostuff/SEMFR_HEIGHT-GESTAGE/working_dataset_TWISIN_20170613_n1159478.RData")

# filter data
dat = dat[which(dat$PARITET_F==1),]  # THIS MUST BE THE FIRST PREGNANCY OF THE MOTHER
dat = dat[which((dat$excl==0)|(is.na(dat$excl))),] #
dat = dat[which(dat$MALDER %in% 18:45),]  # MOTHERS AGE
nrow(dat)

table(dat$status,useNA = "a")
hist(dat$MALDER,breaks=100,col="grey")
range(dat$MALDER)
table(dat$MALDER)
table(dat$MALDER,dat$status)


#########  TWINS vs SINGLETONS
#

##### model with separate groups
mod0 = lm(GRDBS ~ MLANGD + KON + MALDER,data = dat[which(dat$status=="singl"),])
coef(summary(mod0)) # summary(mod0) # 0.1197  beta

mod1 = lm(GRDBS ~ MLANGD + MALDER + factor(KON,levels = c("BG","BB","GG")), data = dat[which(dat$status=="twins"),])
coef(summary(mod1)) # summary(mod1) # 0.2588 beta

#


##### linear interaction model
mod = lm(GRDBS ~ MLANGD * status + MALDER, data=dat)
coef(summary(mod)) # summary(mod)
cfs = coef(summary(mod)); cfs
write.table(cfs,"~/Downloads/TWI-SIN_linear_interaction_coefs.txt",quote=F,sep="\t")


######  quadratic-model-based plot
mod = lm(GRDBS ~ poly(MLANGD,2) * status, data = dat)
summary(mod)
cfs = coef(summary(mod)); cfs
write.table(cfs,"~/Downloads/TWI-SIN_quadratic_interaction_coefs.txt",quote=F,sep="\t")

########################################################################
# Figure 2: Twins-vs-Singletons in SIMPLE LINEAR INTERACTION MODEL
########################################################################

#####   WITH LINEAR MODEL
mod = lm(GRDBS ~ MLANGD * status + MALDER, data = dat)
fun= function(x3,y3) {
df = expand.grid(status=c("twins","singl"),MLANGD=140:196,MALDER=29)
df = df[order(df$status,df$MLANGD),]
prd = predict(mod,df,interval = "confidence",level = 0.95) # predicted values and confidence intervals
df$nw = prd[,1] # nw = new values

tst = df$status=="singl"
plot(df$MLANGD[tst],df$nw[tst],col="cornflowerblue",type="l",ylim=c(min(df$nw)*0.92,max(df$nw)*1.02),
     ylab="child's gestational age at birth (days)",xlab="maternal height (cm)")
points(df$MLANGD[!tst],df$nw[!tst],col="orange",type="l")
grid()
#clean the area for legends
polygon(x=c(159,196,196,159,159),y=c(221,221,245,245,221),border = NA,col="white")

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


###### empirical plot
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



#legend(x = 160,y = 240,legend = c("singletons","twins"),col = c("cornflowerblue","orange"),lwd=2,bty = "n")
# legend
#szs = c(10,1e2,1e3,1e4)
#legend(x = 180,y = 245,legend = szs,pch=1,bty = "n",pt.cex = log(szs,50))

# LEGEND 1
x = x3[1]; dx = 5.0; xs = c(x,x+dx,x+dx,x,x)
y = y3[1]; dy = 2; ys = c(y-dy,y-dy,y+dy,y+dy,y-dy)
polygon(xs,ys,col=ciclr_sin,border = NA)
segments(x,y,x+dx,y,col="cornflowerblue") # SGA
polygon(xs,ys-5,col=ciclr_twi,border = NA)
segments(x,y-5,x+dx,y-5,col="orange")
text(x-1,y+4,"model, 95% CI",pos = 4,cex = 0.7)
text(x+4,y,"singletons",pos = 4,cex = 0.7)
text(x+4,y-5,"twins",pos = 4,cex = 0.7)

# LEGEND 2
x = x3[2]; y = y3[2]; dy = 0.5
text(x-2,y+4,"mean, 95% CI",pos = 4,cex = 0.7)
points(c(x,x-0.2),c(y,y-5),pch=19,cex=1.2,col=c("cornflowerblue","orange"))
segments(x,y-2,x,y+2,col="darkgrey") # SGA
segments(x-0.2,y-5-3,x-0.2,y-5+2,col="darkgrey") # LGA
text(x,y,"singletons",pos = 4,cex = 0.7)
text(x,y-5,"twins",pos = 4,cex = 0.7)

# LEGEND 3
x = x3[3]; y = y3[3]
text(x-1,y-3,"N (observations)",pos = 4,cex = 0.7)
szs = c(1e2,1e3,1e4)
legend(x,y-3,legend = szs,pch=1,bty = "n",pt.cex = log(szs,50))
}
pdf("~/Downloads/Twins_vs_Singletons.pdf",width = 5.2,height = 5.5)
fun(x3=c(158,174,185),y3=c(235,235,242)+3)
dev.off()



########################################################################
# Suppementary Figure: Twins-vs-Singletons in POLYNOMIAL INTERACTION MODEL
########################################################################

mod = lm(GRDBS ~ poly(MLANGD,2) * status + MALDER, data = dat)
fun= function() {
        df = expand.grid(status=c("twins","singl"),MLANGD=140:196,MALDER=29)
        df = df[order(df$status,df$MLANGD),]
        prd = predict(mod,df,interval = "confidence",level = 0.95) # predicted values and confidence intervals
        df$nw = prd[,1] # nw = new values
        
        tst = df$status=="singl"
        plot(df$MLANGD[tst],df$nw[tst],col="cornflowerblue",type="l",ylim=c(min(df$nw)*0.98,max(df$nw)*1.02),
             ylab="child's gestational age at birth (days)",xlab="maternal height (cm)")
        points(df$MLANGD[!tst],df$nw[!tst],col="orange",type="l")
        grid()
        #clean the area for legends
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
        
        
        ###### empirical plot
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
        
        
        
        #legend(x = 160,y = 240,legend = c("singletons","twins"),col = c("cornflowerblue","orange"),lwd=2,bty = "n")
        # legend
        #szs = c(10,1e2,1e3,1e4)
        #legend(x = 180,y = 245,legend = szs,pch=1,bty = "n",pt.cex = log(szs,50))
        
        # LEGEND 1
        x = 155; dx = 5.0; xs = c(x,x+dx,x+dx,x,x)
        y = 230; dy = 2; ys = c(y-dy,y-dy,y+dy,y+dy,y-dy)
        polygon(xs,ys,col=ciclr_sin,border = NA)
        segments(x,y,x+dx,y,col="cornflowerblue") # SGA
        polygon(xs,ys-5,col=ciclr_twi,border = NA)
        segments(x,y-5,x+dx,y-5,col="orange")
        text(x-1,y+4,"model, 95% CI",pos = 4,cex = 0.7)
        text(x+4,y,"singletons",pos = 4,cex = 0.7)
        text(x+4,y-5,"twins",pos = 4,cex = 0.7)
        
        # LEGEND 2
        x = 170; y = 230; dy = 0.5
        text(x-2,y+4,"mean, 95% CI",pos = 4,cex = 0.7)
        points(c(x,x-0.2),c(y,y-5),pch=19,cex=1.2,col=c("cornflowerblue","orange"))
        segments(x,y-2,x,y+2,col="darkgrey") # SGA
        segments(x-0.2,y-5-3,x-0.2,y-5+2,col="darkgrey") # LGA
        text(x,y,"singletons",pos = 4,cex = 0.7)
        text(x,y-5,"twins",pos = 4,cex = 0.7)
        
        # LEGEND 3
        x = 183; y = 235
        text(x-2,y,"observations (n)",pos = 4,cex = 0.7)
        szs = c(1e2,1e3,1e4)
        legend(x,y,legend = szs,pch=1,bty = "n",pt.cex = log(szs,50))
}
fun()









#####
#### alternative method with LGA/SGA
#####

sin = dat[which(dat$status=="singl"),]
sin = sin[which((!is.na(sin$BVIKTBS))&(!is.na(sin$GRDBS))),]  # done already, redundant


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

sin$cat = NA
sin$cat[which(sin$bw_zsc_Marsal<(-1.5))] = "SGA"  ###  or -2
sin$cat[which(sin$bw_zsc_Marsal>1.5)] = "LGA"  ### or +2

table(sin$cat)

sin = sin[which(!is.na(sin$cat)),]
sin$clr = ifelse(sin$cat=="LGA","orange","cornflowerblue")
plot(sin$BVIKTBS ~ sin$GRDBS,col=sin$clr,pch=19,cex=0.5,
     xlab="gestational age at birth (days)",ylab="birth weight (g)")
legend(x = 170,y = 5000,legend = c("LGA","SGA"),pch=19,col=c("orange","cornflowerblue"),cex=0.8,bty = "n")
grid()


coef(summary(lm(sin$GRDBS ~ sin$MLANGD)))
coef(summary(lm(twi$GRDBS ~ twi$MLANGD)))

coef(summary(lm(GRDBS ~ MLANGD, data=sin[which(sin$bw_zsc_Marsal>1.5),])))
coef(summary(lm(GRDBS ~ MLANGD, data=sin[which(sin$bw_zsc_Marsal<(-1.5)),])))


##### model with separate groups

mod0 = lm(GRDBS ~ MLANGD + MALDER + KON, data=sin[which(sin$cat=="SGA"),])
coef(summary(mod0)) # 0.088 beta
mod1 = lm(GRDBS ~ MLANGD + MALDER + KON, data=sin[which(sin$cat=="LGA"),])
coef(summary(mod1)) # 0.161 beta

##### model with linear interaction 
sin$cat = factor(sin$cat,levels = c("SGA","LGA"))
mod2 = lm(GRDBS ~ MLANGD*cat + KON + MALDER,data=sin)
summary(mod2) 
coef(summary(mod2))

table(sin$cat,useNA = "a")


########################################################################
# Figure 3: LGA-vs-SGA in SIMPLE LINEAR INTERACTION MODEL
########################################################################

mod = lm(GRDBS ~ MLANGD * cat + MALDER, data = sin)
coef(summary(mod))
summary(mod)
fun = function(x3,y3) {
        df = expand.grid(cat=c("LGA","SGA"),MLANGD=140:196,MALDER=29)
        df = df[order(df$cat,df$MLANGD),]
        prd = predict(mod,df,interval = "confidence",level = 0.95) # 
        df$nw = prd[,1] # new value (predicted one)
        tst = df$cat=="SGA"
        #points(x=df$MLANGD[tst],y=df$nw[tst],col="cornflowerblue",type="l")
        plot(df$MLANGD[tst],df$nw[tst],col="cornflowerblue",type="l",ylim=c(min(df$nw)*0.98,max(df$nw)*1.01),
             ylab="child's gestational age at birth (days)",xlab="maternal height (cm)")
        grid()
        points(df$MLANGD[!tst],df$nw[!tst],col="orange",type="l")
        xa = x3[1]-3; xb = x3[3]+12
        ya = min(df$nw)*0.98+1; yb = max(y3)+2.5
        polygon(x = c(xa,xb,xb,xa,xa),y=c(ya,ya,yb,yb,ya),col = "white",border = NA)
        
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
        x = x3[1]; dx = 5.0; xs = c(x,x+dx,x+dx,x,x) # x = 152
        y = y3[1]; dy = 0.5; ys = c(y-dy,y-dy,y+dy,y+dy,y-dy) # y = 269
        polygon(xs,ys,col=ciclr_sga,border = NA)
        segments(x,y,x+dx,y,col="cornflowerblue") # SGA
        polygon(xs,ys-1.3,col=ciclr_lga,border = NA)
        segments(x,y-1.3,x+dx,y-1.3,col="orange")
        text(x-2,y+1.5,"model, 95% CI",pos = 4,cex = 0.7)
        text(x+5,y,"SGA",pos = 4,cex = 0.7)
        text(x+5,y-1.3,"LGA",pos = 4,cex = 0.7)
        
        ###### empirical plot
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
        x = x3[2]; dy = 0.5; y = y3[2] #y = 269; x = 168
        text(x-5,y+1.5,"mean, 95% CI",pos = 4,cex = 0.7)
        points(c(x,x-0.2),c(y,y-1.3),pch=19,cex=1.2,col=c("cornflowerblue","orange"))
        segments(x,y-0.5,x,y+0.8,col="darkgrey") # SGA
        segments(x-0.2,y-1.3-0.8,x-0.2,y-1.3+0.5,col="darkgrey") # LGA
        text(x,y,"SGA",pos = 4,cex = 0.7)
        text(x,y-1.3,"LGA",pos = 4,cex = 0.7)
        
        # LEGEND 3
        x = x3[3]; y = y3[3] #y = 271; x = 179
        text(x-1,y+3,"N (observations)",pos = 4,cex = 0.7)
        szs = c(1e2,1e3,1e4)
        legend(x,y+3,legend = szs,pch=1,bty = "n",pt.cex = log(szs,50))
        
}
pdf("~/Downloads/LGA_vs_SGA.pdf",width = 5.2,height = 5.5)
fun(x3=c(159,175,183),y3=c(271,271,269.5))
dev.off()




########################################################################
# Supplemantary Figure XX: LGA-vs-SGA in POLYNOMIAL INTERACTION MODEL
########################################################################

######  quadratic-model-based plot
mod = lm(GRDBS ~ poly(MLANGD,2) * cat, data = sin)
summary(mod)
cfs = coef(summary(mod)); cfs
write.table(cfs,"~/Downloads/LGA_SGA_interaction_coefs.txt",quote=F,sep="\t")

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

###### empirical plot
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




