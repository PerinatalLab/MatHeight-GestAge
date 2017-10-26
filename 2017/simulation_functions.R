# helper functions for simulating gestational age (height paper)
# 2017-09-25  Jonas.B
# location: "~/Dropbox/GIT/SEMFR_HEIGHT-GESTAGE/2017/"
# dependent script: "simulation.Rmd"

library(ggplot2)
library(dplyr)

#####  0  (illustrate 8 combinations of correlation scenarios between three variables)
fun_plotScen = function() {
        #par(mar=c(0,0,0,0)) # mfrow=c(1,1)
        plot(NA,xlim=c(0,16),ylim=c(0,16),xaxt="n",yaxt="n",xlab=NA,ylab=NA)
        offs_sz = 2
        step_sz = 5
        xs = c(0, 1, 2) + offs_sz
        ys = c(0, sqrt(2^2-1^2), 0) + offs_sz 
        
        pcrds = matrix(NA,nc=6,nr=3*3)  # point coordinates
        row_count = 0
        for (i in 0:2) {
                for (j in 0:2) {
                        if (!((i==1)&(j==1))) {
                                # define coordinates of the cluster
                                xss = xs+step_sz*i
                                yss = ys+step_sz*j
                                
                                # draw clusters and connections
                                lines(x = c(xss,xss[1]),y = c(yss,yss[1]),col="grey")
                                points(xss,yss)
                                points(xss,yss,pch=19,col="white",cex=0.8)
                                
                                # record coordinates
                                row_count = row_count + 1
                                pcrds[row_count,] = c(xss,yss)
                                
                                # add text
                                crct = 1 #0.5
                                xss = xss + c(0,0,0)
                                yss = yss + c(-crct,crct,-crct)
                                text(c("MH","FGR","MUS"),x = xss,y=yss,cex = 0.5)
                                
                                # save the center-of-the-figure coordinates 
                        } else {
                                center_coords = c( mean(xs+step_sz*i),mean(range(ys+step_sz*j)) )
                        }
                }
        }
        
        
        # scenario 1
        xs = pcrds[1,1:3]
        ys = pcrds[1,4:6]
        lines(x=xs[c(1,3,2)],y=ys[c(1,3,2)],col="cornflowerblue",lwd=5)
        points(x=xs[c(1,3,2)],y=ys[c(1,3,2)],col="cornflowerblue",pch=19,cex=1.1)
        
        # scenario 2
        xs = pcrds[2,1:3]
        ys = pcrds[2,4:6]
        lines(x=xs[c(2,1,3)],y=ys[c(2,1,3)],col="cornflowerblue",lwd=5)
        points(x=xs[c(2,1,3)],y=ys[c(2,1,3)],col="cornflowerblue",pch=19,cex=1.1)
        
        # scenario 3
        xs = pcrds[3,1:3]
        ys = pcrds[3,4:6]
        lines(x=xs[c(1,2,3)],y=ys[c(1,2,3)],col="cornflowerblue",lwd=5)
        points(x=xs[c(1,2,3)],y=ys[c(1,2,3)],col="cornflowerblue",pch=19,cex=1.1)
        
        # scenario 4
        xs = pcrds[4,1:3]
        ys = pcrds[4,4:6]
        lines(x=c(xs,xs[1]),y=c(ys,ys[1]),col="cornflowerblue",lwd=5)
        points(x=c(xs,xs[1]),y=c(ys,ys[1]),col="cornflowerblue",pch=19,cex=1.1)
        
        # scenario 5
        #  -> all disconnected
        
        # scenario 6
        xs = pcrds[6,1:3]
        ys = pcrds[6,4:6]
        lines(x=xs[1:2],y=ys[1:2],col="cornflowerblue",lwd=5)
        points(x=xs[1:2],y=ys[1:2],col="cornflowerblue",pch=19,cex=1.1)
        
        # scenario 7
        xs = pcrds[7,1:3]
        ys = pcrds[7,4:6]
        lines(x=xs[2:3],y=ys[2:3],col="cornflowerblue",lwd=5)
        points(x=xs[2:3],y=ys[2:3],col="cornflowerblue",pch=19,cex=1.1)
        
        # scenario 8
        xs = pcrds[8,1:3]
        ys = pcrds[8,4:6]
        lines(x=xs[c(1,3)],y=ys[c(1,3)],col="cornflowerblue",lwd=5)
        points(x=xs[c(1,3)],y=ys[c(1,3)],col="cornflowerblue",pch=19,cex=1.1)
        
        ### circle in the middle
        
        x_add = cos( (pi/4)*c(2:0,7:3))
        y_add = sin( (pi/4)*c(2:0,7:3))  #rev((1:8)-1))
        
        x_circ = center_coords[1] + x_add*2 
        y_circ = center_coords[2] + y_add*2
        text(1:8,x = x_circ,y = y_circ,cex=0.9,col="orange")
        #text(labels = "scenarios",x = center_coords[1],y = center_coords[2],cex=0.7,col="orange")
        
}


#####  1  (simulate correlated data in height project)
fun_generateCorData = function(n_ind,mean_fgr,seed,scen) {
        
        # n_ind = number of individuals/pregnancies to simulate
        # mean_fgr = the mid-value of fetal growth rate (cm3/d). ~10
        # seed = reproducibility
        # scen = scenario  (values 1-8)
        
        
        # assumptions:
        # when correlated, dependencies are linear
        # fetal growth rate is constant in time
        
        set.seed(seed)
        
        ###### GENERAL SECTION
        
        ### generate source for bimodal continuous (extreme) fetal growth values:
        xtr_fgr = 0.5 # fraction of all fetal growths that will be considered extreme
        ps1 = seq(from = 1/n_ind,     to = xtr_fgr/2, length.out = floor(n_ind/2))
        ps2 = seq(from = 1-xtr_fgr/2, to = 1-1/n_ind, length.out = floor(n_ind/2))
        fgr_src = qnorm(p = c(ps1,ps2),mean = mean_fgr, sd=  0.5)  #fetal growth rate (g/day)
        fgr_src = sample(fgr_src,size = length(fgr_src),replace = F) # shuffle
        
        # generate source for maternal height values
        mhs_src = rnorm(n_ind,mean=160,sd=10)  # maternal heights
        
        # generate source for maternal uterine size values
        mus_src = rnorm(n_ind,mean=298*10,sd=100)  # maternal uterine size
        
        ###### SPECIFIC SECTIONS
        
        ### Scenario 1
        if (scen==1) {
                # maternal height values
                mhs = sample(mhs_src,replace = F)
                # maternal uterine size values
                mus = sample(mus_src,replace = F)
                # fetal growth rate values
                fgr = sample(fgr_src,replace = F)
        }
        
        ### Scenario 2
        if (scen==2) {
                # generate helper (genetic/confoudning) variable
                hlp = rnorm(n_ind)
                hlp1 = hlp + rnorm(n_ind)  # "offspring"
                hlp2 = hlp + rnorm(n_ind)  # "offspring"
                # maternal height values
                rnk = match(rank(mhs_src),rank(hlp1))
                mhs = mhs_src[order(rnk)]; rm(rnk)
                # maternal uterine size values
                rnk = match(rank(mus_src),rank(hlp2))
                mus = mus_src[order(rnk)]; rm(rnk)
                # fetal growth rate values
                fgr = sample(fgr_src,replace = F)
                ### cleanup
                rm(hlp,hlp1,hlp2)
        }
        
        ### Scenario 3
        if (scen==3) {
                # generate helper (genetic/confoudning) variable
                hlp = rnorm(n_ind) # "source"
                hlp1 = hlp + rnorm(n_ind)  # "offspring"
                hlp2 = hlp + rnorm(n_ind)  # "offspring"
                # maternal height values
                mhs = sample(mhs_src,replace = F)
                # maternal uterine size values
                rnk = match(rank(mus_src),rank(hlp1))
                mus = mus_src[order(rnk)]; rm(rnk)
                # fetal growth rate values
                rnk = match(rank(fgr_src),rank(hlp2))
                fgr = fgr_src[order(rnk)]; rm(rnk)
                ### cleanup
                rm(hlp,hlp1,hlp2)
        }
        
        ### Scenario 4
        if (scen==4) {
                # generate helper (genetic/confoudning) variable
                hlp = rnorm(n_ind) # "source"
                hlp1 = hlp + rnorm(n_ind)  # "offspring"
                hlp2 = hlp + rnorm(n_ind)  # "offspring"
                # maternal height values
                rnk = match(rank(mhs_src),rank(hlp1))
                mhs = mhs_src[order(rnk)]; rm(rnk)
                # maternal uterine size values
                mus = sample(mus_src,replace = F)
                # fetal growth rate values
                rnk = match(rank(fgr_src),rank(hlp2))
                fgr = fgr_src[order(rnk)]; rm(rnk)
                ### cleanup
                rm(hlp,hlp1,hlp2)
        }
        
        ### Scenario 5
        if (scen==5) {
                # generate helper (genetic/confoudning) variable
                hlp = rnorm(n_ind) # "source"
                hlp1 = hlp + rnorm(n_ind)  # "offspring"
                hlp2 = hlp + rnorm(n_ind)  # "offspring"
                hlp3 = hlp + rnorm(n_ind)  # "offspring"
                # maternal height values
                rnk = match(rank(mhs_src),rank(hlp1))
                mhs = mhs_src[order(rnk)]; rm(rnk)
                # maternal uterine size values
                rnk = match(rank(mus_src),rank(hlp2))
                mus = mus_src[order(rnk)]; rm(rnk)
                # fetal growth rate values
                rnk = match(rank(fgr_src),rank(hlp3))
                fgr = fgr_src[order(rnk)]; rm(rnk)
                ### cleanup
                rm(hlp,hlp1,hlp2,hlp3)
        }
        
        ### Scenario 6
        if (scen==6) {
                # generate helper (genetic/confoudning) variable
                hlpA = rnorm(n_ind) # "source"
                hlpB = rnorm(n_ind) # "source"
                hlp1 = hlpA + rnorm(n_ind)  # "offspring"
                hlp2 = hlpA + hlpB + rnorm(n_ind)  # "central offspring"
                hlp3 = hlpB + rnorm(n_ind)  # "offspring"
                # maternal height values
                rnk = match(rank(mhs_src),rank(hlp1))
                mhs = mhs_src[order(rnk)]; rm(rnk)
                # maternal uterine size values
                rnk = match(rank(mus_src),rank(hlp2))
                mus = mus_src[order(rnk)]; rm(rnk)
                # fetal growth rate values
                rnk = match(rank(fgr_src),rank(hlp3))
                fgr = fgr_src[order(rnk)]; rm(rnk)
                ### cleanup
                rm(hlpA,hlpB,hlp1,hlp2,hlp3)
        }
        
        
        ### Scenario 7
        if (scen==7) {
                # generate helper (genetic/confoudning) variable
                hlpA = rnorm(n_ind) # "source"
                hlpB = rnorm(n_ind) # "source"
                hlp1 = hlpA + rnorm(n_ind)  # "offspring"
                hlp2 = hlpA + hlpB + rnorm(n_ind)  # "central offspring"
                hlp3 = hlpB + rnorm(n_ind)  # "offspring"
                # maternal height values
                rnk = match(rank(mhs_src),rank(hlp2))
                mhs = mhs_src[order(rnk)]; rm(rnk)
                # maternal uterine size values
                rnk = match(rank(mus_src),rank(hlp1))
                mus = mus_src[order(rnk)]; rm(rnk)
                # fetal growth rate values
                rnk = match(rank(fgr_src),rank(hlp3))
                fgr = fgr_src[order(rnk)]; rm(rnk)
                ### cleanup
                rm(hlpA,hlpB,hlp1,hlp2,hlp3)
        }
        
        ### Scenario 8
        if (scen==8) {
                # generate helper (genetic/confoudning) variable
                hlpA = rnorm(n_ind) # "source"
                hlpB = rnorm(n_ind) # "source"
                hlp1 = hlpA + rnorm(n_ind)  # "offspring"
                hlp2 = hlpA + hlpB + rnorm(n_ind)  # "central offspring"
                hlp3 = hlpB + rnorm(n_ind)  # "offspring"
                # maternal height values
                rnk = match(rank(mhs_src),rank(hlp1))
                mhs = mhs_src[order(rnk)]; rm(rnk)
                # maternal uterine size values
                rnk = match(rank(mus_src),rank(hlp3))
                mus = mus_src[order(rnk)]; rm(rnk)
                # fetal growth rate values
                rnk = match(rank(fgr_src),rank(hlp2))
                fgr = fgr_src[order(rnk)]; rm(rnk)
                ### cleanup
                rm(hlpA,hlpB,hlp1,hlp2,hlp3)
        }
        
        
        # status (LGA vs SGA)
        sta = as.numeric(fgr>median(fgr))  # 0 = SGA, 1 = LGA
        
        ## assign globally
        fgr<<-fgr
        mus<<-mus
        mhs<<-mhs
        sta<<-sta
}


#####  2  (generate gestational age data, i.e., birth events)
fun_generateGestAge = function(fgr,mus,interactionRisk, baselineRisk) {
        
        
        # fgr (continuous) = fetal growth rate (cm3/d)
        
        # mus (continuous) = maternal uterine size (cm3)
        
        # interactionRisk = TRUE/FALSE; defines whether an interaction
        #       risk to be born on a specific day should be applied
        
        # baselineRisk = TRUE/FALSE; defines whether a baseline risk
        #               to be born on a specific day should be applied
        
        n_ind = length(mus)
        gas = rep(NA,n_ind) # here values of gestational age will be placed
        
        # run simulation going day-by-day through gestational calendar
        count = 0 # day count
        while (any(is.na(gas))) {
                
                count = count + 1 # day count
                fsz = fgr*count  # fetal sizes on this specific day (g)
                
                #  OPTIONAL:  interaction-imposed probability to be born on this day
                if (interactionRisk==TRUE) {
                        pb = 1/(1 + exp(1)^(-log(fsz/mus,1.01) )) # Probability of Birth day
                        fun = function(x) sample(c(F,T),size = 1,prob = c(1-x,x)) # helper 
                        birth_event_1 = unlist(lapply(pb,fun)) # based on fetal size
                } else {
                        birth_event_1 = rep(FALSE, n_ind)
                }
                
                #  OPTIONAL:  baseline probability to be born on this day
                if (baselineRisk==TRUE) {
                        pb = 1 / ( 1 + exp(1)^(log(298/count,1.025)) ) # Probability of Birth
                        birth_event_2 = sample(c(F,T),size=n_ind,prob = c(1-pb,pb),replace = T)
                } else {
                        birth_event_2 = rep(FALSE,n_ind)
                }
                
                # indexes of pregnancies to be delivered today
                ixs = which( is.na(gas) & (birth_event_1 | birth_event_2) )
                if(length(ixs)>0) gas[ixs] = count   # asign birth age (gestAge at birth)
                
        }
        gas
}


#####  3  (visualizeinteraction)
fun_interactionPlots = function(mod,mhs,sta,gas,tit) {
        # mod = linear model (GA ~ MH)
        # mhs = maternal heights
        # sta = status (LGA/SGA; generated from fgr)
        # gas = simulated gestational-age-at-birth values
        # tit = index for plot title
        
        df = expand.grid(sta=c(0,1),mhs=floor(min(mhs)):floor(max(mhs)))
        df = df[order(df$sta,df$mhs),]
        
        # model the gestational age (also confidence intervals)
        prd = predict(mod,df,interval = "confidence",level = 0.95)
        df$nw = prd[,1] # nw = predicted gestAge values
        
        
        # define transparent colors for the scatter plot
        clr = rep(rgb(1,0,0,0.4),length(gas)) # red
        clr[sta==0] = rgb(0,0,0,0.4) # black
        
        #  scatter plot of simulated data
        plot(x = mhs,y=gas,col=clr,#main=paste("interaction pattern",tit),
             ylab="child's gestAge at birth",
             xlab="maternal height")
        grid()
        # add modeled lines
        points(df$mhs[df$sta==0],df$nw[df$sta==0],col="cornflowerblue",type="l") # SGA
        points(df$mhs[df$sta==1],df$nw[df$sta==1],col="orange",type="l") # LGA
        
        
        ##### add confidence intervals for the model curves
        # for SGA
        xs = c(df$mhs[df$sta==0],rev(df$mhs[df$sta==0]),df$mhs[df$sta==0][1])
        ys = c(prd[df$sta==0,2],rev(prd[df$sta==0,3]),prd[df$sta==0,2][1])
        clrcm = col2rgb("cornflowerblue") / 255 # color components
        ciclr_sin = rgb(clrcm[1,1],clrcm[2,1],clrcm[3,1],0.4) # color for conf. intervals
        polygon(xs,ys,col=ciclr_sin,border = NA)
        
        
        # for LGA
        xs = c(df$mhs[df$sta==1],rev(df$mhs[df$sta==1]),df$mhs[df$sta==1][1])
        ys = c(prd[df$sta==1,2],rev(prd[df$sta==1,3]),prd[df$sta==1,2][1])
        clrcm = col2rgb("orange") / 255 # color components
        ciclr_twi = rgb(clrcm[1,1],clrcm[2,1],clrcm[3,1],0.4) # color for conf. intervals
        polygon(xs,ys,col=ciclr_twi,border = NA)
        
}



# gradient histogram
fun_gradHist = function(x,y,BRKS,XLAB,YLAB) {
        
        # x = gest age
        # y = e.g., fetal growth rate and maternal uterine size ratio
        # BRKS = number of breaks in the histogram
        bin_size = round((max(x)-min(x))/BRKS*1.01)  # ***
        h = hist(x,breaks=BRKS,plot = F)
        hlp = data.frame(from=h$breaks,till=lead(h$breaks,1),stringsAsFactors = F)
        hlp = hlp[-nrow(hlp),]
        ctg = rep(NA,length(x))
        for (i in 1:length(h$breaks)) {
                ctg[which( (x>h$breaks[i])&(x<=h$breaks[i+1]))] = i
        }
        hlp1 = data.frame(x,y,ctg,stringsAsFactors = F)
        hlp2 = group_by(hlp1,ctg) %>% summarise(n=n(),col=mean(y)) %>% ungroup()
        hlp2 = as.data.frame(hlp2)
        hlp3 = data.frame(ctg=seq(length(h$mids)),counts=h$counts,mid=h$mids,stringsAsFactors = F)
        dat = merge(hlp2,hlp3,by="ctg",all=T)
        dat = dat[order(dat$ctg),]
        head(dat)
        class(dat$mid)
        dat$counts = as.numeric(dat$counts)
        ggplot(dat, aes(x = mid, y = counts, fill =col)) + 
                geom_bar(stat = "identity",alpha = 0.8,width = bin_size) + xlab(XLAB) + ylab(YLAB) +
                scale_fill_gradient(low="yellow", high="red") +
                theme_minimal() + 
                scale_x_continuous(breaks = seq(220,300,10),labels = seq(220,300,10))
        #scale_x_continuous(breaks = seq(-1,1,0.25), labels = seq(-1,1,0.25)) +
        
        
}



