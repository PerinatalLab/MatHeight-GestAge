

## uterine size and fetal size interaction. assumption: tall fathers = larger fetus
## i expect to observe shorter gestAge when mother is short and father is tall
## data: MOBA
## eseence: interaction between maternal and paternal height and gestAge

.. update the data cleaning
.. should high/low BMI females be excluded?


library(gbm)

mfr = read.table("~/Biostuff/TERESA_PROTEOMICS_PDB540/MoBa_v8_PDB540/Data/PDB540_MFR_410_v8.csv",h=T,sep=",",stringsAsFactors = F)
mfr = mfr[which(mfr$PLURAL==1),]
mfr = mfr[,c("PREG_ID_540","SVLEN_DG")]

q1 = read.table("~/Biostuff/TERESA_PROTEOMICS_PDB540/MoBa_v8_PDB540/Data/PDB540_Skjema1_v8_HeightWeight.csv",h=T,sep=",",stringsAsFactors = F)
q1 = q1[which((q1$AA87>140)&(q1$AA88>150)&(q1$AA85<140)&(q1$AA85>35)),]
q1$bmi = q1$AA85/(q1$AA87/100)^2
m = merge(mfr,q1,by="PREG_ID_540",all=F)

m = m[which((!is.na(m$SVLEN_DG))&(!is.na(m$AA87))&(!is.na(m$AA88))),]

head(m)
plot(m$AA87,m$AA88)
summary(lm(m$AA87~m$AA88))

summary(lm(m$SVLEN_DG~m$AA88+ m$AA87))
summary(lm(m$SVLEN_DG~m$AA88*m$AA87))


gbm1 = gbm(SVLEN_DG ~ AA88 + AA87, data=m,
            var.monotone=c(0,0), # -1: monotone decrease,
            distribution="gaussian", # see the help for other choices
            n.trees=5000, # number of trees
            shrinkage=0.01, # shrinkage or learning rate, 0.05
            # 0.001 to 0.1 usually work
            interaction.depth=3, # 1: additive model, 2: two-way interactions, etc.
            bag.fraction = 0.5, # subsampling fraction, 0.5 is probably best
            train.fraction = 0.5, # fraction of data for training,
            # first train.fraction*N used for training
            n.minobsinnode = 10, # minimum total weight needed in each node
            cv.folds = 3, # do 3-fold cross-validation
            keep.data=FALSE, # keep a copy of the dataset with the object
            verbose=FALSE, # don't print out progress
            n.cores=1) # use only a single core (detecting #cores is
# error-prone, so avoided here)

best.iter <- gbm.perf(gbm1,method="OOB")
print(best.iter)
# check performance using a 50% heldout test set
best.iter <- gbm.perf(gbm1,method="test")
print(best.iter)
# check performance using 5-fold cross-validation
best.iter <- gbm.perf(gbm1,method="cv")
print(best.iter)

# plot the performance # plot variable influence
summary(gbm1,n.trees=1) # based on the first tree
summary(gbm1,n.trees=best.iter) # based on the estimated best number of trees
# compactly print the first and last trees for curiosity
print(pretty.gbm.tree(gbm1,1))
print(pretty.gbm.tree(gbm1,gbm1$n.trees))

plot(gbm1,1,best.iter)
plot(gbm1,2,best.iter)

# contour plot of variables 1 and 2 after "best" iterations
plot(gbm1,2:1,best.iter)
plot(gbm1,2:3,best.iter)
plot(gbm1,c(3,1),best.iter)
plot(gbm1,1:3,best.iter)

sum((m$bmi<18)&(m$AA88>193))
median(m$SVLEN_DG[which((m$bmi<18)&(m$AA88>193))])
