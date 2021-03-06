#--------------------divide into test and training data----------------------
train_size = 0.7

n = nrow(multi.regression.compounds)
#------------collect the data with n*train_size from the dataset------------
perm = sample(n, size = round(n * train_size))

#-------------------training data----------------------------------------
multi.regression.compounds.train <- multi.regression.compounds[perm, ]
preprocessed.y.train <- multi.regression.compounds.train[,c(1)]
multi.regression.x.train <- multi.regression.compounds.train[,-c(1)]
#-----------------------test data----------------------------------------
multi.regression.compounds.test <-multi.regression.compounds[-perm, ]
preprocessed.y.test <- multi.regression.compounds.test[,c(1)]
multi.regression.x.test <- multi.regression.compounds.test[,-c(1)]

#-----------transform into data frame--------------------------
multi.regression.compounds.train <- as.data.frame(multi.regression.compounds.train)

#----------------------MLR regression training--------------------------------------------
compounds.lm <- lm(preprocessed.y~., data=multi.regression.compounds.train)
compounds.lm

summary(compounds.lm)

lm.predicted.y <- predict(compounds.lm)
lm.predicted.y

cor(preprocessed.y.train, lm.predicted.y)
lm.r2 <- cor(preprocessed.y.train, lm.predicted.y)**2
lm.r2

plot(preprocessed.y.train, lm.predicted.y,
     xlab="Observed value",
     ylab="Predicted value", main = "MLR")
abline(a=0, b=1)

#-------------------------MLR test--------------------------------
lm.predicted.y.test <- predict(compounds.lm,newdata = multi.regression.x.test)
lm.predicted.y.test

cor(preprocessed.y.test, lm.predicted.y.test)
lm.r2.test <- cor(preprocessed.y.test, lm.predicted.y.test)**2
lm.r2.test

plot(preprocessed.y.test, lm.predicted.y.test,
     xlab="Observed value",
     ylab="Predicted value", main = "MLR test")
abline(a=0, b=1)


#------------------------PLS training------------------------------------
library(pls)

compounds.plsr <- plsr(preprocessed.y~., data=multi.regression.compounds.train, validation="CV")
summary(compounds.plsr)
plot(compounds.plsr, "validation")

ncomp.onesigma <- selectNcomp(compounds.plsr, method = "randomization", plot = TRUE, ylim = c(.1, 1))
ncomp.onesigma

predict(compounds.plsr)[, , ncomp.onesigma]
plsr.predicted.y <- predict(compounds.plsr)[, , ncomp.onesigma]
plsr.r2 <- cor(multi.regression.compounds.train[,c(1)], plsr.predicted.y)**2
plsr.r2

plot(multi.regression.compounds.train[,c(1)], plsr.predicted.y,
     xlab="Observed value",
     ylab="Predicted value",
     main="PLSR")
abline(a=0, b=1)

compounds.plsr$coefficients[, , ncomp.onesigma]

#-------------------------pls test--------------------------------
pls.predicted.y.test <- predict(compounds.plsr,newdata = multi.regression.x.test)[,, ncomp.onesigma]
cor(preprocessed.y.test, pls.predicted.y.test)
pls.r2.test <- cor(preprocessed.y.test, pls.predicted.y.test)**2
pls.r2.test

plot(preprocessed.y.test, pls.predicted.y.test,
     xlab="Observed value",
     ylab="Predicted value", main = "PLS test")
abline(a=0, b=1)

#--------------------------------------PLS-VIP training-------------------------
library(plsVarSel)

vip.selected <- bve_pls(preprocessed.y.train, multi.regression.x.train, ncomp = ncomp.onesigma, VIP.threshold = 0.8) 
vip.selected

x.vip.train <- multi.regression.x.train[,c(vip.selected$bve.selection)]
x.vip.test <- multi.regression.x.test[,c(vip.selected$bve.selection)]
multi.regression.compounds.train.vip <- cbind.data.frame(preprocessed.y.train,x.vip.train)
multi.regression.compounds.test.vip <- cbind.data.frame(preprocessed.y.test,x.vip.test)

compounds.plsr.vip <- plsr(preprocessed.y.train~., data=multi.regression.compounds.train.vip, validation="CV")
summary(compounds.plsr.vip)
plot(compounds.plsr.vip, "validation")

ncomp.onesigma.vip <- selectNcomp(compounds.plsr.vip, method = "randomization", plot = TRUE, ylim = c(.1, 1))
ncomp.onesigma.vip

predict(compounds.plsr.vip)[, , ncomp.onesigma.vip]
plsr.predicted.y.vip <- predict(compounds.plsr.vip)[, , ncomp.onesigma.vip]
plsr.r2.vip <- cor(preprocessed.y.train, plsr.predicted.y.vip)**2
plsr.r2.vip

plot(preprocessed.y.train, plsr.predicted.y.vip,
     xlab="Observed value",
     ylab="Predicted value",
     main="PLS-VIP")
abline(a=0, b=1)

compounds.plsr.vip$coefficients[, , ncomp.onesigma.vip]

#-------------------------pls-vip test--------------------------------
pls.predicted.y.vip.test <- predict(compounds.plsr.vip,newdata = multi.regression.x.test)[,, ncomp.onesigma.vip]
cor(preprocessed.y.test, pls.predicted.y.vip.test)
pls.r2.vip.test <- cor(preprocessed.y.test, pls.predicted.y.vip.test)**2
pls.r2.vip.test

plot(preprocessed.y.test, pls.predicted.y.vip.test,
     xlab="Observed value",
     ylab="Predicted value", main = "PLS-VIP test")
abline(a=0, b=1)


#-----------------------------Uninformative Variable Elimination in PLS training----------------------
mcuve <- mcuve_pls(preprocessed.y.train, multi.regression.x.train, ncomp = ncomp.onesigma, N = 3)
mcuve

x.uve.train <- multi.regression.x.train[,c(mcuve$mcuve.selection)]
x.uve.test <- multi.regression.x.test[,c(mcuve$mcuve.selection)]
multi.regression.compounds.train.uve <- cbind.data.frame(preprocessed.y.train,x.uve.train)
multi.regression.compounds.test.uve <- cbind.data.frame(preprocessed.y.test,x.uve.test)

compounds.plsr.uve <- plsr(preprocessed.y.train~., data=multi.regression.compounds.train.uve, validation="CV")
summary(compounds.plsr.uve)
plot(compounds.plsr.uve, "validation")

ncomp.onesigma.uve <- selectNcomp(compounds.plsr.uve, method = "randomization", plot = TRUE, ylim = c(.1, 1))
ncomp.onesigma.uve

predict(compounds.plsr.uve)[, , ncomp.onesigma.uve]
plsr.predicted.y.uve <- predict(compounds.plsr.uve)[, , ncomp.onesigma.uve]
plsr.r2.uve <- cor(preprocessed.y.train, plsr.predicted.y.uve)**2
plsr.r2.uve

plot(preprocessed.y.train, plsr.predicted.y.uve,
     xlab="Observed value",
     ylab="Predicted value",
     main="PLS-UVE")
abline(a=0, b=1)

compounds.plsr.uve$coefficients[, , ncomp.onesigma.uve]

#-------------------------pls-uve test--------------------------------
pls.predicted.y.uve.test <- predict(compounds.plsr.uve,newdata = multi.regression.x.test)[,, ncomp.onesigma.uve]
cor(preprocessed.y.test, pls.predicted.y.uve.test)
pls.r2.uve.test <- cor(preprocessed.y.test, pls.predicted.y.uve.test)**2
pls.r2.uve.test

plot(preprocessed.y.test, pls.predicted.y.uve.test,
     xlab="Observed value",
     ylab="Predicted value", main = "PLS-uve test")
abline(a=0, b=1)


#-----------------------------GAPLS training-----------------------------------
library(gaselect)

n = ncol(multi.regression.x.train)
l = cbind(ncomp.onesigma + 10,n)

ctrl <- genAlgControl(populationSize = 100, numGenerations = 100, minVariables = 1,
                      maxVariables = min(l), verbosity = 1)

evaluatorRDCV <- evaluatorPLS(numReplications = 2, innerSegments = 5, outerSegments = 3,
                              numThreads = 1)

# Generate demo-data
set.seed(12345)
X <- as.matrix(multi.regression.x.train)
#View(X)
y <- (preprocessed.y.train)

resultRDCV <- genAlg(y, X, control = ctrl, evaluator = evaluatorRDCV, seed = 123)

RD <- subsets(resultRDCV, 1, names = FALSE)
RD$`1`

x.ga.train <- multi.regression.x.train[,c(RD$`1`)]
x.ga.test <- multi.regression.x.test[,c(RD$`1`)]
multi.regression.compounds.train.ga <- cbind.data.frame(preprocessed.y.train,x.ga.train)
multi.regression.compounds.test.ga <- cbind.data.frame(preprocessed.y.test,x.ga.test)

compounds.plsr.ga <- plsr(preprocessed.y.train~., data=multi.regression.compounds.train.ga, validation="CV")
summary(compounds.plsr.ga)
plot(compounds.plsr.ga, "validation")

ncomp.onesigma.ga <- selectNcomp(compounds.plsr.ga, method = "randomization", plot = TRUE, ylim = c(.1, 1))
ncomp.onesigma.ga

predict(compounds.plsr.ga)[, , ncomp.onesigma.ga]
plsr.predicted.y.ga <- predict(compounds.plsr.ga)[, , ncomp.onesigma.ga]
plsr.r2.ga <- cor(preprocessed.y.train, plsr.predicted.y.ga)**2
plsr.r2.ga

plot(preprocessed.y.train, plsr.predicted.y.ga,
     xlab="Observed value",
     ylab="Predicted value",
     main="PLS-GA")
abline(a=0, b=1)

compounds.plsr.ga$coefficients[, , ncomp.onesigma.ga]

#-------------------------pls-GA test--------------------------------
pls.predicted.y.ga.test <- predict(compounds.plsr.ga,newdata = multi.regression.x.test)[,, ncomp.onesigma.ga]
cor(preprocessed.y.test, pls.predicted.y.ga.test)
pls.r2.ga.test <- cor(preprocessed.y.test, pls.predicted.y.ga.test)**2
pls.r2.ga.test

plot(preprocessed.y.test, pls.predicted.y.ga.test,
     xlab="Observed value",
     ylab="Predicted value", main = "PLS-GA test")
abline(a=0, b=1)

#---------------------summary-----------------------
summary.train <-cbind(lm.r2,plsr.r2,plsr.r2.uve,plsr.r2.vip,plsr.r2.ga)
View(summary.train)
summary.test <- cbind(lm.r2.test,pls.r2.test,pls.r2.uve.test,pls.r2.vip.test,pls.r2.ga.test)
View(summary.test)

#--------------------テストデータと学習データの予測結果を同時にプロット-----------------------
model.plot <- compounds.plsr.vip #どのモデルを使用するか入力
ncomp.plot <- ncomp.onesigma.vip　#どのモデルで最適化したncompを使用するか入力

ttest <- predict(model.plot, newdata = multi.regression.x.test)[,, ncomp.plot]
ttrain <- predict(model.plot)[,, ncomp.plot]

plot(0, 0, type = "n", xlim = c(0, 2), ylim = c(0,2),xlab = "Observed Value", ylab = "Predicted Value")

points(preprocessed.y.test, ttest, col = "orange", pch = 2)
points(preprocessed.y.train, ttrain, col = "darkgray", pch = 3)
abline(a=0, b=1)

r2ttest <- cor(preprocessed.y.test, ttest)^2
r2ttest

r2ttrain <- cor(preprocessed.y.train, ttrain)^2
r2ttrain

