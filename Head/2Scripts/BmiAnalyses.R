###################################
###################################

rm(list=ls(all=TRUE))

Bmi = read.table("../../Body/1Raw/MainTable.txt", header = TRUE, sep = '\t')
Bmi$Sex = Bmi$Sex-1; table(Bmi$Sex) # make it Dummy: 0 - males, 1 - females
Bmi$BmiRatio = Bmi$BMI2/Bmi$BMI1
Bmi$Bmi2Min1 = (Bmi$BMI2 - Bmi$BMI1)
Bmi$mtDNARatio = Bmi$mtDNA2/Bmi$mtDNA1
Bmi$mtDNA2Min1 = Bmi$mtDNA2 - Bmi$mtDNA1

######## pdf
pdf("../../Body/4Figures/Bmi.mtDNA.pdf")

###### 1: BMI decreased after operation (both types of operations, both genders)
summary(Bmi$BMI1)
summary(Bmi$BMI2)
wilcox.test(Bmi$BMI1,Bmi$BMI2, paired = TRUE)
wilcox.test(Bmi[Bmi$Operation == 'LSG',]$BMI1,Bmi[Bmi$Operation == 'LSG',]$BMI2, paired = TRUE)
wilcox.test(Bmi[Bmi$Operation == 'RYGB',]$BMI1,Bmi[Bmi$Operation == 'RYGB',]$BMI2, paired = TRUE)

wilcox.test(Bmi[Bmi$Sex == 0,]$BMI1,Bmi[Bmi$Sex == 0,]$BMI2, paired = TRUE)
wilcox.test(Bmi[Bmi$Sex == 1,]$BMI1,Bmi[Bmi$Sex == 1,]$BMI2, paired = TRUE) # better

# intercept = 10, coefficient of BMI1 is positive => the higher was BMI before the opearaton, the higher it is still after the operation.
a<-lm(Bmi$BMI2 ~ Bmi$BMI1); summary(a)

# gender - doesn't matter
a1<-lm(Bmi$BMI2 ~ Bmi$BMI1 + Bmi$Sex); summary(a1)

# age!!?? the older the higher BMI2. even intercept is getting nonsignificant!!
a2<-lm(Bmi$BMI2 ~ Bmi$BMI1 + Bmi$Age); summary(a2)
a3<-lm(Bmi$BMI2 ~ 0 + Bmi$BMI1 + Bmi$Age); summary(a3)

# directly mtDNA do not correlate with BMI2
a4<-lm(Bmi$BMI2 ~ Bmi$BMI1 + Bmi$mtDNA2 + Bmi$mtDNA1); summary(a4)
a5<-lm(Bmi$BmiRatio ~ Bmi$mtDNARatio); summary(a5)

a5<-lm(Bmi$BMI1 ~ Bmi$mtDNA1 + Bmi$Age); summary(a5)

a6<-lm(Bmi$BMI2/Bmi$BMI1 ~ Bmi$BMI1); summary(a6) # decrease in BMI is more pronounced in case of high BMI1

a7<-lm(Bmi$BMI2/Bmi$BMI1 ~ Bmi$BMI1 + Bmi$Age); summary(a7) # decrease in BMI is more pronounced in case of high BMI1 and less pronounced with Age
a8<-lm(Bmi$BMI2/Bmi$BMI1 ~ Bmi$BMI1 + Bmi$Age + Bmi$mtDNA2 + Bmi$mtDNA1); summary(a8) # no effect of mtDNA
a9<-lm(Bmi$BMI2/Bmi$BMI1 ~ Bmi$BMI1 + Bmi$Age + Bmi$mtDNA2/Bmi$mtDNA1); summary(a9) # no effect of mtDNA

####### 2: mtDNA copies in controls and cases - before and after operations

boxplot(Bmi[Bmi$Operation == 'Controls' & Bmi$Sex == 0,]$mtDNA1,Bmi[Bmi$Operation == 'Controls' & Bmi$Sex == 1,]$mtDNA1, Bmi[Bmi$Operation != 'Controls' & Bmi$Sex == 0,]$mtDNA1,Bmi[Bmi$Operation != 'Controls' & Bmi$Sex == 1,]$mtDNA1, Bmi[Bmi$Operation != 'Controls' & Bmi$Sex == 0,]$mtDNA2,Bmi[Bmi$Operation != 'Controls' & Bmi$Sex == 1,]$mtDNA2, col = c('light blue', 'pink'))
wilcox.test(Bmi[Bmi$Operation != 'Controls' & Bmi$Sex == 1,]$mtDNA1,Bmi[Bmi$Operation != 'Controls' & Bmi$Sex == 1,]$mtDNA2) # 4.829e-05
wilcox.test(Bmi[Bmi$Operation == 'Controls' & Bmi$Sex == 1,]$mtDNA1,Bmi[Bmi$Operation != 'Controls' & Bmi$Sex == 1,]$mtDNA1) # 0.01726
wilcox.test(Bmi[Bmi$Operation == 'Controls' & Bmi$Sex == 1,]$mtDNA1,Bmi[Bmi$Operation != 'Controls' & Bmi$Sex == 1,]$mtDNA2) # 0.6455


####### 3: mtDNA increased after operation (both types of operations) for females but not males

wilcox.test(Bmi[Bmi$Operation != 'Controls' & Bmi$Sex == 1,]$mtDNA1,Bmi[Bmi$Operation != 'Controls' & Bmi$Sex == 1,]$mtDNA2, paired = TRUE)  # 2.187e-05
wilcox.test(Bmi[Bmi$Operation != 'Controls' & Bmi$Sex == 0,]$mtDNA1,Bmi[Bmi$Operation != 'Controls' & Bmi$Sex == 0,]$mtDNA2, paired = TRUE)  # 0.5566

### is there correlation between the strength of changes in BMI and strength of changes in mtDNA? mainly among females!!

cor.test(Bmi$BMI2/Bmi$BMI1,Bmi$mtDNA2/Bmi$mtDNA1, method = 'spearman', alternative = 'less') 
plot(Bmi$BMI2/Bmi$BMI1,Bmi$mtDNA2/Bmi$mtDNA1)

cor.test(Bmi[Bmi$Sex == 1,]$BMI2/Bmi[Bmi$Sex == 1,]$BMI1,Bmi[Bmi$Sex == 1,]$mtDNA2/Bmi[Bmi$Sex == 1,]$mtDNA1, method = 'spearman', alternative = 'less')  # significant
plot(Bmi[Bmi$Sex == 1,]$BMI2/Bmi[Bmi$Sex == 1,]$BMI1,Bmi[Bmi$Sex == 1,]$mtDNA2/Bmi[Bmi$Sex == 1,]$mtDNA1)

######## plot segments
XMin = min(c(Bmi[!is.na(Bmi$mtDNA1),]$mtDNA1,Bmi[!is.na(Bmi$mtDNA2),]$mtDNA2))
XMax = max(c(Bmi[!is.na(Bmi$mtDNA1),]$mtDNA1,Bmi[!is.na(Bmi$mtDNA2),]$mtDNA2))
YMin = min(c(Bmi[!is.na(Bmi$BMI1),]$BMI1,Bmi[!is.na(Bmi$BMI2),]$BMI2))
YMax = max(c(Bmi[!is.na(Bmi$BMI1),]$BMI1,Bmi[!is.na(Bmi$BMI2),]$BMI2))


summary(Bmi$BmiRatio) # all are less than one => segments are going down always
summary(Bmi$mtDNARatio) # in some rare cases mtDNA is decreasing, mark them by red
plot(NA, xlim=c(XMin,XMax), ylim=c(YMin,YMax), xlab='mtDNA (copy numbers)', ylab="BMI", main = '') # , xaxt="n")
# axis(side = 1, at=c(1:13), labels=c(Gene), las = 2) 
for (i in 1:nrow(Bmi))
  { # i = 1
  if (!is.na(Bmi$mtDNA1[i]) & !is.na(Bmi$mtDNA2[i]) & !is.na(Bmi$BMI1[i]) & !is.na(Bmi$BMI2[i]))
    {
    #if (Bmi$mtDNARatio[i] >= 1) 
    #  {
      if (Bmi$Sex[i] == 1) {arrows(Bmi$mtDNA1[i], Bmi$BMI1[i], Bmi$mtDNA2[i], Bmi$BMI2[i], col = 'pink', lwd = 1)}
      if (Bmi$Sex[i] == 0) {arrows(Bmi$mtDNA1[i], Bmi$BMI1[i], Bmi$mtDNA2[i], Bmi$BMI2[i], col = 'blue', lwd = 1)}
      #segments(Bmi$mtDNA1[i], Bmi$BMI1[i], Bmi$mtDNA2[i], Bmi$BMI2[i], col = 'grey', lwd = 1)
    #  }
    # if  (Bmi$mtDNARatio[i] < 1) {segments(Bmi$mtDNA1[i], Bmi$BMI1[i], Bmi$mtDNA2[i], Bmi$BMI2[i], col = 'red', lwd = 1) }
    }
  if (Bmi$Operation[i] == 'Controls') 
    {
    if (Bmi$Sex[i] == 1) {points(Bmi$mtDNA1[i],Bmi$BMI1[i], pch = 16, col = 'pink')}
    if (Bmi$Sex[i] == 0) {points(Bmi$mtDNA1[i],Bmi$BMI1[i], pch = 16, col = 'blue')}
    }
  }

table(Bmi[!is.na(Bmi$mtDNARatio),]$Sex) # 10 32
table(Bmi[Bmi$mtDNARatio < 1 & !is.na(Bmi$mtDNARatio),]$Sex) # 5 4 => amomng red more males. Out of 10, five behave "unnormally"

cor.test(Bmi[Bmi$Sex == 1 & Bmi$Operation != 'Controls',]$BmiRatio,Bmi[Bmi$Sex == 1 & Bmi$Operation != 'Controls',]$mtDNARatio, method = 'spearman', alternative = 'less') # p = 0.04434, cor = 0.3578
plot(Bmi[Bmi$Sex == 1 & Bmi$Operation != 'Controls',]$BmiRatio,Bmi[Bmi$Sex == 1 & Bmi$Operation != 'Controls',]$mtDNARatio) # p = 0.04434, cor = 0.3578

plot(Bmi[Bmi$Sex == 1 & Bmi$Operation != 'Controls',]$Bmi2Min1,Bmi[Bmi$Sex == 1 & Bmi$Operation != 'Controls',]$mtDNA2Min1) # p = 0.04434, cor = 0.3578
cor.test(Bmi[Bmi$Sex == 1 & Bmi$Operation != 'Controls',]$Bmi2Min1,Bmi[Bmi$Sex == 1 & Bmi$Operation != 'Controls',]$mtDNA2Min1, method = 'spearman')
# dev.off()

summary(Bmi[!is.na(Bmi$mtDNARatio) & Bmi$Sex == 0,]$mtDNARatio)
summary(Bmi[!is.na(Bmi$mtDNARatio) & Bmi$Sex == 1,]$mtDNARatio)
boxplot(Bmi[!is.na(Bmi$mtDNARatio) & Bmi$Sex == 0,]$mtDNARatio,Bmi[!is.na(Bmi$mtDNARatio) & Bmi$Sex == 1,]$mtDNARatio, notch = TRUE, names = c('males','females'), ylab = 'mtDNA2 / mtDNA1')
abline(h = 1, col = 'red', lt = 2)

# max in controls 
MaxMtdnaControls = max(Bmi[Bmi$Operation == 'Controls',]$mtDNA1)
MaxMtdnaControls
boxplot(Bmi[!is.na(Bmi$mtDNARatio) & Bmi$Sex == 0,]$Age,Bmi[!is.na(Bmi$mtDNARatio) & Bmi$Sex == 1,]$Age, notch = TRUE, names = c('males','females'), ylab = 'Age')

dev.off()






