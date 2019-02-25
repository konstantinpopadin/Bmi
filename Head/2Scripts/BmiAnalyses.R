###################################
###################################

rm(list=ls(all=TRUE))

Bmi = read.table("../../Body/1Raw/MainTable.txt", header = TRUE, sep = '\t')
Bmi$Sex = Bmi$Sex-1; table(Bmi$Sex) # make it Dummy: 0 - males, 1 - females
Bmi$BmiRatio = Bmi$BMI2/Bmi$BMI1
Bmi$mtDNARatio = Bmi$mtDNA2/Bmi$mtDNA1

###### BMI decreased after operation (both types of operations, both genders)
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


####### mtDNA increased after operation (both types of operations)
summary(Bmi$mtDNA1)
summary(Bmi$mtDNA2)
wilcox.test(Bmi$mtDNA1,Bmi$mtDNA2, paired = TRUE)
wilcox.test(Bmi[Bmi$Operation == 'LSG',]$mtDNA1,Bmi[Bmi$Operation == 'LSG',]$mtDNA2, paired = TRUE)  # 0.011
wilcox.test(Bmi[Bmi$Operation == 'RYGB',]$mtDNA1,Bmi[Bmi$Operation == 'RYGB',]$mtDNA2, paired = TRUE)  # 0.0058

######## pdf

pdf("../../Body/4Figures/Bmi.mtDNA.pdf")

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
    if (Bmi$mtDNARatio[i] >= 1) {segments(Bmi$mtDNA1[i], Bmi$BMI1[i], Bmi$mtDNA2[i], Bmi$BMI2[i], col = 'grey', lwd = 1) }
    if  (Bmi$mtDNARatio[i] < 1) {segments(Bmi$mtDNA1[i], Bmi$BMI1[i], Bmi$mtDNA2[i], Bmi$BMI2[i], col = 'red', lwd = 1) }
    }
  if (Bmi$Operation[i] == 'Controls') 
    {
    points(Bmi$mtDNA1[i],Bmi$BMI1[i], pch = 16, col = 'black')
    }
  }

table(Bmi[!is.na(Bmi$mtDNARatio),]$Sex) # 10 32
table(Bmi[Bmi$mtDNARatio < 1 & !is.na(Bmi$mtDNARatio),]$Sex) # 5 4 => amomng red more males. Out of 10, five behave "unnormally"

summary(Bmi[!is.na(Bmi$mtDNARatio) & Bmi$Sex == 0,]$mtDNARatio)
summary(Bmi[!is.na(Bmi$mtDNARatio) & Bmi$Sex == 1,]$mtDNARatio)
boxplot(Bmi[!is.na(Bmi$mtDNARatio) & Bmi$Sex == 0,]$mtDNARatio,Bmi[!is.na(Bmi$mtDNARatio) & Bmi$Sex == 1,]$mtDNARatio, notch = TRUE, names = c('males','females'), ylab = 'mtDNA2 / mtDNA1')
abline(h = 1, col = 'red', lt = 2)


# max in controls 
MaxMtdnaControls = max(Bmi[Bmi$Operation == 'Controls',]$mtDNA1)
MaxMtdnaControls
boxplot(Bmi[!is.na(Bmi$mtDNARatio) & Bmi$Sex == 0,]$Age,Bmi[!is.na(Bmi$mtDNARatio) & Bmi$Sex == 1,]$Age, notch = TRUE, names = c('males','females'), ylab = 'Age')

dev.off()






