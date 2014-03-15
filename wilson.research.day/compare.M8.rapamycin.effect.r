### need to add standard deviation

rm=(list=ls())

tb = read.csv("Feb8,2011.H2O2effectonLOH.M8star.csv");
tb = read.csv("Feb8,2011.H2O2effectonLOH.M8star.csv", colClasses=c("character",NA, NA, "character", rep("numeric",8 ), NA));
names(tb) = c("Strain", "OD600", "Dilution","Date","H2O2stock", "White", "Black", "halfBlack", "quarterBlack", "ThreeQBlack", "QQBlack", "Other", "Notes")
tb$Rapamycin = 0;
###### some manual curations here
tb$Black[25] = NA #This number is not right!!!!

#############rapamycin table 
tbR = read.csv("M8star.Rap.H2O2.031111.csv"); 
tbR = tbR[1:12, 1:14]
#tbR = read.csv("M8star.Rap.H2O2.031111.csv", colClasses=c("character",NA, NA, "character", rep("numeric",8 ), NA));
names(tbR) = c("Strain", "OD600", "Rapamycin", "Date","H2O2stock", "Dilution","White", "Black", "halfBlack", "quarterBlack", "ThreeQBlack", "QQBlack", "Other", "Notes")
tbR$Rapamycin = 5;

myColNames = c("Strain", "OD600", "Rapamycin", "Date","H2O2stock", "Dilution","White", "Black", "halfBlack", "quarterBlack", "ThreeQBlack", "QQBlack", "Other", "Notes")
tb  = tb[ , myColNames]
tbR = tbR[, myColNames]

tb$H2O2 = tb$H2O2stock/2
tb$tot = tb$White + tb$Black + tb$halfBlack + tb$quarterBlack + tb$ThreeQBlack + tb$QQBlack + tb$Other
tb$Dilution = tb$Dilution / tb$Dilution[1]

tbR$H2O2 = tbR$H2O2stock/2
tbR$tot = tbR$White + tbR$Black + tbR$halfBlack + tbR$quarterBlack + tbR$ThreeQBlack + tbR$QQBlack + tbR$Other
tbR$Dilution = tbR$Dilution / tbR$Dilution[1]

tb = rbind( tb, tbR) #this is the merged table 

######## normalize all data
mycolumns = c("White","Black","halfBlack", "quarterBlack","ThreeQBlack", "QQBlack", "Other","tot"); 
for ( j in mycolumns) {
 tb[,j] = tb[,j] * tb$Dilution
}

###### calculate fractions
tbf = tb; 
tbf$s[tbf$Rapamycin==0] = tbf$tot[tbf$Rapamycin==0] / max(tbf$tot[tbf$Rapamycin==0], na.rm=T)
tbf$s[tbf$Rapamycin==5] = tbf$tot[tbf$Rapamycin==5] / max(tbf$tot[tbf$Rapamycin==5], na.rm=T)

for ( j in 7:13) {
  tbf[, j] = tbf[,j] / tbf$tot
}

summary(lm(tbf$Black ~ tbf$s + tbf$H2O2 + tbf$Rapamycin))
summary(lm(tbf$Black ~   tbf$H2O2 + tbf$Rapamycin))
summary(lm( tbf$halfBlack / tbf$Black ~   tbf$H2O2 + tbf$Rapamycin))
summary(lm( tbf$halfBlack / tbf$Black ~   tbf$H2O2  ))
summary(lm( tbf$halfBlack / tbf$Black ~   tbf$Rapamycin  ))

summary(lm(tbf$s ~   tbf$H2O2 + tbf$Rapamycin))
#Rapamycin affects frequency_black but not vaiability!

#do the plot
tb0 = tbf[tbf$Rapamycin==0, ]
tb5 = tbf[tbf$Rapamycin==5, ]

####### find out means, medians
H2O2 = unique( tb0$H2O2)
tbm0 = data.frame(cbind(H2O2))
for ( i in 1:length(H2O2)) {
  c = H2O2[i]
  tmp = tb0[ tb0$H2O2==c, ]	
  tbm0$s[i] = mean(tmp$s, na.rm=T)
  tbm0$s2[i] = median(tmp$s, na.rm=T)
  tbm0$Black[i] = mean(tmp$Black, na.rm=T) 
  tbm0$Black2[i] = mean(tmp$Black, na.rm=T) 
}

####### find out means, medians
H2O2 = unique( tb5$H2O2)
tbm5 = data.frame(cbind(H2O2))
for ( i in 1:length(H2O2)) {
  c = H2O2[i]
  tmp = tb5[ tb5$H2O2==c, ]	
  tbm5$s[i] = mean(tmp$s, na.rm=T)
  tbm5$s2[i] = median(tmp$s, na.rm=T)
  tbm5$Black[i] = mean(tmp$Black, na.rm=T) 
  tbm5$Black2[i] = mean(tmp$Black, na.rm=T) 
}

myxlim= c(0, 0.17)
plot( tbm0$Black ~ tbm0$H2O2, pch=16, xlab='H2O2', ylab='viability', xlim=myxlim)
lines( tbm0$Black ~ tbm0$H2O2, col='black', lty=2)
points( tbm5$Black ~ tbm5$H2O2, col='red', pch=17)
lines( tbm5$Black ~ tbm5$H2O2, col='red', lty=2)

par(new=T)
plot( tbm0$s2 ~ tbm0$H2O2, pch=18, ty='l', lty=2, xlab='H2O2', ylab='viability', axes=F, col='blue', xlim=myxlim )
lines( tbm5$s2 ~ tbm5$H2O2, col='orange',lty=2)

mylables = c('black control', 'viability control', 'black 5ng/ml', 'viabilty 5ng/ml')
mycolors = c('black','blue','red','orange')
legend( 0.12, 0.7, mylables, lty=2, col=mycolors, pch=c(16, NA, 17, NA))


quit("yes")

###### some manual curations here
tbm$halfBlack[tbm$halfBlack==0 & tbm$H2O2>0] = NA;

###### calculate fractions
tbf = tbm; 
tbf$s = tbf$tot / max(tbf$tot)
for ( j in 3:8) {
  tbf[, j] = tbf[,j] / tbf$tot
}









#pdf("M8Sstar.020811.black.pdf", width=5, height=5)
#plot( tbf$s ~ tbf$H2O2, col="blue", axes=F, xlab='', ylab='', ylim=c(0, 1.2), type='p');
#lines( tbf$s ~ tbf$H2O2, col="blue",lty=2)
#axis( 4, at=pretty(c(0, 1.2)))
#legend ( max(H2O2)*0.7, 0.5, c("viability","black"), col=c("blue","black"), lty=c(2,1), pch=c(1,16) )
#par(new=T)
#plot( tbf$Black ~ tbf$H2O2, pch=16, xlab='H2O2',ylab="black")
#lines(tbf$Black ~ tbf$H2O2)
#title("M8 Met15+/- Feb 8, 2011")
#dev.off()

#tbf$H2O2[tbf$H2O2==0] = 1E-2;
#ith( tbf, plot( tot ~ log10(H2O2), col="blue"));
#par(new=T)
#with( tbf, plot( Black ~ log10(H2O2), pch=16))

	
quit("yes")	
   
   