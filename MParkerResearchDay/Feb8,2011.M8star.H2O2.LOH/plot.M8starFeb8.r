### need to add standard deviation

rm=(list=ls())

tb = read.csv("Feb8,2011.H2O2effectonLOH.M8star.csv");
tb = read.csv("Feb8,2011.H2O2effectonLOH.M8star.csv", colClasses=c("character",NA, NA, "character", rep("numeric",8 ), NA));
names(tb) = c("Strain", "OD600", "Dilution","Date","H2O2stock", "White", "Black", "halfBlack", "quarterBlack", "ThreeQBlack", "QQBlack", "Other", "Notes")

tb$H2O2 = tb$H2O2stock/2
tb$tot = tb$White + tb$Black + tb$halfBlack + tb$quarterBlack + tb$ThreeQBlack + tb$QQBlack + tb$Other
tb.ori = tb; 

tb$Dilution = tb$Dilution / tb$Dilution[1]

###### some manual curations here
#tb$Dilution[c(19,20,21)] = c(10,10,10)  ##this dilution should be 10 times more
tb$Black[25] = NA #This number is not right!!!!

######## normalize all data
mycolumns = c("White","Black","halfBlack", "quarterBlack","ThreeQBlack", "QQBlack", "Other","tot"); 
for ( j in mycolumns) {
 tb[,j] = tb[,j] * tb$Dilution
}

####### find out means
H2O2 = unique( tb$H2O2)
#s = H2O2
tbm = data.frame(cbind(H2O2))
for ( i in 1:length(H2O2)) {
  c = H2O2[i]
  tmp = tb[ tb$H2O2==c, ]	
  tbm$tot[i] = mean(tmp$tot, na.rm=T)
  tbm$White[i] = mean(tmp$White, na.rm=T)
  tbm$Black[i] = mean(tmp$Black, na.rm=T) 
  tbm$halfBlack[i] = mean(tmp$halfBlack, na.rm=T)  
  tbm$quarterBlack[i] = mean(tmp$quarterBlack, na.rm=T)
  tbm$ThreeQBlack[i] = mean(tmp$ThreeQBlack, na.rm=T)  
  tbm$QQBlack[i] = mean(tmp$QQBlack, na.rm=T);  
}

###### some manual curations here
tbm$halfBlack[tbm$halfBlack==0 & tbm$H2O2>0] = NA;

###### calculate fractions
tbf = tbm; 
tbf$s = tbf$tot / max(tbf$tot)
for ( j in 3:8) {
  tbf[, j] = tbf[,j] / tbf$tot
}

pdf("M8Sstar.020811.black.pdf", width=5, height=5)
plot( tbf$s ~ tbf$H2O2, col="blue", axes=F, xlab='', ylab='', ylim=c(0, 1.2), type='p');
lines( tbf$s ~ tbf$H2O2, col="blue",lty=2)
axis( 4, at=pretty(c(0, 1.2)))
legend ( max(H2O2)*0.7, 0.5, c("viability","black"), col=c("blue","black"), lty=c(2,1), pch=c(1,16) )
par(new=T)
plot( tbf$Black ~ tbf$H2O2, pch=16, xlab='H2O2',ylab="black")
lines(tbf$Black ~ tbf$H2O2)
title("M8 Met15+/- Feb 8, 2011")
dev.off()

tbf$H2O2[tbf$H2O2==0] = 1E-2;
with( tbf, plot( tot ~ log10(H2O2), col="blue"));
par(new=T)
with( tbf, plot( Black ~ log10(H2O2), pch=16))

	
quit("yes")	
   
   