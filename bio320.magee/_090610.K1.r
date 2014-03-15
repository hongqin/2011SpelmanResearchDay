
rm(list=ls());
 
require(nlme)
require(survival)
source("/Users/hongqin/lib/R/lifespan.r")
source("/Users/hongqin/lib/R/metropolis.gm.r")

tb = read.csv( "K_rls_082910.1.csv" );

s1 = calculate.s( tb[,1] )
s = s1$s; t = s1$t
g = gnls( s ~  exp( (I/G) *(1 - exp(G* t)) - M*t ), start=list( I=0.003, G=0.16, M=0.01) ) #for debug 
g1 = gnls( s ~  exp( (I/G) *(1 - exp(G* t)) - M*t ), start=list( I=0.003, G=0.16, M=0.01) )
ret1 = optim ( g1$coeff,  llh.GM.single.run, lifespan=tb[,1], lower = c(1E-10, 1E-5, 0), upper = c(0.2, 2, 0.1), method="L-BFGS-B");
ret1g = optim ( g1$coeff[1:2],  llh.G.single.run, lifespan=tb[,1], lower = c(1E-10, 1E-5), upper = c(0.2, 2 ), method="L-BFGS-B");

s2 = calculate.s( tb[,2] )
s = s2$s; t = s2$t
#g2 = gnls( s ~  exp( (I/G) *(1 - exp(G* t)) - M*t ), start=list( I=0.003, G=0.16, M=0.01) )
g2 = gnls( s ~  exp( (I/G) *(1 - exp(G* t)) ), start=list( I=0.003, G=0.16) )
ret2 = optim ( c(g2$coeff, 0),  llh.GM.single.run, lifespan=tb[,2], lower = c(1E-10, 1E-5, 0), upper = c(0.2, 2, 0.1), method="L-BFGS-B");
#ret2g =optim ( g2$coeff[1:2],llh.G.single.run, lifespan=tb[,2], lower = c(1E-10, 1E-5), upper = c(0.2, 2 ), method="L-BFGS-B");
ret2g =optim ( c(0.00003,0.3),llh.G.single.run, lifespan=tb[,2]);
#Metropolis <- function(start=IGM, lifespan=lifespan, R = 1000, sd = c(0.005, 0.05, 0.005))
ret2metro = Metropolis( c(g2$coeff, 0), lifespan=tb[,2] )
ret2m2 = Metropolis( c(0.00001, 0.3, 0), lifespan=tb[,2] )

#pdf("_K1.rls.083110.pdf")
 plot(   s1$s ~ s1$t, xlim=c(0, max(tb[,1:2], na.rm=T)*1.2 ), ylab="viability", xlab="cell divisions", pch=16 )
 points( s2$s ~ s2$t, col="red", pch=17)
mytexts = names(tb)[1:2]
mycol = c("black", "red", "blue", "green")
legend( 0, 0.3, mytexts, col=mycol, pch=16:19)
#dev.off()

t = seq( 0, 100, by=0.1)
sim.s1 = GM.s ( ret1$par, t)
sim.s1g = GM.s ( c(ret1g$par,0), t)
sim.g1 = GM.s( g1$coeff, t)
lines( sim.s1 ~ t )
lines( sim.s1g ~ t, lty=2 )
lines( sim.g1 ~ t, lty=3, col="green" )

sim.s2 = GM.s ( ret2$par, t)
sim.s2g = GM.s ( c(ret2$par,0), t)
lines( sim.s2 ~ t, col="red")
lines( sim.s2g ~ t, col="red", lty=2)

cbind( ret1$value, ret2$value)
fit = data.frame( cbind( ret1$par, ret2$par ))
names(fit) = names(tb)[1:2]
rownames(fit) = c("I", "G", "M")

#for debug
g1 = gnls( s ~  exp( (I/G) *(1 - exp(G* t)) - M*t ), start=list( I=0.003, G=0.16, M=0.01) )

###loop over every column
#for (j in 3:8 ) {
for ( j in 3:length(tb[1,]) ) {
 first = tb[1,j];
 if ( is.na( first ) ) {
 	tmp = c(NA, NA, NA)
  fit = cbind( fit, tmp );  
 } else {
  sfit = calculate.s( tb[,j] )
  s = sfit$s; t = sfit$t

  g$coeff = c(0.05, 0.1, 0.00); 
  doit = function( ) { g = gnls( s ~  exp( (I/G) *(1 - exp(G* t)) - M*t ), start=list( I=0.003, G=0.16, M=0.0) ) };
  try( doit, silent=FALSE); 
  retJ = optim ( g$coeff,  llh.GM.single.run, lifespan=tb[,j], lower = c(1E-10, 1E-5, 0), upper = c(0.2, 2, 0.1), method="L-BFGS-B");
  #retJ = optim ( c(0.05, 0.1, 0.00),  llh.GM.single.run, lifespan=tb[,j], lower = c(1E-10, 1E-5, 0), upper = c(0.2, 2, 0.1), method="L-BFGS-B");
  fit = cbind( fit, retJ$par );
 }
}
head(tb[1,1:8])

names( fit ) = names(tb)
write.csv( fit, "_gmfit.K.RLS.sheet1.090610.csv")


quit("yes")


### partition and fitting by strains
by4743 = c( tb[,1], tb[,2], tb[,3]);
fit4743 = optim ( c(0.01, 0.2),  llh.gompertz.single.run, lifespan=by4743 );  

by4742 = c( tb[,4], tb[,5], tb[,6]);
fit4742 = optim ( c(0.01, 0.2),  llh.gompertz.single.run, lifespan=by4742 );  

#by4741 = c( tb[,7], tb[,8], tb[,9]);
by4741 = c( tb[,7], tb[,8] );
fit4741 = optim ( c(0.01, 0.2),  llh.gompertz.single.run, lifespan=by4741 );  

by4716 = c( tb[,10]);
fit4716 = optim ( c(0.01, 0.2),  llh.gompertz.single.run, lifespan=by4716 );  


#################################
##### LRT to exam whether two data on I and G.
### 1) model H0, same G and same I
### 2) model H1i, same G, different I
### 3) model H1g, different G, same I
### 4) model H2, different G, different I
                                  #    I1       G1       I2         G2
 H0  <- function( rawIG ) { IG <- c(rawIG[1], rawIG[2], rawIG[1], rawIG[2] ) }  #all the same
 H1i <- function( rawIG ) { IG <- c(rawIG[1], rawIG[2], rawIG[3], rawIG[2] ) }  #different I
 H1g <- function( rawIG ) { IG <- c(rawIG[1], rawIG[2], rawIG[1], rawIG[4] ) }  # different G
 H2  <- function( rawIG ) { IG <- c(rawIG[1], rawIG[2], rawIG[3], rawIG[4] ) }  # all different

 Hx.llh.gompertz.single.run <- function( rawIG, model, lifespan1, lifespan2 ) {
   IG = model(rawIG); 
   I1 = IG[1]; G1 = IG[2]; I2 = IG[3]; G2 = IG[4];
   my.data1 = lifespan1[!is.na(lifespan1)];
   my.data2 = lifespan2[!is.na(lifespan2)];
   log_s1 = (I1/G1) *(1 - exp(G1* my.data1))
   log_s2 = (I2/G2) *(1 - exp(G2* my.data2))
   log_m1 = log(I1) +  G1 * my.data1 ; 
   log_m2 = log(I2) +  G2 * my.data2 ; 
   my.lh = sum(log_s1) + sum(log_m1) + sum(log_s2) + sum(log_m2)
   print (IG ); #trace the convergence
   ret = - my.lh # because optim seems to maximize 
 }

## LRT to exam whether BY4742, BY4741 share the same G, I 
llh.H0  = optim( c(0.01,0.2,0.01,0.1), Hx.llh.gompertz.single.run, model=H0,  lifespan1=by4741, lifespan2=by4742 );
llh.H1i = optim( c(0.01,0.2,0.01,0.1), Hx.llh.gompertz.single.run, model=H1i, lifespan1=by4741, lifespan2=by4742 );
llh.H1g = optim( c(0.01,0.2,0.01,0.1), Hx.llh.gompertz.single.run, model=H1g, lifespan1=by4741, lifespan2=by4742 );
llh.H2  = optim( c(0.01,0.2,0.01,0.1), Hx.llh.gompertz.single.run, model=H2,  lifespan1=by4741, lifespan2=by4742 );

cbind(llh.H0$par, llh.H1i$par, llh.H1g$par, llh.H2$par);
LH = c(llh.H0$value, llh.H1i$value, llh.H1g$value, llh.H2$value);
deltaLH =  - LH + llh.H0$value; 
1 - pchisq( 2*deltaLH, df =1 );

## LRT to exam whether BY4743, BY4742 share the same G, I 
llh.H0  = optim( c(0.01,0.2,0.01,0.1), Hx.llh.gompertz.single.run, model=H0,  lifespan1=by4742, lifespan2=by4743 );
llh.H1i = optim( c(0.01,0.2,0.01,0.1), Hx.llh.gompertz.single.run, model=H1i, lifespan1=by4742, lifespan2=by4743 );
llh.H1g = optim( c(0.01,0.2,0.01,0.1), Hx.llh.gompertz.single.run, model=H1g, lifespan1=by4742, lifespan2=by4743 );
llh.H2  = optim( c(0.01,0.2,0.01,0.1), Hx.llh.gompertz.single.run, model=H2,  lifespan1=by4742, lifespan2=by4743 );

cbind(llh.H0$par, llh.H1i$par, llh.H1g$par, llh.H2$par);

LH = c(llh.H0$value, llh.H1i$value, llh.H1g$value, llh.H2$value);
deltaLH =  - LH + llh.H0$value; 
1 - pchisq( 2*deltaLH, df =1 );


#################################
##### LRT to exam THREE data on I and G.
### 1) model H0, same G and same I
### 2) model H1i, same G, different I
### 3) model H1g, different G, same I

 H0  <- function( rawIG ) {
           #    I1       G1       I2         G2      I3           G3
  IG <- c(rawIG[1], rawIG[2], rawIG[1], rawIG[2], rawIG[1], rawIG[2]) }  #all the same

 H1g3  <- function( rawIG ) {
           #    I1       G1       I2         G2      I3           G3
  IG <- c(rawIG[1], rawIG[2], rawIG[1], rawIG[2], rawIG[1], rawIG[6]) } 

 H1i3  <- function( rawIG ) {
           #    I1       G1       I2         G2      I3           G3
  IG <- c(rawIG[1], rawIG[2], rawIG[1], rawIG[2], rawIG[5], rawIG[2]) } 

 H2ig3  <- function( rawIG ) {
           #    I1       G1       I2         G2      I3           G3
  IG <- c(rawIG[1], rawIG[2], rawIG[1], rawIG[2], rawIG[5], rawIG[6]) } 

 H3i2ig3  <- function( rawIG ) {
           #    I1       G1       I2         G2      I3           G3
  IG <- c(rawIG[1], rawIG[2], rawIG[3], rawIG[2], rawIG[5], rawIG[6]) } 

 H6  <- function( rawIG ) {  IG <- rawIG } 


Hx3.llh.gompertz.single.run <- function( rawIG, model, lifespan1, lifespan2, lifespan3 ) {
   IG = model(rawIG); 
   I1 = IG[1]; G1 = IG[2]; I2 = IG[3]; G2 = IG[4]; I3=IG[5]; G3=IG[6];
   my.data1 = lifespan1[!is.na(lifespan1)];
   my.data2 = lifespan2[!is.na(lifespan2)];
   my.data3 = lifespan3[!is.na(lifespan3)];
   log_s1 = (I1/G1) *(1 - exp(G1* my.data1))
   log_s2 = (I2/G2) *(1 - exp(G2* my.data2))
   log_s3 = (I3/G3) *(1 - exp(G3* my.data3))
   log_m1 = log(I1) +  G1 * my.data1 ; 
   log_m2 = log(I2) +  G2 * my.data2 ; 
   log_m3 = log(I3) +  G3 * my.data3 ; 
   my.lh = sum(log_s1) + sum(log_m1) + sum(log_s2) + sum(log_m2) + sum(log_s3) + sum(log_m3)
   print (IG ); #trace the convergence
   ret = - my.lh # because optim seems to maximize 
 }

llh.H0  = optim( c(0.008,0.06,0.008,0.06, 0.01, 0.08), Hx3.llh.gompertz.single.run, model=H0,  lifespan1=by4741, lifespan2=by4742, lifespan3=by4743 );

llh.H1g3  = optim( c(0.008,0.06,0.008,0.06, 0.01, 0.08), Hx3.llh.gompertz.single.run, model=H1g3,  lifespan1=by4741, lifespan2=by4742, lifespan3=by4743 );

llh.H1i3  = optim( c(0.008,0.06,0.008,0.06, 0.01, 0.08), Hx3.llh.gompertz.single.run, model=H1i3,  lifespan1=by4741, lifespan2=by4742, lifespan3=by4743 );

llh.H2ig3  = optim( c(0.008,0.06,0.008,0.06, 0.01, 0.08), Hx3.llh.gompertz.single.run, model=H2ig3,  lifespan1=by4741, lifespan2=by4742, lifespan3=by4743 );

llh.H3i2ig3  = optim( c(0.008,0.06,0.008,0.06, 0.01, 0.08), Hx3.llh.gompertz.single.run, model=H3i2ig3,  lifespan1=by4741, lifespan2=by4742, lifespan3=by4743 );

llh.H6  = optim( c(0.008,0.06,0.008,0.06, 0.01, 0.08), Hx3.llh.gompertz.single.run, model=H6,  lifespan1=by4741, lifespan2=by4742, lifespan3=by4743 );


cbind(llh.H0$par, llh.H1g3$par, llh.H1i3$par, llh.H2ig3$par, llh.H3i2ig3$par, llh.H6$par);
LH=cbind(llh.H0$value, llh.H1g3$value, llh.H1i3$value, llh.H2ig3$value, llh.H3i2ig3$value, llh.H6$value);
deltaLH =  - LH + llh.H0$value; 
1 - pchisq( 2*deltaLH, df =1 );
1 - pchisq( 2*deltaLH, df =3 );


#################################
###### End of LRT    ############
#################################

