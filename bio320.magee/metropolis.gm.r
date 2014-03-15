#  source("/Users/hongqin/lib/R/metropolis.gm.r")

start = c(0.005, 0.1, 0.008); R = 100; sd=c(0.005, 0.05, 0.005) #for debug
Metropolis <- function(start=IGM, lifespan=lifespan, R = 1000, sd = c(0.005, 0.05, 0.005)) {
    parmcount <- length(start)
    #sims <- matrix(NA, nrow=R, ncol = (parmcount+ 1) )
	sims <- matrix(NA, nrow=R, ncol = parmcount)
	p = rep(0, R)
	colnames(sims) <- names(start)

	oldlogalpha <-  - llh.GM.single.run( start, lifespan)
	#sims[1,] <- c( start, oldlogalpha); 
	sims[1, ] = start; 
	p[1] = oldlogalpha; 
    accepts <- 0

    for (i in 2:R) {
	   newlogalpha = NA; y = NA; 
	   while( is.na(newlogalpha) ) {
		   	jump <- rnorm(parmcount, mean=0, sd=sd)
	   		y <- sims[i-1,] + jump
	   		if (y[1] <=0) { y[1]=10^-15;}
	   		if (y[2] <=0) { y[2]=0.001; }
	   		if (y[3] <=0) { y[3]=0; }
        	newlogalpha <- - llh.GM.single.run( y, lifespan); 
			print( c(y, newlogalpha) )
			# c(y, newlogalpha)
       }
  
	  if (log(runif(1)) < newlogalpha - oldlogalpha) { #why not <= ????
	    #sims[i,] <- c( y, newlogalpha)
		sims[i, ] = y; 
		oldlogalpha <- newlogalpha
	    accepts <- accepts + 1
		p[i] = newlogalpha; 
	  } else {
	    sims[i,] <- sims[i-1,]
		p[i] = p[i-1]; 
	  }
   }
    cat('Accepted ',100*accepts/(R-1),'%\n')
    data.frame( cbind(sims, p) );
}
