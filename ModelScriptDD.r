# Script for running delay-difference population model for scallop stocks.
#	The following is an example of a few formulations using data from scallop 
#   stocks in Golden Bay, New Zealand and Georges Bank, Canada.
#	Script reads in data, sets up priors, call WinBUGS to fit the model and 
#   runs a few plotting functions in R 


	# Data Inputs (DDNZ.dat)
	#
	# C = Catch
	# I = Biomass index for commercial size from survey
	# IR = Biomass index for reruit size from survey
	# w.bar = average weight of commercial size scallop
	# w.k = average weight of recruit size scallop

	# Other potential inputs
	#
	# U = Commercial CPUE index
	# I.cv, IR.cv, U.cv = Coefficients of variation for biomass indices
	# CF = annual estimates of average condition factor for the stock
	# l.bar = average shell width of commercial size scallop
	# l.k = average shell width of recruit size scallop
	# N = total numbers of commercial size scallops from survey
	# NR = total numbers of recruit size scallops from survey
	# clappers = total numbers of dead (cluckers) commercial size scallops from survey
	# clappersR =  total numbers of dead (cluckers) recruit size scallops from survey

## R packages
install.packages(c("R2WinBUGS"))

MyDirectory <- 'C:/Users/owner/Dropbox/Public/flatfish' #change to your directory
setwd(MyDirectory) 

# New Zealand scallop data from Golden Bay
SCA7data <- read.csv('data/NZscallopSCA7.csv', skip = 7, header = T, sep = ',')
head(SCA7data)
#these data assume recruit sized scallops will recruit to the commercial size (90 mm) the following year
#the data have been scaled by sampled fraction, area swept by the dredge, and stratum area (but not dredge efficiency).
#instead, catchability q is estimated in the model.


# Georges Bank scallop data 
GBdata <- read.csv('data/GeorgesBankScallop.csv')
head(GBdata)
#these data assume recruit sized scallops will recruit to the commercial size (95 mm) the following year
#the data have been scaled by sampled fraction, area swept by the dredge, and stratum area (but not dredge efficiency).
#instead, catchability q is estimated in the model.


###### Basic Version ######

	DDNZ.dat<-with(subset(SCA7data, Region == 'GB'),list(C=C,I=I,IR=IR,w.bar=w.bar,w.k=w.k))
	
		yrs<-1997:2012	# years of data
		NY<- length(yrs)
		DDNZ.dat$NY<-NY  # number of years


	## Growth

		# predict weight at age from LVB length at age

		Linf <- 144
		K <- 0.4
		t0 <- 0

		# back calculate 1 year to calculate recruit size 
		source("_Rfunctions/age.back.r")
		age.back(90,LVB=list(linf=Linf,k=K,t0=t0),age.dif=-1)

		LWa <- 0.00037
		LWb <- 2.69
		LWa * 100^LWb

		lvb <- Linf * (1 - exp(-K * (1:12 - t0)))
		waa.t <- lvb^LWb * LWa

		# fit weight at age model to calculate rho and alpha
		waa.tm1 <- c(NA, waa.t)
		waa.t <- c(waa.t, NA)
		waa.lm <- lm(waa.t ~ waa.tm1)
		alpha <- coef(waa.lm)[1]
		rho <- coef(waa.lm)[2]

		# Use rho and alpha to calculate Growth terms
		DDNZ.dat$g<-rho + alpha / DDNZ.dat$w.bar # biomass growth multiplier for commercial size
		DDNZ.dat$gR<-rho + alpha / DDNZ.dat$w.k # biomass growth multiplier for recruit size

		# Growth if you have CF, l.bar and l.k
		#DDNZ.dat$g<-getG(CF,l.bar,b=LWb,VB=c(144,0.4,0))$g
		#DDNZ.dat$gR<-getG(CF,l.k,b=LWb,VB=c(144,0.4,0))$g


		DDNZ.dat<-with(DDNZ.dat, list(C=C,I=I,IR=IR,g=g,gR=gR,NY=NY)) # it's important you don't give WinBUGS any information it doesn't need (w.bar, w.k)

	## Priors

		DDNZpriors=list(
			logK=			list(a=5,		b=5,		d="dnorm",		i1=7,	i2=5,	l=1		),		# scaler to total biomass
			r=				list(a=0, 		b=1,		d="dlnorm",		i1=0.2,	i2=0.9,	l=NY	),		# recruit index
			m=				list(a=-1.6, 	b=1,		d="dlnorm",		i1=0.2,	i2=0.3,	l=NY	),		# mortality commercial size
			mR=				list(a=-1.6, 	b=1,		d="dlnorm",		i1=0.2,	i2=0.3,	l=NY	),		# mortality recruits
			q=				list(a=40, 		b=40, 		d="dbeta",		i1=0.2,	i2=0.5,	l=1		),		# catchability for survey commercial size
			qR=				list(a=40, 		b=40, 		d="dbeta",		i1=0.2,	i2=0.5,	l=1		),		# catchability for survey  recruits
			sigma=			list(a=0, 		b=5,		d="dunif",		i1=2,	i2=3,	l=1		),		# process error (SD)
			itau2=			list(a=3, 		b=0.44629,	d="dgamma",		i1=15,	i2=30,	l=1		),		# measurement error for commercial size (precision)
			iepsilon2=		list(a=3, 		b=0.44629,	d="dgamma",		i1=15,	i2=30,	l=1		)		# measurement error for recruits (precision)
		)

		# Check priors by plotting the distributions in R
		#plot(seq(0,5,l=100),dlnorm(seq(0,5,l=100),-1.6,1.4),type='l') 		# m
		#plot(seq(0,1,l=100),dbeta(seq(0,1,l=100),40,40),type='l') 			# q
		#lines(seq(0,1,l=100),dbeta(seq(0,1,l=100),7,15),col="red") 		# qR

	## Runs WinBUGS
	source("_Rfunctions/delayBUGS.r")
	DDNZ.out<-delayBUGS("DDNZ", DDNZ.dat, DDNZpriors, yrs, n = 60000, burn = 30000, thin = 10,debug=F,parameters=c(names(DDNZpriors),'K','P','B','R','mu','Iresid','IRresid','Presid'),wd=MyDirectory,sw='jags')

	## Saves Model Output
	save(list=c("DDNZ.out","DDNZ.dat"),file="DDNZ.Rdata")
	write.csv(DDNZ.out$summary,"PostSumNZGB.csv")

	## Plotting model results

		# plot fits to abundance indices 
		source("_Rfunctions/fit.plt.R")
		fit.plt(DDNZ.out, CI=T,graphic='R')

		# plot the posterior distributions of the estimated parameters
		source("_Rfunctions/post.plt.R")
		post.plt(DDNZ.out,DDNZpriors,years=yrs, graphic='R',nr=3,nc=3,wd=15,multi=F)
		post.plt(DDNZ.out,DDNZpriors,years=yrs, graphic='R',nr=2,nc=3,wd=15,multi=T)

		# plot the biomass estimates for commercial and recruit size scallops
		source("_Rfunctions/biomass.plt.r")
		biomass.plt(DDNZ.out,years=yrs, graphic='R')	

		# plot the expliotion rate and natural survival fraction
		source("_Rfunctions/exploit.plt.r")	
		exploit.plt(DDNZ.out, years=yrs, plt=c('f','m','mR'),graphic='R')

		# plot residuals for the fit and process (process residuals represent the difference between what the dynamics and the data say about biomass)
		# Note: with current version of the data there are some large process reiduals 
		source("_Rfunctions/diag.plt.R")
		diag.plt(DDNZ.out, yrs,graphic='R')

		# add projected biomass to the model output object for various catch senarios (C.p)
		#source("_Rfunctions/projections.r")
		#DDNZ.out<-projections(DDNZ.out,C.p=seq(20,300,20))

		# preform the prediction evaluation procedure:
		# runs model up to years in pe predicting biomass in the nnext year then compares that prediction 
		# to the estimate when the model is fit up to that year
		#source("_Rfunctions/peR.r")
		#peR("DDNZ", DDNZ.dat, DDNZpriors, yrs, pe=2012:2004,n = 60000, burn = 30000, thin = 10, plot=0,lab='NZGB1',debug=F,wd=MyDirectory) # models running
		#peR("DDNZ", DDNZ.dat, DDNZpriors, yrs, pe=2012:2004,run=F, plot=3,graphic='R',lab='NZGB1') # plots results
		

###### Hyperprior m Version ######

	# This version is identical to the basic vesion only that we let WinBUGS estimate the prior on m and mR
	# 
	# Notes:
	# 1. By letting WinBUGS estimate the prior on m you are giving the prior on q more influence because you are giving it less information about m
	# 2. Has trouble estimating both m and mR, changed to try with one m
	# 3. Estimated very high m with current version of the data but fixed the process residual pattern

	DDNZhm.dat<-with(subset(SCA7data, Region == 'GB'),list(C=C,I=I,IR=IR,w.bar=w.bar,w.k=w.k))
	
		yrs<-1997:2012	# years of data
		NY<- length(yrs)
		DDNZhm.dat$NY<-NY  # number of years


	## Growth

		# predict weight at age from LVB length at age
		Linf <- 144
		K <- 0.4
		t0 <- 0

		# back calculate 1 year to calculate recruit size [these parameters produce a range from 65-90 which we think is way too large]
		source("_Rfunctions/age.back.r")
		age.back(90,LVB=list(linf=Linf,k=K,t0=t0),age.dif=-1)

		LWa <- 0.00037
		LWb <- 2.69
		LWa * 100^LWb

		lvb <- Linf * (1 - exp(-K * (1:12 - t0)))
		waa.t <- lvb^LWb * LWa

		# fit weight at age model to calculate rho and alpha
		waa.tm1 <- c(NA, waa.t)
		waa.t <- c(waa.t, NA)
		waa.lm <- lm(waa.t ~ waa.tm1)
		alpha <- coef(waa.lm)[1]
		rho <- coef(waa.lm)[2]

		# Use rho and alpha to calculate Growth terms
		DDNZhm.dat$g<-rho + alpha / DDNZhm.dat$w.bar # biomass growth multiplier for commercial size
		DDNZhm.dat$gR<-rho + alpha / DDNZhm.dat$w.k # biomass growth multiplier for recruit size

		DDNZhm.dat<-with(DDNZhm.dat, list(C=C,I=I,IR=IR,g=g,gR=gR,NY=NY)) # it's important you don't give WinBUGS any information it doesn't need (w.bar, w.k)


	## Priors

		DDNZhmpriors=list(
			logK=			list(a=5,		b=5,		d="dnorm",		i1=7,	i2=5,	l=1		),		# scaler to total biomass
			r=				list(a=0, 		b=1,		d="dlnorm",		i1=0.2,	i2=0.9,	l=NY	),		# recruit index
			m.u=			list(a=0.01, 	b=5,		d="dunif",		i1=0.2,	i2=0.3,	l=1		),		# mean mortality commercial size
			#mR.u=			list(a=0.01, 	b=5,		d="dunif",		i1=0.2,	i2=0.3,	l=1		),		# mean mortality recruits
			m.sd=			list(a=0.01, 	b=5,		d="dunif",		i1=0.2,	i2=0.3,	l=1		),		# sd mortality commercial size
			#mR.sd=			list(a=0.01, 	b=5,		d="dunif",		i1=0.2,	i2=0.3,	l=1		),		# sd mortality recruits
			q=				list(a=40, 		b=40, 		d="dbeta",		i1=0.2,	i2=0.5,	l=1		),		# catchability for survey commercial size
			qR=				list(a=40, 		b=40, 		d="dbeta",		i1=0.2,	i2=0.5,	l=1		),		# catchability for survey  recruits
			sigma=			list(a=0, 		b=5,		d="dunif",		i1=2,	i2=3,	l=1		),		# process error (SD)
			itau2=			list(a=3, 		b=0.44629,	d="dgamma",		i1=15,	i2=30,	l=1		),		# measurement error for commercial size (precision)
			iepsilon2=		list(a=3, 		b=0.44629,	d="dgamma",		i1=15,	i2=30,	l=1		)		# measurement error for recruits (precision)
		)

		# Check priors by plotting the distributions in R
		#plot(seq(0,20,l=100),dgamma(seq(0,20,l=100),3,0.44629),type='l') 		# itau2

	## Runs WinBUGS
	source("_Rfunctions/delayBUGS.r")
	DDNZhm.out<-delayBUGS("DDNZhyperM", DDNZhm.dat, DDNZhmpriors, yrs, n = 60000, burn = 30000, thin = 10,debug=F,parameters=c(names(DDNZhmpriors),'K','P','B','R','mu','m','Iresid','IRresid','Presid'),wd=MyDirectory)

	## Saves Model Output
	save(list=c("DDNZhm.out","DDNZhm.dat"),file="DDNZhm.Rdata")
	write.csv(DDNZhm.out$summary,"PostSumNZGB.csv")

	## Plotting model results

		# plot fits to abundance indices 
		source("_Rfunctions/fit.plt.R")
		fit.plt(DDNZhm.out, CI=T,graphic='R')

		# plot the posterior distributions of the estimated parameters
		source("_Rfunctions/post.plt.R")
		post.plt(DDNZhm.out,DDNZhmpriors,years=yrs, graphic='R',nr=3,nc=4,wd=20,multi=F)
		post.plt(DDNZhm.out,DDNZhmpriors,years=yrs, graphic='R',nr=2,nc=3,wd=15,multi=T)

		# plot the biomass estimates for commercial and recruit size scallops
		source("_Rfunctions/biomass.plt.r")
		biomass.plt(DDNZhm.out,years=yrs, graphic='R')	

		# plot the expliotion rate and natural survival fraction
		source("_Rfunctions/exploit.plt.r")	
		exploit.plt(DDNZhm.out, years=yrs, plt=c('f','m'),graphic='R')

		# plot residuals for the fit and process (process residuals represent the difference between what the dynamics and the data say about biomass)
		source("_Rfunctions/diag.plt.R")
		diag.plt(DDNZhm.out, yrs,graphic='R')

		# add projected biomass to the model output object for various catch senarios (C.p)
		#source("_Rfunctions/projections.r")
		#DDNZhm.out<-projections(DDNZhm.out,C.p=seq(20,300,20))

		# preform the prediction evaluation procedure:
		# runs model up to years in pe predicting biomass in the nnext year then compares that prediction 
		# to the estimate when the model is fit up to that year
		#source("_Rfunctions/peR.r")
		#peR("DDNZhyperM", DDNZhm.dat, DDNZhmpriors, yrs, pe=2012:2004,n = 60000, burn = 30000, thin = 10, plot=0,lab='NZGB1',debug=F,wd=MyDirectory) # models running
		#peR("DDNZhyperM", DDNZhm.dat, DDNZhmpriors, yrs, pe=2012:2004,run=F, plot=3,graphic='R',lab='NZGB1') # plots results
		


###### CV Version using Georges Bank data ######

	# This version incorporates the CVs from the abundance indices to inform (with a prior) the estimated observation error
	# 
	# Notes:
	# 1. 
	# 2. 
	# 3. 

		yrs<-1986:2013	# years of data
		DDGB.dat<-as.list(subset(GBdata,year%in%yrs,c("I","I.cv","IR","IR.cv","g","gR","C","U","U.cv","N","NR","clappers","clappersR")))
	
		NY<- length(yrs)
		DDGB.dat$NY<-NY  # number of years

		# Growth (how we get growth from mean length and condition factor)
		source("_Rfunctions/getG.r")
		getG(GBdata$CF,GBdata$l.bar,b=3,VB=c(149, 0.22, 0.22))

		
		# Set up Priors
		uI=log(DDGB.dat$I.cv^2+1)
		Ip.a=2+(uI/uI)^2
		Ip.b=uI*((uI/uI)^2+1)
		
		uIR=log(DDGB.dat$IR.cv^2+1)
		IRp.a=2+(uIR/uIR)^2
		IRp.b=uIR*((uIR/uIR)^2+1)
		
		uU=log(DDGB.dat$U.cv^2+1)
		Up.a=2+(uU/uU)^2
		Up.b=uU*((uU/uU)^2+1)
		
		DDGBpriors=list(
				logK=			list(a=7,		b=7,		d="dnorm",		i1=8,	i2=10,	l=1		),		# scaler to total biomass
				r=				list(a=0, 		b=1,		d="dlnorm",		i1=0.3,	i2=0.9,	l=NY	),		# recruit index
				q=				list(a=20, 		b=40, 		d="dbeta",		i1=0.2,	i2=0.5,	l=1		),		# catchability for survey fully recruited
				qU=				list(a=0, 		b=1,		d="dunif",		i1=0.4,	i2=0.7,	l=1		),		# catchability for commercial CPUE
				m=				list(a=-2,		b=2,		d="dlnorm",		i1=0.1,	i2=0.3,	l=NY	),		# natural mortality fully recruited
				mR=				list(a=-2,		b=2,		d="dlnorm",		i1=0.2,	i2=0.4,	l=NY	),		# natural mortality  recruits
				S=				list(a=8, 		b=11, 		d="dbeta",		i1=0.5,	i2=0.8,	l=1		),		# clapper dissolution rate
				sigma=			list(a=0, 		b=5,		d="dunif",		i1=2,	i2=3,	l=1		),		# process error (SD)
				sigma.tau=		list(a=2, 		b=1,		d="dgamma",		i1=15,	i2=30,	l=1		),		# measurement error for fully recruited from survey (precision)
				sigma.rho=		list(a=2, 		b=1,		d="dgamma",		i1=15,	i2=30,	l=1		),		# measurement error for recruits from survey (precision)
				sigma.upsilon=	list(a=2, 		b=1,		d="dgamma",		i1=15,	i2=30,	l=1		),		# measurement error for recruits from survey (precision)
				ikappa.tau2=	list(a=3, 		b=0.44629,	d="dgamma",		i1=15,	i2=30,	l=1		),		# measurement error for fully recruited clappers (precision)
				ikappa.rho2=	list(a=3, 		b=0.44629,	d="dgamma",		i1=15,	i2=30,	l=1		),		# measurement error for recruit clappers (precision)
				I.precision=	list(a=Ip.a,	b=Ip.b,		d="dgamma",		i1=15,	i2=30,	l=NY	),		# measurement error for fully recruited from survey (precision)
				IR.precision=	list(a=IRp.a,	b=IRp.b,	d="dgamma",		i1=15,	i2=30,	l=NY	),		# measurement error for recruits from survey (precision)
				U.precision=	list(a=Up.a,	b=Up.b,		d="dgamma",		i1=15,	i2=30,	l=NY	)		# measurement error for recruits from survey (precision)
			)
				
		
		source("_Rfunctions/delayBUGS.r")
		DDGB.out<-delayBUGS("DDGBcv", DDGB.dat, DDGBpriors, 1986:2013, n = 300000, burn =200000, thin = 10,debug=F,add.parameters=c('Imed','Ipred','Irep','IRmed','IRpred','IRrep','sIresid','sIRresid','sPresid','Iresid','IRresid','Presid'),wd=MyDirectory)
	
		save(list=c("DDGB.out","DDGB.dat"),file="DDGB.RData")
		load("DDGB.RData")
		plotsGo<-file.path(MyDirectory,'figures')
		#write.csv(DDGB.out$summary,paste0(plotsGo,"PostSumGBa.csv"))
	
	 	
		source("_Rfunctions/post.plt.R")
		post.plt(DDGB.out,DDGBpriors,years=1986:2013, graphic='pdf',nr=3,nc=4,wd=15,multi=F,path=plotsGo)
		post.plt(DDGB.out,DDGBpriors,years=1986:2013, graphic='pdf',nr=3,nc=4,wd=15,multi=T,path=plotsGo)
	
		source("_Rfunctions/exploit.plt.r")	
		exploit.plt(DDGB.out, years=1986:2013, plt=c('f','m','mR'),graphic='pdf',path=plotsGo)
	
		source("_Rfunctions/fit.plt.R")
		fit.plt(DDGB.out, CI=T,graphic='pdf',path=plotsGo,CV=T)
			
		source("_Rfunctions/diag.plt.R")
		diag.plt(DDGB.out, 1986:2013,graphic='pdf',path=plotsGo)
		
	 	source("_Rfunctions/peR.r")
	 	peR("DDGBcv", DDGB.dat, DDGBpriors, 1986:2013, pe=2013:2000,n = 60000, burn = 30000, thin = 10, plot=0,lab='GBag1',wd=MyDirectory)
	 	peR("DDGBcv", DDGB.dat, DDGBpriors, 1986:2013, pe=2013:2001,run=F, plot=3,proj.from="BUGS",graphic='R',lab='GBag1',path=plotsGo)
			

			
	source("_Rfunctions/projections.r")
	DDGB.out<-projections(DDGB.out,C.p=seq(2000,8000,500)+530)
	source("_Rfunctions/biomass.plt.r")
	biomass.plt(DDGB.out,years=1986:2013, graphic='pdf',TAC=5000+530,path=plotsGo,index=1:24,avg.line=median)	
	 	
		
	source("_Rfunctions/decision.r")
	DtabGBa<-decision("GBa",DDGB.out, mu=0.15,refs=c(12789,4796),post.survey.C=530)
 	#write.csv(DtabGBa,file.path(plotsGo,"Decision1_GB.csv"),row.names=F)

 