


require(nloptr); require(deSolve); require(abind);


#Implement variation in recruitment
qlnormT=function(p,mu,sig) qlnorm(p, meanlog=log(mu/sqrt(1+(sig/mu)^2)), sdlog=sqrt(log(1+(sig/mu)^2)) )
#Indexing function
codify=function(x,cols=1:ncol(x),sep="_") as.matrix(cbind(Index=apply(x[,cols],1,paste0,collapse=sep),x[,-cols]))


#Simulation results processing
errfun=function(x,lout) tryCatch(x, error=function(e) rep(NA,lout))
ODEtrim=function(x,Times,xp=pmax(x[,-1],0)){ if(nrow(xp)==max(Times)) return(xp[Times,]) else return(tail(xp,2)[pmin(Times,1),]); } #If ended prematurely due to explosion, return 2nd to last sln

#Main model
siteODE3=function(t,start,parms,VAR,IDs){
	Ns=NsF=pmax(start,1e-6); mD=rep(unlist(lapply(split(Ns,parms$strat[IDs[[2]]]),sum)), parms$stratRuns[IDs[[1]]])/parms$runs[IDs[[2]]];
	if(parms$sitefeed==1) NsF=mD; gamFeed=parms$gam[IDs[[2]]]*(1-parms$b+parms$b/(1+parms$fg*NsF^GrazingInhibExponent));
	return(list(    parms$grow[IDs[[2]]]*(1+(parms$vartype==2)*(VAR-1))*(Ns*(1-parms$Rr) + mD*parms$Rr)*(parms$Gs[IDs[[2]]]-parms$d*Ns)/(parms$dR*gamFeed+1) - Ns*(parms$mu[IDs[[2]]]+gamFeed)   ))
}

#Set up simulations and implement solver in batches for faster speed
rtesVar4=function(parms,covars,Vars2,ICs=nrow(covars),strat,sitefeed=FALSE,vartype=2,batchsize=1e2,Time=1){
	Nobs=length(strat); stratRuns=unlist(lapply(split(rep(1,length(strat)),strat),sum));
	PARMS=list(strat=strat, stratRuns=stratRuns, d=parms["d"], sitefeed=sitefeed, vartype=vartype,
		fg=parms["fg"]*(parms["fg"]!=1e-5), b=1, runs=rep(stratRuns,stratRuns), Gs=covars[,"L"], Rr=parms["Rr"], dR=parms["RdeltaA"],
		grow = parms["r"]*covars[,"N"], mu = parms["mu"]*covars[,"E"],
		gam = pmax(parms["deltaA"]*covars[,"U"]*exp(-parms["deltaU"]*covars[,"P"]),1e-4) )
	#For TS fit, ICmat assumes ICs already sorted as relevant and absent ICs omitted in all inputs
	if(length(ICs)==1) ICmat=matrix(seq(LNthresh,0.8*UNthresh,length=ICs),Nobs,ICs,byrow=TRUE) else ICmat=t(t(ICs)) 
	ICids=rep(1:ncol(ICmat),length(Vars2)/ncol(ICmat)); NTS=matrix(NA,Nobs,length(Vars2)); #Columns in NTS will loop over all ICs, then over each Var level
	stratL=split(1:Nobs,strat); if(length(stratL)==1) batches=matrix(1,1,1) else batches=suppressWarnings(matrix(1:length(stratL),ncol=min(batchsize/mean(stratRuns),length(stratL)),byrow=TRUE));
	for(n in 1:nrow(batches)){ IDs=list(batches[n,],unlist(stratL[batches[n,]])); 
		for(i in 1:length(Vars2))	NTS[IDs[[2]],i]=errfun(as.vector(t(ODEtrim(lsoda(ICmat[IDs[[2]],ICids[i]],0:Time,siteODE3,parms=PARMS,VAR=rep(Vars2[i],length(IDs[[2]])),IDs=IDs),c(Time+(length(ICs)>0),0)))),length(NTS[IDs[[2]],i])); 
		}
	# if(length(ICs)>1) return(NTS); assign("Nsims",c(Nsims,mean(abs(NTS[,c(TRUE,FALSE)]-NTS[,c(FALSE,TRUE)])>0.05)),.GlobalEnv)
	return(cbind(NTS[,c(TRUE,FALSE)],NTS[,c(FALSE,TRUE)]));
}

#evaluate model results
bests=function(candsFull,targ=NA,strat=NA,NC=ncol(candsFull)/2,COR=TRUE){
	cands=cbind(rowMeans(candsFull[,1:NC]),rowMeans(candsFull[,-(1:NC)]))
	if(is.na(targ[1])) mtch=apply(cands,1,which.max) else mtch=apply(abs(cbind(cands[,1]-targ,cands[,2]-targ)),1,which.min); 
	if(!is.na(strat[1])) mtch=round(unlist(lapply(split(mtch,strat),function(x) rep(mean(x),length(x)))))
	if(is.na(targ[1])) return(t(apply(cbind(mtch,candsFull),1,function(x){ if(x[1]==1) x[2:(NC+1)]; x[-(1:(NC+1))]; })))
	if(COR) return(cor(targ,jitter(cands[cbind(1:nrow(cands),mtch)],1e-5))^2);
	if(!COR) cbind(targ,cands[cbind(1:nrow(cands),mtch)])
}

nrm=function(x,MinScl=min(x,na.rm=TRUE)) (x - MinScl)/(max(x,na.rm=TRUE) - MinScl)
cnt=function(x) x/mean(x,na.rm=TRUE)

#Main likelihood function
LoglikFun=function(parmsTry,parmsFit=parmnames,extras=c(LLvrs=1,plotgive=0,allkelp=0,model=2,waves=0,vartype=1,Sims=0),reps=10){
	print(c(round(parmsTry,5),prevBst=REPORT[which.min(REPORT[,14]),14:17]))
	parms[parmsFit]=parmsTry; fg0=parms["fg"]; parms["fg"]=exp(fg0); parms[!(parmnames%in%parmsFit)]=1e-5; 
	parms[c("fw","deltaU","Rr","b")]=pmax(parms[c("fw","deltaU","Rr","b")],1e-5); parms[c("Rr","b")]=pmin(parms[c("Rr","b")],1); 
	parms[c("b","sigma")]=c(1,0.35*(extras["Sims"]!=999)); if("RV"%in%ls(envir=Environment)) parms["sigma"]=RV; #Constraining parameters fitted in earlier model versions
	
	CovarNms=c("N","L","E","U","R","P"); DatTofit=DATA; DatTofit[,"R"]=0; DatTofit[,"U"]=DatTofit[,"Uq"];
	Erels=cnt(DatTofit[,"E"]); DatTofit[,"E"]=cnt(pmax(Erels*exp(-DatTofit[,"z"]*parms["fw"]/Erels),0.001));
	DatTofit[,"P"]=DatTofit[,"P"]/max(DatTofit[,"P"],na.rm=TRUE) #Standardize predator density across region
	DatTofit[,"L"]=1+parms["Lgrow"]*DatTofit[,"L"]
	Ns=0.98*nrm(DatTofit[,"N"])+0.02; DatTofit[,"N"]=1;
	if("Ngrow"%in%parmsFit) DatTofit[,"N"]=Ns*(exp(parms["Ngrow"]) + 1) / (1 + exp(parms["Ngrow"])*Ns); #Rescaling nitrate limitation recruitment modifier to be between 0 and 1.	
	datTofit=DatTofit[complete.cases(DatTofit[,c("Obs",CovarNms)]),]
	
	NP=1e2; Vars=c(seed=112,simt=8,simICs=2+1*as.numeric(extras["LLvrs"]==1.1),extras["vartype"],sitefeed=NA);
	set.seed(Vars["seed"]); pvars=sort(pnorm(rnorm(reps))); VAR=qlnormT(pvars,1,parms["sigma"]); strat=stratRes=codify(datTofit[,1:2])[,1]; if(extras["model"]!=4) stratRes=NA;
	Preds=(rtesVar4(parms,covars=datTofit[,CovarNms],Vars2=rep(VAR,each=Vars["simICs"]),ICs=Vars["simICs"],strat=strat,sitefeed=extras["model"]==4,vartype=extras["vartype"],batchsize=1e2,Time=Vars["simt"]))

	LLfull=t(apply(cbind(datTofit[,"Obs"],diag(datTofit[,"AdArea"])%*%pmax(Preds,0.01)), 1, function(x) dpois(x[1],x[-1])))
	LLs=-log(rowMeans(LLfull)); LL=sum(LLs);
	pASS=mean((Preds[,-(1:reps)]-Preds[,1:reps])>0.05); OKD=datTofit[,"Obs"]/datTofit[,"AdArea"]; pth=c(1,0.03)[1+(Starts["Region"]>2)];
	stats=round(c(LL=LL, R2state=bests(Preds>pth,OKD>pth), pASS=pASS, R2N=cor(OKD,rowMeans(Preds))^2, R2Nstate=bests(Preds,OKD)),2)
	if(extras["plotgive"]==4){ assign("Rsq",stats,Environment); return(cbind(datTofit,Preds)); }
	parms["fg"]=fg0; assign("REPORT",rbind(REPORT,c(parms,stats)),Environment); return(LL);
	#print(c(round(parms[parmsFit],3),stats[-5],REPORT[which.min(REPORT[,14]),14:17])); return(LL);
}



setwd("/home/vak32/barren/Outs")
#Some parameter names here deviate from manuscript: gamma=Rr, delta_R=RdeltaA, xi_P=deltaU, xi_A=fg, g_L=Lgrow, v=Ngrow.
parmnames=c("r","Rr","d",   "mu","b",   "deltaA","RdeltaA","deltaU","fg",    "fw","Lgrow","Ngrow",  "sigma"); parms=c(0,0,0,0,1,0,0,0,0,1,0,0,0.35); names(parms)=parmnames;
REPORT=matrix(nrow=0,ncol=18) #Track model statistics
extras=c(LLvrs=0,plotgive=0,allkelp=1,model=2,waves=6,vartype=2,Sims=1);


#Full set of models to run. Each NZ fit takes a few days, CA fits take a few weeks
Starts=rbind(as.matrix(expand.grid(Vers=0,allkelp=0,vartype=2,maxfg=0.35,simtype=2,Region=c(2,12),Rep=1,model=c(3:4,6:13))), #Base analysis
    as.matrix(expand.grid(Vers=0,allkelp=0,vartype=2,maxfg=0.35,simtype=2,Region=c(2,12),Rep=11:31,model=14:15)), #Grazing saturation analysis
    as.matrix(expand.grid(Vers=0,allkelp=0,vartype=2,maxfg=c(0.25,0.45),simtype=2,Region=c(2,12),Rep=1,model=c(3:4,6:13))), #Recruitment variation sensitivity analysis
	as.matrix(expand.grid(Vers=0,allkelp=0,vartype=2,maxfg=0.35,simtype=2,Region=12,Rep=1.1,model=c(3:4,6:13))), #Temporal autocorrelation sensitivity analysis (CA only)
	as.matrix(expand.grid(Vers=0,allkelp=0,vartype=2,maxfg=0.35,simtype=2,Region=c(2,12),Rep=11:35,model=c(3:4,6:13)))) #Site effects sensitivity analysis
Starts=Starts0[args,]; #select model
print(Starts)

if(Starts["Region"]==12) Region="CA" else Region="NZ"; ModelNumber=Starts["model"];
Environment=environment()




modelSet=list(m1p=c(1:3),m2p=c(1:3,11),m3p=c(1:3,4,10),m4p=c(1:3,6,7),m5p=c(1:3,6,7,8),m6p=c(1:3,6,7,9),m7p=c(1:3,6,7,9),m8p=c(1:3,11,4,10,6,7),
	m9p=c(1:3,11,4,10,6,7,9),m10p=c(1:3,11,4,10,6,7,9),m11p=c(1:3,11,4,10,6,7,8),m12p=c(1:3,11,4,10,6,7,8,9),m13=c(1:3,11,4,10,6,7,8,9),
	m14=c(1:3,6,7,9),m15=c(1:3,11,4,10,6,7,9),m16=c(1:3,11,4,10,6,7,8,9)) #14-16 are local Type II urchin functional responses
#Model components are environment (E), grazing (U), local behavior feedback (BL), reef-scale behavior feedback (BR), predator avoidance (BP), and local saturation (SL)
#Model numbers correspond to: 3=E, 4=U, 6=U+BL, 7=U+BR, 8=E+U, 9=E+U+BL, 10=E+U+BR, 11=E+U+BP, 12=E+U+BL+BP, 13=E+U+BR+BP, 14=U+SL, 15=E+U+SL, 16=E+U+SL+BP

# if(11%in%pfit) pfit=c(pfit,12)  #Parameter #12 Ngrow is Nitrate limitation; only fitted in preliminary CA analyses

#We explore saturation and behavior as alternative explanations for alternative stable states. Saturation is equivalent to a kelp density exponent of one in the behavioral function:
pfit=unlist(modelSet[[ModelNumber]]); if(!ModelNumber%in%(14:16)) GrazingInhibExponent=2 else GrazingInhibExponent=1;
#Denote whether feedbacks are local (model numbers 2,5) or reef-wide (model numbers 4,6)
extras["model"]=c(2,4)[1+(ModelNumber%in%c(7,10,13))] 
RV=Starts["maxfg"]; 


#To run optimization, we loosely constrain parameters to biologically plausible values,
#with differences in constraints on some parameters among regions reflecting regional differences in mean kelp and urchin densities.
#Note: constraints for fitted proportions can be <0 or >1 to ensure opimizer is not averse to edge cases; LoglikFun constrains all evaluated proportions to be >=0 and <=1
if(Region=="NZ"){ #New Zealand Fits:
	fgnz=3+9*(Starts["model"]>13) #Essentially no limit on saturation models (>13)
	upperBex=c(18,0.3,1.00,   7.000,1.00,   9.00,  12, 0.55, fgnz, 0.8,0.99,6.000,   RV+0.005)
	parmStr= c(8, 0.0,0.50,   1.000,0.995,  2.4,   10, 0.15, 2,    0.5,0.10,5.500,   RV)
	lowerBex=c(6,-2.0,0.15,   0.001,0.99,   0.25,0.2,-1.00,-1,   -3.0,0.00,0.001,   RV-0.005)
	LNthresh=0.25; UNthresh=5.5; #Initial model conditions (kelp density)
	DATA=readRDS("NENZdatFull.rds"); DATA=cbind(DATA[,1:7],AdArea=1,N=0,DATA[,-(1:7)]); #DATA=cbind(DATA[,1:8],N=0,DATA[,-(1:8)]) - for data on PC
	DATA=DATA[DATA[,"z"]>pmax(3-1.67*DATA[,"U"],1),] #Exclude shallow samples where other algae displace kelp and urchins
}


if(Region=="CA"){ #California Fits. For m13, bounds differ for r,d,mu
	upperBex=c(18,2,10,  7.000,3.000,   9,    12.0, 1.50,30,    0.8,0.99,6.000,   RV+0.005)
	parmStr= c(8,0,8,    1.000,0.995,   4.15, 10.0, 0.15, 4,    0.5,0.10,5.500,   RV)
	lowerBex=c(6,-2,2,   0.001,0.990,   0.25, 0.2,-1.00,-1,   -3.0,0.00,0.001,   RV-0.005)
	if(10%in%pfit) pfit=pfit[-which(pfit==10)]; #CA fits with wave attenuation by depth (f_w>0) always outperformed by fits with f_w constrained to 0, but the global optimizer sometimes misses solutions with f_w=0.
	LNthresh=0.025; UNthresh=0.95; #Initial model conditions (kelp density)
	DATA=readRDS("CICAalldat12all.rds")
	#Temporal autocorrelatin sensitivity: omitting consecutive years, removes 40-50% of data
	if(Starts["Rep"]==1.1){
		syu=unique(DATA[,1:2]); cy0=c(0,(head(syu[,2],-1)+1)==syu[-1,2]);
		dsy=cbind(syu,ss=c(0,syu[-1,1]==head(syu[,1],-1)),cy=(cy0*cumsum(cy0)%%2));
		DATA=DATA[codify(DATA[,1:2])%in%codify(dsy[(dsy[,3]*dsy[,4])==0,1:2]),]
	}
}

#Subset out 30% of sites for site-effects analysis
if(Starts["Rep"]>10){ set.seed(Starts["Rep"]); Sites=unique(DATA[,1]); DATA=DATA[DATA[,1]%in%sample(Sites,round(0.70*length(Sites))),]; }




#Implement optimization: global algorithm followed by local optimization starting from solution of global optimizer
#For each model, this step required several days for New Zealand and several weeks for California.
Ntries=1600; Nreps=15; 
g=nloptr(parmsFit=parmnames[pfit],extras=extras, opts=list("algorithm"="NLOPT_GN_DIRECT","xtol_rel"=1e-5,"maxeval"=round(0.8*Ntries),"print_level"=2), x0=parmStr[pfit],eval_f=LoglikFun,lb=lowerBex[pfit],ub=upperBex[pfit],reps=Nreps); 
nloptr(parmsFit=parmnames[pfit],extras=extras, opts=list("algorithm"="NLOPT_LN_COBYLA","xtol_rel"=1e-5,"maxeval"=round(0.2*Ntries),"print_level"=2), x0=g$solution,eval_f=LoglikFun,lb=lowerBex[pfit],ub=upperBex[pfit],reps=Nreps); 








