library(rjags)
library(statmod)
library(coda)
library(MCMCpack)

set.seed(1345221)
ada=20000
iter=200000
nchains=3
thn=2

# define species
spclass <- c("sp","rn","bd","rd")
names(spclass) <- c("Spotted","Ribbon","Bearded","Ringed")
scertclass <- c(3,2,1)
names(scertclass) <- c("guess","likely","pos")

# define age classes
ageclass <- c("p","a")
names(ageclass) <- c("pup","nonpup") 
acertclass <- c(3,2,1)
names(acertclass) <- c("guess","likely","pos")

nspec <- length(spclass)
nage <- length(ageclass)
nscert <- length(scertclass)
nacert <- length(acertclass)
spcats <- nscert*nacert*nage

round1 = data.frame(read.csv("data/IceSeals_CompiledSpeciesID.csv", header=T))
round1key = data.frame(read.csv("data/IceSeals_CompiledSpeciesID_CodeKey.csv", header=TRUE))
round1$detection_id <- subset(round1key, type=="detection_id")$description[match(round1$hotspot_number, subset(round1key, type=="detection_id")$id)]
round1$image_number <- subset(round1key, type=="image_number")$description[match(round1$image_number, subset(round1key, type=="image_number")$id)]

temp <- round1

# remove unknowns and anomolies from round 1 data
temp[which(temp$sp_id_conf %in% c('unk','x')),c("sp_id_conf","age_class_conf")]<-NA

temp$detection_id <- factor(temp$detection_id)
temp <- temp[order(as.numeric(levels(temp$detection_id))[temp$detection_id],as.numeric(temp$obs_id )),]

R=length(unique(temp$detection_id)) #number of seals
num_obs = length(unique(temp$obs_id)) # max observers

# --- Species Matrix ---
temp$SPECIES <- substr(temp$sp_id_conf,1,2)
tempy_it <- matrix(NA, nrow = R, ncol = num_obs)
rownames(tempy_it) <- levels(temp$detection_id)
colnames(tempy_it) <- 1:num_obs
tempy_it[cbind(as.numeric(temp$detection_id), as.numeric(temp$obs_id))] <- temp$SPECIES

# --- Species Confidence Matrix ---
temp$SPECIES_CONF <- substr(temp$sp_id_conf,3,3)
for(c in 1:nscert){
  temp$SPECIES_CONF[which(temp$SPECIES_CONF==scertclass[c])] <- names(scertclass)[c]
}
tempcert_it <- matrix(NA, nrow = R, ncol = num_obs)
tempcert_it[cbind(as.numeric(temp$detection_id), as.numeric(temp$obs_id))] <- temp$SPECIES_CONF

# --- Age Matrix ---
temp$GROSS_AGE <- substr(temp$age_class_conf,1,1)
for(a in 1:nage){
  temp$GROSS_AGE[which(temp$GROSS_AGE==ageclass[a])]<- names(ageclass)[a]
}
temp_age <- matrix(NA, nrow = R, ncol = num_obs)
temp_age[cbind(as.numeric(temp$detection_id), as.numeric(temp$obs_id))] <- temp$GROSS_AGE

# --- Age Confidence Matrix ---
temp$GROSS_AGE_CONF <- substr(temp$age_class_conf,2,2)
for(c in 1:nacert){
  temp$GROSS_AGE_CONF[which(temp$GROSS_AGE_CONF==acertclass[c])] <- names(acertclass)[c]
}
tempcert_age <- matrix(NA, nrow = R, ncol = num_obs)
tempcert_age[cbind(as.numeric(temp$detection_id), as.numeric(temp$obs_id))] <- temp$GROSS_AGE_CONF

# --- Construct Final Observation Matrix (y_cert_it) ---
y_cert_it <- matrix(NA, nrow = R, ncol = num_obs)
rownames(y_cert_it) <- levels(temp$detection_id)
colnames(y_cert_it) <- paste("obs", seq(1:num_obs), sep="")

# Dynamic state generator (accommodates 4 species seamlessly)
sp = rep(spclass, each=spcats)
sp_cert = rep(rep(names(scertclass), each=nacert*nage), times=nspec)
age = rep(rep(names(ageclass), each=nacert), nscert*nspec)
age_cert = rep(names(acertclass), nage*nscert*nspec)

for(i in 1:(nscert*nage*nacert*nspec)){
  matches <- which(tempy_it==sp[i] & tempcert_it==sp_cert[i] & temp_age==age[i] & tempcert_age==age_cert[i])
  y_cert_it[matches] <- i
}

# --- Testing / Dropping Conflicting Positive IDs ---
testsp=matrix(0,ncol=nspec,nrow=nrow(y_cert_it))
for(i in 1:R){
  for(j in 1:nspec){
    for(k in 1:(nage*nacert)){
      testsp[i,j]=testsp[i,j]+any(y_cert_it[i,]==(spcats*(j-1)+(nscert-1)*nage*nacert+k),na.rm=T)
    }
  }
  if(sum(testsp[i,]>0)>1){
    print(c("bad IDs ",rownames(y_cert_it)[i],y_cert_it[i,]))
  }  
}

if(length(which(apply(testsp>0,1,sum)>1))>0){
  y_cert_it=y_cert_it[-which(apply(testsp>0,1,sum)>1),]
}

R=nrow(y_cert_it)

testage=matrix(0,ncol=nage,nrow=nrow(y_cert_it))
for(i in 1:R){
  for(j in 1:nage){
    for(k in 1:nspec){
      for(l in seq(nacert,(nage*nscert)*nacert,(nage*nacert))*(k<(nspec+1))+nspec*(k==(nspec+1))+nacert*(j-1)+spcats*(k-1)){
        testage[i,j]=testage[i,j]+any(y_cert_it[i,]==l,na.rm=T)
      }
    }
  }
  if(sum(testage[i,]>0)>1){
    print(c("bad IDs ",rownames(y_cert_it)[i],y_cert_it[i,]))
  }  
}

if(length(which(apply(testage>0,1,sum)>1))>0){
  y_cert_it=y_cert_it[-which(apply(testage>0,1,sum)>1),]
}

y_cert_it[which(y_cert_it==0)]=NA #remove bad ones for now

R=nrow(y_cert_it)

# --- JAGS Initial Values ---
z1=numeric(nrow(y_cert_it))

init.values<-function(){
  for(i in 1:R) {
    if(any(sp_cert[y_cert_it[i,]]=="pos",na.rm=T))  {
      spec=sp[y_cert_it[i,which(sp_cert[y_cert_it[i,]]=="pos")][1]]
      ztemp=which(spec==spclass)
      z1[i]=(ztemp-1)*nage
      if(any(age_cert[y_cert_it[i,]]=="pos",na.rm=T)) {
        ag=age[y_cert_it[i,which(age_cert[y_cert_it[i,]]=="pos")][1]]
        z1[i]=z1[i]+which(ag==names(ageclass))
      } else {
        z1[i]=z1[i]+sample(c(1,nage),1)
      }
    } else if(any(age_cert[y_cert_it[i,]]=="pos",na.rm=T)) {
      ag=age[y_cert_it[i,which(age_cert[y_cert_it[i,]]=="pos")][1]]    
      z1[i]=which(ag==names(ageclass))+sample(c(0,nage),1)
    } else {
      z1[i]=sample(seq(1,nage*nspec),1)
    }
  }
  return(list(psi=c(rdirichlet(1,rep(1,nage*nspec))),z=z1))
}

obs_indices <- which(!is.na(y_cert_it), arr.ind = TRUE)

y_long <- y_cert_it[obs_indices]

seal_idx <- obs_indices[, "row"]
obs_idx  <- obs_indices[, "col"] 

n_obs_total <- length(y_long)

# --- Run Model ---
data <- list(
  y_long = y_long,           
  seal_idx = seal_idx, 
  obs_idx = obs_idx,
  n = num_obs,
  n_obs_total = n_obs_total,
  R = R, 
  nspec = nspec, 
  nage = nage, 
  nscert = nscert, 
  nacert = nacert, 
  dirchprior = rep(1, nage * nspec)
)

inits <- init.values()
params <- c('psi','alpha','beta','delta','gamma')
mod <- jags.model('code/speciesID_obs.jag',data,inits,n.chains=nchains,n.adapt=ada)
gc()
single.sim <- coda.samples(mod,params,n.iter=iter)

save.image("speciesID_obs.RData")

# --- Process Results ---
psi=single.sim[,(paste("psi[",seq(1:(nspec*nage)),"]",sep=""))]
psi=as.matrix(psi)[seq(1,iter*nchains,thn),]

alpha=single.sim[,(paste0("alpha[",rep(seq(1:num_obs),each=nspec*nspec),",",rep(seq(1:nspec),each=nspec),",",seq(1:nspec),"]"))]
alpha=as.matrix(alpha)[seq(1,iter*nchains,thn),]

beta=single.sim[,(paste0("beta[",rep(seq(1:num_obs),each=nspec*nspec*nscert),",",rep(seq(1:nspec),each=nspec*nscert),",",rep(seq(1:nspec),each=nscert),",",seq(1:nscert),"]"))]
beta=as.matrix(beta)[seq(1,iter*nchains,thn),]

delta=single.sim[,(paste0("delta[",rep(seq(1:num_obs),each=nage*nage),",",rep(seq(1:nage),each=nage),",",seq(1:nage),"]"))]
delta=as.matrix(delta)[seq(1,iter*nchains,thn),]

gamma=single.sim[,(paste0("gamma[",rep(seq(1:num_obs),each=nage*nage*nacert),",",rep(seq(1:nage),each=nage*nacert),",",rep(seq(1:nage),each=nacert),",",seq(1:nacert),"]"))]
gamma=as.matrix(gamma)[seq(1,iter*nchains,thn),]

rm(single.sim)
gc()

p=array(0,dim=c(nspec*nage,spcats*nspec,num_obs,nrow(psi)))
gc()

for(i in 1:num_obs){
  for(l in 1:nspec)  { # true species 
    for(j in 1:nage){ # true age class    
      for(s in 1:nspec){ # observed species
        for(c_s in 1:nscert){ # species confidence
          for(a in 1:nage){  # observed age
            for(c_a in 1:nacert){  # age confidence
              p[(l-1)*nage+j,spcats*(s-1)+nage*nacert*(c_s-1)+nacert*(a-1)+c_a,i,] <- alpha[,paste0("alpha[",i,",",l,",",s,"]")]*beta[,paste0("beta[",i,",",l,",",s,",",c_s,"]")]*delta[,paste0("delta[",i,",",j,",",a,"]")]*gamma[,paste0("gamma[",i,",",j,",",a,",",c_a,"]")]
              gc()
            }
          }
        }
      }
    }    
  }
}

meanp=apply(p,c(1,2,3),mean)
rownames(meanp) <- paste(rep(names(spclass),each=nage),names(ageclass))
colnames(meanp) <- paste(rep(paste(rep(toupper(spclass),each=nscert),names(scertclass)),each=nage*nacert),rep(names(ageclass),times=rep(nacert,nage)),rep(names(acertclass),nspec))

tiff(file="meanp_obs.tiff",width=12,height=6,units="in",res=600,compression="lzw")

par(mfrow=c(2,4),mar=c(2, 1, 0.5, 1)) # Updated to 2x4 grid to fit up to 8 observers
par(oma = c (3, 3, 0, 5))

# Corrected for 4 species: 12 base palettes repeated appropriately
cols <- c(
  rep(
    c(colorRampPalette(c("white", "red"))(100)[seq(30, 100, length = nscert)],
      colorRampPalette(c("white", "green"))(100)[seq(30, 100, length = nscert)],
      colorRampPalette(c("white", "yellow"))(100)[seq(30, 100, length = nscert)],
      colorRampPalette(c("white", "blue"))(100)[seq(30, 100, length = nscert)]), 
    times = rep(nage * nacert, nspec * nscert)
  )
)

dens=seq(15,35,length=nacert)
angs=rep(unique(c(45,360-45,seq(45,360-45,length=nage))),each=nacert)

for(i in 1:num_obs){
  barplot(t(as.matrix(meanp[,,i])),names.arg=paste0(rep(toupper(spclass),each=nage),toupper(substr(names(ageclass),1,1))),xlab=NA,ylab=NA,cex.lab=0.6,cex.axis=0.6,cex.names=0.6,legend.text=F,col =cols)
  barplot(t(as.matrix(meanp[,,i])),add=T,col=1,density=dens,angle=angs,names.arg=paste0(rep(toupper(spclass),each=nage),toupper(substr(names(ageclass),1,1))),xlab=NA,ylab=NA,cex.lab=0.6,cex.axis=0.6,cex.names=0.6,legend.text=F)
}

mtext("True species/age", side = 1, outer = T, line=2)
mtext("Probability of observed species/age and confidence", side = 2, outer = T, line=2)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("right",c(paste(tolower(names(spclass)),"seal"),names(ageclass)), xpd = TRUE, horiz = F, inset = c(0,0), bty = "n", angle=c(rep(45,nspec),unique(c(45,360-45,seq(45,360-45,length=nage)))),density=c(rep(-1,nspec),rep(40,nage)),x = "right", fill=c(cols[seq(nage*nacert,spcats*nspec,spcats)],rep(1,nage)),cex=.75)
dev.off()

p12=array(0,dim=c(nspec,nspec*nscert,num_obs,length(seq(1,iter*nchains,thn))))

for(i in 1:num_obs){
  for(l in 1:nspec){
    for(s in 1:nspec){
      for(c_s in 1:nscert){ 
        p12[l,(s-1)*nscert+c_s,i,]=alpha[,paste0("alpha[",i,",",l,",",s,"]")]*beta[,paste0("beta[",i,",",l,",",s,",",c_s,"]")]
      }
    }
  }
}

dimnames(p12)=list(toupper(spclass),paste0(rep(toupper(spclass),each=nscert),names(scertclass)),paste0("obs",1:num_obs))
save(p12,file="p12_obs.RData")

pspecage=array(0,dim=c(nspec*nage,nspec*nage,num_obs,length(seq(1,iter*nchains,thn))))
for(i in 1:num_obs){
  for(l in 1:nspec){
    tempa=alpha[,rep(paste0("alpha[",i,",",1:nspec,",",l,"]"),each=nage)]
    for(k in 1:nage){
      tempb=delta[,rep(paste0("delta[",i,",",1:nage,",",k,"]"),nage)]
      for(s in 1:nspec){
        for(j in 1:nage){
          pspecage[(s-1)*nage+j,(l-1)*nage+k,i,] = tempa[,(s-1)*nage+j]*tempb[,j]
        }
      }
    }
  }
}

p12mean=apply(p12,c(1,2,3),mean)
tquantiles=function(x){quantile(x,c(.025))}
p12LCI=apply(p12,c(1,2,3),tquantiles)
tquantiles=function(x){quantile(x,c(.975))}
p12UCI=apply(p12,c(1,2,3),tquantiles)

p12summary=array(0,dim=c(nspec*nscert*nspec,3,num_obs))
for(i in 1:num_obs){
  p12summary[,1,i]=c(t(p12mean[,,i]))
  p12summary[,2,i]=c(t(p12LCI[,,i]))
  p12summary[,3,i]=c(t(p12UCI[,,i]))
}
rownames(p12summary) <- paste0(rep(paste0("true",toupper(spclass),":"),each=nspec*nscert),paste0(rep(toupper(spclass),each=nscert),names(scertclass)))
colnames(p12summary) <- c("p","lci","uci")

pspecagemean=apply(pspecage,c(1,2,3),mean)
tquantiles=function(x){quantile(x,c(.025))}
pspecageLCI=apply(pspecage,c(1,2,3),tquantiles)
tquantiles=function(x){quantile(x,c(.975))}
pspecageUCI=apply(pspecage,c(1,2,3),tquantiles)

pspecagesummary=array(0,dim=c(nspec*nage*nspec*nage,3,num_obs))
for(i in 1:num_obs){
  pspecagesummary[,1,i]=c(t(pspecagemean[,,i]))
  pspecagesummary[,2,i]=c(t(pspecageLCI[,,i]))
  pspecagesummary[,3,i]=c(t(pspecageUCI[,,i]))
}

rownames(pspecagemean) <- paste(rep(names(spclass),each=nage),names(ageclass))
colnames(pspecagemean) <- paste(rep(names(spclass),each=nage),names(ageclass))

tiff(file="pspecagemean_obs.tiff",width=12,height=6,units="in",res=600,compression="lzw")
par(mfrow=c(2,4),mar=c(2, 1, 0.5, 1))
par(oma = c (3, 3, 0, 5))

# Corrected for 4 species: 8 colors
cols <- c(colorRampPalette(c("white", "red"))(100)[rep(50,nage)],
          colorRampPalette(c("white", "green"))(100)[rep(50,nage)],
          colorRampPalette(c("white", "yellow"))(100)[rep(50,nage)],
          colorRampPalette(c("white", "blue"))(100)[rep(50,nage)])

for(i in 1:num_obs){
  barplot(t(as.matrix(pspecagemean[,,i])),names.arg=paste0(rep(toupper(spclass),each=nage),toupper(substr(names(ageclass),1,1))),xlab=NA,ylab=NA,cex.lab=.75,cex.axis=0.55,cex.names=0.55,legend.text=F,col = cols)
  barplot(t(as.matrix(pspecagemean[,,i])),names.arg=paste0(rep(toupper(spclass),each=nage),toupper(substr(names(ageclass),1,1))),xlab=NA,ylab=NA,cex.lab=.75,cex.axis=0.55,cex.names=0.55,legend.text=F,add=T,col=1,angle=unique(c(45,360-45,seq(45,360-45,length=nage))),density=rep(25,nspec))
}

mtext("True species/age", side = 1, outer = T, line=2)
mtext("Probability of observed species/age", side = 2, outer = T, line=2)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("right",c(paste(tolower(names(spclass)),"seal"),names(ageclass)), xpd = TRUE, horiz = F, inset = c(0,0), bty = "n", angle=c(rep(45,nspec),unique(c(45,360-45,seq(45,360-45,length=nage)))),density=c(rep(-1,nspec),rep(40,nage)),x = "right", fill=c(cols[seq(nage,nage*nspec,nage)],rep(1,nage)),cex=.75)
dev.off()

pspecagemean_overall=apply(pspecage,c(1,2),mean)
tquantiles=function(x){quantile(x,c(.025))}
pspecageLCI_overall=apply(pspecage,c(1,2),tquantiles)
tquantiles=function(x){quantile(x,c(.975))}
pspecageUCI_overall=apply(pspecage,c(1,2),tquantiles)

pspecagesummary_overall=matrix(0,nrow=nspec*nage*nspec*nage,ncol=3)
pspecagesummary_overall[,1]=c(t(pspecagemean_overall))
pspecagesummary_overall[,2]=c(t(pspecageLCI_overall))
pspecagesummary_overall[,3]=c(t(pspecageUCI_overall))

rownames(pspecagemean_overall) <- paste(rep(names(spclass),each=nage),names(ageclass))
colnames(pspecagemean_overall) <- paste(rep(names(spclass),each=nage),names(ageclass)) # Fixed hard-coded 2-species limitation

tiff(file="pspecagemean_obs.tiff",width=7.25,height=5,units="in",res=600,compression="lzw")
par(mfrow=c(1, 1), mar=c(5, 5, 4, 8))
barplot(t(as.matrix(pspecagemean_overall)),names.arg=paste0(rep(toupper(spclass),each=nage),toupper(substr(names(ageclass),1,1))),xlab=NA,ylab=NA,cex.lab=.75,cex.axis=0.55,cex.names=0.55,legend.text=F,col = cols)
barplot(t(as.matrix(pspecagemean_overall)),names.arg=paste0(rep(toupper(spclass),each=nage),toupper(substr(names(ageclass),1,1))),xlab=NA,ylab=NA,cex.lab=.75,cex.axis=0.55,cex.names=0.55,legend.text=F,add=T,col=1,angle=unique(c(45,360-45,seq(45,360-45,length=nage))),density=rep(25,nspec))
dev.off()

spec_psi=matrix(0,nrow=length(seq(1,iter*nchains,thn)),ncol=nspec)
for(s in 1:nspec){
  spec_psi[,s]=rowSums(psi[,(s-1)*nage+1:nage])
}

Pspec=array(0,dim=c(nspec,nspec,num_obs))
varPspec=array(0,dim=c(nspec,nspec,num_obs))
for(i in 1:num_obs){
  for(s in 1:nspec){
    for(k in 1:nspec){
      tmp=spec_psi[,s]*alpha[,paste0("alpha[",i,",",s,",",k,"]")]/apply(spec_psi*alpha[,paste0("alpha[",i,",",1:nspec,",",k,"]")],1,sum)
      Pspec[s,k,i] = mean(tmp)
      varPspec[s,k,i] = var(tmp)
    }
  }
}

Pspecfull=array(0,dim=c(nspec,nspec*nscert,num_obs))
varPspecfull=array(0,dim=c(nspec,nspec*nscert,num_obs))
for(i in 1:num_obs){
  for(s in 1:nspec){
    for(k in 1:(nspec*nscert)){
      tmp=spec_psi[,s]*p12[s,k,i,]/apply(spec_psi*t(p12[1:nspec,k,i,]),1,sum)
      Pspecfull[s,k,i] = mean(tmp)
      varPspecfull[s,k,i] = var(tmp)
    }
  }
}
rownames(Pspecfull) <- rownames(p12)
colnames(Pspecfull) <- colnames(p12)

age_psi=matrix(0,nrow=length(seq(1,iter*nchains,thn)),ncol=nage)
for(s in 1:nage){
  age_psi[,s]=rowSums(psi[,seq(s,nage*nspec,nage)])
}

Page=array(0,dim=c(nage,nage,num_obs))
varPage=array(0,dim=c(nage,nage,num_obs))
for(i in 1:num_obs){
  for(s in 1:nage){
    for(k in 1:nage){
      tmp=age_psi[,s]*delta[,paste("delta[",i,",",s,",",k,"]",sep="")]/apply(age_psi*delta[,paste("delta[",i,",",1:nage,",",k,"]",sep="")],1,sum)
      Page[s,k,i] = mean(tmp)
      varPage[s,k,i] = var(tmp)
    }
  }
}

P=array(0,dim=c(nage*nspec,nage*nspec,num_obs))
varP=array(0,dim=c(nage*nspec,nage*nspec,num_obs))
for(i in 1:num_obs){
  for(l in 1:nspec){
    tempa=alpha[,rep(paste0("alpha[",i,",",1:nspec,",",l,"]"),each=nage)]
    for(k in 1:nage){
      tempb=delta[,rep(paste0("delta[",i,",",1:nage,",",k,"]"),nspec)]
      for(s in 1:nspec){
        for(j in 1:nage){
          tmp=psi[,(s-1)*nage+j]*tempa[,(s-1)*nage+j]*tempb[,j]/apply(psi*tempa*tempb,1,sum)
          P[(s-1)*nage+j,(l-1)*nage+k,i] = mean(tmp)
          varP[(s-1)*nage+j,(l-1)*nage+k,i] = var(tmp)
        }
      }
    }
  }
}

fullP=array(0,dim=c(nspec*nage,spcats*nspec,num_obs))
fullvarP=array(0,dim=c(nspec*nage,spcats*nspec,num_obs))
for(i in 1:num_obs){
  for(s in 1:(nspec*nage)){
    for(k in 1:(spcats*nspec)){
      tmp=psi[,s]*p[s,k,i,]/apply(psi*t(p[,k,i,]),1,sum)
      fullP[s,k,i] = mean(tmp)
      fullvarP[s,k,i] = var(tmp)
    }
  }
}

rownames(Pspec) <- paste0("true.",tolower(names(spclass)))
colnames(Pspec) <- paste0("obs.",tolower(names(spclass)))

rownames(varPspec) <- paste0("true.",tolower(names(spclass)))
colnames(varPspec) <- paste0("obs.",tolower(names(spclass)))

rownames(P) <- paste0("true.",paste0(rep(tolower(names(spclass)),each=nage),names(ageclass)))
colnames(P) <- paste0("obs.",paste0(rep(tolower(names(spclass)),each=nage),names(ageclass)))

rownames(varP) <- paste0("true.",paste0(rep(tolower(names(spclass)),each=nage),names(ageclass)))
colnames(varP) <- paste0("obs.",paste0(rep(tolower(names(spclass)),each=nage),names(ageclass)))

rownames(fullP) <- paste0("true.",paste0(rep(tolower(names(spclass)),each=nage),names(ageclass)))
colnames(fullP) <- paste0(sp,sp_cert,age,age_cert)

rownames(fullvarP) <- paste0("true.",paste0(rep(tolower(names(spclass)),each=nage),names(ageclass)))
colnames(fullvarP) <- paste0(sp,sp_cert,age,age_cert)