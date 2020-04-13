
# Input:
#   binary areal-average rainfall data
#      (avgimd365.dat)

# Output:
#   Figs. 1,2b,5
#      (Fig1a.pdf,Fig1b.pdf,Fig2b.pdf,Fig5a.pdf,Fig5b.pdf)
#   binary climatological seasonal cycles data
#      (clim120.dat)

rm(list=ls())

lplotfile = FALSE  # Plot figures (TRUE) or save figures (FALSE)

nt = 120     # summer (01Jun-28Sep,JJAS)
tref = 151   # index for 30Apr (fix with leap years beginning 02Jan in input data)
ntall = 365

# Input data:
# Each year (nyrs = 115 in total) is of length ntall = 365 days
# Leap years begin 02Jan so that tref = 151 is index for 30Apr
# and tref+c(1:nt) are indices for 01Jun-28Sep (nt = 120 day season, JJAS)

nyrs = 115
ntall = 365
tref = 151
nt = 120

# Areal-averaged IMD rainfall
prec = readBin('avgimd365.dat','numeric',n=ntall*nyrs*7,size='4',endian='little')
dim(prec) = c(ntall,nyrs,7)

# prec[,,1]   Central India (100%)
# prec[,,2]   All-India
# prec[,,3]   All-India w/o Central India
# prec[,,4:7] Small scale within Central India (10%) (4 regions)


# Function that computes harmonics
fftcomp <- function(var1) {
   nt = length(var1); nm  = floor(nt/2)
   even = FALSE; nm1 = nt/2; if (nm==nm1) even = TRUE
   comp <- array(0,dim=c(nt,1+nm))
   comp[,1] = mean(var1); var1 = var1 - mean(var1)
# Compute FFT
   rfcoef <- array(0,dim=nt)
   rfcoef[1:(nm+1)] = fft(var1)[1:(nm+1)]
   if (even) rfcoef[nm+1] = rfcoef[nm+1]/2.0
# Harmonic coefficients
   x <- c(0:(nt-1)); dim(x) = nt; lam  = 2*pi*x/nt
   amp  = (1.0/nt)*Mod(2*rfcoef[2:(nm+1)])
   phi  = Arg(2*rfcoef[2:(nm+1)])
# Reconstruct for each harmonic
   for (m in 1:nm) comp[,1+m] = amp[m]*cos(m*lam+phi[m])
   return(comp)
   }

# Compute harmonics for each year
nmall = floor(ntall/2)
comp <- array(0,dim=c(ntall,1+nmall,nyrs,7))
for (l in 1:nyrs) for (z in 1:7) comp[,,l,z] = fftcomp(prec[,l,z])


# Correlations for Fig. 1a,1b,2b

# ISMR seasonal means
ismr = apply(prec[tref+c(1:nt),,],c(2,3),mean)
# Sub-seasonal variances
# Only for Central India (100%), All-India, All-India w/o Central India [1:3]
yvars <- array(0,dim=c(nmall-3,nyrs,3))
for (z in 1:3) for (l in 1:nyrs) {
   for (i in 1:(nmall-4)) {
      temp = apply(comp[tref+c(1:nt),(4+i):(1+nmall),l,z],1,sum)
      yvars[i,l,z] = mean((temp-mean(temp))^2)
      }
   temp = comp[tref+c(1:nt),1+nmall,l,z]
   yvars[nmall-3,l,z] = mean((temp-mean(temp))^2)
   }

# Two time periods, 1981-2010 and 1901-2015
iyr <- list(); iyr[[1]] = 81:110; iyr[[2]] = 1:115
ycor <- array(0,dim=c(nmall-3,3,2))
for (zz in 1:2) for (z in 1:3) {
   for (i in 1:(nmall-3)) ycor[i,z,zz] = cor(ismr[iyr[[zz]],z],yvars[i,iyr[[zz]],z])
   ycor[,z,zz] = rev(ycor[,z,zz])
   }
xperall = rev(ntall/c(4:nmall))


if (!lplotfile) pdf('Fig1a.pdf')
plot(xperall,ycor[,1,1],type='l',col='gray',lwd=5,xlab='period (days)',ylab='',main='',cex.lab=1.5,cex.axis=1.5,cex.main=1.8,xlim=c(2,100),ylim=c(-0.25,0.75),xaxt='n',yaxt='n',log='x')
lines(xperall,ycor[,1,2],lwd=5,col='black')
text(35,0.73,'1981-2010',cex=1.5,col='gray',pos=4)
text(35,0.67,'1901-2015',cex=1.5,col='black',pos=4)
text(3.5,0.73,'a) Central',cex=2)
lines(c(1,120),rep(0.183,2),lty=2,lwd=2,col='black')
lines(c(1,120),rep(0.361,2),lty=2,lwd=2,col='gray')
lines(c(1,120),rep(-0.183,2),lty=2,lwd=2,col='black')
axis(1,at=c(2,3,4,6,10,15,30,60,90),col.axis='black',las=1,cex.axis=1.5)
axis(2,at=c(-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7),col.axis='black',las=1,cex.axis=1.5)
if (!lplotfile) dev.off()

if (!lplotfile) pdf('Fig1b.pdf')
plot(xperall,ycor[,2,1],type='l',col='gray',lwd=5,xlab='period (days)',ylab='',main='',cex.lab=1.5,cex.axis=1.5,cex.main=1.8,xlim=c(2,100),ylim=c(-0.25,0.75),xaxt='n',yaxt='n',log='x')
lines(xperall,ycor[,2,2],lwd=5,col='black')
text(35,0.73,'1981-2010',cex=1.5,col='gray',pos=4)
text(35,0.67,'1901-2015',cex=1.5,col='black',pos=4)
text(3.8,0.73,'b) All-India',cex=2)
lines(c(1,120),rep(0.183,2),lty=2,lwd=2,col='black')
lines(c(1,120),rep(0.361,2),lty=2,lwd=2,col='gray')
lines(c(1,120),rep(-0.183,2),lty=2,lwd=2,col='black')
axis(1,at=c(2,3,4,6,10,15,30,60,90),col.axis='black',las=1,cex.axis=1.5)
axis(2,at=c(-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7),col.axis='black',las=1,cex.axis=1.5)
if (!lplotfile) dev.off()

cols = c('black','#FF9933','#138808','#000080')
pdf('Fig2b.pdf')
plot(xperall,ycor[,2,2],type='l',lwd=5,col='#FF9933',xlab='period (days)',ylab='',main='',cex.lab=1.5,cex.axis=1.5,cex.main=1.8,xlim=c(2,100),ylim=c(-0.25,0.75),xaxt='n',yaxt='n',log='x')
lines(xperall,ycor[,1,2],lwd=5,col='black')
lines(xperall,ycor[,3,2],lwd=5,col='#138808')
lines(c(1,120),rep( 0.183,2),lty=2,col='black')
lines(c(1,120),rep(-0.183,2),lty=2,col='black')
text(45,0.49,'All-India',cex=1.5,col='#138808')
text(45,0.44,'without Central',cex=1.5,col='#138808')
text(55,0.14,'Central',cex=1.5,col='black')
text(55,0.07,'All-India',cex=1.5,col='#FF9933')
text(3.8,0.73,'b) correlation',cex=2)
axis(1,at=c(2,3,4,6,10,15,30,60,90),col.axis='black',las=1,cex.axis=1.5)
axis(2,at=c(-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8),col.axis='black',las=1,cex.axis=1.5)
dev.off()


# ISMR climatological seasonal cycles (subtract constant residual also)
clims <- array(0,dim=c(nt,nyrs,7))
clim <- array(0,dim=c(nt,7))
for (z in 1:7) {
   clims[,,z] = apply(comp[tref+c(1:nt),1:4,,z],c(1,3),sum)
   clim[,z] = apply(clims[,,z],1,mean)
   clim1 = mean(prec[tref+c(1:nt),,z])
   clim2 = mean(clim[,z])
   clim[,z] = clim[,z] + clim1-clim2
   }

# Write out climatological seasonal cycles
writeBin(as.vector(clim),'clim120.dat',size='4',endian='little')


# Plot raw data for Fig. 5
# Random wet and dry years
set.seed(8)
idry = which(ismr[,1]<sort(ismr[,1])[ 12]); idry = idry[ceiling(runif(1)*11)]
iwet = which(ismr[,1]>sort(ismr[,1])[104]); iwet = iwet[ceiling(runif(1)*11)]

# 30+ day band
iband = c(1:13)

ylim = c(0,39)
pdf('Fig5a.pdf')
plot(12.5+prec[tref+c(1:nt),idry,1],type='l',lwd=3,col='black',ylim=ylim,cex.lab=1.5,cex.axis=1.5,xaxt='n',yaxt='n',xlab='',ylab='')
for (z in 1:4) lines(prec[tref+c(1:nt),idry,3+z]/5,col='gray')
lines(12.5+apply(comp[tref+c(1:nt),iband,idry,1],1,sum),lwd=3,col='#FF9933')
lines(12.5+clim[,1],lwd=3,col='#138808')
lines(rep(0,nt),lty=2); lines(rep(12.5,nt),lty=2)
text(-3,ylim[2]-(ylim[2]-ylim[1])/40,paste('a) dry year (',1900+idry,')',sep=''),cex=1.5,pos=4)
text(-3,8,'10%',cex=1.5,col='gray',pos=4)
text(-3,10.5,'100%',cex=1.5,col='black',pos=4)
text(75,38,'30+ day periods',col='#FF9933',cex=1.3)
text(75,36,'climatological',col='#138808',cex=1.3)
axis(1,at=c(1,31,62,93,120),labels=c('01Jun','01Jul','01Aug','01Sep','28Sep'),col.axis='black',las=1,cex.axis=1.5)
axis(2,at=c(0,4,8,12.5,17.5,22.5,27.5,32.5,37.5),labels=c('0','20','40','0','5','10','15','20','25'),col.axis='black',las=1,cex.axis=1.5)
dev.off()

pdf('Fig5b.pdf')
plot(12.5+prec[tref+c(1:nt),iwet,1],type='l',lwd=3,col='black',ylim=ylim,cex.lab=1.5,cex.axis=1.5,xaxt='n',yaxt='n',xlab='',ylab='')
for (z in 1:4) lines(prec[tref+c(1:nt),iwet,3+z]/5,col='gray')
lines(12.5+apply(comp[tref+c(1:nt),iband,iwet,1],1,sum),lwd=3,col='#FF9933')
lines(12.5+clim[,1],lwd=3,col='#138808')
lines(rep(0,nt),lty=2); lines(rep(12.5,nt),lty=2)
text(-3,ylim[2]-(ylim[2]-ylim[1])/40,paste('b) wet year (',1900+iwet,')',sep=''),cex=1.5,pos=4)
text(-3,8,'10%',cex=1.5,col='gray',pos=4)
text(-3,10.5,'100%',cex=1.5,col='black',pos=4)
text(100,38,'30+ day periods',col='#FF9933',cex=1.3)
text(100,36,'climatological',col='#138808',cex=1.3)
axis(1,at=c(1,31,62,93,120),labels=c('01Jun','01Jul','01Aug','01Sep','28Sep'),col.axis='black',las=1,cex.axis=1.5)
axis(2,at=c(0,4,8,12.5,17.5,22.5,27.5,32.5,37.5),labels=c('0','20','40','0','5','10','15','20','25'),col.axis='black',las=1,cex.axis=1.5)
dev.off()


