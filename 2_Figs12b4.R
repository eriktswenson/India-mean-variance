
# Input:
#   binary areal-averaged rainfall data
#      (avgimd365.dat)

# Output:
#   Figs. 1, 2b and 4
#      (Fig1a.pdf,Fig1b.pdf,Fig2b.pdf
#       Fig4a.pdf,Fig4b.pdf,Fig4c.pdf,
#       Fig4d.pdf,Fig4e.pdf,Fig4f.pdf)


rm(list=ls())

lplotfile = FALSE  # Plot figures (TRUE) or save figures (FALSE)

ntall = 365   # days in year
nyrs  = 115   # years
nt = 120      # summer (01Jun-28Sep,JJAS)
tref = 151    # index for 30Apr
nreg = 3+4*9  # regions for spatial averaging


# Read in areal-averaged IMD rainfall
prec = readBin('avgimd365.dat','numeric',n=ntall*nyrs*nreg,size='4',endian='little')
dim(prec) = c(ntall,nyrs,nreg)

# prec[,,1]       # Central India
# prec[,,2]       # All-India
# prec[,,3]       # All-India w/o Central India
# prec[,,4:39]    # Sub-regions within Central India


# Function for harmonic decomposition
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

# Number of waves
nmall = floor(ntall/2)

nm = nmall - 3   # all harmonics except 1,2,3
#nm =  150 - 3    # 150 harmonics except 1,2,3 [Saha et al. 2019]

# Corresponding periods (reversed for later plotting)
x1 = rev(365/(3+c(1:nm)))

# Calculate harmonics
comp <- array(0,dim=c(ntall,1+nmall,nyrs,nreg))
for (l in 1:nyrs) for (z in 1:nreg) {
   comp[,,l,z] = fftcomp(prec[,l,z])
   }


# Seasonal means
ymean = apply(prec[tref+c(1:nt),,],c(2,3),mean)

# Consider all years or subset
iyrs1 = 1:nyrs
#iyrs1 = 81:110

# Correlation with Central India or Nino3 (*-1)
#nino3 = matrix(scan('nino3.txt'),nrow=13)
#nino3 = apply(nino3[1+c(6:9),],2,mean) * -1
#ymean[,1] = nino3


yvars <- array(0,dim=c(nm,nyrs,nreg))
ycor  <- array(0,dim=c(nm,nreg))
ycor1 = ycor

for (k in 1:nreg) {
# Calculate cumulative variances
   for (j in 1:nyrs) for (i in 1:nm) {
      temp = apply(comp[tref+c(1:nt),4+c(i:nm),j,k,drop=FALSE],1,sum)
      yvars[i,j,k] = mean((temp-mean(temp))^2)
      }
# Calculate local correlation and fixed Central India correlation
   for (i in 1:nm) {
       ycor[i,k] = cor(ymean[iyrs1,k],yvars[i,iyrs1,k])
      ycor1[i,k] = cor(ymean[iyrs1,1],yvars[i,iyrs1,k])
      }
# Correlations for Fig. 1a, 1b and 2b only with 
# Central India, All-India and All-India w/o Central India
   if (k==3) {
# Correlations over every possible 30-year period
      ycors30 <- array(0,dim=c(nm,3,86))
      for (kkk in 1:86) {
         iyrs = c(1:30)-1+kkk
         for (kk in 1:3) for (i in 1:nm) ycors30[i,kk,kkk] = cor(ymean[iyrs,kk],yvars[i,iyrs,kk])
         }
      ycors30sort = apply(ycors30,c(1,2),sort)
      ycor30med = (ycors30sort[43,,]+ycors30sort[44,,])/2
      ycor30hi  = ycors30sort[65,,]
      ycor30lo  = ycors30sort[22,,]
# Correlations over 1981-2010 period
      iyrs = 81:110
      ycor30 = ycors30[,,81]
      if (!lplotfile) pdf('Fig1a.pdf')
      plot(x1,rev(ycor30med[,1]),type='l',col='blue',lwd=5,xlab='period (days)',ylab='',main='',cex.lab=1.5,cex.axis=1.5,cex.main=1.8,xlim=c(2,100),ylim=c(-0.25,0.75),xaxt='n',yaxt='n',log='x')
      for (z in 1:1000) {
         lines(x1,rev(ycor30med[,1]+(ycor30hi[,1]-ycor30med[,1])*z/1000),lwd=1,col='lightblue')
         lines(x1,rev(ycor30med[,1]+(ycor30lo[,1]-ycor30med[,1])*z/1000),lwd=1,col='lightblue')
         }
      lines(x1,rev(ycor30med[,1]),lwd=5,col='blue')
      lines(x1,rev(ycor30[,1]),lwd=5,col='gray')
      lines(x1,rev(ycor[,1]),lwd=5,col='black')
      text(35,0.73,'1981-2010',cex=1.5,col='gray',pos=4)
      text(35,0.67,'1901-2015',cex=1.5,col='black',pos=4)
      text(25,0.59,'30-year median',cex=1.5,col='blue',pos=4)
      text(25,0.53,'30-year IQR',cex=1.5,col='lightblue',pos=4)
      text(3.5,0.73,'a) Central',cex=2)
      lines(c(1,120),rep(0.183,2),lty=2,lwd=2,col='black')
      lines(c(1,120),rep(0.361,2),lty=2,lwd=2,col='gray')
      lines(c(1,120),rep(-0.183,2),lty=2,lwd=2,col='black')
      axis(1,at=c(2,3,4,6,10,15,30,60,90),col.axis='black',las=1,cex.axis=1.5)
      axis(2,at=c(-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7),col.axis='black',las=1,cex.axis=1.5)
      if (!lplotfile) dev.off()

      if (!lplotfile) pdf('Fig1b.pdf')
      plot(x1,rev(ycor30med[,2]),type='l',col='blue',lwd=5,xlab='period (days)',ylab='',main='',cex.lab=1.5,cex.axis=1.5,cex.main=1.8,xlim=c(2,100),ylim=c(-0.25,0.75),xaxt='n',yaxt='n',log='x')
      for (z in 1:1000) {
         lines(x1,rev(ycor30med[,2]+(ycor30hi[,2]-ycor30med[,2])*z/1000),lwd=1,col='lightblue')
         lines(x1,rev(ycor30med[,2]+(ycor30lo[,2]-ycor30med[,2])*z/1000),lwd=1,col='lightblue')
         }
      lines(x1,rev(ycor30med[,2]),lwd=5,col='blue')
      lines(x1,rev(ycor30[,2]),lwd=5,col='gray')
      lines(x1,rev(ycor[,2]),lwd=5,col='black')
      text(35,0.73,'1981-2010',cex=1.5,col='gray',pos=4)
      text(35,0.67,'1901-2015',cex=1.5,col='black',pos=4)
      text(25,0.59,'30-year median',cex=1.5,col='blue',pos=4)
      text(25,0.53,'30-year IQR',cex=1.5,col='lightblue',pos=4)
      text(3.8,0.73,'b) All-India',cex=2)
      lines(c(1,120),rep(0.183,2),lty=2,lwd=2,col='black')
      lines(c(1,120),rep(0.361,2),lty=2,lwd=2,col='gray')
      lines(c(1,120),rep(-0.183,2),lty=2,lwd=2,col='black')
      axis(1,at=c(2,3,4,6,10,15,30,60,90),col.axis='black',las=1,cex.axis=1.5)
      axis(2,at=c(-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7),col.axis='black',las=1,cex.axis=1.5)
      if (!lplotfile) dev.off()

      cols = c('black','#FF9933','#138808','#000080')
      pdf('Fig2b.pdf')
      plot(x1,rev(ycor[,2]),type='l',lwd=5,col='#FF9933',xlab='period (days)',ylab='',main='',cex.lab=1.5,cex.axis=1.5,cex.main=1.8,xlim=c(2,100),ylim=c(-0.25,0.75),xaxt='n',yaxt='n',log='x')
      lines(x1,rev(ycor[,1]),lwd=5,col='black')
      lines(x1,rev(ycor[,3]),lwd=5,col='#138808')
      lines(c(1,120),rep( 0.183,2),lty=2,col='black')
      lines(c(1,120),rep(-0.183,2),lty=2,col='black')
      text(45,0.49,'All-India',cex=1.5,col='#138808')
      text(45,0.44,'without Central',cex=1.5,col='#138808')
      text(7,0.1,'Central',cex=1.5,col='black')
      text(45,0.1,'All-India',cex=1.5,col='#FF9933')
      text(3.8,0.73,'b) correlation',cex=2)
      axis(1,at=c(2,3,4,6,10,15,30,60,90),col.axis='black',las=1,cex.axis=1.5)
      axis(2,at=c(-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8),col.axis='black',las=1,cex.axis=1.5)
      dev.off()
      }
   cat(paste(round(100*k/nreg),'%\n',sep=''))
   }


# Central India climatological cumulative variance
yvar = apply(yvars[,iyrs1,1],1,mean)

# Average results over 4 sub-regions
ycorave <- array(0,dim=c(nm,9,2))
yvarave <- array(0,dim=c(nm,9))
for (z in 1:9) {
   ii = 3 + 4*(z-1) + c(1:4)
   yvarave[,z]   = apply(yvars[,iyrs1,ii],1,mean)
   ycorave[,z,1] = apply(ycor[,ii],1,mean)
   ycorave[,z,2] = apply(ycor1[,ii],1,mean)
   }


# Fig. 4a, 4b & 4c for cumulative variance

cols = c('black','darkblue','blue','lightblue3','lightblue','bisque2','orange','red1','red3','darkred')
if (!lplotfile) pdf('Fig4a.pdf')
y1 = rev(yvar); ylim = c(0,65)
plot(x1,y1,type='l',lwd=3,xlab='period (days)',ylab='',main='',cex.lab=1.5,cex.axis=1.5,cex.main=1.5,xlim=c(2,150),ylim=ylim,xaxt='n',yaxt='n',log='x')
text(130,y1[length(y1)],paste(100,'%',sep=''),cex=1.5,col=cols[1])
for (z in 1:9) {
   y1 = rev(yvarave[,z])
   lines(x1,y1,lwd=3,col=cols[1+z])
   text(130,y1[length(y1)],paste((1-0.1*z)*100,'%',sep=''),cex=1.5,col=cols[1+z])
   }
y1 = rev(yvarave[,1])
lines(x1,y1,lwd=3)
axis(1,at=c(2,3,4,6,10,15,30,60,90),col.axis='black',las=1,cex.axis=1.5)
axis(2,at=c(0,10,20,30,40,50,60,70),col.axis='black',las=1,cex.axis=1.5)
text(1.75,ylim[2]-(ylim[2]-ylim[1])/40,expression(paste('a) Cumulative variance (mm'^2,'/day'^2,')',sep='')),cex=1.5,pos=4)
if (!lplotfile) dev.off()

if (!lplotfile) pdf('Fig4b.pdf')
y1 = rev(ycor[,1]); ylim = c(-0.25,0.75)
plot(x1,y1,type='l',lwd=3,xlab='period (days)',ylab='',main='',cex.lab=1.5,cex.axis=1.5,cex.main=1.8,xlim=c(2,150),ylim=ylim,xaxt='n',yaxt='n',log='x')
for (z in 1:9) {
   y1 = rev(ycorave[,z,2])
   lines(x1,y1,lwd=3,col=cols[1+z])
   text(130,ylim[2]-z*(ylim[2]-ylim[1])/20,paste((0.1*z)*100,'%',sep=''),cex=1.5,col=cols[11-z])
   }
text(130,ylim[2]-10*(ylim[2]-ylim[1])/20,'100%',cex=1.5,col='black')
y1 = rev(ycor[,1])
lines(x1,y1,lwd=3)
lines(c(1,90),rep(0.183,2),lty=2)
lines(c(1,90),rep(-0.183,2),lty=2)
axis(1,at=c(2,3,4,6,10,15,30,60,90),col.axis='black',las=1,cex.axis=1.5)
axis(2,at=c(-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8),col.axis='black',las=1,cex.axis=1.5)
text(1.75,ylim[2]-(ylim[2]-ylim[1])/40,'b) Central India correlation (cumulative)',cex=1.5,pos=4)
if (!lplotfile) dev.off()

if (!lplotfile) pdf('Fig4c.pdf')
y1 = rev(ycor[,1]); ylim = c(-0.25,0.75)
plot(x1,y1,type='l',lwd=3,xlab='period (days)',ylab='',main='',cex.lab=1.5,cex.axis=1.5,cex.main=1.8,xlim=c(2,150),ylim=ylim,xaxt='n',yaxt='n',log='x')
text(130,y1[length(y1)],paste(100,'%',sep=''),cex=1.5,col=cols[1])
for (z in 1:9) {
   y1 = rev(ycorave[,z,1])
   lines(x1,y1,lwd=3,col=cols[1+z])
   text(130,y1[length(y1)],paste((1-0.1*z)*100,'%',sep=''),cex=1.5,col=cols[1+z])
   }
y1 = rev(ycor[,1])
lines(x1,y1,lwd=3)
lines(c(1,90),rep(0.183,2),lty=2)
lines(c(1,90),rep(-0.183,2),lty=2)
axis(1,at=c(2,3,4,6,10,15,30,60,90),col.axis='black',las=1,cex.axis=1.5)
axis(2,at=c(-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8),col.axis='black',las=1,cex.axis=1.5)
text(1.75,ylim[2]-(ylim[2]-ylim[1])/40,'c) Local correlation (cumulative)',cex=1.5,pos=4)
if (!lplotfile) dev.off()


# Fig. 4d, 4e & 4f for running band variance

n = 39; logbreaks <- array(0,dim=c(2,n))
for (z in 1:n) logbreaks[1,z] = (z-1)*0.1 + 0.5
logbreaks[2,] = logbreaks[1,] + 0.5
breaks = exp(logbreaks)
nharm <- numeric(); mids <- numeric(); iband <- list()
last <- numeric()
for (z in 1:n) {
   mids[z] = exp(mean(logbreaks[,z]))
   iband[[z]] = which(x1>=breaks[1,z]&x1<breaks[2,z])
   nharm[z] = length(iband[[z]])
   last[z] = x1[rev(iband[[z]])[1]]
   }

# Back out individual contributions to variance from cumulative
yvarback = yvars-yvars
yvarback[1,,] = yvars[nm,,]
for (z in 2:nm) yvarback[z,,] = yvars[nm+1-z,,]-yvars[nm+2-z,,]

# Running band variances
yruns <- array(0,dim=c(n,nyrs,nreg))
for (k in 1:nreg) for (j in 1:nyrs) for (z in 1:n) {
   yruns[z,j,k] = sum(yvarback[iband[[z]],j,k])
   }

# Central India climatological running band variance
yrun = apply(yruns[,iyrs1,1],1,mean)

ycorrun <- array(0,dim=c(n,nreg))
ycorrun1 = ycorrun

# Calculate local correlation and fixed Central India correlation
for (k in 1:nreg) for (i in 1:n) {
    ycorrun[i,k] = cor(ymean[iyrs1,k],yruns[i,iyrs1,k])
   ycorrun1[i,k] = cor(ymean[iyrs1,1],yruns[i,iyrs1,k])
   }

# Average results over 4 sub-regions
ycorrunave <- array(0,dim=c(n,9,2))
yrunave <- array(0,dim=c(n,9))
for (z in 1:9) {
   ii = 3 + 4*(z-1) + c(1:4)
   yrunave[,z]   = apply(yruns[,iyrs1,ii],1,mean)
   ycorrunave[,z,1] = apply(ycorrun[,ii],1,mean)
   ycorrunave[,z,2] = apply(ycorrun1[,ii],1,mean)
   }

x1run = mids

ylim = c(0,16)
pdf('Fig4d.pdf')
plot(x1run,yrun,log='x',type='l',lwd=3,ylim=ylim,xlab='period (days)',ylab='',main='',cex.lab=1.5,cex.axis=1.5,cex.main=1.5,xaxt='n',yaxt='n',xlim=c(2,150))
for (zz in 1:9) {
   lines(x1run,yrunave[,zz],lwd=3,col=cols[1+zz])
   text(130,ylim[2]-zz*(ylim[2]-ylim[1])/20,paste((0.1*zz)*100,'%',sep=''),cex=1.5,col=cols[11-zz])
   }
text(130,ylim[2]-10*(ylim[2]-ylim[1])/20,'100%',cex=1.5,col='black')
axis(1,at=c(2,3,4,6,10,15,30,60,90),col.axis='black',las=1,cex.axis=1.5)
axis(2,at=c(0,2,4,6,8,10,12,14,16),col.axis='black',las=1,cex.axis=1.5)
text(1.75,ylim[2]-(ylim[2]-ylim[1])/40,expression(paste('d) Running band variance (mm'^2,'/day'^2,')',sep='')),cex=1.5,pos=4)
dev.off()

pdf('Fig4e.pdf')
y1 = ycorrun[,1]; ylim = c(-0.25,0.75)
plot(x1run,ycorrun1[,1],type='l',lwd=3,xlab='period (days)',ylab='',main='',cex.lab=1.5,cex.axis=1.5,cex.main=1.8,xlim=c(2,150),ylim=ylim,xaxt='n',yaxt='n',log='x')
for (z in 1:9) {
   lines(x1run,ycorrunave[,z,2],lwd=3,col=cols[1+z])
   text(130,ylim[2]-z*(ylim[2]-ylim[1])/20,paste((0.1*z)*100,'%',sep=''),cex=1.5,col=cols[11-z])
   }
text(130,ylim[2]-10*(ylim[2]-ylim[1])/20,'100%',cex=1.5,col='black')
lines(x1run,ycorrun1[,1],lwd=3)
lines(c(1,90),rep(0.183,2),lty=2)
lines(c(1,90),rep(-0.183,2),lty=2)
axis(1,at=c(2,3,4,6,10,15,30,60,90),col.axis='black',las=1,cex.axis=1.5)
axis(2,at=c(-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8),col.axis='black',las=1,cex.axis=1.5)
text(1.75,ylim[2]-(ylim[2]-ylim[1])/40,'e) Central India correlation (running band)',cex=1.5,pos=4)
dev.off()

pdf('Fig4f.pdf')
ylim = c(-0.25,0.75)
plot(x1run,ycorrun[,1],type='l',lwd=3,xlab='period (days)',ylab='',main='',cex.lab=1.5,cex.axis=1.5,cex.main=1.8,xlim=c(2,150),ylim=ylim,xaxt='n',yaxt='n',log='x')
for (z in 1:9) {
   lines(x1run,ycorrunave[,z,1],lwd=3,col=cols[1+z])
   text(130,ylim[2]-z*(ylim[2]-ylim[1])/20,paste((0.1*z)*100,'%',sep=''),cex=1.5,col=cols[11-z])
   }
text(130,ylim[2]-10*(ylim[2]-ylim[1])/20,'100%',cex=1.5,col='black')
lines(x1run,ycorrun[,1],lwd=3)
lines(c(1,90),rep(0.183,2),lty=2)
lines(c(1,90),rep(-0.183,2),lty=2)
axis(1,at=c(2,3,4,6,10,15,30,60,90),col.axis='black',las=1,cex.axis=1.5)
axis(2,at=c(-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8),col.axis='black',las=1,cex.axis=1.5)
text(1.75,ylim[2]-(ylim[2]-ylim[1])/40,'f) Local correlation (running band)',cex=1.5,pos=4)
dev.off()


