
# Input:
#   binary IMD rainfall data for 1901-2015 (must be in same directory)

# Output:
#   binary mask data and GrADS .ctl file for producing Figs. 2a & 3
#      (Fig23.mask.dat,Fig23.mask.ctl)
#         - these files along with GrADS scripts Fig2a.gs & Fig3.gs
#           are contained in 1_GrADS_Figs2a3.tar
#   binary areal-average rainfall data, required for producing Figs. 1 and 2b
#      (avgimd365.dat)  
#   remove stop() statement on line 141 for producing Fig. 4 panels
#      (Fig4a.pdf,Fig4b.pdf,Fig4c.pdf,Fig4d.pdf,Fig4e.pdf,Fig4f.pdf)

rm(list=ls())

ntall = 365 # days in year
nyrs  = 115 # years

# Spatial domain
nlon1 = 135; lon1 = (c(1:nlon1)-1)*0.25 + 66.5
nlat1 = 129; lat1 = (c(1:nlat1)-1)*0.25 +  6.5

ilonsub <- list()
ilatsub <- list()

# All regions for spatial averaging
nreg = 1+4*9+2
# Central India
x1 = 74.5; x2 = 86.5; xdel = x2-x1
y1 = 16.5; y2 = 26.5; ydel = y2-y1
ilonsub[[1]] = which(lon1>=x1&lon1<=x2)
ilatsub[[1]] = which(lat1>=y1&lat1<=y2)
# Central India split up by 2 by 2 quadrants
# that are (100-10*z) % the size
for (z in 1:9) {
   xdel1 = xdel*sqrt(1.0-0.1*z)
   ydel1 = ydel*sqrt(1.0-0.1*z)
   ilonsub[[1+4*(z-1)+1]] = which(lon1>=x1&lon1<(x1+xdel1))
   ilatsub[[1+4*(z-1)+1]] = which(lat1>=y1&lat1<(y1+ydel1))
   ilonsub[[1+4*(z-1)+2]] = which(lon1>=(x2-xdel1)&lon1<=x2)
   ilatsub[[1+4*(z-1)+2]] = ilatsub[[1+4*(z-1)+1]]
   ilonsub[[1+4*(z-1)+3]] = ilonsub[[1+4*(z-1)+1]]
   ilatsub[[1+4*(z-1)+3]] = which(lat1>=(y2-ydel1)&lat1<=y2)
   ilonsub[[1+4*(z-1)+4]] = ilonsub[[1+4*(z-1)+2]]
   ilatsub[[1+4*(z-1)+4]] = ilatsub[[1+4*(z-1)+3]]
# Adjustment of SE box to account for undefined Bay of Bengal
   if ((1.0-0.1*z)<0.25) {
      ilonsub[[1+4*(z-1)+2]] = which(lon1>=((x1+x2)/2)&lon1<=((x1+x2)/2+xdel1))
      ilatsub[[1+4*(z-1)+2]] = which(lat1>=((y1+y2)/2-ydel1)&lat1<=((y1+y2)/2))
      }
   }
# All-India
ilonsub[[nreg-1]] = which(lon1>=66.5&lon1<=100.5)
ilatsub[[nreg-1]] = which(lat1>= 6.5&lat1<=32.5)

# Boxes defining regions
box <- array(-999,dim=c(nlon1,nlat1,nreg))
for (z in 1:(nreg-1)) {
   box[ilonsub[[z]],ilatsub[[z]],z] = 1
   box[ilonsub[[z]],ilatsub[[z]],z] = 1
   }

# Manually define All-India w/o Central India
box[ilonsub[[nreg-1]],ilatsub[[nreg-1]],nreg] = 1
box[ilonsub[[nreg-1]],ilatsub[[nreg-1]],nreg] = 1
box[ilonsub[[1]],ilatsub[[1]],nreg] = -999
box[ilonsub[[1]],ilatsub[[1]],nreg] = -999

# Latitude weighting
g1 = rep(cos(lat1*pi/180),each=nlon1)

# Read in data to get location of Indian land grid points
year = 1901; n = nlon1*nlat1*ntall
filename = paste('ind',year,'_rfp25.grd',sep='')
temp = readBin(filename,'numeric',n=n,size='4',endian='little')
dim(temp) = c(nlon1,nlat1,ntall)

# Find where defined
india <- array(-999,dim=c(nlon1,nlat1))
def = temp[,,1]!=-999
idef = which(def)
india[idef] = 1

# Masks assigned latitude weighting only where defined
mask <- array(0,dim=c(nlon1*nlat1,nreg))
for (z in 1:nreg) {
   idef = which(temp[,,1]!=-999&box[,,z]!=-999)
   mask[idef,z] = rep(g1,each=nlon1)[idef]
   }

# Write out masks for producing Figs. 2a & 3
nbuf = 20
varout <- array(0,dim=c(nbuf+nlon1+nbuf,nbuf+nlat1+nbuf,4*2+4*9+4*2))
varout[(nbuf+1):(nbuf+nlon1),(nbuf+1):(nbuf+nlat1),1:4] = rep(india,4)
varout[(nbuf+1):(nbuf+nlon1),(nbuf+1):(nbuf+nlat1),5:8] = rep(box[,,1],4)
varout[(nbuf+1):(nbuf+nlon1),(nbuf+1):(nbuf+nlat1),9:(8+4*9)] = box[,,2:(1+4*9)]
varout[(nbuf+1):(nbuf+nlon1),(nbuf+1):(nbuf+nlat1),45:48] = rep(box[,,nreg-1],4)
varout[(nbuf+1):(nbuf+nlon1),(nbuf+1):(nbuf+nlat1),49:52] = rep(box[,,nreg],4)
varout[varout==-999] = 0
writeBin(as.vector(varout),'Fig23.mask.dat',size='4',endian='little')
system("echo 'DSET ^Fig23.mask.dat\nTITLE 0.25 degr grid with different region masks\nUNDEF -999\nXDEF 175  LINEAR  61.5 0.25\nYDEF 169  LINEAR  1.5 0.25\nZDEF   1 linear 0 1\nTDEF   1 LINEAR 1jun1979 1yr\nEDEF   13 names india cen 10 20 30 40 50 60 70 80 90 all notcent\nVARS   4\nbox1   0 99\nbox2   0 99\nbox3   0 99\nbox4   0 99\nENDVARS' > Fig23.mask.ctl")
system("unix2dos Fig23.mask.ctl")


# Normalization so that mask sums are areal averages
for (z in 1:nreg) mask[,z] = mask[,z]/sum(mask[,z])

mask = rep(mask,ntall)
dim(mask) = c(nlon1*nlat1,nreg,ntall)

prec <- array(0,dim=c(ntall,nyrs,nreg))
for (l in 1:nyrs) {
   year = 1900 + l
   nt0 = 365; iref = 0
# During leap year, begin 02Jan
   isleap = (l/4)==round(l/4)
   if (isleap) {nt0 = nt0+1; iref = iref+1}

# Read in IMD rainfall 
   n = nlon1*nlat1*nt0
   filename = paste('ind',year,'_rfp25.grd',sep='')
   temp = readBin(filename,'numeric',n=n,size='4',endian='little')
   dim(temp) = c(nlon1*nlat1,nt0)
   temp = temp[,iref + c(1:ntall)]
 
# Take areal averages
   for (z in 1:nreg) prec[,l,z] = apply(temp*mask[,z,],2,sum)
   print(year)
   }

# Write out areal-average IMD rainfall
varout <- array(0,dim=c(ntall,nyrs,7))
varout[,,1] = prec[,,1]       # Central India
varout[,,2] = prec[,,nreg-1]  # All-India
varout[,,3] = prec[,,nreg]    # All-India w/o Central India
varout[,,4:7] = prec[,,34:37] # Small scale within Central India
writeBin(as.vector(varout),'avgimd365.dat',size='4',endian='little')


stop()


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

# Calculate harmonics
comp <- array(0,dim=c(ntall,1+nmall,nyrs,nreg))
for (l in 1:nyrs) {
   year = 1900 + l
   for (z in 1:nreg) {
      comp[,,l,z] = fftcomp(prec[,l,z])
      }
   print(year)
   }

nt = 120   # summer (01Jun-28Sep,JJAS)
tref = 151 # index for 30Apr

# Seasonal means
ymean = apply(prec[tref+c(1:nt),,],c(2,3),mean)

y <- array(0,dim=c(nmall-3,nyrs,nreg))
ycor <- array(0,dim=c(nmall-3,nreg))
ycor1 <- array(0,dim=c(nmall-3,9))
yvar1 = ycor1
ycor2 = ycor1
yvar1yr <- array(0,dim=c(nmall-3,9,nyrs))

# Cumulative variances
for (k in 1:nreg) for (j in 1:nyrs) for (i in 1:(nmall-4)) {
   temp = apply(comp[tref+c(1:nt),(3+i):nmall,j,k,drop=FALSE],1,sum)
   y[i,j,k] = mean((temp-mean(temp))^2)
   }
# Cumulative variance local correlation
for (k in 1:nreg) for (i in 1:(nmall-4)) {
   ycor[i,k] = cor(ymean[,k],y[i,,k])
   }

# Climatological cumulative variance
yvar = apply(y[,,1],1,mean)

# Cumulative variance Central India correlation
for (z in 1:9) for (i in 1:(nmall-4)) {
   ii = 1 + 4*(z-1) + c(1:4)
   temp = apply(y[i,,ii],1,mean)
   yvar1[i,z] = mean(temp)
   yvar1yr[i,z,] = temp
   ycor1[i,z] = cor(ymean[,1],temp)
   for (zz in 1:4) ycor2[i,z] = ycor2[i,z] + cor(ymean[,1],y[i,,ii[zz]])/4
   }


cols = c('black','darkblue','blue','lightblue3','lightblue','bisque2','orange','red1','red3','darkred')

# periods
x1 = rev(365/c(4:nmall))

# Fig. 4a, 4b & 4c for cumulative variance
pdf('Fig4a.pdf')
y1 = rev(yvar); ylim = c(0,65)
plot(x1,y1,type='l',lwd=3,xlab='period (days)',ylab='',main='',cex.lab=1.5,cex.axis=1.5,cex.main=1.5,xlim=c(2,150),ylim=ylim,xaxt='n',yaxt='n',log='x')
text(130,y1[length(y1)],paste(100,'%',sep=''),cex=1.5,col=cols[1])
for (z in 1:9) {
   y1 = rev(yvar1[,z])
   lines(x1,y1,lwd=3,col=cols[1+z])
   text(130,y1[length(y1)],paste((1-0.1*z)*100,'%',sep=''),cex=1.5,col=cols[1+z])
   }
y1 = rev(yvar1[,1])
lines(x1,y1,lwd=3)
axis(1,at=c(2,3,4,6,10,15,30,60,90),col.axis='black',las=1,cex.axis=1.5)
axis(2,at=c(0,10,20,30,40,50,60,70),col.axis='black',las=1,cex.axis=1.5)
text(1.75,ylim[2]-(ylim[2]-ylim[1])/40,expression(paste('a) Cumulative variance (mm'^2,'/day'^2,')',sep='')),cex=1.5,pos=4)
dev.off()

pdf('Fig4b.pdf')
y1 = rev(ycor[,1]); ylim = c(-0.25,0.75)
plot(x1,y1,type='l',lwd=3,xlab='period (days)',ylab='',main='',cex.lab=1.5,cex.axis=1.5,cex.main=1.8,xlim=c(2,150),ylim=ylim,xaxt='n',yaxt='n',log='x')
for (z in 1:9) {
   y1 = rev(ycor2[,z])
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
text(1.75,ylim[2]-(ylim[2]-ylim[1])/40,'b) Central India correlation',cex=1.5,pos=4)
dev.off()

pdf('Fig4c.pdf')
y1 = rev(ycor[,1]); ylim = c(-0.25,0.75)
plot(x1,y1,type='l',lwd=3,xlab='period (days)',ylab='',main='',cex.lab=1.5,cex.axis=1.5,cex.main=1.8,xlim=c(2,150),ylim=ylim,xaxt='n',yaxt='n',log='x')
text(130,y1[length(y1)],paste(100,'%',sep=''),cex=1.5,col=cols[1])
for (z in 1:9) {
   ii = 1 + 4*(z-1) + c(1:4)
   y1 = rev(apply(ycor[,ii],1,mean,na.rm=TRUE))
   lines(x1,y1,lwd=3,col=cols[1+z])
   text(130,y1[length(y1)],paste((1-0.1*z)*100,'%',sep=''),cex=1.5,col=cols[1+z])
   }
y1 = rev(ycor[,1])
lines(x1,y1,lwd=3)
lines(c(1,90),rep(0.183,2),lty=2)
lines(c(1,90),rep(-0.183,2),lty=2)
axis(1,at=c(2,3,4,6,10,15,30,60,90),col.axis='black',las=1,cex.axis=1.5)
axis(2,at=c(-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8),col.axis='black',las=1,cex.axis=1.5)
text(1.75,ylim[2]-(ylim[2]-ylim[1])/40,'c) Local correlation',cex=1.5,pos=4)
dev.off()


# Fig. 4d, 4e & 4f for running band variance

var1 <- numeric(); var1yr <- array(0,dim=c(nmall-3,nyrs))
var2 <- array(0,dim=c(nmall-3,9)); var2yr <- array(0,dim=c(nmall-3,9,nyrs))
var1[1] = rev(yvar)[1]
yvaryr = y[,,1]
var1yr[1,] = apply(yvaryr,2,rev)[1,]
# Climatological variance calculation
for (zz in 1:9) var2[1,zz] = rev(yvar1[,zz])[1]
# Back out individual contributions to variance from cumulative
for (z in 2:(nmall-3)) {
   var1[z] = rev(yvar)[z]-rev(yvar)[z-1]
   for (zz in 1:9) var2[z,zz] = rev(yvar1[,zz])[z]-rev(yvar1[,zz])[z-1]
   }
# Variance calculation for each year separately
for (zz in 1:9) var2yr[1,zz,] = apply(yvar1yr[,zz,],2,rev)[1,]
# Back out individual contributions to variance from cumulative
for (z in 2:(nmall-3)) {
   var1yr[z,] = apply(yvaryr,2,rev)[z,] - apply(yvaryr,2,rev)[z-1,]
   for (zz in 1:9) var2yr[z,zz,] = apply(yvar1yr[,zz,],2,rev)[z,] - apply(yvar1yr[,zz,],2,rev)[z-1,]
   }


# Running band window
n = 39; logbreaks <- array(0,dim=c(2,n))
for (z in 1:n) logbreaks[1,z] = (z-1)*0.1 + 0.5
logbreaks[2,] = logbreaks[1,] + 0.5
breaks = exp(logbreaks)
x1 = rev(365/c(4:nmall))
nharm <- numeric(); mids <- numeric(); iband <- list()
last <- numeric()
for (z in 1:n) {
   mids[z] = exp(mean(logbreaks[,z]))
   iband[[z]] = which(x1>=breaks[1,z]&x1<breaks[2,z])
   nharm[z] = length(iband[[z]])
   last[z] = x1[rev(iband[[z]])[1]]
   }

# Back out individual contributions to variance from cumulative
yback = y-y
yback[1,,] = y[nmall-3,,]
for (z in 2:(nmall-3)) yback[z,,] = y[nmall-3+1-z,,]-y[nmall-3+2-z,,]

# Running band variances
yrun <- array(0,dim=c(n,nyrs,nreg))
for (k in 1:nreg) for (j in 1:nyrs) for (z in 1:n) {
   yrun[z,j,k] = sum(yback[iband[[z]],j,k])
   }

ycorrun <- array(0,dim=c(n,nreg))
ycorrun1 <- array(0,dim=c(n,9))

# Running band variance local correlation
for (k in 1:nreg) for (i in 1:n) {
   ycorrun[i,k] = cor(ymean[,k],yrun[i,,k])
   }
# Running band variance CI correlation
for (z in 1:9) for (i in 1:n) {
   ii = 1 + 4*(z-1) + c(1:4)
   for (zz in 1:4) ycorrun1[i,z] = ycorrun1[i,z] + cor(ymean[,1],yrun[i,,ii[zz]])/4
   }

# Climatological running band variance calculation
var3 <- numeric()
var4 <- array(0,dim=c(n,9))
for (z in 1:n) {
   var3[z] = sum(var1[iband[[z]]])
   for (zz in 1:9) var4[z,zz] = sum(var2[iband[[z]],zz])
   }

# Running band variance calculation for each year separately
var3yr <- array(0,dim=c(n,nyrs))
var4yr <- array(0,dim=c(n,9,nyrs))
for (z in 1:n) {
   var3yr[z,] = apply(var1yr[iband[[z]],,drop=FALSE],2,sum)
   for (zz in 1:9) var4yr[z,zz,] = apply(var2yr[iband[[z]],zz,,drop=FALSE],3,sum)
   }

# Check
#var3 == apply(var3yr,1,mean) 
#var4 == apply(var4yr,c(1,2),mean)

x1run = mids

ylim = c(0,16)
pdf('Fig4d.pdf')
plot(x1run,var3,log='x',type='l',lwd=3,ylim=ylim,xlab='period (days)',ylab='',main='',cex.lab=1.5,cex.axis=1.5,cex.main=1.5,xaxt='n',yaxt='n',xlim=c(2,150))
for (zz in 1:9) {
   lines(x1run,var4[,zz],lwd=3,col=cols[1+zz])
   text(130,ylim[2]-zz*(ylim[2]-ylim[1])/20,paste((0.1*zz)*100,'%',sep=''),cex=1.5,col=cols[11-zz])
   }
text(130,ylim[2]-10*(ylim[2]-ylim[1])/20,'100%',cex=1.5,col='black')
axis(1,at=c(2,3,4,6,10,15,30,60,90),col.axis='black',las=1,cex.axis=1.5)
axis(2,at=c(0,2,4,6,8,10,12,14,16),col.axis='black',las=1,cex.axis=1.5)
text(1.75,ylim[2]-(ylim[2]-ylim[1])/40,expression(paste('d) Running band variance (mm'^2,'/day'^2,')',sep='')),cex=1.5,pos=4)
dev.off()

pdf('Fig4e.pdf')
y1 = ycorrun[,1]; ylim = c(-0.25,0.75)
plot(x1run,y1,type='l',lwd=3,xlab='period (days)',ylab='',main='',cex.lab=1.5,cex.axis=1.5,cex.main=1.8,xlim=c(2,150),ylim=ylim,xaxt='n',yaxt='n',log='x')
for (z in 1:9) {
   y1 = ycorrun1[,z]
   lines(x1run,y1,lwd=3,col=cols[1+z])
   text(130,ylim[2]-z*(ylim[2]-ylim[1])/20,paste((0.1*z)*100,'%',sep=''),cex=1.5,col=cols[11-z])
   }
text(130,ylim[2]-10*(ylim[2]-ylim[1])/20,'100%',cex=1.5,col='black')
y1 = ycorrun[,1]
lines(x1run,y1,lwd=3)
lines(c(1,90),rep(0.183,2),lty=2)
lines(c(1,90),rep(-0.183,2),lty=2)
axis(1,at=c(2,3,4,6,10,15,30,60,90),col.axis='black',las=1,cex.axis=1.5)
axis(2,at=c(-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8),col.axis='black',las=1,cex.axis=1.5)
text(1.75,ylim[2]-(ylim[2]-ylim[1])/40,'e) Central India correlation',cex=1.5,pos=4)
dev.off()

pdf('Fig4f.pdf')
y1 = ycorrun[,1]; ylim = c(-0.25,0.75)
plot(x1run,y1,type='l',lwd=3,xlab='period (days)',ylab='',main='',cex.lab=1.5,cex.axis=1.5,cex.main=1.8,xlim=c(2,150),ylim=ylim,xaxt='n',yaxt='n',log='x')
for (z in 1:9) {
   ii = 1 + 4*(z-1) + c(1:4)
   y1 = apply(ycorrun[,ii],1,mean)
   lines(x1run,y1,lwd=3,col=cols[1+z])
   text(130,ylim[2]-z*(ylim[2]-ylim[1])/20,paste((0.1*z)*100,'%',sep=''),cex=1.5,col=cols[11-z])
   }
text(130,ylim[2]-10*(ylim[2]-ylim[1])/20,'100%',cex=1.5,col='black')
y1 = ycorrun[,1]
lines(x1run,y1,lwd=3)
lines(c(1,90),rep(0.183,2),lty=2)
lines(c(1,90),rep(-0.183,2),lty=2)
axis(1,at=c(2,3,4,6,10,15,30,60,90),col.axis='black',las=1,cex.axis=1.5)
axis(2,at=c(-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8),col.axis='black',las=1,cex.axis=1.5)
text(1.75,ylim[2]-(ylim[2]-ylim[1])/40,'f) Local correlation',cex=1.5,pos=4)
dev.off()


