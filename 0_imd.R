
# Input:
#   binary IMD rainfall data for 1901-2015 (must be in same directory)

# Output:
#   binary mask data and GrADS .ctl file for producing Figs. 2a & 3
#      (Fig23.mask.dat,Fig23.mask.ctl)
#         - these files along with GrADS scripts Fig2a.gs & Fig3.gs
#           are contained in 1_GrADS_Figs2a3.tar
#   binary areal-average rainfall data, required for producing Figs. 1,2b,4
#      (avgimd365.dat)  

rm(list=ls())

ntall = 365 # days in year
nyrs  = 115 # years

# Spatial domain
nlon1 = 135; lon1 = (c(1:nlon1)-1)*0.25 + 66.5
nlat1 = 129; lat1 = (c(1:nlat1)-1)*0.25 +  6.5

ilonsub <- list()
ilatsub <- list()

# All regions for spatial averaging
nreg = 3+4*9
# Central India
x1 = 74.5; x2 = 86.5; xdel = x2-x1
y1 = 16.5; y2 = 26.5; ydel = y2-y1
ilonsub[[1]] = which(lon1>=x1&lon1<=x2)
ilatsub[[1]] = which(lat1>=y1&lat1<=y2)
# All-India
ilonsub[[2]] = which(lon1>=66.5&lon1<=100.5)
ilatsub[[2]] = which(lat1>= 6.5&lat1<=32.5)
# All-India (later remove Central India)
ilonsub[[3]] = which(lon1>=66.5&lon1<=100.5)
ilatsub[[3]] = which(lat1>= 6.5&lat1<=32.5)
# Central India split up by 2 by 2 quadrants
# that are (100-10*z) % the size
for (z in 1:9) {
   xdel1 = xdel*sqrt(1.0-0.1*z)
   ydel1 = ydel*sqrt(1.0-0.1*z)
   ilonsub[[3+4*(z-1)+1]] = which(lon1>=x1&lon1<(x1+xdel1))
   ilatsub[[3+4*(z-1)+1]] = which(lat1>=y1&lat1<(y1+ydel1))
   ilonsub[[3+4*(z-1)+2]] = which(lon1>=(x2-xdel1)&lon1<=x2)
   ilatsub[[3+4*(z-1)+2]] = ilatsub[[3+4*(z-1)+1]]
   ilonsub[[3+4*(z-1)+3]] = ilonsub[[3+4*(z-1)+1]]
   ilatsub[[3+4*(z-1)+3]] = which(lat1>=(y2-ydel1)&lat1<=y2)
   ilonsub[[3+4*(z-1)+4]] = ilonsub[[3+4*(z-1)+2]]
   ilatsub[[3+4*(z-1)+4]] = ilatsub[[3+4*(z-1)+3]]
# Adjustment of SE box to account for undefined Bay of Bengal
   if ((1.0-0.1*z)<0.25) {
      ilonsub[[3+4*(z-1)+2]] = which(lon1>=((x1+x2)/2)&lon1<=((x1+x2)/2+xdel1))
      ilatsub[[3+4*(z-1)+2]] = which(lat1>=((y1+y2)/2-ydel1)&lat1<=((y1+y2)/2))
      }
   }
# All-India
ilonsub[[2]] = which(lon1>=66.5&lon1<=100.5)
ilatsub[[2]] = which(lat1>= 6.5&lat1<=32.5)

# Boxes defining regions
box <- array(-999,dim=c(nlon1,nlat1,nreg))
for (z in 1:nreg) {
   box[ilonsub[[z]],ilatsub[[z]],z] = 1
   box[ilonsub[[z]],ilatsub[[z]],z] = 1
   }

# Remove Central India for All-India without Central India
box[ilonsub[[1]],ilatsub[[1]],3] = -999
box[ilonsub[[1]],ilatsub[[1]],3] = -999

# Latitude weighting
g1 = rep(cos(lat1*pi/180),each=nlon1)

# Read in data to get location of Indian land grid points
filedir = '/project/mjo/monsoon_heating/precip_india/'
year = 1901; ntall = 365; n = nlon1*nlat1*ntall
filename = paste(filedir,'ind',year,'_rfp25.grd',sep='')
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
varout[(nbuf+1):(nbuf+nlon1),(nbuf+1):(nbuf+nlat1),9:(8+4*9)] = box[,,4:(3+4*9)]
varout[(nbuf+1):(nbuf+nlon1),(nbuf+1):(nbuf+nlat1),45:48] = rep(box[,,2],4)
varout[(nbuf+1):(nbuf+nlon1),(nbuf+1):(nbuf+nlat1),49:52] = rep(box[,,3],4)
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
   filename = paste(filedir,'ind',year,'_rfp25.grd',sep='')
   temp = readBin(filename,'numeric',n=n,size='4',endian='little')
   dim(temp) = c(nlon1*nlat1,nt0)
   temp = temp[,iref + c(1:ntall)]
 
# Take areal averages
   for (z in 1:nreg) prec[,l,z] = apply(temp*mask[,z,],2,sum)
   print(year)
   }

# Write out areal-average IMD rainfall
# prec[,,1]       # Central India
# prec[,,2]       # All-India
# prec[,,3]       # All-India w/o Central India
# prec[,,4:39]    # Sub-regions within Central India
writeBin(as.vector(prec),'avgimd365.dat',size='4',endian='little')


