#Computes all terms of the latent heating contribution of the LWA budget
#Written by Joonsuk M. Kang
#Modified by Giorgio M. Sarro

import numpy as np
from netCDF4 import Dataset
from lwabudget import *
import sys
import os
#================================#
#Load directory where the data is stored
loaddir='/project2/tas1/gmsarro/LWA/'
#Specify where to save the data
savedir='/project2/tas1/gmsarro/LWA/'

#Specify the value of latent heating (L) and the relaxation zonal wind (U)
Llist=np.arange(0.0)
Ulist= np.array([1.0], dtype=float)
#Specify the length of your dataset: ##Attention: If the datset is too long, the computation might crash. 
max_lenghth= 10000

sname='LP_%s_2.0_0.1_%s.nc'%(str(np.round(Llist[0],2)),str(np.round(Ulist[0],2)))  
#(A) VARIABLES================================#
#1 load full variables, including precipitation (P)
read=Dataset(loaddir+'N128_%s_2.0_0.1_%s.3d.nc'%(str(np.round(Llist[0],2)),str(np.round(Ulist[0],2)))) 
qdat=read.variables['q1'][:,:,:].data
pdat=read.variables['P'][:,:,:].data
xs=read.variables['x'][:].data
ys=read.variables['y'][:].data
read.close()
#2 load reference state
read=Dataset(loaddir+'N128_%s_2.0_0.1_%s.qref1_2.nc'%(str(np.round(Llist[0],2)),str(np.round(Ulist[0],2)))  )
Qref=read.variables['qref1'][:,:].data
read.close()
print('variables loaded')
#(B) EDDY VARIABLES================================#
L=float(Llist[l])
#(C) BUTGET CALCULATION================================#
# prepare coordinates
times=np.linspace(0,max_lenghth,max_lenghth,endpoint=False)[:]
dt=times[1]-times[0]
dx=xs[1]-xs[0]
dy=ys[1]-ys[0]

LP=LH(pdat,qdat,Qref,L,dx,dy,filt=False)
print('budget calculated')
#Save the result       
#==========================================================#
        
os.system('rm -f %s%s'%(savedir,sname))
write = Dataset(savedir+sname,'w')
write.createDimension('time',size = len(times))
write.createDimension('latitude',size = len(ys))
write.createDimension('longitude',size = len(xs))
#
time = write.createVariable('time','f4',dimensions=['time'])
latitude = write.createVariable('latitude','f4',dimensions=['latitude'])
longitude = write.createVariable('longitude','f4',dimensions=['longitude'])
#
term1 = write.createVariable('LH','f4',dimensions=['time','latitude','longitude'])
#
longitude[:]=xs[:]
latitude[:]=ys[:]
time[:]=times
#
term1[:,:,:]=LP[:,:,:]

write.close()
print("output saved ! done !")
