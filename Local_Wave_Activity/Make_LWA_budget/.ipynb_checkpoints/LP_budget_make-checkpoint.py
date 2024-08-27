import numpy as np
from netCDF4 import Dataset
from lwabudget import *
import sys
#================================#
#runname,t1,t2=sys.argv[1],int(sys.argv[2]),int(sys.argv[3])
loaddir='/mnt/winds/data2/gmsarro/Rossbypalloza_project_22/LWA/'
savedir='/mnt/winds/data2/gmsarro/Rossbypalloza_project_22/LWA/'

Llist=np.array([0.15,0.5], dtype=float)
Ulist= np.array([1.0,1.1,2.0,2.1], dtype=float)

#Llist=np.arange(0,0.55,0.05)
#Ulist=np.arange(1,2.5,0.1)
for u in range(len(Ulist)):
    for l in range(len(Llist)):

        sname='L_eddy_LP_%s_2.0_0.1_%s.nc'%(str(np.round(Llist[l],2)),str(np.round(Ulist[u],2)))  
        #(A) VARIABLES================================#
        #1 load full variables
        read=Dataset(loaddir+'L_eddy_N128_%s_2.0_0.1_%s.3d.nc'%(str(np.round(Llist[l],2)),str(np.round(Ulist[u],2)))) 
        qdat=read.variables['q1'][:int(40000/8),:,:].data
        pdat=read.variables['Removed_P'][:int(40000/8),:,:].data
        xs=read.variables['x'][:].data
        ys=read.variables['y'][:].data
        read.close()
        #2 load reference state
        read=Dataset(loaddir+'L_eddy_N128_%s_2.0_0.1_%s.qref1_2.nc'%(str(np.round(Llist[l],2)),str(np.round(Ulist[u],2)))  )
        Qref=read.variables['qref1'][:int(40000/8),:].data
        read.close()
        print('variables loaded')
        #(B) EDDY VARIABLES================================#
        L=float(Llist[l])
        #(C) BUTGET CALCULATION================================#
        # prepare coordinates
        times=np.linspace(0,10000,10000,endpoint=False)[:int(40000/8)]
        dt=times[1]-times[0]
        dx=xs[1]-xs[0]
        dy=ys[1]-ys[0]
        Ld=1.

        LP=LH(pdat,qdat,Qref,L,dx,dy,filt=False)



        print('budget calculated')
        #==========================================================#
        import os
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
