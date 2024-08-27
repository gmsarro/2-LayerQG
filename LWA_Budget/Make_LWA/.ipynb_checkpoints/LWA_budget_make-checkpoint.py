import numpy as np
from netCDF4 import Dataset
from lwabudget import *
import sys
#================================#
loaddir='/project2/tas1/gmsarro/LWA/'
savedir='/project2/tas1/gmsarro/LWA/'
Llist=np.arange(0.0,0.35,0.05)
Ulist= np.array([1.0,1.2,1.8,2.4], dtype=float) #np.arange(2.3,2.5,0.1)
for u in range(len(Ulist)):
    for l in range(len(Llist)):
        sname='LH1_%s_2.0_0.1_%s.nc'%(str(np.round(Llist[l],2)),str(np.round(Ulist[u],2)))  
        #(A) VARIABLES================================#
        #1 load full variables
        read=Dataset(loaddir+'N128_%s_2.0_0.1_%s.3d.nc'%(str(np.round(Llist[l],2)),str(np.round(Ulist[u],2)))) 
        qdat=read.variables['q1'][:int(40000/8),:,:].data
        vdat=read.variables['v1'][:int(40000/8),:,:].data
        udat=read.variables['u1'][:int(40000/8),:,:].data
        tdat=read.variables['tau'][:int(40000/8),:,:].data
        xs=read.variables['x'][:].data
        ys=read.variables['y'][:].data
        read.close()
        #2 load reference state
        read=Dataset(loaddir+'N128_%s_2.0_0.1_%s.qref1_2.nc'%(str(np.round(Llist[l],2)),str(np.round(Ulist[u],2)))  )
        Qref=read.variables['qref1'][:int(40000/8),:].data
        read.close()
        read=Dataset(loaddir+'N128_%s_2.0_0.1_%s.uref1_2.nc'%(str(np.round(Llist[l],2)),str(np.round(Ulist[u],2)))  )
        Uref=read.variables['uref1'][:int(40000/8),:].data
        read.close()
        read=Dataset(loaddir+'N128_%s_2.0_0.1_%s.tref1_2.nc'%(str(np.round(Llist[l],2)),str(np.round(Ulist[u],2)))  )
        Tref=read.variables['tref1'][:int(40000/8),:].data
        read.close()
        #3 load LWA
        read=Dataset(loaddir+'N128_%s_2.0_0.1_%s.wac1_2.nc'%(str(np.round(Llist[l],2)),str(np.round(Ulist[u],2)))  )
        LWAC=read.variables['wac1'][:int(40000/8),:,:].data
        read.close()
        read=Dataset(loaddir+'N128_%s_2.0_0.1_%s.waa1_2.nc'%(str(np.round(Llist[l],2)),str(np.round(Ulist[u],2)))  )
        LWAA=read.variables['waa1'][:int(40000/8),:,:].data
        read.close()
        LWA=LWAA[:,:,:]+LWAC[:,:,:]
        print('variables loaded')
        #(B) EDDY VARIABLES================================#
        qe=qdat[:,:,:]-Qref[:,:,np.newaxis] # at y' =0
        ue=udat[:,:,:]-Uref[:,:,np.newaxis] # at y' =0
        ve=vdat[:,:,:] # at y'=0
        te=tdat[:,:,:]-Tref[:,:,np.newaxis] # at y' =0
        print('eddy variables calculated')
        #(C) BUTGET CALCULATION================================#
        # prepare coordinates
        times=np.linspace(0,10000,10000,endpoint=False)[:int(40000/8)]
        dt=times[1]-times[0]
        dx=xs[1]-xs[0]
        dy=ys[1]-ys[0]
        Ld=1.

        LWAtend=lwatend(LWA,dt)
        Urefadv=urefadv(LWA,Uref,dx,filt=False)
        ueqeadv=ueadv(qdat,Qref,udat,Uref,dx,dy,filt=False)
        EF_x=eddyflux_x(ue,ve,dx,filt=False)
        EF_y=eddyflux_y(ue,ve,dy,filt=False)
        EF_z=eddyflux_z(ve,te,Ld,filt=False)
        EF=eddyflux(ve,qe,filt=False)


        '''
        def run_mean(a,npt):
            # running mean along the last axis, a is 3d array
            a_pad=np.pad(a,((0,0),(0,0),(npt,npt)),mode='wrap')
            xn=np.shape(a)[2]
            out=np.zeros(np.shape(a))
            for x in range(xn):
                out[:,:,x]=np.mean(a_pad[:,:,x:x+2*npt+1],axis=2)
            return out
        N_PT=2 # half-width of running mean
        LWAtend=run_mean(LWAtend,N_PT)
        Urefadv=run_mean(Urefadv,N_PT)
        ueqeadv=run_mean(ueqeadv,N_PT)
        EF_x=run_mean(EF_x,N_PT)
        EF_y=run_mean(EF_y,N_PT)
        EF_z=run_mean(EF_z,N_PT)
        '''

        # 7 Residual
        RHS=Urefadv+ueqeadv+EF_x+EF_y+EF_z
        #RHS=Urefadv+ueqeadv+EF
        RES=LWAtend-RHS

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
        term1 = write.createVariable('lwatend','f4',dimensions=['time','latitude','longitude'])
        term2 = write.createVariable('urefadv','f4',dimensions=['time','latitude','longitude'])
        term3 = write.createVariable('ueqeadv','f4',dimensions=['time','latitude','longitude'])
        term4 = write.createVariable('ef_x','f4',dimensions=['time','latitude','longitude'])
        term5 = write.createVariable('ef_y','f4',dimensions=['time','latitude','longitude'])
        term6 = write.createVariable('ef_z','f4',dimensions=['time','latitude','longitude'])
        term7 = write.createVariable('res','f4',dimensions=['time','latitude','longitude'])
        #
        longitude[:]=xs[:]
        latitude[:]=ys[:]
        time[:]=times
        #
        term1[:,:,:]=LWAtend[:,:,:]
        term2[:,:,:]=Urefadv[:,:,:]
        term3[:,:,:]=ueqeadv[:,:,:]
        term4[:,:,:]=EF_x[:,:,:]
        term5[:,:,:]=EF_y[:,:,:]
        term6[:,:,:]=EF_z[:,:,:]
        term7[:,:,:]=RES[:,:,:]

        write.close()
        print("output saved ! done !" )
