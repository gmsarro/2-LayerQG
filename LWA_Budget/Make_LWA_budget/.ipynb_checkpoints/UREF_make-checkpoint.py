import numpy as np
from netCDF4 import Dataset
import os,sys
#------------------------#
# calculate U_REF from Q_REF, all variables are nondimensionalized
datadir='/mnt/winds/data2/gmsarro/Rossbypalloza_project_22/LWA/' #input directory containing q_ref and zonal mean U
savedir=datadir # output directory
savedir='/mnt/winds/data2/gmsarro/Rossbypalloza_project_22/LWA/' #output directory

Llist= np.array([0.15,0.5], dtype=float) #np.arange(0.5,0.55,0.05)
Ulist= np.array([1.0, 1.1, 2, 2.1], dtype=float)#np.arange(2.4,2.5,0.1)
for u in range(len(Ulist)):
    for l in range(len(Llist)):
        sname_u='L_eddy_N128_%s_2.0_0.1_%s.uref1_2.nc'%(str(np.round(Llist[l],2)),str(np.round(Ulist[u],2)))  # output name U
        sname_t='L_eddy_N128_%s_2.0_0.1_%s.tref1_2.nc'%(str(np.round(Llist[l],2)),str(np.round(Ulist[u],2)))  # output name T

        #read zonal-mean QGPV for comparison
        readqm=Dataset(datadir+'L_eddy_N128_%s_2.0_0.1_%s.3d.nc'%(str(np.round(Llist[l],2)),str(np.round(Ulist[u],2)))) 
        qm=np.mean(readqm.variables['q1'][:,:,:].data,axis=2)

        # read Q_REF
        readq=Dataset(datadir+'L_eddy_N128_%s_2.0_0.1_%s.qref1_2.nc'%(str(np.round(Llist[l],2)),str(np.round(Ulist[u],2)))) 
        qref=readq.variables['qref1'][:,:].data
        tn,yn=np.shape(qref)

        # read zonal-mean U and T (used as boundary condition)
        readu=Dataset(datadir+'L_eddy_N128_%s_2.0_0.1_%s.nc'%(str(np.round(Llist[l],2)),str(np.round(Ulist[u],2))))
        um=readu.variables['zu1'][:,:].data
        umb=readu.variables['zu2'][:,:].data
        tm=readu.variables['ztau'][:,:].data
        ys=readu.variables['y'][:].data # nondimensional latitude

        # other dynamic parameters
        beta=0.2 # nondimensional beta
        Ld=1 # deformation radius

        # parameters for numerical solver
        maxerr=1E-6 # maximum tolerated error
        maxIT=10000 # maximum number of iteration
        relax=1.9 # relaxation parameter
        dy=ys[1]-ys[0] # grids should be evenly spaced
        AC=np.array([1/dy**2,-2/dy**2,1/dy**2]) # array for second-derivative

        #
        # calculate y gradient of Q_REF
        qref_y=np.zeros((tn,yn))
        qref_y[:,1:-1]=(qref[:,2:]-qref[:,:-2])/(2*dy)

        #========================================================================================#
        # solve for U_REF using SOR method
        uref=np.zeros((tn,yn))+um[:,:] # use um as the initial guess. This would make convergence faster if U_REF is close to zonal mean. Also um at the meridional boundaries are used.

        # below process solves 
        #$d^2 U_REF/ d y^2 - U_REF/Ld^2 - Beta + d Q_REF/ d y =0.


        # initial values
        nIT=0 #number of iteration
        err=1E+5 #initial error

        # start iteration
        while nIT < maxIT and err > maxerr:
            utemp=np.zeros((tn,yn)) # temporary array
            utemp[:,:]=uref[:,:] # upate with previous solution
            for y in range(1,yn-1):
                RS=(AC[0]*uref[:,y-1]+AC[1]*uref[:,y]+AC[2]*uref[:,y+1])-beta+qref_y[:,y]-uref[:,y]/Ld**2 + umb[:,y]/Ld**2
                uref[:,y]=uref[:,y]-relax*RS/(AC[1]-1/Ld**2)
            err=np.max(np.abs(uref[:,:]-utemp[:,:]))
            nIT+=1

        if nIT == maxIT:	print('Not fully converged')
        else:	print('converged at %i th iteration'%nIT)


        #========================================================================================#
        # numerically integrate to calculate T_REF
        tref=np.zeros((tn,yn)) # we will start with boundary condition being zero at the ends
        ushear=(uref[:,1:]+uref[:,:-1])*0.5 # yn-1 array defines shears between gridpoints
        for y in range(yn-1):
            tref[:,y+1]=tref[:,0]-np.sum(ushear[:,:y+1]*dy,axis=1)
        # tref boundary condition is updated to conserve temperature
        for t in range(tn):
            offset=np.mean(tm[t,:]-tref[t,:])
            tref[t,:]+=offset

        #========================================================================================#
        os.system('rm -f %s'%(savedir+sname_u))
        write=Dataset(savedir+sname_u,'w')
        write.createDimension('time',size = tn)
        write.createDimension('latitude',size = yn)
        var_u = write.createVariable('uref1','f4',dimensions=['time','latitude'])
        time = write.createVariable('time','f4',dimensions=['time'])
        latitude = write.createVariable('y','f4',dimensions=['latitude'])
        var_u[:,:]=uref[:,:]
        latitude=ys[:]
        time[:]=np.arange(tn)
        write.close()
        print('output saved!')
        #======================================================================================#
        #========================================================================================#
        os.system('rm -f %s'%(savedir+sname_t))
        write=Dataset(savedir+sname_t,'w')
        write.createDimension('time',size = tn)
        write.createDimension('latitude',size = yn)
        var_t = write.createVariable('tref1','f4',dimensions=['time','latitude'])
        time = write.createVariable('time','f4',dimensions=['time'])
        latitude = write.createVariable('y','f4',dimensions=['latitude'])
        var_t[:,:]=tref[:,:]
        latitude=ys[:]
        time[:]=np.arange(tn)
        write.close()
        print('output saved!')
        #======================================================================================#

