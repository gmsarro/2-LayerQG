#_____________________________________#
# Script for making dimensional QGPV  #
#_____________________________________#
import numpy as np
import sys
from netCDF4 import Dataset
import subprocess
import os
import xarray as xr

C= 2.0
Er= 0.1
U_1 = 0.5
L_arr= np.array([0.15,0.5])
U_arr = np.array([1,1.1,2,2.1])

for U_1 in U_arr:
    for L in L_arr:

        Lstr= str(L)
        Cstr= str(C)
        Estr= str(Er)
        Ustr= str(U_1)

# load namelist and 3d output #
#        dir_in='/project2/tas1/gmsarro/LWA/N128_%s_%s_%s_%s.3d.nc'%(Lstr,Cstr,Estr,Ustr)
#namelist='namelist.npz' # namelist output
#npz=np.load(dir_in+namelist) #reading namelist output
#for v in npz.files:	exec("%s=npz[v]"%v) # reading namelist output
#        file_name='.3d.nc' # 3d output

# Set dimensional variable
#U_DIM = 30 # dimensional maximum westerly
#Ld_DIM = 1E+6 # dimensional deformation radius
#T_DIM = Ld_DIM/(U_1*U_DIM) # dimensional time
#F0_DIM = 1E-4 # dimensional Coriolis factor

# Read Nondimensional Coordinates and QGPV
#read=Dataset(dir_in+file_name,'r')
#x_nondim=read.variables['x'][:].data
#y_nondim=read.variables['y'][:].data
#q1_nondim=read.variables['q1'][:,:,:].data
#read.close()

# Give dimensions
#x_dim = np.copy(x_nondim[:])#Ld_DIM*x_nondim[:] # dimensionalize longitude
#y_dim = Ld_DIM*y_nondim[:] # dimensionalize latitude
#t_dim = (T_DIM)*st*np.arange(np.shape(q1_nondim)[0]) # dimensionalize time
#t_dim /= (24*3600) # into days
#q1_dim = (U_1*U_DIM)/(Ld_DIM)*q1_nondim[:,:,:] # dimensionalize QGPV
        dir_in='/mnt/winds/data2/gmsarro/Rossbypalloza_project_22/LWA/L_eddy_N128_%s_%s_%s_%s.3d.nc'%(Lstr,Cstr,Estr,Ustr)
        file_name='.3d.nc'
# Write output
        dir_out = '/mnt/winds/data2/gmsarro/Rossbypalloza_project_22/LWA/L_eddy_N128_%s_%s_%s_%s.'%(Lstr,Cstr,Estr,Ustr) # output directory
        out_name = 'QGPV.nc' # output name
#write = Dataset(dir_out+out_name,'w')
#write.createDimension('time',size = len(t_dim))
#write.createDimension('latitude',size = len(y_dim))
#write.createDimension('longitude',size = len(x_dim))
#qgpv = write.createVariable('qgpv1','f4',dimensions=['time','latitude','longitude'])
#time = write.createVariable('time','f4',dimensions=['time'])
#latitude = write.createVariable('latitude','f4',dimensions=['latitude'])
#longitude = write.createVariable('longitude','f4',dimensions=['longitude'])
#qgpv.setncatts({'title': out_name})
#qgpv[:,:,:]=q1_dim[:,:,:]
#longitude[:]=x_dim[:]
#latitude[:]=y_dim[:]
#time[:]=t_dim
#write.close()

# Make link for fortran code
        dir_fortran ='/mnt/winds/data2/gmsarro/Rossbypalloza_project_22/LWA/run_LWA/'
    
##        os.system('rm %swac1.nc'%(dir_out))
##        os.system('rm %swaa1.nc'%(dir_out))
##        os.system('rm %sqref1.nc'%(dir_out))    
        
#        os.system('cp %s %sQGPV.nc'%(dir_in,dir_fortran)) #Move the file you are interested in
        
        ds = xr.open_dataset(dir_in)
        ds_cropped = ds 
        q1_mean = np.nanmean(ds_cropped['q1'], axis=2)

        # Find zonal mean location of maximum and minimum
        location_max = np.argmax(q1_mean, axis=1)
        location_min = np.argmin(q1_mean, axis=1)

        # Find zonal mean value of maximum and minimum
        value_max = np.max(q1_mean, axis=1)
        value_min = np.min(q1_mean, axis=1)

        # Replace values in the original q1 variable
        for i in range(ds.sizes['time']):
            ds_cropped['q1'][i, location_max[i]:, :] = value_max[i]
            ds_cropped['q1'][i, :location_min[i], :] = value_min[i]

        ds_cropped = xr.Dataset({'q1': ds_cropped['q1']})
        ds_cropped.to_netcdf(os.path.join(dir_fortran, out_name))
        ds.close()
        ds_cropped.close()
        os.system('chmod 744 QGPV.nc')
#        os.system('module load netcdf-fortran')
#        os.system('gfortran -c -g -fcheck=all $NC_INC rp2.f90 -o rp2.o')
#        os.system('gfortran -o rp2 rp2.o $NC_LIB')
#        os.system('./rp2')
        
            
        os.system('ifort -c rp2.f90 -I/opt/netcdf4/intel/include -L/opt/netcdf4/intel/lib -lnetcdf -lnetcdff -I/opt/hdf4/intel/include -L/opt/hdf4/intel/lib -lmhdf -ldf -ljpeg -lz -o rp2.o')
        os.system('ifort -O3 "./rp2.o" -mkl -I/opt/hdf4/intel/include -L/opt/hdf4/intel/lib -lmfhdf -ldf -ljpeg -lz -L/opt/netcdf4/intel/lib -lnetcdf -lnetcdff -o "./rp2"')
        os.system('rm *.o')
        os.system('./rp2')
        

#        os.system('gfortran -c -g -fcheck=all $NC_INC rp4.f90 -o rp4.o')
#        os.system('gfortran -o rp4 rp4.o $NC_LIB')
 #       os.system('./rp4')
        
        os.system('ifort -c rp4.f90 -I/opt/netcdf4/intel/include -L/opt/netcdf4/intel/lib -lnetcdf -lnetcdff -I/opt/hdf4/intel/include -L/opt/hdf4/intel/lib -lmhdf -ldf -ljpeg -lz -o rp4.o')
        os.system('ifort -O3 "./rp4.o" -mkl -I/opt/hdf4/intel/include -L/opt/hdf4/intel/lib -lmfhdf -ldf -ljpeg -lz -L/opt/netcdf4/intel/lib -lnetcdf -lnetcdff -o "./rp4"')
        os.system('rm *.o')
        os.system('./rp4')


#_____________________________________#
#remove temporary codes and move outputs
        os.system('rm -f rp2.o')
        os.system('rm -f rp4.o')
        os.system('rm -f rp2')
        os.system('rm -f rp4')
        os.system('rm -f QGPV.nc')
#
        os.system('mv %swac1.nc %swac1_2.nc'%(dir_fortran,dir_out))
        os.system('mv %swaa1.nc %swaa1_2.nc'%(dir_fortran,dir_out))
        os.system('mv %sqref1.nc %sqref1_2.nc'%(dir_fortran,dir_out))


        os.system('rm -f wac1.nc')
        os.system('rm -f waa1.nc')
        os.system('rm -f qref1.nc')

#_____________________________________#

