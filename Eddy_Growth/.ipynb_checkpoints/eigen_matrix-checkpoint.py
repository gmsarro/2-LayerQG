#Computes the fastest eddy growth rates in a 2-layer QG model.
#This is a modified version of the Phillips model.
#Written by Giorgio M. Sarro


#This is an eigensolver of a matrix built by inputting wave perturbatios in the form A2(y)e^(i(kx-ct)) and A1(y)e^(i(kx-ct)) into the lower and upper layer PV equations.
#The wave amplitudes (A1 and A2) are unknown, so they are optimized to reach the largest possible growth rate

import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import netCDF4 as nc
from scipy.optimize import root
from scipy.optimize import fsolve

L = 1  # Rossby radius of deformation
#==========================================================#
# Load your 3D model data: upper layer PV (q1), lower layer PV (q2), upper layer wind (u1) and lower layer wind (u2)
#Take the time and zonal mean
data_path = "/mnt/winds/data2/gmsarro/Rossbypalloza_project_22/LWA/N128_0.0_2.0_0.1_1.0.3d.nc"
with nc.Dataset(data_path, 'r') as nc_file:
    q1_mean = np.mean(np.mean(nc_file.variables['q1'][:], axis=2), axis=0)
    q2_mean = np.mean(np.mean(nc_file.variables['q2'][:], axis=2), axis=0)
    u1_mean_1 = np.mean(np.mean(nc_file.variables['u1'][:], axis=2), axis=0)
    u2_mean_1 = np.mean(np.mean(nc_file.variables['u2'][:], axis=2), axis=0)
    y = nc_file.variables['y'][:]
length_of_y = int(np.size(y))
#==========================================================#
#Calculate the meridional gradient of PV in both layers
dq1_dy = np.gradient(q1_mean, y, axis=0)
dq2_dy = np.gradient(q2_mean, y, axis=0)

#==========================================================#
#These are our zonal wavenumbers (from 0 to 5, it can be modified according to your domain size)'
k_value = np.linspace(0, 5, 1000)
#k=1 Loop through k
largest_imaginary_eigenvalues = []  # To store the largest imaginary eigenvalues for each k
optimal_d_A1_dy_values = []  # To store the optimal d_A1_dy arrays for each k
optimal_d_A2_dy_values = []  
#==========================================================#

#Let's build our matrix.
for k in k_value:
    # Define a function that computes the objective to maximize (the negative imaginary part of the eigenvalue)
    def objective(d_A):
        d_A1_dy = d_A[:length_of_y]
        d_A2_dy = d_A[length_of_y:]    
        main_matrix = np.zeros((length_of_y*2, length_of_y*2))    
        for column in range(length_of_y*2):
            for row in range(length_of_y*2):
                if column < length_of_y and row < length_of_y:
                    if column == row:
                        main_matrix[column, row] = u1_mean_1[column] * k * (-k ** 2 + d_A1_dy[column] - 1 / L ** 2) + k * dq1_dy[column]
                if column >= length_of_y and row < length_of_y:
                    if column - length_of_y == row:
                        main_matrix[column, row] = u1_mean_1[column - length_of_y] * k / L ** 2
                if column < length_of_y and row >= length_of_y:
                    if column == row - length_of_y:
                        main_matrix[column, row] = u2_mean_1[column] * k / L ** 2
                if column >= length_of_y and row >= length_of_y:
                    if column - length_of_y == row - length_of_y:
                        main_matrix[column, row] = u2_mean_1[column - length_of_y] * k * (-k ** 2 + d_A2_dy[column - length_of_y] - 1 / L ** 2) + k * dq2_dy[column - length_of_y]
#==========================================================#
        #Let's solve the eigenvalues
        
        eigenvalues, _ = np.linalg.eig(main_matrix)
        index_of_max_imaginary = np.argmax(np.imag(eigenvalues))
        largest_imaginary_eigenvalue = eigenvalues[index_of_max_imaginary].imag

        # We want to maximize the negative imaginary part, so we return the negative value
        return -largest_imaginary_eigenvalue

    # Define the constraints for integration
    # Initial guess for d_A1_dy and d_A2_dy
    initial_guess = np.zeros(length_of_y*2)
    initial_guess[:length_of_y] = (np.random.rand(length_of_y)-0.5)*20
    initial_guess[length_of_y:] = (np.random.rand(length_of_y)-0.5)*10


    # Run the optimization
    result = minimize(objective, initial_guess, method="SLSQP")

    # Extract the optimized values
    optimal_d_A1_dy = result.x[:length_of_y]
    optimal_d_A2_dy = result.x[length_of_y:]
    largest_imaginary_eigenvalue = result.fun  # Negative of the largest imaginary eigenvalue

    # Save the results
    largest_imaginary_eigenvalues.append(-largest_imaginary_eigenvalue)
    optimal_d_A1_dy_values.append(optimal_d_A1_dy)
    optimal_d_A2_dy_values.append(optimal_d_A2_dy)

#==========================================================#
#Save the growth rates as a function of wavenumber and the meridional derivative of the wave amplitude as a netCDF file

output_file = "eigenvalue_results.nc"
with nc.Dataset(output_file, 'w', format='NETCDF4_CLASSIC') as nc_file:
    # Create dimensions for wavenumber
    nc_file.createDimension('k_dim', len(k_value))
    # Create dimension for 'y' 
    nc_file.createDimension('y_dim', length_of_y)
    # Create variables to store the data
    nc_file.createVariable('k', 'f8', ('k_dim',))
    nc_file.createVariable('largest_imaginary_eigenvalues', 'f8', ('k_dim',))
    nc_file.createVariable('optimal_d_A1_dy_values', 'f8', ('k_dim', 'y_dim'))
    nc_file.createVariable('optimal_d_A2_dy_values', 'f8', ('k_dim', 'y_dim'))

    # Save data into variables
    nc_file['k'][:] = k_value
    nc_file['largest_imaginary_eigenvalues'][:] = largest_imaginary_eigenvalues
    nc_file['optimal_d_A1_dy_values'][:] = optimal_d_A1_dy_values
    nc_file['optimal_d_A2_dy_values'][:] = optimal_d_A2_dy_values


