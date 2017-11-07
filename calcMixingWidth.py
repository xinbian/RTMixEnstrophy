import numpy as np
import h5py

file_name = 'tests_single_new.h5'
H5File = h5py.File(file_name, 'r')
#rho_l = 1.0/3.0
#rho_h = 1.0
skip = 10

L_z = 0.8
Nz = 512
dz = L_z/Nz

for x in range(500, 19920, 500):
    rho = H5File['/Fields/Prho/' + str(x).zfill(6)]
    rho = np.array(rho)
    rho = np.transpose(rho)
    rho_l = np.average(rho[:,:,skip], axis=(0,1))
    rho_h = np.average(rho[:,:,Nz-skip], axis=(0,1))
    rho = (rho - rho_l) / (rho_h - rho_l)
    rho_z = np.average(rho, axis=(0,1) )

    for i in range(np.size(rho_z)):
        if rho_z[i] <= 0.5:
            rho_z[i] = 2 * rho_z[i]
        else:
            rho_z[i] = 2 * (1 - rho_z[i])

    width = np.sum(rho_z) * dz
    print('At step ', x ,', width is ',width)
