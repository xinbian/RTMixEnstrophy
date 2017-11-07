import numpy as np
import h5py

def Create_matrix_fd(num):
    v1 = np.concatenate(([5.0, 2.0/11.0], np.ones(num-4)/3.0, [2.0/11.0]))
    B = np.diagflat(v1, 1) + np.diagflat(np.ones(num))
    v1 = np.concatenate(([2.0/11.0], np.ones(num-4)/3.0, [2.0/11.0, 5.0]))
    B = B + np.diagflat(v1, -1)
    A = np.zeros((num, num))
    A[0,:] = np.concatenate(([-197.0/60,-5.0/12,5.0,-5.0/3,5.0/12,-1.0/20],np.zeros(num-6)))
    A[1,:] = np.concatenate(([-20.0/33,-35.0/132,34.0/33,-7.0/33,2.0/33,-1.0/132],np.zeros(num-6)))
    A[num-2,:] = np.concatenate((np.zeros(num-6),[1.0/132,-2.0/33,7.0/33,-34.0/33,35.0/132,20.0/33]))
    A[num-1,:] = np.concatenate((np.zeros(num-6),[1.0/20,-5.0/12,5.0/3,-5.0,5.0/12,197.0/60]))
    for ii in range(num-4):
        A[ii+2,:]=np.concatenate((np.zeros(ii),[-1.0/36,-7.0/9,0.0,7.0/9,1.0/36],np.zeros(num-ii-5)))
    return A, B

def Create_matrix_fd2(num):
    v1 = np.concatenate(([-1.0], np.zeros(num-2), [1] ))
    A = np.diagflat(v1, 0)
    v1 = np.concatenate(([1], 0.5* np.ones(num-2) ))
    A += np.diagflat(v1,1)
    v1 = np.concatenate(( -0.5* np.ones(num-2), [-1] ))
    A += np.diagflat(v1,-1)
    return A

def Create_matrix_filtering(num):
    B = np.diagflat(0.4*np.ones(num-1),1) + np.diagflat(0.4*np.ones(num-1),-1) + np.diagflat(np.ones(num),0)
    A = np.zeros((num,num))
    A[0,:] = np.concatenate(([317.0/320,73.0/160,-9.0/64,3.0/16,-9.0/64,9.0/160,-3.0/320] , np.zeros(num-7)))
    A[1,:]=np.concatenate(([129.0/320,157.0/160,143.0/320,-1.0/16,3.0/64,-3.0/160,1.0/320] , np.zeros(num-7)))
    A[2,:]=np.concatenate(([-1.0/320,67.0/160,61.0/64,37.0/80,-3.0/64,3.0/160,-1.0/320] , np.zeros(num-7)))
    A[num-3,:]=np.concatenate((np.zeros(num-7),[-1.0/320,3.0/160,-3.0/64,37.0/80,61.0/64,67.0/160,-1.0/320]))
    A[num-2,:]=np.concatenate((np.zeros(num-7),[1.0/320,-3.0/160,3.0/64,-1.0/16,143.0/320,157.0/160,129.0/320]))
    A[num-1,:]=np.concatenate((np.zeros(num-7),[-3.0/320,9.0/160,-9.0/64,3.0/16,-9.0/64,73.0/160,317.0/320]))
    for ii in range(num-6):
        A[ii+3,:] = np.concatenate((np.zeros(ii),[1.0/320,-3.0/160,143.0/320,15.0/16,143.0/320,-3.0/160,1.0/320],np.zeros(num-7-ii)))
    return A, B

def calc_x_deri(rho, CFD_x):
    nx, ny, nz = rho.shape
    result = np.zeros((nx, ny, nz))
    for j in range(ny):
        for k in range(nz):
            result[:, j, k] = np.dot(CFD_x, rho[:, j, k])
    return result

def calc_y_deri(rho, CFD_y):
    nx, ny, nz = rho.shape
    result = np.zeros((nx, ny, nz))
    for i in range(nx):
        for k in range(nz):
            result[i, :, k] = np.dot(CFD_y, rho[i, :, k])
    return result

def calc_z_deri(rho, CFD_z):
    nx, ny, nz = rho.shape
    result = np.zeros((nx, ny, nz))
    for i in range(nx):
        for j in range(ny):
            result[i, j, :] = np.dot(CFD_z, rho[i, j, :])
    return result

if __name__ == "__main__":

    file_name = "tests_multi_ablative_lowerq.h5"
    H5File = h5py.File(file_name,'r')
    
    steps = [12000, 13500, 15000, 17500, 19000, 20500]
    
    nx = 256
    ny = 256
    nz = 512
    Lx = 3.2
    Ly = 3.2
    Lz = 6.4
    dx = Lx/nx
    dy = Ly/ny
    dz = Lz/nz

    
    #A, B = Create_matrix_fd(nx)
    #A = np.array(A)
    #B = np.array(B)
    #CFD_x = np.linalg.solve(B, A)
    #CFD_x = CFD_x / dx

    #CFD_y = CFD_x

    #A, B = Create_matrix_fd(nz)
    #A = np.array(A)
    #B = np.array(B)
    #CFD_z = np.linalg.solve(B, A)
    #CFD_z = CFD_z / dz

    CFD_x = Create_matrix_fd2(nx) / dx
    CFD_y = CFD_x
    CFD_z = Create_matrix_fd2(nz) / dz


    #A, B = Create_matrix_filtering(nz)
    #A = np.array(A)
    #B = np.array(B)
    #CFD_filter_z = np.linalg.solve(B, A)
    #CFD_filter_z = CFD_filter_z

    
    #for x in steps:
    for x in range(1000, 27000, 1000):
        Vx = H5File['/Fields/PVx/' + str(x).zfill(6)]
        Vy = H5File['/Fields/PVy/' + str(x).zfill(6)]
        Vz = H5File['/Fields/PVz/' + str(x).zfill(6)]
        Vx = np.array(Vx)
        Vy = np.array(Vy)
        Vz = np.array(Vz)
        Vx = Vx.transpose()
        Vy = Vy.transpose()
        Vz = Vz.transpose()
        
        dyVx = calc_y_deri(Vx, CFD_y)
        dzVx = calc_z_deri(Vx, CFD_z)
        dxVy = calc_x_deri(Vy, CFD_x)
        dzVy = calc_z_deri(Vy, CFD_z)
        dxVz = calc_x_deri(Vz, CFD_x)
        dyVz = calc_y_deri(Vz, CFD_y)

        Wx = dyVz - dzVy
        Wy = dzVx - dxVz
        Wz = dxVy - dyVx

        Enstrophy = 0.5*(np.multiply(Wx,Wx)+np.multiply(Wy,Wy)+np.multiply(Wz,Wz))
        sum_x = np.sum(Enstrophy)*dx*dy*dz
        print("At step ", x, ", enstrophy is",sum_x)
        
        #hahax = np.gradient(Vx, dx ,dy, dz)
        #hahay = np.gradient(Vy, dx ,dy, dz)
        #hahaz = np.gradient(Vz, dx ,dy, dz)
        #hahax = np.gradient(Vx, dx ,dy, dz)
        #hahay = np.gradient(Vy, dx ,dy, dz)
        #hahaz = np.gradient(Vz, dx ,dy, dz)
        #hahax = np.array(hahax)
        #hahay = np.array(hahay)
        #hahaz = np.array(hahaz)
        #Wx = hahaz[1,:,:,:] - hahay[2,:,:,:]
        #Wy = hahax[2,:,:,:] - hahaz[0,:,:,:]
        #Wz = hahay[0,:,:,:] - hahax[1,:,:,:]

        #Enstrophy = 0.5*(np.multiply(Wx,Wx)+np.multiply(Wy,Wy)+np.multiply(Wz,Wz))
        
        #sum_x = np.sum(Enstrophy)*dx*dy*dz
        #print("At step ", x, ", enstrophy is",sum_x)
        
