# -*- coding: utf-8 -*-
"""
To calculate 3D mixing layer width and enstrophy,
mixing layer defination:
https://search.proquest.com/docview/194682903?pq-origsite=gscholar
@author: Xin
"""
import pylab
import matplotlib.pyplot as plt
import h5py
import numpy as np
import os
import os.path

curfilePath = os.path.abspath(__file__)
curDir = os.path.abspath(os.path.join(curfilePath,os.pardir))
parentDir = os.path.abspath(os.path.join(curDir,os.pardir)) 


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



#not used 
#rhoL=1
#rhoH=1.0833

#specify inout parameters here
gamma=5.0/3.0
g=1.0
inFile="tests_single_new.h5"
Lz=0.8
#input done

mylist = [parentDir,'/',inFile]
delimiter = ''
filepath = delimiter.join(mylist)
#nz enlarged only 
variable = ['PVx','PVy','PVz','PPress', 'Prho']
h5file = h5py.File(filepath,'r+')
#read dataset dimensions
mylist = ['Fields/','Prho','/','002000']
filepath = delimiter.join(mylist)
databk = h5file.get(filepath)
m1 = np.array(databk)
nz=m1.shape[0]
ny=m1.shape[1]
nx=m1.shape[2]

dz=dy=dx=Lz/nz
specout = 500
skip = 10
seq = 0
step = []
for i in range(39):
    step.append(str((i+1)*specout).zfill(6))
#initialize time dependent mixing layer width, KE, PE, enstropy
h = np.zeros(len(step))
ie = np.zeros(len(step))
ke = np.zeros(len(step))
pe = np.zeros(len(step))
vorx = np.zeros((nz, ny, nx))
vory = np.zeros((nz, ny, nx))
vorz = np.zeros((nz, ny, nx))
enstropy = np.zeros(len(step))  

CFD_x = Create_matrix_fd2(nx) / dx
CFD_y = CFD_x
CFD_z = Create_matrix_fd2(nz) / dz


#calculate mixing layer width
for istep in step:
	delimiter = ''
	mylist = ['Fields/',variable[4],'/',istep]
	filepath = delimiter.join(mylist)
	databk = h5file.get(filepath)
	rho = np.array(databk)
	x=np.zeros((nz, ny, nx))
	xMean=np.zeros(nz)
	rhoL = rho[skip,:, :].mean()
	rhoH = rho[nz-skip,:, :].mean()
	x=(rho-rhoL)/(rhoH-rhoL)
	xMean=x.reshape(nz, ny*nx).mean(axis=1)
	for i in range(nz):
		h[seq] = h[seq] + 2*min(xMean[i], 1-xMean[i]) 
		#potential energy
		pe[seq] = pe[seq] + np.sum(rho[i, :, :] * g * i * dz)
 	#kinetic energy
	mylist = ['Fields/',variable[0],'/',istep]
	filepath = delimiter.join(mylist)
	databk = h5file.get(filepath)
	vx = np.array(databk)
	mylist = ['Fields/',variable[1],'/',istep]
	filepath = delimiter.join(mylist)
	databk = h5file.get(filepath)
	vy = np.array(databk)
	mylist = ['Fields/',variable[2],'/',istep]
	filepath = delimiter.join(mylist)
	databk = h5file.get(filepath)
	vz = np.array(databk)
	ke[seq] = (0.5 * (vx**2 + vy**2 + vz**2) * rho).mean()
	#internal energy
	mylist = ['Fields/',variable[3],'/',istep]
	filepath = delimiter.join(mylist)
	databk = h5file.get(filepath)
	press = np.array(databk)
	ie[seq] = ((1/(gamma-1))*press).mean() 
	#enstrophy
	if nx == 1:
		vorx = np.gradient(vz, dz, axis=1) - np.gradient(vy, dz, axis=0)
	else:
		vorx = np.gradient(vz, dz, axis=1) - np.gradient(vy, dz, axis=0)
		vory = np.gradient(vx, dz, axis=0) - np.gradient(vz, dz, axis=2)
		vorz = np.gradient(vy, dz, axis=2) - np.gradient(vx, dz, axis=1)
	enstropy[seq] += (vorx**2+vory**2+vorz**2).mean() 
	seq += 1

        Vx = H5File['/Fields/PVx/' + istep]
        Vy = H5File['/Fields/PVy/' + istep]
        Vz = H5File['/Fields/PVz/' + istep]
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
	
	#test omega
#	mylist = ['Fields/','PomegaX','/',istep]
#	filepath = delimiter.join(mylist)
#	h5file.create_dataset(filepath,data=vorx)


#nomalize pe/calcuate losed pe to this time
pe = pe/(nx*ny*nz)
#normalize h
h=h*dz
#output
all_data = np.column_stack((np.asarray(step),h, ke, pe[0]-pe, ie-ie[0], enstropy, sum_x)
np.savetxt('savedMixAndEnstro', all_data,delimiter='\t',fmt='%s')

h5file.close()


#plt.plot(np.asarray(step),h)
#plt.title('mixing layer width vs time step')
#plt.savefig('mix.eps', format='eps', dpi=1000)
#plt.show()
#plt.plot(np.asarray(step),ke , label='KE')
#plt.plot(np.asarray(step),pe[0]-pe, label='released PE')
#plt.plot(np.asarray(step),ie-ie[0], label='increased IE')
##plt.plot(np.asarray(step),ke+pe,label='KE+PE')
#plt.plot(np.asarray(step),pe[0]-pe-(ie-ie[0]), label='released PE -increased IE')
#plt.title('energy vs time step')
#pylab.legend(loc='best')
#plt.savefig('energy.eps', format='eps', dpi=1000)
#plt.show()
#plt.plot(np.asarray(step), enstropy)
#plt.title('enstrophy vs time step')
#plt.savefig('enstropy.eps', format='eps', dpi=1000)
#plt.show()
#
#f = open('output.d','w')
#for zz_ref in range(nz):
# f.write("%4s\t%10s\n" % (zz_ref, np.mean(m1[zz_ref,:,:])))
#f.close()
    
    
    
