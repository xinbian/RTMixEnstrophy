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
from CFDmodule import *

curfilePath = os.path.abspath(__file__)
curDir = os.path.abspath(os.path.join(curfilePath,os.pardir))
parentDir = os.path.abspath(os.path.join(curDir,os.pardir)) 

#specify inout parameters here
gamma=5.0/3.0
g=1.0
inFile="tests_single.h5"
Lz=3.2
waveLen = 0.4
CFDmethod = False
outPut = True
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
if nx == 1:
	dx=1.0
specout = 2000
skip = 10
seq = 0
step = []
for i in range(600000/2000):
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
sum_x = np.zeros(len(step))  

if CFDmethod:
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
	rhoL = 1.0
	rhoH = 1.0833
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
	ke[seq] = 0.5 * np.sum((vx**2 + vy**2 + vz**2) * rho)*dx*dy*dz
	#internal energy
	mylist = ['Fields/',variable[3],'/',istep]
	filepath = delimiter.join(mylist)
	databk = h5file.get(filepath)
	press = np.array(databk)
	ie[seq] =np.sum((1/(gamma-1))*press)*dx*dy*dz
	#enstrophy
	if CFDmethod:
		#CFD method to calc enstrophy
	        vx = vx.transpose()
	        vy = vy.transpose()
	        vz = vz.transpose()
	        
	        dyVx = calc_y_deri(vx, CFD_y)
	        dzVx = calc_z_deri(vx, CFD_z)
	        dxVy = calc_x_deri(vy, CFD_x)
	        dzVy = calc_z_deri(vy, CFD_z)
	        dxVz = calc_x_deri(vz, CFD_x)
	        dyVz = calc_y_deri(vz, CFD_y)
	        Wx = dyVz - dzVy
	        Wy = dzVx - dxVz
	        Wz = dxVy - dyVx
	        EnstrophyCFD = 0.5*(np.multiply(Wx,Wx)+np.multiply(Wy,Wy)+np.multiply(Wz,Wz))
	        sum_x[seq] = np.sum(EnstrophyCFD)*dx*dy*dz
	else:
		if nx == 1:
			vorx = np.gradient(vz, dz, axis=1) - np.gradient(vy, dz, axis=0)
		else:
			vorx = np.gradient(vz, dz, axis=1) - np.gradient(vy, dz, axis=0)
			vory = np.gradient(vx, dz, axis=0) - np.gradient(vz, dz, axis=2)
			vorz = np.gradient(vy, dz, axis=2) - np.gradient(vx, dz, axis=1)
		enstropy[seq] = 0.5*np.sum(vorx**2+vory**2+vorz**2)*dx*dy*dz 
	
	seq += 1

#nomalize pe/calcuate losed pe to this time
pe = pe*dx*dy*dz
#normalize h
h=h*dz/waveLen
#output
all_data = np.column_stack((np.asarray(step),h, ke, pe[0]-pe, ie-ie[0], enstropy, sum_x))
np.savetxt('savedMixAndEnstro', all_data,delimiter='\t',fmt='%s')

h5file.close()

if outPut:
	plt.plot(np.asarray(step),h)
	plt.title('mixing layer width vs time step')
	plt.savefig('mix.eps', format='eps', dpi=1000)
	plt.show()
	plt.plot(np.asarray(step),ke , label='KE')
	plt.plot(np.asarray(step),pe[0]-pe, label='released PE')
	plt.plot(np.asarray(step),ie-ie[0], label='increased IE')
	plt.plot(np.asarray(step),ke+pe,label='KE+PE')
	plt.plot(np.asarray(step),pe[0]-pe-(ie-ie[0]), label='released PE -increased IE')
	plt.title('energy vs time step')
	pylab.legend(loc='best')
	plt.savefig('energy.eps', format='eps', dpi=1000)
	plt.show()
	plt.plot(np.asarray(step), enstropy)
	plt.title('enstrophy vs time step')
	plt.savefig('enstropy.eps', format='eps', dpi=1000)
	plt.show()

#f = open('output.d','w')
#for zz_ref in range(nz):
# f.write("%4s\t%10s\n" % (zz_ref, np.mean(m1[zz_ref,:,:])))
#f.close()
    
    
    
