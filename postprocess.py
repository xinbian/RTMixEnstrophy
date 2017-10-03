# -*- coding: utf-8 -*-
"""
To calculate 3D mixing layer width and enstrophy,
mixing layer defination:
https://search.proquest.com/docview/194682903?pq-origsite=gscholar
@author: Xin
"""

import h5py
import numpy as np
import os
import os.path

curfilePath = os.path.abspath(__file__)
curDir = os.path.abspath(os.path.join(curfilePath,os.pardir))
parentDir = os.path.abspath(os.path.join(curDir,os.pardir)) 

#not used 
#rhoL=1
#rhoH=1.0833

#specify inout parameters here
g=1.0
inFile="tests_single_new.h5"
rhoH=1.0
rhoL=0.1111
Lz=3.2
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

dz=Lz/nz
specout = 1000
seq = 0
step = []
for i in range(141):
    step.append(str((i+1)*specout).zfill(6))
#initialize time dependent mixing layer width, KE, PE, enstropy
h = np.zeros(len(step))
ke = np.zeros(len(step))
pe = np.zeros(len(step))
vorx = np.zeros((nz, ny, nx))
vory = np.zeros((nz, ny, nx))
vorz = np.zeros((nz, ny, nx))
enstropy = np.zeros(len(step))  
#calculate mixing layer width
for istep in step:
	delimiter = ''
	mylist = ['Fields/',variable[4],'/',istep]
	filepath = delimiter.join(mylist)
	databk = h5file.get(filepath)
	rho = np.array(databk)
	x=np.zeros((nz, ny, nx))
	xMean=np.zeros(nz)
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
	#enstrophy
	if nx == 1:
		vorx += np.gradient(vz, dz, axis=1) - np.gradient(vy, dz, axis=0)
	else:
		vorx += np.gradient(vz, dz, axis=1) - np.gradient(vy, dz, axis=0)
		vory += np.gradient(vx, dz, axis=0) - np.gradient(vz, dz, axis=2)
		vorz += np.gradient(vy, dz, axis=2) - np.gradient(vx, dz, axis=1)
	enstropy[seq] += (vorx**2+vory**2+vorz**2).mean() 
	seq += 1
	
	#test omega
	mylist = ['Fields/','PomegaX','/',istep]
	filepath = delimiter.join(mylist)
	h5file.create_dataset(filepath,data=vorx)


#nomalize pe
pe = pe/(nx*ny*nz)
#output
all_data = np.column_stack((np.asarray(step),h, ke, pe, ke + pe, enstropy))
np.savetxt('savedMixAndEnstro', all_data,delimiter='\t',fmt='%s')

h5file.close()

#f = open('output.d','w')
#for zz_ref in range(nz):
# f.write("%4s\t%10s\n" % (zz_ref, np.mean(m1[zz_ref,:,:])))
#f.close()
    
    
    
