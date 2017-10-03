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
istep='000001'
inFile="tests_single_new.h5"
rhoH=1.0833
rhoL=1.0
#input done

mylist = [parentDir,'/',inFile]
delimiter = ''
filepath = delimiter.join(mylist)
#nz enlarged only 
variable = ['PVx','PVy','PVz','PPress', 'Prho']
h5file = h5py.File(filepath,'r')
#read dataset dimensions
mylist = ['Fields/','Prho','/',istep]
filepath = delimiter.join(mylist)
databk = h5file.get(filepath)
m1 = np.array(databk)
nz=m1.shape[0]
ny=m1.shape[1]
nx=m1.shape[2]
dz=3.2/nz
specout = 1000
h = np.zeros(len(steps))
seq = 0
step = []
for i in range(721):
    step.append(str((i+1)*specout).zfill(6))

#calculate mixing layer width
for istep in steps:
	delimiter = ''
	mylist = ['Fields/',variable[4],'/',istep]
	filepath = delimiter.join(mylist)
	databk = h5file.get(filepath)
	rho = np.array(databk)
	x=np.zeros(nz, ny, nx)
	xMean=np.zeros(nz)
	x=(rho-rhoL)/(rhoH-rhoL)
	xMean=x.reshape(nz, ny*nx).mean(axis=1)
	for i in range(nz):
		h[seq] = h[seq] + 2*min(xMean(i), 1-xMean(i)) 
	seq += 1


all_data = np.column_stack(np.asarray(steps),h)
np.savetxt('savedMixAndEnstro', all_data,delimiter='\t',fmt='%s')

h5file.close()

#f = open('output.d','w')
#for zz_ref in range(nz):
# f.write("%4s\t%10s\n" % (zz_ref, np.mean(m1[zz_ref,:,:])))
#f.close()
    
    
    
