# -*- coding: utf-8 -*-
"""
To calculate enstrophy behind bubble and spike,
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
from readStep import *


#specify inout parameters here
g=1.0
Lz=3.2
##maximum average height (grid point)
#average over how much percent of total height
winPercent = 0.025

rhoH = 1.0833
rhoL = 1.0
waveLen = 0.4
#input done

rd = rhoL/rhoH
k=2*np.pi/waveLen
rho_inc = (rhoH + rhoL)/2.0

#nz enlarged only 
variable = ['PVx','PVy','PVz','PPress', 'Prho']

winPoint =  int(winPercent*nz)

delimiter = ''
dz=dy=dx=Lz/nz
if nx == 1:
	dx = 1.0

FieldPoint = h5file.get('Fields').values()
totalsteps, specout = get_LatestTime(FieldPoint)

for i in range(totalsteps/specout):
    step.append(str((i+1)*specout).zfill(6))


vorx = np.zeros((nz, ny, nx))
vory = np.zeros((nz, ny, nx))
vorz = np.zeros((nz, ny, nx))


bub_loc_all = np.zeros(len(step))
bub_loc_all_ori = np.zeros(len(step))
sp_loc_all = np.zeros(len(step))
bub_velo_all = np.zeros(len(step))
bub_velo_all_aver = np.zeros(len(step))
sp_velo_all = np.zeros(len(step))
bub_velo_all_ori = np.zeros(len(step))
ensbub = np.zeros(len(step))
enspk = np.zeros(len(step))
ensbub2 = np.zeros(len(step))
enspk2 = np.zeros(len(step))


i=0


test = []
#calculate bubble and spike location
for istep in step:
    mylist = ['Fields/', variable[4], '/', istep]
    filepath = delimiter.join(mylist)
    databk = h5file.get(filepath)
    rho_data = np.array(databk)
    if nx == 1:
	    m1 = (rho_data[:, ny/2-1, 0] + rho_data[:, ny/2, 0] )/2
    else:
	    m1 = (rho_data[:, ny/2-1, nx/2-1] + rho_data[:, ny/2, nx/2] 
	  + rho_data[:, ny/2-1, nx/2] + rho_data[:, ny/2, nx/2-1])/4.0
    m2 = rho_data[:, 0, 0]
    m1_filter=m1.copy();    
    m2_filter=m2.copy();    
    for jstep in range(2,nz-3):
        m1_filter[jstep]=(m1[jstep-2]+m1[jstep-1]+m1[jstep]+m1[jstep+1]+m1[jstep+2])/5;
        m2_filter[jstep]=(m2[jstep-2]+m2[jstep-1]+m2[jstep]+m2[jstep+1]+m2[jstep+2])/5;

    m1_grad = high_order_gradient(m1_filter,dx,6)
    m2_grad = high_order_gradient(m2_filter,dx,6)

    sp_loc = np.argmax(m1_grad)
    bub_loc = np.argmax(m2_grad)


    sp_loc_all[i] = sp_loc
    bub_loc_all[i] = bub_loc

    print 'finish', 100*float(istep)/totalstep,'%'
    
    #obtain velocity 
    
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
    
    
    if nx == 1:
    	    m1 = (vz[:, ny/2-1, 0] + vz[:, ny/2, 0] )/2
    else:
    	    m1 = (vz[:, ny/2-1, nx/2-1] + vz[:, ny/2, nx/2] 
    	  + vz[:, ny/2-1, nx/2] + vz[:, ny/2, nx/2-1])/4.0

    bub_velo_all[i] = m1[bub_loc]
    sp_velo_all[i] = m1[sp_loc]
    
    #cacluate vorticity
    if nx == 1:
		vorx = np.gradient(vz, dz, axis=1) - np.gradient(vy, dz, axis=0)
    else:
		vorx = np.gradient(vz, dz, axis=1) - np.gradient(vy, dz, axis=0)
		vory = np.gradient(vx, dz, axis=0) - np.gradient(vz, dz, axis=2)
		vorz = np.gradient(vy, dz, axis=2) - np.gradient(vx, dz, axis=1)


    #set avearging area
    bubRegion = int(bub_loc_all[i] - winPoint)
    spkRegion = int(sp_loc_all[i] + winPoint)
    
    #calcuate vorticity behind bub 
    for j in range(ny/2):
    
        for k in reversed(range(nz)):
            
            if (rho_data[k, j, 0] > rho_inc) and (rho_data[k-1, j, 0] < rho_inc):
                bub_inc_loc = k
                break
  
        if bub_inc_loc > bubRegion:
            ensbub[i] = ensbub[i] + np.sum(vorx[bubRegion:bub_inc_loc, j, :]**2
                     + vory[bubRegion:bub_inc_loc, j, :]**2
                     + vorz[bubRegion:bub_inc_loc, j, :]**2)*dx*dy*dz 
                  
                  
                 
                  
                  
    #calcuate vorticity behind spike             
    for j in reversed(range(ny/2)):
    
        for k in range(nz):
            
            if (rho_data[k, j, 0] < rho_inc) and (rho_data[k+1, j, 0] > rho_inc):
                spk_inc_loc = k
                break
    
        if spk_inc_loc < spkRegion:
            enspk[i] = enspk[i] + np.sum(vorx[spk_inc_loc:spkRegion, j, :]**2
                     + vory[spk_inc_loc:spkRegion, j, :]**2
                     + vorz[spk_inc_loc:spkRegion, j, :]**2)*dx*dy*dz 
        test.append(spk_inc_loc)
                 
    ensbub2[i] = ensbub2[i] + np.sum(vorx[bubRegion:int(bub_loc_all[i]), :, :]**2
                     + vory[bubRegion:int(bub_loc_all[i]), :, :]**2
                     + vorz[bubRegion:int(bub_loc_all[i]), :, :]**2)*dx*dy*dz 
    enspk2[i] = enspk2[i] + np.sum(vorx[int(sp_loc_all[i]):spkRegion, :, :]**2
                     + vory[int(sp_loc_all[i]):spkRegion, :, :]**2
                     + vorz[int(sp_loc_all[i]):spkRegion, :, :]**2)*dx*dy*dz 
       


    i = i + 1

    
ensbub = 2 * ensbub
enspk = 2 * enspk

bubspeed = np.sqrt(1/(3*np.pi)+rd/(1-rd)*ensbub/(4*np.pi*k*g))
plt.plot(ensbub, label='bub curve')
plt.plot(ensbub2, label='bub square')
plt.plot(enspk, label='spk curve')
plt.plot(enspk2, label='spk square')
pylab.legend(loc='best')
pylab.savefig('ensaftertip')
plt.show()

plt.plot(bubspeed)


h5file.close()

np.savetxt('bubspeed', bubspeed, delimiter=',')   
np.savetxt('bub', ensbub, delimiter=',')     
np.savetxt('bubsq', ensbub2, delimiter=',')     
np.savetxt('spk', enspk, delimiter=',')     
np.savetxt('spksq', enspk2, delimiter=',')     
 
