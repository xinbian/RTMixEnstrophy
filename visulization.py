# -*- coding: utf-8 -*-
"""
To calculate 3D mixing layer width and enstrophy,
mixing layer defination:
https://search.proquest.com/docview/194682903?pq-origsite=gscholar
@author: Xin
"""
import matplotlib.pyplot as plt
import h5py
import numpy as np
import os
import os.path
import copy
from CFDmodule import *
from readStep import *



#specify inout parameters here
gamma=5.0/3.0
g=1.0
inFile="tests_single_new.h5"
Ly=0.4
Lz=3.2
CFDmethod = False
rho_l=0.4
rho_h=1.6
rho_inc = (rho_h + rho_l)/2.0
winPercent = 0.05
winPoint =  int(winPercent*nz)
#input done
variable = ['PVx','PVy','PVz','PPress', 'Prho']
dz=dy=dx=Lz/nz

if nx == 1:
    dx=1.0


horizon_lim = (0, Ly-dy)
vert_lim = (0, Lz-dz)
extent=horizon_lim+vert_lim


FieldPoint = h5file.get('Fields').values()
totalsteps, specout = get_LatestTime(FieldPoint)

step = []
#totalstep=4000
#specout=2000
for i in range(totalstep/specout):
    step.append(str((i+1)*specout).zfill(6))


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
    print istep
    mylist = ['Fields/','Prho','/',istep]
    filepath = delimiter.join(mylist)
    databk = h5file.get(filepath)
    rho = np.array(databk)
    rho = np.transpose(rho)
    rho = np.transpose(np.reshape(rho, (ny,nz)))
    rho_mark = copy.deepcopy(rho)
    
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
            rho_mark[bubRegion:bub_inc_loc, j] = np.nan
            rho_mark[bubRegion:bub_inc_loc, ny-j-1] = np.nan
            

                                   
                   
    from PIL import Image  
    import matplotlib                


    fig_pi = plt.figure(1, figsize = (5*1.214,5.0), dpi=200)
    plt.rc('font', family='serif', size=8)
    plt.imshow(rho, origin='none', extent=horizon_lim+vert_lim, aspect=1,
               cmap='coolwarm',vmin=rho_l,vmax=rho_h)
    #plt.colorbar()
    plt.savefig('ori.jpg')
 


    cmap = matplotlib.cm.coolwarm
    cmap.set_bad('black',1.)
    fig_mark = plt.figure(1, figsize = (5*1.214,5.0), dpi=200)
    plt.rc('font', family='serif', size=8)
    plt.imshow(rho_mark, origin='none', extent=horizon_lim+vert_lim, aspect=1,
               cmap=cmap,vmin=rho_l,vmax=rho_h)
    #plt.colorbar()
    plt.savefig('mark.jpg')
   


    ori = Image.open("ori.jpg")
    mark = Image.open("mark.jpg")

    ori = ori.convert("RGBA")
    mark = mark.convert("RGBA")
    new_img = Image.blend(ori, mark, 0.5)
    mylist = ['new', istep]
    filepath = delimiter.join(mylist)
    new_img.save(filepath,"PNG")
    
    

    
    

