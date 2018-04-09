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

curfilePath = os.path.abspath(__file__)
curDir = os.path.abspath(os.path.join(curfilePath,os.pardir))
parentDir = os.path.abspath(os.path.join(curDir,os.pardir)) 

def high_order_gradient(fx,dx,order):
    length=len(fx)
    fxgrad = np.zeros(length)
    if order==4: 
        for i in range(2):
            fxgrad[i]=(-25*fx[i]+48*fx[i+1]-36*fx[i+2]+16*fx[i+3]-3*fx[i+4])/(12*dx)
        for i in range(2,length-2):
            fxgrad[i]=(-fx[i+2]+fx[i+1]*8-fx[i-1]*8+fx[i-2])/(12*dx)
        for i in range(length-2,length):
            fxgrad[i]=(25*fx[i]-48*fx[i-1]+36*fx[i-2]-16*fx[i-3]+3*fx[i-4])/(12*dx)
    if order==6:
        for i in range(3):
            fxgrad[i]=(-49/20*fx[i]+6*fx[i+1]-15/2*fx[i+2]+20/3*fx[i+3]-15/4*fx[i+4]+6/5*fx[i+5]-1/6*fx[i+6])/(dx)
        for i in range(3,length-3):
            fxgrad[i]=(fx[i+3]-9*fx[i+2]+45*fx[i+1]-45*fx[i-1]+9*fx[i-2]-fx[i-3])/(60*dx)
        for i in range(length-3,length):
            fxgrad[i]=(49/20*fx[i]-6*fx[i-1]+15/2*fx[i-2]-20/3*fx[i-3]+15/4*fx[i-4]-6/5*fx[i-5]+1/6*fx[i-6])/(dx)

        
    return fxgrad

inFile="tests_single_new.h5"
mylist = [parentDir,'/',inFile]
delimiter = ''
filepath = delimiter.join(mylist)
#nz enlarged only 
h5file = h5py.File(filepath,'r+')

variable = ['PVx','PVy','PVz','PPress', 'Prho']
mylist = ['Fields/','Prho','/','600000']
filepath = delimiter.join(mylist)
databk = h5file.get(filepath)
rho = np.array(databk)

nz=rho.shape[0]
ny=rho.shape[1]
nx=rho.shape[2]
#specify inout parameters here
gamma=5.0/3.0
g=1.0
inFile="tests_single_new.h5"
Ly=0.4
Lz=3.2
CFDmethod = False
dz=dy=dx=Lz/nz
rho_l=0.4
rho_h=1.6
rho_inc = (rho_h + rho_l)/2.0
winPercent = 0.1
winPoint =  int(winPercent*nz)
#input done

mylist = [parentDir,'/',inFile]
delimiter = ''
filepath = delimiter.join(mylist)

rho = np.transpose(rho)
rho = np.transpose(np.reshape(rho, (ny,nz)))
rho_mark = copy.deepcopy(rho)
horizon_lim = (0, Ly-dy)
vert_lim = (0, Lz-dz)
extent=horizon_lim+vert_lim



dx=1.0
specout = 1000
step = []
totalstep=839626
for i in range(totalstep/specout):
    step.append(str((i+1)*specout).zfill(6))
step = ['133000']

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
plt.show()


cmap = matplotlib.cm.coolwarm
cmap.set_bad('black',1.)
fig_mark = plt.figure(1, figsize = (5*1.214,5.0), dpi=200)
plt.rc('font', family='serif', size=8)
plt.imshow(rho_mark, origin='none', extent=horizon_lim+vert_lim, aspect=1,
           cmap=cmap,vmin=rho_l,vmax=rho_h)
#plt.colorbar()
plt.savefig('mark.jpg')
plt.show()


ori = Image.open("ori.jpg")
mark = Image.open("mark.jpg")

ori = ori.convert("RGBA")
mark = mark.convert("RGBA")
new_img = Image.blend(ori, mark, 0.5)
new_img.save("new.png","PNG")

plt.show()
