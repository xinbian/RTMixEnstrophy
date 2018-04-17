# -*- coding: utf-8 -*-
"""
To calculate 3D mixing layer width and enstrophy,
mixing layer defination:
https://search.proquest.com/docview/194682903?pq-origsite=gscholar
@author: Xin
"""
import matplotlib
matplotlib.use('Agg')
import pylab
import matplotlib.pyplot as plt
import matplotlib
import h5py
import numpy as np
import os
import os.path
from CFDmodule import *
from readStep import *

#curfilePath = os.path.abspath(__file__)
#curDir = os.path.abspath(os.path.join(curfilePath,os.pardir))
#parentDir = os.path.abspath(os.path.join(curDir,os.pardir)) 

#specify inout parameters here
gamma=5.0/3.0
g=1.0
Lz=3.2
waveLen = 0.4
mu =1.13137E-4
#select one if calculate enstrophy
CFDmethod = False
npCalGrad = False
showPerct = False
calcEnergy = False
calcMix = False
outPut = False
rhoH = 1.0833
rhoL = 1.0
skip = 20
#####input done

dz=dy=dx=Lz/nz

if nx == 1:
    dx=1.0

FieldPoint = h5file.get('Fields').values()
totalsteps, specout = get_LatestTime(FieldPoint)

seq = 0
step = []
for i in range(totalsteps/specout):
    step.append(str((i+1)*specout).zfill(6))

#initialize time dependent mixing layer width, KE, PE, enstropy
rho_avr=np.zeros(nz, dtype=np.float64)   
h = np.zeros(len(step))
ie = np.zeros(len(step))
ke = np.zeros(len(step))
pe = np.zeros(len(step))
enstropy = np.zeros(len(step))  
sum_x = np.zeros(len(step))  
dissRate = np.zeros(len(step))

bub_loc_all = np.zeros(len(step))
bub_loc_all_ori = np.zeros(len(step))
sp_loc_all = np.zeros(len(step))
bub_velo_all = np.zeros(len(step))
bub_velo_all_aver = np.zeros(len(step))
sp_velo_all = np.zeros(len(step))
bub_velo_all_ori = np.zeros(len(step))



if CFDmethod == True:
	CFD_x = Create_matrix_fd2(nx) / dx
	CFD_y = Create_matrix_fd2(ny) / dy
	CFD_z = Create_matrix_fd2(nz) / dz


#calculate mixing layer width
for istep in step:
	if showPerct == True:
		print "doing", float(istep)/totalsteps*100, "%"
	if calcEnergy == True:
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
		ke[seq] = 0.5 * np.sum((vx**2 + vy**2 + vz**2) * rho)*dx*dy*dz
		#internal energy
		mylist = ['Fields/',variable[3],'/',istep]
		filepath = delimiter.join(mylist)
		databk = h5file.get(filepath)
		press = np.array(databk)
		ie[seq] =np.sum((1/(gamma-1))*press)*dx*dy*dz
	#enstrophy
	if CFDmethod == True:
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
		
		dxVx = calc_x_deri(vx, CFD_x)
		dyVy = calc_y_deri(vy, CFD_y)
		dzVz = calc_z_deri(vz, CFD_z)
		
	        Wx = dyVz - dzVy
	        Wy = dzVx - dxVz
	        Wz = dxVy - dyVx
		
		Sxx = dxVx
		Sxy = (dxVy + dyVx)/2.0
		Sxz = (dzVx + dxVz)/2.0
		Syy = dyVy
		Syz = (dzVy + dyVz)/2.0
		Szz = dzVz
		
		if nx == 1:
			ndim = 2.0
		else:
			ndim = 3.0
			
	        EnstrophyCFD = 0.5*(np.multiply(Wx,Wx)+np.multiply(Wy,Wy)+np.multiply(Wz,Wz))
	        sum_x[seq] = np.sum(EnstrophyCFD)*dx*dy*dz
		dissRate[seq] = 2*mu*dx*dy*dz*np.sum((Sxx**2 + Sxy**2 + Sxz**2 + Sxy**2 + Syy**2 + Syz**2 + \
				Sxz**2 + Syz**2 + Szz**2 - (Sxx**2 + Syy**2 + Szz**2)/ndim))
	if npCalGrad == True:
		dyVx = np.gradient(vx, dy, axis=1)
	        dzVx = np.gradient(vx, dz, axis=0)
	        dxVy = np.gradient(vy, dx, axis=2)
	        dzVy = np.gradient(vy, dz, axis=0)
	        dxVz = np.gradient(vz, dx, axis=2)
	        dyVz = np.gradient(vz, dy, axis=1)
		
		dxVx = np.gradient(vx, dx, axis=2)
		dyVy = np.gradient(vy, dy, axis=1)
		dzVz = np.gradient(vz, dz, axis=0)
	
		Wx = dyVz - dzVy
		Wy = dzVx - dxVz
		Wz = dxVy - dyVx
		
		Sxx = dxVx
		Sxy = (dxVy + dyVx)/2.0
		Sxz = (dzVx + dxVz)/2.0
		Syy = dyVy
		Syz = (dzVy + dyVz)/2.0
		Szz = dzVz
		
		
		if nx == 1:
			ndim = 2.0
		else:
			ndim = 3.0
		
		enstropy[seq] = 0.5*np.sum(Wx**2+Wy**2+Wz**2)*dx*dy*dz 


        mylist = ['Fields/', 'Prho', '/', istep]
        filepath = delimiter.join(mylist)
        databk = h5file.get(filepath)
        np_data = np.array(databk)
        if nx == 1:
    	    m1 = (np_data[:, ny/2-1, 0] + np_data[:, ny/2, 0] )/2
        else:
    	    m1 = (np_data[:, ny/2-1, nx/2-1] + np_data[:, ny/2, nx/2] 
    	  + np_data[:, ny/2-1, nx/2] + np_data[:, ny/2, nx/2-1])/4.0
        m2 = np_data[:, 0, 0]
        m1_filter=m1.copy();    
        m2_filter=m2.copy();    
        for jstep in range(2,nz-3):
            m1_filter[jstep]=(m1[jstep-2]+m1[jstep-1]+m1[jstep]+m1[jstep+1]+m1[jstep+2])/5;
            m2_filter[jstep]=(m2[jstep-2]+m2[jstep-1]+m2[jstep]+m2[jstep+1]+m2[jstep+2])/5;
    
    
        m1_grad = high_order_gradient(m1_filter,dx,6)
        m2_grad = high_order_gradient(m2_filter,dx,6)
    
    
        sp_loc = np.argmax(m1_grad)
        bub_loc = np.argmax(m2_grad)
    
        sp_loc_all[seq] = sp_loc
        bub_loc_all[seq] = bub_loc
    
    
        mylist = ['Fields/', 'PVz', '/', istep]
        filepath = delimiter.join(mylist)
        databk = h5file.get(filepath)
        np_data = np.array(databk)
      
        #consider 2D/3D case
        if nx == 1:
    	    m1 = (np_data[:, ny/2-1, 0] + np_data[:, ny/2, 0] )/2
        else:
    	    m1 = (np_data[:, ny/2-1, nx/2-1] + np_data[:, ny/2, nx/2] 
    	  + np_data[:, ny/2-1, nx/2] + np_data[:, ny/2, nx/2-1])/4.0
        m2 = np_data[:, 0, 0]
        sp_velo = m1[sp_loc]
        bub_velo = m2[bub_loc]
        bub_velo_all[seq] = bub_velo
        sp_velo_all[seq] = sp_velo
    
    	
        seq += 1
    
#normalize pe/calculate lost pe to this time
pe = pe*dx*dy*dz

#normalize h
h=h*dz/waveLen
#output
all_data = np.column_stack((np.asarray(step),h, ke, pe[0]-pe, ie-ie[0], enstropy, sum_x, dissRate))
np.savetxt('savedMixAndEnstro', all_data,delimiter='\t',fmt='%s')

all_data = np.column_stack((bub_loc_all,bub_velo_all,sp_loc_all,sp_velo_all))
np.savetxt('saved_bub_velo',all_data,delimiter='\t',fmt='%s')

plt.plot(bub_velo_all)
plt.plot(sp_velo_all)
plt.savefig('vel.eps', format='eps', dpi=300)
plt.clf()

plt.plot(bub_loc_all)
plt.plot(sp_loc_all)
plt.savefig('loc.eps', format='eps', dpi=300)
plt.clf()

#mixing layer
if calcMix == True:	
	f=open('rho1d.txt','w')
	
	step=[]
	for i in range(totalsteps/specout+1):
	    step.append(str((i+1)*specout).zfill(6))
	
	
	for ii in range(1,totalsteps+1):
	    f.write("%i \n" % ii)
	    if ii==1:
	        rho_avr[0:nz/2-1]=rhoL
	        rho_avr[nz/2:]=rhoH
	        np.savetxt(f,rho_avr)
	    if ii%specout==0:            
	                  rho_avr=np.zeros(nz, dtype=np.float64)       
	                  delimiter = ''
	                  mylist = ['Fields/','Prho','/',step[ii/specout-1]]
	                  filepath = delimiter.join(mylist)
	                  databk = h5file.get(filepath)
	                  np_data = np.array(databk)
	                  m1=np_data
	                  for i in range(nz):
	                      rho_avr[i]=np.mean(m1[i,:,:])
	                  np.savetxt(f,rho_avr)
            
	f.close()
h5file.close()



if outPut == True:
	plt.plot(np.asarray(step),h)
	plt.title('mixing layer width vs time step')
	plt.savefig('mix.eps', format='eps', dpi=1000)
	plt.show()
	plt.plot(np.asarray(step),ke , label='KE')
	plt.plot(np.asarray(step),pe[0]-pe, label='released PE')
	plt.plot(np.asarray(step),ie-ie[0], label='increased IE')
	#plt.plot(np.asarray(step),ke+pe,label='KE+PE')
	plt.plot(np.asarray(step),pe[0]-pe-(ie-ie[0]), label='released PE -increased IE')
	plt.title('energy vs time step')
	pylab.legend(loc='best')
	plt.savefig('energy.eps', format='eps', dpi=1000)
	plt.show()
	plt.plot(np.asarray(step), enstropy)
	plt.title('enstrophy vs time step')
	plt.savefig('enstropy.eps', format='eps', dpi=1000)
	plt.show()
	   
	    
	  
