# -*- coding: utf-8 -*-
"""
module for reading time step
@author: Xin
"""
import h5py
import numpy as np
import os
import os.path

curfilePath = os.path.abspath(__file__)
curDir = os.path.abspath(os.path.join(curfilePath,os.pardir))
parentDir = os.path.abspath(os.path.join(curDir,os.pardir)) 

inFile='tests_single_new.h5'
mylist = [parentDir,'/',inFile]
delimiter = ''
filepath = delimiter.join(mylist)
h5file = h5py.File(filepath,'r')

#read dataset dimensions
mylist = ['Fields/','Prho','/','002000']
filepath = delimiter.join(mylist)
databk = h5file.get(filepath)
m1 = np.array(databk)
nz=m1.shape[0]
ny=m1.shape[1]
nx=m1.shape[2]



def get_timestepstr(dset):

    return os.path.split(dset.name)[1]

def get_LatestTime(ChkPoints):

	maxtstep = 0
	for grp in ChkPoints:
		dsets = grp.values()
		i = 0
		specout = 0
		for dset in dsets:
			dTimestep = int(get_timestepstr(dset))
			if dTimestep > maxtstep:
			    maxtstep = dTimestep
			    
			if i == 2:
			    specout = dTimestep
			if i == 3:
			    specout = dTimestep - specout
			i += 1
	return maxtstep, specout

