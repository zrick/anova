#!/opt/local/bin/python 
import numpy as np 
from netCDF4 import Dataset 

pi2=2.*np.pi

nx=50
ny=50
nz=50
nt=50

fnc=Dataset('rand_{}_{}_{}_{}.nc'.format(nx,ny,nz,nt),'w',format='NETCDF4') 


fnc.createDimension('x',nx) 
fnc.createDimension('y',ny)
fnc.createDimension('z',nz) 
fnc.createDimension('t',nt) 

x=np.arange(0,pi2,pi2/nx) 
y=np.arange(0,pi2,pi2/ny) 
z=np.arange(0,pi2,pi2/nz) 
t=np.arange(0,pi2,pi2/nt) 
print(len(x)) 

vx=fnc.createVariable('x','f4',dimensions=('x')) 
vy=fnc.createVariable('y','f4',dimensions=('y')) 
vz=fnc.createVariable('z','f4',dimensions=('z')) 
vt=fnc.createVariable('t','f4',dimensions=('t')) 
vd=fnc.createVariable('data','f4',dimensions=('x','y','z','t')) 

fnc.variables['x'][:]=x 
fnc.variables['y'][:]=y 
fnc.variables['z'][:]=z 
fnc.variables['t'][:]=t 
var=0
min=0
max=0
for it in range(nt): 
    print('IT:',it) 
    for ix in range(nx):   
        rf=np.random.rand(ny,nz) 
        fnc.variables['data'][ix,:,:,it]=rf
        var=var+np.var(rf)  
        min=np.amin(rf) if np.amin(rf)<min else min
        max=np.amax(rf) if np.amax(rf)>max else max
print(min,max,var/(it*ix))

fnc.close()
