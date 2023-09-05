#!/usr/bin/env python3
import sys, os, datetime, re, shutil, netCDF4
from numpy import *
from scipy.special import logsumexp
from mpi4py import MPI

def write_netcdf(filename, lon, lat, c, e, b, o):
   os.system('rm -f '+filename)
   nc = netCDF4.Dataset(filename, 'w', format='NETCDF3_CLASSIC')
   nc.Conventions = 'COARDS'
   
   ny, nx = c.shape
   nc.createDimension('x', nx)
   nc.createDimension('y', ny)
   var = nc.createVariable('lon', 'f4', ('y','x'))
   var[:,:] = lon.astype(float32)
   var.units = 'degreesE'
   var = nc.createVariable('lat', 'f4', ('y','x'))
   var[:,:] = lat.astype(float32)
   var.units = 'degreesN'
   
   var = nc.createVariable('c', 'f4', ('y','x'))
   var[...] = c.astype(float32)
   var.units = 'mm'
   var.longname = 'Control'
   var.coordinates = 'lon lat'
   var = nc.createVariable('e', 'f4', ('y','x'))
   var[...] = e.astype(float32)
   var.units = 'mm'
   var.longname = 'Ensemble mean'
   var.coordinates = 'lon lat'
   var = nc.createVariable('b', 'f4', ('y','x'))
   var[...] = b.astype(float32)
   var.units = 'mm'
   var.longname = 'Barycenter'
   var.coordinates = 'lon lat'
   var = nc.createVariable('o', 'f4', ('y','x'))
   var[...] = o.astype(float32)
   var.units = 'mm'
   var.longname = 'Observation'
   var.coordinates = 'lon lat'

area = 'Kyushu'; analysis_date = '202007030600'
iexp = 100; nmember0 = 21
#experiment = 'Kyushu02km'; start = 6*60; end = 15*60
experiment = 'Fugaku05km'; start = 9*60; end = 18*60

epsilon = 0.0001 # enstropic regulation
phi = 0.99999 #phi = tau/(tau+epsilon) #tau: mass regulation
q = 2.
amin = 1.e-30 # minimum value of a to avoid nan
stopping = 1.e-4
niteration = 10000

# Directories
root_dir = '/home/leduc/project/ot4ens'
input_dir = root_dir
input_file = input_dir+'/'+area+'%2.2d'%iexp+'_'+'%2.2d'%(start/60)+'%2.2d'%(end/60)+'.nc'
output_dir = input_dir
output_file = output_dir+'/'+area+'%2.2d'%iexp+'mean_'+'%2.2d'%(start/60)+'%2.2d'%(end/60)+'.nc'

# Work dir
comm = MPI.COMM_WORLD
nprocessor = comm.Get_size()
myid = comm.Get_rank()
work_dir = root_dir+'/../work'
current_dir = os.getcwd()
os.chdir(work_dir)
if myid == 0: os.system('rm -rf *')

# Member
nmember = 0
for imember in range(1,nmember0):
   if (imember-1)%nprocessor != myid: continue
   nmember += 1
member = list(range(nmember))
kmember = 0
for imember in range(1,nmember0):
   if (imember-1)%nprocessor != myid: continue
   member[kmember] = '%4.4d'%imember
   kmember += 1
pass

# Read
nc = netCDF4.Dataset(input_file, 'r')
lon = nc.variables['lon'][:]
lat = nc.variables['lat'][:]
rain = nc.variables['rain'][:]
nc.close()
ny,nx = lon.shape
w = ones((nmember))/(nmember0-1)
logw = zeros((nmember,ny,nx))
for imember in range(nmember): logw[imember] = log(w[imember])
tmp = zeros((ny,nx))
c = zeros((ny,nx)) # control
o = zeros((ny,nx)) # obs
b = ones((2,ny,nx)) # barycenter
e = zeros((ny,nx)) # ensemble mean
a = zeros((nmember,ny,nx))
kmember = 0
for imember in range(1,nmember0):
   if (imember-1)%nprocessor != myid: continue
   a[kmember] = rain[imember]
   kmember += 1
c[:,:] = rain[0]
o[:,:] = rain[-1]
# Ensemble mean
e[:,:] = a.sum(axis=0)
comm.Allreduce(e, tmp, MPI.SUM)
e[:,:] = tmp/(nmember0-1)

# Cost
nxy=max(nx,ny)
Cx = zeros((nx,nx))
Cy = zeros((ny,ny))
for i in range(nx):
   for j in range(i+1,nx):
      Cx[i,j] = abs((i-j)/(nxy-1))**q
      Cx[j,i] = Cx[i,j]
Cx /= epsilon
for i in range(ny):
   for j in range(i+1,ny):
      Cy[i,j] = abs((i-j)/(nxy-1))**q
      Cy[j,i] = Cy[i,j]
Cy /= epsilon

# Sinkhorn
a += amin
loga = log(a)
f = zeros((nmember,ny,nx))
g = zeros((nmember,ny,nx))
if myid == 0: tmp = zeros((nprocessor,ny,nx))
for it in range(niteration):
   # Kv or f
   for imember in range(nmember):
      for j in range(ny):
         g[imember,j] = logsumexp(-Cx+g[imember,j,None,:],axis=1)
      for i in range(nx):
         g[imember,:,i] = logsumexp(-Cy+g[imember,None,:,i],axis=1)
      f[imember] = phi*(loga[imember]-g[imember])
   # K.T*u or g
   for imember in range(nmember):
      for j in range(ny):
         f[imember,j] = logsumexp(-Cx+f[imember,j,:,None],axis=0)
      for i in range(nx):
         f[imember,:,i] = logsumexp(-Cy+f[imember,:,None,i],axis=0)
   b[0] = b[1]
   b[1] = logsumexp((1-phi)*f+logw,axis=0)
   if myid == 0:
      tmp[0] = b[1]
      for iprocessor in range(1,nprocessor):
         k = comm.recv(source=iprocessor, tag=0)
         tmp[k] = comm.recv(source=iprocessor, tag=1)
      b[1] = 1/(1-phi)*logsumexp(tmp,axis=0)
   else:
      comm.send(myid,dest=0,tag=0)
      comm.send(b[1],dest=0,tag=1)
   b[1] = comm.bcast(b[1],root=0)
   if it > 0: ratio = abs(b[1]-b[0]).max()/max(abs(b[1]).max(),abs(b[0]).max())
   if it > 0 and myid == 0: print('barycenter:', it, ratio, exp(b[1].min()), exp(b[1].max()))
   if it > 0 and ratio < stopping: break
   for imember in range(nmember): g[imember] = phi*(b[1]-f[imember])
b[:,:,:] = exp(b)

if myid == 0:
   b -= amin
   invalid = b < 0.; b[invalid] = 0.
   write_netcdf(output_file, lon, lat, c, e, b[1], o)
comm.Barrier()

