#!/usr/bin/env python3
import sys, os, datetime, re, shutil, netCDF4
from numpy import *
from scipy.special import logsumexp
from mpi4py import MPI

def write_netcdf(filename, c, e, b, o):
   os.system('rm -f '+filename)
   nc = netCDF4.Dataset(filename, 'w', format='NETCDF3_CLASSIC')
   nc.Conventions = 'COARDS'
   nc.createDimension('time', ntime)
   var = nc.createVariable('c', 'f4', ('time'))
   var[...] = c.astype(float32)
   var.units = 'mm'
   var.longname = 'Control'
   var = nc.createVariable('e', 'f4', ('time'))
   var[...] = e.astype(float32)
   var.units = 'mm'
   var.longname = 'Ensemble mean'
   var = nc.createVariable('b', 'f4', ('time'))
   var[...] = b.astype(float32)
   var.units = 'mm'
   var.longname = 'Barycenter'
   var = nc.createVariable('o', 'f4', ('time'))
   var[...] = o.astype(float32)
   var.units = 'mm'
   var.longname = 'Observation'

area = 'Ichifusa'; 
iexp = 2; nmember0 = 1001
#iexp = 4; nmember0 = 101
#iexp = 100; nmember0 = 21
#experiment = 'Kyushu02km'
experiment = 'Fugaku05km'; analysis_date = '202007030900'

epsilon = 0.0001 # enstropic regulation
phi = 0.99999 #phi = tau/(tau+epsilon) tau: mass regulation
q = 2
amin = 1.e-30 # minimum value of a to avoid nan
stopping = 1.e-5
niteration = 10000

# Directories
root_dir = '/home/leduc/project/ot4ens'
input_dir = root_dir
input_file = input_dir+'/'+area+'%2.2d'%iexp+'.nc'
output_dir = input_dir
output_file = output_dir+'/'+area+'%2.2d'%iexp+'mean'+'.nc'

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
rain = nc.variables['rain'][:]
nc.close()
n,ntime = rain.shape
w = ones((nmember))/(nmember0-1)
logw = zeros((nmember,ntime))
for imember in range(nmember): logw[imember] = log(w[imember])
tmp = zeros((ntime))
c = zeros((ntime)) # control
o = zeros((ntime)) # obs
b = ones((2,ntime)) # barycenter
e = zeros((ntime)) # ensemble mean
a = zeros((nmember,ntime))
kmember = 0
for imember in range(1,nmember0):
   if (imember-1)%nprocessor != myid: continue
   a[kmember] = rain[imember]
   kmember += 1
c[:] = rain[0]
o[:] = rain[-1]
# Ensemble mean
e[:] = a.sum(axis=0)
comm.Allreduce(e, tmp, MPI.SUM)
e[:] = tmp/(nmember0-1)

# Cost
C = zeros((ntime,ntime))
for i in range(ntime):
   for j in range(i+1,ntime):
      C[i,j] = abs((i-j)/(ntime-1))**q
      C[j,i] = C[i,j]
C /= epsilon

# Sinkhorn
a += amin
loga = log(a)
f = zeros((nmember,ntime))
g = zeros((nmember,ntime))
if myid == 0: tmp = zeros((nprocessor,ntime))
for it in range(niteration):
   # Kv or f
   for imember in range(nmember):
      f[imember] = phi*(loga[imember]-logsumexp(-C+g[imember,None,:],axis=1))
   # Ku or g
   for imember in range(nmember):
      g[imember] = logsumexp(-C+f[imember,:,None],axis=0)
   b[0] = b[1]
   b[1] = logsumexp((1-phi)*g+logw,axis=0)
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
   for imember in range(nmember): g[imember] = phi*(b[1]-g[imember])
b[:,:] = exp(b)
   
if myid == 0:
   b -= amin
   invalid = b < 0.; b[invalid] = 0.
   write_netcdf(output_file, c, e, b[1], o)
comm.Barrier()
