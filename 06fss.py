#!/usr/bin/env python3
import sys, os, datetime, re, shutil, netCDF4
from mri_mapproj import mapproj
from verification import verification
from nhmlib import *
from numpy import *
from mpi4py import MPI

missing_value = -9999.

def write_field(filename, lon, lat, field):
   os.system('rm -f '+filename)
   nc = netCDF4.Dataset(filename, 'w', format='NETCDF3_CLASSIC')
   nx = len(lon)
   ny = len(lat)
   nc.createDimension('x', nx)
   nc.createDimension('y', ny)
   var = nc.createVariable('x', 'f4', ('x'))
   var[:] = lon.astype(float32)
   var.units = 'degreesE'
   var = nc.createVariable('y', 'f4', ('y'))
   var[:] = lat.astype(float32)
   var.units = 'degreesN'
   var = nc.createVariable('field', 'f4', ('y','x'))
   var[...] = field.astype(float32)
   var.missing_value = missing_value
   nc.close()

def write_netcdf(filename, time, timescale, scale, threshold):
   os.system('rm -f '+filename)
   nc = netCDF4.Dataset(filename, 'w', format='NETCDF3_CLASSIC')
   nc.Conventions = 'COARDS'
   
   ntime, ntimescale, nscale, nthreshold = pf2.shape
   nc.createDimension('threshold', nthreshold)
   nc.createDimension('scale', nscale)
   nc.createDimension('timescale', ntimescale)
   nc.createDimension('time', ntime)
   
   var = nc.createVariable('threshold', 'f4', ('threshold'))
   var[:] = array(threshold).astype(float32)
   var = nc.createVariable('scale', 'f4', ('scale'))
   var[:] = scale.astype(float32)
   var = nc.createVariable('timescale', 'f4', ('timescale'))
   var[:] = timescale.astype(float32)
   var = nc.createVariable('time', 'i2', ('time'))
   var[:] = array(time).astype(int16)
   
   var = nc.createVariable('nfss', 'i4', ('time'))
   var[:] = nfss.astype(int32)
   var = nc.createVariable('npo', 'i4', ('time','timescale'))
   var[:] = npo.astype(int32)
   var = nc.createVariable('po', 'f4', ('time','timescale','threshold'))
   var[:] = po.astype(float32)
   var = nc.createVariable('pfpo', 'f4', ('time','timescale','scale','threshold'))
   var[:] = pfpo.astype(float32)
   var = nc.createVariable('pf2', 'f4', ('time','timescale','scale','threshold'))
   var[:] = pf2.astype(float32)

area = 'Kyushu'; variable = 'e'; analysis_date = '202007030600'
iexp = 100; nmember0 = 21
#iexp = 2; nmember0 = 1001
#experiment = 'Kyushu02km'; start = 6*60; end = 15*60
experiment = 'Fugaku05km'; resolution = 5.; start = 12*60; end = 18*60; period = 3*60
scale = array([30.,45])
threshold = [1.,5.,10.,25.,50.,75.,100.]
timescale = array([1,]).astype(int32)
time = list(range(start,end+1,period))
ntime = len(time)

# Directories
root_dir = '/home/leduc/project/uot4ens'
nhm_dir = root_dir+'/../../nhm'
const_dir = nhm_dir+'/exp/'+experiment+'/const'
input_dir = root_dir
output_dir = input_dir
output_file = output_dir+'/'+area+'%2.2d'%iexp+'fss'+variable+'_'+'%2.2d'%((start-period)/60)+'%2.2d'%(end/60)+'.nc'
# Work dir
work_dir = root_dir+'/../tmp'
comm = MPI.COMM_WORLD
nprocessor = comm.Get_size()
myid = comm.Get_rank()
current_dir = os.getcwd()
os.chdir(work_dir)
if myid == 0: os.system('rm -rf *')
comm.Barrier()

# Inventory
input_file = input_dir+'/'+area+'%2.2d'%iexp+'mean_'+'%2.2d'%((start-period)/60)+'%2.2d'%(start/60)+'.nc'
nc = netCDF4.Dataset(input_file, 'r')
lon = nc.variables['lon'][:]
lat = nc.variables['lat'][:]
nc.close()
ny,nx = lon.shape
domain = ones((ny,nx))
fcst = zeros((ntime,ny,nx))
obs = zeros((ntime,ny,nx))
valid = ones((ntime,ny,nx),dtype=bool)

# Read
for itime in range(ntime):
   if itime%nprocessor != myid: continue
   input_file = input_dir+'/'+area+'%2.2d'%iexp+'mean_'+'%2.2d'%((time[itime]-period)/60)+'%2.2d'%(time[itime]/60)+'.nc'
   if not os.path.isfile(input_file): sys.exit(0)
   nc = netCDF4.Dataset(input_file, 'r')
   fcst[itime] = nc.variables[variable][:]
   obs[itime] = nc.variables['o'][:]
   nc.close()
   invalid = fcst[itime] < 0.
   fcst[itime,invalid] = 0.
   invalid = obs[itime] < 0.
   obs[itime,invalid] = 0.
pass

# Variables
nthreshold = len(threshold)
nscale = len(scale)
ntimescale = len(timescale)
pfpo = zeros((ntime,ntimescale,nscale,nthreshold))
pf2 = zeros((ntime,ntimescale,nscale,nthreshold))
po = zeros((ntime,ntimescale,nthreshold))
npo = zeros((ntime,ntimescale)).astype(int32)
nfss = zeros((ntime)).astype(int32)
valid_int = valid.astype(int32)

# Run
for itime in range(ntime):
   if itime%nprocessor != myid: continue
   nfss[itime] = sum(valid[itime])
   for ithreshold in range(nthreshold):
      print('Thread '+str(myid)+' computes FSS at '+str(time[itime])+' with threshold '+str(threshold[ithreshold])+'mm')
      # Observation fraction of entire domain
      npo[itime], po[itime,:,ithreshold] = verification.compute_po(itime+1, threshold[ithreshold], timescale, valid_int, obs)
      # Observation fraction with scales
      po_tmp = verification.compute_fraction(itime+1, resolution, threshold[ithreshold], timescale, scale, valid_int, obs)
      po2_tmp = po_tmp[:,:,valid[itime]]**2
      po2_tmp = sum(po2_tmp,axis=2)
      # FSS
      pf_tmp = verification.compute_fraction(itime+1, resolution, threshold[ithreshold], timescale, scale, valid_int, fcst)
      pf2_tmp = pf_tmp[:,:,valid[itime]]**2
      pf2[itime,:,:,ithreshold] = sum(pf2_tmp,axis=2) + po2_tmp
      pfpo_tmp = pf_tmp[:,:,valid[itime]]*po_tmp[:,:,valid[itime]]
      pfpo[itime,:,:,ithreshold] = sum(pfpo_tmp,axis=2)
      #print(itime, pf_tmp.min(), pf_tmp.max())
comm.Barrier()

# Gather data
for itime in range(ntime):
   if itime%nprocessor != myid: continue
   if myid == 0:
      for iprocessor in range(1,nprocessor):
         if itime+iprocessor > ntime-1: break
         i = comm.recv(source=iprocessor, tag=0)
         nfss[i] = comm.recv(source=iprocessor, tag=1)
         npo[i] = comm.recv(source=iprocessor, tag=2)
         po[i] = comm.recv(source=iprocessor, tag=3)
         pfpo[i] = comm.recv(source=iprocessor, tag=4)
         pf2[i] = comm.recv(source=iprocessor, tag=5)
   else:
      comm.send(itime,dest=0,tag=0)
      comm.send(nfss[itime],dest=0,tag=1)
      comm.send(npo[itime],dest=0,tag=2)
      comm.send(po[itime],dest=0,tag=3)
      comm.send(pfpo[itime],dest=0,tag=4)
      comm.send(pf2[itime],dest=0,tag=5)
comm.Barrier()

# Write daily files for using in bootstraping
if myid == 0:
   print('Write file '+output_file)
   #if not os.path.isdir(output_dir): os.system('mkdir -p '+output_dir)
   write_netcdf(output_file, time, timescale*period, scale, threshold)
comm.Barrier()
os.chdir(current_dir)
