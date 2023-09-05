#!/usr/bin/env python3
import sys, os, datetime, re, shutil, random, netCDF4
from scipy.stats import scoreatpercentile, percentileofscore
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

def write_netcdf(filename, quantile, experiment, timescale, scale, threshold):
   os.system('rm -f '+filename)
   nc = netCDF4.Dataset(filename, 'w', format='NETCDF3_CLASSIC')
   nc.Conventions = 'COARDS'
   #nc.lag_hours = lag_hours
   
   nquantile, nexp, nexp, ntimescale, nscale, nthreshold = fss_quantile.shape
   nc.createDimension('threshold', nthreshold)
   nc.createDimension('scale', nscale)
   nc.createDimension('timescale', ntimescale)
   nc.createDimension('experiment', nexp)
   nc.createDimension('quantile', nquantile)
   
   var = nc.createVariable('threshold', 'f4', ('threshold'))
   var[:] = array(threshold).astype(float32)
   var = nc.createVariable('scale', 'f4', ('scale'))
   var[:] = scale.astype(float32)
   var = nc.createVariable('timescale', 'f4', ('timescale'))
   var[:] = timescale.astype(float32)
   var = nc.createVariable('experiment', 'i1', ('experiment'))
   for iexp in range(nexp): var.setncattr('experiment'+'%2.2d'%iexp, experiment[iexp])
   var[:] = arange(nexp, dtype=int8)
   var = nc.createVariable('quantile', 'f4', ('quantile'))
   var[:] = quantile.astype(float32)
   
   var = nc.createVariable('fss', 'f4', ('experiment','timescale','scale','threshold'))
   var[...] = fss[-1].astype(float32)
   var.missing_value = missing_value
   var = nc.createVariable('fss_quantile', 'f4', ('quantile','experiment','experiment','timescale','scale','threshold'))
   var[...] = fss_quantile.astype(float32)
   var.missing_value = missing_value
   var = nc.createVariable('fss_confidence', 'f4', ('experiment','experiment','timescale','scale','threshold'))
   var[...] = fss_confidence.astype(float32)
   var.missing_value = missing_value

experiment = 'Kyushu'; iexp = 100
name = ['c', 'b']
start = 12*60; end = 18*60; period = 3*60
quantile = array([2.5,5.0,50.0,95.0,97.5])
nquantile = len(quantile)

# Directories
root_dir = '/home/leduc/project/uot4ens'
input_dir = root_dir
output_dir = root_dir
output_file = output_dir+'/'+experiment+'%2.2d'%iexp+'fss'+name[1]+'vs'+name[0]+'_'+'%2.2d'%((start-period)/60)+'%2.2d'%(end/60)+'.nc'
work_dir = root_dir+'/../tmp'

# Work dir
comm = MPI.COMM_WORLD
nprocessor = comm.Get_size()
myid = comm.Get_rank()
work_dir = work_dir+'/tmp'+'%2.2d'%myid
os.system('mkdir -p '+work_dir)
os.system('rm -f '+work_dir+'/*')
current_dir = os.getcwd()
os.chdir(work_dir)
os.system('rm -rf *')

# Find dates
nexp = len(name)
tail = '_'+'%2.2d'%((start-period)/60)+'%2.2d'%(end/60)+'.nc'
nfile = 1
nbootstrap = 1000*nfile
if nfile == 1: nbootstrap = 0

# Inventory
nc = netCDF4.Dataset(input_dir+'/'+experiment+'%2.2d'%iexp+'fss'+name[0]+tail, 'r')
timescale = nc.variables['timescale'][:]
scale = nc.variables['scale'][:]
threshold = nc.variables['threshold'][:]
ntimescale = len(timescale)
nscale = len(scale)
nthreshold = len(threshold)
nc.close()

# Set istart, iend
istart = list(range(nexp))
iend = list(range(nexp))
for i in range(nexp):
   nc = netCDF4.Dataset(input_dir+'/'+experiment+'%2.2d'%iexp+'fss'+name[i]+tail, 'r')
   time = nc.variables['time'][:]
   ntime = len(time)
   for itime in range(ntime):
      if time[itime] == start: istart[i] = itime
      if time[itime] == end: iend[i] = itime
   nc.close()
pass

# File variables
pfpo_file = zeros((nexp,nfile,ntimescale,nscale,nthreshold))
pf2_file = zeros((nexp,nfile,ntimescale,nscale,nthreshold))
nfss_file = zeros((nexp,nfile)).astype(int32)
for j in range(nexp):
   for i in range(nfile):
      #print input_dir+'/'+experiment[j]+'/forecast/'+experiment[j]+iexp[j]+'/'+input_file[i]
      nc = netCDF4.Dataset(input_dir+'/'+experiment+'%2.2d'%iexp+'fss'+name[j]+tail, 'r')
      data = nc.variables['pfpo'][:]
      pfpo_file[j,i] = sum(data[istart[j]:iend[j]+1],axis=0)
      data = nc.variables['pf2'][:]
      pf2_file[j,i] = sum(data[istart[j]:iend[j]+1],axis=0)
      data = nc.variables['nfss'][:]
      nfss_file[j,i] = sum(data[istart[j]:iend[j]+1],axis=0)
      nc.close()
#sys.exit(0)

# Variables
pfpo = zeros((nexp,ntimescale,nscale,nthreshold))
pf2 = zeros((nexp,ntimescale,nscale,nthreshold))
fss = zeros((nbootstrap+1,nexp,ntimescale,nscale,nthreshold)) # -1 is for the normal calculation
fss_confidence = zeros((nexp,nexp,ntimescale,nscale,nthreshold))
fss_quantile = zeros((nquantile,nexp,nexp,ntimescale,nscale,nthreshold))

# Bootstraping
for i in range(nbootstrap+1):
   if i%nprocessor != myid: continue
   #print 'Bootstrap: ', i
   pfpo[...] = 0.
   pf2[...] = 0.
   nfss = 0
   # FSS
   for j in range(nfile):
      if i == nbootstrap: k = j
      else: k = random.randint(0,nfile-1)
      pfpo += pfpo_file[:,k]
      pf2 += pf2_file[:,k]
      nfss += nfss_file[:,k]
   for iexp in range(nexp):
      invalid = pf2[iexp] <= 0.
      valid = logical_not(invalid)
      fss[i,iexp,valid] = 2*pfpo[iexp,valid]/pf2[iexp,valid]
      fss[i,iexp,invalid] = missing_value
      #print iexp, (pf2[iexp,0,0]-2*pfpo[iexp,0,0])/nfss[iexp], fss[i,iexp,0,0]
   # Gather data
   if myid == 0:
      for iprocessor in range(1,nprocessor):
         if i+iprocessor > nbootstrap: break
         j = comm.recv(source=iprocessor, tag=0)
         fss[j] = comm.recv(source=iprocessor, tag=3)
   else:
      comm.send(i,dest=0,tag=0)
      comm.send(fss[i],dest=0,tag=3)
comm.Barrier()

# Confidence interval
if myid == 0:
   for itimescale in range(ntimescale):
      for iscale in range(nscale):
         for ithreshold in range(nthreshold):
            # FSS
            for j in range(nexp):
               fss1 = fss[:-1,j,itimescale,iscale,ithreshold]
               for k in range(nexp):
                  if k == j: continue
                  fss2 = fss[:-1,k,itimescale,iscale,ithreshold]
                  diff = fss2 - fss1
                  valid = logical_and(fss2>missing_value+1, fss1>missing_value+1)
                  nvalid = sum(valid)
                  for iquantile in range(nquantile):
                     if nvalid > 2: fss_quantile[iquantile,j,k,itimescale,iscale,ithreshold] = scoreatpercentile(diff[valid],quantile[iquantile])
                     else: fss_quantile[iquantile,j,k,itimescale,iscale,ithreshold] = missing_value
                  if nvalid > 2: fss_confidence[j,k,itimescale,iscale,ithreshold] = 100. - percentileofscore(diff[valid],0.,kind='weak')
                  else: fss_confidence[j,k,itimescale,iscale,ithreshold] = missing_value
   print('Write file '+output_file)
   write_netcdf(output_file, quantile, name, timescale, scale, threshold)
#sys.exit(0)
comm.Barrier()
os.chdir(current_dir)
#shutil.rmtree(work_dir)
