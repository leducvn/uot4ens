#!/usr/bin/env python3
import sys, os, datetime, re, shutil, netCDF4
from numpy import *
from scipy.special import logsumexp

def plot_histogram(hour, a, e, cumulative=False, figtype=''):
   nmember, ntime = a.shape
   nexp, ntime = e.shape
   left = 0.; right = ntime-1
   bottom = 0.; top = a.max()+0.01
   dx = '5'; dy = '5'; y_text = ''
   if cumulative: dy = '50'
   
   os.system('gmtset PAPER_MEDIA a3')
   os.system('gmtset X_ORIGIN 3c')
   os.system('gmtset Y_ORIGIN 3c')
   os.system('gmtset BASEMAP_TYPE PLAIN')
   os.system('gmtset TICK_LENGTH 0')
   os.system('gmtset ANNOT_FONT_SIZE 14p')
   os.system('gmtset LABEL_FONT_SIZE 16p')
   os.system('gmtset LABEL_OFFSET 0.1c')
   os.system('gmtset MEASURE_UNIT inch')
   os.system('psbasemap -R'+str(left)+'/'+str(right)+'/'+str(bottom)+'/'+str(top)+' -JX10i/3i -Ba'+dx+'g'+dx+':"Time":/a'+dy+'g'+dy+':"'+y_text+'":WeSn -Gwhite -K > graph.ps')
   
   for imember in range(nmember):
      graph_file=open('graph.txt', 'w')
      for itime in range(ntime):
         x = itime
         graph_file.write(str(x)+' '+str(a[imember,itime])+'\n')
      graph_file.close()
      #os.system('psxy graph.txt -R -J -O -K -SC0.12 -Gblack >> graph.ps')
      os.system('psxy graph.txt -R -J -O -K -Wfat,black >> graph.ps')
   # Mean
   for iexp in range(nexp):
      graph_file=open('graph.txt', 'w')
      for itime in range(ntime):
         x = itime
         graph_file.write(str(x)+' '+str(e[iexp,itime])+'\n')
      graph_file.close()
      os.system('psxy graph.txt -R -J -O -K -S'+symbol[iexp]+'0.12 -G'+color[iexp]+' >> graph.ps')
      os.system('psxy graph.txt -R -J -O -K -Wthickest,'+color[iexp]+' >> graph.ps')

   lgd_file = open('graph.lgd','w')
   lgd_file.write('C black\n')
   #lgd_file.write('S 0.13i C 0.13i darkgreen 0.3p 0.25i Members\n')
   for iexp in range(nexp):
      lgd_file.write('S 0.13i '+symbol[iexp]+' 0.13i '+color[iexp]+' 0.3p 0.25i @~e@~='+str(epsilon[iexp])+', @~t@~='+str(tau)+'\n')
   #lgd_file.write('S 0.13i A 0.13i brown 0.3p 0.25i Best member\n')
   lgd_file.close()
   os.system('pslegend -R -J -O -K -C0.05i/0.05i -D'+str(left)+'/'+str(top)+'/2.5i/1.0i/LT graph.lgd >> graph.ps')
   
   os.system('echo "'+str(right)+' '+str(top)+' 36 0 4 RT ('+figtype+')" | pstext -R -J -O -K -Dj0.1/0.1 -Gblack -N >> graph.ps')
   if cumulative: os.system('echo "'+str(0.5*(left+right))+' '+str(top)+' 20 0 4 CB Cumulative rainfall" | pstext -R -J -O -K -Dj0./0.1 -Gblack -N >> graph.ps')
   else: os.system('echo "'+str(0.5*(left+right))+' '+str(top)+' 20 0 4 CB Regularized Gaussian-Hellinger Barycenter" | pstext -R -J -O -K -Dj0./0.1 -Gblack -N >> graph.ps')
   #os.system('echo "'+str(0.5*(left+right))+' '+str(top)+' 16 0 4 CB '+description+'" | pstext -R -J -O -K -Dj0./0.1 -Gblack -N >> graph.ps')
   #os.system('echo "'+str(left)+' '+str(bottom)+' 18 0 4 LT Initial date: '+initial_date+'" | pstext -R -J -O -K -Dj0./0.3 -Gblack -N >> graph.ps')
   #os.system('echo "'+str(right)+' '+str(bottom)+' 18 0 4 RT Forecast range: '+'%2.2d'%(time[-1]/60)+'h" | pstext -R -J -O -K -Dj0./0.3 -Gblack -N >> graph.ps')
   
   os.system('ps2raster -Tf graph.ps')
   os.system('convert -background white -flatten -trim +repage +rotate 90 '+density+' graph.pdf graph.'+format)

def write_netcdf(filename, c, e, b, o):
   os.system('rm -f '+filename)
   nc = netCDF4.Dataset(filename, 'w', format='NETCDF3_CLASSIC')
   nc.Conventions = 'COARDS'
   nc.createDimension('time', ntime)
   var = nc.createVariable('c', 'f4', ('time'))
   var[...] = c.astype(float32)
   var = nc.createVariable('e', 'f4', ('time'))
   var[...] = e.astype(float32)
   var = nc.createVariable('b', 'f4', ('time'))
   var[...] = b.astype(float32)
   var = nc.createVariable('o', 'f4', ('time'))
   var[...] = o.astype(float32)

# Sinkhorn log: stable
sigma0 = 5; mean0 = 30; factor0 = 1.
#sigma1 = 5; mean1 = 70; factor1 = 1.
sigma1 = 8; mean1 = 60; factor1 = 1.3
nmember = 2; ntime = 101
format = 'png'; density = '-density 100'
symbol = ['C','S','T','A']
color = ['blue','red','darkgreen','darkbrown']

epsilon = [0.0001, 0.001, 0.01, 0.1] # enstropic regulation (does not work if < 1e-4)
nexp = len(epsilon)
tau = 10.# mass regulation
q = 2.
amin = 1.e-30 # minimum value of a to avoid nan
stopping = 1.e-5
niteration = 10000

# Directories
root_dir = '/home/leduc/project/ot4ens'
input_dir = root_dir
output_dir = input_dir
# Work dir
work_dir = root_dir+'/../tmp'
current_dir = os.getcwd()
os.chdir(work_dir)
os.system('rm -rf *')

w = ones((nmember))/nmember
w[0] = 0.5; w[1] = 1.-w[0]
logw = zeros((nmember,ntime))
for imember in range(nmember): logw[imember] = log(w[imember])
tmp = zeros((ntime))
b = ones((2,ntime)) # barycenter
e = zeros((nexp,ntime)) # ensemble mean
a = zeros((nmember,ntime))
time = arange(ntime)
a[0,:] = factor0*exp(-0.5*(time-mean0)**2/sigma0**2)/(sqrt(2*pi)*sigma0)
a[1,:] = factor1*exp(-0.5*(time-mean1)**2/sigma1**2)/(sqrt(2*pi)*sigma1)

# Cost
C = zeros((ntime,ntime))
D = zeros((ntime,ntime))
for i in range(ntime):
   for j in range(i+1,ntime):
      D[i,j] = abs((i-j)/(ntime-1))**q
      D[j,i] = D[i,j]
pass

# Sinkhorn
a += amin
loga = log(a)
f = zeros((nmember,ntime))
g = zeros((nmember,ntime))
for iexp in range(nexp):
   b[:,:] = 1.; f[:,:] = 0.; g[:,:] = 0.
   phi = tau/(tau+epsilon[iexp])
   C[:,:] = D/epsilon[iexp]
   for it in range(niteration):
      # Kv or f
      for imember in range(nmember):
         f[imember] = phi*(loga[imember]-logsumexp(-C+g[imember,None,:],axis=1))
      # Ku or g
      for imember in range(nmember):
         g[imember] = logsumexp(-C+f[imember,:,None],axis=0)
      b[0] = b[1]
      b[1] = 1/(1-phi)*logsumexp((1-phi)*g+logw,axis=0)
      if it > 0:
         ratio = abs(b[1]-b[0]).max()/min(abs(b[1]).max(),abs(b[0]).max())
         print('barycenter:', it, ratio, exp(b[1].min()), exp(b[1].max()))
      if it > 0 and ratio < stopping: break
      for imember in range(nmember): g[imember] = phi*(b[1]-g[imember])
   e[iexp,:] = exp(b[1])
e -= amin
invalid = e < 0.; e[invalid] = 0.
#write_netcdf(output_file, c, e, b[1], o)
output_file = output_dir+'/Figure04a.'+format
plot_histogram(time, 500*a, 500*e, cumulative=False, figtype='a')
shutil.move('graph.'+format, output_file)
