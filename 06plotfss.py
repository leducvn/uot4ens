#!/usr/bin/env python3
import sys, os, datetime, re, shutil, netCDF4
from numpy import *

missing_value = -9999.

def plot_fss(experiment, threshold, fss):
   nexp, nthreshold = fss.shape
   bottom = 0.; top = 1.
   left = 0.; right = nthreshold
   dx = '1'; dy = '0.2'
   
   os.system('gmtset PAPER_MEDIA a3')
   os.system('gmtset X_ORIGIN 3c')
   os.system('gmtset Y_ORIGIN 3c')
   os.system('gmtset BASEMAP_TYPE PLAIN')
   os.system('gmtset TICK_LENGTH 0')
   os.system('gmtset ANNOT_FONT_SIZE 18p')
   os.system('gmtset LABEL_FONT_SIZE 20p')
   os.system('gmtset LABEL_OFFSET 0.1c')
   os.system('gmtset MEASURE_UNIT inch')
   os.system('psbasemap -R'+str(left)+'/'+str(right)+'/'+str(bottom)+'/'+str(top)+' -JX7i/5i -Bg'+\
              dx+':"":/a'+dy+'g'+dy+':"FSS":WeSn -Gwhite -K > graph.ps')
   for i in range(nthreshold):
      if threshold[i] < 1.: os.system('echo "'+str(i+0.5)+' '+str(bottom)+' 18 0 4 CT '+'%3.1f'%float(threshold[i])+'" | pstext -R -J -O -K -Dj0./0.1 -Gblack -N >> graph.ps')
      else: os.system('echo "'+str(i+0.5)+' '+str(bottom)+' 20 0 4 CT '+'%2.2d'%int(threshold[i])+'" | pstext -R -J -O -K -Dj0./0.1 -Gblack -N >> graph.ps')
   os.system('echo "'+str(0.5*(left+right))+' '+str(bottom)+' 20 0 4 CT Threshold (mm)" | pstext -R -J -O -K -Dj0./0.3 -Gblack -N >> graph.ps')
   
   for j in range(nexp):
      graph_file=open('graph.txt', 'w')
      for i in range(nthreshold):
         graph_file.write(str(i+0.5)+' '+str(fss[j,i])+'\n')
      graph_file.close()
      os.system('psxy graph.txt -R -J -O -K -S'+symbol[j]+'0.1 -G'+color[j]+' >> graph.ps')
      os.system('psxy graph.txt -R -J -O -K -Wthickest,'+color[j]+' >> graph.ps')
   pass
   
   lgd_file = open('graph.lgd','w')
   lgd_file.write('C blue\n')
   for j in range(nexp):
      lgd_file.write('S 0.1i '+symbol[j]+' 0.2i '+color[j]+' 0.20p 0.3i '+experiment[j]+'\n')
   lgd_file.close()
   os.system('pslegend -R -J -O -K -C0.05i/0.05i -D'+str(right)+'/'+str(top)+'/1.8i/2.0i/RT graph.lgd >> graph.ps')
   
   #os.system('echo "'+str(right)+' '+str(bottom)+' 36 0 4 RB ('+figtype+')" | pstext -R -J -O -K -Dj0.1/0.1 -Gblack -N >> graph.ps')
   os.system('echo "'+str(left)+' '+str(bottom)+' 18 0 9 LB '+variable_desc+'" | pstext -R -J -O -K -Dj0.1/0.3 -Gblue -N >> graph.ps')
   os.system('echo "'+str(left)+' '+str(bottom)+' 18 0 9 LB Spatial scale: 15 km" | pstext -R -J -O -K -Dj0.1/0.1 -Gblue -N >> graph.ps')
   os.system('ps2raster -Tf graph.ps')
   os.system('convert -background white -flatten -trim +repage +rotate 90 '+density+' graph.pdf graph.'+format)

format = 'png'
density = '-density 100'
figtype = 'b'
experiment = 'Kyushu'; iexp = 100
name = ['c', 'b', 'e']; fcstname = ['MEPS DET', 'MEPS UOT', 'MEPS ENS']
start_date = '202007030600'; end_date = '202007030600'
interval = 3*60; start = 9*60; end = 18*60
color = ['red', 'blue', 'orange', 'green', 'magenta', 'cyan', 'darkyellow']
symbol = ['C', 'S', 'T', 'A', 'T', 'A', 'I', 'D', 'N', 'H']
nexp = len(name)

# Directories
root_dir = '/home/leduc/project/uot4ens'
input_dir = root_dir
input_file = []
for i in range(nexp):
   data_file = input_dir+'/'+experiment+'%2.2d'%iexp+'fss'+name[i]+'vs'+name[0]+'_'+'%2.2d'%(start/60)+'%2.2d'%(end/60)+'.nc'
   input_file.append(data_file)
   if i > 0 and not os.path.isfile(input_file[i]): sys.exit(0)
output_dir = input_dir
output_file = output_dir+'/Figure06.'+format
#if os.path.isfile(output_file): sys.exit(0)
work_dir = root_dir+'/../tmp'
current_dir = os.getcwd()
os.chdir(work_dir)

# Inventory
nc = netCDF4.Dataset(input_file[1])
nc.set_auto_mask(False)
timescale = nc.variables['timescale'][:]
scale = nc.variables['scale'][:]
threshold = nc.variables['threshold'][:]
ntimescale = len(timescale)
nscale = len(scale)
nthreshold = len(threshold)
nc.close()

# Read
fss = zeros((nexp,nthreshold))
for iexp in range(1,nexp):
   nc = netCDF4.Dataset(input_file[iexp])
   nc.set_auto_mask(False)
   if iexp == 1: fss[:2] = nc.variables['fss'][:][:,0,0,:]
   else: fss[iexp] = nc.variables['fss'][:][1,0,0]
   nc.close()
print(fss)

# Plot fss
print('Plot file '+output_file)
if interval < 60: variable = 'R'+'%2.2d'%interval+'M  '
else: variable = 'R'+'%2.2d'%(interval/60)+'H  '
variable_desc = variable+'%2.2d'%(start/60)+'-'+'%2.2d'%(end/60)+'h'
plot_fss(fcstname, threshold, fss)
shutil.move('graph.'+format, output_file)
os.chdir(current_dir)
#shutil.rmtree(work_dir)
