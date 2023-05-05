#!/usr/bin/env python3
import sys, os, datetime, re, shutil, netCDF4
import scipy.stats
from numpy import *

def plot_histogram(hour, rain, cumulative=False, figtype=''):
   nmember, ntime = rain.shape
   nmember -= 2
   width = 0.5 # box width
   left = 3.; right = ntime-1
   #bottom = 1000.; top = 0.
   ensrain = rain[1:-2,:]
   #for itime in range(ntime):
      #bottom = min(bottom,scipy.stats.scoreatpercentile(ensrain[:,itime],2.5))
      #top = max(top,scipy.stats.scoreatpercentile(ensrain[:,itime],97.5))
   #bottom = min(bottom,min(rain[0]),min(rain[-1])); top = max(top,max(rain[0]),max(rain[-1]))
   #print(bottom,top)
   #bottom -= 1; top += 1
   bottom = 0.; top = rain.max()+0.1
   dx = '1'; dy = '10'; y_text = 'Rainrate (mm/hour)'
   if cumulative: dy = '50'; y_text = 'Accumulated Rainfall (mm)'
   
   os.system('gmtset PAPER_MEDIA a3')
   os.system('gmtset X_ORIGIN 3c')
   os.system('gmtset Y_ORIGIN 3c')
   os.system('gmtset BASEMAP_TYPE PLAIN')
   os.system('gmtset TICK_LENGTH 0')
   os.system('gmtset ANNOT_FONT_SIZE 14p')
   os.system('gmtset LABEL_FONT_SIZE 16p')
   os.system('gmtset LABEL_OFFSET 0.1c')
   os.system('gmtset MEASURE_UNIT inch')
   os.system('psbasemap -R'+str(left)+'/'+str(right)+'/'+str(bottom)+'/'+str(top)+' -JX10i/3i -Bg'+dx+'::/a'+dy+'g'+dy+':"'+y_text+'":WeSn -Gwhite -K > graph.ps')
   for itime in range(3,ntime):
      os.system('echo "'+str(itime)+' '+str(bottom)+' 16 0 4 CT '+hour[itime]+'" | pstext -R -J -O -K -Dj0./0.1 -Gblack -N >> graph.ps')
   os.system('echo "'+str(0.5*(right+left))+' '+str(bottom)+' 20 0 4 CT Japanese Time (hour)" | pstext -R -J -O -K -Dj0./0.3 -Gblack -N >> graph.ps')
   
   # Members
   for imember in range(1,nmember):
      graph_file=open('graph.txt', 'w')
      for itime in range(ntime):
         x = itime
         graph_file.write(str(x)+' '+str(rain[imember,itime])+'\n')
      graph_file.close()
      #os.system('psxy graph.txt -R -J -O -K -SC0.12 -Ggray >> graph.ps')
      os.system('psxy graph.txt -R -J -O -K -Wthin,darkgray >> graph.ps')
   # Barycenter
   graph_file=open('graph.txt', 'w')
   for itime in range(ntime):
      x = itime
      graph_file.write(str(x)+' '+str(rain[0,itime])+'\n')
   graph_file.close()
   os.system('psxy graph.txt -R -J -O -K -SS0.12 -Gred >> graph.ps')
   os.system('psxy graph.txt -R -J -O -K -Wthicker,red >> graph.ps')
   # Mean
   graph_file=open('graph.txt', 'w')
   for itime in range(ntime):
      x = itime
      graph_file.write(str(x)+' '+str(rain[-1,itime])+'\n')
   graph_file.close()
   os.system('psxy graph.txt -R -J -O -K -ST0.12 -Gblue >> graph.ps')
   os.system('psxy graph.txt -R -J -O -K -Wthicker,blue >> graph.ps')
   # Obs
   graph_file=open('graph.txt', 'w')
   for itime in range(ntime):
      x = itime
      graph_file.write(str(x)+' '+str(rain[-2,itime])+'\n')
   graph_file.close()
   os.system('psxy graph.txt -R -J -O -K -SA0.12 -Gdarkgreen >> graph.ps')
   os.system('psxy graph.txt -R -J -O -K -Wthicker,darkgreen >> graph.ps')
   
   lgd_file = open('graph.lgd','w')
   lgd_file.write('C black\n')
   lgd_file.write('S 0.13i A 0.13i darkgreen 0.3p 0.25i Observation\n')
   lgd_file.write('S 0.13i T 0.13i blue 0.3p 0.25i Ensemble mean\n')
   lgd_file.write('S 0.13i S 0.13i red 0.3p 0.25i GH Barycenter\n')
   #lgd_file.write('S 0.13i A 0.13i brown 0.3p 0.25i Best member\n')
   lgd_file.close()
   os.system('pslegend -R -J -O -K -C0.05i/0.05i -D'+str(left)+'/'+str(top)+'/2.1i/1.0i/LT graph.lgd >> graph.ps')
   
   os.system('echo "'+str(right)+' '+str(bottom)+' 36 0 4 RB ('+figtype+')" | pstext -R -J -O -K -Dj0.1/0.1 -Gblack -N >> graph.ps')
   if cumulative: os.system('echo "'+str(0.5*(left+right))+' '+str(top)+' 20 0 4 CB Accumulated rainfall over '+area+'" | pstext -R -J -O -K -Dj0./0.4 -Gblack -N >> graph.ps')
   else: os.system('echo "'+str(0.5*(left+right))+' '+str(top)+' 20 0 4 CB Time series of rainfall over '+area+'" | pstext -R -J -O -K -Dj0./0.4 -Gblack -N >> graph.ps')
   os.system('echo "'+str(0.5*(left+right))+' '+str(top)+' 16 0 4 CB '+description+'" | pstext -R -J -O -K -Dj0./0.1 -Gblack -N >> graph.ps')
   #os.system('echo "'+str(left)+' '+str(bottom)+' 18 0 4 LT Initial date: '+initial_date+'" | pstext -R -J -O -K -Dj0./0.3 -Gblack -N >> graph.ps')
   #os.system('echo "'+str(right)+' '+str(bottom)+' 18 0 4 RT Forecast range: '+'%2.2d'%(time[-1]/60)+'h" | pstext -R -J -O -K -Dj0./0.3 -Gblack -N >> graph.ps')
   
   os.system('ps2raster -Tf graph.ps')
   os.system('convert -background white -flatten -trim +repage +rotate 90 '+density+' graph.pdf graph.'+format)

iexp = 100
analysis_date = '202007030600'
#experiment = 'Kyushu02km'; nmember = 101; start = 1*60; end = 24*60; period = 1*60
experiment = 'Fugaku05km'; nmember = 21; start = 1*60; end = 27*60; period = 1*60
fortran_file = 'fort'
format = 'png'; density = '-density 100'
#description = 'Experiment: LETKF100; Ensemble members: '+str(nmember-1)
description = 'Experiment: MEPS; Ensemble members: '+str(nmember-1)
area = 'Ichifusa'

# Directories
root_dir = '/home/leduc/project/ot4ens'
input_dir = root_dir
input_file1 = input_dir+'/'+area+'%2.2d'%iexp+'.nc'
input_file2 = input_dir+'/'+area+'%2.2d'%iexp+'mean.nc'
output_dir = input_dir
# Work dir
work_dir = root_dir+'/../tmp'
current_dir = os.getcwd()
os.chdir(work_dir)
os.system('rm -rf *')

# Dates
time = list(range(start-period,end+1,period))
ntime = len(time)
hour = []
date = datetime.datetime(int(analysis_date[:4]), int(analysis_date[4:6]), int(analysis_date[6:8]), int(analysis_date[8:10]), int(analysis_date[10:12]))
date += datetime.timedelta(hours=9)
initial_date = date.strftime('%Y')+'/'+date.strftime('%m')+'/'+date.strftime('%d')+' - '+date.strftime('%H')+':'+date.strftime('%M')+' JST'
for itime in range(ntime):
   valid_date = date + datetime.timedelta(minutes=int(time[itime]))
   if valid_date.strftime('%H') == '00': hour.append(valid_date.strftime('%b')+' '+valid_date.strftime('%d'))
   else: hour.append(valid_date.strftime('%H'))
print(time)

# Member
member = list(range(nmember))
for imember in range(nmember): member[imember] = '%4.4d'%imember
member.append('obs')
member.append('ens')

# Read
histogram = zeros((nmember+2,ntime))
cumulative = zeros((nmember+2,ntime))
nc = netCDF4.Dataset(input_file1, 'r')
histogram[:nmember+1,:] = nc.variables['rain'][:]
nc.close()
nc = netCDF4.Dataset(input_file2, 'r')
histogram[nmember+1] = nc.variables['e'][:]
histogram[0] = nc.variables['b'][:]
nc.close()
for imember in range(nmember+2):
   cumulative[imember] = histogram[imember]
   for itime in range(1,ntime): cumulative[imember,itime] += cumulative[imember,itime-1]
#output_file = output_dir+'/Figure06a.'+format
#plot_histogram(hour, histogram, cumulative=False, figtype='a')
#shutil.move('graph.'+format, output_file)
output_file = output_dir+'/Figure06a.'+format
plot_histogram(hour, cumulative, cumulative=True, figtype='a')
shutil.move('graph.'+format, output_file)

