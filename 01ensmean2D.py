#!/usr/bin/env python3
import sys, os, shutil, datetime, re, netCDF4
from mri_mapproj import mapproj
from nhmlib import *
from numpy import *
from mpi4py import MPI

missing_value = -9999.
capitals={'Tokyo'        : ('139.581', '35.615', 'CT'),
          'Taipei'       : ('121.63', '25.07', 'CT'),
          'Seoul'        : ('127.00', '37.50', 'CT'),
          'Beijing'      : ('116.43', '39.92', 'CT')}
cities={'Osaka'          : ('135.502', '34.694', 'CT'),
        'Nagoya'         : ('136.906', '35.181', 'CT'),
        'Sapporo'        : ('141.347', '43.064', 'CT'),
        'Fukuoka'        : ('130.418', '33.606', 'CT'),}

def write_field(filename, lon, lat, field):
   os.system('rm -f '+filename)
   nc = netCDF4.Dataset(filename, 'w',format='NETCDF3_CLASSIC')
   nlon = len(lon)
   nlat = len(lat)
   nc.createDimension('lon', nlon)
   nc.createDimension('lat', nlat)
   var = nc.createVariable('lon', 'f4', ('lon'))
   var[:] = lon.astype(float32)
   var.units = 'degreesE'
   var = nc.createVariable('lat', 'f4', ('lat'))
   var[:] = lat.astype(float32)
   var.units = 'degreesN'
   
   var = nc.createVariable('field', 'f4', ('lat','lon'))
   var[:,:] = field.astype(float32)
   var.missing_value = missing_value
   nc.close()

def plot_mer(area):
   factor = 1.1
   slon = str(first_lon)
   slat = str(first_lat)
   elon = str(last_lon)
   elat = str(last_lat)
   dlon = float(elon) - float(slon)
   dlat = float(elat) - float(slat)
   if dlat >= dlon:
      xscale = factor*6.4/dlon
      yscale = factor*9.8/dlat
      portrait = True
   else:
      xscale = factor*6.4/dlat
      yscale = factor*9.8/dlon
      portrait = False
   scale = min(xscale, yscale)
   region = ' -R'+str(slon)+'/'+str(elon)+'/'+str(slat)+'/'+str(elat)
   projection = ' -Jm'+str(scale)+'i'
   if area == 'whole': boundary = ' -Ba2g2/a2g2wesn'
   else: boundary = ' -Ba1g1/a1g1wesn'
   
   # Basemap
   os.system('gmtset PAPER_MEDIA a3')
   os.system('gmtset X_ORIGIN 1c')
   os.system('gmtset Y_ORIGIN 2c')
   os.system('gmtset BASEMAP_TYPE PLAIN')
   os.system('gmtset TICK_LENGTH 0')
   os.system('gmtset MEASURE_UNIT inch')
   if portrait: os.system('psbasemap'+region+projection+boundary+' -P -K > map.ps')
   else: os.system('psbasemap'+region+projection+boundary+' -K > map.ps')
   #os.system('pscoast -R -J -B -Dl -Glightbrown -Slightblue -N1/faint/black -Wfaint,black -O -K >> map.ps')
   
   # Color bar
   if period <= 6*60:
      #hue={0.:'255 255 255', 1.:'184 184 255', 5.:'127 153 229', 10.:'127 180 210', 15.:'127 210 180', 20.:'200 255 150', \
           #25.:'150 255 120', 30.:'100 250 100', 35.:'50 240 50', 40.:'255 255 180', 45.:'255 255 130', 50.:'255 255 80', \
           #55.:'255 240 50', 60.:'255 210 150', 65.:'255 180 120', 70.:'255 150 90', 75.:'255 100 50', 80.:'255 75 40', \
           #85.:'255 50 30', 90.:'255 25 19', 100.:''}
      hue={0.:'255 255 255', 0.4:'135 206 255', 1.:'0 255 255', 5.:'0 210 0', 10.:'255 255 0', 20.:'255 192 0', \
           50.:'255 0 0', 100.:''}
   elif period <= 12*60:
      #hue={0.:'255 255 255', 5.:'184 184 255', 20.:'127 153 229', 30.:'127 180 210', 40.:'127 210 180', 50.:'200 255 150', \
           #60.:'150 255 120', 70.:'100 250 100', 80.:'50 240 50', 90.:'255 255 180', 100.:'255 255 130', 110.:'255 255 80', \
           #120.:'255 240 50', 130.:'255 210 150', 140.:'255 180 120', 150.:'255 150 90', 160.:'255 100 50', 170.:'255 75 40', \
           #180.:'255 50 30', 190.:'255 25 19', 200.:''}
      hue={0.:'255 255 255', 1.:'135 206 255', 5.:'0 255 255', 10.:'0 210 0', 20.:'255 255 0', 50.:'255 192 0', \
           100.:'255 0 0', 200.:''}
   else:
      hue={0.:'255 255 255', 10.:'184 184 255', 40.:'127 153 229', 60.:'127 180 210', 80.:'127 210 180', 100.:'200 255 150', \
           120.:'150 255 120', 140.:'100 250 100', 160.:'50 240 50', 180.:'255 255 180', 200.:'255 255 130', 220.:'255 255 80', \
           240.:'255 240 50', 260.:'255 210 150', 280.:'255 180 120', 300.:'255 150 90', 320.:'255 100 50', 340.:'255 75 40', \
           360.:'255 50 30', 380.:'255 25 19', 400.:''}
   cpt_file = open('map.cpt','w')
   cpt_file.write('#COLOR_MODEL = +RGB\n')
   #threshold = hue.keys(); threshold.sort()
   threshold = sorted(hue)
   for i in range(len(threshold)-1):
      from_value = str(threshold[i])
      to_value = str(threshold[i+1])
      cpt_file.write(from_value+' '+hue[threshold[i]]+' '+to_value+' '+hue[threshold[i]]+'\n')
   cpt_file.write('B 255 255 255\n')
   if period <= 12*60: cpt_file.write('F 255 192 255\n')
   else: cpt_file.write('F 250 0 0\n')
   cpt_file.write('N -\n')
   cpt_file.close()
   f = os.popen('echo "'+str(last_lon)+' '+str(last_lat)+'" | mapproject -R -J')
   xy = re.split('\s', f.readline())
   xpos = float(xy[0]) + 0.1
   ypos = 0.5*float(xy[1]) - 0.15
   length = float(xy[1]) - 0.3
   bar_scale = ' -D'+str(xpos)+'i/'+str(ypos)+'i/'+str(length)+'i/0.2i'
   os.system('psscale -Cmap.cpt -O -K'+bar_scale+' -B:"":/:mm: -L -Ef0.3i >> map.ps')
   
   # Rain
   if area != 'whole': os.system('grdsample rain.nc -I0.01 -Q -Grain.nc')
   os.system('grdsample rain.nc -I0.01 -Q -Grain.nc')
   os.system('grdimage rain.nc -Cmap.cpt -R -J -O -K -Q >> map.ps')
   
   # Cities
   for capital in capitals.keys():
      lon, lat, justify = capitals[capital]
      os.system('echo "'+lon+' '+lat+' 13 0 4 '+justify+' '+capital+'" | pstext -R -J -O -K -Dj0.1/0.1 -Gred >> map.ps')
      os.system('echo "'+lon+' '+lat+'" | psxy -R -J -O -K -Sa0.12 -Gred >> map.ps')
   #for city in cities.keys():
      #lon, lat, justify = cities[city]
      #os.system('echo "'+lon+' '+lat+' 13 0 4 '+justify+' '+city+'" | pstext -R -J -O -K -Dj0.1/0.1 -Gred >> map.ps')
      #os.system('echo "'+lon+' '+lat+'" | psxy -R -J -O -K -Sc0.07 -Gred >> map.ps')
   if area == 'whole': os.system('pscoast -R -J'+boundary+' -Di -N1/thin/brown -Wthin,brown -O -K >> map.ps')
   else: os.system('pscoast -R -J'+boundary+' -Dh -N1/thin/brown -Wthin,brown -O -K >> map.ps')
   
   f = os.popen('echo "'+str(last_lon)+' '+str(last_lat)+'" | mapproject -R -J')
   xy = re.split('\s', f.readline())
   f = os.popen('echo "0 '+xy[1]+'" | mapproject -R -J -I')
   ll = re.split('\s', f.readline())
   if period >= 60: desc = 'Obs R'+'%2.2d'%(period/60)+'h'
   else: 'Obs R'+'%2.2d'%period+'min'
   os.system('echo "'+ll[0]+' '+ll[1]+' 28 0 4 LT GSMaP R'+desc+'" | pstext -R -J -O -K -Dj0.1/0.1 -Gblue -N >> map.ps')
   f = os.popen('echo "'+xy[0]+' 0" | mapproject -R -J -I')
   ll = re.split('\s', f.readline())
   os.system('echo "'+ll[0]+' '+ll[1]+' 16 0 4 RB '+valid_date+'" | pstext -R -J -O -K -Dj0.1/0.1 -Gblue -N >> map.ps')
   os.system('ps2raster -Tf map.ps')
   os.system('convert -background white -flatten -trim +repage +rotate 90 map.pdf map.gif')

def plot_lmn(area, fcsttype, figtype):
   scale = 9
   region = ' -R'+str(slon0)+'/'+str(slat0)+'/'+str(elon0)+'/'+str(elat0)+'r'
   projection = ' -JL'+str(stdlon)+'/'+str(stdlat)+'/30/60/'+str(scale)+'i'
   if area == 'whole': boundary = ' -Ba2g2/a2g2wesn'
   else: boundary = ' -Ba1g1/a1g1wesn'
   
   # Basemap
   os.system('gmtset PAPER_MEDIA a3')
   os.system('gmtset X_ORIGIN 1c')
   os.system('gmtset Y_ORIGIN 2c')
   os.system('gmtset BASEMAP_TYPE PLAIN')
   os.system('gmtset TICK_LENGTH 0')
   os.system('gmtset ANNOT_FONT_SIZE 20p')
   os.system('gmtset LABEL_FONT_SIZE 20p')
   os.system('gmtset MEASURE_UNIT inch')
   os.system('psbasemap'+region+projection+boundary+' -K > map.ps')
   #os.system('pscoast -R -J -B -Dl -Glightbrown -Slightblue -N1/faint/black -Wfaint,black -O -K >> map.ps')
   
   # Color bar
   if period <= 6*60:
      #hue={0.:'255 255 255', 0.4:'135 206 255', 1.:'0 255 255', 5.:'0 210 0', 10.:'255 255 0', 20.:'255 192 0', \
           #50.:'255 0 0', 100.:''}
      hue={0.:'255 255 255', 1.:'119 170 221', 5.:'153 221 255', 10.:'68 187 153', 20.:'187 204 51', 50.:'238 221 136', \
           100.:'238 136 102', 200.:'255 170 187'}
   elif period <= 12*60:
      #hue={0.:'255 255 255', 1.:'135 206 255', 5.:'0 255 255', 10.:'0 210 0', 20.:'255 255 0', 50.:'255 192 0', \
           #100.:'255 0 0', 200.:''}
      hue={0.:'255 255 255', 1.:'119 170 221', 5.:'153 221 255', 10.:'68 187 153', 20.:'187 204 51', 50.:'238 221 136', \
           100.:'238 136 102', 200.:'255 170 187'}
   else:
      hue={0.:'255 255 255', 10.:'184 184 255', 40.:'127 153 229', 60.:'127 180 210', 80.:'127 210 180', 100.:'200 255 150', \
           120.:'150 255 120', 140.:'100 250 100', 160.:'50 240 50', 180.:'255 255 180', 200.:'255 255 130', 220.:'255 255 80', \
           240.:'255 240 50', 260.:'255 210 150', 280.:'255 180 120', 300.:'255 150 90', 320.:'255 100 50', 340.:'255 75 40', \
           360.:'255 50 30', 380.:'255 25 19', 400.:''}
   cpt_file = open('map.cpt','w')
   cpt_file.write('#COLOR_MODEL = +RGB\n')
   #threshold = hue.keys(); threshold.sort()
   threshold = sorted(hue)
   for i in range(len(threshold)-1):
      from_value = str(threshold[i])
      to_value = str(threshold[i+1])
      cpt_file.write(from_value+' '+hue[threshold[i]]+' '+to_value+' '+hue[threshold[i]]+'\n')
   cpt_file.write('B '+hue[threshold[0]]+'\n')
   cpt_file.write('F '+hue[threshold[-1]]+'\n')
   cpt_file.write('N -\n')
   cpt_file.close()
   f = os.popen('echo "'+str(elon0)+' '+str(elat0)+'" | mapproject -R -J')
   xy = re.split('\s', f.readline())
   xpos = float(xy[0]) + 0.1
   ypos = 0.5*float(xy[1]) - 0.15
   length = float(xy[1]) - 0.3
   bar_scale = ' -D'+str(xpos)+'i/'+str(ypos)+'i/'+str(length)+'i/0.2i'
   os.system('psscale -Cmap.cpt -O -K'+bar_scale+' -B:"":/:mm: -L -Ef0.3i >> map.ps')
   
   # Rain
   if area != 'whole': os.system('grdsample rain.nc -I0.01 -Q -Grain.nc')
   os.system('grdsample rain.nc -I0.01 -Q -Grain.nc')
   os.system('grdimage rain.nc -Cmap.cpt -R -J -O -K -Q >> map.ps')
   
   # Cities
   for capital in capitals.keys():
      lon, lat, justify = capitals[capital]
      os.system('echo "'+lon+' '+lat+' 13 0 4 '+justify+' '+capital+'" | pstext -R -J -O -K -Dj0.1/0.1 -Gred >> map.ps')
      os.system('echo "'+lon+' '+lat+'" | psxy -R -J -O -K -Sa0.12 -Gred >> map.ps')
   #for city in cities.keys():
      #lon, lat, justify = cities[city]
      #os.system('echo "'+lon+' '+lat+' 13 0 4 '+justify+' '+city+'" | pstext -R -J -O -K -Dj0.1/0.1 -Gred >> map.ps')
      #os.system('echo "'+lon+' '+lat+'" | psxy -R -J -O -K -Sc0.07 -Gred >> map.ps')
   if area == 'whole': os.system('pscoast -R -J'+boundary+' -Di -N1/thin/brown -Wthin,brown -O -K >> map.ps')
   else: os.system('pscoast -R -J'+boundary+' -Dh -N1/thin/brown -Wthin,brown -O -K >> map.ps')
   
   os.system('echo "'+str(slon0)+' '+str(slat0)+' 36 0 4 LB ('+figtype+')" | pstext -R -J -O -K -Dj0.1/0.1 -Gblue -N >> map.ps')
   f = os.popen('echo "'+str(elon0)+' '+str(elat0)+'" | mapproject -R -J')
   xy = re.split('\s', f.readline())
   f = os.popen('echo "0 '+xy[1]+'" | mapproject -R -J -I')
   ll = re.split('\s', f.readline())
   if period >= 60: desc = fcsttype+' '+'%2.2d'%(period/60)+'h'
   else: desc = fcsttype+' '+'%2.2d'%period+'min'
   os.system('echo "'+ll[0]+' '+ll[1]+' 36 0 4 LT '+desc+'" | pstext -R -J -O -K -Dj0.1/0.1 -Gblue -N >> map.ps')
   #f = os.popen('echo "'+xy[0]+' '+xy[1]+'" | mapproject -R -J -I')
   #ll = re.split('\s', f.readline())
   #os.system('echo "'+ll[0]+' '+ll[1]+' 28 0 9 RT '+str(fcsttime/60.)+'h Forecast" | pstext -R -J -O -K -Dj0.1/0.1 -Gblue -N >> map.ps')
   f = os.popen('echo "'+xy[0]+' 0" | mapproject -R -J -I')
   ll = re.split('\s', f.readline())
   os.system('echo "'+ll[0]+' '+ll[1]+' 28 0 9 RB Valid: '+valid_date+'" | pstext -R -J -O -K -Dj0.1/0.1 -Gblue -N >> map.ps')
   os.system('ps2raster -Tf map.ps')
   os.system('convert -background white -flatten -trim +repage +rotate 90 '+density+' map.pdf map.'+format)

area = 'Kyushu'; analysis_date = '202007030600'
#iexp = 2; nmember0 = 1001
iexp = 100; nmember0 = 21
#experiment = 'Kyushu02km'; start = 6*60; end = 15*60
experiment = 'Fugaku05km'; start = 9*60; end = 18*60

format = 'png'
density = '-density 100'
fortran_file = 'fort'
figtype = ['a', 'b', 'c']
fcsttype = ['OBS', 'DET', 'ENS']
period = end-start
if period >= 60: variable = 'rain'+'%2.2d'%(period/60)+'h'
else: variable = 'rain'+'%2.2d'%(period/60)+'m'

# Directories
root_dir = '/home/leduc/project/uot4ens'
nhm_dir = root_dir+'/../../nhm'
const_dir = nhm_dir+'/exp/'+experiment+'/const'
input_dir = root_dir
input_file = input_dir+'/'+area+'%2.2d'%iexp+'_'+'%2.2d'%(start/60)+'%2.2d'%(end/60)+'.nc'
output_dir = input_dir
# Work dir
comm = MPI.COMM_WORLD
nprocessor = comm.Get_size()
myid = comm.Get_rank()
work_dir = root_dir+'/../tmp'
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
tmp = zeros((ny,nx))
c = zeros((ny,nx)) # control
o = zeros((ny,nx)) # obs
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

# Inventory
nc = netCDF4.Dataset(const_dir+'/parameters.nc','r')
nx0 = nc.nx; ny0 = nc.ny
projection = nc.projection
resolution = nc.resolution
xi = nc.xi; xj = nc.xj
xlon = nc.xlon; xlat = nc.xlat
stdlon = nc.slon; stdlat = nc.slat
nc.close()

# Domain
slon0 = 127.8; elon0 = 133.8
slat0 = 29.6; elat0 = 34.6
if myid == 0: print('Domain: ', slon0, elon0, slat0, elat0)
lon2D = zeros((1,1)); lat2D = zeros((1,1))
lon2D[0,0] = slon0; lat2D[0,0] = slat0
x2D, y2D = mapproj.latlon2xy(projection, 1000*resolution, 1000*resolution, xi, xj, xlon, xlat, stdlon, stdlat, lon2D, lat2D)
obsx = x2D[0,0] - 1; imin = int(round(obsx))
obsy = ny0+1 - y2D[0,0] - 1; jmin = int(round(obsy))
lon2D[0,0] = elon0; lat2D[0,0] = elat0
x2D, y2D = mapproj.latlon2xy(projection, 1000*resolution, 1000*resolution, xi, xj, xlon, xlat, stdlon, stdlat, lon2D, lat2D)
obsx = x2D[0,0] - 1; imax = int(round(obsx))
obsy = ny0+1 - y2D[0,0] - 1; jmax = int(round(obsy))

# Latlon grid
slon = lon.min(); elon = lon.max()
slat = lat.min(); elat = lat.max()
dlon = dlat = (lon[0,-1]-lon[0,0])/(nx-1)
nlon = int(floor((elon-slon)/dlon)) + 1
nlat = int(floor((elat-slat)/dlat)) + 1
elon = slon + (nlon-1)*dlon
elat = slat + (nlat-1)*dlat
lon = linspace(slon,elon,nlon)
lat = linspace(slat,elat,nlat)
if myid == 0:
   print('Latlon grid:')
   print(nlon, nlat)
   print(slon, elon, slat, elat, dlon, dlat)
lon2D = zeros((nlon,nlat))
lat2D = zeros((nlon,nlat))
for i in range(nlon): lat2D[i,:] = lat[:]
for j in range(nlat): lon2D[:,j] = lon[:]
xlon, ylat = mapproj.latlon2xy(projection, 1000*resolution, 1000*resolution, xi, xj, xlon, xlat, stdlon, stdlat, lon2D, lat2D)
#print xlon.min(), xlon.max(), ylat.min(), ylat.max()
ylat = ny0+1 - ylat
x = linspace(1,nx0,nx0)
y = linspace(1,ny0,ny0)

# Interpolate
tmp = zeros((ny0,nx0))
valid = ones((ny0,nx0),dtype=int32)
date = datetime.datetime(int(analysis_date[:4]), int(analysis_date[4:6]), int(analysis_date[6:8]), int(analysis_date[8:10]), 0)
valid_date = date + datetime.timedelta(minutes=(end+9*60))
minute = valid_date.strftime('%M')
hour   = valid_date.strftime('%H')
day    = valid_date.strftime('%d')
month  = valid_date.strftime('%m')
year   = valid_date.strftime('%Y')
valid_date = year+'/'+month+'/'+day+' '+hour+':'+minute+' JST'
for k in range(3):
   if myid != 0: continue
   if k == 0: tmp[jmin:jmax+1,imin:imax+1] = o[:,:]
   elif k == 1: tmp[jmin:jmax+1,imin:imax+1] = c[:,:]
   else: tmp[jmin:jmax+1,imin:imax+1] = e[:,:]
   invalid, rain = mapproj.nearest_itpl(x, y, valid, tmp, xlon.T, ylat.T)
   invalid = invalid == 1
   rain[invalid] = 0.
   write_field('rain.nc', lon, lat, rain)
   output_file = output_dir+'/Figure01'+figtype[k]+'.'+format
   if projection == 'MER ': plot_mer(area)
   elif projection == 'LMN ': plot_lmn(area, fcsttype[k], figtype[k])
   shutil.move('map.'+format, output_file)
for i in range(nprocessor):
   if i%nprocessor != myid: continue
   os.chdir(current_dir)
   #os.system('rm -rf '+work_dir)
