import netCDF4
from matplotlib import pylab

nc = netCDF4.Dataset('coords_CF_ORCA12_GO6-2.nc')
# read the latitudes with various staggerings

latt = nc.variables['latt'][0:10, 0:10]
latu = nc.variables['latu'][0:10, 0:10]
latv = nc.variables['latv'][0:10, 0:10]
latw = nc.variables['latw'][0:10, 0:10]

lont = nc.variables['lont'][0:10, 0:10]
lonu = nc.variables['lonu'][0:10, 0:10]
lonv = nc.variables['lonv'][0:10, 0:10]
lonw = nc.variables['lonw'][0:10, 0:10]

pylab.plot(lonw, latw, 'mo')
pylab.plot(lont, latt, 'bx')
pylab.plot(lonu, latu, 'b>')
pylab.plot(lonv, latv, 'g^')


pylab.show()
