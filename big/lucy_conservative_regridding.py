import iris
import iris.quickplot as qplt
import iris.plot as iplt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mplcolors
import cf_units
from iris.experimental.regrid_conservative import regrid_conservative_via_esmpy
import datetime
import os

def make_cube(grid_res=0.1,datavalues='zero'):
    """Build a cube at required resolution, filled with zeros"""
    grid_values = np.arange(0.,10.+grid_res,grid_res)

    lat_lon_coord_system = iris.coord_systems.GeogCS(semi_major_axis =
                                               iris.fileformats.pp.EARTH_RADIUS)
    lon_coord = iris.coords.DimCoord(grid_values,
                                     standard_name='longitude',
                                     units='degrees',
                                     coord_system=lat_lon_coord_system)
    lat_coord = iris.coords.DimCoord(grid_values,
                                     standard_name='latitude',
                                     units='degrees',
                                     coord_system=lat_lon_coord_system)

    if isinstance(datavalues,float):
        print 'here'
        data = np.zeros((len(grid_values), len(grid_values)))
        data[:] = datavalues
    elif datavalues == 'zero':
        data = np.zeros((len(grid_values), len(grid_values)))
    elif datavalues == 'random':
        data = np.random.random((len(grid_values), len(grid_values)))
    else:
        raise IOError('Unknown value for datavalues')

    cube = iris.cube.Cube(data,long_name='Emissions',units='ug/m2/s')
    cube.add_dim_coord(lat_coord,0)
    cube.add_dim_coord(lon_coord,1)
    cube.coord('latitude').guess_bounds()
    cube.coord('longitude').guess_bounds()
    return cube

def get_cube_from_file(filename):
    cube = iris.load_cube(filename)
    coord_names = [ coord.name() for coord in cube.coords()]

    if 'longitude' in coord_names:
        cube.coord('latitude').guess_bounds()
        cube.coord('longitude').guess_bounds()
    else:
        cube.coord('grid_latitude').guess_bounds()
        cube.coord('grid_longitude').guess_bounds()

    if 'time' in coord_names:    
        cube.remove_coord('time')
        
    if cube.data.max() < 1.e-6:
        cube.units = 'kg/m2/s'
    else:    
        cube.units = 'ug/m2/s'

    return cube

def emission_shapes(data):
    """Add some random shapes to represent emissions"""
    data[10:20,10:20] = 3.
    data[30,30] = 5.
    data[39:42,39:42] = 8.
    data[50:60,50] = 9.
    data[73:74,63:69] = 7.
    data[80:95,83:99] = 2.
    data[86:89,91:94] = 1.
    data[10:20:2,81:91:2] = 6.
    data[73:91,11:17:2] = 4.
    return data

def emission_total(cube):
    """Calculate total emission in kg/s"""

    grid_areas = iris.analysis.cartography.area_weights(cube)
    coord_names = [ coord.name() for coord in cube.coords()]

    if 'longitude' in coord_names:
        total_cube = cube.collapsed(['longitude', 'latitude'],
                                    iris.analysis.SUM, weights=grid_areas)
    else:
        total_cube = cube.collapsed(['grid_longitude', 'grid_latitude'],
                                    iris.analysis.SUM, weights=grid_areas)
        
    total_cube.units *= cf_units.Unit('m2')
    total_cube.convert_units("kg/s")
    return total_cube.data

def plot(hires_cube, lowres_cube, regridding='', deltadt=None):
    """Plot cubes"""

    if hires_cube.data.max() < 1.e-6:
        levels = np.logspace(-16,-9,8)
    else:
        levels = np.arange(0.,10.5, 0.5)
    norm = mplcolors.BoundaryNorm(levels, 256)

    plt.figure(figsize=(11,4))
    
    ax = plt.gcf().add_subplot(121) 
    cf = iplt.pcolormesh(hires_cube, norm=norm)
    plt.gca().coastlines()
    total = "%.4f" % emission_total(hires_cube)
    plt.title('Hires emissions:\n'+total+'kg/s')
    
    ax = plt.gcf().add_subplot(122)
    iplt.pcolormesh(lowres_cube, norm=norm)
    plt.gca().coastlines()
    total = "%.3f" % emission_total(lowres_cube)
    title = 'Regridded emissions:\n'+regridding
    title += '\nTime: '+str(deltadt)+ 's'
    title += '\n '+total+'kg/s'
    plt.title(title)
    #plt.colorbar(cf, orientation='vertical',
    #             shrink=0.5, format = '%.1e')
    #plt.show()
    plt.savefig(regridding+'.png')
    plt.close()

def regrid(hires_cube, regridding, rotate=False):

    lowres_cube = None

    startdt = datetime.datetime.now()
    
    if regridding == 'iris.regrid.linear':
        lowres_cube = hires_cube.regrid(lowres_grid,iris.analysis.Linear())
    elif regridding == 'iris.regrid.linear.nan':
        lowres_cube = hires_cube.regrid(lowres_grid,iris.analysis.Linear('nan'))
    elif regridding ==  'iris.regrid.AreaWeighted':
        if rotate:
            print 'Doesnt work with changing coordinate systems'
        else:    
            lowres_cube = hires_cube.regrid(lowres_grid,iris.analysis.AreaWeighted())
    elif regridding == 'iris.experimental.area_weighted':
        if rotate:
            print 'Doesnt work with changing coordinate systems'
        else:    
            lowres_cube = iris.experimental.regrid.regrid_area_weighted_rectilinear_src_and_grid(
                             hires_cube, lowres_grid)
    elif regridding == 'iris.esmpy':
        lowres_cube = regrid_conservative_via_esmpy(hires_cube, lowres_grid)
    elif regridding == 'tidl.pp_regrid' :  
        lowres_cube = get_cube_from_file('lowres_ppregrid.pp')
    elif regridding == 'tidl.pp_regrid_avg' :
        lowres_cube = get_cube_from_file('lowres_ppregridavg.pp')

    enddt = datetime.datetime.now()
    deltadt = (enddt - startdt).total_seconds()

    return lowres_cube, deltadt
        
if __name__ == '__main__':

    #hires_cube = make_cube(0.1,'zero')
    #hires_cube.data = emission_shapes(hires_cube.data)
    
    #lowres_grid = make_cube(0.05,'zero')
    #dirs = [ '0.05','0.1','0.2','0.5','CAMS_hilev','CAMS_sfc','EMEP_to_hiresll','NAEIarea',
    #         'NAEIareapt','NAEIareapt_to_lowresll']
    #for path in dirs:
    #    os.chdir('/home/h06/apdl/python/play/'+path)

    hires_cube = get_cube_from_file('hires.pp')
    lowres_grid = get_cube_from_file('lowres_grid_input.pp')

    #Make even lower resolution
    lowres_grid = lowres_grid[::3,::3]
    lowres_grid.coord('grid_latitude').bounds=None
    lowres_grid.coord('grid_latitude').guess_bounds()
    lowres_grid.coord('grid_longitude').bounds=None
    lowres_grid.coord('grid_longitude').guess_bounds()

    
    #Save to pp so can call tidl routines
    #iris.save(hires_cube,'hires.pp')
    #iris.save(lowres_grid,'lowres_grid.pp')
    #stop

    regrid_options = ['iris.regrid.linear', 'iris.regrid.AreaWeighted',
                      'iris.experimental.area_weighted', 'iris.esmpy',
                      'iris.regrid.linear.nan','tidl.pp_regrid', 'tidl.pp_regrid_avg', ]
    #regrid_options = ['iris.regrid.linear.nan']

    for regridding in regrid_options:

        print regridding

        lowres_cube, deltadt = regrid(hires_cube, regridding, rotate=True)
        if lowres_cube is not None:
            plot(hires_cube, lowres_cube, regridding, deltadt)

    #regridding = 'iris.regrid.linear'
    #lowres_cube = hires_cube.regrid(lowres_grid,iris.analysis.Linear())

    #regridding = 'iris.regrid.AreaWeighted'
    #lowres_cube = hires_cube.regrid(lowres_grid,iris.analysis.AreaWeighted())

    #regridding = 'iris.experimental.area_weighted'
    #lowres_cube = iris.experimental.regrid.regrid_area_weighted_rectilinear_src_and_grid(
    #                 hires_cube, lowres_grid)

    #regridding = 'iris.esmpy'
    #lowres_cube = regrid_conservative_via_esmpy(hires_cube, lowres_grid)

    #regridding = 'tidl.pp_regrid'
    #lowres_cube = get_cube_from_file('lowres_ppregrid.pp')
    
    #regridding = 'tidl.pp_regrid_avg'
    #lowres_cube = get_cube_from_file('lowres_ppregridavg.pp')

    #plot(hires_cube, lowres_cube, regridding)

    
    #!!! Try point sources
    #!!! Try changing resolution.
