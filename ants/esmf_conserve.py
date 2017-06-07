import iris
import ants.decomposition as decomp
from ants.regrid import _ESMFRegridder
import ants.utils.cube

iris.FUTURE.netcdf_promote = True

def regrid_op(cube1, cube2):
    regridder = _ESMFRegridder(cube1, cube2)
    return regridder(cube1)

#filename = 'ESACCI-LC-L4-LCCS-Map-300m-P5Y-2010-v1.6.1-inlandwater.nc'
#varname = 'lccs_class'
#src_cube = iris.load(filename, iris.Constraint(cube_func = lambda c: c.var_name == varname))[0]

src_filename = 'src.nc'
src_varname = 'cellData'
src_cube = iris.load(src_filename, iris.Constraint(cube_func = lambda c: c.var_name == src_varname))[0]
ants.utils.cube.set_crs(src_cube)
dst_filename = 'dst.nc'
dst_varname = 'cellData'
dst_cube = iris.load(dst_filename, iris.Constraint(cube_func = lambda c: c.var_name == dst_varname))[0]
ants.utils.cube.set_crs(dst_cube)

#splitter = decomp.MosaicBySplit(src_cube, (2, 4))
cube = decomp.decompose(regrid_op, src_cube, dst_cube, split=(2, 4))
#print(list(splitter()))
#print(cube.data.dtype)

