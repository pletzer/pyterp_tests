import iris

filename = 'ESACCI-LC-L4-LCCS-Map-300m-P5Y-2010-v1.6.1-inlandwater.nc'
varname = 'lccs_class'
cube = iris.load(filename, iris.Constraint(cube_func = lambda c: c.var_name == varname))[0]
print(cube)


