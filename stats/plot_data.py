

def plotData(dataId):
	from matplotlib import pylab
	xtypep = c_int()
	fillValuePtr = c_void_p()
	gridId = c_int()
	coordIds = (c_int * ndims)()

	dataPtr = POINTER(c_double)()
	latPtr = POINTER(c_double)()
	lonPtr = POINTER(c_double)()

	ier = pycf.nccf.nccf_get_data_pointer(dataId, byref(xtypep),
                                          byref(dataPtr), byref(fillValuePtr))
	assert(ier == pycf.NC_NOERR)
	ier = pycf.nccf.nccf_get_data_pointer(dataId, byref(xtypep),
                                          byref(dataPtr), byref(fillValuePtr))
	assert(ier == pycf.NC_NOERR)
	ier = pycf.nccf.nccf_inq_data_gridid(dataId, byref(gridId))
	assert(ier == pycf.NC_NOERR)
	ier = pycf.nccf.nccf_inq_grid_coordids (gridId, coordIds)
	assert(ier == pycf.NC_NOERR)
	ier = pycf.nccf.nccf_get_coord_data_pointer(coordIds[0], byref(latPtr))
	assert(ier == pycf.NC_NOERR)
	ier = pycf.nccf.nccf_get_coord_data_pointer(coordIds[1], byref(lonPtr))
	assert(ier == pycf.NC_NOERR)

	ntot, dims = inquireDataSizes(dataId)
	data = numpy.ctypeslib.as_array(dataPtr, shape=tuple(dims))
	lats = numpy.ctypeslib.as_array(latPtr, shape=tuple(dims))
	lons = numpy.ctypeslib.as_array(lonPtr, shape=tuple(dims))

	pylab.pcolor(lons, lats, data)
	pylab.show()

