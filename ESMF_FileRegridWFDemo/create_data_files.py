import os

import numpy as np
from ocgis import RequestDataset
from ocgis.spatial.grid_splitter import GridSplitter
from ocgis.test.base import create_exact_field


# ------------------------------------------------------------------------------
# Customize these path locations.

WD = os.getcwd()
DATADIR = '/media/benkoziol/Extra Drive 1/data/bekozi-work/i39-tripole-grids'
OUTDIR = os.path.join(WD, 'data')

# ------------------------------------------------------------------------------
# These globals can likely remain intact.

TRIPOLE_PATH = os.path.join(DATADIR, 'coords_CF_ORCA12_GO6-2.nc')
DST_FILENAME = os.path.join(DATADIR, 'dst.nc')
DECOMP = (10, 10)
DMAP_TRIPOLE = {'x': {'variable': 'lont', 'dimension': ['x'], 'bounds': 'lont_bounds'},
                'y': {'variable': 'latt', 'dimension': ['y'], 'bounds': 'latt_bounds'}}
DMAP_RECT = {'x': {'variable': 'lonMid', 'dimension': ['dim1_0']}, 'y': {'variable': 'latMid', 'dimension': ['dim0_0']}}

SRC_TEMPLATE = os.path.join(OUTDIR, 'src_split_{}.nc')
DST_TEMPLATE = os.path.join(OUTDIR, 'dst_split_{}.nc')
INDEX_FILENAME = os.path.join(OUTDIR, '01-split-index.nc')
WGT_TEMPLATE = 'weights_{}.nc'


if __name__ == '__main__':
    # --------------------------------------------------------------------------
    # Get fields and remove unneeded variables.

    src_rd = RequestDataset(TRIPOLE_PATH, dimension_map=DMAP_TRIPOLE, grid_abstraction='point')
    dst_rd = RequestDataset(DST_FILENAME, dimension_map=DMAP_RECT, grid_abstraction='point')

    src_field = src_rd.get()
    dst_field = dst_rd.get()

    dst_to_keep = ['lonMid_bnds', 'lonMid', 'latMid_bnds', 'latMid', 'latitude_longitude']
    dst_to_pop = set(dst_field.keys()).difference(set(dst_to_keep))
    for d in dst_to_pop:
        dst_field.pop(d)

    src_to_keep = ['latt', 'lont', 'latt_bounds', 'lont_bounds', 'latitude_longitude']
    src_to_pop = set(src_field.keys()).difference(set(src_to_keep))
    for s in src_to_pop:
        src_field.pop(s)

    # --------------------------------------------------------------------------
    # Create exact field values for testing.

    src_efield = create_exact_field(src_field.grid, 'exact', ntime=3, fill_data_var=True)
    dst_efield = create_exact_field(dst_field.grid, 'exact', ntime=3, fill_data_var=False)
    assert np.all(dst_efield['exact'].get_value() == dst_efield['exact'].fill_value)

    # --------------------------------------------------------------------------
    # Execute grid splitting routine.

    if not os.path.exists(OUTDIR):
        os.mkdir(OUTDIR)

    gs = GridSplitter(src_efield.grid, dst_efield.grid, DECOMP, check_contains=False, allow_masked=True)
    gs.write_subsets(SRC_TEMPLATE, DST_TEMPLATE, WGT_TEMPLATE, INDEX_FILENAME)
