import os

from create_data_files import *
from ocgis.spatial.grid_splitter import GridSplitter

index_path = os.path.join(WD, INDEX_FILENAME)
dst_netcdf = os.path.join(WD, DST_FILENAME)

GridSplitter.insert_weighted(index_path, WD, dst_netcdf)
