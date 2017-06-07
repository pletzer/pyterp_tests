import os

import numpy as np
from create_data_files import *
from ocgis import RequestDataset
from ocgis.util.helpers import create_exact_field_value

field = RequestDataset(os.path.join(WD, DST_FILENAME)).get()
lat_random = np.random.randint(0, high=field.grid.dimensions[0].size, size=20)
lon_random = np.random.randint(0, high=field.grid.dimensions[1].size, size=20)
dv = field.data_variables[0]
abs_errors = []
total = 20 * 20 * 3
ctr = 1
for tidx in range(field.time.shape[0]):
    for latr in lat_random:
        for lonr in lon_random:
            if ctr % 100 == 0 or ctr == 1:
                print('Processing {} of {}...'.format(ctr, total))
            sub = dv[tidx, latr, lonr]
            target_value = sub.get_value().flatten()[0]
            exact = create_exact_field_value(sub.parent.grid.x.get_value()[0], sub.parent.grid.y.get_value()[0])
            exact = exact + 10 * (tidx + 1)
            abs_errors.append(np.abs(exact - target_value))
            ctr += 1
abs_errors = np.array(abs_errors)
print('Minimum Absolute Relative Error', abs_errors.min())
print('   Mean Absolute Relative Error', abs_errors.mean())
print('Maximum Absolute Relative Error', abs_errors.max())
print('                    Sample Size', abs_errors.size)
