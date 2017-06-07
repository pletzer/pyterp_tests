# README - ESMF_FileRegridWFDemo

## License Information

**Earth System Modeling Framework (ESMF)**

* Copyright 2002-2017, University Corporation for Atmospheric Research, Massachusetts Institute of Technology, Geophysical Fluid Dynamics Laboratory, University of Michigan, National Centers for Environmental Prediction, Los Alamos National Laboratory, Argonne National Laboratory, NASA Goddard Space Flight Center.

* Licensed under the University of Illinois-NCSA License.

## Requirements

Please see the [Step-by-Step](#step-by-step) section for a full tutorial on running this external demo.

* ESMF must be built and installed for this testing application to run properly.
* The ESMF library must be built with NetCDF support.
* If the ESMF library has already been built on the system you are using, you should only have to set the `ESMFMKFILE` environment variable. This environment variable should point to the esmf.mk file in the *installed* version of the ESMF library.
* OCGIS must be installed on your system. See the [Step-by-Step](#step-by-step) section for instructions.

Some useful commands you can issue in this directory:

`make`       - build the application
`make dust`  - clean out .o and .mod files
`make clean` - clean the directory and start over

Please contact esmf_support@list.woc.noaa.gov with any questions or problems.

## Tutorial

Paths and variable names are hardcoded into the source files. If you wish to change the paths to match the location of your data, modify the source files and replace the paths/names with your requirements.

* Python variables are located at the top of `create_data_files.py`.
* ESMF variables are located at the top of `ESMF_FileRegridWFDemo.F90`.
* Both must be adjusted if filenames, etc. change.

### Limitations and Recommendations

* Works with structured grids and spherical coordinate systems only.
* Works with one and only one data variable.
* Works with three-dimensional data only with dimension order `(time, latitude, longitude)`. The `time` dimension *must* be named `time`.
* It is highly recommended that any external data used in this procedure closely follows the example dataset metadata created in this demo.

### Step-by-Step

* Build the demo code.

```
export ESMFMKFILE=<path>/esmf.mk
git clone git://git.code.sf.net/p/esmf/external_demos esmf-external_demos
cd esmf-external_demos/ESMF_FileRegridWFDemo
make clean
make
```

* Install OCGIS from the repository master branch. Full installation documentation available here: [http://ocgis.readthedocs.io/en/latest/install.html](http://ocgis.readthedocs.io/en/latest/install.html).

```
conda install -c conda-forge ocgis mpi4py
conda remove ocgis
git clone https://github.com/NCPP/ocgis.git
cd ocgis
python setup.py install
```

* Create input data files then split them into smaller source and destination files.

```
# Creates example datasets in the ./data directory.
python create_data_files.py
  *or*
mpirun -n <nprocs> python create_data_files.py
```

* Run the regridder and weighting utility.

```
# Generates regridding weights and applies the sparse matrix.
mpirun -n <procs> ./ESMF_FileRegridWFDemo
```

* Insert the weighted data into the master destination file.

```
# Takes data from regridded variable and inserts it into the master 
# destination file.
python insert_weighted_data.py
```

* Optionally, run validation to print absolute relative errors using an exact field comparison.

```
python validate.py
```
