# Comparison of regridding packages used for climate/weather data

## Overview

Regridding, that is the process of mapping fields from a source grid to another destination grid with different resolution,
is widely used in pre- and post-processing tools. Here we provide a comparison of regridding tools in terms of:

 * Accuracy
 * Numerical Performance
 * And their ability to handle masking (invalid data)

 We will use 2D structured grids for input data and typically uniform latitude/longitude grids as target grids. 

The regridding tools we considered are:

 * Iris built-in regridding
 * libcf 1.6.9
 * sigrid
 * ESMF 7.x

 Their capabilities are summarized below:

|               |  grid type    |   bilinear?   | conservative? |  store weights? |
|---------------|---------------|--------------|----------------|-----------------|
|  iris         |  uniform     |    yes       |     yes        |     no          | 
|  libcf        |  structured  |    yes        |    no         |     yes         |  
| sigrid        |  structured  |    no         |    yes        |     yes         |
| ESMF          |  structured  |    yes        |    yes        |     yes         |

Note that some tools have capabilities we have not tested. For instance, libcf supportd multilinear 
interpolation in n-dimensions and ESMF supports interpolation from and onto unstructured grids in 
2D and 3D. 


## Prequisites

 * Anaconda python 2.7 with numpy 1.11.1 or later
 * Iris 1.10.0-DEV. Install with `conda install -c scitools iris`
 * ESMF/ESMPy 7.0.0 or later [https://www.earthsystemcog.org/projects/esmpy/]
     
Typical build instructions on Linux

```cd esmf; export ESMF_DIR=$PWD; export ESMF_COMM=mpich2; export ESMF_INSTALL_PREFIX=<path-to-esmf>```

```make; make install```

```cd src/addon/ESMPy/; python setup.py build --ESMFMKFILE=<path-to-esmf>/lib/libO/<platform>/esmf.mk install```
     
 * libcf/pycf 1.6.5 or later. Install with `pip install pycf`
 * sigrid 0.1.0 or later. `git clone https://github.com/pletzer/sigrid && cd sigrid && python setup.py install`

## Results

Shown are execution times for the computation of the interpolation weights and the times it takes to apply 
the weights to the fields (eval). 

### Polar grid to uniform grid

Source grid: 400 x 800

Target grid: 200 x 400


#### Conservative (cell centered field)


|               |  weight sec   | eval sec     | error       |
| ------------- |---------------|--------------|-------------|
| sigrid        |               |              |             |
| ESMF/ESMPy    |  7.86         |   0.003      | -0.0013     |


### Rotated pole grid to uniform grid

Source grid: 1281 x 2560

Target grid: 641 x 1281


#### Bilinear (nodal field)

|               |  weight sec   | eval sec     | error       |
| ------------- |---------------|--------------|-------------|
| libcf/pycf    |  12.9         |  0.028       |   2.8e-06   |
| ESMF/ESMPy    |  234          |  0.53        |   9.3e-06   |

#### Conservative (cell centered field)

|               |  weight sec   | eval sec     | error       |
| ------------- |---------------|--------------|-------------|
| sigrid        |               |              |             |
| ESMF/ESMPy    |               |              |             |

### ORCA grid 

Source grid: 3606 x 4322

Target grid: 720 x 1440


#### Bilinear (nodal field)

|               |  weight sec   | eval sec     |
| ------------- |---------------|--------------|
| libcf/pycf    |  537          |  0.025       |
| ESMF/ESMPy    |  229          |  0.024       |

Source grid: 3606 x 4322

Target grid: 3601 x 7201

#### Bilinear (nodal field)

|               |  weight sec   | eval sec     |
| ------------- |---------------|--------------|
| libcf/pycf    |  1.3e+04      |  0. 38       |
| ESMF/ESMPy    | 1.82e+03      |  0.40        |



