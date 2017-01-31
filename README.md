# pyterp_tests: Python interpolation tests

## Overview

This repository contains curvilinear interpolation tests. 

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

