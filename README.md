# pyterp_tests: Python interpolation tests

## Overview

This repository contains curvilinear interpolation tests. 

## Prequisites

 * Anaconda python 2.7 with numpy 1.11.1 or later
 * Iris 1.10.0-DEV. Install with `conda install -c scitools iris`
 * ESMF/ESMPy 7.0.0 or later [https://www.earthsystemcog.org/projects/esmpy/]
 * libcf/pycf 1.6.1 or later. Install with `pip install pycf`
 * sigrid 0.1.0 or later. `git clone https://github.com/pletzer/sigrid && cd sigrid && python setup.py install`

## Results

### Polar grid to uniform grid

Source grid: 

Target grid: 


#### Conservative (cell centered field)


|               |  weight sec   | eval sec     | error       |
| ------------- |---------------|--------------|-------------|
| sigrid        |               |              |             |
| ESMF/ESMPy    |               |              |             |


### Rotated pole grid to uniform grid

Source grid: 

Target grid: 


#### Bilinear (nodal field)

|               |  weight sec   | eval sec     | error       |
| ------------- |---------------|--------------|-------------|
| libcf/pycf    |               |              |             |
| ESMF/ESMPy    |               |              |             |     

#### Conservative (cell centered field)

|               |  weight sec   | eval sec     | error       |
| ------------- |---------------|--------------|-------------|
| sigrid        |               |              |             |
| ESMF/ESMPy    |               |              |             |

