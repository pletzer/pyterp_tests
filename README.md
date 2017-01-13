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

Shown are execution times for the computation of the interpolation weights and the times it takes to apply 
the weights to the fiels (eval). 

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
| libcf/pycf    |  9.41         |  0.0333      |   4.73e-06  |
| ESMF/ESMPy    |  100.0        |  0.0216      |   9.26e-06  |     

#### Conservative (cell centered field)

|               |  weight sec   | eval sec     | error       |
| ------------- |---------------|--------------|-------------|
| sigrid        |               |              |             |
| ESMF/ESMPy    |               |              |             |

