# pyterp_tests: Python interpolation tests

## Overview

## Prequisites

## How to run the tests

## Results

### Uniform nodal 2D field (double)

Source grid: 401 x 801

Target grid: 1201 x 2401

|               | time sec      | memory  MB   | error       |
| ------------- |---------------|--------------|-------------|
| Iris          |  2.46         |    1533      | 9.25e-06    |
| libcf/pycf    |  8.64+0.07    |     962      | 9.25e-06    |
| ESMF/ESMPy    |               |              |

### Nodal 2D field with a vertical axis

|               | uniform       | rectilinear  | curvilinear |
| ------------- |---------------|--------------|-------------|
| Iris          |               |              |             |
| libcf         |               |              |             |
| ESMPy         |               |              |             |

### Cell 2D field (conservative)

|               | uniform       | rectilinear  | curvilinear |
| ------------- |---------------|--------------|-------------|
| Iris          |               |              |             |
| libcf         |               |              |             |
| ESMPy         |               |              |             |

### Cell 2D field with a vertical axis (conservative)

|               | uniform       | rectilinear  | curvilinear |
| ------------- |---------------|--------------|-------------|
| Iris          |               |              |             |
| libcf         |               |              |             |
| ESMPy         |               |              |             |
