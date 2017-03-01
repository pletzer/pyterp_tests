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

 The considered interpolation methods are bilinear and conservative. Bilinear is suitable for nodal 
 data whereas conservative should be applied to cell centred data, for which there is a need to 
 conserve global quantities sucha as mass, energy, etc. None of these methods are appropriate for 
 vector fields such as the velocity on a Arakawa C grid for which the components are staggered.

 The capabilities of the regridding tools are summarized below:

|               |  grid type    |   bilinear?   | conservative? |  store weights? |
|---------------|---------------|--------------|----------------|-----------------|
|  iris         |  uniform     |    yes       |     yes        |     no          | 
|  libcf        |  structured  |    yes        |    no         |     yes         |  
| sigrid        |  structured  |    no         |    yes        |     yes         |
| ESMF          |  structured  |    yes        |    yes        |     yes         |

Note that some tools have capabilities we have not tested. For instance, libcf supports linear 
interpolation in n-dimensions while ESMF supports interpolation from and onto unstructured grids in 
2D and 3D. 

## Accuracy

### Bilinear

#### No masking

#### Masking

### Conservative

## Performance

### Uniform to uniform grid

### Rotated pole to uniform grid

### Tripolar grid to uniform grid

## Summary and recommendations



