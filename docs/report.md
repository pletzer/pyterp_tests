# Comparison of regridding packages used for climate/weather data

Alex Pletzer (NIWA/NeSI), Chris Scott (NeSI) and Jamie Kettleborough (UK Met Office)

13 March 2017

## Overview

Regridding, that is the process of mapping fields from one grid to another,
is widely used in pre- and post-processing tools. Here we provide a comparison of existing and newly developed regridding tools. We 
will look at:

 * Accuracy
 * Numerical performance (execution time and memory footprint)
 * Ability to handle masking (invalid data)

Our emphasis will be on general, 2D structured source grids. Uniform latitude/longitude grids will be used as target grids (most regridding packages support any type of target grids). 

The regridding tools we consider are:

 * **iris** 1.10.0-DEV's built-in regridding. Iris is a community driven Python library for analyzing earth science data sets. Iris can be installed with `conda install -c scitools iris`.
 * **libcf** 1.6.9 and its Python interface. Libcf was developed at UCAR to address the need to produce and read CF compliant files. The library also supports regriding. Libcf can be installed with `pip install pycf`. 
 * **sigrid** is a github project that computes the intersection of structured grids: `git clone https://github.com/pletzer/sigrid && cd sigrid && python setup.py install` 
 * The Earth System Modeling Framework **ESMF** 7.0 is a high performance software for building coupled earth modeling applications. ESMF includes a regridding class with a Python callable interface: https://www.earthsystemcog.org/projects/esmf/. ESMF can be installed with 
 `conda install -c conda-forge esmpy=7.0.0`.

Two regridding methods are considered: _bilinear_ and _conservative_. Bilinear is suitable for nodal 
 data whereas conservative should be applied to cell centred data. Conservative regridding 
 conserves global quantities such as mass, energy, etc. None of these methods are appropriate for 
 vector fields for which the components are staggered (eg Arakawa C-grid staggering).

 The capabilities of the considered regridding tools are summarized below:

|               |  grid type    |   bilinear?   | conservative? |  stores weights? |
|---------------|---------------|--------------|----------------|-----------------|
|  iris         |  rectilinear  |    yes       |     yes        |     no          | 
|  libcf        |  structured  |    yes        |    no         |     yes         |  
| sigrid        |  structured  |    no         |    yes        |     yes         |
| ESMF          |  structured  |    yes        |    yes        |     yes         |

Note that some tools have capabilities we have not tested. For instance, libcf supports 
interpolation in n-dimensions while ESMF supports interpolation from and onto unstructured grids in 
2D and 3D. Also note that structured grids (regular topology/irregular points) are more general than rectilinear grids
(regular topology and points). A rectilinear grid is a structured grid but not vice versa.

## Accuracy

### Bilinear regridding accuracy

A simple sinusoidal field on a uniform source grid with invalid data inside a quarter of a disk is chosen. Shown are 
the source field values on nodes (large spheres) and their interpolated values on the much finer destination grid (small spheres).
Invalid points are shown as grey cubes. 

For ESMF, destination points falling inside partially valid source cells, i.e. cells that have at least one invalid node, are not interpolated. In the case of libcf on the other hand, destination points will be interpolated if these fall within a valid triangular subcell. As such, libcf will give a smoother transition from valid to invalid regions.


ESMF bilinear with masked data | libcf bilinear with masked data
:-----------------------------:|:--------------------------------:
![alt text](https://github.com/pletzer/pyterp_tests/blob/master/masking/vis_esmf1_dst.png "ESMF bilinear regridding of masked field")                        | ![alt text](https://github.com/pletzer/pyterp_tests/blob/master/masking/vis_libcf1_dst.png "libcf bilinear regridding of masked field")


### Conservative regridding accuracy

To evaluate the accuracy of conservative regridding, we chose a polar source grid with a sine wave imprinted on it. The 
Cartesian destination grid was chosen to be entirely contained within the source grid as shown below. The picture to the right
shows the difference between the regridded ESMF and sigrid fields, which is of order 1.e-10 or smaller.

Source grid (black) and destination field/grid | Difference between ESMF and sigrid regridding
:---------------------------------------------:|:---------------------------------------------:
![alt text](https://github.com/pletzer/pyterp_tests/blob/master/polar/srcAndDst.png "Src and dst grids") | ![alt text](https://github.com/pletzer/pyterp_tests/blob/master/polar/esmfMinusSigrid.png "Difference between ESMF and sigrid")


## Performance

### Uniform to uniform grid

The execution time required to regrid a source field on a uniform grid to a uniform destination grid is measured. Both the source and destination grid resolutions are varied.

Conservative regridding with ESMF is one order of magnitude slower than sigrid. Bilinear ESMF regridding is 2x faster than conservative ESMF. Libcf is nearly an order of magnitude faster than bilinear ESMF and bilinear iris is another two aroders of magnitudes faster. 

A surprising result is that conservative iris and ESMF regridding run neck to neck - there is no performance penalty in using ESMF, which is more general than iris's regridding. Recall that iris's conservative regridding is restricted to uniform source and destination grids.

| Comparing the cost of different regridding schemes |
|:--------------------------------------------------:|
|![alt text](https://github.com/pletzer/pyterp_tests/blob/master/uniform/run.png "comparing the execution times of different regridding methods and packages")|


### Rotated pole to uniform grid

Rotated pole grids are widely used in regional models as a way to circumvent the problem of poles. These grids share many 
characteristics of curvilinear grids, including significant bending of the grid lines and the presence of poles, which tend to 
be problematic for many regridding tools. Shown below are the execution times split between the computation of interpolation 
weights (wgts) and the application of the weights to the source field (eval) for bilinear and conservative regridding.

The evaluation step is typically orders of magnituds faster than the computation of weights,
indicating the need to reuse the weights whenever possible. 

Libcf is faster than ESMF for high resolution (total number of source and destination cells larger than 1e8) but this advantage mostly disappears at very high resolutions. Conservative interpolation is only a facto 2x slower than ESMF bilinear at very high resolution. 

Source (green) and destination (red) grids     | Regridding execution times
:---------------------------------------------:|:---------------------------------------------:
![alt text](https://github.com/pletzer/pyterp_tests/blob/master/rotated_pole/rotated_pole_grid.png "rotated pole to uniform grid") | ![alt text](https://github.com/pletzer/pyterp_tests/blob/master/rotated_pole/run.png "execution times for rotated pole to uniform regridding")


### Tripolar grid to uniform grid

The source is a tripolar grid of fixed resolution 3606 x 4322 with the destination grid's resolution being varied. Shown are the execution
times for computing the conservative interpolation weights (solid lines) and the time to evaluate the conservative 
interpolation (dashed lines). Running in parallel only moderately reduces the overall wall clock time. 

Also shown is the peak memory consumption returned by the SLURM scheduler on Pan (magenta) and the value obtained by eye balling the 
Unix `top` command (cyan) while the application ran. As such the cyan curve has much lower sampling frequency and should be regarded 
as a time average of the cyan curve. From the data required to store the source grid, its coordinates and the fields we the meory footprint
to be 0.5GB (compare to 10GB observed). Thus, conservative regridding requires 20-50x more memory than would be required to store the data only.

| ESMF conservative serial vs MPI 4 processor execution                |  ESMF conservative memory consumption (serial)                                                           
|:--------------------------------------------------:|:------------------------------------------------------------------------------:
|![alt text](https://github.com/pletzer/pyterp_tests/blob/master/big/run_conserve.png "tripolar to uniform conservative regridding") | ![alt text](https://github.com/pletzer/pyterp_tests/blob/master/big/memory.png "Memory consumption (serial)")

## Summary and recommendations

 * The method of regridding depends on the field. Cell centred fields with bounds arrays require conservative interpolations as they associate average values across cells. Field values attached to cell nodes should use bilinear interpolation.

 * Regardless of whether bilinear interpolation or conservative interpolation is used, the computation of the interpolation weights is typically orders of magnitudes slower than the evaluation of the interpolation. Hence, we recommend to split the weight computation from the interpolation proper. This will allow users to amortize the cost of computing the interpolation weights, particularly in scenarios where the field is time dependent or has a vertical axis. The same interpolation weights can also be used for other fields provided they share the same staggering and grids. 

 * Only one package (ESMF) supports bilinear and conservative interpolation on general structured grids. ESMF also supports regridding on 
 unstructured grids in up to 3 dimensions and second order accutate interpolation for nodal data. Second order accurate conservative interpolation will become available in version 7.1.

 * In contrast to libcf, ESMF bilinear will not interpolate data in cells which have invalid nodes. 

 * There appears some initial overhead in computing the weights with ESMF, which appears to depend mostly on the source grid resolution with little dependence on the target grid resolution. Libcf, which uses a Newton algorithm to locate the target cells, has more favourable scaling for 
 low resolution source grid but this advantage almost vanishes as the target grid resolution exceeds the source grid resolution.

 * Our recommendation is to implement ESMF conservative and bilinear into iris with support for general structured grids (both source and destination). Ideally, regridding should be a class taking source and destination grids at construction. Given the cost of computing the interpolation weights, 
 it would be advantageous to have methods to store the weights to disk and load the weights from disk. Finally there should be an evaluation step 
 which fills the field with interpolated values.


