**Under Development**

![myfempy_logo](figures\logoB.png)

# Welcome to myfempy\'s web documentation!

Copyright © Antonio Vinicius G. Campos 2022. Processo INPI BR512022001484-0

# About

**Myfempy** is a python package based on finite element method for
scientific analysis. The code is open source and *intended for
educational and scientific purposes only, not recommended to commercial
use*. You can help us by contributing with a donation on the main
project page, read the support options. **If you use myfempy in your
research, the developers would be grateful if you could cite in your
work.**

# Installation

## To install myfempy manually in your directory, following the steps

1.  Clone/ Download the main code \[latest version\] from
    [github/myfempy/main](https://github.com/easycae-3d/myfempy/)
2.  Unzip the pack in your preferred location
3.  In the **myfempy-main** folder, open a terminal and enter with the
    command:

``` bash
>> python -m pip install --upgrade pip

>> pip install .

or

>> python -m pip install --upgrade build

>> python -m build
```

**Note: is recommend to create a virtual environment previously the
installation of** myfempy\*\* and dependencies packs. You can use the
[virtualenv](https://virtualenv.pypa.io/en/latest/) or [conda
environments](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html)\*\*

# Dependencies

**Myfempy** can be used in systems based on Linux, MacOS and Windows.
**Myfempy** requires Python 3.

## Installation prerequisites, required to build **myfempy**

You can use either of two python development environments to run myfempy

-   [Python 3.x](https://www.python.org/) - *Python is a programming
    language that lets you work quickly and integrate systems more
    effectively.*
-   [Anaconda](https://www.anaconda.com/) - *Anaconda offers the easiest
    way to perform Python/R data science and machine learning on a
    single machine.*

## Python packages required for using **myfempy**

The following python packages are required to run myfempy. Before to
install myfempy-main, install this packages. Check if they are already
installed on your machine

-   [numpy](https://numpy.org/) - The fundamental package for scientific
    computing with Python
-   [cython](https://cython.org/) - Cython is a language that makes
    writing C extensions for Python as easy as Python itself
-   [scipy](https://scipy.org/) - Fundamental algorithms for scientific
    computing in Python
-   [vedo](https://vedo.embl.es/) - A python module for scientific
    analysis and visualization of эd objects
-   [vtk](https://pypi.org/project/vtk/)(optional) - VTK is an
    open-source toolkit for 3D computer graphics, image processing, and
    visualization
-   try

``` bash
>> pip install numpy, cython, scipy, vedo
```

### Outhers prerequisites

-   [gmsh/External Generator Mesh](https://gmsh.info/) - Gmsh is an open
    source 3D finite element mesh generator with a built-in CAD engine
    and post-processor. *Notes: 1 - Gmsh is NOT part of myfempy
    projects; 2 - Is Needed install Gmsh manually*
-   try

``` bash
>> pip install --upgrade gmsh
```

-   [gmsh PyPi](https://pypi.org/project/gmsh/)

# Tutorial

A **Basic Tutorial** is available
[here](https://myfempy.readthedocs.io/en/1.dev9/tutorial.html).

Many **Examples** are available
[here](https://github.com/easycae-3d/myfempy/tree/master/examples).

# Documentation

The myfempy is documented using Sphinx under **doc**. The myfempy\'s
documents versions can be found in html, pdf or epub.

The **Web Documentation** is available on \[Read the
Docs\](<https://myfempy.readthedocs.io/>).

To compile the documentation use *sphinx* in the **doc** folder. Do,

``` bash
>> make html {in the root folder where the index.rst file is} --> This command generates *.html* files

>> make latexpdf # [optional] todo doc pdf
```

# Release

The all release versions is available
[here](https://github.com/easycae-3d/myfempy/releases)

# Features

The *main myfempy features* are available here:

-   [Features
    List](https://github.com/easycae-3d/myfempy/blob/main/docs/Myfempy%20Features.pdf)

# License

**myfempy** is published under the [GPLv3
license](https://www.gnu.org/licenses/gpl-3.0.en.html). See the
[myfempy/LICENSE](https://github.com/easycae-3d/myfempy/blob/main/LICENSE.txt).

```{=html}
<!-- ## >> Acknowledgment -->
```
# Citing

Have you found this software useful for your research? Star the project
and cite it as:

-   APA:

``` bash
Antonio Vinicius Garcia Campos. (2022). myfempy (1.5.1). Zenodo. https://doi.org/10.5281/zenodo.6958796
```

-   BibTex:

``` bash
@software{antonio_vinicius_garcia_campos_2022_6958796,
author       = {Antonio Vinicius Garcia Campos},
title        = {myfempy},
month        = aug,
year         = 2022,
publisher    = {Zenodo},
version      = {1.5.1},
doi          = {10.5281/zenodo.6958796},
url          = {https://doi.org/10.5281/zenodo.6958796}
}
```

# References

-   [Myfempy](https://myfempy.readthedocs.io/) - *A python package for
    scientific analysis based on finite element method.*
-   [FEM](https://en.wikipedia.org/wiki/Finite_element_method) - *The
    finite element method (FEM) is a popular method for numerically
    solving differential equations arising in engineering and
    mathematical modeling.*
-   [Solid Mechanics](https://en.wikipedia.org/wiki/Solid_mechanics)
    -*Solid mechanics, also known as mechanics of solids, is the branch
    of continuum mechanics that studies the behavior of solid materials,
    especially their motion and deformation under the action of forces,
    temperature changes, phase changes, and other external or internal
    agents.*
-   [PDE](https://en.wikipedia.org/wiki/Partial_differential_equation)
    \- *In mathematics, a partial differential equation (PDE) is an
    equation which imposes relations between the various partial
    derivatives of a multivariable function.*

------------------------------------------------------------------------

# Changelog

The changelog is available
[here](https://github.com/easycae-3d/myfempy/wiki/Changelog)

# Project tree structure

``` bash
/myfempy
|   version.py
|   __init__.py
|
+---core
|       assembler.py
|       solver.py
|       solverset.py
|       staticlinear.py
|       vibralinear.py
|       __init__.py
|
+---felib
|   |   crossec.py
|   |   felemset.py
|   |   materset.py
|   |   physicset.py
|   |   quadrature.py
|   |   __init__.py
|   |
|   +---fluid
|   |       __init__.py
|   |
|   +---fsi
|   |       __init__.py
|   |
|   +---materials
|   |       axial.py
|   |       lumped.py
|   |       planestrain.py
|   |       planestress.py
|   |       solid.py
|   |       __init__.py
|   |
|   +---physics
|   |       force2node.py
|   |       getnode.py
|   |       loadsconstr.py
|   |       __init__.py
|   |
|   \---struct
|           beam21.py
|           frame21.py
|           frame22.py
|           plane31.py
|           plane41.py
|           solid41.py
|           solid81.py
|           spring21.py
|           truss21.py
|           __init__.py
|
+---io
|       filters.py
|       iomsh.py
|       iovtk.py
|       __init__.py
|
+---mesh
|       genmesh.py
|       gmsh.py
|       legacy.py
|       __init__.py
|
+---plots
|       meshquality.py
|       physics.py
|       plotmesh.py
|       plotxy.py
|       postplot.py
|       prevplot.py
|       __init__.py
|
+---postprc
|       displcalc.py
|       postcomp.py
|       postset.py
|       __init__.py
|
\---tools
        logo.png
        logo.txt
        path.py
        tools.py
        __init__.py
```