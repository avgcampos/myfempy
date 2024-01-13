
# Under Development

<!-- ![logo_v1](https://user-images.githubusercontent.com/54820276/159730160-871ac41a-958a-4398-a014-506619c4cb56.png) -->

<!-- ![logoB](https://user-images.githubusercontent.com/54820276/181291975-38e2a37e-be4f-40a3-b333-09724c61827b.png) -->

![myfempy_logo](docs/assets/logo.png)

Copyright Antonio Vinicius G. Campos and 3D EasyCAE, 2022

[![Python Versions]()]()
[![Documentation Status](https://readthedocs.org/projects/myfempy/badge/?version=latest)](https://myfempy.readthedocs.io/en/latest/?badge=latest)
[![PyPI]()]() [![conda]()]()
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6958796.svg)](https://doi.org/10.5281/zenodo.6958796)
[![Quality Gate Status](https://sonarcloud.io/api/project_badges/measure?project=avgcampos_myfempy&metric=alert_status)](https://sonarcloud.io/summary/new_code?id=avgcampos_myfempy)
[![Downloads]()
[![lics](https://img.shields.io/badge/license-GPL-blue.svg)](https://en.wikipedia.org/wiki/GNU_General_Public_License)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg?style=flat&logo=appveyor)]([https://github.com/psf/black](https://github.com/psf/black))


-----------

## About

**Myfempy** is a python package based on finite element method to multiphysics analysis. The code is open source and *intended for educational and scientific purposes only, not recommended to commercial use*. You can help us by contributing with a donation on the main project page, send us a email [3deasycaebr.contato@gmail.com]. **If you use myfempy in your research, the  developers would be grateful if you could cite in your work.**

## Installation

### To install myfempy manually in your directory, following the steps

1. Clone/ Download the main code [latest version] from [github/myfempy/main](https://github.com/easycae-3d/myfempy/)

2. Unzip the pack in your preferred location

3. In the **myfempy-main** folder, open a terminal and enter with the command:

```bash

>> python -m pip install --upgrade pip

>> pip install .

or

>> python -m pip install --upgrade build

>> python -m build

```

**Note: is recommend to create a virtual environment previously the installation of **myfempy** and dependencies packs. You can use the [virtualenv](https://virtualenv.pypa.io/en/latest/) or [conda environments](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html)**

## Dependencies

**Myfempy** can be used in systems based on Linux, MacOS and Windows. **Myfempy** requires Python 3.

### Installation prerequisites, required to build **myfempy**

You can use either of two python development environments to run myfempy

- [Python 3.x](https://www.python.org/) - *Python is a programming language that lets you work quickly and integrate systems more effectively.*
- [Anaconda](https://www.anaconda.com/) - *Anaconda offers the easiest way to perform Python/R data science and machine learning on a single machine.*

### Python packages required for using **myfempy**

The following python packages are required to run myfempy. Before to install myfempy-main, install this packages. Check if they are already installed on your machine

- [numpy](https://numpy.org/) - The fundamental package for scientific computing with Python
- [cython](https://cython.org/) - Cython is a language that makes writing C extensions for Python as easy as Python itself
- [scipy](https://scipy.org/) - Fundamental algorithms for scientific computing in Python
- [numdifftools](https://pypi.org/project/numdifftools/) - Solves automatic numerical differentiation problems in one or more variables.
- [vedo](https://vedo.embl.es/) - A python module for scientific analysis and visualization of эd objects
- [vtk](https://pypi.org/project/vtk/)(optional) - VTK is an open-source toolkit for 3D computer graphics, image processing, and visualization


- try

```bash

>> pip install numpy

```

#### Outhers prerequisites

- [gmsh/External Generator Mesh](https://gmsh.info/) - Gmsh is an open source 3D finite element mesh generator with a built-in CAD engine and post-processor. *Notes: 1 - Gmsh is NOT part of myfempy projects;  2 - Is Needed install Gmsh manually*

- try

```bash

>> pip install --upgrade gmsh

```

- [gmsh PyPi](https://pypi.org/project/gmsh/)

## Tutorial

A **Basic Tutorial** is available [here](https://myfempy.readthedocs.io/en/1.dev9/tutorial.html).

Many **Examples** are available [here](https://github.com/easycae-3d/myfempy/tree/master/examples).

## Documentation

The myfempy is documented using Mkdocs under `docs/`. The myfempy's documents versions can be found in html, pdf or epub.

The **Web Documentation** is available on [Read the Docs](https://myfempy.readthedocs.io/).

The **User's Manual**(pdf) is available on [manual_myfempy](https://github.com/easycae-3d/myfempy/blob/main/doc/manual%20myfempy.pdf).

To compile the documentation use *sphinx* in the **\doc** folder. Do,

```bash

>> make html {in the root folder where the index.rst file is} --> This command generates *.html* files

```

## Release

The all release versions is available [here](https://github.com/easycae-3d/myfempy/releases)

## Features

The *main myfempy features* are available here:

- [Features List](https://github.com/easycae-3d/myfempy/blob/main/docs/Myfempy%20Features.pdf)

## License

**myfempy** is published under the [GPLv3 license](https://www.gnu.org/licenses/gpl-3.0.en.html). See the [myfempy/LICENSE](https://github.com/easycae-3d/myfempy/blob/main/LICENSE.txt).

<!-- ## >> Acknowledgment -->

## Citing

Have you found this software useful for your research? Star the project and cite it as:

- APA:

```bash

Antonio Vinicius Garcia Campos. (2022). myfempy (1.5.1). Zenodo. https://doi.org/10.5281/zenodo.6958796

```

- BibTex:

```bash

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
  
## References

- [Myfempy](https://myfempy.readthedocs.io/) - *A python package for scientific analysis based on finite element method.*

- [FEM](https://en.wikipedia.org/wiki/Finite_element_method) - *The finite element method (FEM) is a popular method for numerically solving differential equations arising in engineering and mathematical modeling.*

- [Solid Mechanics](https://en.wikipedia.org/wiki/Solid_mechanics) - *Solid mechanics, also known as mechanics of solids, is the branch of continuum mechanics that studies the behavior of solid materials, especially their motion and deformation under the action of forces, temperature changes, phase changes, and other external or internal agents.*

- [PDE](https://en.wikipedia.org/wiki/Partial_differential_equation) - *In mathematics, a partial differential equation (PDE) is an equation which imposes relations between the various partial derivatives of a multivariable function.*

-----------

## Changelog

The changelog is available [here](https://github.com/easycae-3d/myfempy/wiki/Changelog)

# Project tree structure

```bash

/myfempy
│
├───core
│      alglin.py
│      harmonicforced.py
│      mesh.py
│      modallinear.py
│      solver.py
│      staticlinear.py
│      utils.py
│   
│
├───felib
│   │
│   ├───elements
│   │     element.py
│   │     plane.py
│   │     platekc.py
│   │     solid.py
│   │
│   ├───geometry
│   │     geometry.py
│   │     rectangle.py
│   │     thickness.py
│   │   
│   │
│   ├───material
│   │     material.py
│   │     planestrain.py
│   │     planestress.py
│   │     solid.py
│   │
│   ├───physic
│   │     bcstruct.py
│   │     loadstruct.py
│   │     structural.py
│   │
│   ├───shapes
│         hexa8.py
│         line.py
│         line2.py
│         quad4.py
│         shape.py
│         tetr4.py
│         tria3.py
│
├───io
│      iocsv.py
│      iogmsh.py
│      iovtk.py
│   
│
├───misc
│      logo.png
│      logo.txt
│      path.py
│      utils.py
│
├───plots
│     meshquality.py
│     physics.py
│     plotmesh.py
│     plotxy.py
│     postplot.py

```
