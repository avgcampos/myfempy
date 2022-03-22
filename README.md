# myfempy project
![logo](https://user-images.githubusercontent.com/54820276/159380676-a35b17f3-e9cc-4ee4-b8b1-0a7c4c3fb0d6.svg)


[![lics](https://img.shields.io/badge/license-GPL-blue.svg)](https://en.wikipedia.org/wiki/GNU_General_Public_License)

[![DOI](https://zenodo.org/badge/462762513.svg)](https://zenodo.org/badge/latestdoi/462762513)

<!--

[![Downloads](https://pepy.tech/badge/vedo)](https://pepy.tech/project/vedo)

[![CircleCI](https://circleci.com/gh/marcomusy/vedo.svg?style=svg)](https://circleci.com/gh/marcomusy/vedo)
 
-->


## About
**Myfempy** is a python package based on finite element method for scientific analysis. The code is open source and *intended for educational and scientific purposes only, not recommended to commercial use*. You can help us by contributing with a donation on the main project page, read the support options. **If you use myfempy in your research, the  developers would be grateful if you could cite in your work.**

**Myfempy** uses command lines to solve simulations.

- Author: Antonio Vinicius Garcia Campos [GitHub page](https://github.com/antonio-vinicius-garcia-campos)
- Version: 1.0.0 dev/ mar 2022

Copyright © Antonio Vinicius G. Campos and 3D EasyCAE, 2022

## Installation
### To install myfempy manually in your directory, following the steps

1. Download the main code from [github/myfempy/main](https://github.com/easycae-3d/myfempy/tree/main)

2. Unzip the pack in your preferred location

3. In the **myfempy-main** folder, open a terminal and enter with the command:

```bash
pip install .
```

*Note: is recommend to create a new virtual environment previously the installation of **myfempy** and dependencies packs. You can use the [virtualenv](https://virtualenv.pypa.io/en/latest/)* 

## Requirements

**Myfempy** can be used in systems based on Linux, MacOS and Windows. **Myfempy** requires Python 3.

### Installation prerequisites, required to build **myfempy**:
- [Python 3.x](https://www.python.org/) - *Python is a programming language that lets you work quickly and integrate systems more effectively.*
- [Anaconda](https://www.anaconda.com/) - *Anaconda offers the easiest way to perform Python/R data science and machine learning on a single machine.*
------------
#### Outhers prerequisites
- [Gmsh](https://gmsh.info/) - Gmsh is an open source 3D finite element mesh generator with a built-in CAD engine and post-processor. *Note: needed install manually*

### Python packages required for using **myfempy**:
- [numpy](https://numpy.org/) - The fundamental package for scientific computing with Python
- [scipy](https://scipy.org/) - Fundamental algorithms for scientific computing in Python
- [vedo](https://vedo.embl.es/) - A python module for scientific analysis and visualization of эd objects


```bash
pip install numpy, scipy, vedo
```

## Documentation
The main documentation is available here: [myfempy web doc](https://myfempy.readthedocs.io/).

The myfempy GitHub page/download page is available here: [myfempy github page](https://github.com/easycae-3d/myfempy/).

The **User's Manual** is available here: [myfempy/docs](https://github.com/easycae-3d/myfempy/blob/master/docs/Users_Manual.pdf).

Many examples are available here: [myfempy/examples](https://github.com/easycae-3d/myfempy/tree/master/examples).

## Release

The version up to date is **myfempy v1.0.0 dev**. Go to *Features List/Version History* to visualization all versions releses.

## Features

[Features List](https://docs.google.com/spreadsheets/d/1k9kiXk2PPuUvcsiukAni005zQc-IOCmP2r-Z6B02304/edit?usp=sharing)


## License

**myfempy** is published under the [GPLv3 license](https://en.wikipedia.org/wiki/GNU_General_Public_License)

<!-- ## >> Acknowledgment -->

## Citing

Have you found this software useful for your research? Star the project and cite it as:

Vinicius Campos. (2022). easycae-3d/myfempy: beta (v1.0.1). Zenodo. https://doi.org/10.5281/zenodo.6376522


## References

- [Myfempy](https://myfempy.readthedocs.io/) - *A python package for scientific analysis based on finite element method.* 

- [FEM](https://en.wikipedia.org/wiki/Finite_element_method) - *The finite element method (FEM) is a popular method for numerically solving differential equations arising in engineering and mathematical modeling.*

- [Solid Mechanics](https://en.wikipedia.org/wiki/Solid_mechanics) - *Solid mechanics, also known as mechanics of solids, is the branch of continuum mechanics that studies the behavior of solid materials, especially their motion and deformation under the action of forces, temperature changes, phase changes, and other external or internal agents.*

- [PDE](https://en.wikipedia.org/wiki/Partial_differential_equation) - *In mathematics, a partial differential equation (PDE) is an equation which imposes relations between the various partial derivatives of a multivariable function.*

-----------
# project tree structure
```bash
/myfempy
|--/bin
|	gui.py
|	plotter.py
|
|
|
|--/felib
|	|--/fluid
|	|
|	|
|	|
|	|--/fsi
|	|
|	|
|	|
|	|--/struct
|		beam21.py
|		frame22.py
|		frame23.py
|		plane32.py
|		plane42.py
|		solid83.py
|		spring20.py
|		truss22.py
|		.py
|		.py
|		.py
|
|	integrat.py
|	material.py
|	postproc.py
|
|
|
|--/help
|	help.py
|	version.py
|
|
|
|--/io
|	filters.py
|	ioctrl.py
|	miscel.py
|
|
|
|--/mesh
|	meshgen.py
|
|
|
|--/solver
|		assembly.py
|		bcloads.py
|		plotter.py
|		solverset.py
|		static.py
|		vibra.py
|
```
