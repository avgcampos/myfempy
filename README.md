![logo_v1](https://user-images.githubusercontent.com/54820276/159730160-871ac41a-958a-4398-a014-506619c4cb56.png)

Copyright © Antonio Vinicius G. Campos and 3D EasyCAE, 2022

|  | Badges |
| --- | :---: |
| Documentation | [![Documentation Status](https://readthedocs.org/projects/myfempy/badge/?version=latest)](https://myfempy.readthedocs.io/en/latest/?badge=latest) |
| Packages | [![PyPI]()]() [![conda]()]() |
| Paper | [![DOI](https://zenodo.org/badge/462762513.svg)](https://zenodo.org/badge/latestdoi/462762513) |
| Code quality | [![Quality Gate Status](https://sonarcloud.io/api/project_badges/measure?project=easycae-3d_myfempy&metric=alert_status)](https://sonarcloud.io/summary/new_code?id=easycae-3d_myfempy) |
| Downloads | [![Downloads]() |
| License | [![lics](https://img.shields.io/badge/license-GPL-blue.svg)](https://en.wikipedia.org/wiki/GNU_General_Public_License) |

-----------

## About
**Myfempy** is a python package based on finite element method for scientific analysis. The code is open source and *intended for educational and scientific purposes only, not recommended to commercial use*. You can help us by contributing with a donation on the main project page, read the support options. **If you use myfempy in your research, the  developers would be grateful if you could cite in your work.**


## Installation
### To install myfempy manually in your directory, following the steps

1. Download the main code from [github/myfempy/main](https://github.com/easycae-3d/myfempy/tree/main)

2. Unzip the pack in your preferred location

3. In the **myfempy-main** folder, open a terminal and enter with the command:

```bash
pip install .
```

*Note: is recommend to create a new virtual environment previously the installation of **myfempy** and dependencies packs. You can use the [virtualenv](https://virtualenv.pypa.io/en/latest/)* 

## Dependencies


**Myfempy** can be used in systems based on Linux, MacOS and Windows. **Myfempy** requires Python 3.


### Installation prerequisites, required to build **myfempy**:
- [Python 3.x](https://www.python.org/) - *Python is a programming language that lets you work quickly and integrate systems more effectively.*
- [Anaconda](https://www.anaconda.com/) - *Anaconda offers the easiest way to perform Python/R data science and machine learning on a single machine.*


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

The project is documented using Sphinx under docs/. Built version can be found from Read the Docs. Here are direct links to additional resources:

The main documentation is available [here](https://myfempy.readthedocs.io/).

The myfempy GitHub page/download page is available [here](https://github.com/easycae-3d/myfempy/).

The **User's Manual [PT-BR]** is available [here](https://github.com/easycae-3d/myfempy/blob/master/docs/Users_Manual.pdf).

Many examples are available [here](https://github.com/easycae-3d/myfempy/tree/master/examples).

## Release

The version up to date is available here: 
[myfempy releases](https://github.com/easycae-3d/myfempy/releases)

Go to *Features List/Version History* to visualization all versions releses.

## Features

[Features List](https://docs.google.com/spreadsheets/d/1k9kiXk2PPuUvcsiukAni005zQc-IOCmP2r-Z6B02304/edit?usp=sharing)


## License

**myfempy** is published under the [GPLv3 license](https://en.wikipedia.org/wiki/GNU_General_Public_License)

<!-- ## >> Acknowledgment -->

## Citing

Have you found this software useful for your research? Star the project and cite it as:

- APA:
  ```
  Antonio Vinicius Garcia Campos. (2022). easycae-3d/myfempy: beta (v1.0.1). Zenodo. https://doi.org/10.5281/zenodo.6376522
  ```

- BibTex:
  ```bash
   @software{Campos_easycae-3d_myfempy_beta_2022,
            author = {Antonio Vinicius Garcia Campos},
            title = {easycae-3d/myfempy: beta},
            version = {v1.0.1},
            url = {https://github.com/easycae-3d/myfempy/},
            doi = {10.5281/zenodo.6376522},
            month = {3},
            year = {2022}
            }
  ```
  
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
