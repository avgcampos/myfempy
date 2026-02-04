# Welcome to myfempy's online documentation

**Under Development**

The **myfempy** project is under development, updates and code modifications may occur in future versions without prior notice from the developers.

![myfempy_logo](assets\logo2.png)

Copyright © Antonio Vinicius G. Campos 2022. Processo INPI BR512022001484-0

## About

**myfempy** is a python package based on finite element method to multiphysics analysis. The code is open source and *intended for educational and scientific purposes only, not recommended to commercial use*. The name **myfempy** is an acronym for **M**ultiph**Y**sics **F**inite **E**lements **M**odule to **PY**thon. You can help us by contributing with the main project, send us a mensage [Github Discussions](https://github.com/avgcampos/myfempy/discussions/10). **If you use myfempy in your research, the  developers would be grateful if you could cite in your work.**

![home](assets\sim.png)

## Installation

To install myfempy manually in your directory, following the steps

1.  Clone/ Download the main code latest version from
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
installation of myfempy and dependencies packs. You can use the
[virtualenv](https://virtualenv.pypa.io/en/latest/) or [conda
environments](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html)**

## Dependencies

**Myfempy** can be used in systems based on Linux, MacOS and Windows, and requires Python 3 to running.

Installation prerequisites, required to build **myfempy**

You can use either of two python development environments to run myfempy

-   [Python 3.x](https://www.python.org/) - *Python is a programming
    language that lets you work quickly and integrate systems more
    effectively.*
-   [Anaconda](https://www.anaconda.com/) - *Anaconda offers the easiest
    way to perform Python/R data science and machine learning on a
    single machine.*

Basic python packages required for using **myfempy**

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
-   [vtk](https://pypi.org/project/vtk/) - VTK is an
    open-source toolkit for 3D computer graphics, image processing, and
    visualization

try

``` bash
>> pip install numpy, cython, scipy, vedo
```

Outhers prerequisites

-   [gmsh/External Generator Mesh](https://gmsh.info/) - Gmsh is an open
    source 3D finite element mesh generator with a built-in CAD engine
    and post-processor. *Notes: 1 - Gmsh is NOT part of myfempy
    projects; 2 - Is Needed install Gmsh manually*

try

``` bash
>> pip install --upgrade gmsh
```

-   [gmsh PyPi](https://pypi.org/project/gmsh/)

## Tutorial

An **User's Guide** is available
[here](user_manual.md).

Many **Examples** are available
[here](https://github.com/avgcampos/myfempy/tree/main/examples).

## Documentation

The myfempy's web documents can be found [here](https://myfempy.readthedocs.io/en/latest/). The myfempy is documented using the [Mkdocs](https://www.mkdocs.org/) under **docs** folder. 

## Releases/ Versions

The changelog is available
[here](https://github.com/avgcampos/myfempy/releases)

## License

**myfempy** is published under the [GPLv3
license](https://www.gnu.org/licenses/gpl-3.0.en.html). See the
[myfempy/LICENSE](https://github.com/avgcampos/myfempy/blob/main/LICENSE).

## Citing

Have you found this software useful for your research? Star the project
and cite it as:

-   APA:

``` bash
Campos, A. V. G. (2022). myfempy. Zenodo. https://doi.org/10.5281/zenodo.15756128
```

-   BibTex:

``` bash
@software{campos_2022_15756128,
  author       = {Campos, Antonio Vinicius Garcia},
  title        = {myfempy},
  month        = jul,
  year         = 2022,
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.15756128},
  url          = {https://doi.org/10.5281/zenodo.15756128},
}
```

## References

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
    *In mathematics, a partial differential equation (PDE) is an
    equation which imposes relations between the various partial
    derivatives of a multivariable function.*