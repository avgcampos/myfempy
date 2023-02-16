.. myfempy documentation master file, created by
   sphinx-quickstart on Mon Mar 14 21:53:28 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.
   
.. only:: html

**Under Development**

.. figure:: logoB.png

Welcome to myfempy's web documentation!
=======================================

Copyright © Antonio Vinicius G. Campos and 3D EasyCAE, 2022

About
=====

**Myfempy** is a python package based on finite element method for
scientific analysis. The code is open source and *intended for
educational and scientific purposes only, not recommended to commercial
use*. You can help us by contributing with a donation on the main
project page, read the support options. **If you use myfempy in your
research, the developers would be grateful if you could cite in your
work.**

Installation
============

To install myfempy manually in your directory, following the steps
------------------------------------------------------------------

1. Clone/ Download the main code [latest version] from
   `github/myfempy/main <https://github.com/easycae-3d/myfempy/>`__

2. Unzip the pack in your preferred location

3. In the **myfempy-main** folder, open a terminal and enter with the
   command:

.. code:: bash


   >> python -m pip install --upgrade pip

   >> pip install .

   or

   >> python -m pip install --upgrade build

   >> python -m build

**Note: is recommend to create a virtual environment previously the
installation of** myfempy*\* and dependencies packs. You can use the
`virtualenv <https://virtualenv.pypa.io/en/latest/>`__ or `conda
environments <https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html>`__\ \*\*

Dependencies
============

**Myfempy** can be used in systems based on Linux, MacOS and Windows.
**Myfempy** requires Python 3.

Installation prerequisites, required to build **myfempy**
---------------------------------------------------------

You can use either of two python development environments to run myfempy

-  `Python 3.x <https://www.python.org/>`__ - *Python is a programming
   language that lets you work quickly and integrate systems more
   effectively.*
-  `Anaconda <https://www.anaconda.com/>`__ - *Anaconda offers the
   easiest way to perform Python/R data science and machine learning on
   a single machine.*

Python packages required for using **myfempy**
----------------------------------------------

The following python packages are required to run myfempy. Before to
install myfempy-main, install this packages. Check if they are already
installed on your machine

-  `numpy <https://numpy.org/>`__ - The fundamental package for
   scientific computing with Python

-  `cython <https://cython.org/>`__ - Cython is a language that makes
   writing C extensions for Python as easy as Python itself

-  `scipy <https://scipy.org/>`__ - Fundamental algorithms for
   scientific computing in Python

-  `vedo <https://vedo.embl.es/>`__ - A python module for scientific
   analysis and visualization of эd objects

-  `vtk <https://pypi.org/project/vtk/>`__\ (optional) - VTK is an
   open-source toolkit for 3D computer graphics, image processing, and
   visualization

-  try

.. code:: bash


   >> pip install numpy, cython, scipy, vedo

Outhers prerequisites
~~~~~~~~~~~~~~~~~~~~~

-  `gmsh/External Generator Mesh <https://gmsh.info/>`__ - Gmsh is an
   open source 3D finite element mesh generator with a built-in CAD
   engine and post-processor. *Notes: 1 - Gmsh is NOT part of myfempy
   projects; 2 - Is Needed install Gmsh manually*

-  try

.. code:: bash


   >> pip install --upgrade gmsh

-  `gmsh PyPi <https://pypi.org/project/gmsh/>`__

Tutorial
========

A **Basic Tutorial** is available
`here <https://myfempy.readthedocs.io/en/1.dev9/tutorial.html>`__.

Many **Examples** are available
`here <https://github.com/easycae-3d/myfempy/tree/master/examples>`__.

Documentation
=============

*The myfempy is documented using Sphinx under ``docs/``. The downloads
built versions can be found in pdf, html or epub.*

The **Web Documentation** is available on `Read the
Docs <https://myfempy.readthedocs.io/>`__.

The **User’s Manual [pdf]** is available
`here <https://github.com/easycae-3d/myfempy/blob/main/docs/Users%20Manual.pdf>`__.

Release
=======

The all release versions is available
`here <https://github.com/easycae-3d/myfempy/releases>`__

Features
========

The *main myfempy features* are available here:

-  `Features
   List <https://github.com/easycae-3d/myfempy/blob/main/docs/Myfempy%20Features.pdf>`__

License
=======

**myfempy** is published under the `GPLv3
license <https://www.gnu.org/licenses/gpl-3.0.en.html>`__. See the
`myfempy/LICENSE <https://github.com/easycae-3d/myfempy/blob/main/LICENSE.txt>`__.

.. raw:: html

   <!-- ## >> Acknowledgment -->

Citing
======

Have you found this software useful for your research? Star the project
and cite it as:

-  APA:

.. code:: bash


   Antonio Vinicius Garcia Campos. (2022). myfempy (1.5.1). Zenodo. https://doi.org/10.5281/zenodo.6958796

-  BibTex:

.. code:: bash


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

References
==========

-  `Myfempy <https://myfempy.readthedocs.io/>`__ - *A python package for
   scientific analysis based on finite element method.*

-  `FEM <https://en.wikipedia.org/wiki/Finite_element_method>`__ - *The
   finite element method (FEM) is a popular method for numerically
   solving differential equations arising in engineering and
   mathematical modeling.*

-  `Solid Mechanics <https://en.wikipedia.org/wiki/Solid_mechanics>`__ -
   *Solid mechanics, also known as mechanics of solids, is the branch of
   continuum mechanics that studies the behavior of solid materials,
   especially their motion and deformation under the action of forces,
   temperature changes, phase changes, and other external or internal
   agents.*

-  `PDE <https://en.wikipedia.org/wiki/Partial_differential_equation>`__
   - *In mathematics, a partial differential equation (PDE) is an
   equation which imposes relations between the various partial
   derivatives of a multivariable function.*

--------------

Changelog
=========

The changelog is available
`here <https://github.com/easycae-3d/myfempy/wiki/Changelog>`__

Project tree structure
======================

.. code:: bash


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

--------------

.. Indices and tables
.. ==================

.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`

.. toctree::
   :hidden:
   :maxdepth: 2
   
   Quick Start <self>
   tutorial
   documentation
   Source Code <https://github.com/easycae-3d/myfempy>