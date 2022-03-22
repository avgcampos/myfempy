.. myfempy documentation master file, created by
   sphinx-quickstart on Mon Mar 14 21:53:28 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.
   
.. only:: html

Welcome to myfempy's web documentation!
=======================================

About
----------------

**Myfempy** is a python package based on finite element method for scientific analysis. The code is open source and *intended for educational and scientific purposes only, not recommended to commercial use*. You can help us by contributing with a donation on the main project page, read the support options. **If you use myfempy in your research, the  developers would be grateful if you could cite in your work.**

**Myfempy** uses command lines to solve simulations.

- Author: Vinicius Campos `GitHub page <https://github.com/antonio-vinicius-garcia-campos>`_
- Version: 1.0.0 dev/ mar 2022

Copyright © Antonio Vinicius G. Campos and 3D EasyCAE, 2022

Installation
----------------

	To install myfempy manually in your directory, following the steps

	1. Download the main code from `github/myfempy/main <https://github.com/easycae-3d/myfempy/tree/main>`_

	2. Unzip the pack in your preferred location

	3. In the **myfempy-main** folder, open a terminal and enter with the command:

	- try: ::
	
		pip install .
	

Note: is recommend to create a new virtual environment previously the installation of **myfempy** and dependencies packs. You can use the `virtualenv <https://virtualenv.pypa.io/en/latest/>`_  

Requirements
----------------

**Myfempy** can be used in systems based on Linux, MacOS and Windows. **Myfempy** requires Python 3.

Installation prerequisites, required to build **myfempy**:

- `Python 3.x <https://www.python.org/>`_ - *Python is a programming language that lets you work quickly and integrate systems more effectively.*

- `Anaconda <https://www.anaconda.com/>`_ - *Anaconda offers the easiest way to perform Python/R data science and machine learning on a single machine.*
------------
Outhers prerequisites

- `Gmsh <https://gmsh.info/>`_ - Gmsh is an open source 3D finite element mesh generator with a built-in CAD engine and post-processor. *Note: needed install manually*

Python packages required for using **myfempy**:

- `numpy <https://numpy.org/>`_ - The fundamental package for scientific computing with Python

- `scipy <https://scipy.org/>`_ - Fundamental algorithms for scientific computing in Python

- `vedo <https://vedo.embl.es/>`_ - A python module for scientific analysis and visualization of эd objects

- try: ::

	pip install numpy, scipy, vedo


Documentation
----------------

The main documentation is available here: `myfempy web doc <https://myfempy.readthedocs.io/>`_

The myfempy GitHub page/download page is available here: `myfempy github page <https://github.com/easycae-3d/myfempy/>`_

The **User's Manual** is available here: `myfempy/docs <https://github.com/easycae-3d/myfempy/blob/main/docs/Users_Manual.pdf>`_

Many examples are available here: `myfempy/examples <https://github.com/easycae-3d/myfempy/tree/main/examples>`_

Release
----------------

The version up to date is **myfempy v1.0.0 dev**. Go to *Features List/Version History* to visualization all versions releses.

Features
----------------

`Features List <https://docs.google.com/spreadsheets/d/1k9kiXk2PPuUvcsiukAni005zQc-IOCmP2r-Z6B02304/edit?usp=sharing>`_


License
----------------

**myfempy** is published under the `GPLv3 license <https://en.wikipedia.org/wiki/GNU_General_Public_License>`_

Citing
----------------

Have you found this software useful for your research? Star the project and cite it as:

- BibTeX::

	@software{campos_2022,
	  author       = {Antonio Vinicius Garcia Campos,
					  },
	  title        = {myfempy: 1.0.0 dev},
	  month        = mar,
	  year         = 2022,
	  publisher    = {Zenodo},
	  version      = {},
	  doi          = {},
	  url          = {}
	}

References
----------------

- `Myfempy <https://myfempy.readthedocs.io/>`_ - *A python package for scientific analysis based on finite element method.* 

- `FEM <https://en.wikipedia.org/wiki/Finite_element_method>`_ - *The finite element method (FEM) is a popular method for numerically solving differential equations arising in engineering and mathematical modeling.*

- `Solid Mechanics <https://en.wikipedia.org/wiki/Solid_mechanics>`_ - *Solid mechanics, also known as mechanics of solids, is the branch of continuum mechanics that studies the behavior of solid materials, especially their motion and deformation under the action of forces, temperature changes, phase changes, and other external or internal agents.*

- `PDE <https://en.wikipedia.org/wiki/Partial_differential_equation>`_ - *In mathematics, a partial differential equation (PDE) is an equation which imposes relations between the various partial derivatives of a multivariable function.*


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. toctree::
   :hidden:
   :maxdepth: 3
   
   About <self>
   documentation
   tutorial