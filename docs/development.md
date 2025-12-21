# Development Stages

The main goal of this repository is to develop tests and configure a new version of the **myfempy** project (new version) based entirely on object-oriented programming (OOP).

This new version will allow for an expansion of the current **myfempy** version to high-performance computing (HPC), as well as a large data reading capacity, which is not currently possible due to its class complexity.

A future advancement will be to send data so that the solver executes the solution in parallel (multi-core) and also perform multi-physics analysis (Fluid-Structure interaction and Thermal-Structure interaction).

With the goal of keeping the project clean and clear throughout this development journey, the guidelines of Clean Architecture and the Bridge Design Pattern will be used for OOP code.

This new version allows the inclusion of modules and code written in C/Cython, as well as the possibility of using parallel processing for high-performance computing.

1. [ADDED] Study the characteristics of a FEM code with OOP

2. [ADDED] Extend the functionalities with OOP for the **myfempy** project

3. [ADDED] Generate a UML map of the classes in the code layers, making these maps available in the **User's Guide** documentation

4. [ADDED] Use the bridge design pattern to write the main classes of the system, as well as their feature() and method() functions

5. [ADDED] Develop an OOP code for the **myfempy** project

6. [ADDED] Develop an internal mesh generator using gmsh with .msh1 and .vtk mesh readers

7. [ADDED] Develop the output of results using .vtk files

8. [ ] Test the program to solve the multi-physics problem

9. [ADDED] Implementation of higher-order (quadritic) elements

10. [ADDED] Implementation of routines for multicore execution, e.g., Matrix assembly

11. [ADDED] Implement functionalities of the _myfempy/core_ code in C/Cython

12. [ ] Advanced solvers using the [PETSc](https://petsc.org/release/) and [SLEPc](https://slepc.upv.es/) packages

Implement the following analysis solutions:

1. Solvers:

	1. Steady State Linear

		- [ADDED] Direct
		- [ADDED] Iterative
		- [ADDED] Eigen (Modal)
		- [ADDED] Harmonic (FRF)
		- [ ] Harmonic (Modal)
		- [ADDED] Cyclic Symmetry
		- [ ] Buckling
		- [ADDED] Phononic Crystal 2D (Plane Wave Elastic in Periodic Micro Cell)

	2. Transient Response

		- [ ] Algo. Newmark

	3. Non-Linear

		- [ ] Algo. Newton-Raphson

2. Mechanical behavior of material:

	- [ADDED] Plane Stress
	- [ADDED] Plane Strain
	- [ADDED] Solid Elastic
	- [ADDED] Euler-Bernouilli Space Beam
	- [DEPRECATED] Timoshenko Space Beam
	- [ ] Plate Kirchhoff
	- [DEPRECATED] Plate Reissner-Mindlin
	- [ADDED] Homogenized Elastic Tensor: Heterogeneous Material Micro Base Cell
	- [ ] Large Displacement
	- [ ] Plasticity
	- [ADDED] Heat Plane
	- [ ] Fluid Flow Plane
	- [ ] Acustic Plane

3. Finite Elements Library:

	- [ADDED] line2 - linha 2 nós
	- [ADDED] line3 - linha 3 nós
	- [ADDED] tria3 - triagular 3 nós
	- [ADDED] tria6 - triagular 6 nós
	- [ADDED] quad4 - quadrilateral 4 nós
	- [ADDED] quad8 - quadrilateral 8 nós
	- [ADDED] tetr4 - tetraedro 4 nós
	- [DEPRECATED] tetr10 - tetraedro 10 nós
	- [ADDED] hexa8 - hexaedro 8 nós
	- [DEPRECATED] hexa20 - hexaedro 20 nós

4. Coupling:

	1. Multi-material

		- [ADDED] Multi Material Interface (multi E, v, ...)

	2. Multi-physical Analysis

		- [ADDED] Thermal-Structure interaction
		- [ ] Fluid-Structure interaction
		- [ ] Acoustic-Structure interaction

After the new implementation, a _merge_ was performed into the **myfempy_dev** repository (official repository for development and testing), and after passing all tests, it was uploaded to the project's main **myfempy** repository.
