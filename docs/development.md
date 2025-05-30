# Etapas de desenvolvimento

Este repositório tem como objetivo principal desenvolver testes e configurar uma nova versão do projeto **myfempy** (new version) baseando-se totalmente em programação orientada a objetos (OOP).
Esta nova versão permitirá uma expansão da atual versão do **myfempy** para alta performance computacional (HPC), além de uma capacidade grande de leitura de dados, o que não é possível atualmente devido a sua complexidade de classes.
Um avanço futuro, será de enviar dados para que o solver execute a solução em paralelo (multi-core) e também realizar a análise em multi-física ([Fluid-Structure interaction](https://en.wikipedia.org/wiki/Fluid%E2%80%93structure_interaction) 
e [Thermal-Structure interaction](https://en.wikipedia.org/wiki/Thermal_stress)).

Com o objetivo de deixar o projeto limpo e claro nesta jornada de desenvolvimento, será utilizado as normativas da _clean archicteture_ e da _bridge design pattern_ para códigos com OOP.

Esta nova versão permite a inclusão de módulos e códigos escritos em C/Cython além da possibilidade de utilizar o processamento em paralelo para a computação de alta performance.

1. [X] Estudar as caracteristicas de um código FEM com OOP

2. [X] Extender as funcionalidades com OOP para o projeto **myfempy**

3. [X] Gerar um mapa UML das clases nas layers do código, disponibilizando estes mapas na documentação do **User's Guide**

4. [X] Usar _bridge design pattern_ para escrerver as principais classes do sistema, assim como as suas _feacture()_ e _method()_

5. [X] Desenvolver um código OOP para o projeto **myfempy**

6. [X] Desenvolver um gerador de malha interno utilizando o gmsh com leitor de malha .msh1 e .vtk

7. [X] Desenvolver a saida dos resultados por meio de arquivos .vtk

8. [X] Testar o programa para resolver o problema multi-físicos

9. [X] Implementação de elementos de alta ordem (quadritico)

10. [X] Implementação de rotinas para execução em multicore, exe: Montagem das matrizes

11. [X] Implementar funcionalidades do código do _core_ em C/Cython

12. [ ] Solvers avançados utilizando os pacotes [PETSc](https://petsc.org/release/) e [SLEPc](https://slepc.upv.es/)

Implementar as seguintes soluções de análise:

	1. Solver:

		1. Permanente linear

			- [X] Direct
			- [X] Iterative
			- [X] Dynamic Modal (Eigen)
			- [X] Dynamic FRF (Direct)
			- [ ] Dynamic FRF (Modal)
            - [X] Cyclic Symmetry
			- [ ] Buckling (Eigen)
			- [X] Phononic Crystal 2D (Plane Wave Elastic in Periodic Micro Cell)

		2. Transiênte linear

			- [ ] Algo. Newmark

		3. Non Linear

			- [ ] Algo. Newton Raphson

    2. Comportamento mecânico de material:

        - [X] Plane Stress
        - [X] Plane Strain
        - [X] Solid Elastic
		- [X] Euler-Bernouilli Space Beam
		- [ ] Timoshenko Space Beam
        - [ ] Plate Kirchhoff
        - [ ] Plate Reissner-Mindlin
        - [X] Homogenized Elastic Tensor/ Heterogeneous Material Micro Base Cell
		- [ ] Grande Deslocamento
		- [ ] Plasticidade
		- [X] Heat Plane
		- [ ] Fluid Flow Plane
		- [ ] Acustic Plane

	3. Elementos isoparametricos:

		- [X] line2 - linha 2 nós
		- [X] line3 - linha 3 nós
		- [X] tria3 - triagular 3 nós
		- [X] tria6 - triagular 6 nós
		- [X] quad4 - quadrilateral 4 nós
		- [X] quad8 - quadrilateral 8 nós
		- [X] tetr4 - tetraedro 4 nós
		- [ ] tetr10 - tetraedro 10 nós
		- [X] hexa8 - hexaedro 8 nós
		- [ ] hexa20 - hexaedro 20 nós

	4. Acoplamento MMulti-físico:

		1. Multi-material

			- [X] Interface de multi-material (multi E, v, ...)

		2. Acoplamento multi-físico

			- [X] Thermal-Structure interaction
			- [ ] Fluid-Structure interaction
			- [ ] Acoustic-Structure interaction

Após a nova implementação, foi feito um _merge_ no repositório **myfempy_dev** (repositório oficial para desenvolvimento e testes), e após passar por todos os teste foi feito o upload para
o repositório principal **myfempy** do projeto.
