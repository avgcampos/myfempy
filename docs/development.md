# Etapas de desenvolvimento

Este repositório tem como objetivo principal desenvolver testes e configurar uma nova versão do projeto **myfempy** (new version) baseando-se totalmente em programação orientada a objetos (OOP). Esta nova versão permitirá uma expansão da atual versão do **myfempy** para alta performance computacional (HPC), além de uma capacidade grande de leitura de dados, o que não é possível atualmente devido a sua complexidade de classes. Um avanço futuro, será de enviar dados para que o solver execute a solução em paralelo (multi-core) e também realizar a análise em multi-física ([Fluid-Structure interaction](https://en.wikipedia.org/wiki/Fluidstructure_interaction) e [Thermal-Structure interaction](https://en.wikipedia.org/wiki/Thermal_stress)).

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

10. [ ] Implementação de rotinas para execução em multicore, exe: Montagem das matrizes

11. [ ] Solvers avançados utilizando os pacotes [PETSc](https://petsc.org/release/) e [SLEPc](https://slepc.upv.es/)

Implementar as seguintes soluções de análise:

	1. Estrutural estático:
		1. Linear
			- [X] Elastico
		2. Non-linear
			- [ ] Grande deslocamento
			- [ ] Plasticidade
			- [ ] Estabilidade (flambagem)
	2. Estrutural dinamica linear:
		1. Domínio frequência
			- [X] Modal
			- [X] Força harmônica (direto)
			- [ ] Força harmônica (modal)
		2. Transiênte linear
			- [ ] Algo. Newmark
	3. Domínio termico e fluido:
		1. Condução
			- [X] Regime permanente
			- [ ] Regime transiente
		2. Escoamento
			- [ ] Laminar permanente
			- [ ] Laminar transiente
		3. Acústico
			- [ ] Modal
	4. Multi-domínio:
		1. Multi-material
			- [X] Interface de multi-material (multi E, v, ...)
		2. Acoplamento multi-físico
			- [ ] Fluid-Structure interaction
			- [X] Thermal-Structure interaction
			- [ ] Acoustic-Structure interaction
    5. Comportamento mecânico e material:
        - [X] Plane stress
        - [X] Plane strain
        - [ ] Solid isotropic
        - [ ] Plate Kirchhoff
        - [ ] Plate Reissner-Mindlin
        - [X] Cyclic symmetry boundary
		- [X] Pressure load application
        - [ ] Homogeinização/ micro escala (tensor)

Desenvolver os seguintes tipos de elementos isoparametrico:

    - [ ] line2 - linha 2 nós
    - [ ] line3 - linha 3 nós
    - [X] tria3 - triagular 3 nós
    - [X] tria6 - triagular 6 nós
    - [X] quad4 - quadrilateral 4 nós
    - [X] quad8 - quadrilateral 8 nós
    - [ ] tetr4 - tetraedro 4 nós
    - [ ] tetr10 - tetraedro 10 nós
    - [ ] hexa8 - hexaedro 8 nós
    - [ ] hexa20 - hexaedro 20 nós



Após a nova implementação, foi feito um _merge_ no repositório **myfempy_dev** (repositório oficial para desenvolvimento e testes), e após passar por todos os teste foi feito o upload para o repositório principal **myfempy** do projeto.