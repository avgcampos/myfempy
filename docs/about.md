# mymfempy

```bash
#                     __                                
#  _ __ ___   _   _  / _|  ___  _ __ ___   _ __   _   _ 
# | '_ ` _ \ | | | || |_  / _ \| '_ ` _ \ | '_ \ | | | |
# | | | | | || |_| ||  _||  __/| | | | | || |_) || |_| |
# |_| |_| |_| \__, ||_|   \___||_| |_| |_|| .__/  \__, |
#             |___/                       |_|     |___/ 

# Copyright © Antonio Vinicius G. Campos, 2023. All rights reserved. 

```

## Motivação

Este repositório tem como objetivo principal desenvolver testes e configurar uma nova versão do projeto **myfempy** (new version) baseando-se totalmente em programação orientada a objetos (OOP). Esta nova versão permitirá uma expansão da atual versão do **myfempy** para alta performance computacional (HPC), além de uma capacidade grande de leitura de dados, o que não é possível atualmente devido a sua complexidade de classes. Um avanço futuro, será de enviar dados para que o solver execute a solução em paralelo (multi-core) e também realizar a análise em multi-física (FSI, FTI, TMI).

Esta nova versão permite a inclusão de módulos e códigos escritos em C/Cython e os solvers avançados utilizando os pacotes [PETSc](https://petsc.org/release/) e [SLEPc](https://slepc.upv.es/), além da possibilidade de utilizar o [JAX](https://jax.readthedocs.io/en/latest/index.html) para computação de alta performance.

Após a nova implementação, foi feito um _merge_ no repositório **myfempy_dev** (repositório oficial para desenvolvimento e testes), e após passar por todos os teste será feito o upload para o repositório principal **myfempy** do projeto.

Com o objetivo de deixar o projeto limpo e claro nesta jornada de desenvolvimento, será utilizado as normativas da _clean archicteture_ e da _bridge design pattern_ para códigos com OOP.

## Etapas de desenvolvimento

1. Estudar as caracteristicas de um código FEM com OOP;

2. Extender as funcionalidades com OOP para o projeto **myfempy**;

3. Gerar um mapa UML das clases nas layers do código, disponibilizando estes mapas na documentação do **User's Manual**;

4. Usar _bridge design pattern_ para escrerver as principais classes do sistema, assim como as suas _feacture()_ e _method()_;

5. Desenvolver um código OOP para o projeto **myfempy**;

6. Testar o programa para resolver o problema de elasticidade bidimensinal de placa (Modelo de Mindlin)/ casca (Plate+PlaneStress) retangular;

7. Implementar as seguintes soluções de análise:
	
	1. Estático:
		1. Linear
			- [X] Elastico
		2. Non-linear
			- [ ] Grande deslocamento
			- [ ] Plasticidade
			- [ ] Estabilidade (flambagem)
			- [ ] Contato sem atrito??
	2. Dinamica linear:
		1. Domínio frequência
			- [X] Modal
			- [X] Força harmônica
		2. Transiênte linear
			- [ ] Algo. implicito
			- [ ] Algo. explicito?
	3. Multi-domínio:
		1. Multi-material 	
			- [X] Interface de multi-material (multi E, v, ...)
		2. Acoplamento multi-físico
			- [ ] FSI
			- [ ] TMI
			- [ ] Kratos/ openfoam?
    4. Comportamento mecânico e material:
        - [X] Plane stress
        - [X] Plane strain
        - [X] Solid isotropic
        - [ ] Plate Kirchhoff
        - [ ] Plate Reissner-Mindlin
        - [ ] Cyclic symmetry boundary
        - [ ] Homogeinização/ micro escala (tensor)

8. Validar a análise com a implementação seguintes tipos de elementos isoparametrico:

    - [X] tria3 - triagular 3 nós
    - [ ] tria6 - triagular 6 nós
    - [X] quad4 - quadrilateral 4 nós
    - [ ] quad8 - quadrilateral 8 nós
    - [X] tetr4 - tetraedro 4 nós
    - [ ] tetr10 - tetraedro 10 nós
    - [X] hexa8 - hexaedro 8 nós
    - [ ] hexa20 - hexaedro 20 nós

9. Desenvolver um gerador de malha interno utilizando o gmsh (pygmsh?) com leitor de malha .msh1 e .vtk;

10. Desenvolver a saida dos resultados por meio de arquivos .vtk;

# Regras para realizar upload {*commit*} no projeto myfempy

Este repositório destina-se ao gerenciamento das versões de desenvolvimento do projeto **myfempy**. Para que se tenha um código claro e limpo, as seguintes regras são aplicadas para manter um padrão geral do projeto.

# 1 - Estilos gerais utilizados no projeto

## 1.1 - Estilo dos códigos *.py*

O padrão utilizado na sintax da escrita dos códigos é o [PEP8](https://peps.python.org/pep-0008).

### 1.1.1 - Padrão de criação de classes *class AbcAbc*

```python

# the file name is abcabc.py

class AbcAbc: # class name

    def __init__(self, x:int, y:float, z:str):  # constructor method in python

        # add parameter in class
        self.x = x
        self.y = y
        self.z = z

```

### 1.1.2 - Padrão de criação de funções *def abc_abc*

```python

# the file name is abcabc.py

def abc_abc(x:int, y:float): # function

    soma_x = x + x
    div_y = x/y

    return div_y

```

## 1.2 - ```__doc__``` de classes e funções

Toda classe e função deve conter sua documentação interna {```__doc__```}.

### 1.2.1 - Procedimento

1. Utilizar a referência [PEP8](https://peps.python.org/pep-0008).

2. Apresentar em uma única linha o que aquela classe ou função executa.

3. Use comentarios para informações adicionais.

4. Todas as entradas e saídas e a tipagem de dados deve ser indicados.

5. O idioma padrão é o inglês técnico.

Exemplo:

```python

# the file name is abcabc.py
# '__doc__' to a short documentation, try AbcAbc.__doc__ ....

__doc__ = """ 
This is a short documentation of class AbcAbc...
"""
class AbcAbc: # class name
    """
    Abc: class

    This class do...
    """    

    def __init__(self, x:int, y:float, z:str):  # constructor method in python
        """
        This class do...

        Args:
            x (int): int number
            y (float): float number
            z (str): string name
        """        

        # add parameter in class
        self.x = x
        self.y = y
        self.z = z

    def abc_abc(self): # class function

        """
        This function do...

        Returns:
            div_y (float): _description_
        """        

        soma_x = self.x + self.x
        div_y = self.x/soma_x

        return div_y
```

## 1.3 - Documentação online (Read The Docs)

O padrão utilizado para escrever a documentação oficial do projeto é no formato [Markdown](https://www.markdownguide.org/basic-syntax).

Atualmente o projeto utiliza o pacote [mkdocs](https://www.mkdocs.org/) para gerar a documentação online {*.html*}, que é hospedada no site [GitHub Page](https://avgcampos.github.io/myfempy/).

O *sphinx* converte arquivos *.rst* {reStruturedText} para *.html/pdf* de uma forma prática. Também é possível converter arquivos *.md* {markdown} em *.rst*, uma opção interessante é [mdToRst](https://github.com/kata198/mdToRst) ou de forma onlline com o [vertopal](https://www.vertopal.com/).

### 1.3.1 - Procedimento de atualização da documentação online

1. Para editar a documentação online oficial do projeto basta atualizar os arquivos *.md* na pasta **\docs**.

2. Gere uma documentação automática a partir do arquivo de configuração *mkdocs.yml* com o comando:

```bash

>> mkdocs build  # contruir os arquivos em \site\index.html

>> mmkdocs gh-deploy # deploy para o endereço de hospedagem noo github pages

```

Veja a seção [1.2 - ```__doc__``` de classes e funções]

6.Verificar se não a **bugs** na documentação.

# 2 - Testes, bugs, e documentação

## 2.1 - Utilizando o line_profiler para verificar performance do código, linha-linha

```python

@profile
def slow_function(a, b, c):
    ...
	
```

```bash

>> kernprof -l -v 'script_to_profile.py'

>> python -m line_profiler script_to_profile.py.lprof

```

## 2.2 - Teste de verificação manual, antes realizar o commit para o _main_

Executar na pasta principal do projeto \myfempy

```bash

>> black nome_do_script.py/<pasta>  # biblioteca que formata o código de acordo com a PEP 8


>> flake8 --extend-ignore=E501{line too long (82>79 characters) erro} nome_do_script.py/<pasta> #framework que checa o estilo e qualidade de código Python

>> isort nome_do_script.py/<pasta> # biblioteca Python para formatar automaticamente as importações de acordo com a PEP 8

>> interrogate -vv nome_do_script.py/<pasta> # framework para verificar a ausência de documentações (docstring) no código.

```

Estes sistemas também podem ser verificados de forma automatica via *Makefile*. Faça,

```bash
>> make install # instala o pacote com poetry

>> make format  # verifica a formação com black e isort

>> make lint    # verifica a o estilo (PEP 08) com flake8 e interrogate

>> make test    # test arquivos de itegração do pacote com pytest

>> make sec     # verifica vulnerabilidades com pip-audit
```

https://medium.com/gbtech/aprimorando-qualidade-de-c%C3%B3digo-python-com-black-flake8-isort-e-interrogate-d5c089121357

# 3 - Git

Utilize commits para subir atualizações no repositório do github (https://github.com/easycae-3d/simple_myfempy_oop).

Faça,

```bash
git status        # verifica as modificações no sistema

git add <file>... # to update what will be committed

git commit -m 'Commit menssage' # commit

git push         # push para a _main_ do projeto
```

Após realizar o envio das modificações, faça realeases quando necessário, utilize tags das versões como marcadores e gere os changelogs para documentar as modificações ao longo do projeto.

# 4 - Instalação do pacote

Faça,

```bash

pip install .   # instala o pacote myfempy do pyproject.toml

```

