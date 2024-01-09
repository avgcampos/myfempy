# Contribution

## Regras para realizar upload {*commit*} no projeto myfempy

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

## 1.3 - Documentação online

O padrão utilizado para escrever a documentação oficial do projeto é no formato [Markdown](https://www.markdownguide.org/basic-syntax).

Atualmente o projeto utiliza o pacote [mkdocs](https://www.mkdocs.org/) para gerar a documentação online {*.html*}, que é hospedada no site [GitHub Page](https://avgcampos.github.io/myfempy/).

O *sphinx* converte arquivos *.rst* {reStruturedText} para *.html/pdf* de uma forma prática. Também é possível converter arquivos *.md* {markdown} em *.rst*, uma opção interessante é [mdToRst](https://github.com/kata198/mdToRst) ou de forma onlline com o [vertopal](https://www.vertopal.com/).

### 1.3.1 - Procedimento de atualização da documentação online

1. Para editar a documentação online oficial do projeto basta atualizar os arquivos *.md* na pasta **\docs**.

2. Gere uma documentação automática a partir do arquivo de configuração *mkdocs.yml* com o comando:

```bash

>> mkdocs build  # constroe os arquivos em \site

>> mkdocs gh-deploy # deploy para o endereço de hospedagem no github pages

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