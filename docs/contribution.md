# General Styles Used in the Project

**Rules for Uploading {*commit*} to the myfempy Project**

This repository is intended for managing the development versions of the **myfempy** project. To maintain clear and clean code, the following rules are applied to ensure a general project standard.

## Style of *.py* Codes

The syntax used for writing the code is [PEP8](https://peps.python.org/pep-0008).

### Class Creation Pattern *class AbcAbc*

```python

# the file name is abcabc.py

class AbcAbc: # class name

def __init__(self, x:int, y:float, z:str): # constructor method in python

# add parameter in class
self.x = x
self.y = y
self.z = z

```

### Function Creation Pattern *def abc_abc*

```python

# the file name is abcabc.py

def abc_abc(x:int, y:float): # function

sum_x = x + x

div_y = x/y

return div_y

```

## ```__doc__``` for classes and functions

Every class and function must contain its internal documentation {```__doc__```}.

### Procedure

1. Use the reference [PEP8](https://peps.python.org/pep-0008)

2. Present in a single line what that class or function does

3. Use comments for additional information

4. All inputs and outputs and data types must be indicated

5. The default language is technical English

Example:

```python

# the file name is abcabc.py
# '__doc__' to a short documentation, try AbcAbc.__doc__ ....

__doc__ = """ This is a short documentation of class AbcAbc...
"""
class AbcAbc: # class name
"""
Abc: class

This class does...

"""

def __init__(self, x:int, y:float, z:str): # constructor method in python
"""
This class do...

Args:

x (int): int number

y (float): float number
z (str): string name

""

# add parameter in class
self.x = x
self.y = y
self.z = z

def abc_abc(self): # class function

""
This function do...

Returns:

div_y (float): _description_

""

soma_x = self.x + self.x

div_y = self.x/soma_x

return div_y

```

## Online Documentation

The standard used to write the official project documentation is in [Markdown](https://www.markdownguide.org/basic-syntax) format.

Currently, the project uses the [mkdocs](https://www.mkdocs.org/) package to generate online documentation {*.html*}, which is hosted on the [GitHub Page](https://avgcampos.github.io/myfempy/) website.

*sphinx* converts *.rst* {reStructuredText} files to *.html/pdf* in a practical way. It is also possible to convert *.md* {markdown} files to *.rst*; an interesting option is [mdToRst](https://github.com/kata198/mdToRst) or online with [vertopal](https://www.vertopal.com/).

### Procedure for updating the online documentation

1. To edit the official online project documentation, simply update the *.md* files in the **\docs** folder.

2. Generate automatic documentation from the *mkdocs.yml* configuration file using the command:

```bash

>> mkdocs build # builds the files in \site

>> mkdocs gh-deploy # deploys to the hosting address on GitHub Pages

or

>> make doc # builds and generates local documentation

```

See the [```__doc__``` section on classes and functions]

6. Check for **bugs** in the documentation.

# Tests, bugs, and documentation

## Using line_profiler to check code performance, line by line

```python

@profile
def slow_function(a, b, c):

...

```

```bash

>> kernprof -l -v 'script_to_profile.py'

>> python -m line_profiler script_to_profile.py.lprof

```

## Manual verification test, before committing to _main_

Run in the project's main folder \myfempy

```bash

>> black script_name.py/<folder> # library that formats the code according to PEP 8

>> flake8 --extend-ignore=E501{line too long (82>79 characters) error} script_name.py/<folder> # framework that checks the style and quality of Python code

>> `isort script_name.py/<folder>` # Python library to automatically format imports according to PEP 8

`>> interrogate -vv script_name.py/<folder>` # framework to check for missing documentation (docstring) in the code.

```

These systems can also be checked automatically via *Makefile*. Do the following:

```bash
>> make install # installs the package with poetry

>> make format # checks the formatting with black and isort

>> make lint # checks the style (PEP 08) with flake8 and interrogate

>> make test # tests the package integration files with pytest

>> make sec # checks for vulnerabilities with pip-audit

```

https://medium.com/gbtech/aprimorando-qualidade-de-c%C3%B3digo-python-com-black-flake8-isort-e-interrogate-d5c089121357

# Git

Do the following:

```bash
git status # checks the system modifications

git add <file>... # to update what