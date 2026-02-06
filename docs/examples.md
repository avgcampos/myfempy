
[TOC]

# Examples and Benchmarks

Many examples are provided for various types of problems that the project currently addresses. More examples exploring new functionalities will be added in the future. To verify the accuracy of myfempy, several Benchmarks Reports (e.g. analytics, ansys, etc.) are available below.

But first, follow these steps:

**Install myfempy**

1. Clone/ Download the main code [latest version] from [github/myfempy/main](https://github.com/avgcampos/myfempy)

2. Unzip the pack in your preferred location

3. Install the myfempy in a isolated Python environments (with virtualenv or conda)

4. Change to myfempy/docs/gallery folder

5. Run any file example

**Run this Commands on Terminal (win, linux, ...)**

```bash

>> git clone https://github.com/avgcampos/myfempy   # install git is needed before this command

>> python -m pip install --upgrade pip               # update pip
	
>> pip install .                                     # install myfempy-main

>> cd myfempy/examples/gallery                       # change to the examples gallery directory

>> python run_first_test.py                           # run any file to check installation

```

## Example List

Many Examples are available to download [here](https://github.com/avgcampos/myfempy/tree/main/examples).

### Static Linear Analysis

#### run_first_test.py
You can use this code to first test required to check myfempy install
[Relatório de Benchmark (PT-BR PDF)](benchs/estatico1.pdf)

```python linenums="1"
--8<-- "examples/run_first_test.py"
```

#### static_linear_beam3d_LoganBook.py
[Relatório de Benchmark (PT-BR PDF)](benchs/estatico2.pdf)

```python linenums="1"
--8<-- "examples/static_linear_beam3d_LoganBook.py"
```

#### myfempy_nafems_le1.py
[Relatório de Benchmark (PT-BR PDF)](benchs/estatico3.pdf)

```python linenums="1"
--8<-- "examples/myfempy_nafems_le1.py"
```

#### patchtest_quad4.py
[Relatório de Benchmark (PT-BR PDF)](benchs/patch_test.pdf)

```python linenums="1"
--8<-- "examples/patchtest_quad4.py"
```

#### myfempyVSsu2_StaticLinear.py
Analysis to convergence of myfempy solver with [su2](https://su2code.github.io/tutorials/Linear_Elasticity)

Resposta do MYFEMPY
![Relatório de Benchmark (PT-BR PDF)](gallery/myfempy_su2_test_displ.png)
```python linenums="1"
--8<-- "examples/myfempyVSsu2_StaticLinear.py"
```

### Vibration Modal Linear Analysis

#### vibration_plane.py
[Relatório de Benchmark (PT-BR PDF)](benchs/vibracoes1.pdf)

```python linenums="1"
--8<-- "examples/vibration_plane.py"
```
#### vibration_beam.py

```python linenums="1"
--8<-- "examples/vibration_beam.py"
```

#### vibration_solid.py
[Relatório de Benchmark (PT-BR PDF)](benchs/vibracoes2.pdf)

```python linenums="1"
--8<-- "examples/vibration_solid.py"
```

#### tuning_fork.py
[Relatório de Benchmark (PT-BR PDF)](benchs/vibracoes3.pdf)

```python linenums="1"
--8<-- "examples/tuning_fork.py"
```

### Vibration Harmonic Linear Analysis

#### vibration_frf.py
```python linenums="1"
--8<-- "examples/gallery/vibration_frf.py"
```

### Heat Flow

#### heat_StSd.py
```python linenums="1"
--8<-- "examples/gallery/heat_StSd.py"
```

### Cyclic Symmetric

#### cyclic_symmetric_plane.py
```python linenums="1"
--8<-- "examples/gallery/cyclic_symmetric_plane.py"
```

### Thermal-Mechanics

#### coup_struct_heat_StSd.py
```python linenums="1"
--8<-- "examples/gallery/coup_struct_heat_StSd.py"
```

#### sim_tira_bimetalica.py
```python linenums="1"
--8<-- "examples/gallery/sim_tira_bimetalica.py"
```

### Multi-Scale

#### myfempy_multi_escala_Al_core_square_hole.py
```python linenums="1"
--8<-- "examples/gallery/myfempy_multi_escala_Al_core_square_hole.py"
```