
# Get Started

```bash
#                     __                                
#  _ __ ___   _   _  / _|  ___  _ __ ___   _ __   _   _ 
# | '_ ` _ \ | | | || |_  / _ \| '_ ` _ \ | '_ \ | | | |
# | | | | | || |_| ||  _||  __/| | | | | || |_) || |_| |
# |_| |_| |_| \__, ||_|   \___||_| |_| |_|| .__/  \__, |
#             |___/                       |_|     |___/ 

# Copyright © Antonio Vinicius G. Campos and 3D EasyCAE, 2023. All rights reserved. 

```

# Install myfempy

1. Clone/ Download the main code [latest version] from [github/myfempy/main](https://github.com/easycae-3d/myfempy/)

2. Unzip the pack in your preferred location

3. Install the myfempy in a isolated Python environments (with virtualenv or conda)

4. Change to myfempy-main/examples folder

5. Run any file example

**Commands on terminal**

```bash

>> git clone https://github.com/easycae-3d/myfempy   # install git is needed before this command

>> python -m pip install --upgrade pip               # update pip

>> pip install numpy, cython, scipy, vedo            # install packages required
	
>> pip install .                                     # install myfempy-main

>> cd myfempy/examples

>> python ex_first_test.py                           # run any file to test installation

```

## Example List

### Static Linear Analysis

1. <ex_first_test.py> 							-- First test required to test install myfempy
2. <ex_ex_su2_Linear Elasticity_Case.py>		-- Analysis to convergence of myfempy solver with su2 code [https://su2code.github.io/tutorials/Linear_Elasticity]
3. <ex_quad4_patchtest.py>						[dev]
4. <ex_manual_mesh.py>							[dev]
5. <ex_beam_meshlegacy.py> 						[dev]								
6. <ex_plane2d_meshlegacy.py>					[dev]
7. <ex_solid_meshlegacy.py>						[dev]
8. <ex_frame2d_gmsh.py>							[dev]
9. <ex_frame3d_gmsh.py>							[dev]
10. <ex_plane2d_gmsh.py>						[dev]
11. <ex_plane2d_fromcad.py>						[dev]
12. <ex_solid_fromcad.py>						[dev]

### Modal Linear Analysis

1. <ex_beam_modal.py> 							[dev]
2. <ex_plane2d_modal.py>						[dev]
3. <ex_lumped_modal.py> 						[dev]

### Harmonic Dynamic Linear Analysis

1. <ex_beam_dyharm.py> 							[dev]
2. <ex_plane2d_dyharm.py>						[dev]