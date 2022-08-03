# -*- coding: utf-8 -*-
"""
========================================================================
~~~ MODULO DE SIMULACAO ESTRUTURAL PELO METODO DOS ELEMENTOS FINITOS ~~~
       	                    __                                
       	 _ __ ___   _   _  / _|  ___  _ __ ___   _ __   _   _ 
       	| '_ ` _ \ | | | || |_  / _ \| '_ ` _ \ | '_ \ | | | |
       	| | | | | || |_| ||  _||  __/| | | | | || |_) || |_| |
       	|_| |_| |_| \__, ||_|   \___||_| |_| |_|| .__/  \__, |
       	            |___/                       |_|     |___/ 

~~~      Mechanical studY with Finite Element Method in PYthon       ~~~
~~~                PROGRAMA DE ANÁLISE COMPUTACIONAL                 ~~~
~~~              copyright @ 2022, all rights reserved               ~~~
========================================================================
"""
import numpy as np
from myfempy.felib.felemset import get_elemset


def get_stress(modelinfo, U, num_elm):

    if modelinfo['tabmat'][int(modelinfo['inci'][num_elm, 2])-1, -1] == 2:
        from myfempy.felib.materials.axial import Tensor
        return Tensor(modelinfo, U, num_elm)

    elif modelinfo['tabmat'][int(modelinfo['inci'][num_elm, 2])-1, -1] == 3:
        from myfempy.felib.materials.planestress import Tensor
        return Tensor(modelinfo, U, num_elm)

    elif modelinfo['tabmat'][int(modelinfo['inci'][num_elm, 2])-1, -1] == 4:
        pass

    elif modelinfo['tabmat'][int(modelinfo['inci'][num_elm, 2])-1, -1] == 5:
        from myfempy.felib.materials.solid import Tensor
        return Tensor(modelinfo, U, num_elm)

    elif modelinfo['tabmat'][int(modelinfo['inci'][num_elm, 2])-1, -1] == 11:
        print('Not implemented yet')


def get_displ(modelinfo, U):

    element = get_elemset(int(modelinfo['elemid'][0]))

    dofdef = element.elemset()

    if dofdef['def'] == 'struct 1D':
        if dofdef['dofs'] == ['ux', 'uy']:
            from myfempy.postprc.displcalc import Deformation
            dis = Deformation(modelinfo)
            return dis.ux(U)

        elif dofdef['dofs'] == ['uy', 'rz']:
            from myfempy.postprc.displcalc import Deformation
            dis = Deformation(modelinfo)
            return dis.uy_rz(U)

        elif dofdef['dofs'] == ['ux', 'uy', 'rz']:
            from myfempy.postprc.displcalc import Deformation
            dis = Deformation(modelinfo)
            return dis.ux_uy_rz(U)

        elif dofdef['dofs'] == ['ux', 'uy', 'uz', 'rx', 'ry', 'rz']:
            from myfempy.postprc.displcalc import Deformation
            dis = Deformation(modelinfo)
            return dis.ux_uy_uz_rx_ry_rz(U)

    elif dofdef['def'] == 'struct 2D':
        from myfempy.postprc.displcalc import Deformation
        dis = Deformation(modelinfo)
        return dis.ux_uy(U)

    elif dofdef['def'] == 'struct 3D':
        from myfempy.postprc.displcalc import Deformation
        dis = Deformation(modelinfo)
        return dis.ux_uy_uz(U)
