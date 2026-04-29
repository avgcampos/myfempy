
# >> python setup_wrap_cy_pyx.py build_ext --inplace


__docformat__ = "google"

__doc__ = """

==========================================================================
                            __                                
         _ __ ___   _   _  / _|  ___  _ __ ___   _ __   _   _ 
        | '_ ` _ \ | | | || |_  / _ \| '_ ` _ \ | '_ \ | | | |
        | | | | | || |_| ||  _||  __/| | | | | || |_) || |_| |
        |_| |_| |_| \__, ||_|   \___||_| |_| |_|| .__/  \__, |
                    |___/                       |_|     |___/ 
        myfempy -- MultiphYsics Finite Element Module to PYthon    
                    COMPUTATIONAL ANALYSIS PROGRAM                   
        Copyright (C) 2022-2026 Antonio Vinicius Garcia Campos        
==========================================================================
This Python file is part of myfempy project.

myfempy is a python package based on finite element method to multiphysics
analysis. The code is open source and *intended for educational and scientific
purposes only, not recommended to commercial use. The name myfempy is an acronym
for MultiphYsics Finite Elements Module to PYthon. You can help us by contributing
with the main project, send us a mensage on https://github.com/avgcampos/myfempy/discussions/10
If you use myfempy in your research, the  developers would be grateful if you 
could cite in your work.
																		
The code is written by Antonio Vinicius Garcia Campos.                                  
																		
A github repository, with the most up to date version of the code,      
can be found here: https://github.com/avgcampos/myfempy.                 
																		
The code is open source and intended for educational and scientific     
purposes only. If you use myfempy in your research, the developers      
would be grateful if you could cite this. The myfempy project is published
under the GPLv3, see the myfempy LICENSE on
https://github.com/avgcampos/myfempy/blob/main/LICENSE.
																		
Disclaimer:                                                             
The authors reserve all rights but do not guarantee that the code is    
free from errors. Furthermore, the authors shall not be liable in any   
event caused by the use of the program.

"""


import platform
from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext as _build_ext
from Cython.Build import cythonize
import numpy

include_dirs = [numpy.get_include()]

if platform.system() == 'Windows':
    compile_args = ['/openmp']
    link_args = []
elif platform.system() == 'Linux':
    compile_args = ['-fopenmp', '-static', '-static-libgcc', '-static-libstdc++']
    link_args = ['-fopenmp', '-static-libgcc', '-static-libstdc++']
else: # MAC-OS
    compile_args = []
    link_args = []

extension_kwargs = dict(
    include_dirs=include_dirs,
    extra_compile_args=compile_args,
    extra_link_args=link_args,
    language='c',
    )

class build_ext(_build_ext):
    def build_extensions(self):
        super().build_extensions()

    def initialize_options(self):
        super().initialize_options()
        if self.distribution.ext_modules == None:
            self.distribution.ext_modules = []

        extensions=[
            Extension("*", sources=[".myfempy/core/solver/assemblersymm_cython.pyx"], **extension_kwargs),
            Extension("*", sources=[".myfempy/core/solver/assemblerfull_cython.pyx"], **extension_kwargs),
            Extension("*", sources=[".myfempy/core/shapes/line2_tasks.pyx"], **extension_kwargs),
            Extension("*", sources=[".myfempy/core/shapes/line3_tasks.pyx"], **extension_kwargs),
            Extension("*", sources=[".myfempy/core/shapes/tria3_tasks.pyx"], **extension_kwargs),
            Extension("*", sources=[".myfempy/core/shapes/tria6_tasks.pyx"], **extension_kwargs),
            Extension("*", sources=[".myfempy/core/shapes/quad4_tasks.pyx"], **extension_kwargs),
            Extension("*", sources=[".myfempy/core/shapes/quad8_tasks.pyx"], **extension_kwargs),
            Extension("*", sources=[".myfempy/core/shapes/tetr4_tasks.pyx"], **extension_kwargs),
            Extension("*", sources=[".myfempy/core/shapes/hexa8_tasks.pyx"], **extension_kwargs),
        ]
        self.distribution.ext_modules.extend(cythonize(extensions))