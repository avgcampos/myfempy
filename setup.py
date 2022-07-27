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
# SETUP SYSTEM & PIP INSTALL
from setuptools import setup, find_packages
 
try:
    VERSIONFILE = "myfempy/version.py"
    verstrline = open(VERSIONFILE,"rt").read()
    verstr = verstrline.split("=")[1].replace("\n", "").replace("'", "")
except:
    verstr = "unknown"


setup(
    python_requires = '>=3',
    include_package_data = True,
    name = 'myfempy',
    version = verstr,
    license = 'GNU',
    license_files = ['LICENSE'], 
    
    author = 'Campos, A. V. G.',
    maintainer = "Campos, A. V. G. & 3D EasyCAE",
    maintainer_email = '3deasycaebr.contato@gmail.com',
    
    description = "A python package for scientific analysis based on finite element method",
    long_description = "The myfempy is a python based on finite element method for scientific analysis. The code is open source and intended for educational and scientific purposes only, not recommended to commercial use. You can help us by contributing with a donation on the main project page, read the support options. If you use myfempy in your research, the  developers would be grateful if you could cite in your work",
    long_description_content_type='text/markdown',
    
    url = 'https://myfempy.readthedocs.io/en/latest/',
    download_url = 'https://github.com/easycae-3d/myfempy',
    keywords = ['Finite Element', 'Mechanics', 'Python Package'],
    
    packages = find_packages(),
    
    # packages = [
        # "myfempy",
        # "myfempy.core",
        # "myfempy.felib",
        # "myfempy.felib.fluid",
        # "myfempy.felib.fsi",
        # "myfempy.felib.materials",
        # "myfempy.felib.physics",
        # "myfempy.felib.struct",
        # "myfempy.io",
        # "myfempy.mesh",
        # "myfempy.plots",
        # "myfempy.postprc",
        # "myfempy.tools",
    # ],
    
    package_dir = {
        "myfempy": "myfempy",
    },
    
    install_requires = [
        "numpy",
        "scipy",
        "matplotlib",
        "gmsh",
        "vedo",
    ],
    
    # include_package_data=True,
    zip_safe = False,
        
    classifiers = [
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Programming Language :: Python",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Education",
        "Topic :: Scientific/Engineering :: Physics",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: POSIX :: Linux",
        "Operating System :: MacOS :: MacOSX",
    ],
)














