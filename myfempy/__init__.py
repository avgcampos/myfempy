#!/usr/bin/env python
# ==========================================================================#
#  This Python file is part of myfempy project                             #
#                                                                          #
#  The code is written by A. V. G. Campos                                  #
#                                                                          #
#  A github repository, with the most up to date version of the code,      #
#  can be found here:                                                      #
#     https://github.com/easycae-3d/myfempy                                #
#                                                                          #
#  The code is open source and intended for educational and scientific     #
#  purposes only. If you use myfempy in your research, the developers      #
#  would be grateful if you could cite this.                               #
#                                                                          #
#  Disclaimer:                                                             #
#  The authors reserve all rights but do not guarantee that the code is    #
#  free from errors. Furthermore, the authors shall not be liable in any   #
#  event caused by the use of the program.                                 #
# ==========================================================================#
"""
myfempy -- MultiphYsics Finite Element Method with PYthon
Copyright (C) 2023 Antonio Vinicius Garcia Campos
"""
# ==========================================================
from __future__ import absolute_import
import pkg_resources

# solvers
from .core.staticlinear import StaticLinear
from .core.modallinear import ModalLinear
from .core.harmonicforced import HarmonicForced

# fea analysis
from .setup.feaanalysis import FEANewAnalysis

# ==========================================================
# __author__ = pkg_resources.get_distribution('myfempy').authors
# __copyright__ = "Copyright @ 2023, Antonio Vinicius Garcia Campos"
# __license__ = pkg_resources.get_distribution('myfempy').license
# __website__ = pkg_resources.get_distribution('myfempy').homepage
__version__ = 'dev' #pkg_resources.get_distribution('myfempy').version
# ==========================================================
__all__ = ["FEANewAnalysis", "StaticLinear", "ModalLinear", "HarmonicForced", "__version__"]
# ==========================================================
