from myfempy.core.material.material import Material
__docformat__ = "google"

__doc__ = """
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


class NullMat(Material):
    """
    Null Material Class <ConcreteClassService>
    """

    def getMaterialSet():
        matset = {
            "mat": "null",
            "type": "null",
        }
        return matset

    def getElasticTensor(tabmat, inci, element_number, Model=None):
        return 0.0

    def getElementStrain(Model, U, ptg, element_number):
        epsilon = 0.0
        strain = [
            0.00000000000,
        ]
        return epsilon, strain

    def getTitleStrain():
        title = [
            "NULL",
        ]
        return title

    def getElementStress(Model, epsilon, element_number):
        sigma = 0.0
        stress = [
            0.00000000000,
        ]
        return sigma, stress

    def getTitleStress():
        title = [
            "NULL",
        ]
        return title

    def getStrainEnergyDensity(sigma, epsilon, elemvol):
        return 0.0

    def getTitleCompliance():
        title = ["NULL"]
        return title

    def getFailureCriteria(sigma):
        return 0.0

    def getTitleFoS():
        title = ["NULL"]
        return title
