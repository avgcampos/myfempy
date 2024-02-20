from abc import ABC, abstractmethod


def getMaterial(set_material):
    if set_material['mat'] == 'axial':
        pass
        
    elif set_material['mat'] == 'planestress':
        if set_material['type'] == 'isotropic':
            from myfempy.core.material.planestress import PlaneStressIsotropic
            return PlaneStressIsotropic
        else:
            pass
                
    elif set_material['mat'] == 'planestrain':
        if set_material['type'] == 'isotropic':
            from myfempy.core.material.planestrain import PlaneStrainIsotropic
            return PlaneStrainIsotropic
        else:
            pass
        
    elif set_material['mat'] == 'axisymmetric':
        pass
        
    elif set_material['mat'] == 'solid':
        if set_material['type'] == 'isotropic':
            from myfempy.core.material.solid import SolidIsotropic
            return SolidIsotropic
        else:
            pass
    
    
       
    
class Material(ABC):
    '''Material API Class <ClassService>'''
    
    def __init__(self) -> None:
        pass

    @abstractmethod
    def MaterialSet():
        pass
    
    @abstractmethod
    def getElasticTensor():
        pass
    
    @abstractmethod
    def getElementStrain():
        pass
    
    @abstractmethod
    def getElementStress():
        pass

    @abstractmethod
    def setTitleStrain():
        pass

    @abstractmethod
    def setTitleStress():
        pass