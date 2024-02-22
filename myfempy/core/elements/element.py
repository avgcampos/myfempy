from abc import ABC, abstractmethod

import numpy as np

def setElement(set_element):
    if set_element['type'] == 'plane':
        from myfempy.core.elements.plane import Plane
        return Plane
    
    elif set_element['type'] == 'platekc':
        from myfempy.core.elements.platekc import PlateKC
        return PlateKC
    
    elif set_element['type'] == 'solid':
        from myfempy.core.elements.solid import Solid
        return Solid
    
    else:
        pass

class Element(ABC):
    '''Element API Class <ClassService>'''
                        
    @abstractmethod
    def getElementSet():
        pass

    @abstractmethod
    def getElementShape():
        pass

    @abstractmethod
    def getElementMesh():
        pass

    @abstractmethod
    def getElementMaterial():
        pass

    @abstractmethod
    def getElementGeometry():
        pass

    @abstractmethod
    def getB():
        pass

    @abstractmethod
    def getH():
        pass
           
    @abstractmethod
    def getStifLinearMat():
        pass
    
    @abstractmethod
    def getStifNonLinMat():
        pass
    
    @abstractmethod
    def getMassConsistentMat():
        pass
    
    @abstractmethod
    def getMassLumpedMat():
        pass
    
    @abstractmethod
    def getDamperConsistentMat():
        pass
    
    @abstractmethod
    def getElementDeformation():
        pass

    @abstractmethod
    def getTitleDeformation():
        pass

    @abstractmethod
    def getElementVolume():
        pass