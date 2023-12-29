from abc import ABC, abstractmethod

import numpy as np

def getElement(set_element):
    if set_element['type'] == 'plane':
        from myfempy.felib.elements.plane import Plane
        return Plane
    
    elif set_element['type'] == 'platekc':
        from myfempy.felib.elements.platekc import PlateKC
        return PlateKC
    
    elif set_element['type'] == 'solid':
        from myfempy.felib.elements.solid import Solid
        return Solid
    
    else:
        pass

class Element(ABC):
    '''Element API Class <ClassService>'''
                        
    @abstractmethod
    def getElementSet():
        pass

    @abstractmethod
    def setElementShape():
        pass

    @abstractmethod
    def setElementMesh():
        pass

    @abstractmethod
    def setElementMaterial():
        pass

    @abstractmethod
    def setElementGeometry():
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
    def setTitleDeformation():
        pass

    @abstractmethod
    def getElementVolume():
        pass