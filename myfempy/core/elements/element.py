from abc import ABC, abstractmethod

import numpy as np


def setElement(set_element):
    if set_element["type"] == "structplane":
        from myfempy.core.elements.structPlane import StructuralPlane
        return StructuralPlane

    elif set_element["type"] == "structplate":
        from myfempy.core.elements.structPlate import StructuralPlate
        return StructuralPlate

    elif set_element["type"] == "structsolid":
        from myfempy.core.elements.structSolid import StructuralSolid
        return StructuralSolid
    
    elif set_element["type"] == "heatplane":
        from myfempy.core.elements.heatPlane import HeatPlane
        return HeatPlane

    else:
        pass


class Element(ABC):
    """Element API Class <ClassService>"""

    @abstractmethod
    def getElementSet():
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
    def getUpdateMatrix():
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

    def get_side_fcapp(set_side):
        side = {
            "0 1": "0",  # quad4
            "1 0": "0",  # quad4
            "1 2": "1",  # quad4
            "2 1": "1",  # quad4
            "2 3": "2",  # quad4
            "3 2": "2",  # quad4
            "3 0": "3",  # quad4
            "0 3": "3",  # quad4
            "2 0": "2",  # tria3
            "0 2": "2",  # tria3
        }
        return side[set_side]