from abc import ABC, abstractmethod

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
    