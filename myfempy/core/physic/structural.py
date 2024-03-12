from abc import ABC, abstractmethod

class Structural(ABC):
    '''Load API Class <ClassService>'''
    
    @abstractmethod 
    def getForceApply():
        pass
 
    @abstractmethod
    def setLoadDof():
        pass
    
    @abstractmethod
    def getBCApply():
        pass 

    @abstractmethod
    def setBCDof():
        pass
