from abc import ABC, abstractmethod

class Structural(ABC):
    '''Load API Class <ClassService>'''
    
    @abstractmethod 
    def getForceApply():
        pass
    
    @abstractmethod
    def getBCApply():
        pass 