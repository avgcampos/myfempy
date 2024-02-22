from abc import ABC, abstractmethod

class Solver(ABC):
    
    """
     Solver API Class <ClassService>
    """
        
    @abstractmethod
    def runSolve():
        
        """
        Solve Solve FE System: Model + Physics
        """
        pass
    
    @abstractmethod
    def getMatrixAssembler():
        pass
    
    
    @abstractmethod
    def getLoadAssembler():
        pass
    
    
    @abstractmethod
    def getConstrains():
        pass