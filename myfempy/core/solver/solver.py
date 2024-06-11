from abc import ABC, abstractmethod


class Solver(ABC):

    """
    Solver API Class <ClassService>
    """

    @abstractmethod
    def runSolve():
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

    @abstractmethod
    def getDirichletNH():
        pass
