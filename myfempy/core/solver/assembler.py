from abc import ABC, abstractmethod


class Assembler(ABC):
    """
    Assembler API Class <ClassService>
    """

    @abstractmethod
    def getGlobalMatrixAssembler():
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

    @abstractmethod
    def getRotationMatrix():
        pass

    @abstractmethod
    def getSaveAssemblerFile():
        pass