from abc import ABC, abstractmethod


class Shape(ABC):
    """Shape API Class <ClassService>"""

    @abstractmethod
    def getShapeSet():
        pass

    @abstractmethod
    def getShapeFunctions():
        pass

    @abstractmethod
    def getDiffShapeFuntion():
        pass

    @abstractmethod
    def getJacobian():
        pass

    @abstractmethod
    def getinvJacobi():
        pass

    @abstractmethod
    def getdetJacobi():
        pass

    @abstractmethod
    def getNodeList():
        pass

    @abstractmethod
    def getNodeCoord():
        pass

    @abstractmethod
    def getLocKey():
        pass

    @abstractmethod
    def getNodesLines():
        pass

    @abstractmethod
    def getIsoParaSide():
        pass

    @abstractmethod
    def getEdgeLength():
        pass

    @abstractmethod
    def getSideAxis():
        pass