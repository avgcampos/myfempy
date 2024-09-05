from abc import ABC, abstractmethod


def setShape(set_shape):
    # if set_shape["type"] == "line2":
    #     from myfempy.core.shapes.line2 import Line2
    #     return Line2
    
    if set_shape["type"] == "tria3":
        from myfempy.core.shapes.tria3 import Tria3
        return Tria3
    
    elif set_shape["type"] == "tria6":
        from myfempy.core.shapes.tria6 import Tria6
        return Tria6

    elif set_shape["type"] == "quad4":
        from myfempy.core.shapes.quad4 import Quad4
        return Quad4

    elif set_shape["type"] == "quad8":
        from myfempy.core.shapes.quad8 import Quad8
        return Quad8

    # elif set_shape["type"] == "hexa8":
    #     from myfempy.core.shapes.hexa8 import Hexa8

    #     return Hexa8

    # elif set_shape["type"] == "tetr4":
    #     from myfempy.core.shapes.tetr4 import Tetra4
        # return Tetra4

    else:
        pass


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
