from abc import ABC, abstractmethod


def setGeometry(set_geometry):
    if set_geometry["geo"] == "thickness":
        from myfempy.core.geometry.thickness import Thickness

        return Thickness
    else:
        pass


class Geometry(ABC):
    """Geoemtry API Class <ClassService>"""

    def __init__(self) -> None:
        pass

    @abstractmethod
    def GeometrySet():
        pass

    @abstractmethod
    def getSectionProp():
        pass

    @abstractmethod
    def getCGCoord():
        pass
