from abc import ABC, abstractmethod


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
