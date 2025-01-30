from abc import ABC, abstractmethod


class Thermal(ABC):
    """Load API Class <ClassService>"""

    @abstractmethod
    def getForceApply():
        pass

    @abstractmethod
    def getBCApply():
        pass

    @abstractmethod
    def getUpdateMatrix():
        pass

    @abstractmethod
    def getUpdateLoad():
        pass
