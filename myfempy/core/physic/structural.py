from abc import ABC, abstractmethod


class Structural(ABC):
    """Load API Class <ClassService>"""

    @abstractmethod
    def getLoadApply():
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
