from abc import ABC, abstractmethod


class Material(ABC):
    """Material API Class <ClassService>"""

    def __init__(self) -> None:
        pass

    @abstractmethod
    def getMaterialSet():
        pass

    @abstractmethod
    def getElasticTensor():
        pass
