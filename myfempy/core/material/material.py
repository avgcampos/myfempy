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

    # @abstractmethod
    # def getElementStrain():
    #     pass

    # @abstractmethod
    # def getElementStress():
    #     pass

    # @abstractmethod
    # def getFailureCriteria():
    #     pass

    # @abstractmethod
    # def getTitleStrain():
    #     pass

    # @abstractmethod
    # def getTitleStress():
    #     pass

    # @abstractmethod
    # def getTitleFoS():
    #     pass
