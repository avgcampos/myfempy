from myfempy.io.iogmsh import get_gmsh_geo, get_gmsh_msh

def setElement(set_element):
    if set_element["TYPE"] == "structbeam":
        from myfempy.core.elements.structBeam import StructuralBeam
        return StructuralBeam

    elif set_element["TYPE"] == "structplane":
        from myfempy.core.elements.structPlane import StructuralPlane
        return StructuralPlane

    # TO DO
    elif set_element["TYPE"] == "structplate":
        pass

    elif set_element["TYPE"] == "structsolid":
        from myfempy.core.elements.structSolid import StructuralSolid
        return StructuralSolid

    elif set_element["TYPE"] == "heatplane":
        from myfempy.core.elements.heatPlane import HeatPlane
        return HeatPlane
    
    elif set_element["TYPE"] == "heatsolid":
        from myfempy.core.elements.heatSolid import HeatSolid
        return HeatSolid
    
    elif set_element["TYPE"] == "userelement":
        NewClass = set_element["CLASS"]
        return NewClass

    else:
        pass


def setShape(set_shape):
    if set_shape["SHAPE"] == "line2":
        from myfempy.core.shapes.line2 import Line2
        return Line2

    elif set_shape["SHAPE"] == "line3":
        from myfempy.core.shapes.line3 import Line3
        return Line3

    elif set_shape["SHAPE"] == "tria3":
        from myfempy.core.shapes.tria3 import Tria3
        return Tria3

    elif set_shape["SHAPE"] == "tria6":
        from myfempy.core.shapes.tria6 import Tria6
        return Tria6

    elif set_shape["SHAPE"] == "quad4":
        from myfempy.core.shapes.quad4 import Quad4
        return Quad4

    elif set_shape["SHAPE"] == "quad8":
        from myfempy.core.shapes.quad8 import Quad8
        return Quad8

    elif set_shape["SHAPE"] == "hexa8":
        from myfempy.core.shapes.hexa8 import Hexa8
        return Hexa8

    elif set_shape["SHAPE"] == "tetr4":
        from myfempy.core.shapes.tetr4 import Tetra4
        return Tetra4
    
    elif set_shape["SHAPE"] == "usershape":
        NewClass = set_shape["CLASS"]
        return NewClass

    else:
        pass


def setMesh(set_mesh):
    if set_mesh["TYPE"] == "manual":
        from myfempy.core.mesh.mesh import MeshMANUAL
        return MeshMANUAL

    elif set_mesh["TYPE"] == "legacy":
        if set_mesh["SHAPE"] == "quad4":
            from myfempy.core.mesh.legacyquad4 import LegacyQuad4
            return LegacyQuad4

        elif set_mesh["SHAPE"] == "tria3":
            from myfempy.core.mesh.legacytria3 import LegacyTria3
            return LegacyTria3

        elif set_mesh["SHAPE"] == "line2":
            from myfempy.core.mesh.legacyline2 import LegacyLine2
            return LegacyLine2

        else:
            pass

    elif set_mesh["TYPE"] == "gmsh":
        if "meshimport" in set_mesh.keys():
            pass
        else:
            filename = set_mesh["user_path"] + "/" + set_mesh["filename"]
            get_gmsh_geo(filename, set_mesh)
            get_gmsh_msh(filename, set_mesh)
        from myfempy.core.mesh.gmsh import MeshGmsh
        return MeshGmsh
    else:
        pass


def setMaterial(set_material):
    if set_material["MAT"] == "uniaxialstress":
        if set_material["TYPE"] == "isotropic":
            from myfempy.core.material.uniaxialstress import UniAxialStress
            return UniAxialStress
        
        elif set_material["TYPE"] == "usermaterial":
            NewMatClass = set_material["CLASS"]
            return NewMatClass    
    
        else:
            pass

    elif set_material["MAT"] == "planestress":
        if set_material["TYPE"] == "isotropic":
            from myfempy.core.material.planestress import PlaneStress
            return PlaneStress
        
        elif set_material["TYPE"] == "usermaterial":
            NewMatClass = set_material["CLASS"]
            return NewMatClass
            
        else:
            pass


    elif set_material["MAT"] == "planestrain":
        if set_material["TYPE"] == "isotropic":
            from myfempy.core.material.planestrain import PlaneStrain
            return PlaneStrain
        
        elif set_material["TYPE"] == "usermaterial":
            NewMatClass = set_material["CLASS"]
            return NewMatClass            
        
        else:
            pass

    elif set_material["MAT"] == "solidelastic":
        if set_material["TYPE"] == "isotropic":
            from myfempy.core.material.solidelastic import SolidElastic
            return SolidElastic
        
        elif set_material["TYPE"] == "usermaterial":
            NewMatClass = set_material["CLASS"]
            return NewMatClass            
        
        else:
            pass

    elif set_material["MAT"] == "heatplane":
        if set_material["TYPE"] == "isotropic":
            from myfempy.core.material.heatplane import HeatPlane
            return HeatPlane
        
        elif set_material["TYPE"] == "usermaterial":
            NewMatClass = set_material["CLASS"]
            return NewMatClass            
        
        else:
            pass

    elif set_material["MAT"] == "heatsolid":
        if set_material["TYPE"] == "isotropic":
            from myfempy.core.material.heatsolid import HeatSolid
            return HeatSolid
        
        elif set_material["TYPE"] == "usermaterial":
            NewMatClass = set_material["CLASS"]
            return NewMatClass            
        
        else:
            pass
    
    else:
        pass


def setGeometry(set_geometry):
    if set_geometry["GEO"] == "thickness":
        from myfempy.core.geometry.thickness import Thickness
        return Thickness

    elif set_geometry["GEO"] == "rectangle":
        from myfempy.core.geometry.rectangle import Rectangle
        return Rectangle

    elif set_geometry["GEO"] == "rectangle_tube":
        from myfempy.core.geometry.rectangle_tube import RectangleTube
        return RectangleTube

    elif set_geometry["GEO"] == "circle":
        from myfempy.core.geometry.circle import Circle
        return Circle

    elif set_geometry["GEO"] == "circle_tube":
        from myfempy.core.geometry.circle_tube import CircleTube
        return CircleTube

    elif set_geometry["GEO"] == "isection":
        from myfempy.core.geometry.isection import ISection
        return ISection

    elif set_geometry["GEO"] == "userdefined":
        from myfempy.core.geometry.userdefined import UserDefined
        return UserDefined

    elif set_geometry["TYPE"] == "usergeometry":
        NewClass = set_geometry["CLASS"]
        return NewClass

    else:
        pass


def setDomain(set_domain):
    if set_domain["DOMAIN"] == "structural":
        from myfempy.core.physic.bcstruct import BoundCondStruct
        from myfempy.core.physic.loadstruct import LoadStructural
        return LoadStructural, BoundCondStruct

    elif set_domain["DOMAIN"] == "thermal":
        from myfempy.core.physic.bcthermal import BoundCondThermal
        from myfempy.core.physic.loadthermal import LoadThermal
        return LoadThermal, BoundCondThermal

    elif set_domain["DOMAIN"] == "fluidflow":
        pass
    
    else:
        pass
    
def setCoupling(set_coupling):
    # themal structural iteration
    if set_coupling["TYPE"] == "thermalstress":
        from myfempy.core.physic.thermstructcoup import \
            ThermalStructuralCoupling
        from myfempy.core.physic.bcstruct import BoundCondStruct
        return ThermalStructuralCoupling, BoundCondStruct

    # fluidflow structural iteration
    # elif physicdata["COUPLING"]['TYPE'] == "fsi":

    # acoustic structural iteration
    # elif physicdata["COUPLING"]['TYPE'] == "asi":

    else:
        pass
    
def setPoints2NumericalIntegration(type_shape):
    shape = {
            "line2": 2,
            "line3": 4,
            "tria3": 1,
            "tria6": 7,
            "quad4": 2,
            "quad8": 3,
            "tetr4": 4,
            "tetr10": 5,
            "hexa8": 2,
            "hexa20": 3,
        }
    return shape[type_shape]