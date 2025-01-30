from myfempy.io.iogmsh import get_gmsh_geo, get_gmsh_msh

def setElement(set_element):
    if set_element["type"] == "structbeam":
        from myfempy.core.elements.structBeam import StructuralBeam
        return StructuralBeam
    
    elif set_element["type"] == "structplane":
        from myfempy.core.elements.structPlane import StructuralPlane
        return StructuralPlane

    # elif set_element["type"] == "structplate":
    #     from myfempy.core.elements.structPlate import StructuralPlate
    #     return StructuralPlate

    elif set_element["type"] == "structsolid":
        from myfempy.core.elements.structSolid import StructuralSolid
        return StructuralSolid
    
    elif set_element["type"] == "heatplane":
        from myfempy.core.elements.heatPlane import HeatPlane
        return HeatPlane

    else:
        pass
    

def setShape(set_shape):
    if set_shape["type"] == "line2":
        from myfempy.core.shapes.line2 import Line2
        return Line2
    
    elif set_shape["type"] == "line3":
        from myfempy.core.shapes.line3 import Line3
        return Line3
    
    elif set_shape["type"] == "tria3":
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

    elif set_shape["type"] == "tetr4":
        from myfempy.core.shapes.tetr4 import Tetra4
        return Tetra4

    else:
        pass


def setMesh(set_mesh):
    if set_mesh["TYPE"] == "add":
        from myfempy.core.mesh.mesh import MeshADD
        return MeshADD

    elif set_mesh["TYPE"] == "legacy":
        if set_mesh["shape"] == "quad4":
            from myfempy.core.mesh.legacyquad4 import LegacyQuad4
            return LegacyQuad4
        
        elif set_mesh["shape"] == "tria3":
            from myfempy.core.mesh.legacytria3 import LegacyTria3
            return LegacyTria3
        
        elif set_mesh["shape"] == "line2":
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
    if set_material["mat"] == "uniaxialstress":
        if set_material["type"] == "isotropic":
            from myfempy.core.material.uniaxialstress import UniAxialStressIsotropic
            
            return UniAxialStressIsotropic

    elif set_material["mat"] == "planestress":
        if set_material["type"] == "isotropic":
            from myfempy.core.material.planestress import PlaneStressIsotropic

            return PlaneStressIsotropic
        else:
            pass

    elif set_material["mat"] == "planestrain":
        if set_material["type"] == "isotropic":
            from myfempy.core.material.planestrain import PlaneStrainIsotropic

            return PlaneStrainIsotropic
        else:
            pass

    elif set_material["mat"] == "axisymmetric":
        pass

    elif set_material["mat"] == "solid":
        if set_material["type"] == "isotropic":
            from myfempy.core.material.solid import SolidIsotropic

            return SolidIsotropic
        else:
            pass
        
    elif set_material["mat"] == 'heatplane':
        if set_material["type"] == "isotropic":
            from myfempy.core.material.heatplane import HeatPlaneIsotropic
            
            return HeatPlaneIsotropic
            
    else:
        pass
    

def setGeometry(set_geometry):
    if set_geometry["geo"] == "thickness":
        from myfempy.core.geometry.thickness import Thickness
        return Thickness
   
    elif set_geometry["geo"] == "rectangle":
        from myfempy.core.geometry.rectangle import Rectangle
        return Rectangle 
    
    elif set_geometry["geo"] == "rectangle_tube":
        from myfempy.core.geometry.rectangle_tube import RectangleTube
        return RectangleTube 
    
    elif set_geometry["geo"] == "circle":
        from myfempy.core.geometry.circle import Circle
        return Circle 
    
    elif set_geometry["geo"] == "circle_tube":
        from myfempy.core.geometry.circle_tube import CircleTube
        return CircleTube 
    
    elif set_geometry["geo"] == "isection":
        from myfempy.core.geometry.isection import ISection
        return ISection 
    
    elif set_geometry["geo"] == "userdefined":
        from myfempy.core.geometry.userdefined import UserDefined
        return UserDefined 
    
    else:
        pass