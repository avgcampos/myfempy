from __future__ import annotations

from myfempy.core.shapes.shape import Shape
from myfempy.core.shapes.tria3_tasks import ShapeFunctions, DiffShapeFuntion, Jacobian, invJacobi, detJacobi, NodeList, NodeCoord, LocKey

class Tria3(Shape):
    '''Triangular 3-Node Shape Class <ConcreteClassService>'''
    
    def getShapeSet():
        shapeset = {
            "def": "3-nodes_conec 1-first_order",
            "key": "tria3",
            "id": 31,
            "nodes": ["i", "j", "k"],
        }
        return shapeset
    
    
    def getShapeFunctions(r_coord, nodedof):
        return ShapeFunctions(r_coord, nodedof)


    def getDiffShapeFuntion(r_coord, nodedof):
        return DiffShapeFuntion(r_coord, nodedof)
        
         
    def getJacobian(r_coord, element_coord):  
        return Jacobian(r_coord, element_coord)
        

    def getinvJacobi(r_coord, element_coord, nodedof):
        return invJacobi(r_coord, element_coord, nodedof)
    
    
    def getdetJacobi(r_coord, element_coord):
        return detJacobi(r_coord, element_coord)
    
  
    def getNodeList(inci, element_number):            
        return NodeList(inci, element_number)
             
    def getNodeCoord(coord, node_list):
        return NodeCoord(coord, node_list)
    
 
    def getLocKey(node_list, nodedof):
        return LocKey(node_list, nodedof)