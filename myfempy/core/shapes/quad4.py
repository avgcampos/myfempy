from __future__ import annotations

from myfempy.core.shapes.shape import Shape
from myfempy.core.shapes.quad4_tasks import ShapeFunctions, DiffShapeFuntion, Jacobian, invJacobi, detJacobi, NodeList, NodeCoord, LocKey

class Quad4(Shape):
    '''Quadrilateral 4-Node Shape Class <ConcreteClassService>'''
    
    def getShapeSet():
        shapeset = {
            "def": "4-nodes_conec 1-first_order",
            "key": "quad4",
            "id": 41,
            "nodes": ["i", "j", "k", "l"],
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
    