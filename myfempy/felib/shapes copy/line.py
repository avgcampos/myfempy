from __future__ import annotations

import numdifftools as nd
import numpy as np
from shape import Shape


class Line(Shape):
    '''Line Shape Class <ConcreteClassService>'''
    
    def getShapeFunctions(poly_order, r_coord):
        try: 
            
            if poly_order == '1':
                N = np.zeros((1,2))
                N[0, 0] = 0.5*(1.0-r_coord)
                N[0, 1] = 0.5*(1.0+r_coord)
            
            elif poly_order == '2':
                N = np.zeros((1,3))
                N[0, 0] = 0.5*(1.0-r_coord)-0.5*(1.0-r_coord*r_coord) 
                N[0, 1] = 1-r_coord*r_coord
                N[0, 2] = 0.5*(1.0+r_coord)-0.5*(1.0-r_coord*r_coord)
            
            return N
        
        except:
            print('Error Line Shape Function')
        
    def getDiffShapeFuntion(shape_function, r_coord, n):
        try:
            
            dN = nd.Gradient(shape_function, n = n)
            
            return dN(r_coord)
    
        except:
            print('Error Line DiffShape Function')
            
    def getJacobian(diff_shape_function, element_coord):    
        try:
            
            return np.dot(diff_shape_function, element_coord)
        
        except:
            print('Error Line Jacobian Function')
            
        