# -*- coding: utf-8 -*-
"""
Created on Sat May 16 12:05:12 2020
@author: ANTONIO VINICIUS GARCIA CAMPOS
_______________________________________________________________________________
 ~~~~~~~~~~          PROPRIEDADES DE SEÇÕES TRANSVERSAIS             ~~~~~~~~~~

> FUNCIONALIDADES
--- ENTRADAS: 
--- SAIDA: 

===============================================================================

> ATUALIZACOES DA VERSAO:

_______________________________________________________________________________
"""

import numpy as np

def sect_prop(typSection,dimSection):
    b = dimSection[0]
    h = dimSection[1]
    t = dimSection[2]
    d = dimSection[3]
    
    if typSection == 'rectangle':
        A = b*h      
        Izz = (1/12)*b*h**3
        Iyy = (1/12)*h*b**3 
        Jxx = Iyy+Izz
               
        y_max = h*0.5
        y_min = -h*0.5
        z_max = b*0.5
        z_min = -b*0.5
        r_max = y_max
        CGcoord = np.array([y_max,y_min,z_max,z_min,r_max])
        
        
    elif typSection == 'rectangle_tube':
        A = b*h - ((b-2*t)*(h-2*t))
        Izz = (1/12)*(b*h**3) - (1/12)*((b-2*t)*(h-2*t)**3)
        Iyy = (1/12)*(h*b**3) - (1/12)*((h-2*t)*(b-2*t)**3)
        Jxx = Iyy+Izz
        
        y_max = h*0.5
        y_min = -h*0.5
        z_max = b*0.5
        z_min = -b*0.5
        r_max = y_max
        CGcoord = np.array([y_max,y_min,z_max,z_min,r_max])

        
    elif typSection == 'circle':
        A = (1/4)*np.pi*d**2    
        Izz = (1/64)*np.pi*d**4    
        Iyy = Izz
        Jxx = Iyy + Izz
        
        y_max = d*0.5
        y_min = -d*0.5
        z_max = d*0.5
        z_min = -d*0.5
        r_max = y_max
        CGcoord = np.array([y_max,y_min,z_max,z_min,r_max])
                
    elif typSection == 'circle_tube':
        A = (1/4)*np.pi*(d**2 - (d-2*t)**2)
        Izz = (1/64)*np.pi*(d**4 - (d-2*t)**4)
        Iyy = Izz
        Jxx = Iyy + Izz
        
        y_max = d*0.5
        y_min = -d*0.5
        z_max = d*0.5
        z_min = -d*0.5
        r_max = y_max
        CGcoord = np.array([y_max,y_min,z_max,z_min,r_max])
        
        
    elif typSection == 'I_section':
        A = 2*b*d + t*(h-2*d)
        Izz = (b*h**3)/12 - ((b-t)*(h-2*d)**3)/12 
        Iyy = ((h-2*d)*t**3)/12 + 2*(d*b**3)/12
        Jxx = Iyy + Izz
        
        y_max = h*0.5
        y_min = -h*0.5
        z_max = b*0.5
        z_min = -b*0.5
        r_max = y_max
        CGcoord = np.array([y_max,y_min,z_max,z_min,r_max])
        
    elif typSection == 'spring':
        A = 1
        Izz = 1
        Iyy = 1
        Jxx = 1
        
        y_max = d*0.5
        y_min = -d*0.5
        z_max = d*0.5
        z_min = -d*0.5
        r_max = y_max
        CGcoord = np.array([y_max,y_min,z_max,z_min,r_max])
        
    return A,Izz,Iyy,Jxx,CGcoord

def area_coord(gmshgeo_file,typSection,dimSection,numPts):
    b = dimSection[0]
    h = dimSection[1]
    t = dimSection[2]
    d = dimSection[3]
    
    if typSection == 'circle':
        with open(gmshgeo_file,'w') as file_object:
            file_object.write('// GMSH GEOMETRY FILE FROM MYFEMPY *SECTION CIRCLE*\n')
            file_object.write('SetFactory("OpenCASCADE");\n')
            file_object.write('Circle(1) = {0, 0, 0,'+str(d)+', 0, 2*Pi};\n')
            
            file_object.write('Curve Loop(1) = {1};\n')
            file_object.write('Plane Surface(1) = {1};\n')
            file_object.write('Characteristic Length {1} = '+str(d/numPts)+';\n')
    
    elif typSection == 'circle_tube':
        with open(gmshgeo_file,'w') as file_object:
            file_object.write('// GMSH GEOMETRY FILE FROM MYFEMPY *SECTION CIRCLE TUBE*\n')
            file_object.write('SetFactory("OpenCASCADE");\n')
            file_object.write('Circle(1) = {0, 0, 0,'+str(d)+', 0, 2*Pi};\n')
            file_object.write('Circle(2) = {0, 0, 0,'+str(d-2*t)+', 0, 2*Pi};\n')
            
            file_object.write('Curve Loop(1) = {1};\n')
            file_object.write('Curve Loop(2) = {2};\n')
            file_object.write('Plane Surface(1) = {1,2};\n')
            file_object.write('Characteristic Length {1,2} = '+str(d/numPts)+';\n')
    
    if (typSection != 'circle')and(typSection != 'circle_tube'):
        
        if typSection == 'rectangle':
            numponlist = 4
            numlinlist = 4
            numlooplist = 1
            pointlist = np.array([[-b/2,-h/2,0],[b/2,-h/2,0],[b/2,h/2,0],[-b/2,h/2,0]])
            linelist = np.array([[1,2],[2,3],[3,4],[4,1]])
            plaLin_list = ['1,2,3,4']
            line_list = ['1,2,3,4']
            title = '// GMSH GEOMETRY FILE FROM MYFEMPY *SECTION RECTANGLE*\n'
            
        elif typSection == 'rectangle_tube':
            numponlist = 8
            numlinlist = 8
            numlooplist = 2
            pointlist = np.array([[-b/2,-h/2,0],[b/2,-h/2,0],[b/2,h/2,0],[-b/2,h/2,0],\
                                  [-b/2+t,-h/2+t,0],[b/2-t,-h/2+t,0],[b/2-t,h/2-t,0],[-b/2+t,h/2-t,0]])
            linelist = np.array([[1,2],[2,3],[3,4],[4,1],[5,6],[6,7],[7,8],[8,5]])
            plaLin_list = ['1,2,3,4','5,6,7,8']
            line_list = ['1,2,3,4','5,6,7,8']
            title = '// GMSH GEOMETRY FILE FROM MYFEMPY *SECTION RECTANGLE TUBE*\n'
        
        
        elif typSection == 'I_section':
            numponlist = 12
            numlinlist = 12
            numlooplist = 1
            pointlist = np.array([[-b/2,-h/2+d,0],[-b/2,-h/2,0],[b/2,-h/2,0],[b/2,-h/2+d,0],[t/2,-h/2+d,0],[t/2,h/2-d,0],\
                                  [b/2,h/2-d,0],[b/2,h/2,0],[-b/2,h/2,0],[-b/2,h/2-d,0],[-t/2,h/2-d,0],[-t/2,-h/2+d,0]])
            linelist = np.array([[1,2],[2,3],[3,4],[4,5],[5,6],[6,7],[7,8],[8,9],[9,10],[10,11],[11,12],[12,1]])
            plaLin_list = ['1,2,3,4,5,6,7,8,9,10,11,12']
            line_list = ['1,2,3,4,5,6,7,8,9,10,11,12']
            title = '// GMSH GEOMETRY FILE FROM MYFEMPY *I SECTION*\n'
    
    
        
        with open(gmshgeo_file,'w') as file_object:
            file_object.write(title)
            file_object.write('SetFactory("OpenCASCADE");\n')
            for i in range(0,numponlist):
                file_object.write('Point('+str(i+1)+') = {'+str(pointlist[i,0])+','+str(pointlist[i,1])+','+str(pointlist[i,2])+'};'+'\n')
        
            for i in range(0,numlinlist):
                file_object.write('Line('+str(i+1)+') = {'+str(linelist[i,0])+','+str(linelist[i,1])+'};\n')
            
            for i in range(0,numlooplist):
                file_object.write('Curve Loop('+str(i+1)+') = {'+plaLin_list[i]+'};\n')
                
            if numlooplist == 2:
                file_object.write('Plane Surface(1) = {1,2};\n')
            else: 
                file_object.write('Plane Surface(1) = {1};\n')
            
            for i in range(0,numlooplist):
                file_object.write('Transfinite Curve {'+line_list[i]+'} = '+str(numPts)+' Using Progression 1;\n')
            if numlooplist == 2:
                file_object.write('Transfinite Surface {1,2};\n')
            else: 
                file_object.write('Transfinite Surface {1};\n')

            
            
        
        
    