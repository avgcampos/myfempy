# -*- coding: utf-8 -*-
"""
========================================================================
~~~ MODULO DE SIMULACAO ESTRUTURAL PELO METODO DOS ELEMENTOS FINITOS ~~~
       	                    __                                
       	 _ __ ___   _   _  / _|  ___  _ __ ___   _ __   _   _ 
       	| '_ ` _ \ | | | || |_  / _ \| '_ ` _ \ | '_ \ | | | |
       	| | | | | || |_| ||  _||  __/| | | | | || |_) || |_| |
       	|_| |_| |_| \__, ||_|   \___||_| |_| |_|| .__/  \__, |
       	            |___/                       |_|     |___/ 

~~~      Mechanical studY with Finite Element Method in PYthon       ~~~
~~~                PROGRAMA DE ANÁLISE COMPUTACIONAL                 ~~~
~~~              copyright @ 2022, all rights reserved               ~~~
========================================================================
"""

#%% READ INPUTDATA FROM USER PATH
import sys
import os
# import imp
import numpy as np
# from colorama import Fore, Back, Style

# from myfempy.felib.felib import fe_key
# from myfempy.felib.material import mat_def, mat_beh
# from myfempy.felib.crossec import sect_prop, sec_def


#%% 

def gmsh_key(meshtype):
    
    l = {'line2':'-1',
         'tria3': '-2',
         'quad4':'-2',
         'hexa8':'-3',
         'tetr4':'-3',
        }
    
    return l[meshtype]


# gera a malha no GMSH, inserir arquivo de geometria .geo --> run gmsh 
def get_gmsh_msh(meshdata):
    # cmd = 'gmsh geo_1d_beam.geo -1 -o mesh_out_1d_beam.msh1'
    
    cmd = "gmsh"+" "+(meshdata['GMSH']['filename']+'.geo')+" "+gmsh_key(meshdata['GMSH']['meshconfig']['mesh'])+" -o "+(meshdata['GMSH']['filename']+'.msh1') #+'.msh1'
        
    os.system("echo GENERATING MESH FROM EXTERNAL GMSH")
    os.system(cmd)
    os.system("echo MESHING IS DONE")
    os.system("echo SAVING AND EXIT")


#-----------------------------------------------------------------------------#
def get_gmsh_geo(meshdata):
    with open((meshdata['GMSH']['filename']+'.geo'),'w') as file_object:
        
        file_object.write('// GMSH GEOMETRY FILE FROM MYFEMPY\n')
        file_object.write('SetFactory("OpenCASCADE");\n')
        if 'pointlist' in meshdata['GMSH'].keys():
            
        
            # if 'linelist' in meshdata['GMSH'].keys():
            numlinlist = len(meshdata['GMSH']['linelist'])
        
            line_list = ""
            for i in range(numlinlist):
                line_list += str(i+1)+','
            line_list = line_list[0:-1]
            
            if 'planelist' in meshdata['GMSH'].keys():
                numplalistP = len(meshdata['GMSH']['planelist'])
                planes = ""
                for i in range(numplalistP):
                    planes += str(i+1)+','
                planes = planes[0:-1]
            else:
                pass
            
            numpnt = len(meshdata['GMSH']['pointlist'])
            for i in range(0,numpnt):
                file_object.write('Point('+str(i+1)+') = {'+str(meshdata['GMSH']['pointlist'][i][0])+','+str(meshdata['GMSH']['pointlist'][i][1])+','+str(meshdata['GMSH']['pointlist'][i][2])+','+str(meshdata['GMSH']['meshconfig']['sizeelement'])+'};'+'\n')
            
                
            if 'arc' in meshdata['GMSH'].keys():
                numincl = len(meshdata['GMSH']['arc'])
                for inl in range(numincl):
                    d = meshdata['GMSH']['arc'][inl][0]
                    cx = meshdata['GMSH']['arc'][inl][1][0]
                    cy = meshdata['GMSH']['arc'][inl][1][1]
                    cz = meshdata['GMSH']['arc'][inl][1][2]
                    arc0 = meshdata['GMSH']['arc'][inl][2][0]
                    arc1 = meshdata['GMSH']['arc'][inl][2][1]
                    file_object.write('Circle('+str(numlinlist+inl+1)+') = {'+str(cx)+','+str(cy)+','+str(cz)+','+str(d)+','+arc0+','+arc1+'};\n')
                    
                    # file_object.write('Circle('+str(numlinlist+inl+1)+') = {'+str(meshdata['GMSH']['arc'][inl][0])+','+str(meshdata['GMSH']['arc'][inl][1])+','+str(meshdata['GMSH']['arc'][inl][2])+'};\n')
    
                    # if meshdata['GMSH']['curve'][inl][0] == 'hole':
                    #     d = meshdata['GMSH']['curve'][inl][2]
                    #     cx = meshdata['GMSH']['curve'][inl][3][0]
                    #     cy = meshdata['GMSH']['curve'][inl][3][1]
                    #     cz = meshdata['GMSH']['curve'][inl][3][2]
                    #     file_object.write('Circle('+str(numlinlist+inl+1)+') = {'+str(cx)+','+str(cy)+','+str(cz)+','+str(d)+',0,2*Pi};\n')
                    
                    # elif meshdata['GMSH']['curve'][inl][0] == 'arc':
            else:
                numincl = 0
                
            for i in range(0,numlinlist):
                file_object.write('Line('+str(i+1)+') = {'+str(meshdata['GMSH']['linelist'][i][0])+','+str(meshdata['GMSH']['linelist'][i][1])+'};\n')
                   
                    # numplalistP += 1
                    # plane_list  += ','+str(numplalistP+inl+1)
                
            # if 'hole' in meshdata['GEOMETRY'].keys():
            #     numincl = len(meshdata['GEOMETRY']['hole'])
            #     for inl in range(numincl):
            #         d = meshdata['GEOMETRY']['hole'][inl][0]
            #         sx1 = meshdata['GEOMETRY']['hole'][inl][1][0]+1.5*d
            #         sy1 = meshdata['GEOMETRY']['hole'][inl][1][1]+1.5*d
            #         sx2 = meshdata['GEOMETRY']['hole'][inl][1][0]-1.5*d
            #         sy2 = meshdata['GEOMETRY']['hole'][inl][1][1]+1.5*d
            #         sx3 = meshdata['GEOMETRY']['hole'][inl][1][0]-1.5*d
            #         sy3 = meshdata['GEOMETRY']['hole'][inl][1][1]-1.5*d
            #         sx4 = meshdata['GEOMETRY']['hole'][inl][1][0]+1.5*d
            #         sy4 = meshdata['GEOMETRY']['hole'][inl][1][1]-1.5*d
                    
            #         cx1 = meshdata['GEOMETRY']['hole'][inl][1][0]+d*np.cos(np.pi/4)
            #         cy1 = meshdata['GEOMETRY']['hole'][inl][1][1]+d*np.cos(np.pi/4)
            #         cx2 = meshdata['GEOMETRY']['hole'][inl][1][0]-d*np.cos(np.pi/4)
            #         cy2 = meshdata['GEOMETRY']['hole'][inl][1][1]+d*np.sin(np.pi/4)
            #         cx3 = meshdata['GEOMETRY']['hole'][inl][1][0]-d*np.cos(np.pi/4)
            #         cy3 = meshdata['GEOMETRY']['hole'][inl][1][1]-d*np.sin(np.pi/4)
            #         cx4 = meshdata['GEOMETRY']['hole'][inl][1][0]+d*np.cos(np.pi/4)
            #         cy4 = meshdata['GEOMETRY']['hole'][inl][1][1]-d*np.sin(np.pi/4)
                    
            #         file_object.write('Point('+str(numpnt+inl+1)+') = {'+str(meshdata['GEOMETRY']['hole'][inl][1][0])+','+str(meshdata['GEOMETRY']['hole'][inl][1][1])+','+str(meshdata['GEOMETRY']['hole'][inl][1][2])+'};\n')
            #         file_object.write('Point('+str(numpnt+inl+2)+') = {'+str(cx1)+','+str(cy1)+','+str(meshdata['GEOMETRY']['hole'][inl][1][2])+'};\n')
            #         file_object.write('Point('+str(numpnt+inl+3)+') = {'+str(cx2)+','+str(cy2)+','+str(meshdata['GEOMETRY']['hole'][inl][1][2])+'};\n')
            #         file_object.write('Point('+str(numpnt+inl+4)+') = {'+str(cx3)+','+str(cy3)+','+str(meshdata['GEOMETRY']['hole'][inl][1][2])+'};\n')
            #         file_object.write('Point('+str(numpnt+inl+5)+') = {'+str(cx4)+','+str(cy4)+','+str(meshdata['GEOMETRY']['hole'][inl][1][2])+'};\n')
            #         file_object.write('Point('+str(numpnt+inl+6)+') = {'+str(sx1)+','+str(sy1)+','+str(meshdata['GEOMETRY']['hole'][inl][1][2])+'};\n')
            #         file_object.write('Point('+str(numpnt+inl+7)+') = {'+str(sx2)+','+str(sy2)+','+str(meshdata['GEOMETRY']['hole'][inl][1][2])+'};\n')
            #         file_object.write('Point('+str(numpnt+inl+8)+') = {'+str(sx3)+','+str(sy3)+','+str(meshdata['GEOMETRY']['hole'][inl][1][2])+'};\n')
            #         file_object.write('Point('+str(numpnt+inl+9)+') = {'+str(sx4)+','+str(sy4)+','+str(meshdata['GEOMETRY']['hole'][inl][1][2])+'};\n')
                    
            #         file_object.write('Circle('+str(numlinlist+inl+1)+') = {'+str(numpnt+inl+2)+','+str(numpnt+inl+1)+','+str(numpnt+inl+3)+'};\n')
            #         file_object.write('Circle('+str(numlinlist+inl+2)+') = {'+str(numpnt+inl+3)+','+str(numpnt+inl+1)+','+str(numpnt+inl+4)+'};\n')
            #         file_object.write('Circle('+str(numlinlist+inl+3)+') = {'+str(numpnt+inl+4)+','+str(numpnt+inl+1)+','+str(numpnt+inl+5)+'};\n')
            #         file_object.write('Circle('+str(numlinlist+inl+4)+') = {'+str(numpnt+inl+5)+','+str(numpnt+inl+1)+','+str(numpnt+inl+2)+'};\n')
                    
            #         file_object.write('Line('+str(numlinlist+inl+5)+') = {'+str(numpnt+inl+6)+','+str(numpnt+inl+7)+'};\n')
            #         file_object.write('Line('+str(numlinlist+inl+6)+') = {'+str(numpnt+inl+7)+','+str(numpnt+inl+8)+'};\n')
            #         file_object.write('Line('+str(numlinlist+inl+7)+') = {'+str(numpnt+inl+8)+','+str(numpnt+inl+9)+'};\n')
            #         file_object.write('Line('+str(numlinlist+inl+8)+') = {'+str(numpnt+inl+9)+','+str(numpnt+inl+6)+'};\n')
            #         file_object.write('Line('+str(numlinlist+inl+9)+') = {'+str(numpnt+inl+2)+','+str(numpnt+inl+6)+'};\n')
            #         file_object.write('Line('+str(numlinlist+inl+10)+') = {'+str(numpnt+inl+3)+','+str(numpnt+inl+7)+'};\n')
            #         file_object.write('Line('+str(numlinlist+inl+11)+') = {'+str(numpnt+inl+4)+','+str(numpnt+inl+8)+'};\n')
            #         file_object.write('Line('+str(numlinlist+inl+12)+') = {'+str(numpnt+inl+5)+','+str(numpnt+inl+9)+'};\n')
            #         plane_list  += ','+str(numplalistP+inl+1)
                
        #---------------------------
        if meshdata['GMSH']['meshconfig']['mesh'] == 'line2':           
            
            if 'numbernodes' in meshdata['GMSH']['meshconfig'].keys():
                file_object.write('Transfinite Curve {'+line_list+'} = '+str(meshdata['GMSH']['meshconfig']['numbernodes'])+' Using Progression 1;\n')
                
            # file_object.write('Characteristic Length {:} = '+str(meshdata['GMSH']['meshconfig']['sizeelement'])+';\n')
            
            elif 'sizeelement' in meshdata['GMSH']['meshconfig'].keys():
                for i in range(0,numpnt):
                    file_object.write('Point('+str(i+1)+') = {'+str(meshdata['GMSH']['pointlist'][i][0])+','+str(meshdata['GMSH']['pointlist'][i][1])+','+str(meshdata['GMSH']['pointlist'][i][2])+','+str(meshdata['GMSH']['meshconfig']['sizeelement'])+'};'+'\n')
            
            else:
                pass
        #---------------------------
        elif meshdata['GMSH']['meshconfig']['mesh'] == 'tria3':
            
            if 'cadimport' in meshdata['GMSH'].keys():
                file_object.write('Merge "'+meshdata['GMSH']['cadimport']['object']+'";\n')   
            
            else:
                npl = 0
                phl = 0
                for i in range(0,numplalistP):
                    npl+=1
                    file_object.write('Curve Loop('+str(npl)+') = {'+(str(meshdata['GMSH']['planelist'][i][:]))[1:-1]+'};\n')
                    plane_hole = ""
                    plane_hole += str(npl)+','
                    for inl in range(numincl):
                        if meshdata['GMSH']['arc'][inl][0] == i+1:
                            npl+=1
                            plane_hole += str(npl)+','
                            phl=meshdata['GMSH']['arc'][inl][1]
                            file_object.write('Curve Loop('+str(i+1+inl+1)+') = {'+str(numlinlist+inl+1)+'};\n')
                                                
                    if 'arc' in meshdata['GMSH'].keys():
                        file_object.write('Plane Surface('+str(i+1)+') = {'+plane_hole[0:-1]+'};\n')
                        
                    else:
                        file_object.write('Plane Surface('+str(i+1)+') = {'+str(npl)+'};\n')

                            
            file_object.write('Characteristic Length {:} = '+str(meshdata['GMSH']['meshconfig']['sizeelement'])+';\n')

            if meshdata['GMSH']['meshconfig']['meshmap'] == True:
                file_object.write('//FACE MAPPING \n')
                # file_object.write('Transfinite Curve {:} = '+str(meshdata['GMSH']['meshconfig']['numbernodes'])+' Using Progression 1;\n')
                file_object.write('Transfinite Surface {:};\n')

            else:
                pass
            
            
            file_object.write('// MESH CONFIGURATION\n')
            file_object.write('Mesh.CharacteristicLengthExtendFromBoundary = 1;\n')
            file_object.write('Mesh.CharacteristicLengthMin = 0;\n')
            file_object.write('Mesh.CharacteristicLengthMax = '+str(meshdata['GMSH']['meshconfig']['sizeelement'])+';\n')
            file_object.write('Mesh.CharacteristicLengthFromPoints = 1;\n')
            # file_object.write('Mesh.CharacteristicLengthFromCurvature = 1;\n')
            # file_object.write('Mesh.MinimumCurvePoints = 12;\n')
            # file_object.write('Mesh.MinimumElementsPerTwoPi = 12;\n')
            # file_object.write('Mesh.CharacteristicLengthExtendFromBoundary = 1;\n')
            file_object.write('Mesh.Optimize = 1;\n')
            file_object.write('Mesh.HighOrderOptimize = 0;\n')
            file_object.write('Mesh.Algorithm = 8;\n')
            file_object.write('Mesh.ElementOrder = 1;\n')

                
        #---------------------------
        elif meshdata['GMSH']['meshconfig']['mesh'] == 'quad4':
            
            if 'cadimport' in meshdata['GMSH'].keys():
                file_object.write('Merge "'+meshdata['GMSH']['cadimport']['object']+'";\n')   
            
            else:
                npl = 0
                phl = 0
                for i in range(0,numplalistP):
                    npl+=1
                    file_object.write('Curve Loop('+str(npl)+') = {'+(str(meshdata['GMSH']['planelist'][i][:]))[1:-1]+'};\n')
                    plane_hole = ""
                    plane_hole += str(npl)+','
                    for inl in range(numincl):
                        if meshdata['GMSH']['arc'][inl][0] == i+1:
                            npl+=1
                            plane_hole += str(npl)+','
                            phl=meshdata['GMSH']['arc'][inl][1]
                            file_object.write('Curve Loop('+str(i+1+inl+1)+') = {'+str(numlinlist+inl+1)+'};\n')
                                                
                    if 'arc' in meshdata['GMSH'].keys():
                        file_object.write('Plane Surface('+str(i+1)+') = {'+plane_hole[0:-1]+'};\n')
                        
                    else:
                        file_object.write('Plane Surface('+str(i+1)+') = {'+str(npl)+'};\n')

                            
            file_object.write('Characteristic Length {:} = '+str(meshdata['GMSH']['meshconfig']['sizeelement'])+';\n')

            if meshdata['GMSH']['meshconfig']['meshmap'] == True:
                file_object.write('//FACE MAPPING \n')
                # file_object.write('Transfinite Curve {:} = '+str(meshdata['GMSH']['meshconfig']['numbernodes'])+' Using Progression 1;\n')
                file_object.write('Transfinite Surface {:};\n')

            else:
                pass
            
                file_object.write('Recombine Surface {:};\n')
                file_object.write('// MESH CONFIGURATION\n')
                file_object.write('Mesh.RecombinationAlgorithm = 1;\n')
                file_object.write('Mesh.RecombineAll = 1;\n')
                file_object.write('Mesh.SubdivisionAlgorithm = 1;\n')
                file_object.write('Mesh.CharacteristicLengthExtendFromBoundary = 1;\n')
                file_object.write('Mesh.CharacteristicLengthMin = 0;\n')
                file_object.write('Mesh.CharacteristicLengthMax = '+str(meshdata['GMSH']['meshconfig']['sizeelement'])+';\n')
                file_object.write('Mesh.CharacteristicLengthFromPoints = 1;\n')
                # file_object.write('Mesh.CharacteristicLengthFromCurvature = 1;\n')
                # file_object.write('Mesh.MinimumCurvePoints = 12;\n')
                # file_object.write('Mesh.MinimumElementsPerTwoPi = 12;\n')
                # file_object.write('Mesh.CharacteristicLengthExtendFromBoundary = 1;\n')
                file_object.write('Mesh.Optimize = 1;\n')
                file_object.write('Mesh.HighOrderOptimize = 0;\n')
                file_object.write('Mesh.Algorithm = 8;\n')
                file_object.write('Mesh.ElementOrder = 1;\n')
            
            
        #---------------------------    
        elif meshdata['GMSH']['meshconfig']['mesh'] == 'hexa8':
            thck = 1
            if 'cadimport' in meshdata['GMSH'].keys():
                file_object.write('Merge "'+meshdata['GMSH']['cadimport']['object']+'";\n')   
            
            else:
                npl = 0
                phl = 0
                for i in range(0,numplalistP):
                    npl+=1
                    file_object.write('Curve Loop('+str(npl)+') = {'+(str(meshdata['GMSH']['planelist'][i][:]))[1:-1]+'};\n')
                    plane_hole = ""
                    plane_hole += str(npl)+','
                    for inl in range(numincl):
                        if meshdata['GMSH']['arc'][inl][0] == i+1:
                            npl+=1
                            plane_hole += str(npl)+','
                            phl=meshdata['GMSH']['arc'][inl][1]
                            file_object.write('Curve Loop('+str(i+1+inl+1)+') = {'+str(numlinlist+inl+1)+'};\n')
                                                
                    if 'arc' in meshdata['GMSH'].keys():
                        file_object.write('Plane Surface('+str(i+1)+') = {'+plane_hole[0:-1]+'};\n')
                        
                    else:
                        file_object.write('Plane Surface('+str(i+1)+') = {'+str(npl)+'};\n')

                thck = meshdata['GMSH']['meshconfig']['extrude']
                file_object.write('Extrude {0, 0, '+str(float(thck))+'} {Surface{:};Layers{'+str(int(float(thck)/float(meshdata['GMSH']['meshconfig']['sizeelement'])))+'};Recombine;};\n')
 
            file_object.write('Characteristic Length {:} = '+str(meshdata['GMSH']['meshconfig']['sizeelement'])+';\n')

            if meshdata['GMSH']['meshconfig']['meshmap'] == True:
                file_object.write('//FACE MAPPING \n')
                # file_object.write('Transfinite Curve {:} = '+str(meshdata['GMSH']['meshconfig']['numbernodes'])+' Using Progression 1;\n')
                file_object.write('Transfinite Surface {:};\n')

            else:
                pass

                file_object.write('Recombine Surface {:};\n')
                file_object.write('// MESH CONFIGURATION\n')
                file_object.write('Mesh.Algorithm =8;\n')
                file_object.write('Mesh.Algorithm3D = 7;\n')
                file_object.write('Mesh.CharacteristicLengthMin = 0;\n')
                file_object.write('Mesh.CharacteristicLengthMax = '+str(meshdata['GMSH']['meshconfig']['sizeelement'])+';\n')
                file_object.write('Mesh.ElementOrder = 1;\n')
                # file_object.write('Mesh.CharacteristicLengthFromPoints = 1;\n')
                # file_object.write('Mesh.CharacteristicLengthFromCurvature = 1;\n')
                # file_object.write('Mesh.MinimumCurvePoints = 12;\n')
                # file_object.write('Mesh.MinimumElementsPerTwoPi = 12;\n')
                # file_object.write('Mesh.CharacteristicLengthExtendFromBoundary = 1;\n')
                file_object.write('Mesh.Optimize = 1;\n')
                file_object.write('Mesh.HighOrderOptimize = 0;\n')
                file_object.write('Mesh.RecombinationAlgorithm = 0;\n')
                file_object.write('Mesh.SubdivisionAlgorithm = 2;\n')
                file_object.write('Mesh.RecombineAll = 1;\n')
                
        else:
            print("input erro: mesh_cfg don't defined")
                    
            


# #-----------------------------------------------------------------------------#
# def get_gmsh_geo_cs(gmshgeo_file,sec_def,dimSection,numPts):
#     b = dimSection[0]
#     h = dimSection[1]
#     t = dimSection[2]
#     d = dimSection[3]
    
#     if sec_def == 'circle':
#         with open(gmshgeo_file,'w') as file_object:
#             file_object.write('// GMSH GEOMETRY FILE FROM MYFEMPY *SECTION CIRCLE*\n')
#             file_object.write('SetFactory("OpenCASCADE");\n')
#             file_object.write('Circle(1) = {0, 0, 0,'+str(d)+', 0, 2*Pi};\n')
            
#             file_object.write('Curve Loop(1) = {1};\n')
#             file_object.write('Plane Surface(1) = {1};\n')
#             file_object.write('Characteristic Length {1} = '+str(d/numPts)+';\n')
    
#     elif sec_def == 'circle_tube':
#         with open(gmshgeo_file,'w') as file_object:
#             file_object.write('// GMSH GEOMETRY FILE FROM MYFEMPY *SECTION CIRCLE TUBE*\n')
#             file_object.write('SetFactory("OpenCASCADE");\n')
#             file_object.write('Circle(1) = {0, 0, 0,'+str(d)+', 0, 2*Pi};\n')
#             file_object.write('Circle(2) = {0, 0, 0,'+str(d-2*t)+', 0, 2*Pi};\n')
            
#             file_object.write('Curve Loop(1) = {1};\n')
#             file_object.write('Curve Loop(2) = {2};\n')
#             file_object.write('Plane Surface(1) = {1,2};\n')
#             file_object.write('Characteristic Length {1,2} = '+str(d/numPts)+';\n')
    
#     if (sec_def != 'circle')and(sec_def != 'circle_tube'):
        
#         if sec_def == 'rectangle':
#             numponlist = 4
#             numlinlist = 4
#             numlooplist = 1
#             pointlist = np.array([[-b/2,-h/2,0],[b/2,-h/2,0],[b/2,h/2,0],[-b/2,h/2,0]])
#             linelist = np.array([[1,2],[2,3],[3,4],[4,1]])
#             plaLin_list = ['1,2,3,4']
#             line_list = ['1,2,3,4']
#             title = '// GMSH GEOMETRY FILE FROM MYFEMPY *SECTION RECTANGLE*\n'
            
#         elif sec_def == 'rectangle_tube':
#             numponlist = 8
#             numlinlist = 8
#             numlooplist = 2
#             pointlist = np.array([[-b/2,-h/2,0],[b/2,-h/2,0],[b/2,h/2,0],[-b/2,h/2,0],\
#                                   [-b/2+t,-h/2+t,0],[b/2-t,-h/2+t,0],[b/2-t,h/2-t,0],[-b/2+t,h/2-t,0]])
#             linelist = np.array([[1,2],[2,3],[3,4],[4,1],[5,6],[6,7],[7,8],[8,5]])
#             plaLin_list = ['1,2,3,4','5,6,7,8']
#             line_list = ['1,2,3,4','5,6,7,8']
#             title = '// GMSH GEOMETRY FILE FROM MYFEMPY *SECTION RECTANGLE TUBE*\n'
        
        
#         elif sec_def == 'I_section':
#             numponlist = 12
#             numlinlist = 12
#             numlooplist = 1
#             pointlist = np.array([[-b/2,-h/2+d,0],[-b/2,-h/2,0],[b/2,-h/2,0],[b/2,-h/2+d,0],[t/2,-h/2+d,0],[t/2,h/2-d,0],\
#                                   [b/2,h/2-d,0],[b/2,h/2,0],[-b/2,h/2,0],[-b/2,h/2-d,0],[-t/2,h/2-d,0],[-t/2,-h/2+d,0]])
#             linelist = np.array([[1,2],[2,3],[3,4],[4,5],[5,6],[6,7],[7,8],[8,9],[9,10],[10,11],[11,12],[12,1]])
#             plaLin_list = ['1,2,3,4,5,6,7,8,9,10,11,12']
#             line_list = ['1,2,3,4,5,6,7,8,9,10,11,12']
#             title = '// GMSH GEOMETRY FILE FROM MYFEMPY *I SECTION*\n'
    
    
        
#         with open(gmshgeo_file,'w') as file_object:
#             file_object.write(title)
#             file_object.write('SetFactory("OpenCASCADE");\n')
#             for i in range(0,numponlist):
#                 file_object.write('Point('+str(i+1)+') = {'+str(pointlist[i,0])+','+str(pointlist[i,1])+','+str(pointlist[i,2])+'};'+'\n')
        
#             for i in range(0,numlinlist):
#                 file_object.write('Line('+str(i+1)+') = {'+str(linelist[i,0])+','+str(linelist[i,1])+'};\n')
            
#             for i in range(0,numlooplist):
#                 file_object.write('Curve Loop('+str(i+1)+') = {'+plaLin_list[i]+'};\n')
                
#             if numlooplist == 2:
#                 file_object.write('Plane Surface(1) = {1,2};\n')
#             else: 
#                 file_object.write('Plane Surface(1) = {1};\n')
            
#             for i in range(0,numlooplist):
#                 file_object.write('Transfinite Curve {'+line_list[i]+'} = '+str(numPts)+' Using Progression 1;\n')
#             if numlooplist == 2:
#                 file_object.write('Transfinite Surface {1,2};\n')
#             else: 
#                 file_object.write('Transfinite Surface {1};\n')

            
            

