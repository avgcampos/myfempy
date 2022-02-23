"""   
_______________________________________________________________________________
~~~~~~ MODULO DE SIMULACAO ESTRUTURAL PELO METODO DOS ELEMENTOS FINITOS ~~~~~~~
~~~~~~                   PROGRAMA DE PROPOSITO GERAL                    ~~~~~~~  
~~~~~~                  copyright all rights reserved                   ~~~~~~~    
===============================================================================
BSD 3-Clause License
Copyright (c) 2021, 3D EasyCAE
All rights reserved.
"""

#%% READ INPUTDATA FROM USER PATH
import sys
import os
import imp
import numpy as np
from colorama import Fore, Back, Style

#%% 
# gera a malha no GMSH, inserir arquivo de geometria .geo --> run gmsh 
def gmsh_meshing(file_geo,name_out_mesh,geo_select):
    # cmd = 'gmsh geo_1d_beam.geo -1 -o mesh_out_1d_beam.msh1'
    if geo_select == 'beam':
        cmd = "gmsh"+" "+file_geo+" -1 -o "+name_out_mesh
        
    if geo_select == 'plate':
        cmd = "gmsh"+" "+file_geo+" -2 -o "+name_out_mesh
    
    if geo_select == 'solid':
        cmd = "gmsh"+" "+file_geo+" -3 -o "+name_out_mesh
    
    os.system("echo GENERATING MESH FILE WITH GMSH")
    os.system(cmd)
    os.system("echo MESHING IS DONE")
    os.system("echo SAVING AND EXIT")
              
#%%
def gmsh_geo2msh(gmshgeo_file,pointlist,linelist,planelist,meshconfig,propgeonewlist):
    with open(gmshgeo_file,'w') as file_object:
        numpnt = len(pointlist)
        numlinlist = len(linelist)
        numplalistP = len(planelist)
        # numplalistL = (planelist).shape[1]
        
        line_list = ""
        for i in range(numlinlist):
            line_list += str(i+1)+','
        line_list = line_list[0:-1]
        
        plane_list = ""
        for i in range(numplalistP):
            plane_list += str(i+1)+','
        plane_list = plane_list[0:-1]
        
        for i in range(numplalistP):
            numplalistL = len(planelist[i])
            exec(f'plaLin_list_{i} = ""')
            for j in range(numplalistL):
                exec(f"plaLin_list_{i} += planelist[i][j]+','")
        
        plaLin_list = []
        for i in range(numplalistP):
            exec(f"plaLin_list_{i} = plaLin_list_{i}[0:-1]")
            exec(f"plaLin_list.append(plaLin_list_{i})")       
                
        file_object.write('// GMSH GEOMETRY FILE FROM MYFEMPY\n')
        for i in range(0,numpnt):
            file_object.write('Point('+str(i+1)+') = {'+str(pointlist[i,1])+','+str(pointlist[i,2])+','+str(pointlist[i,3])+'};'+'\n')
        
        for i in range(0,numlinlist):
            file_object.write('Line('+str(i+1)+') = {'+str(linelist[i,1])+','+str(linelist[i,2])+'};\n')
        
        if meshconfig[0,1] == 'line2node':
            if meshconfig[0,2] == 'sizeelement':
                file_object.write('Characteristic Length {'+line_list+'} = '+meshconfig[0,3]+';\n')
            elif meshconfig[0,2] == 'numbernodes':
                file_object.write('Transfinite Curve {'+line_list+'} = '+meshconfig[0,3]+' Using Progression 1;\n')
            else:
                print("input erro: mesh_cfg don't defined")
        
        elif meshconfig[0,1] == 'tria3node':
            for i in range(0,numplalistP):
                file_object.write('Curve Loop('+str(i+1)+') = {'+plaLin_list[i]+'};\n')
            file_object.write('Plane Surface(1) = {'+plane_list+'};\n')
            if meshconfig[0,2] == 'sizeelement':
                file_object.write('Characteristic Length {'+line_list+'} = '+meshconfig[0,3]+';\n')
            elif meshconfig[0,2] == 'numbernodes':
                file_object.write('Transfinite Curve {'+line_list+'} = '+meshconfig[0,3]+' Using Progression 1;\n')
            else:
                print("input erro: mesh_cfg don't defined")
            if meshconfig[0,4] == 'facemap':
                file_object.write('//FACE MAPPING \n')
                file_object.write('Transfinite Surface {1};\n')
            elif meshconfig[0,4] == 'nofacemap':
                file_object.write('//NOT FACE MAPPING \n')
            else:
                print("input erro: mesh_cfg don't defined")
        
        elif meshconfig[0,1] == 'quad4node':
            for i in range(0,numplalistP):
                file_object.write('Curve Loop('+str(i+1)+') = {'+plaLin_list[i]+'};\n')
            file_object.write('Plane Surface(1) = {'+plane_list+'};\n')
            if meshconfig[0,2] == 'sizeelement':
                file_object.write('Characteristic Length {'+line_list+'} = '+meshconfig[0,3]+';\n')
            elif meshconfig[0,2] == 'numbernodes':
                file_object.write('Transfinite Curve {'+line_list+'} = '+meshconfig[0,3]+' Using Progression 1;\n')
            else:
                print("input erro: mesh_cfg don't defined")
            if meshconfig[0,4] == 'facemap':
                file_object.write('//FACE MAPPING \n')
                file_object.write('Transfinite Surface {1};\n')
            elif meshconfig[0,4] == 'nofacemap':
                file_object.write('//NOT FACE MAPPING \n')
            else:
                print("input erro: mesh_cfg don't defined")
            file_object.write('Mesh.RecombineAll = 1;\n')
            
        elif meshconfig[0,1] == 'hexa8node':
            thck = propgeonewlist[0,5]
            for i in range(0,numplalistP):
                file_object.write('Curve Loop('+str(i+1)+') = {'+plaLin_list[i]+'};\n')
            file_object.write('Plane Surface(1) = {'+plane_list+'};\n')
            file_object.write('Characteristic Length {'+line_list+'} = '+meshconfig[0,3]+';\n')
            if meshconfig[0,4] == 'extrude':
                file_object.write('//FACE MAPPING AND EXTRUDE SURFACE\n')
                file_object.write('Transfinite Surface {1};\n')
                file_object.write('Recombine Surface {1};\n')
                file_object.write('Extrude {0, 0, '+thck+'} {Surface{1};Layers{'+str(int(float(thck)/float(meshconfig[0,3])))+'};Recombine;};\n')
            elif meshconfig[0,4] == 'nofacemap':
                file_object.write('//NOT FACE MAPPING \n')
            else:
                print("input erro: mesh_cfg don't defined")
        else:
            print("input erro: mesh_cfg don't defined")
            
            
#%% IMPORT MESH FROM GMSH, VERSION NEW
#version old 
# # type_elm == 1 -> 2-node line
def import_mesh_line(file_imp):
    with open(file_imp,'r') as file_object:
        file_object.readline()
        NNOD = int(file_object.readline())
        coord = np.zeros((NNOD,4))
        for ii in range(0,NNOD):
            line = file_object.readline()
            lineaux = line.split()
            contstr = lineaux[0:4]
            coord[ii,:] = [float(contstr[0]),float(contstr[1]),float(contstr[2]),float(contstr[3])]
        
        file_object.readline()
        file_object.readline()
        NELM = int(file_object.readline())
        conec_elm = np.zeros((NELM,9))
        prop_elm = np.zeros((NELM,3))
        contelm=0
        for kk in range(0,NELM):
            line = file_object.readline()
            lineaux = line.split()
            type_elm = int(lineaux[1])     
            if type_elm == 1:
                contelm+=1
                prop_elm[contelm-1,:] = [type_elm,float(lineaux[3]),float(lineaux[3])]
                conec_elm[contelm-1,:] = [contelm,float(lineaux[5]),float(lineaux[6]),0,0,0,0,0,0]
            
    conec_elm = conec_elm[:contelm,:]
    prop_elm = prop_elm[:contelm,:]
    inci = np.concatenate((conec_elm[:,0][:, np.newaxis],prop_elm,conec_elm[:,1:]),axis=1)
    return coord, inci, type_elm


# type_elm == 2 -> 3-node triangle
def import_mesh_tria(file_imp):
    with open(file_imp,'r') as file_object:
        file_object.readline()
        NNOD = int(file_object.readline())
        coord = np.zeros((NNOD,4))
        for ii in range(0,NNOD):
            line = file_object.readline()
            lineaux = line.split()
            contstr = lineaux[0:4]
            coord[ii,:] = [float(contstr[0]),float(contstr[1]),float(contstr[2]),float(contstr[3])]
        
        file_object.readline()
        file_object.readline()
        NELM = int(file_object.readline())
        conec_elm = np.zeros((NELM,9))
        prop_elm = np.zeros((NELM,3))
        contelm=0
        for kk in range(0,NELM):
            line = file_object.readline()
            lineaux = line.split()
            type_elm = int(lineaux[1])
            if type_elm == 2:
                contelm+=1
                prop_elm[contelm-1,:] = [type_elm,float(lineaux[3]),float(lineaux[3])]
                conec_elm[contelm-1,:] = [contelm,float(lineaux[5]),float(lineaux[6]),float(lineaux[7]),0,0,0,0,0]


    conec_elm = conec_elm[:contelm,:]
    prop_elm = prop_elm[:contelm,:]
    inci = np.concatenate((conec_elm[:,0][:, np.newaxis],prop_elm,conec_elm[:,1:]),axis=1)      
    return coord, inci, type_elm


# type_elm == 3 -> 4-node quadrilateral
def import_mesh_quad(file_imp):
    with open(file_imp,'r') as file_object:
        file_object.readline()
        NNOD = int(file_object.readline())
        coord = np.zeros((NNOD,4))
        for ii in range(0,NNOD):
            line = file_object.readline()
            lineaux = line.split()
            contstr = lineaux[0:4]
            coord[ii,:] = [float(contstr[0]),float(contstr[1]),float(contstr[2]),float(contstr[3])]
        
        file_object.readline()
        file_object.readline()
        NELM = int(file_object.readline())
        conec_elm = np.zeros((NELM,9))
        prop_elm = np.zeros((NELM,3))
        contelm=0
        for kk in range(0,NELM):
            line = file_object.readline()
            lineaux = line.split()
            type_elm = int(lineaux[1])
            if type_elm == 3:
                contelm+=1
                prop_elm[contelm-1,:] = [type_elm,float(lineaux[3]),float(lineaux[3])]
                conec_elm[contelm-1,:] = [contelm,float(lineaux[5]),float(lineaux[6]),float(lineaux[7]),float(lineaux[8]),0,0,0,0]

    conec_elm = conec_elm[:contelm,:]
    prop_elm = prop_elm[:contelm,:]
    inci = np.concatenate((conec_elm[:,0][:, np.newaxis],prop_elm,conec_elm[:,1:]),axis=1)
    return coord, inci, type_elm

# type_elm == 5 -> 8-node hexaedro
def import_mesh_hexa(file_imp):
    with open(file_imp,'r') as file_object:
        file_object.readline()
        NNOD = int(file_object.readline())
        coord = np.zeros((NNOD,4))
        for ii in range(0,NNOD):
            line = file_object.readline()
            lineaux = line.split()
            contstr = lineaux[0:4]
            coord[ii,:] = [float(contstr[0]),float(contstr[1]),float(contstr[2]),float(contstr[3])]
        
        file_object.readline()
        file_object.readline()
        NELM = int(file_object.readline())
        conec_elm = np.zeros((NELM,9))
        prop_elm = np.zeros((NELM,3))
        contelm=0
        for kk in range(0,NELM):
            line = file_object.readline()
            lineaux = line.split()
            type_elm = int(lineaux[1])
            if type_elm == 5:
                contelm+=1
                prop_elm[contelm-1,:] = [type_elm,float(lineaux[3]),float(lineaux[3])]
                conec_elm[contelm-1,:] = [contelm,float(lineaux[5]),float(lineaux[6]),float(lineaux[7]),float(lineaux[8]),float(lineaux[9]),float(lineaux[10]),float(lineaux[11]),float(lineaux[12])]

    conec_elm = conec_elm[:contelm,:]
    prop_elm = prop_elm[:contelm,:]
    inci = np.concatenate((conec_elm[:,0][:, np.newaxis],prop_elm,conec_elm[:,1:]),axis=1)
    return coord, inci, type_elm




#%% 
def meshgen_data(usrlog,name_out_mesh):
    
    if usrlog.mesh_typ == "line2node":
        from myfempy.solve.read_mesh import import_mesh_line
        coord, inci, type_elm = import_mesh_line(name_out_mesh)
    elif usrlog.mesh_typ == "tria3node":
        from myfempy.solve.read_mesh import import_mesh_tria
        coord, inci, type_elm  = import_mesh_tria(name_out_mesh)
    elif usrlog.mesh_typ == "quad4node": 
        from myfempy.solve.read_mesh import import_mesh_quad
        coord, inci, type_elm = import_mesh_quad(name_out_mesh)
    elif usrlog.mesh_typ == "hexa8node":
        from myfempy.solve.read_mesh import import_mesh_hexa
        coord, inci, type_elm = import_mesh_hexa(name_out_mesh)
    else:
        print("input erro: mesh_typ don't defined")
    
    if usrlog.mod_typ == "beam":
        if usrlog.mod_opt == "spring20":
            IdDOF = 'spring20'
            ngdl = 2        # num gdl por no
            datamesh = [ngdl,len(coord),len(inci),ngdl*len(coord),IdDOF]
        elif usrlog.mod_opt == "truss22":
            IdDOF = 'truss22'
            ngdl = 2        # num gdl por no
            datamesh = [ngdl,len(coord),len(inci),ngdl*len(coord),IdDOF]
        elif usrlog.mod_opt == "beam21":
            IdDOF = "beam21"
            ngdl = 2        # num gdl por no
            datamesh = [ngdl,len(coord),len(inci),ngdl*len(coord),IdDOF]
        elif usrlog.mod_opt == "frame22":
            IdDOF = "frame22"
            ngdl = 3        # num gdl por no
            datamesh = [ngdl,len(coord),len(inci),ngdl*len(coord),IdDOF]
        elif usrlog.mod_opt == "frame23":
            IdDOF = "frame23"
            ngdl = 6        # num gdl por no
            datamesh = [ngdl,len(coord),len(inci),ngdl*len(coord),IdDOF]
        else:
            print("input erro: fin_elm_stf don't defined")
    elif usrlog.mod_typ == "plate":
        if usrlog.mod_opt == "plane32":
            IdDOF = "plane32"
            ngdl = 2        # num gdl por no
            datamesh = [ngdl,len(coord),len(inci),ngdl*len(coord),IdDOF]
        elif usrlog.mod_opt == "plane42":
            IdDOF = "plane42"
            ngdl = 2        # num gdl por no
            datamesh = [ngdl,len(coord),len(inci),ngdl*len(coord),IdDOF]
        else:
            print("input erro: fin_elm_stf don't defined")
    elif usrlog.mod_typ == "solid":
        if usrlog.mod_opt == "solid83":
            IdDOF = "solid83"
            ngdl = 3        # num gdl por no
            datamesh = [ngdl,len(coord),len(inci),ngdl*len(coord),IdDOF]
    else:
        print("input erro: fin_elm_stf don't defined")
    
    return datamesh, coord, inci, type_elm