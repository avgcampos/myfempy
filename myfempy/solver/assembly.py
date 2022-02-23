# -*- coding: utf-8 -*-
"""
Created on Sat May 16 12:05:12 2020
@author: ANTONIO VINICIUS GARCIA CAMPOS
@version: v20
_______________________________________________________________________________
 ~~~~~~~~~~        MATRIZ DE RIGIDEZ DOS ELEMENTOS FINITOS           ~~~~~~~~~~

ESTE MODULO DETERMINA A MATRIZ DE RIGIDEZ DO ELEMENTO FINITO E FAZ A MONTAGEM 
DA RIGIDEZ GLOBAL

> FUNCIONALIDADES
--- ENTRADAS: datamesh,inci,coord,tabmat,tabgeo
--- SAIDA: mtKG   
_______________________________________________________________________________
"""
import sys
import numpy as np

def stifness_matrix(usrlog,datamesh,coord,inci,tabmat,tabgeo,forcenodeaply):
    # TIPO DE ELEMENTO FINITO DA ANALISE   
    
    if usrlog.mod_typ == "beam":
        if usrlog.mod_opt == "spring20":
            from myfempy.lib.beam_fem_lib import spring20_stif
            typeMechMat = [usrlog.mat_opt,usrlog.mat_def]
            KG = spring20_stif(datamesh,inci,coord,typeMechMat,tabmat,tabgeo)
        elif usrlog.mod_opt == "truss22":
            from myfempy.lib.beam_fem_lib import truss22_stif
            typeMechMat = [usrlog.mat_opt,usrlog.mat_def]
            KG = truss22_stif(datamesh,inci,coord,typeMechMat,tabmat,tabgeo)
        elif usrlog.mod_opt == "beam21":
            from myfempy.lib.beam_fem_lib import beam21_stif
            typeMechMat = [usrlog.mat_opt,usrlog.mat_def]
            KG = beam21_stif(datamesh,inci,coord,typeMechMat,tabmat,tabgeo)
        elif usrlog.mod_opt == "frame22":
            from myfempy.lib.beam_fem_lib import frame22_stif
            typeMechMat = [usrlog.mat_opt,usrlog.mat_def]
            KG = frame22_stif(datamesh,inci,coord,typeMechMat,tabmat,tabgeo)
        elif usrlog.mod_opt == "frame23":
            from myfempy.lib.beam_fem_lib import frame23_stif
            typeMechMat = [usrlog.mat_opt,usrlog.mat_def]
            KG = frame23_stif(datamesh,inci,coord,typeMechMat,tabmat,tabgeo)
        else:
            print("input erro: fin_elm_stf don't defined")
    elif usrlog.mod_typ == "plate":
        if usrlog.mod_opt == "plane32":
            from myfempy.lib.plate_fem_lib import plane32_stif
            typeMechMat = [usrlog.mat_opt,usrlog.mat_def]
            KG = plane32_stif(datamesh,inci,coord,typeMechMat,tabmat,tabgeo)
        elif usrlog.mod_opt == "plane42":
            from myfempy.lib.plate_fem_lib import plane42_stif
            typeMechMat = [usrlog.mat_opt,usrlog.mat_def]
            npp = 4 #num. points Gauss
            KG = plane42_stif(datamesh,inci,coord,typeMechMat,tabmat,tabgeo,npp)
        else:
            print("input erro: fin_elm_stf don't defined")
    elif usrlog.mod_typ == "solid":
        if usrlog.mod_opt == "solid83":
            from myfempy.lib.solid_fem_lib import solid83_stif
            typeMechMat = [usrlog.mat_opt,usrlog.mat_def]
            npp = 8 #num. points Gauss
            KG = solid83_stif(datamesh,inci,coord,typeMechMat,tabmat,tabgeo,npp)
    else:
        print("input erro: fin_elm_stf don't defined")
        
    if np.any(forcenodeaply[np.where(forcenodeaply[:,1]>=6),:][0]) == True:
        springadd = forcenodeaply[np.where(forcenodeaply[:,1]==8.0),:][0]

        for ii in range(len(springadd)):
            noi = int(springadd[ii,0])
            loc = np.array([datamesh[0]*noi-2,datamesh[0]*noi-1])
            k0d =  springadd[ii,2]*np.array([[1,-1],[-1,1]])
            KG[np.ix_(loc,loc)] += k0d
        
    return KG


def mass_matrix(usrlog,datamesh,coord,inci,tabmat,tabgeo,forcenodeaply):
    # TIPO DE ELEMENTO FINITO DA ANALISE
    if usrlog.mod_typ == "beam":
        if usrlog.mod_opt == "spring20":
            from myfempy.lib.beam_fem_lib import spring20_mass
            typeMechMat = [usrlog.mat_opt,usrlog.mat_def]
            MG = spring20_mass(datamesh,inci,coord,typeMechMat,tabmat,tabgeo)
        if usrlog.mod_opt == "beam21":
            from myfempy.lib.beam_fem_lib import beam21_mass
            typeMechMat = [usrlog.mat_opt,usrlog.mat_def]
            MG = beam21_mass(datamesh,inci,coord,typeMechMat,tabmat,tabgeo)
        elif usrlog.mod_opt == "frame22":
            from myfempy.lib.beam_fem_lib import frame22_mass
            typeMechMat = [usrlog.mat_opt,usrlog.mat_def]
            MG = frame22_mass(datamesh,inci,coord,typeMechMat,tabmat,tabgeo)
        elif usrlog.mod_opt == "frame23":
            from myfempy.lib.beam_fem_lib import frame23_mass
            typeMechMat = [usrlog.mat_opt,usrlog.mat_def]
            MG = frame23_mass(datamesh,inci,coord,typeMechMat,tabmat,tabgeo)
        else:
            print("input erro: fin_elm_stf don't defined")
    elif usrlog.mod_typ == "plate":
        if usrlog.mod_opt == "plane32":
            from myfempy.lib.plate_fem_lib import plane32_mass
            typeMechMat = [usrlog.mat_opt,usrlog.mat_def]
            MG = plane32_mass(datamesh,inci,coord,typeMechMat,tabmat,tabgeo)
        elif usrlog.mod_opt == "plane42":
            from myfempy.lib.plate_fem_lib import plane42_mass
            typeMechMat = [usrlog.mat_opt,usrlog.mat_def]
            npp = 4 #num. points Gauss
            MG = plane42_mass(datamesh,inci,coord,typeMechMat,tabmat,tabgeo,npp)
        else:
            print("input erro: fin_elm_stf don't defined")
    elif usrlog.mod_typ == "solid":
        if usrlog.mod_opt == "solid83":
            from myfempy.lib.solid_fem_lib import solid83_mass
            typeMechMat = [usrlog.mat_opt,usrlog.mat_def]
            npp = 8 #num. points Gauss
            MG = solid83_mass(datamesh,inci,coord,typeMechMat,tabmat,tabgeo,npp)
    else:
        print("input erro: fin_elm_stf don't defined")
    
    if np.any(forcenodeaply[np.where(forcenodeaply[:,1]>=6),:][0]) == True:
        massadd = forcenodeaply[np.where(forcenodeaply[:,1]==7.0),:][0]
        for ii in range(len(massadd)):
            noi = int(massadd[ii,0])
            loc = np.array([datamesh[0]*noi-2,datamesh[0]*noi-1])
            m0d =  massadd[ii,2]*np.eye(2)
            MG[np.ix_(loc,loc)] += m0d

    
    return MG

#%%----------------------------------------------------------------------------
# Calculo da carga distribuida equivalente nodal
def force_edge(datamesh,inci,coord,force_value,force_dirc,node_list_fc,dir_fc,tabgeo):
    
    if (dir_fc == 'y_x')or(dir_fc == 'z_x'):
        coord_fc = 1
    elif (dir_fc == 'x_y')or(dir_fc == 'z_y'):
        coord_fc = 2
    elif (dir_fc == 'x_z')or(dir_fc == 'y_z'): 
        coord_fc = 3            
            
    if datamesh[4] == 'plane32':
        elmlist = np.zeros((1))
        for ii in range(len(node_list_fc)):
            elm2list = inci[(np.asarray(np.where(inci[:,4:]==node_list_fc[ii])))[0][:],0]
            elmlist = np.append(elmlist,[elm2list[0]],axis=0)
        
        elmlist = np.unique(elmlist)
        elmlist = elmlist[2::][::]
        
        force_value_vector = np.zeros((datamesh[0]*len(node_list_fc),1))
        for i in range(len(elmlist)):
            noi = int(inci[int(elmlist[i]-1),4])
            noj = int(inci[int(elmlist[i]-1),5])
            nok = int(inci[int(elmlist[i]-1),6])
            
            if np.any([noi==node_list_fc[:]]):
                no1=noi
                if np.any([noj==node_list_fc[:]]):
                    no2=noj
                elif np.any([nok==node_list_fc[:]]):
                    no2=nok
            elif np.any([noj==node_list_fc[:]]):
                no1=noj
                if np.any([nok==node_list_fc[:]]):
                    no2=nok
            
            L = abs(coord[no1-1,coord_fc] - coord[no2-1,coord_fc])
            
            no1dof = np.where(no1==node_list_fc)[0][0]
            no2dof = np.where(no2==node_list_fc)[0][0]
            loc = np.array([datamesh[0]*no1dof,datamesh[0]*no1dof+1,datamesh[0]*no2dof,datamesh[0]*no2dof+1])
            tck = tabgeo[int(inci[int(elmlist[i]-1),3]-1),4]
            if force_dirc == 'fx':
                force_value_vector[loc,0] += np.array([force_value*tck*L/2,0,force_value*tck*L/2,0])
                fc_type_dof = np.array([1,0])
            elif force_dirc == 'fy':
                force_value_vector[loc,0] += np.array([0,force_value*tck*L/2,0,force_value*tck*L/2])
                fc_type_dof = np.array([0,2])
            

    elif datamesh[4] == 'plane42':
        elmlist = np.zeros((1))
        for ii in range(len(node_list_fc)):
            elm2list = inci[(np.asarray(np.where(inci[:,4:]==node_list_fc[ii])))[0][:],0]
            elmlist = np.append(elmlist,[elm2list[0]],axis=0)
        
        elmlist = np.unique(elmlist)
        elmlist = elmlist[1::][::]
        
        force_value_vector = np.zeros((datamesh[0]*len(node_list_fc),1))
        for i in range(len(elmlist)):
            noi = int(inci[int(elmlist[i]-1),4])
            noj = int(inci[int(elmlist[i]-1),5])
            nok = int(inci[int(elmlist[i]-1),6])
            nol = int(inci[int(elmlist[i]-1),7])
            
            if np.any([noi==node_list_fc[:]]):
                no1=noi
                if np.any([noj==node_list_fc[:]]):
                    no2=noj
                elif np.any([nol==node_list_fc[:]]):
                    no2=nol
            elif np.any([noj==node_list_fc[:]]):
                no1=noj
                if np.any([nok==node_list_fc[:]]):
                    no2=nok
            elif np.any([nok==node_list_fc[:]]):
                no1=nok
                if np.any([nol==node_list_fc[:]]):
                    no2=nol
            
            L = abs(coord[no1-1,coord_fc] - coord[no2-1,coord_fc])
            
            no1dof = np.where(no1==node_list_fc)[0][0]
            no2dof = np.where(no2==node_list_fc)[0][0]
            loc = np.array([datamesh[0]*no1dof,datamesh[0]*no1dof+1,datamesh[0]*no2dof,datamesh[0]*no2dof+1])
            tck = tabgeo[int(inci[int(elmlist[i]-1),3]-1),4]
            if force_dirc == 'fx':
                force_value_vector[loc,0] += np.array([force_value*tck*L/2,0,force_value*tck*L/2,0])
                fc_type_dof = np.array([1,0])
            elif force_dirc == 'fy':
                force_value_vector[loc,0] += np.array([0,force_value*tck*L/2,0,force_value*tck*L/2])
                fc_type_dof = np.array([0,2])
            
            
    return force_value_vector,fc_type_dof

# Calculo da carga distribuida equivalente nodal
def force_surf(datamesh,inci,coord,force_value,force_dirc,node_list_fc,dir_fc):
    if dir_fc == 'x':
        coord_fc = [3,2]
    elif dir_fc == 'y':
        coord_fc = [1,3]
    elif dir_fc == 'z':
        coord_fc = [1,2]

    elmlist = np.zeros((1))
    for ii in range(len(node_list_fc)):
        elm2list = inci[(np.asarray(np.where(inci[:,4:]==node_list_fc[ii])))[0][:],0]
        elmlist = np.append(elmlist,[elm2list[0]],axis=0)
        
    elmlist = np.unique(elmlist)
    elmlist = elmlist[1::][::]
    
    force_value_vector = np.zeros((datamesh[0]*len(node_list_fc),1))
    for i in range(len(elmlist)):
        noi = int(inci[int(elmlist[i]-1),4])
        noj = int(inci[int(elmlist[i]-1),5])
        nok = int(inci[int(elmlist[i]-1),6])
        nol = int(inci[int(elmlist[i]-1),7])
        nom = int(inci[int(elmlist[i]-1),8])
        non = int(inci[int(elmlist[i]-1),9])
        noo = int(inci[int(elmlist[i]-1),10])
        nop = int(inci[int(elmlist[i]-1),11])

        if np.any([noi==node_list_fc[:]]):
            no1=noi
            if np.any([noj==node_list_fc[:]]):
                no2=noj
                if np.any([nok==node_list_fc[:]]):
                    no3=nok
                    no4=nol
            elif np.any([nom==node_list_fc[:]]):
                no2=nom
                if np.any([nop==node_list_fc[:]]):
                    no3=nop
                    no4=nol
                elif np.any([non==node_list_fc[:]]):
                    no3=non
                    no4=noj
        elif np.any([nom==node_list_fc[:]]):
            no1=nom
            if np.any([non==node_list_fc[:]]):
                no2=non
                if np.any([noo==node_list_fc[:]]):
                    no3=non
                    no4=nop
        elif np.any([non==node_list_fc[:]]):
            no1=non
            if np.any([noo==node_list_fc[:]]):
                no2=noo
                no3=nok
                no4=noj
        elif np.any([noo==node_list_fc[:]]):
            no1=non
            if np.any([nok==node_list_fc[:]]):
                no2=nok
                no3=nol
                no4=nop
        
        
        
        coordX_max = max(coord[no1-1,coord_fc[0]],coord[no2-1,coord_fc[0]],coord[no3-1,coord_fc[0]],coord[no4-1,coord_fc[0]])
        coordX_min = min(coord[no1-1,coord_fc[0]],coord[no2-1,coord_fc[0]],coord[no3-1,coord_fc[0]],coord[no4-1,coord_fc[0]])
        coordY_max = max(coord[no1-1,coord_fc[1]],coord[no2-1,coord_fc[1]],coord[no3-1,coord_fc[1]],coord[no4-1,coord_fc[1]])
        coordY_min = min(coord[no1-1,coord_fc[1]],coord[no2-1,coord_fc[1]],coord[no3-1,coord_fc[1]],coord[no4-1,coord_fc[1]])
        
        coordX = abs(coordX_max-coordX_min)
        coordY = abs(coordY_max-coordY_min)
        
        A = coordX*coordY
        
        no1dof = np.where(no1==node_list_fc)[0][0]
        no2dof = np.where(no2==node_list_fc)[0][0]
        no3dof = np.where(no3==node_list_fc)[0][0]
        no4dof = np.where(no4==node_list_fc)[0][0]
        loc = np.array([datamesh[0]*no1dof,datamesh[0]*no1dof+1,datamesh[0]*no1dof+2,\
                        datamesh[0]*no2dof,datamesh[0]*no2dof+1,datamesh[0]*no2dof+2,\
                        datamesh[0]*no3dof,datamesh[0]*no3dof+1,datamesh[0]*no3dof+2,\
                        datamesh[0]*no4dof,datamesh[0]*no4dof+1,datamesh[0]*no4dof+2,])
            
        if force_dirc == 'fx':
            force_value_vector[loc,0] += np.array([force_value*A/4,0.0,0.0,force_value*A/4,0.0,0.0,force_value*A/4,0.0,0.0,force_value*A/4,0.0,0.0])
            fc_type_dof = np.array([1,0,0])
        elif force_dirc == 'fy':
            force_value_vector[loc,0] += np.array([0.0,force_value*A/4,0.0,0.0,force_value*A/4,0.0,0.0,force_value*A/4,0.0,0.0,force_value*A/4,0.0])
            fc_type_dof = np.array([0,2,0])
        elif force_dirc == 'fz':
            force_value_vector[loc,0] += np.array([0.0,0.0,force_value*A/4,0.0,0.0,force_value*A/4,0.0,0.0,force_value*A/4,0.0,0.0,force_value*A/4])
            fc_type_dof = np.array([0,0,3])
            
    return force_value_vector,fc_type_dof

def force_line(datamesh,inci,coord,force_value,force_dirc,dir_fc,node_list_fc,line_fc):
              
    inci_fc = inci[np.where(inci[:,3]==line_fc),:][0]
    
    if dir_fc == 'x':
        coord_fc = 1
    elif dir_fc == 'y':
        coord_fc = 2
    elif dir_fc == 'z':
        coord_fc = 3
    
    if datamesh[4] == 'beam21':
        elmlist = np.zeros((1))
        for ii in range(len(node_list_fc)):
            elm2list = inci_fc[(np.asarray(np.where(inci_fc[:,4:]==node_list_fc[ii])))[0][:],0]
            if elm2list.size != 0: 
                elmlist = np.append(elmlist,[elm2list[0]],axis=0)                
        
        elmlist = np.unique(elmlist)
        elmlist = elmlist[1::][::]
        
        force_value_vector = np.zeros((datamesh[0]*len(node_list_fc),1))
        for i in range(len(elmlist)):
            no1 = int(inci[int(elmlist[i]-1),4])
            no2 = int(inci[int(elmlist[i]-1),5])
            
            no1dof = np.asarray(np.where(no1==node_list_fc))
            no2dof = np.asarray(np.where(no2==node_list_fc))
            if no1dof.size != 0:
                if no2dof.size != 0:
                    L = abs(coord[no1-1,coord_fc] - coord[no2-1,coord_fc])
        
                    loc = np.array([datamesh[0]*no1dof[0][0],datamesh[0]*no1dof[0][0]+1,datamesh[0]*no2dof[0][0],datamesh[0]*no2dof[0][0]+1])
                    force_value_vector[loc,0] += np.array([force_value*L/2,(force_value*L**2)/12,force_value*L/2,(-1*force_value*L**2)/12])
                    fc_type_dof = np.array([2,6])
            
    elif datamesh[4] == 'frame22':
        elmlist = np.zeros((1))
        for ii in range(len(node_list_fc)):
            elm2list = inci_fc[(np.asarray(np.where(inci_fc[:,4:]==node_list_fc[ii])))[0][:],0]
            if elm2list.size != 0: 
                elmlist = np.append(elmlist,[elm2list[0]],axis=0)
                
        elmlist = np.unique(elmlist)
        elmlist = elmlist[1::][::]
        
        force_value_vector = np.zeros((datamesh[0]*len(node_list_fc),1))
        for i in range(len(elmlist)):
            no1 = int(inci[int(elmlist[i]-1),4])
            no2 = int(inci[int(elmlist[i]-1),5])
            
            no1dof = np.asarray(np.where(no1==node_list_fc))
            no2dof = np.asarray(np.where(no2==node_list_fc))
            if no1dof.size != 0:
                if no2dof.size != 0:
                    L = abs(coord[no1-1,coord_fc] - coord[no2-1,coord_fc])
        
                    loc = np.array([datamesh[0]*no1dof[0][0],datamesh[0]*no1dof[0][0]+1,datamesh[0]*no1dof[0][0]+2,datamesh[0]*no2dof[0][0],datamesh[0]*no2dof[0][0]+1,datamesh[0]*no2dof[0][0]+2])
                    if force_dirc == 'fx':
                        force_value_vector[loc,0] += np.array([(force_value*L**2)/2,0.0,0.0,(force_value*L**2)/2,0.0,0.0])
                        fc_type_dof = np.array([1,0,0])
                    elif force_dirc == 'fy':
                        force_value_vector[loc,0] += np.array([0,force_value*L/2,(force_value*L**2)/12,0,force_value*L/2,(-1*force_value*L**2)/12])
                        fc_type_dof = np.array([0,2,6])
            
        
    elif datamesh[4] == 'frame23':
        elmlist = np.zeros((1))
        for ii in range(len(node_list_fc)):
            elm2list = inci_fc[(np.asarray(np.where(inci_fc[:,4:]==node_list_fc[ii])))[0][:],0]
            if elm2list.size != 0: 
                elmlist = np.append(elmlist,[elm2list[0]],axis=0)
        
        elmlist = np.unique(elmlist)
        elmlist = elmlist[1::][::]
        
        force_value_vector = np.zeros((datamesh[0]*len(node_list_fc),1))
        for i in range(len(elmlist)):
            no1 = int(inci[int(elmlist[i]-1),4])
            no2 = int(inci[int(elmlist[i]-1),5])
            
            no1dof = np.asarray(np.where(no1==node_list_fc))
            no2dof = np.asarray(np.where(no2==node_list_fc))
            if no1dof.size != 0:
                if no2dof.size != 0:
                    L = abs(coord[no1-1,coord_fc] - coord[no2-1,coord_fc])
        
                    loc = np.array([datamesh[0]*no1dof[0][0],datamesh[0]*no1dof[0][0]+1,datamesh[0]*no1dof[0][0]+2,datamesh[0]*no1dof[0][0]+3,datamesh[0]*no1dof[0][0]+4,datamesh[0]*no1dof[0][0]+5,\
                                   datamesh[0]*no2dof[0][0],datamesh[0]*no2dof[0][0]+1,datamesh[0]*no2dof[0][0]+2,datamesh[0]*no2dof[0][0]+3,datamesh[0]*no2dof[0][0]+4,datamesh[0]*no2dof[0][0]+5])
                    if force_dirc == 'fx':
                        force_value_vector[loc,0] += np.array([(force_value*L**2)/2, 0, 0, 0, 0, 0, (force_value*L**2)/2, 0, 0, 0, 0, 0])
                        fc_type_dof = np.array([1,0,0,0,0,0])
                    elif force_dirc == 'fy':
                        force_value_vector[loc,0] += np.array([0, force_value*L/2, 0, (force_value*L**2)/12, 0, 0, 0, force_value*L/2, 0, (-1*force_value*L**2)/12, 0, 0])
                        fc_type_dof = np.array([0,2,0,4,0,0])
                    elif force_dirc == 'fz':
                        force_value_vector[loc,0] += np.array([0, 0, force_value*L/2, 0, (force_value*L**2)/12, 0, 0, 0, force_value*L/2, 0, (-1*force_value*L**2)/120, 0])
                        fc_type_dof = np.array([0,0,3,0,5,0])
  
    
    return force_value_vector,fc_type_dof

#------------------------------------------------------------------------------
#%% busca de nós específicos na malha X
def search_edgeX(edge_coordX,coord,erro):
    dif = abs(edge_coordX*np.ones_like(coord[:,1]) - coord[:,1])
    node_posi=(np.asarray(np.where(dif<erro)))[0][:]
    node=coord[node_posi,0].astype(int)
    return node

# busca de nós específicos na malha Y
def search_edgeY(edge_coordY,coord,erro):
    dif = abs(edge_coordY*np.ones_like(coord[:,2]) - coord[:,2])
    node_posi=(np.asarray(np.where(dif<erro)))[0][:]
    node=coord[node_posi,0].astype(int)
    return node

# busca de nós específicos na malha Z
def search_edgeZ(edge_coordZ,coord,erro):
    dif = abs(edge_coordZ*np.ones_like(coord[:,3]) - coord[:,3])
    node_posi=(np.asarray(np.where(dif<erro)))[0][:]
    node=coord[node_posi,0].astype(int)
    return node

# busca de nós específicos na malha XYZ
def search_surfXY(orthg_coordZ,coord,erro):
    dif = abs(orthg_coordZ*np.ones_like(coord[:,3]) - coord[:,3])
    node_posi=(np.asarray(np.where(dif<erro)))[0][:]
    node=coord[node_posi,0].astype(int)
    return node

# busca de nós específicos na malha XYZ
def search_surfYZ(orthg_coordX,coord,erro):
    dif = abs(orthg_coordX*np.ones_like(coord[:,1]) - coord[:,1])
    node_posi=(np.asarray(np.where(dif<erro)))[0][:]
    node=coord[node_posi,0].astype(int)
    return node

# busca de nós específicos na malha XYZ
def search_surfZX(orthg_coordY,coord,erro):
    dif = abs(orthg_coordY*np.ones_like(coord[:,2]) - coord[:,2])
    node_posi=(np.asarray(np.where(dif<erro)))[0][:]
    node=coord[node_posi,0].astype(int)
    return node

# # busca de nós específicos na malha XY
# def search_nodeXY(node_coordX,node_coordY,coord,erro):
#     difx = abs(node_coordX*np.ones_like(coord[:,1]) - coord[:,1])
#     dify = abs(node_coordY*np.ones_like(coord[:,2]) - coord[:,2])
#     node_posi=(np.asarray(np.where((difx<erro)&(dify<erro))))[0][:]
#     node=coord[node_posi,0].astype(int)
#     return node

# busca de nós específicos na malha XYZ
def search_nodeXYZ(node_coordX,node_coordY,node_coordZ,coord,erro):
    difx = abs(node_coordX*np.ones_like(coord[:,1]) - coord[:,1])
    dify = abs(node_coordY*np.ones_like(coord[:,2]) - coord[:,2])
    difz = abs(node_coordZ*np.ones_like(coord[:,3]) - coord[:,3])
    node_posi=(np.asarray(np.where((difx<erro)&(dify<erro)&(difz<erro))))[0][:]
    node=coord[node_posi,0].astype(int)
    return node


def model_physics(usrlog,datamesh,coord,inci,forcelistfull,boncdlist,tabgeo):
    
    forcenodeaply = np.zeros((1,4))
    for k in range(int(usrlog.solver_end)):

        forcelist = forcelistfull[np.where(forcelistfull[:,8]==str(k+1))[0],:]
        
        numfrclist = len(forcelist)
        node_list_fc = np.array([])
        for i in range(numfrclist):
            if eval('usrlog.force_opt_'+str(i)) == "linex":
                coord_0 = float(forcelist[i,6])
                coord_1 = float(forcelist[i,7])
                node_list_fc = coord[np.where((coord[:,1]>=coord_0)&(coord[:,1]<=coord_1)),0][0]
                dir_fc = 'x'
            
            elif eval('usrlog.force_opt_'+str(i)) == "liney":
                coord_0 = float(forcelist[i,6])
                coord_1 = float(forcelist[i,7])
                node_list_fc = coord[np.where((coord[:,2]>=coord_0)&(coord[:,2]<=coord_1)),0][0]
                dir_fc = 'y'
            
            elif eval('usrlog.force_opt_'+str(i)) == "linez":
                coord_0 = float(forcelist[i,6])
                coord_1 = float(forcelist[i,7])
                node_list_fc = coord[np.where((coord[:,3]>=coord_0)&(coord[:,3]<=coord_1)),0][0]
                dir_fc = 'z'
           
            elif eval('usrlog.force_opt_'+str(i)) == "edgex":
                # sys.path.append('../lib')
                # from domain_physics import search_edgeX
                edge_coordX = float(forcelist[i,5])
                            
                if float(forcelist[i,6]) == 999:
                    dir_fc = 'x_y'
                    coord_fc=(coord[np.where(coord[:,3]==float(forcelist[i,7])),:])[0]
                elif float(forcelist[i,7]) == 999:
                    dir_fc = 'x_z'
                    coord_fc=(coord[np.where(coord[:,2]==float(forcelist[i,6])),:])[0]
                    
                node_list_fc = search_edgeX(edge_coordX,coord_fc,2E-3)
            
            elif eval('usrlog.force_opt_'+str(i)) == "edgey":
                # sys.path.append('../lib')
                # from domain_physics import search_edgeY
                edge_coordY = float(forcelist[i,6])
                
                if float(forcelist[i,5]) == 999:
                    dir_fc = 'y_x'
                    coord_fc=(coord[np.where(coord[:,3]==float(forcelist[i,7])),:])[0]
                elif float(forcelist[i,7]) == 999:
                    dir_fc = 'y_z'
                    coord_fc=(coord[np.where(coord[:,1]==float(forcelist[i,5])),:])[0]
                
                node_list_fc = search_edgeY(edge_coordY,coord_fc,2E-3)
                
            elif eval('usrlog.force_opt_'+str(i)) == "edgez":
                # sys.path.append('../lib')
                # from domain_physics import search_edgeZ
                edge_coordZ = float(forcelist[i,7])
                
                if float(forcelist[i,5]) == 999:
                    dir_fc = 'z_x'
                    coord_fc=(coord[np.where(coord[:,2]==float(forcelist[i,6])),:])[0]
                elif float(forcelist[i,6]) == 999:
                    dir_fc = 'z_x'
                    coord_fc=(coord[np.where(coord[:,1]==float(forcelist[i,5])),:])[0]
                
                node_list_fc = search_edgeZ(edge_coordZ,coord_fc,2E-3)      
        
                
            elif eval('usrlog.force_opt_'+str(i)) == "surfxy":
                # sys.path.append('../lib')
                # from domain_physics import search_surfXY
                orthg_coordZ = float(forcelist[i,7])
                node_list_fc = search_surfXY(orthg_coordZ,coord,2E-3)
                dir_fc = 'z'
                
            elif eval('usrlog.force_opt_'+str(i)) == "surfyz":
                # sys.path.append('../lib')
                # from domain_physics import search_surfYZ
                orthg_coordX = float(forcelist[i,5])
                node_list_fc = search_surfYZ(orthg_coordX,coord,2E-3)
                dir_fc = 'x'
            
            elif eval('usrlog.force_opt_'+str(i)) == "surfzx":
                # sys.path.append('../lib')
                # from domain_physics import search_surfZX
                orthg_coordY = float(forcelist[i,6])
                node_list_fc = search_surfZX(orthg_coordY,coord,2E-3)
                dir_fc = 'y'
            
            elif eval('usrlog.force_opt_'+str(i)) == "point":
                # sys.path.append('../lib')
                # from domain_physics import search_nodeXYZ
                node_coordX = float(forcelist[i,5])
                node_coordY = float(forcelist[i,6])
                node_coordZ = float(forcelist[i,7])
                node_list_fc = search_nodeXYZ(node_coordX,node_coordY,node_coordZ,coord,2E-3)
            else:
                print("input erro: force_opt don't defined")
        
            
            if eval('usrlog.force_typ_'+str(i)) == "forcenode":
                force_value_vector = np.ones_like(node_list_fc)*float(forcelist[i,3])
                
                if eval('usrlog.force_opt_dir'+str(i)) == 'fx':
                    fc_type_dof = 1*np.ones_like(node_list_fc)
                elif eval('usrlog.force_opt_dir'+str(i)) == 'fy':
                    fc_type_dof = 2*np.ones_like(node_list_fc)
                elif eval('usrlog.force_opt_dir'+str(i)) == 'fz':
                    fc_type_dof = 3*np.ones_like(node_list_fc)
                elif eval('usrlog.force_opt_dir'+str(i)) == 'tx':
                    fc_type_dof = 4*np.ones_like(node_list_fc)
                elif eval('usrlog.force_opt_dir'+str(i)) == 'ty':
                    fc_type_dof = 5*np.ones_like(node_list_fc)
                elif eval('usrlog.force_opt_dir'+str(i)) == 'tz':
                    fc_type_dof = 6*np.ones_like(node_list_fc)
                elif eval('usrlog.force_opt_dir'+str(i)) == 'massadd':
                    fc_type_dof = 7*np.ones_like(node_list_fc)
                elif eval('usrlog.force_opt_dir'+str(i)) == 'spring2gd':
                    fc_type_dof = 8*np.ones_like(node_list_fc)
                elif eval('usrlog.force_opt_dir'+str(i)) == 'spring2gd':
                    fc_type_dof = 9*np.ones_like(node_list_fc)
                else:
                    print("input erro: force_opt_dir don't defined")
                    
            
                for j in range(len(node_list_fc)):
                    fcdef = np.array([[int(node_list_fc[j]),fc_type_dof[j],force_value_vector[j],int(eval('usrlog.force_step'+str(k)))]])
                    forcenodeaply = np.append(forcenodeaply,fcdef,axis=0)
        
            elif  eval('usrlog.force_typ_'+str(i)) == "forceline":
                   # sys.path.append('../lib')
                   # from domain_physics import force_line
                   
                   line_fc = int(float(forcelist[i,5]))
                   
                   force_value = float(forcelist[i,3])
                   force_dirc = forcelist[i,2]
                   
                   force_value_vector,fc_type_dof = force_line(datamesh,inci,coord,force_value,force_dirc,dir_fc,node_list_fc,line_fc)
        
                   fc_type_dof = np.tile(fc_type_dof, len(node_list_fc))
                    
                   node_list_fc = np.repeat(node_list_fc, datamesh[0], axis=0)
                                    
                   for j in range(len(node_list_fc)):
                       fcdef = np.array([[int(node_list_fc[j]),fc_type_dof[j],force_value_vector[j,0],int(eval('usrlog.force_step'+str(k)))]])
                       forcenodeaply = np.append(forcenodeaply,fcdef,axis=0)                
                       
            elif eval('usrlog.force_typ_'+str(i)) == "forceedge":
                # sys.path.append('../lib')
                # from domain_physics import force_edge
                
                force_value = float(forcelist[i,3])
                force_dirc = forcelist[i,2]
                
                force_value_vector,fc_type_dof = force_edge(datamesh,inci,coord,force_value,force_dirc,node_list_fc,dir_fc,tabgeo)
        
                fc_type_dof = np.tile(fc_type_dof, len(node_list_fc))
                
                node_list_fc = np.repeat(node_list_fc, datamesh[0], axis=0)
                                
                for j in range(len(node_list_fc)):
                    fcdef = np.array([[int(node_list_fc[j]),fc_type_dof[j],force_value_vector[j,0],int(eval('usrlog.force_step'+str(k)))]])
                    forcenodeaply = np.append(forcenodeaply,fcdef,axis=0)
                
                # sys.exit()
            elif eval('usrlog.force_typ_'+str(i)) == "forcesurf":
                # sys.path.append('../lib')
                # from domain_physics import force_surf
                
                force_value = float(forcelist[i,3])
                force_dirc = forcelist[i,2]        
                
                force_value_vector,fc_type_dof = force_surf(datamesh,inci,coord,force_value,force_dirc,node_list_fc,dir_fc)
                
                fc_type_dof = np.tile(fc_type_dof, len(node_list_fc))
                
                node_list_fc = np.repeat(node_list_fc, datamesh[0], axis=0)
                                
                for j in range(len(node_list_fc)):
                    fcdef = np.array([[int(node_list_fc[j]),fc_type_dof[j],force_value_vector[j,0],int(eval('usrlog.force_step'+str(k)))]])
                    forcenodeaply = np.append(forcenodeaply,fcdef,axis=0)
            
                
            
            # elif eval('usrlog.force_typ_'+str(i)) == "gravity":
            
            # elif eval('usrlog.force_typ_'+str(i)) == "pressure":
                
            else:
                print("input erro: force_typ don't defined")
                
    # CONDICCOES DE CONTORNO
    numbdclist = len(boncdlist)
    boncdnodeaply = np.zeros((1,2))
    node_list_bc = np.array([])
    for i in range(numbdclist):
        if eval('usrlog.bc_typ_'+str(i)) == "fixed":
            if eval('usrlog.bc_opt_'+str(i)) == "edgex":
                # sys.path.append('../lib')
                # from domain_physics import search_edgeX
                edge_coordX = float(boncdlist[i,4])
                
                if float(boncdlist[i,5]) == 999:
                    coord_bc=(coord[np.where(coord[:,3]==float(boncdlist[i,6])),:])[0]
                elif float(boncdlist[i,6]) == 999:
                    coord_bc=(coord[np.where(coord[:,2]==float(boncdlist[i,5])),:])[0]
                
                node_list_bc = search_edgeX(edge_coordX,coord_bc,2E-3)
            
            elif eval('usrlog.bc_opt_'+str(i)) == "edgey":
                # sys.path.append('../lib')
                # from domain_physics import search_edgeY
                edge_coordY = float(boncdlist[i,5])
                
                if float(boncdlist[i,4]) == 999:
                    coord_bc=(coord[np.where(coord[:,3]==float(boncdlist[i,6])),:])[0]
                elif float(boncdlist[i,6]) == 999:
                    coord_bc=(coord[np.where(coord[:,1]==float(boncdlist[i,4])),:])[0]
                
                node_list_bc = search_edgeY(edge_coordY,coord_bc,2E-3)
            
            elif eval('usrlog.bc_opt_'+str(i)) == "edgey":
                # sys.path.append('../lib')
                # from domain_physics import search_edgeZ
                edge_coordZ = float(boncdlist[i,6])
                
                if float(boncdlist[i,4]) == 999:
                    coord_bc=(coord[np.where(coord[:,2]==float(boncdlist[i,5])),:])[0]
                elif float(boncdlist[i,5]) == 999:
                    coord_bc=(coord[np.where(coord[:,1]==float(boncdlist[i,4])),:])[0]
                
                node_list_bc = search_edgeZ(edge_coordZ,coord_bc,2E-3)
            
            elif eval('usrlog.bc_opt_'+str(i)) == "surfxy":
                # sys.path.append('../lib')
                # from domain_physics import search_surfXY
                orthg_coordZ = float(boncdlist[i,6])
                node_list_bc = search_surfXY(orthg_coordZ,coord,2E-3)
            
            elif eval('usrlog.bc_opt_'+str(i)) == "surfyz":
                # sys.path.append('../lib')
                # from domain_physics import search_surfYZ
                orthg_coordX = float(boncdlist[i,4])
                node_list_bc = search_surfYZ(orthg_coordX,coord,2E-3)
            
            elif eval('usrlog.bc_opt_'+str(i)) == "surfzx":
                # sys.path.append('../lib')
                # from domain_physics import search_surfZX
                orthg_coordY = float(boncdlist[i,5])
                node_list_bc = search_surfZX(orthg_coordY,coord,2E-3)
            
            elif eval('usrlog.bc_opt_'+str(i)) == "point":
                # sys.path.append('../lib')
                # from domain_physics import search_nodeXYZ
                node_coordX = float(boncdlist[i,4])
                node_coordY = float(boncdlist[i,5]) 
                node_coordZ = float(boncdlist[i,6])
                node_list_bc = search_nodeXYZ(node_coordX,node_coordY,node_coordZ,coord,2E-3)
            else:
                print("input erro: bc_opt don't defined")
            
    
            nnodelist_bc = len(node_list_bc)
            for j in range(nnodelist_bc):            
                if eval('usrlog.bc_opt_dir'+str(i)) == 'ux':
                    bcdef = np.array([[1,int(node_list_bc[j])]])
                    boncdnodeaply = np.append(boncdnodeaply,bcdef,axis=0)
                elif eval('usrlog.bc_opt_dir'+str(i)) == 'uy':
                    bcdef = np.array([[2,int(node_list_bc[j])]])
                    boncdnodeaply = np.append(boncdnodeaply,bcdef,axis=0)
                elif eval('usrlog.bc_opt_dir'+str(i)) == 'uz':
                    bcdef = np.array([[3,int(node_list_bc[j])]])
                    boncdnodeaply = np.append(boncdnodeaply,bcdef,axis=0)
                elif eval('usrlog.bc_opt_dir'+str(i)) == 'rx':
                    bcdef = np.array([[4,int(node_list_bc[j])]])
                    boncdnodeaply = np.append(boncdnodeaply,bcdef,axis=0)
                elif eval('usrlog.bc_opt_dir'+str(i)) == 'ry':
                    bcdef = np.array([[5,int(node_list_bc[j])]])
                    boncdnodeaply = np.append(boncdnodeaply,bcdef,axis=0)
                elif eval('usrlog.bc_opt_dir'+str(i)) == 'rz':
                    bcdef = np.array([[6,int(node_list_bc[j])]])
                    boncdnodeaply = np.append(boncdnodeaply,bcdef,axis=0)
                elif eval('usrlog.bc_opt_dir'+str(i)) == 'full':
                    bcdef = np.array([[0,int(node_list_bc[j])]])
                    boncdnodeaply = np.append(boncdnodeaply,bcdef,axis=0)
                else:
                    print("input erro: bc_opt_dir don't defined")
        # elif force_typ == "forceedge":
        else:
            print("input erro: bc_typ don't defined")
    
    forcenodeaply = forcenodeaply[1::][::]
    boncdnodeaply = boncdnodeaply[1::][::]
    return forcenodeaply, boncdnodeaply