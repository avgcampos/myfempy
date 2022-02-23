# -*- coding: utf-8 -*-
"""
Created on Sat May 16 12:05:12 2020
@author: ANTONIO VINICIUS GARCIA CAMPOS
@version: v15
_______________________________________________________________________________
 ~~~~~~~~~~ MODULO DE SIMULACAO PELO METODO DOS ELEMENTOS FINITOS ~~~~~~~~~~

ESTE MODULO CONSTROE RESULTADOS PARA POS PROCESSAMENTO DA ANALISE 
- ENTRADAS: 
- SAIDA:     
===============================================================================

> ATUALIZACOES DA VERSAO:
_______________________________________________________________________________
"""
import numpy as np
import matplotlib.pyplot as plt
import sys
import re

def save_data(path_user,usr_analysi_name,file_dir,file_save,datamesh,inci,coord,type_elm,data_display,Udef,meshDefU,scale_mesh,data_result_save,data_result_view,data_title,title_display,typeData,scale_bar,title_scale_bar,screen_on,title_win,save_screen,screen_mode_shape,ModeNumb):
    from myfempy.solve.view_mesh import view_results
    from myfempy.solve.read_mesh import export_mesh, animation_mesh
    from myfempy.setup.myfempy_welcome import myfempy_ctrllVer
    myfempy_version = myfempy_ctrllVer()
    
    export_mesh(file_dir,type_elm,datamesh[1],datamesh[2],meshDefU,inci,data_result_view,data_title,typeData)
    view_results(file_dir,file_save,scale_bar,title_scale_bar,screen_on,save_screen,usr_analysi_name,title_win,myfempy_version,title_display,data_display,scala_view=0)
    
    export_mesh(file_dir,type_elm,datamesh[1],datamesh[2],meshDefU,inci,data_result_save,data_title,typeData)
    
    if screen_mode_shape == 'false':
        print(' ')
    else:
        time_out = 30
        animation_mesh(path_user,usr_analysi_name,type_elm,inci,coord,datamesh,Udef,data_result_save,data_title,typeData,time_out,scale_mesh,ModeNumb)

def beam_output(path_user,usr_analysi_name,usrlog,graphout,datamesh,inci,coord,tabgeo,CGcoord,tabmat,type_elm,KG,U,U2,U3,screen_on,save_screen,scale_mesh,datastep,hist_dof_list):
    
    screen_mode_shape = 'n'
     
    y_max = CGcoord[0]
    y_min = CGcoord[1]
    z_max = CGcoord[2]
    z_min = CGcoord[3]
    r_max = CGcoord[4]
    
    typeMechMat = [usrlog.mat_opt,usrlog.mat_def]
    
    max_displ = np.zeros((1))
    max_stress = np.zeros((1))
    max_strain = np.zeros((1))
    
    ModeNumb = int(datastep)
    if usrlog.mod_opt == "spring20":
        from myfempy.lib.beam_posproc import mesh_def_axial_rigid
        
        meshDefU,Udef,Umag = mesh_def_axial_rigid(datamesh,coord,U,scale_mesh)
        
        if np.any([[graphout == 'displ{ux}'],[graphout == 'displ{uy}'],[graphout == 'displ{mag}']]):
            screen_mode_shape = 'true'
            file_dir=path_user+'/'+usr_analysi_name+'_result_view_'+datastep+'.vtk'
            file_save=path_user+'/'+usr_analysi_name+'_posproc_'+datastep+'.png'
            
            result_list = np.concatenate((meshDefU,Umag),axis=1)
            
            if np.any(graphout == 'displ{ux}'):
                title_scale_bar = 'DISPL X'
                data_result_view = np.reshape(Udef[:,0],(-1,1))
            
            if np.any(graphout == 'displ{uy}'):
                title_scale_bar = 'DISPL Y'
                data_result_view = np.reshape(Udef[:,1],(-1,1))
                
            if np.any(graphout == 'displ{mag}'):
                title_scale_bar = 'DISPL MAG'
                data_result_view = Umag 
            
            data_result_save = np.concatenate((Udef,Umag),axis=1)   
            data_title = ['DISPL_X','DISPL_Y','DISPL_Z','DISPL_MAG']
            title_win = 'NODAL DISPLACEMENTS '+datastep
            scale_bar = 'true'
            title_scale_bar = 'DISPLACEMENT'
            title_display = ["MAX DISPL. U: ","MID DISPL. U: ","MIN DISPL. U: ",""]
            typeData='avr'
            
            data_display = np.array([[max(abs(data_result_view[:,0]))],[sum(data_result_view[:,0])/datamesh[1]],[min(data_result_view[:,0])],[0],[usrlog.mod_opt]])
            save_data(path_user,usr_analysi_name,file_dir,file_save,datamesh,inci,coord,type_elm,data_display,Udef,meshDefU,scale_mesh,data_result_save,data_result_view,data_title,title_display,typeData,scale_bar,title_scale_bar,screen_on,title_win,save_screen,screen_mode_shape,ModeNumb)
            screen_mode_shape = 'false'
            
        if np.any(graphout == 'mode'):
            
                result_list = U2
                
                ModeNumb = int(U2[0])
                Wn_rad = round(U2[1],4)
                Wn_hz = round(U2[2],4)
                
                screen_mode_shape = 'true'
                file_dir=path_user+'/'+usr_analysi_name+'_mode_shape_'+str(Wn_hz)+'_Hz_'+datastep+'.vtk'
                file_save=path_user+'/'+usr_analysi_name+'_mode_shape_'+str(Wn_hz)+'_Hz_'+datastep+'.png'
                
                data_result_save = Umag
                data_result_view = data_result_save
                data_title = ['MODE'+str(ModeNumb)+'_'+str(Wn_hz)+'_Hz']
                scale_bar = 'no'
                title_scale_bar = 'DISPL'
                title_display = ["MODE NUMB.: ","FREQ. RAD/S: ","FREQ. HZ: ","MAX DISPL.: "]
                title_win = 'MODE SHAPE'
                typeData='avr'
                
                data_display = np.array([[ModeNumb],[Wn_rad],[Wn_hz],[round(max(abs(data_result_view[:,0])),4)],[usrlog.mod_opt]])
                save_data(path_user,usr_analysi_name,file_dir,file_save,datamesh,inci,coord,type_elm,data_display,Udef,meshDefU,scale_mesh,data_result_save,data_result_view,data_title,title_display,typeData,scale_bar,title_scale_bar,screen_on,title_win,save_screen,screen_mode_shape,ModeNumb)
                screen_mode_shape = 'false'
    

    elif usrlog.mod_opt == "truss22":
        from myfempy.lib.beam_posproc import mesh_def_axial_rigid,react_suport,stress_axial_only_elm
        
        meshDefU,Udef,Umag = mesh_def_axial_rigid(datamesh,coord,U,scale_mesh)
        strs_elm = stress_axial_only_elm(datamesh,U,inci,coord,typeMechMat,tabmat)
        stress_list = strs_elm
        if np.any([[graphout == 'displ{ux}'],[graphout == 'displ{uy}'],[graphout == 'displ{mag}']]):
            screen_mode_shape = 'true'
            file_dir=path_user+'/'+usr_analysi_name+'_result_view_'+datastep+'.vtk'
            file_save=path_user+'/'+usr_analysi_name+'_posproc_'+datastep+'.png'
            
            result_list = np.concatenate((meshDefU,Umag),axis=1)
        
            if np.any(graphout == 'displ{ux}'):
                title_scale_bar = 'DISPL X'
                data_result_view = np.reshape(Udef[:,0],(-1,1))
            
            if np.any(graphout == 'displ{uy}'):
                title_scale_bar = 'DISPL Y'
                data_result_view = np.reshape(Udef[:,1],(-1,1))
                
            if np.any(graphout == 'displ{mag}'):
                title_scale_bar = 'DISPL MAG'
                data_result_view = Umag 
                
            data_result_save = np.concatenate((Udef,Umag),axis=1)   
            data_title = ['DISPL_X','DISPL_Y','DISPL_Z','DISPL_MAG']    
            title_win = 'NODAL DISPLACEMENTS'+datastep
            scale_bar = 'true'
            title_display = ["MAX DISPL. U: ","MID DISPL. U: ","MIN DISPL. U: ",""]
            typeData='avr'
            
            data_display = np.array([[max(abs(data_result_view[:,0]))],[sum(data_result_view[:,0])/datamesh[1]],[min(data_result_view[:,0])],[0],[usrlog.mod_opt]])
            save_data(path_user,usr_analysi_name,file_dir,file_save,datamesh,inci,coord,type_elm,data_display,Udef,meshDefU,scale_mesh,data_result_save,data_result_view,data_title,title_display,typeData,scale_bar,title_scale_bar,screen_on,title_win,save_screen,screen_mode_shape,ModeNumb)
            screen_mode_shape = 'false'
            
            if np.any(graphout == 'steptracker'):
                from myfempy.lib.domain_physics import search_nodeXYZ
                
                for jj in range(len(hist_dof_list)):
                    node_coordX = float(hist_dof_list[jj,0])
                    node_coordY = float(hist_dof_list[jj,1])
                    node_coordZ = float(hist_dof_list[jj,2])
                    hist_node_list = search_nodeXYZ(node_coordX,node_coordY,node_coordZ,coord,2E-3)
                    
                    max_displ = Umag[hist_node_list-1,0]
            
        if np.any(graphout == 'stress_elem{max}'):
            file_dir=path_user+'/'+usr_analysi_name+'_result_view_'+datastep+'_.vtk'
            file_save=path_user+'/'+usr_analysi_name+'_posproc_'+datastep+'_.png'
                        
            data_result_save  = strs_elm
            data_result_view = data_result_save
            data_title = ['STRESS XX']
            title_win = 'TRESS AXIAL'
            scale_bar = 'true'
            title_scale_bar = 'STRESS MAX'
            title_display = ["MAX DISPL. U: ","MAX STRESS: ","MID STRESS: ","MIN STRESS.: "]
            typeData='elm'
        
            data_display = np.array([[max(abs(data_result_view[:,0]))],[sum(data_result_view[:,0])/datamesh[1]],[min(data_result_view[:,0])],[0],[usrlog.mod_opt]])
            save_data(path_user,usr_analysi_name,file_dir,file_save,datamesh,inci,coord,type_elm,data_display,Udef,meshDefU,scale_mesh,data_result_save,data_result_view,data_title,title_display,typeData,scale_bar,title_scale_bar,screen_on,title_win,save_screen,screen_mode_shape,ModeNumb)
         
        if np.any(graphout == 'beam_extforc'):
            
            Fr,Mr = react_suport(datamesh,coord,KG,U)
            # print('FORCA {0:.4f} | MOMENTO {0:.4f}\n'.format(Fr,Mr))
            print('FORCA \n'+str(Fr))
            print('\n')
            print('MOMENTO \n'+str(Mr))
            data_result = np.zeros([datamesh[2],1])
            result_name = 'DEFORMATION'
            title_display = ["MAX DISPL. U: ","MAX STRESS: ","MID STRESS: ","MIN STRESS.: "]
            typeData='elm'
                
            
    elif usrlog.mod_opt == "beam21":
            from myfempy.lib.beam_posproc import mesh_def_bending_2d,intn_force_bending_only,react_suport,stress_bending_only_elm
            
            meshDefU,meshRotZ,Udef,Umag = mesh_def_bending_2d(datamesh,coord,U,scale_mesh)
            strs_elm_max = stress_bending_only_elm(datamesh,U,inci,coord,typeMechMat,tabmat,tabgeo,y_max)
            strs_elm_min = stress_bending_only_elm(datamesh,U,inci,coord,typeMechMat,tabmat,tabgeo,y_min)
            stress_list = np.concatenate((strs_elm_max,strs_elm_min),axis=1)                        
            
            if np.any([[graphout == 'displ{ux}'],[graphout == 'displ{uy}'],[graphout == 'displ{mag}']]):
                screen_mode_shape = 'true'
                file_dir=path_user+'/'+usr_analysi_name+'_displacement_'+datastep+'.vtk'
                file_save=path_user+'/'+usr_analysi_name+'_displacement_'+datastep+'.png'
                
                result_list = np.concatenate((meshDefU,Umag),axis=1)
            
                if np.any(graphout == 'displ{ux}'):
                    title_scale_bar = 'DISPL X'
                    data_result_view = np.reshape(Udef[:,0],(-1,1))
                
                if np.any(graphout == 'displ{uy}'):
                    title_scale_bar = 'DISPL Y'
                    data_result_view = np.reshape(Udef[:,1],(-1,1))
                    
                if np.any(graphout == 'displ{mag}'):
                    title_scale_bar = 'DISPL MAG'
                    data_result_view = Umag 
                
                data_result_save = np.concatenate((Udef,Umag),axis=1)   
                data_title = ['DISPL_X','DISPL_Y','DISPL_Z','DISPL_MAG']
                title_win = 'NODAL DISPLACEMENTS '+datastep
                scale_bar = 'true'
                title_display = ["MAX DISPL. U: ","MID DISPL. U: ","MIN DISPL. U: ",""]
                typeData='avr'
            
                data_display = np.array([[max(abs(data_result_view[:,0]))],[sum(data_result_view[:,0])/datamesh[1]],[min(data_result_view[:,0])],[0],[usrlog.mod_opt]])
                save_data(path_user,usr_analysi_name,file_dir,file_save,datamesh,inci,coord,type_elm,data_display,Udef,meshDefU,scale_mesh,data_result_save,data_result_view,data_title,title_display,typeData,scale_bar,title_scale_bar,screen_on,title_win,save_screen,screen_mode_shape,ModeNumb)
                screen_mode_shape = 'false'
            
                if np.any(graphout == 'steptracker'):
                    from myfempy.lib.domain_physics import search_nodeXYZ
                    
                    for jj in range(len(hist_dof_list)):
                        node_coordX = float(hist_dof_list[jj,0])
                        node_coordY = float(hist_dof_list[jj,1])
                        node_coordZ = float(hist_dof_list[jj,2])
                        hist_node_list = search_nodeXYZ(node_coordX,node_coordY,node_coordZ,coord,2E-3)
                        
                        max_displ = Umag[hist_node_list-1,0]
            
            if np.any([[graphout == 'stress_elem{max}'],[graphout == 'stress_elem{min}']]):           
                
                data_result_save = np.concatenate((strs_elm_max,strs_elm_min),axis=1)
                
                if np.any(graphout == 'stress_elem{max}'):
                    title_scale_bar = 'STRESS XX MAX'
                    data_result_view = strs_elm_max
                
                if np.any(graphout == 'stress_elem{min}'):
                    title_scale_bar = 'STRESS XX MIN'
                    data_result_view = strs_elm_min
                
                
                data_title = ['STRXX_MAX','STRXX_MIN']
                title_win = 'STRESS BENDING'
                scale_bar = 'true'
                title_display = ["MAX DISPL. U: ","MAX STRESS: ","MID STRESS: ","MIN STRESS.: "]
                typeData='elm'
                
                
                file_dir=path_user+'/'+usr_analysi_name+'_stress_flexural_'+datastep+'.vtk'
                file_save=path_user+'/'+usr_analysi_name+'_stress_flexural_'+datastep+'.png'     
                
                data_display = np.array([[max(abs(U))[0]],[max(data_result_view[:,0])],[sum(data_result_view[:,0])/datamesh[1]],[min(data_result_view[:,0])],[usrlog.mod_opt]])
                save_data(path_user,usr_analysi_name,file_dir,file_save,datamesh,inci,coord,type_elm,data_display,Udef,meshDefU,scale_mesh,data_result_save,data_result_view,data_title,title_display,typeData,scale_bar,title_scale_bar,screen_on,title_win,save_screen,screen_mode_shape,ModeNumb)
                

            if np.any(graphout == 'beam_intforc'):
                from myfempy.solve.view_mesh import plotterXY_resultgraph
                
                domL,Vy,Mz = intn_force_bending_only(datamesh,coord,inci,typeMechMat,tabmat,tabgeo,U)
                                
                Mz_max_idx = np.argmax(abs(Mz), axis=0)
                Mz_min_idx = np.argmin(abs(Mz), axis=0)
                Vy_max_idx = np.argmax(abs(Vy), axis=0)
                Vy_min_idx = np.argmin(abs(Vy), axis=0)
                
                Mz_max = Mz[Mz_max_idx,0]
                Mz_min = Mz[Mz_min_idx,0]
                Vy_max = Vy[Vy_max_idx,0]
                Vy_min = Vy[Vy_min_idx,0]
                
                coordX_Mz_max = domL[Mz_max_idx,0]
                coordX_Mz_min = domL[Mz_min_idx,0]
                coordX_Vy_max = domL[Vy_max_idx,0]
                coordX_Vy_min = domL[Vy_min_idx,0]
                
                Xlabel = 'L_beam'
                Ylabel1 = 'VY'
                Ylabel2 = 'MZ'
                plot_name = 'momento_fletor'
                X = domL
                numPoints = datamesh[1]
                
                file_save1=path_user+'/'+usr_analysi_name+plot_name+'_plotter_Vy_'+datastep+'.png'
                file_save2=path_user+'/'+usr_analysi_name+plot_name+'_plotter_Mz_'+datastep+'.png'

                plotterXY_resultgraph(file_save1,save_screen,screen_on,numPoints,X,Vy,Xlabel,Ylabel1)
                plotterXY_resultgraph(file_save2,save_screen,screen_on,numPoints,X,Mz,Xlabel,Ylabel2)
             
            if np.any(graphout == 'beam_extforc'):
                Fr,Mr = react_suport(datamesh,coord,KG,U)
                # print('FORCA {0:.4f} | MOMENTO {0:.4f}\n'.format(Fr,Mr))
                print('FORCA \n'+str(Fr))
                print('\n')
                print('MOMENTO \n'+str(Mr))
                data_result = np.zeros([datamesh[2],1])
                result_name = 'DEFORMATION'
                title_display = ["MAX DISPL. U: ","MAX STRESS: ","MID STRESS: ","MIN STRESS.: "]
                typeData='elm'
                
            if np.any(graphout == 'mode'):
                
                result_list = U2
                
                ModeNumb = int(U2[0])
                Wn_rad = round(U2[1],4)
                Wn_hz = round(U2[2],4)
                
                screen_mode_shape = 'true'
                file_dir=path_user+'/'+usr_analysi_name+'_mode_shape_'+str(Wn_hz)+'_Hz_'+datastep+'.vtk'
                file_save=path_user+'/'+usr_analysi_name+'_mode_shape_'+str(Wn_hz)+'_Hz_'+datastep+'.png'
                
                data_result_save = Umag
                data_result_view = data_result_save
                data_title = ['MODE'+str(ModeNumb)+'_'+str(Wn_hz)+'_Hz']
                scale_bar = 'false'
                title_scale_bar = 'DISPL'
                title_display = ["MODE NUMB.: ","FREQ. RAD/S: ","FREQ. HZ: ","MAX DISPL.: "]
                title_win = 'MODE SHAPE'
                typeData='avr'
                
                data_display = np.array([[ModeNumb],[Wn_rad],[Wn_hz],[round(max(abs(data_result_view[:,0])),4)],[usrlog.mod_opt]])
                save_data(path_user,usr_analysi_name,file_dir,file_save,datamesh,inci,coord,type_elm,data_display,Udef,meshDefU,scale_mesh,data_result_save,data_result_view,data_title,title_display,typeData,scale_bar,title_scale_bar,screen_on,title_win,save_screen,screen_mode_shape,ModeNumb)
                screen_mode_shape = 'false'
   
            if np.any(graphout == 'frf'):
                from myfempy.lib.domain_physics import search_nodeXYZ
                
                for jj in range(len(hist_dof_list)):
                    node_coordX = float(hist_dof_list[jj,0])
                    node_coordY = float(hist_dof_list[jj,1])
                    node_coordZ = float(hist_dof_list[jj,2])
                    hist_node_list = search_nodeXYZ(node_coordX,node_coordY,node_coordZ,coord,2E-3)
                    node_dof_hist = datamesh[0]*hist_node_list-1
                    
                    from myfempy.solve.view_mesh import plotterXY_resultgraph
            
                    frf = U3[node_dof_hist,:]
                
                    Xlabel = 'FREQUENCY [Hz]'
                    Ylabel = 'FRF NODE '+str(node_dof_hist)
                    plot_name = 'Response in Frequency'
                    
                    X = np.ravel(U2)
                    Y = np.log(abs(np.ravel(frf.todense())))
                    numPoints = U2.shape[0]
                    
                    file_save=path_user+'/'+usr_analysi_name+plot_name+'_plotter_FRF_'+datastep+str(jj)+'.png'
                    plotterXY_resultgraph(file_save,save_screen,screen_on,numPoints,X,Y,Xlabel,Ylabel)
   
    
    elif usrlog.mod_opt == "frame22":
            from myfempy.lib.beam_posproc import mesh_def_flexural_2d, intn_force_flexural_full_2d,stress_flexural_full_2d_elm
            
            meshDefU,meshRotZ,Udef,Umag = mesh_def_flexural_2d(datamesh,coord,U,scale_mesh) 
            strs_elm_max = stress_flexural_full_2d_elm(datamesh,U,inci,coord,typeMechMat,tabmat,tabgeo,y_max)
            strs_elm_min = stress_flexural_full_2d_elm(datamesh,U,inci,coord,typeMechMat,tabmat,tabgeo,y_min)
            stress_list = np.concatenate((strs_elm_max,strs_elm_min),axis=1)            
            if np.any([[graphout == 'displ{ux}'],[graphout == 'displ{uy}'],[graphout == 'displ{mag}']]):
                screen_mode_shape = 'true'
                file_dir=path_user+'/'+usr_analysi_name+'_displacement_'+datastep+'.vtk'
                file_save=path_user+'/'+usr_analysi_name+'_displacement_'+datastep+'.png'
                
                result_list = np.concatenate((meshDefU,Umag),axis=1)
            
                if np.any(graphout == 'displ{ux}'):
                    title_scale_bar = 'DISPL X'
                    data_result_view = np.reshape(Udef[:,0],(-1,1))
                
                if np.any(graphout == 'displ{uy}'):
                    title_scale_bar = 'DISPL Y'
                    data_result_view = np.reshape(Udef[:,1],(-1,1))
                    
                if np.any(graphout == 'displ{mag}'):
                    title_scale_bar = 'DISPL MAG'
                    data_result_view = Umag 
                
                data_result_save = np.concatenate((Udef,Umag),axis=1)   
                data_title = ['DISPL_X','DISPL_Y','DISPL_Z','DISPL_MAG']
                scale_bar = 'true'
                title_scale_bar = 'STRESS'
                result_name = 'DISPLACEMENT'
                title_win = 'NODAL DISPLACEMENTS'+datastep
                title_display = ["MAX DISPL. U: ","MID DISPL. U: ","MIN DISPL. U: ",""]
                typeData='avr'
                
                data_display = np.array([[max(abs(data_result_view[:,0]))],[sum(data_result_view[:,0])/datamesh[1]],[min(data_result_view[:,0])],[0],[usrlog.mod_opt]])
                save_data(path_user,usr_analysi_name,file_dir,file_save,datamesh,inci,coord,type_elm,data_display,Udef,meshDefU,scale_mesh,data_result_save,data_result_view,data_title,title_display,typeData,scale_bar,title_scale_bar,screen_on,title_win,save_screen,screen_mode_shape,ModeNumb)
                screen_mode_shape = 'false'
                
                if np.any(graphout == 'steptracker'):
                    from myfempy.lib.domain_physics import search_nodeXYZ
                    
                    for jj in range(len(hist_dof_list)):
                        node_coordX = float(hist_dof_list[jj,0])
                        node_coordY = float(hist_dof_list[jj,1])
                        node_coordZ = float(hist_dof_list[jj,2])
                        hist_node_list = search_nodeXYZ(node_coordX,node_coordY,node_coordZ,coord,2E-3)
                        
                        max_displ = Umag[hist_node_list-1,0]
                

            if np.any([[graphout == 'stress_elem{max}'],[graphout == 'stress_elem{min}']]):                
                
                data_result_save = np.concatenate((strs_elm_max,strs_elm_min),axis=1)
                
                if np.any(graphout == 'stress_elem{max}'):
                    title_scale_bar = 'STRESS XX MAX'
                    data_result_view = strs_elm_max
                
                if np.any(graphout == 'stress_elem{min}'):
                    title_scale_bar = 'STRESS XX MIN'
                    data_result_view = strs_elm_min
                
                data_title = ['STRXX_MAX','STRXX_MIN']
                title_win = 'STRESS FLEXURAL'
                scale_bar = 'true'
                title_display = ["MAX DISPL. U: ","MAX STRESS: ","MID STRESS: ","MIN STRESS.: "]
                typeData='elm'
                
                
                file_dir=path_user+'/'+usr_analysi_name+'_stress_flexural_'+datastep+'.vtk'
                file_save=path_user+'/'+usr_analysi_name+'_stress_flexural_'+datastep+'.png'     
                
                data_display = np.array([[max(abs(U))[0]],[max(data_result_view[:,0])],[sum(data_result_view[:,0])/datamesh[1]],[min(data_result_view[:,0])],[usrlog.mod_opt]])
                save_data(path_user,usr_analysi_name,file_dir,file_save,datamesh,inci,coord,type_elm,data_display,Udef,meshDefU,scale_mesh,data_result_save,data_result_view,data_title,title_display,typeData,scale_bar,title_scale_bar,screen_on,title_win,save_screen,screen_mode_shape,ModeNumb)
                
            if np.any(graphout == 'beam_intforc'):
                from myfempy.solve.view_mesh import plotterXY_resultgraph
                
                domL,Nx,Vy,Mz = intn_force_flexural_full_2d(datamesh,coord,inci,typeMechMat,tabmat,tabgeo,U)
                
                # data_result = np.concatenate((np.ones((datamesh[1],1)),2*np.ones((datamesh[1],1)),3*np.ones((datamesh[1],1))),axis=1)
                # data_title = ['Nx','Vy','Mz']
                
                # title_win = 'BEAM FORCE INTERNAL'
                # scale_bar = 'no'
                # title_scale_bar = ''
                # title_display = ["MAX DISPL. U: ","MAX STRESS: ","MID STRESS: ","MIN STRESS.: "]
                # typeData='avr'
                
                
                # file_dir1=path_user+'/'+usr_analysi_name+'_Nx.vtk'
                # file_save1=path_user+'/'+usr_analysi_name+'_Nx.png'   
                # file_dir2=path_user+'/'+usr_analysi_name+'_Vy.vtk'
                # file_save2=path_user+'/'+usr_analysi_name+'_Vy.png' 
                # file_dir1=path_user+'/'+usr_analysi_name+'_Mz.vtk'
                # file_save1=path_user+'/'+usr_analysi_name+'_Mz.png' 
                
                # data_display = np.array([[1],[1],[1],[1],[usrlog.mod_opt]])
                

                # save_data(path_user,usr_analysi_name,file_dir1,file_save2,datamesh,inci,coord,type_elm,data_display,Udef,meshBeamNx,scale_mesh,data_result,data_title,title_display,typeData,scale_bar,title_scale_bar,screen_on,title_win,save_screen,screen_mode_shape,ModeNumb)
                # save_data(path_user,usr_analysi_name,file_dir1,file_save2,datamesh,inci,coord,type_elm,data_display,Udef,meshBeamVy,scale_mesh,data_result,data_title,title_display,typeData,scale_bar,title_scale_bar,screen_on,title_win,save_screen,screen_mode_shape,ModeNumb)                
                # save_data(path_user,usr_analysi_name,file_dir1,file_save2,datamesh,inci,coord,type_elm,data_display,Udef,meshBeamMz,scale_mesh,data_result,data_title,title_display,typeData,scale_bar,title_scale_bar,screen_on,title_win,save_screen,screen_mode_shape,ModeNumb)


                Mz_max_idx = np.argmax(abs(Mz), axis=0)
                Mz_min_idx = np.argmin(abs(Mz), axis=0)
                Vy_max_idx = np.argmax(abs(Vy), axis=0)
                Vy_min_idx = np.argmin(abs(Vy), axis=0)
                
                Mz_max = Mz[Mz_max_idx,0]
                Mz_min = Mz[Mz_min_idx,0]
                Vy_max = Vy[Vy_max_idx,0]
                Vy_min = Vy[Vy_min_idx,0]
                
                coordX_Mz_max = domL[Mz_max_idx,0]
                coordX_Mz_min = domL[Mz_min_idx,0]
                coordX_Vy_max = domL[Vy_max_idx,0]
                coordX_Vy_min = domL[Vy_min_idx,0]
                
                Xlabel = 'L_beam'
                Ylabel1 = 'NX'
                Ylabel2 = 'VY'
                Ylabel3 = 'MZ'
                plot_name = 'momento_fletor'
                X = domL
                numPoints = datamesh[1]
                
                file_save1=path_user+'/'+usr_analysi_name+plot_name+'_plotter_Nx_'+datastep+'.png'
                file_save2=path_user+'/'+usr_analysi_name+plot_name+'_plotter_Vy_'+datastep+'.png'
                file_save3=path_user+'/'+usr_analysi_name+plot_name+'_plotter_Mz_'+datastep+'.png'

                plotterXY_resultgraph(file_save1,save_screen,screen_on,numPoints,X,Nx,Xlabel,Ylabel1)
                plotterXY_resultgraph(file_save2,save_screen,screen_on,numPoints,X,Vy,Xlabel,Ylabel2)
                plotterXY_resultgraph(file_save2,save_screen,screen_on,numPoints,X,Mz,Xlabel,Ylabel3)
                
             
            if np.any(graphout == 'beam_extforc'):
                Fr,Mr = react_suport(datamesh,coord,KG,U)
                # print('FORCA {0:.4f} | MOMENTO {0:.4f}\n'.format(Fr,Mr))
                print('FORCA \n'+str(Fr))
                print('\n')
                print('MOMENTO \n'+str(Mr))
                data_result = np.zeros([datamesh[2],1])
                result_name = 'DEFORMATION'
                title_display = ["MAX DISPL. U: ","MAX STRESS: ","MID STRESS: ","MIN STRESS.: "]
                typeData='elm'
                
            if np.any(graphout == 'mode'):
                
                ModeNumb = int(U2[0])
                Wn_rad = round(U2[1],4)
                Wn_hz = round(U2[2],4)
                
                screen_mode_shape = 'true'
                file_dir=path_user+'/'+usr_analysi_name+'_mode_shape_'+str(Wn_hz)+'_Hz_'+datastep+'.vtk'
                file_save=path_user+'/'+usr_analysi_name+'_mode_shape_'+str(Wn_hz)+'_Hz_'+datastep+'.png'
                
                data_result_save = Umag
                data_result_view = data_result_save
                data_title = ['MODE'+str(ModeNumb)+'_'+str(Wn_hz)+'_Hz']
                scale_bar = 'false'
                title_scale_bar = 'DISPL'
                title_display = ["MODE NUMB.: ","FREQ. RAD/S: ","FREQ. HZ: ","MAX DISPL.: "]
                title_win = 'MODE SHAPE'
                typeData='avr'
                
                data_display = np.array([[ModeNumb],[Wn_rad],[Wn_hz],[round(max(abs(data_result_view[:,0])),4)],[usrlog.mod_opt]])
                save_data(path_user,usr_analysi_name,file_dir,file_save,datamesh,inci,coord,type_elm,data_display,Udef,meshDefU,scale_mesh,data_result_save,data_result_view,data_title,title_display,typeData,scale_bar,title_scale_bar,screen_on,title_win,save_screen,screen_mode_shape,ModeNumb)
                screen_mode_shape = 'false'
                
            if np.any(graphout == 'frf'):
                from myfempy.lib.domain_physics import search_nodeXYZ
                
                for jj in range(len(hist_dof_list)):
                    node_coordX = float(hist_dof_list[jj,0])
                    node_coordY = float(hist_dof_list[jj,1])
                    node_coordZ = float(hist_dof_list[jj,2])
                    hist_node_list = search_nodeXYZ(node_coordX,node_coordY,node_coordZ,coord,2E-3)
                    node_dof_hist = datamesh[0]*hist_node_list-1
                    
                    from myfempy.solve.view_mesh import plotterXY_resultgraph
            
                    frf = U3[node_dof_hist,:]
                
                    Xlabel = 'FREQUENCY [Hz]'
                    Ylabel = 'FRF NODE '+str(node_dof_hist)
                    plot_name = 'Response in Frequency'
                    
                    X = np.ravel(U2)
                    Y = np.log(abs(np.ravel(frf.todense())))
                    numPoints = U2.shape[0]
                    
                    file_save=path_user+'/'+usr_analysi_name+plot_name+'_plotter_FRF_'+datastep+str(jj)+'.png'
                    plotterXY_resultgraph(file_save,save_screen,screen_on,numPoints,X,Y,Xlabel,Ylabel)   
            
            
    elif usrlog.mod_opt == "frame23":
            from myfempy.lib.beam_posproc import mesh_def_flexural_3d,intn_force_flexural_full_3d,stress_flexural_full_3d_elm
            
            meshDefU,meshRotZ,Udef,Umag = mesh_def_flexural_3d(datamesh,coord,U,scale_mesh)
            strs_elm_max = stress_flexural_full_3d_elm(datamesh,U,inci,coord,typeMechMat,tabmat,tabgeo,y_max,z_max,r_max)
            strs_elm_min = stress_flexural_full_3d_elm(datamesh,U,inci,coord,typeMechMat,tabmat,tabgeo,y_min,z_min,r_max)
            stress_list = np.concatenate((strs_elm_max,strs_elm_min),axis=1)
                        
            if np.any([[graphout == 'displ{ux}'],[graphout == 'displ{uy}'],[graphout == 'displ{mag}']]):
                
                file_dir=path_user+'/'+usr_analysi_name+'_displacement_'+datastep+'.vtk'
                file_save=path_user+'/'+usr_analysi_name+'_displacement_'+datastep+'.png'
                
                result_list = np.concatenate((meshDefU,Umag),axis=1)
                            
                if np.any(graphout == 'displ{ux}'):
                    title_scale_bar = 'DISPL X'
                    data_result_view = np.reshape(Udef[:,0],(-1,1))
                
                if np.any(graphout == 'displ{uy}'):
                    title_scale_bar = 'DISPL Y'
                    data_result_view = np.reshape(Udef[:,1],(-1,1))
                    
                if np.any(graphout == 'displ{mag}'):
                    title_scale_bar = 'DISPL MAG'
                    data_result_view = Umag 
                
                data_result_save = np.concatenate((Udef,Umag),axis=1)   
                data_title = ['DISPL_X','DISPL_Y','DISPL_Z','DISPL_MAG']

                scale_bar = 'true'
                title_scale_bar = 'DISPLA'
                result_name = 'DISPLACEMENT'
                title_win = 'NODAL DISPLACEMENTS'+datastep
                title_display = ["MAX DISPL. U: ","MID DISPL. U: ","MIN DISPL. U: ",""]
                typeData='avr'
                
                data_display = np.array([[max(abs(data_result_view[:,0]))],[sum(data_result_view[:,0])/datamesh[1]],[min(data_result_view[:,0])],[0],[usrlog.mod_opt]])
                save_data(path_user,usr_analysi_name,file_dir,file_save,datamesh,inci,coord,type_elm,data_display,Udef,meshDefU,scale_mesh,data_result_save,data_result_view,data_title,title_display,typeData,scale_bar,title_scale_bar,screen_on,title_win,save_screen,screen_mode_shape,ModeNumb)
                screen_mode_shape = 'false'
                
                if np.any(graphout == 'steptracker'):
                    from myfempy.lib.domain_physics import search_nodeXYZ
                    
                    for jj in range(len(hist_dof_list)):
                        node_coordX = float(hist_dof_list[jj,0])
                        node_coordY = float(hist_dof_list[jj,1])
                        node_coordZ = float(hist_dof_list[jj,2])
                        hist_node_list = search_nodeXYZ(node_coordX,node_coordY,node_coordZ,coord,2E-3)
                        
                        max_displ = Umag[hist_node_list-1,0]
                
                
            if np.any([[graphout == 'stress_elem{max}'],[graphout == 'stress_elem{min}']]):                
                
                data_result_save = np.concatenate((strs_elm_max,strs_elm_min),axis=1)
                
                if np.any(graphout == 'stress_elem{max}'):
                    title_scale_bar = 'STRESS XX MAX'
                    data_result_view = strs_elm_max
                
                if np.any(graphout == 'stress_elem{min}'):
                    title_scale_bar = 'STRESS XX MIN'
                    data_result_view = strs_elm_min
                
                data_title = ['STRXX_MAX','STRXX_MIN']
                title_win = 'STRESS FLEXURAL'
                scale_bar = 'true'
                title_display = ["MAX DISPL. U: ","MAX STRESS: ","MID STRESS: ","MIN STRESS.: "]
                typeData='elm'
                
                file_dir=path_user+'/'+usr_analysi_name+'_stress_flexural_'+datastep+'.vtk'
                file_save=path_user+'/'+usr_analysi_name+'_stress_flexural_'+datastep+'.png'     
                
                data_display = np.array([[max(abs(U))[0]],[max(data_result_view[:,0])],[sum(data_result_view[:,0])/datamesh[1]],[min(data_result_view[:,0])],[usrlog.mod_opt]])
                save_data(path_user,usr_analysi_name,file_dir,file_save,datamesh,inci,coord,type_elm,data_display,Udef,meshDefU,scale_mesh,data_result_save,data_result_view,data_title,title_display,typeData,scale_bar,title_scale_bar,screen_on,title_win,save_screen,screen_mode_shape,ModeNumb)

                
            if np.any(graphout == 'beam_intforc'):
                from myfempy.solve.view_mesh import plotterXY_resultgraph
                
                domL,Nx,Vy,Vz,Tx,My,Mz = intn_force_flexural_full_3d(datamesh,coord,inci,typeMechMat,tabmat,tabgeo,U)
                
                Mz_max_idx = np.argmax(abs(Mz), axis=0)
                Mz_min_idx = np.argmin(abs(Mz), axis=0)
                Vy_max_idx = np.argmax(abs(Vy), axis=0)
                Vy_min_idx = np.argmin(abs(Vy), axis=0)
                
                Mz_max = Mz[Mz_max_idx,0]
                Mz_min = Mz[Mz_min_idx,0]
                Vy_max = Vy[Vy_max_idx,0]
                Vy_min = Vy[Vy_min_idx,0]
                
                coordX_Mz_max = domL[Mz_max_idx,0]
                coordX_Mz_min = domL[Mz_min_idx,0]
                coordX_Vy_max = domL[Vy_max_idx,0]
                coordX_Vy_min = domL[Vy_min_idx,0]
                
                Xlabel = 'L_beam'
                Ylabel1 = 'NX'
                Ylabel2 = 'VY'
                Ylabel3 = 'VZ'
                Ylabel4 = 'TX'
                Ylabel5 = 'MY'
                Ylabel6 = 'MZ'
                plot_name = 'momento_fletor'
                X = domL
                numPoints = datamesh[1]
                
                file_save1=path_user+'/'+usr_analysi_name+plot_name+'_plotter_Nx_'+datastep+'.png'
                file_save2=path_user+'/'+usr_analysi_name+plot_name+'_plotter_Vy_'+datastep+'.png'
                file_save3=path_user+'/'+usr_analysi_name+plot_name+'_plotter_Vz_'+datastep+'.png'
                file_save4=path_user+'/'+usr_analysi_name+plot_name+'_plotter_Tx_'+datastep+'.png'
                file_save5=path_user+'/'+usr_analysi_name+plot_name+'_plotter_My_'+datastep+'.png'
                file_save6=path_user+'/'+usr_analysi_name+plot_name+'_plotter_Mz_'+datastep+'.png'

                plotterXY_resultgraph(file_save1,save_screen,screen_on,numPoints,X,Nx,Xlabel,Ylabel1)
                plotterXY_resultgraph(file_save2,save_screen,screen_on,numPoints,X,Vy,Xlabel,Ylabel2)
                plotterXY_resultgraph(file_save3,save_screen,screen_on,numPoints,X,Vz,Xlabel,Ylabel3)
                plotterXY_resultgraph(file_save4,save_screen,screen_on,numPoints,X,Tx,Xlabel,Ylabel4)
                plotterXY_resultgraph(file_save5,save_screen,screen_on,numPoints,X,My,Xlabel,Ylabel5)
                plotterXY_resultgraph(file_save6,save_screen,screen_on,numPoints,X,Mz,Xlabel,Ylabel6)  
                
            if np.any(graphout == 'beam_extforc'):
                Fr,Mr = react_suport(datamesh,coord,KG,U)
                # print('FORCA {0:.4f} | MOMENTO {0:.4f}\n'.format(Fr,Mr))
                print('FORCA \n'+str(Fr))
                print('\n')
                print('MOMENTO \n'+str(Mr))
                data_result = np.zeros([datamesh[2],1])
                result_name = 'DEFORMATION'
                title_display = ["MAX DISPL. U: ","MAX STRESS: ","MID STRESS: ","MIN STRESS.: "]
                typeData='elm'
                
            if np.any(graphout == 'mode'):
                
                ModeNumb = int(U2[0])
                Wn_rad = round(U2[1],4)
                Wn_hz = round(U2[2],4)
                
                screen_mode_shape = 'true'
                file_dir=path_user+'/'+usr_analysi_name+'_mode_shape_'+str(Wn_hz)+'_Hz_'+datastep+'.vtk'
                file_save=path_user+'/'+usr_analysi_name+'_mode_shape_'+str(Wn_hz)+'_Hz_'+datastep+'.png'
                
                data_result_save = Umag
                data_result_view = data_result_save
                data_title = ['MODE'+str(ModeNumb)+'_'+str(Wn_hz)+'_Hz']
                scale_bar = 'false'
                title_scale_bar = 'DISPL'
                title_display = ["MODE NUMB.: ","FREQ. RAD/S: ","FREQ. HZ: ","MAX DISPL.: "]
                title_win = 'MODE SHAPE'
                typeData='avr'
                
                data_display = np.array([[ModeNumb],[Wn_rad],[Wn_hz],[round(max(abs(data_result_view [:,0])),4)],[usrlog.mod_opt]])
                save_data(path_user,usr_analysi_name,file_dir,file_save,datamesh,inci,coord,type_elm,data_display,Udef,meshDefU,scale_mesh,data_result_save,data_result_view,data_title,title_display,typeData,scale_bar,title_scale_bar,screen_on,title_win,save_screen,screen_mode_shape,ModeNumb)
                screen_mode_shape = 'false'
                
            if np.any(graphout == 'frf'):
                from myfempy.lib.domain_physics import search_nodeXYZ
                
                for jj in range(len(hist_dof_list)):
                    node_coordX = float(hist_dof_list[jj,0])
                    node_coordY = float(hist_dof_list[jj,1])
                    node_coordZ = float(hist_dof_list[jj,2])
                    hist_node_list = search_nodeXYZ(node_coordX,node_coordY,node_coordZ,coord,2E-3)
                    node_dof_hist = datamesh[0]*hist_node_list-1
                    
                    from myfempy.solve.view_mesh import plotterXY_resultgraph
            
                    frf = U3[node_dof_hist,:]
                
                    Xlabel = 'FREQUENCY [Hz]'
                    Ylabel = 'FRF NODE '+str(node_dof_hist)
                    plot_name = 'Response in Frequency'
                    
                    X = np.ravel(U2)
                    Y = np.log(abs(np.ravel(frf.todense())))
                    numPoints = U2.shape[0]
                    
                    file_save=path_user+'/'+usr_analysi_name+plot_name+'_plotter_FRF_'+datastep+str(jj)+'.png'
                    plotterXY_resultgraph(file_save,save_screen,screen_on,numPoints,X,Y,Xlabel,Ylabel) 
            
    else:
        print("input erro: elm_opt don't defined")

    return max_displ, max_stress, max_strain, result_list, stress_list
    
    
def plate_output(path_user,usr_analysi_name,usrlog,graphout,datamesh,inci,coord,tabmat,type_elm,KG,U,U2,U3,screen_on,save_screen,scale_mesh,datastep,hist_dof_list):
    
    screen_mode_shape = 'false'
    
    typeMechMat = [usrlog.mat_opt,usrlog.mat_def]
    
    max_displ = np.zeros((1))
    max_stress = np.zeros((1))
    max_strain = np.zeros((1))
    
    ModeNumb = int(datastep)
    if usrlog.mod_opt == "plane32":
        from myfempy.lib.plate_posproc import mesh_def_plane_2d, stress_vm_plane_t3,stress_avr_plane
        
        meshDefU, Udef, Umag = mesh_def_plane_2d(datamesh,coord,U,scale_mesh)
        stress_list, strain_list = stress_vm_plane_t3(U,inci,coord,typeMechMat,tabmat,datamesh)
                
        if np.any([[graphout == 'displ{ux}'],[graphout == 'displ{uy}'],[graphout == 'displ{mag}']]):
            
            result_list = np.concatenate((meshDefU,Umag),axis=1)
            
            if np.any(graphout == 'displ{ux}'):
                title_scale_bar = 'DISPL X'
                data_result_view = np.reshape(Udef[:,0],(-1,1))
            
            if np.any(graphout == 'displ{uy}'):
                title_scale_bar = 'DISPL Y'
                data_result_view = np.reshape(Udef[:,1],(-1,1))
                
            if np.any(graphout == 'displ{mag}'):
                title_scale_bar = 'DISPL MAG'
                data_result_view = Umag 
             
            data_result_save = np.concatenate((Udef,Umag),axis=1)    
            data_title = ['DISPL_X','DISPL_Y','DISPL_Z','DISPL_MAG']
            screen_mode_shape = 'true'
            file_dir=path_user+'/'+usr_analysi_name+'_displacement_'+datastep+'.vtk'
            file_save=path_user+'/'+usr_analysi_name+'_displacement_'+datastep+'.png'
                            
            scale_bar = 'true'
            title_win = 'NODAL DISPLACEMENTS '+datastep
            title_display = ["MAX DISPL. U: ","MID DISPL. U: ","MIN DISPL. U: ",""]
            typeData='avr'
            
            data_display = np.array([[max(abs(data_result_view[:,0]))],[sum(data_result_view[:,0])/datamesh[1]],[min(data_result_view[:,0])],[0],[usrlog.mod_opt]])
            
            save_data(path_user,usr_analysi_name,file_dir,file_save,datamesh,inci,coord,type_elm,data_display,Udef,meshDefU,scale_mesh,data_result_save,data_result_view,data_title,title_display,typeData,scale_bar,title_scale_bar,screen_on,title_win,save_screen,screen_mode_shape,ModeNumb)  
            screen_mode_shape = 'false'
            
            if np.any(graphout == 'steptracker'):
                from myfempy.lib.domain_physics import search_nodeXYZ
                
                for jj in range(len(hist_dof_list)):
                    node_coordX = float(hist_dof_list[jj,0])
                    node_coordY = float(hist_dof_list[jj,1])
                    node_coordZ = float(hist_dof_list[jj,2])
                    hist_node_list = search_nodeXYZ(node_coordX,node_coordY,node_coordZ,coord,2E-3)
                    
                    max_displ = Umag[hist_node_list-1,0]
            
        if np.any(graphout == 'stress_elem{svm}'):

            data_result_view = stress_list
            data_result_save = np.concatenate((stress_list, strain_list),axis=1)
                                              
            data_title = ['STRESS_EQV','STRAIN_EQV']
            
            file_dir=path_user+'/'+usr_analysi_name+'_stress_vm_'+datastep+'.vtk'
            file_save=path_user+'/'+usr_analysi_name+'_stress_vm_'+datastep+'.png'
            
            title_win = 'STRESS EQV VON MISES'
            scale_bar = 'true'
            title_scale_bar = 'STRESS MAIN ELEM'
            title_display = ["MAX DISPL. U: ","MAX STRESS: ","MID STRESS: ","MIN STRESS.: "]
            typeData='elm'
            data_display = np.array([[max(abs(U))[0]],[max(data_result_view[:,0])],[sum(data_result_view[:,0])/datamesh[1]],[min(data_result_view[:,0])],[usrlog.mod_opt]])
            
            save_data(path_user,usr_analysi_name,file_dir,file_save,datamesh,inci,coord,type_elm,data_display,Udef,meshDefU,scale_mesh,data_result_save,data_result_view,data_title,title_display,typeData,scale_bar,title_scale_bar,screen_on,title_win,save_screen,screen_mode_shape,ModeNumb)

            
        if np.any(graphout == 'stress_node{svm}'):
            
            strs_avr = stress_avr_plane(stress_list,datamesh,inci)
            strn_avr = stress_avr_plane(strain_list,datamesh,inci)
            
            file_dir=path_user+'/'+usr_analysi_name+'_stress_vm_avr_'+datastep+'.vtk'
            file_save=path_user+'/'+usr_analysi_name+'_stress_vm_avr_'+datastep+'.png'
   
            data_result_view = strs_avr
            data_result_save = np.concatenate((strs_avr,strn_avr),axis=1)
            
            data_title = ['STRESS_EQV','STRAIN_EQV']
            
            title_win = 'STRESS EQV VON MISES'
            scale_bar = 'true'
            title_scale_bar = 'STRESS AVR'
            title_display = ["MAX DISPL. U: ","MAX STRESS: ","MID STRESS: ","MIN STRESS.: "]
            typeData='avr' 
            data_display = np.array([[max(abs(U))[0]],[max(data_result_view[:,0])],[sum(data_result_view[:,0])/datamesh[0]],[min(data_result_view[:,0])],[usrlog.mod_opt]])
                     
            save_data(path_user,usr_analysi_name,file_dir,file_save,datamesh,inci,coord,type_elm,data_display,Udef,meshDefU,scale_mesh,data_result_save,data_result_view,data_title,title_display,typeData,scale_bar,title_scale_bar,screen_on,title_win,save_screen,screen_mode_shape,ModeNumb)
            
            if np.any(graphout == 'steptracker'):
                from myfempy.lib.domain_physics import search_nodeXYZ
                
                for jj in range(len(hist_dof_list)):
                    node_coordX = float(hist_dof_list[jj,0])
                    node_coordY = float(hist_dof_list[jj,1])
                    node_coordZ = float(hist_dof_list[jj,2])
                    hist_node_list = search_nodeXYZ(node_coordX,node_coordY,node_coordZ,coord,2E-3)
                    
                    max_stress = strs_avr[hist_node_list-1,0]
                    max_strain = strn_avr[hist_node_list-1,0]
            
            
        if np.any(graphout == 'mode'):
            
            result_list = U2
            
            ModeNumb = int(U2[0])
            Wn_rad = round(U2[1],4)
            Wn_hz = round(U2[2],4)
            
            screen_mode_shape = 'true'
            file_dir=path_user+'/'+usr_analysi_name+'_mode_shape_'+str(Wn_hz)+'_Hz_'+datastep+'.vtk'
            file_save=path_user+'/'+usr_analysi_name+'_mode_shape_'+str(Wn_hz)+'_Hz_'+datastep+'.png'
            
            data_result_save = Umag
            data_result_view = data_result_save
            data_title = ['MODE'+str(ModeNumb)+'_'+str(Wn_hz)+'_Hz']
            scale_bar = 'false'
            title_scale_bar = 'DISPL'
            title_display = ["MODE NUMB.: ","FREQ. RAD/S: ","FREQ. HZ: ","MAX DISPL.: "]
            title_win = 'MODE SHAPE'
            typeData='avr'
            
            data_display = np.array([[ModeNumb],[Wn_rad],[Wn_hz],[round(max(abs(data_result_view[:,0])),4)],[usrlog.mod_opt]])
            save_data(path_user,usr_analysi_name,file_dir,file_save,datamesh,inci,coord,type_elm,data_display,Udef,meshDefU,scale_mesh,data_result_save,data_result_view,data_title,title_display,typeData,scale_bar,title_scale_bar,screen_on,title_win,save_screen,screen_mode_shape,ModeNumb)
            screen_mode_shape = 'false'
    
        
        if np.any(graphout == 'frf'):
            from myfempy.lib.domain_physics import search_nodeXYZ
            
            for jj in range( len(hist_dof_list)):
                node_coordX = float(hist_dof_list[jj,0])
                node_coordY = float(hist_dof_list[jj,1])
                node_coordZ = float(hist_dof_list[jj,2])
                hist_node_list = search_nodeXYZ(node_coordX,node_coordY,node_coordZ,coord,2E-3)
                
                if np.any(graphout == 'frf{ux}'):
                    node_dof_hist = datamesh[0]*hist_node_list-2
                
                if np.any(graphout == 'frf{uy}'):
                    node_dof_hist = datamesh[0]*hist_node_list-1              
                                
                from myfempy.solve.view_mesh import plotterXY_resultgraph
        
                frf = U3[node_dof_hist,:]
            
                Xlabel = 'FREQUENCY [Hz]'
                Ylabel = 'FRF NODE '+str(node_dof_hist)
                plot_name = 'Response in Frequency'
                
                X = np.ravel(U2)
                Y = np.log(abs(np.ravel(frf.todense())))
                numPoints = U2.shape[0]
                
                file_save=path_user+'/'+usr_analysi_name+plot_name+'_plotter_FRF_'+datastep+str(jj)+'.png'
                plotterXY_resultgraph(file_save,save_screen,screen_on,numPoints,X,Y,Xlabel,Ylabel) 


    elif usrlog.mod_opt == "plane42":
         from myfempy.lib.plate_posproc import mesh_def_plane_2d, stress_vm_plane_q4, stress_avr_plane
         
         meshDefU,Udef,Umag = mesh_def_plane_2d(datamesh,coord,U,scale_mesh)
         stress_list, strain_list = stress_vm_plane_q4(U,datamesh,inci,coord,typeMechMat,tabmat,npp=1)
                           
         if np.any([[graphout == 'displ{ux}'],[graphout == 'displ{uy}'],[graphout == 'displ{mag}']]):
            
            result_list = np.concatenate((meshDefU,Umag),axis=1)
            
            if np.any(graphout == 'displ{ux}'):
                title_scale_bar = 'DISPL X'
                data_result_view = np.reshape(Udef[:,0],(-1,1))
            
            if np.any(graphout == 'displ{uy}'):
                title_scale_bar = 'DISPL Y'
                data_result_view = np.reshape(Udef[:,1],(-1,1))
                
            if np.any(graphout == 'displ{mag}'):
                title_scale_bar = 'DISPL MAG'
                data_result_view = Umag 
            
            data_result_save = np.concatenate((Udef,Umag),axis=1)
            data_title = ['DISPL_X','DISPL_Y','DISPL_Z','DISPL_MAG']         
            screen_mode_shape = 'true'
            file_dir=path_user+'/'+usr_analysi_name+'_displacement_'+datastep+'.vtk'
            file_save=path_user+'/'+usr_analysi_name+'_displacement_'+datastep+'.png'
                            
            scale_bar = 'true'
            title_win = 'NODAL DISPLACEMENTS '+datastep
            title_display = ["MAX DISPL. U: ","MID DISPL. U: ","MIN DISPL. U: ",""]
            typeData='avr'
           
            data_display = np.array([[max(abs(data_result_view[:,0]))],[sum(data_result_view[:,0])/datamesh[1]],[min(data_result_view[:,0])],[0],[usrlog.mod_opt]])
            save_data(path_user,usr_analysi_name,file_dir,file_save,datamesh,inci,coord,type_elm,data_display,Udef,meshDefU,scale_mesh,data_result_save,data_result_view,data_title,title_display,typeData,scale_bar,title_scale_bar,screen_on,title_win,save_screen,screen_mode_shape,ModeNumb)
            screen_mode_shape = 'false'
            
            if np.any(graphout == 'steptracker'):
               from myfempy.lib.domain_physics import search_nodeXYZ
               
               for jj in range(len(hist_dof_list)):
                   node_coordX = float(hist_dof_list[jj,0])
                   node_coordY = float(hist_dof_list[jj,1])
                   node_coordZ = float(hist_dof_list[jj,2])
                   hist_node_list = search_nodeXYZ(node_coordX,node_coordY,node_coordZ,coord,2E-3)
                   
                   max_displ = Umag[hist_node_list-1,0]
             
         if np.any([[graphout == 'stress_elem{sxx}'],[graphout == 'stress_elem{syy}'],[graphout == 'stress_elem{sxy}'],[graphout == 'stress_elem{svm}']]):
             
            file_dir=path_user+'/'+usr_analysi_name+'_stress_vm_'+datastep+'.vtk'
            file_save=path_user+'/'+usr_analysi_name+'_stress_vm_'+datastep+'.png'
             
            if np.any(graphout == 'stress_elem{sxx}'):
               title_scale_bar = 'STRESS XX'
               data_result_view = np.reshape(stress_list[:,1],(-1,1))
            
            if np.any(graphout == 'stress_elem{syy}'):
                title_scale_bar = 'STRESS YY'
                data_result_view = np.reshape(stress_list[:,2],(-1,1))
                
            if np.any(graphout == 'stress_elem{sxy}'):
                title_scale_bar = 'STRESS XY'
                data_result_view = np.reshape(stress_list[:,3],(-1,1))
                
            if np.any(graphout == 'stress_elem{svm}'):
                title_scale_bar = 'STRESS VM'
                data_result_view = np.reshape(stress_list[:,0],(-1,1))
            
                
            data_result_save = np.concatenate((stress_list,strain_list),axis=1)
            data_title = ['STRESS_EQV','STRESS_XX','STRESS_YY','STRESS_XY','STRAIN_EQV','STRAIN_XX','STRAIN_YY','STRAIN_XY']
           
            title_win = 'STRESS ELEM EQV VON MISES'
            scale_bar = 'true'
            title_scale_bar = 'STRESS MAIN ELM'
            title_display = ["MAX DISPL. U: ","MAX STRESS: ","MID STRESS: ","MIN STRESS.: "]
            typeData='elm'
            data_display = np.array([[max(abs(U))[0]],[max(data_result_view[:,0])],[sum(data_result_view[:,0])/datamesh[1]],[min(data_result_view[:,0])],[usrlog.mod_opt]])
            
            save_data(path_user,usr_analysi_name,file_dir,file_save,datamesh,inci,coord,type_elm,data_display,Udef,meshDefU,scale_mesh,data_result_save,data_result_view,data_title,title_display,typeData,scale_bar,title_scale_bar,screen_on,title_win,save_screen,screen_mode_shape,ModeNumb)
             
             
         if np.any([[graphout == 'stress_node{sxx}'],[graphout == 'stress_node{syy}'],[graphout == 'stress_node{sxy}'],[graphout == 'stress_node{svm}']]):
            strs_avr_vm = stress_avr_plane(np.reshape(stress_list[:,0],(-1,1)),datamesh,inci)
            strs_avr_xx = stress_avr_plane(np.reshape(stress_list[:,1],(-1,1)),datamesh,inci)
            strs_avr_yy = stress_avr_plane(np.reshape(stress_list[:,2],(-1,1)),datamesh,inci)
            strs_avr_xy = stress_avr_plane(np.reshape(stress_list[:,3],(-1,1)),datamesh,inci)
            
            strn_avr_vm = stress_avr_plane(np.reshape(strain_list[:,0],(-1,1)),datamesh,inci)
            strn_avr_xx = stress_avr_plane(np.reshape(strain_list[:,1],(-1,1)),datamesh,inci)
            strn_avr_yy = stress_avr_plane(np.reshape(strain_list[:,2],(-1,1)),datamesh,inci)
            strn_avr_xy = stress_avr_plane(np.reshape(strain_list[:,3],(-1,1)),datamesh,inci)
            
            max_stress = max(strs_avr_vm)
            max_strain = max(strn_avr_vm)
            
            file_dir=path_user+'/'+usr_analysi_name+'_stress_vm_avr_'+datastep+'.vtk'
            file_save=path_user+'/'+usr_analysi_name+'_stress_vm_avr_'+datastep+'.png'
            
            if np.any(graphout == 'stress_node{sxx}'):
                title_scale_bar = 'STRESS XX'
                data_result_view = strs_avr_xx
            
            if np.any(graphout == 'stress_node{syy}'):
                title_scale_bar = 'STRESS YY'
                data_result_view = strs_avr_yy
                
                
            if np.any(graphout == 'stress_node{sxy}'):
                title_scale_bar = 'STRESS XY'
                data_result_view = strs_avr_xy
                
            if np.any(graphout == 'stress_node{svm}'):
                title_scale_bar = 'STRESS VM'
                data_result_view = strs_avr_vm
            
            
            
            data_result_save = np.concatenate((strs_avr_vm, strs_avr_xx, strs_avr_yy, strs_avr_xy, strn_avr_vm, strn_avr_xx, strn_avr_yy, strn_avr_xy),axis=1)
            data_title = ['STRESS_EQV','STRESS_XX','STRESS_YY','STRESS_XY','STRAIN_EQV','STRAIN_XX','STRAIN_YY','STRAIN_XY']
            
            title_win = 'STRESS NODE EQV VON MISES'
            scale_bar = 'true'
            title_display = ["MAX DISPL. U: ","MAX STRESS: ","MID STRESS: ","MIN STRESS.: "]
            typeData='avr'
            data_display = np.array([[max(abs(U))[0]],[max(data_result_view[:,0])],[sum(data_result_view[:,0])/datamesh[0]],[min(data_result_view[:,0])],[usrlog.mod_opt]])
        
            save_data(path_user,usr_analysi_name,file_dir,file_save,datamesh,inci,coord,type_elm,data_display,Udef,meshDefU,scale_mesh,data_result_save,data_result_view,data_title,title_display,typeData,scale_bar,title_scale_bar,screen_on,title_win,save_screen,screen_mode_shape,ModeNumb)
            
            if np.any(graphout == 'steptracker'):
                from myfempy.lib.domain_physics import search_nodeXYZ
                
                for jj in range(len(hist_dof_list)):
                    node_coordX = float(hist_dof_list[jj,0])
                    node_coordY = float(hist_dof_list[jj,1])
                    node_coordZ = float(hist_dof_list[jj,2])
                    hist_node_list = search_nodeXYZ(node_coordX,node_coordY,node_coordZ,coord,2E-3)
                    
                    max_stress = strs_avr_vm[hist_node_list-1,0]
                    max_strain = strn_avr_vm[hist_node_list-1,0]
            
         if np.any(graphout == 'mode'):
             
           result_list = U2  
           
           ModeNumb = int(U2[0])
           Wn_rad = round(U2[1],4)
           Wn_hz = round(U2[2],4)
           
           screen_mode_shape = 'true'
           file_dir=path_user+'/'+usr_analysi_name+'_mode_shape_'+str(Wn_hz)+'_Hz_'+datastep+'.vtk'
           file_save=path_user+'/'+usr_analysi_name+'_mode_shape_'+str(Wn_hz)+'_Hz_'+datastep+'.png'
           
           data_result_save = Umag
           data_result_view = data_result_save
           data_title = ['MODE'+str(ModeNumb)+'_'+str(Wn_hz)+'_Hz']
           scale_bar = 'false'
           title_scale_bar = 'DISPL'
           title_display = ["MODE NUMB.: ","FREQ. RAD/S: ","FREQ. HZ: ","MAX DISPL.: "]
           title_win = 'MODE SHAPE'
           typeData='avr'
           
           data_display = np.array([[ModeNumb],[Wn_rad],[Wn_hz],[round(max(abs(data_result_view[:,0])),4)],[usrlog.mod_opt]])
           save_data(path_user,usr_analysi_name,file_dir,file_save,datamesh,inci,coord,type_elm,data_display,Udef,meshDefU,scale_mesh,data_result_save,data_result_view,data_title,title_display,typeData,scale_bar,title_scale_bar,screen_on,title_win,save_screen,screen_mode_shape,ModeNumb)
           screen_mode_shape = 'false'
       
         if np.any(graphout == 'frf'):
            from myfempy.lib.domain_physics import search_nodeXYZ
            
            for jj in range(len(hist_dof_list)):
                node_coordX = float(hist_dof_list[jj,0])
                node_coordY = float(hist_dof_list[jj,1])
                node_coordZ = float(hist_dof_list[jj,2])
                hist_node_list = search_nodeXYZ(node_coordX,node_coordY,node_coordZ,coord,2E-3)
                
                if np.any(graphout == 'frf{ux}'):
                    node_dof_hist = datamesh[0]*hist_node_list-2
                
                if np.any(graphout == 'frf{uy}'):
                    node_dof_hist = datamesh[0]*hist_node_list-1 
                
                from myfempy.solve.view_mesh import plotterXY_resultgraph
        
                frf = U3[node_dof_hist,:]
            
                Xlabel = 'FREQUENCY [Hz]'
                Ylabel = 'FRF NODE '+str(node_dof_hist)
                plot_name = 'Response in Frequency'
                
                X = np.ravel(U2)
                Y = np.log(abs(np.ravel(frf.todense())))
                numPoints = U2.shape[0]
                
                file_save=path_user+'/'+usr_analysi_name+plot_name+'_plotter_FRF_'+datastep+str(jj)+'.png'
                plotterXY_resultgraph(file_save,save_screen,screen_on,numPoints,X,Y,Xlabel,Ylabel) 

    return max_displ, max_stress, max_strain, result_list, stress_list, strain_list

def solid_output(path_user,usr_analysi_name,usrlog,graphout,datamesh,inci,coord,tabmat,type_elm,KG,U,U2,U3,screen_on,save_screen,scale_mesh,datastep,hist_dof_list):
    
    screen_mode_shape = 'false'
    
    typeMechMat = [usrlog.mat_opt,usrlog.mat_def]
    
    max_displ = np.zeros((1))
    max_stress = np.zeros((1))
    max_strain = np.zeros((1))
    
    ModeNumb = int(datastep)
    if usrlog.mod_opt == "solid83":
        from myfempy.lib.solid_posproc import  mesh_def_solid_3d, stress_vm_plane_h8, stress_avr_solid
        
        meshDefU, Udef, Umag = mesh_def_solid_3d(datamesh,coord,U,scale_mesh)

        stress_list, strain_list = stress_vm_plane_h8(U,datamesh,inci,coord,typeMechMat,tabmat,npp=1)
        
        if np.any([[graphout == 'displ{ux}'],[graphout == 'displ{uy}'],[graphout == 'displ{uz}'],[graphout == 'displ{mag}']]):
            
            result_list = np.concatenate((meshDefU,Umag),axis=1)
            
            if np.any(graphout == 'displ{ux}'):
                title_scale_bar = 'DISPL X'
                data_result_view = np.reshape(Udef[:,0],(-1,1))
            
            if np.any(graphout == 'displ{uy}'):
                title_scale_bar = 'DISPL Y'
                data_result_view = np.reshape(Udef[:,1],(-1,1))
                
            if np.any(graphout == 'displ{uz}'):
                title_scale_bar = 'DISPL Z'
                data_result_view = np.reshape(Udef[:,2],(-1,1))
                
            if np.any(graphout == 'displ{mag}'):
                title_scale_bar = 'DISPL MAG'
                data_result_view = Umag 
            
            data_result_save = np.concatenate((Udef,Umag),axis=1)
            data_title = ['DISPL_X','DISPL_Y','DISPL_Z','DISPL_MAG']
            screen_mode_shape = 'true'
            file_dir=path_user+'/'+usr_analysi_name+'_displacement_'+datastep+'.vtk'
            file_save=path_user+'/'+usr_analysi_name+'_displacement_'+datastep+'.png'
                            
            scale_bar = 'true'
            
            title_display = ["MAX DISPL. U: ","MID DISPL. U: ","MIN DISPL. U: ",""]
            title_win = 'NODAL DISPLACEMENTS '+datastep
            typeData='avr'
            
            data_display = np.array([[max(abs(data_result_view[:,0]))],[sum(data_result_view[:,0])/datamesh[1]],[min(data_result_view[:,0])],[0],[usrlog.mod_opt]])
            save_data(path_user,usr_analysi_name,file_dir,file_save,datamesh,inci,coord,type_elm,data_display,Udef,meshDefU,scale_mesh,data_result_save,data_result_view,data_title,title_display,typeData,scale_bar,title_scale_bar,screen_on,title_win,save_screen,screen_mode_shape,ModeNumb)
            screen_mode_shape = 'false'
            
            if np.any(graphout == 'steptracker'):
                from myfempy.lib.domain_physics import search_nodeXYZ
                
                for jj in range(len(hist_dof_list)):
                    node_coordX = float(hist_dof_list[jj,0])
                    node_coordY = float(hist_dof_list[jj,1])
                    node_coordZ = float(hist_dof_list[jj,2])
                    hist_node_list = search_nodeXYZ(node_coordX,node_coordY,node_coordZ,coord,2E-3)
                    
                    max_displ = Umag[hist_node_list-1,0]
                    
        if np.any([[graphout == 'stress_elem{sxx}'],[graphout == 'stress_elem{syy}'],[graphout == 'stress_elem{szz}'],[graphout == 'stress_elem{sxy}'],[graphout == 'stress_elem{syz}'],[graphout == 'stress_elem{szx}'],[graphout == 'stress_elem{svm}']]):
                        
            file_dir=path_user+'/'+usr_analysi_name+'_stress_vm_'+datastep+'.vtk'
            file_save=path_user+'/'+usr_analysi_name+'_stress_vm_'+datastep+'.png'
            
            if np.any(graphout == 'stress_elem{sxx}'):
                title_scale_bar = 'STRESS XX'
                data_result_view = np.reshape(stress_list[:,1],(-1,1))
            
            if np.any(graphout == 'stress_elem{syy}'):
                title_scale_bar = 'STRESS YY'
                data_result_view = np.reshape(stress_list[:,2],(-1,1))
                
            if np.any(graphout == 'stress_elem{szz}'):
                title_scale_bar = 'STRESS ZZ'
                data_result_view = np.reshape(stress_list[:,3],(-1,1))
                
            if np.any(graphout == 'stress_elem{sxy}'):
                title_scale_bar = 'STRESS XY'
                data_result_view = np.reshape(stress_list[:,4],(-1,1))
                
            if np.any(graphout == 'stress_elem{syz}'):
                title_scale_bar = 'STRESS YZ'
                data_result_view = np.reshape(stress_list[:,5],(-1,1))
                
            if np.any(graphout == 'stress_elem{szx}'):
                title_scale_bar = 'STRESS ZX'
                data_result_view = np.reshape(stress_list[:,6],(-1,1))
                
            if np.any(graphout == 'stress_elem{svm}'):
                title_scale_bar = 'STRESS VM'
                data_result_view = np.reshape(stress_list[:,0],(-1,1))
            
            
            data_result_save = np.concatenate((stress_list,strain_list),axis=1)
            data_title = ['STRESS_EQV','STRESS_XX','STRESS_YY','STRESS_ZZ','STRESS_XY','STRESS_YZ','STRESS_ZX','STRAIN_EQV','STRAIN_XX','STRAIN_YY','STRAIN_ZZ','STRAIN_XY','STRAIN_YZ','STRAIN_ZX']
            
            title_win = 'STRESS ELEM'
            scale_bar = 'true'
            title_display = ["MAX DISPL. U: ","MAX STRESS: ","MID STRESS: ","MIN STRESS.: "]
            typeData='elm'
            data_display = np.array([[max(abs(U))[0]],[max(data_result_view[:,0])],[sum(data_result_view[:,0])/datamesh[1]],[min(data_result_view[:,0])],[usrlog.mod_opt]])
            
            if np.any(graphout == 'steptracker'):
                print("output erro: steptracker is don't applied in stress_elem function")
            
            save_data(path_user,usr_analysi_name,file_dir,file_save,datamesh,inci,coord,type_elm,data_display,Udef,meshDefU,scale_mesh,data_result_save,data_result_view,data_title,title_display,typeData,scale_bar,title_scale_bar,screen_on,title_win,save_screen,screen_mode_shape,ModeNumb)

            
        if np.any([[graphout == 'stress_node{sxx}'],[graphout == 'stress_node{syy}'],[graphout == 'stress_node{szz}'],[graphout == 'stress_node{sxy}'],[graphout == 'stress_node{syz}'],[graphout == 'stress_node{szx}'],[graphout == 'stress_node{svm}']]):
            
            strs_avr_vm = stress_avr_solid(np.reshape(stress_list[:,0],(-1,1)),datamesh,inci)
            strs_avr_xx = stress_avr_solid(np.reshape(stress_list[:,1],(-1,1)),datamesh,inci)
            strs_avr_yy = stress_avr_solid(np.reshape(stress_list[:,2],(-1,1)),datamesh,inci)
            strs_avr_zz = stress_avr_solid(np.reshape(stress_list[:,3],(-1,1)),datamesh,inci)
            strs_avr_xy = stress_avr_solid(np.reshape(stress_list[:,4],(-1,1)),datamesh,inci)
            strs_avr_yz = stress_avr_solid(np.reshape(stress_list[:,5],(-1,1)),datamesh,inci)
            strs_avr_zx = stress_avr_solid(np.reshape(stress_list[:,6],(-1,1)),datamesh,inci)
            
            strn_avr_vm = stress_avr_solid(np.reshape(strain_list[:,0],(-1,1)),datamesh,inci)
            strn_avr_xx = stress_avr_solid(np.reshape(strain_list[:,1],(-1,1)),datamesh,inci)
            strn_avr_yy = stress_avr_solid(np.reshape(strain_list[:,2],(-1,1)),datamesh,inci)
            strn_avr_zz = stress_avr_solid(np.reshape(strain_list[:,3],(-1,1)),datamesh,inci)
            strn_avr_xy = stress_avr_solid(np.reshape(strain_list[:,4],(-1,1)),datamesh,inci)
            strn_avr_yz = stress_avr_solid(np.reshape(strain_list[:,5],(-1,1)),datamesh,inci)
            strn_avr_zx = stress_avr_solid(np.reshape(strain_list[:,6],(-1,1)),datamesh,inci)
                                        
            file_dir=path_user+'/'+usr_analysi_name+'_stress_vm_avr_'+datastep+'.vtk'
            file_save=path_user+'/'+usr_analysi_name+'_stress_vm_avr_'+datastep+'.png'
            
            if np.any(graphout == 'stress_node{sxx}'):
                title_scale_bar = 'STRESS XX'
                data_result_view = strs_avr_xx
            
            if np.any(graphout == 'stress_node{syy}'):
                title_scale_bar = 'STRESS YY'
                data_result_view = strs_avr_yy
                
            if np.any(graphout == 'stress_node{szz}'):
                title_scale_bar = 'STRESS ZZ'
                data_result_view = strs_avr_zz
                
            if np.any(graphout == 'stress_node{sxy}'):
                title_scale_bar = 'STRESS XY'
                data_result_view = strs_avr_xy
                
            if np.any(graphout == 'stress_node{syz}'):
                title_scale_bar = 'STRESS YZ'
                data_result_view = strs_avr_yz 
                
            if np.any(graphout == 'stress_node{szx}'):
                title_scale_bar = 'STRESS ZX'
                data_result_view = strs_avr_zx 
                
            if np.any(graphout == 'stress_node{svm}'):
                title_scale_bar = 'STRESS VM'
                data_result_view = strs_avr_vm
            
            
            data_result_save = np.concatenate((strs_avr_vm, strs_avr_xx, strs_avr_yy, strs_avr_zz, strs_avr_xy, strs_avr_yz, strs_avr_zx, strn_avr_vm, strn_avr_xx, strn_avr_yy, strn_avr_zz, strn_avr_xy, strn_avr_yz, strn_avr_zx),axis=1)
            data_title = ['STRESS_EQV','STRESS_XX','STRESS_YY','STRESS_ZZ','STRESS_XY','STRESS_YZ','STRESS_ZX','STRAIN_EQV','STRAIN_XX','STRAIN_YY','STRAIN_ZZ','STRAIN_XY','STRAIN_YZ','STRAIN_ZX']
            
            title_win = 'STRESS NODE'
            scale_bar = 'true'
            title_display = ["MAX DISPL. U: ","MAX STRESS: ","MID STRESS: ","MIN STRESS.: "]
            typeData='avr'
            data_display = np.array([[max(abs(U))[0]],[max(data_result_view[:,0])],[sum(data_result_view[:,0])/datamesh[0]],[min(data_result_view[:,0])],[usrlog.mod_opt]])
            
            if np.any(graphout == 'steptracker'):
                from myfempy.lib.domain_physics import search_nodeXYZ
                
                for jj in range(len(hist_dof_list)):
                    node_coordX = float(hist_dof_list[jj,0])
                    node_coordY = float(hist_dof_list[jj,1])
                    node_coordZ = float(hist_dof_list[jj,2])
                    hist_node_list = search_nodeXYZ(node_coordX,node_coordY,node_coordZ,coord,2E-3)
                    
                    max_stress = strs_avr_vm[hist_node_list-1,0]
                    max_strain = strn_avr_vm[hist_node_list-1,0]
            
            save_data(path_user,usr_analysi_name,file_dir,file_save,datamesh,inci,coord,type_elm,data_display,Udef,meshDefU,scale_mesh,data_result_save,data_result_view,data_title,title_display,typeData,scale_bar,title_scale_bar,screen_on,title_win,save_screen,screen_mode_shape,ModeNumb)
            
        if np.any(graphout == 'mode'):
            
           result_list = U2 #np.concatenate((meshDefU,Umag),axis=1)
            
           ModeNumb = int(U2[0])
           Wn_rad = round(U2[1],4)
           Wn_hz = round(U2[2],4)
           
           screen_mode_shape = 'true'
           file_dir=path_user+'/'+usr_analysi_name+'_mode_shape_'+str(Wn_hz)+'_Hz_'+datastep+'.vtk'
           file_save=path_user+'/'+usr_analysi_name+'_mode_shape_'+str(Wn_hz)+'_Hz_'+datastep+'.png'
           
           data_result_save = Umag
           data_result_view = data_result_save
           data_title = ['MODE'+str(ModeNumb)+'_'+str(Wn_hz)+'_Hz']
           scale_bar = 'false'
           title_scale_bar = 'DISPL'
           title_display = ["MODE NUMB.: ","FREQ. RAD/S: ","FREQ. HZ: ","MAX DISPL.: "]
           title_win = 'MODE SHAPE'
           typeData='avr'
           
           data_display = np.array([[ModeNumb],[Wn_rad],[Wn_hz],[round(max(abs(data_result_view[:,0])),4)],[usrlog.mod_opt]])
           save_data(path_user,usr_analysi_name,file_dir,file_save,datamesh,inci,coord,type_elm,data_display,Udef,meshDefU,scale_mesh,data_result_save,data_result_view,data_title,title_display,typeData,scale_bar,title_scale_bar,screen_on,title_win,save_screen,screen_mode_shape,ModeNumb)
           screen_mode_shape = 'false'
                                     
        
        if np.any([[graphout == 'frf{ux}'],[graphout == 'frf{uy}'],[graphout == 'frf{uz}']]):
            from myfempy.lib.domain_physics import search_nodeXYZ
            
            result_list = np.zeros((np.shape(U3)[0],len(hist_dof_list)))
            
            for jj in range(len(hist_dof_list)):
                node_coordX = float(hist_dof_list[jj,0])
                node_coordY = float(hist_dof_list[jj,1])
                node_coordZ = float(hist_dof_list[jj,2])
                hist_node_list = search_nodeXYZ(node_coordX,node_coordY,node_coordZ,coord,2E-3)
                
                if np.any(graphout == 'frf{ux}'):
                    node_dof_hist = datamesh[0]*hist_node_list-3
                
                if np.any(graphout == 'frf{uy}'):
                    node_dof_hist = datamesh[0]*hist_node_list-2
                    
                if np.any(graphout == 'frf{uz}'):
                    node_dof_hist = datamesh[0]*hist_node_list-1

                from myfempy.solve.view_mesh import plotterXY_resultgraph
        
                frf = U3[node_dof_hist,:]
            
                Xlabel = 'FREQUENCY [Hz]'
                Ylabel = 'FRF NODE '+str(node_dof_hist)
                plot_name = 'RESPONSE IN FREQUENCY'
                
                X = np.ravel(U2)
                Y = np.log(abs(np.ravel(frf.todense())))
                numPoints = U2.shape[0]
                
                file_save=path_user+'/'+usr_analysi_name+plot_name+'_plotter_FRF_'+datastep+str(jj)+'.png'
                plotterXY_resultgraph(file_save,save_screen,screen_on,numPoints,X,Y,Xlabel,Ylabel) 

    return max_displ, max_stress, max_strain, result_list


def myfempy_plotstep(x,y,xlabel,ylabel,fignumb):
    
    # plt.gcf().set_size_inches(8, 6)               
    plt.figure(fignumb)
    plt.plot(x,y,'-or')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid(True)
    plt.show()
    plt.pause(0.1)
    

def graphic_prosproc(path_user,usr_analysi_name,usrlog,graphout,fileout,time_list,datamesh,inci,coord,tabgeo,CGcoord,tabmat,type_elm,KG,U,U2,U3,screen_on,save_screen):
    scale_mesh = 1
    hist_dof = np.zeros((1,3))
    hist_dof_list = np.zeros((1,3))
    for jj in range(1,len(graphout)):
        straux = graphout[jj]
        if straux.count('scale{') == 1:
            idx0 = straux.find('{')
            idx1 = straux.find('}')
            scale_mesh = float(straux[idx0+1:idx1])
        
        elif straux.count('point{') == 1:
            hist_list = re.findall(r"[-+]?\d*\.\d+|\d+", straux)
            hist_dof[0,0] = float(hist_list[0])
            hist_dof[0,1] = float(hist_list[1])
            hist_dof[0,2] = float(hist_list[2])
            hist_dof_list = np.append(hist_dof_list,hist_dof,axis=0)
    
    hist_dof_list = hist_dof_list[1::][::]  
  
    max_displ_list = np.zeros((1))
    max_stress_list = np.zeros((1))
    max_strain_list = np.zeros((1))
    displ_list = np.zeros((1))
    natfreq_list = np.zeros((1))
    stress_list = np.zeros((1))
    strain_list = np.zeros((1))
    max_displ = np.zeros((1))
    max_stress = np.zeros((1))
    max_strain = np.zeros((1))
    
    # if usrlog.solver_opt == 'step':
        
        
    # elif usrlog.solver_opt == 'mode':
    #     print()
        
    # elif usrlog.solver_opt == 'fqrp':
    #     print()
        
    # elif usrlog.solver_opt == 'time':
    #     print()
    
    step = np.zeros((1))
    
    plt.close('all')
    plt.style.use('classic')
    for ii in range(0,U.shape[1]):
        datastep = str(ii)
        Ustep = U[:,ii]
        U2step = U2[ii,:]
            
        if usrlog.mod_typ == "beam":
            max_displ, max_stress, max_strain, result_list, stress_list = beam_output(path_user,usr_analysi_name,usrlog,graphout,datamesh,inci,coord,tabgeo,CGcoord,tabmat,type_elm,KG,Ustep,U2step,U3,screen_on,save_screen,scale_mesh,datastep,hist_dof_list)
                        
        if usrlog.mod_typ == "plate":
            max_displ, max_stress, max_strain, displ_list, stress_list, strain_list, natfreq_list = plate_output(path_user,usr_analysi_name,usrlog,graphout,datamesh,inci,coord,tabmat,type_elm,KG,Ustep,U2step,U3,screen_on,save_screen,scale_mesh,datastep,hist_dof_list)

        if usrlog.mod_typ == "solid":
            max_displ, max_stress, max_strain, result_list = solid_output(path_user,usr_analysi_name,usrlog,graphout,datamesh,inci,coord,tabmat,type_elm,KG,Ustep,U2step,U3,screen_on,save_screen,scale_mesh,datastep,hist_dof_list)
        
        max_displ_list = np.append(max_displ_list,max_displ,axis=0)
        max_stress_list = np.append(max_stress_list,max_stress,axis=0)
        max_strain_list = np.append(max_strain_list,max_strain,axis=0)
        step = np.append(step,np.array([ii+1]),axis=0)
        
        if np.any(graphout == 'steptracker'):
            if np.any(graphout == 'displ'):
                myfempy_plotstep(step,max_displ_list,'Steps','Displacement',1)
                
            if np.any(graphout == 'stress_avr'):
                myfempy_plotstep(max_strain_list,max_stress_list,'Strain','Stress',2)
    
    
    # SAVE OUTPUT FILE
    from myfempy.setup.myfempy_preproc import white_output_logofile
    output_file = 'myfempy_output_logfile.txt'
    white_output_logofile(path_user,usr_analysi_name,output_file,fileout,time_list,result_list,datamesh,inci,coord,tabgeo,tabmat)
    