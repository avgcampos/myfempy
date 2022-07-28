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

import numpy as np
import os
import scipy.sparse as sp
import matplotlib.pyplot as plt
from myfempy.plots.postplot import postproc_plot
from myfempy.felib.felemset import get_elemset
from myfempy.postprc.postset import get_stress, get_displ
from myfempy.io.iovtk import convert_to_vtk
# from myfempy.felib.geometry.getnode import search_nodexyz
from myfempy.plots.plotxy import tracker_plot
from myfempy.io.iovtk import convert_to_vtk
from myfempy.tools.tools import print_console

class PostComputer:
    
    def __init__(self, modelinfo):
        
        self.modelinfo = modelinfo
        self.dofe = modelinfo['nodecon'][0]*modelinfo['nodedof'][0]
        self.nodecon = modelinfo['nodecon'][0]
        self.fulldof = modelinfo["nodedof"][0]*len(modelinfo["coord"])
        self.nodedof = modelinfo['nodedof'][0]
        self.nelem = len(modelinfo["inci"])
        self.nnode =len(modelinfo["coord"])
        self.inci = modelinfo['inci']
        self.coord = modelinfo['coord']
        self.tabmat = modelinfo['tabmat']
        self.tabgeo = modelinfo['tabgeo']
        self.ntensor = modelinfo['ntensor'][0]
    
    
    def displ(self, U, scale):
        
        meshdefU, result = get_displ(self.modelinfo, U, scale)
              
        return meshdefU, result
    
    
    def stress(self, U):
            
        # nelem = len(modelinfo["inci"])
        stress_list = np.zeros((self.nelem,self.modelinfo['ntensor'][0]+1),dtype=float)
        strain_list = np.zeros((self.nelem,self.modelinfo['ntensor'][0]+1),dtype=float)
        for ee in range(self.nelem):
            
            tensor = get_stress(self.modelinfo, U, ee)
            
            epsilon, strain, tistrn = tensor.strain()
            stress, tistrs = tensor.stress(epsilon)
                      
            stress_list[ee,:] = stress
            strain_list[ee,:] = strain
                            
        result = np.concatenate((stress_list, strain_list), axis=1)
        title = np.concatenate((tistrs, tistrn), axis=0)
        return result, title
    
    
    
    def intforces(self, U):
        # nelem = len(modelinfo["inci"])
            
        element = get_elemset(int(self.modelinfo['elemid'][0]))
        setelement = element(self.modelinfo)
        lines = [[None]*2]
        if len(self.modelinfo['regions']) == 0:
            for ll in range(len(self.modelinfo['inci'])):
                lines.append([self.modelinfo['inci'][ll,0], self.modelinfo['inci'][ll,[4,5]].tolist()])
            lines = lines[1::][::]
        else:
            lines = self.modelinfo['regions'][1][1]
            
        
        ifb, title = setelement.intforces(U, lines)
            
        return ifb, title
            
    
    def elem2nodes_filter(self):
    
        # inci = modelinfo['inci']
        # nnode = modelinfo['coord'].shape[0]
        # nelem = modelinfo['inci'].shape[0]
        # # fulldof = modelinfo["nodedof"][0]*len(modelinfo["coord"])
        # dofe = modelinfo['nodecon'][0]*modelinfo['nodedof'][0]
        ith = np.zeros((self.nelem*(self.dofe*self.dofe)),dtype=int)
        jth = np.zeros((self.nelem*(self.dofe*self.dofe)),dtype=int)
        val = np.zeros((self.nelem*(self.dofe*self.dofe)),dtype=int)
        
        q0 = 0
        for i in range(self.nnode):
            
            elmlist = self.modelinfo['inci'][(np.asarray(np.where(self.modelinfo['inci'][:,4:]==i+1)))[0][:],0]
            q1 = elmlist.size
            ith[q0:q1+q0] = i
            jth[q0:q1+q0] = elmlist-1
            val[q0:q1+q0] = elmlist
            q0=q1+q0 
            
        S = sp.csc_matrix((val, (ith, jth)), shape=(self.nnode, self.nelem))
        
        return S
    
    
    #-----------------------------------------------------------------------------#
    def results_average(self, results_elm):
        
        S =  PostComputer.elem2nodes_filter(self)
        
        results_avr = np.zeros((self.nnode), dtype=float)
        for mm in range(self.nnode):
            results_avr[mm] = np.mean(results_elm[(S[mm,:].nonzero())[1]])
        return results_avr
    
    
    # #-----------------------------------------------------------------------------#
    # def gen_animation(path_user,usr_analysi_name,type_elm,inci,coord,datamesh,Udef,data_result,data_title,typeData,time_out,scale_mesh,ModeNumb):
    #     from myfempy.setup.myfempy_preproc import create_user_dir
    #     path_user = create_user_dir(path_user,'/anime_file')
    #     # scale_mesh_frac0 = np.linspace(-time_out,time_out,time_out)
    #     scale_mesh_frac0 = np.concatenate((np.linspace(-time_out,time_out,int(time_out/2)),np.linspace(time_out,-time_out,int(time_out/2))),axis=0)
    #     for itime in range(0,time_out):
    #         scale_mesh_frac1 = 0.01*scale_mesh*scale_mesh_frac0[itime]
    #         file_dir=path_user+'/'+usr_analysi_name+'_mode_shape_anime_'+str(ModeNumb)+'_'+(str(itime))+'.vtk'
    #         meshDefU_anime = np.concatenate((coord[:,0][:, np.newaxis],np.add(coord[:,1:],scale_mesh_frac1*Udef)),axis=1)
    #         export_mesh(file_dir,type_elm,datamesh[1],datamesh[2],meshDefU_anime,inci,data_result,data_title,typeData)
            
                   
    
    
    def save_vtk(self, postporc_result, postprocset):
        
        path = os.getcwd()
        
        plotdata = dict()
        
        plotdata['inci'] = self.modelinfo['inci']
        plotdata['nodecon'] = self.modelinfo['nodecon']
        plotdata['elemid']  = self.modelinfo['elemid']
        plotdata['coord'] = self.modelinfo['coord']
        
        plotdata['displ_POINT_DATA_val'] = []
        plotdata['displ_POINT_DATA_name'] = []
        plotdata['displ_POINT_DATA_title'] = []
        
        plotdata['stress_CELL_DATA_val'] = []
        plotdata['stress_CELL_DATA_name'] = []
        plotdata['stress_CELL_DATA_title'] = []
        
        plotdata['stress_POINT_DATA_val'] = []
        plotdata['stress_POINT_DATA_name'] = []
        plotdata['stress_POINT_DATA_title'] = []        
        
        plotdata['modes_POINT_DATA'] = []
        
        if 'displ' in postprocset['COMPUTER'].keys():
            for st in range(len(postporc_result['displ'])):
    
                plotdata['filename'] = (path+'/'+postprocset["PLOTSET"]['filename']+'_results_step-'+str(st+1)) #(postprocset["PLOTSET"]['filename']+'_'+str(st+1)+'_'+postporc_result['solution'][pp]['name'])
                           
                # plotdata['defor_POINT_DATA_val'] = postporc_result['displ'][st]['val'][:,0]
                plotdata['displ_POINT_DATA_val'] = postporc_result['displ'][st]['val'][:,1:]
                # plotdata['displ_POINT_DATA_name'] = postporc_result['displ'][st]['name']
                plotdata['displ_POINT_DATA_title'] = postporc_result['displ'][st]['title']
                
                if postprocset['COMPUTER']['stress'] == True:
                    plotdata['stress_CELL_DATA_val'] = postporc_result['stress_elm'][st]['val']
                    # plotdata['stress_CELL_DATA_name'] = postporc_result['stress_elm'][st]['name']
                    plotdata['stress_CELL_DATA_title'] = postporc_result['stress_elm'][st]['title']
                    
                if postprocset['COMPUTER']['average'] == True:
                    plotdata['stress_POINT_DATA_val'] = postporc_result['stress_avr'][st]['val']
                    # plotdata['stress_POINT_DATA_name'] = postporc_result['stress_avr'][st]['name']
                    plotdata['stress_POINT_DATA_title'] = postporc_result['stress_avr'][st]['title']
                    
                convert_to_vtk(plotdata)
                    
        
        if 'eig' in postprocset['COMPUTER'].keys():
            # plotdata['defor_POINT_DATA_val'] = postporc_result['displ'][st]['val'][:,0]
            
            plotdata['filename'] = (path+'/'+postprocset["PLOTSET"]['filename']+'_results_modes')
            
            plotdata['modes_POINT_DATA'] = postporc_result['modes']#[st]['val'][:,1:]
            # plotdata['modes_POINT_DATA_name'] = postporc_result['modes'][st]['name']
            # plotdata['modes_POINT_DATA_title'] = postporc_result['modes'][st]['title']
            
            convert_to_vtk(plotdata)
            
            
    

    #-----------------------------------------------------------------------------#
    # @profile
    def main(self, postprocset):
        
        print_console('post')

        postporc_result = dict()
        postporc_result['displ'] = []
        postporc_result['stress_avr'] = []
        postporc_result['stress_elm'] = []
        postporc_result['balance'] = []
        postporc_result['modes'] = []
        # postporc_result['set'] = {'avr': postprocset["COMPUTER"]["average"]}
        # scale = postprocset["SCALE"]
        
        scale = (postprocset["SCALE"]/100)*(1/(max(abs(postprocset["SOLUTION"][:,0]))))
        
        result_disp = np.zeros((postprocset["SOLUTION"].shape[1],self.modelinfo['coord'].shape[0],4))
        meshdefU = np.zeros((postprocset["SOLUTION"].shape[1],self.modelinfo['coord'].shape[0],self.modelinfo['coord'].shape[1]))
        for ns in range(postprocset["SOLUTION"].shape[1]):
            meshdefU[ns,:,:], result_disp[ns,:,:] = PostComputer.displ(self, postprocset["SOLUTION"][:,ns], scale)
        
        postporc_result['defcoord'] = meshdefU
        
        if 'displ' in postprocset['COMPUTER'].keys():
            if postprocset['COMPUTER']['displ'] == True:
                for st in range(postprocset["SOLUTION"].shape[1]):
                    # postporc_result['displ'] = {'val':result_disp, 'title':title}
                    postporc_result['displ'].append({'val':result_disp[st][:][:], 'title':'DISPLACEMENT', 'avr': True})
            else:
                pass
                
        if 'stress' in postprocset['COMPUTER'].keys():
            if postprocset['COMPUTER']['stress'] == True:
                result_stress = np.zeros((postprocset["SOLUTION"].shape[1],self.modelinfo['inci'].shape[0],2*self.modelinfo["ntensor"][0]+2))
                # result_avr = np.zeros((postprocset["SOLUTION"].shape[1],modelinfo['coord'].shape[0],modelinfo["ntensor"][0]))
                for ns in range(postprocset["SOLUTION"].shape[1]):
                    result_stress[ns,:,:], title = PostComputer.stress(self, postprocset["SOLUTION"][:,ns])
                    
                # if postprocset["COMPUTER"]["average"] == True:
                results_avr = np.zeros((postprocset["SOLUTION"].shape[1],self.modelinfo['coord'].shape[0],2*self.modelinfo["ntensor"][0]+2))
                for ns in range(postprocset["SOLUTION"].shape[1]):
                    for nt in range(2*self.modelinfo["ntensor"][0]+2):
                        results_avr[ns,:,nt] = PostComputer.results_average(self, result_stress[ns,:,nt])
                result_stress_avr = results_avr   
                      
                for st in range(postprocset["SOLUTION"].shape[1]): 
                    postporc_result['stress_elm'].append({'val':result_stress[st][:][:], 'title':title, 'avr': False})
                    postporc_result['stress_avr'].append({'val':result_stress_avr[st][:][:], 'title':title, 'avr': True})
        
        
        if 'balance' in postprocset['COMPUTER'].keys():
            if postprocset['COMPUTER']['balance'] == True:
                # result_ifb = dict()
                # result_ifb['intforces'] = [] #np.zeros((postprocset["SOLUTION"].shape[1],modelinfo['coord'].shape[0],modelinfo['nodedof'][0]+1))
                for ns in range(postprocset["SOLUTION"].shape[1]):
                    ifb, title = PostComputer.intforces(self, postprocset["SOLUTION"][:,ns])
                    postporc_result['balance'].append({'val': [ifb['le'], ifb['val']], 'title': title})
            else:
                pass
            
        if 'eig' in postprocset['COMPUTER'].keys():
            if postprocset['COMPUTER']['eig'] == True:
                for st in range(postprocset["SOLUTION"].shape[1]):
                    # postporc_result['displ'] = {'val':result_disp, 'title':title}
                    postporc_result['modes'].append({'val':result_disp[st][:][:], 'title':('MODE_'+str(st+1)+'-'+'FREQ_'+str(round(postprocset["FREQ"][st][2],2))+'Hz'), 'avr': True})
            else:
                pass
                        
        PostComputer.save_vtk(self, postporc_result, postprocset) # save in vtk file
               
        # data_result_save = {'disp':displ_result_save,'stress':stress_result_save}
        return postporc_result

    
    
