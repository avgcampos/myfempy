#!/usr/bin/env python
"""
Post-Process Calculator
"""
__author__ = "Antonio Vinicius Garcia Campos"
__copyright__ = "Copyright @ 2022, Antonio Vinicius Garcia Campos"
__credits__ = ["Antonio Vinicius Garcia Campos", "3D EasyCAE"]
__license__ = "GPL"
__status__ = "Development"
import numpy as np
import os
import scipy.sparse as sp
from myfempy.felib.felemset import get_elemset
from myfempy.postprc.postset import get_stress, get_displ
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
        self.nnode = len(modelinfo["coord"])
        self.inci = modelinfo['inci']
        self.coord = modelinfo['coord']
        self.tabmat = modelinfo['tabmat']
        self.tabgeo = modelinfo['tabgeo']
        self.ntensor = modelinfo['ntensor'][0]

    def displ(self, U):
        result = get_displ(self.modelinfo, U)
        return result

    def stress(self, U):
        stress_list = np.zeros(
            (self.nelem, self.modelinfo['ntensor'][0]+1), dtype=float)
        strain_list = np.zeros(
            (self.nelem, self.modelinfo['ntensor'][0]+1), dtype=float)
        for ee in range(self.nelem):
            tensor = get_stress(self.modelinfo, U, ee)
            epsilon, strain, tistrn = tensor.strain()
            stress, tistrs = tensor.stress(epsilon)
            stress_list[ee, :] = stress
            strain_list[ee, :] = strain
        result = np.concatenate((stress_list, strain_list), axis=1)
        title = np.concatenate((tistrs, tistrn), axis=0)
        return result, title

    def intforces(self, U):
        element = get_elemset(int(self.modelinfo['elemid'][0]))
        setelement = element(self.modelinfo)
        lines = [[None]*2]
        if len(self.modelinfo['regions']) == 0:
            for ll in range(len(self.modelinfo['inci'])):
                lines.append([self.modelinfo['inci'][ll, 0],
                             self.modelinfo['inci'][ll, [4, 5]].tolist()])
            lines = lines[1::][::]
        else:
            lines = self.modelinfo['regions'][1][1]
        ifb, title = setelement.intforces(U, lines)
        return ifb, title

    def elem2nodes_filter(self):
        ith = np.zeros((self.nelem*(self.dofe*self.dofe)), dtype=int)
        jth = np.zeros((self.nelem*(self.dofe*self.dofe)), dtype=int)
        val = np.zeros((self.nelem*(self.dofe*self.dofe)), dtype=int)
        q0 = 0
        for i in range(self.nnode):
            elmlist = self.modelinfo['inci'][(np.asarray(
                np.where(self.modelinfo['inci'][:, 4:] == i+1)))[0][:], 0]
            q1 = elmlist.size
            ith[q0:q1+q0] = i
            jth[q0:q1+q0] = elmlist-1
            val[q0:q1+q0] = elmlist
            q0 = q1+q0
        S = sp.csc_matrix((val, (ith, jth)), shape=(self.nnode, self.nelem))
        return S

    def results_average(self, results_elm):
        S = PostComputer.elem2nodes_filter(self)
        results_avr = np.zeros((self.nnode), dtype=float)
        for mm in range(self.nnode):
            results_avr[mm] = np.mean(results_elm[(S[mm, :].nonzero())[1]])
        return results_avr

    def save_vtk(self, postporc_result, postprocset):
        path = os.getcwd()
        plotdata = dict()
        plotdata['inci'] = self.modelinfo['inci']
        plotdata['nodecon'] = self.modelinfo['nodecon']
        plotdata['elemid'] = self.modelinfo['elemid']
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
                plotdata['filename'] = (
                    path+'/'+postprocset["PLOTSET"]['filename']+'_results_step-'+str(st+1))
                plotdata['displ_POINT_DATA_val'] = postporc_result['displ'][st]['val'][:, 1:]
                plotdata['displ_POINT_DATA_title'] = postporc_result['displ'][st]['title']
                if 'stress' in postprocset['COMPUTER'].keys():
                    if postprocset['COMPUTER']['stress'] == True:
                        plotdata['stress_CELL_DATA_val'] = postporc_result['stress_elm'][st]['val']
                        plotdata['stress_CELL_DATA_title'] = postporc_result['stress_elm'][st]['title']
                if 'average' in postprocset['COMPUTER'].keys():
                    if postprocset['COMPUTER']['average'] == True:
                        plotdata['stress_POINT_DATA_val'] = postporc_result['stress_avr'][st]['val']
                        plotdata['stress_POINT_DATA_title'] = postporc_result['stress_avr'][st]['title']
                convert_to_vtk(plotdata)
        if 'eig' in postprocset['COMPUTER'].keys():
            plotdata['filename'] = (
                path+'/'+postprocset["PLOTSET"]['filename']+'_results_modes')
            plotdata['modes_POINT_DATA'] = postporc_result['modes']
            convert_to_vtk(plotdata)

    def main(self, postprocset):
        print_console('post')
        postporc_result = dict()
        postporc_result['displ'] = []
        postporc_result['stress_avr'] = []
        postporc_result['stress_elm'] = []
        postporc_result['balance'] = []
        postporc_result['modes'] = []
        postporc_result['frf'] = []
        result_disp = np.zeros(
            (postprocset["SOLUTION"].shape[1], self.modelinfo['coord'].shape[0], 4))
        for ns in range(postprocset["SOLUTION"].shape[1]):
            result_disp[ns, :, :] = PostComputer.displ(
                self, postprocset["SOLUTION"][:, ns])
        if 'displ' in postprocset['COMPUTER'].keys():
            if postprocset['COMPUTER']['displ'] == True:
                for st in range(postprocset["SOLUTION"].shape[1]):
                    # postporc_result['displ'] = {'val':result_disp, 'title':title}
                    postporc_result['displ'].append(
                        {'val': result_disp[st][:][:], 'title': 'DISPLACEMENT', 'avr': True})
            else:
                pass
        if 'stress' in postprocset['COMPUTER'].keys():
            if postprocset['COMPUTER']['stress'] == True:
                result_stress = np.zeros(
                    (postprocset["SOLUTION"].shape[1], self.modelinfo['inci'].shape[0], 2*self.modelinfo["ntensor"][0]+2))
                for ns in range(postprocset["SOLUTION"].shape[1]):
                    result_stress[ns, :, :], title = PostComputer.stress(
                        self, postprocset["SOLUTION"][:, ns])
                results_avr = np.zeros(
                    (postprocset["SOLUTION"].shape[1], self.modelinfo['coord'].shape[0], 2*self.modelinfo["ntensor"][0]+2))
                for ns in range(postprocset["SOLUTION"].shape[1]):
                    for nt in range(2*self.modelinfo["ntensor"][0]+2):
                        results_avr[ns, :, nt] = PostComputer.results_average(
                            self, result_stress[ns, :, nt])
                result_stress_avr = results_avr
                for st in range(postprocset["SOLUTION"].shape[1]):
                    postporc_result['stress_elm'].append(
                        {'val': result_stress[st][:][:], 'title': title, 'avr': False})
                    postporc_result['stress_avr'].append(
                        {'val': result_stress_avr[st][:][:], 'title': title, 'avr': True})
        if 'balance' in postprocset['COMPUTER'].keys():
            if postprocset['COMPUTER']['balance'] == True:
                for ns in range(postprocset["SOLUTION"].shape[1]):
                    ifb, title = PostComputer.intforces(
                        self, postprocset["SOLUTION"][:, ns])
                    postporc_result['balance'].append(
                        {'val': [ifb['le'], ifb['val']], 'title': title})
            else:
                pass
        if 'eig' in postprocset['COMPUTER'].keys():
            if postprocset['COMPUTER']['eig'] == True:
                for st in range(postprocset["SOLUTION"].shape[1]):
                    # postporc_result['displ'] = {'val':result_disp, 'title':title}
                    postporc_result['modes'].append({'val': result_disp[st][:][:], 'title': (
                        'MODE_'+str(st+1)+'-'+'FREQ_'+str(round(postprocset["FREQ"][st][2], 2))+'Hz'), 'avr': True})
            else:
                pass
        if 'frf' in postprocset['COMPUTER'].keys():
            if postprocset['COMPUTER']['frf'] == True:
                postporc_result['frf'].append(
                    {'val': postprocset["SOLUTION"], 'freqlog': postprocset["FREQ"]})
            else:
                pass
        PostComputer.save_vtk(self, postporc_result,
                              postprocset)  # save in vtk file
        return postporc_result
