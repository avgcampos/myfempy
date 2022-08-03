#!/usr/bin/env python
"""
Plotter Post Process
"""
__author__ = "Antonio Vinicius Garcia Campos"
__copyright__ = "Copyright @ 2022, Antonio Vinicius Garcia Campos"
__credits__ = ["Antonio Vinicius Garcia Campos", "3D EasyCAE"]
__license__ = "GPL"
__status__ = "Development"
import numpy as np
from myfempy.felib.physics.getnode import search_nodexyz
from myfempy.plots.plotxy import plot_forces, tracker_plot, frf_plot
from myfempy.plots.plotmesh import post_show_mesh


def postproc_plot(postprocset, postporc_result, modelinfo):
    plotset = dict()
    if "TRACKER" in postprocset.keys():
        if postprocset["TRACKER"]['show'] == True:
            # for pp in range(len(postporc_result['solution'])):
            for st in range(len(plotset[postprocset["TRACKER"]['result2plot']])):
                plotset['step'] = st+1
                plotset['val_list'] = plotset[postprocset["TRACKER"]
                                              ['result2plot']][st]
                plotset['fignumb'] = 99
                plotset['rstl'] = [0, modelinfo["ntensor"][0]+1]
                tracker_plot(postprocset, plotset, modelinfo['coord'])
    if "PLOTSET" in postprocset.keys():
        if postprocset["PLOTSET"]['show'] == True:
            if 'step' in postprocset["PLOTSET"].keys():
                step = int(postprocset["PLOTSET"]['step'])
            else:
                step = int(1)
            plotset['text_plot'] = ('DISPL'+' step: '+str(step))
            plotset['step'] = step-1
            file2plot = (postprocset["PLOTSET"]
                         ['filename']+'_results_step-'+str(step))
            if 'edge' in postprocset["PLOTSET"].keys():
                plotset['edge'] = postprocset["PLOTSET"]['edge']
            else:
                plotset['edge'] = False
            if 'average' in postprocset["PLOTSET"]['result2plot'].keys():
                if postprocset["COMPUTER"]['average'] == True:
                    plotset['apply'] = 'points'
                else:
                    plotset['apply'] = 'cells'
            else:
                plotset['apply'] = 'cells'
            if 'intforces' in postprocset["PLOTSET"]['result2plot'].keys():
                if 'beam' in postprocset["PLOTSET"].keys():
                    nbeam = postprocset["PLOTSET"]['beam']
                else:
                    nbeam = [1]
                lenx = np.around(
                    postporc_result['balance'][step-1]['val'][0], decimals=3)
                leny = np.around(
                    postporc_result['balance'][step-1]['val'][1], decimals=3)
                xlabel = 'lenght ---> x'
                ylabel = postporc_result['balance'][step-1]['title']
                size = len(ylabel)
                plot_forces(lenx, leny, xlabel, ylabel, size, nbeam)
            if 'displ' in postprocset["PLOTSET"]['result2plot'].keys():
                post_show_mesh(file2plot, plotset)
            if 'frf' in postprocset["PLOTSET"]['result2plot'].keys():
                node_coordX = float(
                    postprocset["PLOTSET"]['result2plot']['frf']['point']['x'])
                node_coordY = float(
                    postprocset["PLOTSET"]['result2plot']['frf']['point']['y'])
                node_coordZ = float(
                    postprocset["PLOTSET"]['result2plot']['frf']['point']['z'])
                hist_node = search_nodexyz(
                    node_coordX, node_coordY, node_coordZ, modelinfo['coord'], 2E-3)
                hist_node = hist_node[0]
                plotset['fignumb'] = 3
                plotset['rstl'] = modelinfo['nodedof'][0]*hist_node - \
                    (modelinfo['nodedof'][0]-postprocset["PLOTSET"]
                     ['result2plot']['frf']['dof'])
                plotset['val_y'] = postporc_result['frf'][0]['val']
                plotset['val_x'] = postporc_result['frf'][0]['freqlog']
                frf_plot(plotset, hist_node)
            else:
                pass
    else:
        pass
