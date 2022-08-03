#!/usr/bin/env python
"""
Plot XY
"""
__author__ = "Antonio Vinicius Garcia Campos"
__copyright__ = "Copyright @ 2022, Antonio Vinicius Garcia Campos"
__credits__ = ["Antonio Vinicius Garcia Campos", "3D EasyCAE"]
__license__ = "GPL"
__status__ = "Development"

import numpy as np
import matplotlib.pyplot as plt
from myfempy.felib.physics.getnode import search_nodexyz


def plot(x, y, xlabel, ylabel, fignumb):
    plt.gcf().set_size_inches(10, 8)
    plt.plot(x, y, '-sm')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid(True)
    plt.draw()
    plt.pause(0.001)


def tracker_plot(postprocset, plotset, coord):
    if 'point' in postprocset["TRACKER"].keys():
        node_coordX = float(postprocset["TRACKER"]['point']['x'])
        node_coordY = float(postprocset["TRACKER"]['point']['y'])
        node_coordZ = float(postprocset["TRACKER"]['point']['z'])
        hist_node = search_nodexyz(
            node_coordX, node_coordY, node_coordZ, coord, 2E-3)
        hist_node = hist_node[0]
        if postprocset["TRACKER"]['result2plot'] == 'stress':
            val_Y = plotset['val_list'][hist_node-1, plotset['rstl'][0]]
            val_X = plotset['val_list'][hist_node-1, plotset['rstl'][1]]
            xlabel = 'STRAIN'
            ylabel = ('STRESS NODE: '+str(hist_node))
        elif postprocset["TRACKER"]['result2plot'] == 'displ':
            val_Y = plotset['val_list'][hist_node-1, 0]
            val_X = plotset['step']
            xlabel = 'STEP'
            ylabel = ('DISPL NODE: '+str(hist_node))
    elif 'max' in postprocset["TRACKER"].keys():
        if postprocset["TRACKER"]['result2plot'] == 'stress':
            val_Y = max(abs(plotset['val_list'][:, 0]))
            val_X = max(abs(plotset['val_list'][:, 4]))
            xlabel = 'STRAIN'
            ylabel = ('STRESS VM MAX')
        elif postprocset["TRACKER"]['result2plot'] == 'displ':
            val_Y = max(abs(plotset['val_list'][:, 0]))
            val_X = plotset['step']
            xlabel = 'STEP'
            ylabel = ('DISPL MAG MAX')
    elif 'min' in postprocset["TRACKER"].keys():
        if postprocset["TRACKER"]['result2plot'] == 'stress':
            val_Y = min(abs(plotset['val_list'][:, 0]))
            val_X = min(abs(plotset['val_list'][:, 4]))
            xlabel = 'STRAIN'
            ylabel = ('STRESS VM MIN')
        elif postprocset["TRACKER"]['result2plot'] == 'displ':
            val_Y = min(abs(plotset['val_list'][:, 0]))
            val_X = plotset['step']
            xlabel = 'STEP'
            ylabel = ('DISPL MAG MIN')
    else:
        val_X = 0
        val_Y = 0
        xlabel = 'erro'
        ylabel = 'erro'
    plot(val_X, val_Y, xlabel, ylabel, plotset['fignumb'])


def plot_forces(lenx, leny, xlabel, yl, size, nbeam):
    plt.gcf().set_size_inches(16, 8)
    plt.subplots_adjust(hspace=0.5)
    cont = 1
    for ff in range(len(leny)):
        for bb in range(len(nbeam)):
            x = lenx[:, nbeam[bb]-1]
            y = leny[ff][:, nbeam[bb]-1]
            ylabel = (yl[ff]+'_beam_'+str(nbeam[bb]))

            plt.subplot(size, len(nbeam), cont)
            plt.plot(x, y, '-sc')
            plt.xlabel(xlabel)
            plt.ylabel(ylabel)
            plt.grid(True)
            cont += 1
    plt.show()


def frf_plot(plotset, hist_node):
    # np.log(abs(np.ravel(plotset['val_y'][plotset['rstl'],:])))
    val_Y = np.log(abs(plotset['val_y'][plotset['rstl'], :]))
    val_X = plotset['val_x']
    xlabel = 'FREQ HZ'
    ylabel = 'DISPL'
    plot(val_X, val_Y, xlabel, ylabel, plotset['fignumb'])
