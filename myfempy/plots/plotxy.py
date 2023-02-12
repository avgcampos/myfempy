#!/usr/bin/env python
from myfempy.felib.physics.getnode import search_nodexyz
from myfempy.io.iocsv import writer2csv
import matplotlib.pyplot as plt
import numpy as np

__doc__ = """
Plot XY
"""


def plot(x, y, xlabel, ylabel, fignumb):
    """_summary_

    Arguments:
        x -- _description_
        y -- _description_
        xlabel -- _description_
        ylabel -- _description_
        fignumb -- _description_
    """
    plt.ion()
    plt.gcf().set_size_inches(10, 8)
    plt.plot(x, y, "-sm")
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid(True)
    plt.draw()
    plt.pause(0.001)


def tracker_plot(postprocset, plotset, coord, nodedof):
    """_summary_

    Arguments:
        postprocset -- _description_
        plotset -- _description_
        coord -- _description_
    """
    if "displ" in postprocset["TRACKER"]['data'].keys():
        if "point" in postprocset["TRACKER"]['data']['displ'].keys():
            node_coordX = float(postprocset["TRACKER"]['data']['displ']["point"]["x"])
            node_coordY = float(postprocset["TRACKER"]['data']['displ']["point"]["y"])
            node_coordZ = float(postprocset["TRACKER"]['data']['displ']["point"]["z"])
            hist_node = search_nodexyz(node_coordX, node_coordY, node_coordZ, coord, 1e-6)
            hist_node = hist_node[0]
            # plotset["rstl"] = nodedof[0]*hist_node - (nodedof[0]-postprocset["TRACKER"]["data"]["displ"]["dof"])
            val_Y = plotset["val_list"][hist_node-1, postprocset["TRACKER"]["data"]["displ"]["dof"]]
            val_X = plotset["step"]
            xlabel = "STEP"
            ylabel = "DISPL NODE: " + str(hist_node)
        elif "max" in postprocset["TRACKER"]['data']['displ'].keys():
            val_Y = max(abs(plotset["val_list"][:, 0]))
            val_X = plotset["step"]
            xlabel = "STEP"
            ylabel = "DISPL MAG MAX"
        elif "min" in postprocset["TRACKER"]['data']['displ'].keys():
            val_Y = min(abs(plotset["val_list"][:, 0]))
            val_X = plotset["step"]
            xlabel = "STEP"
            ylabel = "DISPL MAG MIN"
        else:
            val_X = 0
            val_Y = 0
            xlabel = "erro"
            ylabel = "erro"
        # plot(val_X, val_Y, xlabel, ylabel, plotset["fignumb"])
    return val_X, val_Y, xlabel, ylabel


def forces_plot(lenx, leny, xlabel, yl, size, nbeam):
    """_summary_

    Arguments:
        lenx -- _description_
        leny -- _description_
        xlabel -- _description_
        yl -- _description_
        size -- _description_
        nbeam -- _description_
    """
    plt.gcf().set_size_inches(16, 8)
    plt.subplots_adjust(hspace=0.5)
    cont = 1
    for ff in range(len(leny)):
        for bb in range(len(nbeam)):
            val_X = lenx[:, nbeam[bb] - 1]
            val_Y = leny[ff][:, nbeam[bb] - 1]
            ylabel = yl[ff] + "_beam_" + str(nbeam[bb])

            plt.subplot(size, len(nbeam), cont)
            plt.plot(val_X, val_Y, "-sc")
            plt.xlabel(xlabel)
            plt.ylabel(ylabel)
            plt.grid(True)
            cont += 1
            writer2csv_file('BEAM-BALANCE_'+str(nbeam[bb])+'_csvfile.txt',[val_X,val_Y],[xlabel,ylabel])
    plt.show()

def frf_plot(plotset, hist_node):
    """_summary_

    Arguments:
        plotset -- _description_
        hist_node -- _description_
    """
    val_Y = plotset["val_y"] #20*np.log((abs(plotset["val_y"][plotset["rstl"], :]))/10E-12)
    val_X = plotset["val_x"]
    xlabel = "FREQ HZ"
    ylabel = "dB 20*LOG[DISPL/10E-12 m]"
    plot(val_X, val_Y, xlabel, ylabel, plotset["fignumb"])
    writer2csv_file('FRF_csvfile.txt',[val_X,val_Y],['Hz','dB(20LOG(U))'])
