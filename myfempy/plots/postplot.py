#!/usr/bin/env python
import numpy as np

from myfempy.core.utilities import search_nodexyz
from myfempy.io.iocsv import write2log, writer2csv
from myfempy.plots.plotmesh import post_show_mesh
from myfempy.plots.plotxy import forces_plot, frf_plot, tracker_plot


__docformat__ = "google"

__doc__ = """

==========================================================================
                            __                                
         _ __ ___   _   _  / _|  ___  _ __ ___   _ __   _   _ 
        | '_ ` _ \ | | | || |_  / _ \| '_ ` _ \ | '_ \ | | | |
        | | | | | || |_| ||  _||  __/| | | | | || |_) || |_| |
        |_| |_| |_| \__, ||_|   \___||_| |_| |_|| .__/  \__, |
                    |___/                       |_|     |___/ 
        myfempy -- MultiphYsics Finite Element Module to PYthon    
                    COMPUTATIONAL ANALYSIS PROGRAM                   
        Copyright (C) 2022-2026 Antonio Vinicius Garcia Campos        
==========================================================================
This Python file is part of myfempy project.

myfempy is a python package based on finite element method to multiphysics
analysis. The code is open source and *intended for educational and scientific
purposes only, not recommended to commercial use. The name myfempy is an acronym
for MultiphYsics Finite Elements Module to PYthon. You can help us by contributing
with the main project, send us a mensage on https://github.com/avgcampos/myfempy/discussions/10
If you use myfempy in your research, the  developers would be grateful if you 
could cite in your work.
																		
The code is written by Antonio Vinicius Garcia Campos.                                  
																		
A github repository, with the most up to date version of the code,      
can be found here: https://github.com/avgcampos/myfempy.                 
																		
The code is open source and intended for educational and scientific     
purposes only. If you use myfempy in your research, the developers      
would be grateful if you could cite this. The myfempy project is published
under the GPLv3, see the myfempy LICENSE on
https://github.com/avgcampos/myfempy/blob/main/LICENSE.
																		
Disclaimer:                                                             
The authors reserve all rights but do not guarantee that the code is    
free from errors. Furthermore, the authors shall not be liable in any   
event caused by the use of the program.

"""

def postproc_plot(postprocset: dict, postporc_result: dict, modelinfo: dict):
    plotset = dict()
    hist_X = []
    hist_Y = []

    if "TRACKER" in postprocset.keys():
        if postprocset["TRACKER"]["show"]:
            for st in range(len(postporc_result["displ"])):
                plotset["step"] = st + 1
                plotset["val_list"] = postporc_result["displ"][st]["val"]
                plotset["fignumb"] = 99
                plotset["rstl"] = [0, modelinfo["ntensor"][0] + 1]
                val_X, val_Y, xlabel, ylabel = tracker_plot(
                    postprocset, plotset, modelinfo["coord"], modelinfo["nodedof"]
                )
                hist_X.append(val_X)
                hist_Y.append(val_Y)

            # hist_X = hist_X[1::][::]
            # hist_Y = hist_Y[1::][::]
            writer2csv("TRACKER_csvfile.txt", [hist_X, hist_Y], [xlabel, ylabel])

    if "PLOTSET" in postprocset.keys():
        if postprocset["PLOTSET"]["show"]:
            if "step" in postprocset["PLOTSET"].keys():
                step = int(postprocset["PLOTSET"]["step"])
            else:
                step = int(1)
            plotset["text_plot"] = "DISPL" + " step: " + str(step)
            plotset["step"] = step - 1
            file2plot = (
                postprocset["PLOTSET"]["filename"] + "_results_step-" + str(step)
            )

            if "edge" in postprocset["PLOTSET"].keys():
                plotset["edge"] = postprocset["PLOTSET"]["edge"]
            else:
                plotset["edge"] = False

            if "average" in postprocset["PLOTSET"]["data"].keys():
                if postprocset["COMPUTER"]["average"]:
                    plotset["apply"] = "points"
                else:
                    plotset["apply"] = "cells"
            else:
                plotset["apply"] = "cells"

            if "displ" in postprocset["PLOTSET"]["data"].keys():
                post_show_mesh(file2plot, plotset)

            if "intforces_plot" in postprocset["PLOTSET"]["data"].keys():
                if "beam" in postprocset["PLOTSET"]["data"]["intforces"].keys():
                    nbeam = postprocset["PLOTSET"]["data"]["intforces"]["beam"]
                else:
                    nbeam = [1]
                lenx = np.around(
                    postporc_result["intforces_plot"][step - 1]["val"][0], decimals=3
                )
                leny = np.around(
                    postporc_result["intforces_plot"][step - 1]["val"][1], decimals=3
                )
                xlabel = "lenght ---> x"
                ylabel = postporc_result["intforces_plot"][step - 1]["title"]
                size = len(ylabel)
                forces_plot(lenx, leny, xlabel, ylabel, size, nbeam)

            if "frf" in postprocset["PLOTSET"]["data"].keys():
                node_coordX = float(postprocset["PLOTSET"]["data"]["frf"]["point"]["x"])
                node_coordY = float(postprocset["PLOTSET"]["data"]["frf"]["point"]["y"])
                node_coordZ = float(postprocset["PLOTSET"]["data"]["frf"]["point"]["z"])
                hist_node = search_nodexyz(
                    node_coordX, node_coordY, node_coordZ, modelinfo["coord"], 2e-3
                )
                hist_node = hist_node[0]
                plotset["fignumb"] = 3
                plotset["filename"] = postprocset["PLOTSET"]["filename"]
                plotset["rstl"] = modelinfo["nodedof"][0] * hist_node - (
                    modelinfo["nodedof"][0]
                    - postprocset["PLOTSET"]["data"]["frf"]["dof"]
                )
                plotset["val_y"] = 20 * np.log(
                    (abs(postporc_result["frf"][0]["val"][plotset["rstl"], :])) / 10e-12
                )
                plotset["val_x"] = postporc_result["frf"][0]["freqlog"]
                frf_plot(plotset, hist_node)
            else:
                pass

    if "OUTPUT" in postprocset.keys():
        write2log(
            "OUTPUT_logfile.txt",
            postprocset["OUTPUT"],
            modelinfo,
            postprocset["SOLUTION"]["solvestatus"],
        )

    else:
        pass
