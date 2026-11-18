from __future__ import annotations

import numpy as np

from myfempy.core.utilities import (gauss_points, results_average,
                                    search_nodexyz)
from myfempy.io.iocsv import write2log, writer2csv
from myfempy.io.iovtk import convert_to_vtk

__docformat__ = "google"

__doc__ = """
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


class setPostProcess:
    """PostProcess Class <ClassOrder>"""

    def getCompute(self, postprocset):
        # print_console("post")
        postprocdata = dict()
        postprocdata["SOLUTION"] = []
        SOLUTION = postprocset["SOLVERDATA"]["solution"]["U"]

        postporc_result = dict()

        if "structural" in postprocset["COMPUTE"].keys():
            result_solu = np.zeros((SOLUTION.shape[1], self.model.coord.shape[0], 3))

            for st in range(SOLUTION.shape[1]):
                result_solu[st, :, :], sol_title = setPostProcess.__displ(
                    self, SOLUTION[:, st]
                )

            if (
                "displ" in postprocset["COMPUTE"]["structural"].keys()
                and postprocset["COMPUTE"]["structural"]["displ"] == True
            ):
                postporc_result["solution"] = []
                for st in range(SOLUTION.shape[1]):
                    postporc_result["solution"].append(
                        {
                            "val": result_solu[st][:][:],
                            "title": sol_title,
                            "avr": True,
                        }
                    )

                    postprocdata["SOLUTION"].append(
                        {sol_title: result_solu[st, :, :], "STEP": st + 1}
                    )

                if (
                    "stress" in postprocset["COMPUTE"]["structural"].keys()
                    and postprocset["COMPUTE"]["structural"]["stress"] == True
                ):
                    
                    result_stress = np.zeros(
                        (
                            SOLUTION.shape[1],
                            self.model.inci.shape[0],
                            2 * self.model.modelinfo["tensor"] + 4,
                        )
                    )

                    for st in range(SOLUTION.shape[1]):
                        result_stress[st, :, :], title = setPostProcess.__stress(
                            self, SOLUTION[:, st]
                        )

                    postprocdata["STRESS"] = {}
                    for setpost in range(result_stress.shape[2]):
                        post_list_step = []
                        for st in range(SOLUTION.shape[1]):
                            post_list_step.append({'ELEMVAL': result_stress[st, :, setpost], "STEP": st + 1})
                        postprocdata["STRESS"][title[setpost]] = post_list_step

                    
                    postporc_result["stress_elm"] = []
                    for st in range(SOLUTION.shape[1]):
                        postporc_result["stress_elm"].append(
                            {"val": result_stress[st][:][:], "title": title, "avr": False}
                        )

            if (
                "modes" in postprocset["COMPUTE"]["structural"].keys()
                and postprocset["COMPUTE"]["structural"]["modes"] == True
            ):
                FREQUENCY = postprocset["SOLVERDATA"]["solution"]["FREQ"]
                postporc_result["solution"] = []
                for st in range(SOLUTION.shape[1]):
                    postporc_result["solution"].append(
                        {
                            "val": result_solu[st][:][:],
                            "title": (
                                "MODE_"
                                + str(st + 1)
                                + "_"
                                + "FREQ_"
                                + str(np.round(FREQUENCY[st][2], 3))
                                + "Hz"
                            ),
                            "avr": True,
                        }
                    )

                    postprocdata["SOLUTION"].append(
                        {
                            "VAL": result_solu[st, :, :],
                            "MODE": st + 1,
                            "FREQ": np.round(FREQUENCY[st][2], 3),
                        }
                    )

            if (
                "freq" in postprocset["COMPUTE"]["structural"].keys()
                and postprocset["COMPUTE"]["structural"]["freq"] == True
            ):
                FREQUENCY = postprocset["SOLVERDATA"]["solution"]["FREQ"]
                postporc_result["solution"] = []
                for st in range(FREQUENCY.shape[0]):
                    postporc_result["solution"].append(
                        {
                            "val": result_solu[st][:][:],
                            "title": (
                                "RESPONSE_"
                                + str(st)
                                + "_FREQ_"
                                + ",".join(map(str, np.round(FREQUENCY[st][:], 3).flatten())) #str(np.round(FREQUENCY[st][:], 3))
                                + "Hz"
                            ),
                            # "freq": FREQUENCY[st][:].tolist(),
                            "avr": True,
                        }
                    )
                    
                    postprocdata["SOLUTION"].append(
                        {
                            "VAL": result_solu[st, :, :],
                            "FREQ": FREQUENCY[st][:].tolist(),
                        }
                    )
            else:
                pass

        if "thermal" in postprocset["COMPUTE"].keys():
            result_solu = np.zeros((SOLUTION.shape[1], self.model.coord.shape[0], 1))
            for st in range(SOLUTION.shape[1]):
                result_solu[st, :, :], sol_title = setPostProcess.__displ(
                    self, SOLUTION[:, st]
                )
            if (
                "temp" in postprocset["COMPUTE"]["thermal"].keys()
                and postprocset["COMPUTE"]["thermal"]["temp"] == True
            ):
                postporc_result["solution"] = []
                for st in range(SOLUTION.shape[1]):
                    postporc_result["solution"].append(
                        {
                            "val": result_solu[st][:][:],
                            "title": sol_title,
                            "avr": True,
                        }
                    )

                    postprocdata["SOLUTION"].append(
                        {sol_title: result_solu[st, :, :], "STEP": st + 1}
                    )

                if (
                    "heatflux" in postprocset["COMPUTE"]["thermal"].keys()
                    and postprocset["COMPUTE"]["thermal"]["heatflux"] == True
                ):
                    result_stress = np.zeros(
                        (
                            SOLUTION.shape[1],
                            self.model.inci.shape[0],
                            2 * self.model.modelinfo["tensor"] + 2,
                        )
                    )

                    for st in range(SOLUTION.shape[1]):
                        result_stress[st, :, :], title = setPostProcess.__heatflux(
                            self, SOLUTION[:, st]
                        )

                    postporc_result["stress_elm"] = []
                    for st in range(SOLUTION.shape[1]):
                        postporc_result["stress_elm"].append(
                            {"val": result_stress[st][:][:], "title": title, "avr": False}
                        )

                        for setpost in range(result_stress.shape[2]):
                            postprocdata[title[setpost]] = result_stress[st, :, setpost]

            else:
                pass
        else:
            pass

        # save in vtk file
        setPostProcess.__tovtkplot(self, postprocset, postporc_result)

        return postprocdata

    def getPlotCSV(self, postprocset, postporc_result):
        plotset = dict()
        hist_X = []
        hist_Y = []
        plotset["nodedof"] = self.model.modelinfo["nodedof"]
        if "data" in postprocset["PLOT"].keys():
            hist_Y = (postprocset["PLOT"]["data"]["y_data"]).tolist()
            hist_X = (postprocset["PLOT"]["data"]["x_data"]).tolist()
            xlabel = postprocset["PLOT"]["data"]["x_label"]
            ylabel = postprocset["PLOT"]["data"]["y_label"]
        else:
            max_steps = len(postporc_result["compute"]["SOLUTION"])
            for st in range(max_steps):
                plotset["step"] = st
                val_X, val_Y, xlabel, ylabel = setPostProcess.__hist_value(
                    postporc_result["compute"], postprocset, plotset, self.model.coord
                )
                hist_X.append(val_X)
                hist_Y.append([val_Y])
        filename = (
            str(self.path)
            + "/"
            + postprocset["POSTSET"]["filename"]
            + "_csv-file.txt"
        )
        writer2csv(filename, [hist_X, hist_Y], [xlabel, ylabel])

    def getLog(self, postprocset, postporc_result):
        filename = (
            str(self.path)
            + "/"
            + postprocset["POSTSET"]["filename"]
            + "_solver-log.txt"
        )
        postprocset["REPORT"]["timesolver"] = self.timenow
        write2log(
            self.model, self.physic, postprocset["REPORT"], postporc_result, filename
        )

    def __tovtkplot(self, postprocset, postporc_result):
        # path = os.getcwd()
        plotdata = dict()

        plotdata["inci"] = self.model.inci
        plotdata["nodecon"] = self.model.modelinfo["nodecon"]
        plotdata["filename"] = (
            str(self.path)
            + "/"
            + postprocset["POSTSET"]["filename"]
            + "_myfempy_postprocess"
        )
        plotdata["coord"] = self.model.coord
        plotdata["material_CELL_DATA_val"] = (
            (np.array([self.model.inci[:, 2]])).T
        ).astype(int)
        plotdata["material_CELL_DATA_title"] = ["TOPOLOGY"]

        if "structural" in postprocset["COMPUTE"].keys():
            if "displ" in postprocset["COMPUTE"]["structural"].keys():

                plotdata["displ_POINT_DATA_val"] = []
                plotdata["displ_POINT_DATA_title"] = []
                plotdata["stress_CELL_DATA_val"] = []
                plotdata["stress_CELL_DATA_title"] = []
                # ... any field in future

                for st in range(len(postporc_result["solution"])):
                    plotdata["filename"] = (
                        str(self.path)
                        + "/"
                        + postprocset["POSTSET"]["filename"]
                        + "_myfempy_postprocess_step-"
                        + str(st)
                    )

                    plotdata["displ_POINT_DATA_val"] = postporc_result["solution"][st][
                        "val"
                    ][:, :]
                    plotdata["displ_POINT_DATA_title"] = postporc_result["solution"][
                        st
                    ]["title"]

                    if (
                        "stress" in postprocset["COMPUTE"]["structural"].keys()
                        and postprocset["COMPUTE"]["structural"]["stress"] == True
                    ):
                        plotdata["stress_CELL_DATA_val"] = postporc_result[
                            "stress_elm"
                        ][st]["val"]
                        plotdata["stress_CELL_DATA_title"] = postporc_result[
                            "stress_elm"
                        ][st]["title"]
                    convert_to_vtk(plotdata)

            elif "modes" in postprocset["COMPUTE"]["structural"].keys():
                plotdata["modes_POINT_DATA"] = []
                plotdata["modes_POINT_DATA"] = postporc_result["solution"]
                convert_to_vtk(plotdata)

            elif "freq" in postprocset["COMPUTE"]["structural"].keys():
                plotdata["freq_POINT_DATA"] = []
                plotdata["freq_POINT_DATA"] = postporc_result["solution"]
                convert_to_vtk(plotdata)

        elif "thermal" in postprocset["COMPUTE"].keys():

            if "temp" in postprocset["COMPUTE"]["thermal"].keys():

                plotdata["temp_POINT_DATA_val"] = []
                plotdata["temp_POINT_DATA_title"] = []
                plotdata["heatflux_CELL_DATA_val"] = []
                plotdata["heatflux_CELL_DATA_title"] = []
                # ... any field in future

                for st in range(len(postporc_result["solution"])):
                    plotdata["filename"] = (
                        str(self.path)
                        + "/"
                        + postprocset["POSTSET"]["filename"]
                        + "_post_process_step-"
                        + str(st)
                    )

                    plotdata["temp_POINT_DATA_val"] = postporc_result["solution"][st][
                        "val"
                    ][:, :]
                    plotdata["temp_POINT_DATA_title"] = postporc_result["solution"][st][
                        "title"
                    ]

                    if (
                        "heatflux" in postprocset["COMPUTE"]["thermal"].keys()
                        and postprocset["COMPUTE"]["thermal"]["heatflux"] == True
                    ):
                        plotdata["stress_CELL_DATA_val"] = postporc_result[
                            "stress_elm"
                        ][st]["val"]
                        plotdata["stress_CELL_DATA_title"] = postporc_result[
                            "stress_elm"
                        ][st]["title"]

                    convert_to_vtk(plotdata)

            else:
                pass

    def __displ(self, U):
        disp = self.model.element.getElementDeformation(U, self.model.modelinfo)
        title = self.model.element.setTitleDeformation()
        return disp, title

    def __stress(self, U):
        stress_list = np.zeros(
            (self.model.modelinfo["nelem"], self.model.modelinfo["tensor"] + 1),
            dtype=float,
        )

        strain_list = np.zeros(
            (self.model.modelinfo["nelem"], self.model.modelinfo["tensor"] + 1),
            dtype=float,
        )

        compliance_list = np.zeros((self.model.modelinfo["nelem"], 1), dtype=float)

        factor_of_safety = np.zeros((self.model.modelinfo["nelem"], 1), dtype=float)

        pt, wt = gauss_points(self.model.modelinfo["shape"], 1)

        for ee in range(self.model.modelinfo["nelem"]):

            epsilon, strain = self.model.material.getElementStrain(
                self.model, U, pt[0], ee
            )

            sigma, stress = self.model.material.getElementStress(
                self.model, epsilon, ee
            )

            strain_energy = self.model.material.getStrainEnergyDensity(
                sigma, epsilon, self.model.elemvol[ee]
            )

            FoS = self.model.material.getFailureCriteria(stress)

            stress_list[ee, :] = stress
            strain_list[ee, :] = strain
            compliance_list[ee, :] = strain_energy
            factor_of_safety[ee, :] = FoS

        result = np.concatenate(
            (stress_list, strain_list, compliance_list, factor_of_safety), axis=1
        )

        tistrs = self.model.material.getTitleStress()
        tistrn = self.model.material.getTitleStrain()
        ticomp = self.model.material.getTitleCompliance()
        tifos = self.model.material.getTitleFoS()
        title = np.concatenate((tistrs, tistrn, ticomp, tifos), axis=0)
        return result, title

    def __heatflux(self, U):
        stress_list = np.zeros(
            (self.model.modelinfo["nelem"], self.model.modelinfo["tensor"] + 1),
            dtype=float,
        )

        strain_list = np.zeros(
            (self.model.modelinfo["nelem"], self.model.modelinfo["tensor"] + 1),
            dtype=float,
        )

        pt, wt = gauss_points(self.model.modelinfo["shape"], 1)

        for ee in range(self.model.modelinfo["nelem"]):

            epsilon, strain = self.model.material.getElementGradTemp(
                self.model, U, pt[0], ee
            )

            sigma, stress = self.model.material.getElementHeatFlux(
                self.model, epsilon, ee
            )

            stress_list[ee, :] = stress
            strain_list[ee, :] = strain

        result = np.concatenate((stress_list, strain_list), axis=1)

        tistrs = self.model.material.getTitleHeatFlux()
        tistrn = self.model.material.getTitleGradTemp()
        title = np.concatenate((tistrs, tistrn), axis=0)
        return result, title

    def __hist_value(postporc_result, postprocset, plotset, coord):
        if "u" in postprocset["PLOT"].keys():
            node_coordX = float(postprocset["PLOT"]["u"]["x"])
            node_coordY = float(postprocset["PLOT"]["u"]["y"])
            node_coordZ = float(postprocset["PLOT"]["u"]["z"])
            hist_node = search_nodexyz(
                node_coordX, node_coordY, node_coordZ, coord, 1e-6
            )
            val_Y = postporc_result["SOLUTION"][plotset["step"]][postprocset["PLOT"]["u"]["val"]][
                hist_node[0] - 1, postprocset["PLOT"]["u"]["dof"] - 1
            ]
            val_X = plotset["step"] + 1
            xlabel = "STEP"
            ylabel = "U NODE: " + str(hist_node)

        elif "freq" in postprocset["PLOT"].keys():
            node_coordX = float(postprocset["PLOT"]["freq"]["x"])
            node_coordY = float(postprocset["PLOT"]["freq"]["y"])
            node_coordZ = float(postprocset["PLOT"]["freq"]["z"])
            hist_node = search_nodexyz(
                node_coordX, node_coordY, node_coordZ, coord, 1e-6
            )
            val_Y = abs(
                postporc_result["SOLUTION"][plotset["step"]]["VAL"][
                    hist_node[0] - 1, postprocset["PLOT"]["freq"]["dof"] - 1
                ]
            )
            val_X = postporc_result["SOLUTION"][plotset["step"]]["FREQ"]
            xlabel = "FREQUENCY [Hz]"
            ylabel = "RESPONSE [mag(U)]"

        elif "frf" in postprocset["PLOT"].keys():
            node_coordX = float(postprocset["PLOT"]["frf"]["x"])
            node_coordY = float(postprocset["PLOT"]["frf"]["y"])
            node_coordZ = float(postprocset["PLOT"]["frf"]["z"])
            hist_node = search_nodexyz(
                node_coordX, node_coordY, node_coordZ, coord, 1e-6
            )
            val_Y = 20 * np.log10(
                (
                    abs(
                        postporc_result["SOLUTION"][plotset["step"]]["VAL"][
                            hist_node[0] - 1,
                            postprocset["PLOT"]["frf"]["dof"] - 1,
                        ]
                    )
                )
            )
            val_X = postporc_result["SOLUTION"][plotset["step"]]["FREQ"]
            xlabel = "FREQUENCY [Hz]"
            ylabel = "RESPONSE [dB]"

        elif "tc" in postprocset["PLOT"].keys():
            node_coordX_in = float(postprocset["PLOT"]["tc"]['indof']["x"])
            node_coordY_in = float(postprocset["PLOT"]["tc"]['indof']["y"])
            node_coordZ_in = float(postprocset["PLOT"]["tc"]['indof']["z"])
            hist_node_in = search_nodexyz(
                node_coordX_in, node_coordY_in, node_coordZ_in, coord, 1e-6
            )
            node_coordX_out = float(postprocset["PLOT"]["tc"]['outdof']["x"])
            node_coordY_out = float(postprocset["PLOT"]["tc"]['outdof']["y"])
            node_coordZ_out = float(postprocset["PLOT"]["tc"]['outdof']["z"])
            hist_node_out = search_nodexyz(
                node_coordX_out, node_coordY_out, node_coordZ_out, coord, 1e-6
            )
            H_in = postporc_result["SOLUTION"][plotset["step"]]["VAL"][
                    hist_node_in[0] - 1, postprocset["PLOT"]["tc"]['indof']["dof"] - 1]
            
            H_out = postporc_result["SOLUTION"][plotset["step"]]["VAL"][
                    hist_node_out[0] - 1, postprocset["PLOT"]["tc"]['outdof']["dof"] - 1]
            
            val_Y = 20*np.log10(np.abs(H_out  / (H_in + 1e-30)))
            val_X = postporc_result["SOLUTION"][plotset["step"]]["FREQ"]
            xlabel = "FREQUENCY [Hz]"
            ylabel = "TRANSMISSIBILITY [dB]"

        else:
            val_X = 0.0
            val_Y = 0.0
            xlabel = "NONE"
            ylabel = "NONE"
        return val_X, val_Y, xlabel, ylabel
