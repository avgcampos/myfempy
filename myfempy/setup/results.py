from __future__ import annotations

import os
from abc import ABC, abstractmethod

import numpy as np

from myfempy.core.utilities import (gauss_points, results_average,
                                    search_nodexyz)
from myfempy.io.iocsv import write2log, writer2csv
from myfempy.io.iovtk import convert_to_vtk


class setPostProcess(ABC):
    """PostProcess Class <ClassOrder>"""

    def getCompute(self, postprocset):
        """_summary_

        Arguments:
            postprocset -- _description_

        Returns:
            _description_
        """
        # print_console("post")

        postprocdata = dict()
        postprocdata["SOLUTION"] = []
        postprocdata["solverstatus"] = postprocset["SOLVERDATA"]["solverstatus"]
        SOLUTION = postprocset["SOLVERDATA"]["solution"]["U"]

        postporc_result = dict()

        if "structural" in postprocset["COMPUTER"].keys():
            result_solu = np.zeros(
                (SOLUTION.shape[1], self.model.coord.shape[0], 3)
            )
            for ns in range(SOLUTION.shape[1]):
                result_solu[ns, :, :], sol_title = setPostProcess.__displ(
                    self, SOLUTION[:, ns]
                )
            if "displ" in postprocset["COMPUTER"]["structural"].keys():
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
                "stress" in postprocset["COMPUTER"]["structural"].keys()
                and postprocset["COMPUTER"]["structural"]["stress"] == True
            ):
                result_stress = np.zeros(
                    (
                        SOLUTION.shape[1],
                        self.model.inci.shape[0],
                        2 * self.model.modelinfo["tensor"] + 4,
                    )
                )

                for st in range(SOLUTION.shape[1]):
                    result_stress[ns, :, :], title = setPostProcess.__stress(
                        self, SOLUTION[:, ns]
                    )

                    for setpost in range(result_stress.shape[2]):
                        postprocdata[title[setpost]] = result_stress[st, :, setpost]

                postporc_result["stress_elm"] = []
                for st in range(SOLUTION.shape[1]):
                    postporc_result["stress_elm"].append(
                        {"val": result_stress[st][:][:], "title": title, "avr": False}
                    )
            else:
                pass

            if "modes" in postprocset["COMPUTER"]["structural"].keys():
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
                "frf" in postprocset["COMPUTER"]["structural"].keys()
                and postprocset["COMPUTER"]["structural"]["frf"] == True
            ):
                FREQUENCY = postprocset["SOLVERDATA"]["solution"]["FREQ"]
                postporc_result["solution"] = []
                for st in range(SOLUTION.shape[1]):
                    postporc_result["solution"].append(
                        {
                            "val": result_solu[st][:][:],
                            "title": (
                                "RESPONSE_"
                                + str(st)
                                + "_FREQ_"
                                + str(np.round(FREQUENCY[st], 3))
                                + "Hz"
                            ),
                            "freq": FREQUENCY[st],
                            "avr": True,
                        }
                    )

                    postprocdata["SOLUTION"].append(
                        {
                            "VAL": result_solu[st, :, :],
                            "FREQ": np.round(FREQUENCY[st], 3),
                        }
                    )
            else:
                pass

        if "thermal" in postprocset["COMPUTER"].keys():
            result_solu = np.zeros(
                (SOLUTION.shape[1], self.model.coord.shape[0], 1)
            )
            for st in range(SOLUTION.shape[1]):
                result_solu[st, :, :], sol_title = setPostProcess.__displ(
                    self, SOLUTION[:, st]
                )
            if "temp" in postprocset["COMPUTER"]["thermal"].keys():
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
                "heatflux" in postprocset["COMPUTER"]["thermal"].keys()
                and postprocset["COMPUTER"]["thermal"]["heatflux"] == True
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

        setPostProcess.__tovtkplot(
            self, postprocset, postporc_result
        )  # save in vtk file

        return postprocdata

    def getTracker(self, postprocset, postporc_result):
        plotset = dict()
        hist_X = []
        hist_Y = []
        plotset["nodedof"] = self.model.modelinfo["nodedof"]

        max_steps = len(postporc_result["SOLUTION"])

        # if 'frf' in postprocset["TRACKER"].keys():
        #     max_steps = len(postporc_result["frf"])

        for st in range(max_steps):
            plotset["step"] = st
            # plotset["val_list"] = postporc_result['solution'][st]['val']
            # plotset["fignumb"] = 99
            # plotset["rstl"] = [0, modelinfo["ntensor"][0] + 1]
            val_X, val_Y, xlabel, ylabel = setPostProcess.__tracker_value(
                postporc_result, postprocset, plotset, self.model.coord
            )
            hist_X.append(val_X)
            hist_Y.append(val_Y)

        # hist_X = hist_X[1::][::]
        # hist_Y = hist_Y[1::][::]
        # path = os.getcwd()
        filename = (
            str(self.path)
            + "/"
            + postprocset["PLOTSET"]["filename"]
            + "_myfempy_csv-file.txt"
        )
        writer2csv(filename, [hist_X, hist_Y], [xlabel, ylabel])

    def getLog(self, postprocset, postporc_result):
        # path = os.getcwd()
        filename = (
            str(self.path)
            + "/"
            + postprocset["PLOTSET"]["filename"]
            + "_myfempy_solver-status.txt"
        )
        write2log(self.model, self.physic, postprocset["OUTPUT"], postporc_result, filename)

    def __tovtkplot(self, postprocset, postporc_result):
        """_summary_

        Arguments:
            postporc_result -- _description_
            postprocset -- _description_
        """
        # path = os.getcwd()
        plotdata = dict()

        plotdata["inci"] = self.model.inci
        plotdata["nodecon"] = self.model.modelinfo["nodecon"]
        plotdata["filename"] = (
            str(self.path) + "/" + postprocset["PLOTSET"]["filename"] + "_post_process"
        )
        plotdata["coord"] = self.model.coord
        plotdata["material_CELL_DATA_val"] = (
            (np.array([(self.model.inci[:, 2])])).T
        ).astype(int)
        plotdata["material_CELL_DATA_title"] = ["Material_Set"]

        # if "structural" in postprocset["COMPUTER"].keys():
        #     if "displ" in postprocset["COMPUTER"]["structural"].keys():
        if "structural" in postprocset["COMPUTER"].keys():

            if "displ" in postprocset["COMPUTER"]["structural"].keys():

                plotdata["displ_POINT_DATA_val"] = []
                plotdata["displ_POINT_DATA_title"] = []
                plotdata["stress_CELL_DATA_val"] = []
                plotdata["stress_CELL_DATA_title"] = []
                # ... any field in future

                for st in range(len(postporc_result["solution"])):
                    plotdata["filename"] = (
                        str(self.path)
                        + "/"
                        + postprocset["PLOTSET"]["filename"]
                        + "_post_process_step-"
                        + str(st)
                    )

                    plotdata["displ_POINT_DATA_val"] = postporc_result["solution"][st][
                        "val"
                    ][:, :]
                    plotdata["displ_POINT_DATA_title"] = postporc_result["solution"][
                        st
                    ]["title"]

                    if (
                        "stress" in postprocset["COMPUTER"]["structural"].keys()
                        and postprocset["COMPUTER"]["structural"]["stress"] == True
                    ):
                        plotdata["stress_CELL_DATA_val"] = postporc_result[
                            "stress_elm"
                        ][st]["val"]
                        plotdata["stress_CELL_DATA_title"] = postporc_result[
                            "stress_elm"
                        ][st]["title"]
                    convert_to_vtk(plotdata)

            elif "modes" in postprocset["COMPUTER"]["structural"].keys():
                plotdata["modes_POINT_DATA"] = []
                plotdata["modes_POINT_DATA"] = postporc_result["solution"]
                convert_to_vtk(plotdata)

            elif "frf" in postprocset["COMPUTER"]["structural"].keys():
                plotdata["frf_POINT_DATA"] = []
                plotdata["frf_POINT_DATA"] = postporc_result["solution"]
                convert_to_vtk(plotdata)

        elif "thermal" in postprocset["COMPUTER"].keys():

            if "temp" in postprocset["COMPUTER"]["thermal"].keys():

                plotdata["temp_POINT_DATA_val"] = []
                plotdata["temp_POINT_DATA_title"] = []
                plotdata["heatflux_CELL_DATA_val"] = []
                plotdata["heatflux_CELL_DATA_title"] = []
                # ... any field in future

                for st in range(len(postporc_result["solution"])):
                    plotdata["filename"] = (
                        str(self.path)
                        + "/"
                        + postprocset["PLOTSET"]["filename"]
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
                        "heatflux" in postprocset["COMPUTER"]["thermal"].keys()
                        and postprocset["COMPUTER"]["thermal"]["heatflux"] == True
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
        """_summary_

        Arguments:
            U -- _description_

        Returns:
            _description_
        """
        disp = self.model.element.getElementDeformation(U, self.model.modelinfo)
        title = self.model.element.setTitleDeformation()
        return disp, title

    def __stress(self, U):
        """_summary_

        Arguments:
            U -- _description_

        Returns:
            _description_
        """
        stress_list = np.zeros(
            (self.model.modelinfo["nelem"], self.model.modelinfo["tensor"] + 1), dtype=float
        )

        strain_list = np.zeros(
            (self.model.modelinfo["nelem"], self.model.modelinfo["tensor"] + 1), dtype=float
        )

        compliance_list = np.zeros((self.model.modelinfo["nelem"], 1), dtype=float)

        factor_of_safety = np.zeros((self.model.modelinfo["nelem"], 1), dtype=float)

        pt, wt = gauss_points(self.model.modelinfo["type_shape"], 1)

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
        """_summary_

        Arguments:
            U -- _description_

        Returns:
            _description_
        """
        stress_list = np.zeros(
            (self.model.modelinfo["nelem"], self.model.modelinfo["tensor"] + 1), dtype=float
        )

        strain_list = np.zeros(
            (self.model.modelinfo["nelem"], self.model.modelinfo["tensor"] + 1), dtype=float
        )

        pt, wt = gauss_points(self.model.modelinfo["type_shape"], 1)

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

    def __tracker_value(postporc_result, postprocset, plotset, coord):
        """trancker plot function"""

        if "point" in postprocset["TRACKER"].keys():
            node_coordX = float(postprocset["TRACKER"]["point"]["x"])
            node_coordY = float(postprocset["TRACKER"]["point"]["y"])
            node_coordZ = float(postprocset["TRACKER"]["point"]["z"])
            hist_node = search_nodexyz(
                node_coordX, node_coordY, node_coordZ, coord, 1e-3
            )
            # hist_node = hist_node[0]-1
            # plotset["rstl"] = nodedof[0]*hist_node - (nodedof[0]-postprocset["TRACKER"]["data"]["displ"]["dof"])
            val_Y = postporc_result["SOLUTION"][plotset["step"]]["val"][
                hist_node[0] - 1, postprocset["TRACKER"]["point"]["dof"]
            ]
            val_X = plotset["step"] + 1
            xlabel = "STEP"
            ylabel = "DISPL NODE: " + str(hist_node)

        elif "data" in postprocset["TRACKER"].keys():
            # node_coordX = float(postprocset["TRACKER"]["point"]["x"])
            # node_coordY = float(postprocset["TRACKER"]["point"]["y"])
            # node_coordZ = float(postprocset["TRACKER"]["point"]["z"])
            # hist_node = search_nodexyz(node_coordX, node_coordY, node_coordZ, coord, 1e-3)
            # hist_node = hist_node[0]-1
            # plotset["rstl"] = nodedof[0]*hist_node - (nodedof[0]-postprocset["TRACKER"]["data"]["displ"]["dof"])
            val_Y = postprocset["TRACKER"]["data"][
                "y_data"
            ]  # plotset["val_list"][hist_node, postprocset["TRACKER"]["point"]["dof"]]
            val_X = postprocset["TRACKER"]["data"]["x_data"]  # plotset["step"]
            xlabel = postprocset["TRACKER"]["data"]["x_label"]
            ylabel = postprocset["TRACKER"]["data"]["y_label"]

        # elif "max" in postprocset["TRACKER"].keys():
        #     val_Y = max(abs(plotset["val_list"][:, 0]))
        #     val_X = plotset["step"]
        #     xlabel = "STEP"
        #     ylabel = "DISPL MAG MAX"

        # elif "min" in postprocset["TRACKER"].keys():
        #     val_Y = min(abs(plotset["val_list"][:, 0]))
        #     val_X = plotset["step"]
        #     xlabel = "STEP"
        #     ylabel = "DISPL MAG MIN"

        elif "frf" in postprocset["TRACKER"].keys():
            node_coordX = float(postprocset["TRACKER"]["frf"]["x"])
            node_coordY = float(postprocset["TRACKER"]["frf"]["y"])
            node_coordZ = float(postprocset["TRACKER"]["frf"]["z"])
            hist_node = search_nodexyz(
                node_coordX, node_coordY, node_coordZ, coord, 1e-6
            )
            rstl = plotset["nodedof"] * (hist_node[0] - 1) - (
                plotset["nodedof"] - postprocset["TRACKER"]["frf"]["dof"] - 1
            )
            # val_Y = 20*np.log((abs(plotset["val_list"][rstl, :]))/10E-12)
            # val_Y = 20*np.log((abs(postporc_result['frf'][plotset["step"]]['val'][rstl, :]))/10E-12)
            val_Y = 20 * np.log(
                (
                    abs(
                        postporc_result["SOLUTION"][plotset["step"]]["VAL"][
                            hist_node[0] - 1, postprocset["TRACKER"]["frf"]["dof"]
                        ]
                    )
                )
                # / 10e-12
            )
            val_X = postporc_result["SOLUTION"][plotset["step"]]["FREQ"]
            xlabel = "FREQUENCY RESPONSE [Hz]"
            ylabel = "FREQUENCY RESPONSE [dB]"

        else:
            val_X = 0
            val_Y = 0
            xlabel = "erro"
            ylabel = "erro"
        # plot(val_X, val_Y, xlabel, ylabel, plotset["fignumb"])
        return val_X, val_Y, xlabel, ylabel
