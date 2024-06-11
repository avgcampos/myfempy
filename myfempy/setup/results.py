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

    # def __init__(self) -> None:

    # self.domain = Model.domain
    # self.coord = Model.coord
    # self.inci = Model.inci
    # self.ntensor = Model.ntensor
    # self.nnode = Model.nnode
    # self.nelem = Model.nelem
    # self.dofe = Model.dofe
    # self.nodecon = Model.nodecon
    # self.elemid = Model.elemid
    # self.modelinfo
    def getCompute(self, postprocset):
        """_summary_

        Arguments:
            postprocset -- _description_

        Returns:
            _description_
        """
        # print_console("post")
        postporc_result = dict()
        postporc_result["solution"] = []
        postporc_result["stress_avr"] = []
        postporc_result["stress_elm"] = []
        postporc_result["strain_energy_elm"] = []

        postporc_result["balance"] = []
        # postporc_result["modes"] = []
        postporc_result["frf"] = []
        postporc_result["solverstatus"] = postprocset["SOLVERDATA"]["solverstatus"]

        if self.modelinfo["domain"] == "structural":
            SOLUTION = postprocset["SOLVERDATA"]["solution"]["U"]
            result_disp = np.zeros(
                (SOLUTION.shape[1], self.modelinfo["coord"].shape[0], 4)
            )
            for ns in range(SOLUTION.shape[1]):
                result_disp[ns, :, :] = setPostProcess.__displ(self, SOLUTION[:, ns])

        if "structural" in postprocset["COMPUTER"].keys():
            if "displ" in postprocset["COMPUTER"]["structural"].keys():
                for st in range(SOLUTION.shape[1]):
                    postporc_result["solution"].append(
                        {
                            "val": result_disp[st][:][:],
                            "title": "DEFORMATION",
                            "avr": True,
                        }
                    )

            if (
                "stress" in postprocset["COMPUTER"]["structural"].keys()
                and postprocset["COMPUTER"]["structural"]["stress"] == True
            ):
                result_stress = np.zeros(
                    (
                        SOLUTION.shape[1],
                        self.modelinfo["inci"].shape[0],
                        2 * self.modelinfo["tensor"] + 4,
                    )
                )

                for ns in range(SOLUTION.shape[1]):
                    result_stress[ns, :, :], title = setPostProcess.__stress(
                        self, SOLUTION[:, ns]
                    )

                # if "average" in postprocset["COMPUTER"]["structural"].keys() and postprocset["COMPUTER"]["structural"]["average"]==True:

                #     results_avr = np.zeros(
                #             (
                #             SOLUTION.shape[1],
                #             self.modelinfo['coord'].shape[0],
                #             2 * self.modelinfo['tensor'] + 3))

                #     for ns in range(SOLUTION.shape[1]):
                #         for nt in range(2 * self.modelinfo['tensor'] + 2):
                #             # results_avr[ns, :, nt] = PostComputer.results_average(self, result_stress[ns, :, nt])
                #             # results_avr[ns, :, nt] = results_average(self, result_stress[ns, :, nt])
                #             results_avr[ns, :, nt] = results_average(
                #                 result_stress[ns, :, nt],
                #                 self.modelinfo['nnode'],
                #                 self.modelinfo['nelem'],
                #                 self.modelinfo['nodedof'],
                #                 self.modelinfo['inci'],
                #             )

                # results_avr = results_avr

                for st in range(SOLUTION.shape[1]):
                    postporc_result["stress_elm"].append(
                        {"val": result_stress[st][:][:], "title": title, "avr": False}
                    )

                    # if "average" in postprocset["COMPUTER"]["structural"].keys() and postprocset["COMPUTER"]["structural"]["average"]==True:
                    #     postporc_result["stress_avr"].append(
                    #         {
                    #             "val": results_avr[st][:][:],
                    #             "title": title,
                    #             "avr": True,
                    #         }
                    #     )
            else:
                pass

        if "dynamic" in postprocset["COMPUTER"].keys():
            if "modes" in postprocset["COMPUTER"]["dynamic"].keys():
                FREQUENCY = postprocset["SOLVERDATA"]["solution"]["FREQ"]
                for st in range(SOLUTION.shape[1]):
                    postporc_result["solution"].append(
                        {
                            "val": result_disp[st][:][:],
                            "title": (
                                "MODE_"
                                + str(st + 1)
                                + "-"
                                + "FREQ_"
                                + str(np.round(FREQUENCY[st][2], 3))
                                + "Hz"
                            ),
                            "avr": True,
                        }
                    )

            if "frf" in postprocset["COMPUTER"]["dynamic"].keys():
                FREQUENCY = postprocset["SOLVERDATA"]["solution"]["FREQ"]
                for st in range(SOLUTION.shape[1]):
                    postporc_result["solution"].append(
                        {
                            "val": result_disp[st][:][:],
                            "title": (
                                "FREQUENCY_RESPONSE_"
                                + str(np.round(FREQUENCY[st], 3))
                                + "Hz"
                            ),
                            "freq": FREQUENCY[st],
                            "avr": True,
                        }
                    )
            else:
                pass

        else:
            pass

        setPostProcess.__tovtkplot(
            self, postprocset, postporc_result
        )  # save in vtk file

        return postporc_result

    def getTracker(self, postprocset, postporc_result):
        plotset = dict()
        hist_X = []
        hist_Y = []
        plotset["nodedof"] = self.modelinfo["nodedof"]

        max_steps = len(postporc_result["solution"])

        # if 'frf' in postprocset["TRACKER"].keys():
        #     max_steps = len(postporc_result["frf"])

        for st in range(max_steps):
            plotset["step"] = st
            # plotset["val_list"] = postporc_result['solution'][st]['val']
            # plotset["fignumb"] = 99
            # plotset["rstl"] = [0, modelinfo["ntensor"][0] + 1]
            val_X, val_Y, xlabel, ylabel = setPostProcess.__tracker_value(
                postporc_result, postprocset, plotset, self.modelinfo["coord"]
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
        write2log(filename, postprocset["OUTPUT"], self.modelinfo, postporc_result)

    def __tovtkplot(self, postprocset, postporc_result):
        """_summary_

        Arguments:
            postporc_result -- _description_
            postprocset -- _description_
        """
        # path = os.getcwd()
        plotdata = dict()
        plotdata["inci"] = self.modelinfo["inci"]
        plotdata["nodecon"] = self.modelinfo["nodecon"]
        plotdata["filename"] = (
            str(self.path) + "/" + postprocset["PLOTSET"]["filename"] + "_post_process"
        )
        plotdata["coord"] = self.modelinfo["coord"]

        plotdata["displ_POINT_DATA_val"] = []
        plotdata["displ_POINT_DATA_title"] = []

        plotdata["stress_POINT_DATA_val"] = []
        plotdata["stress_POINT_DATA_title"] = []

        plotdata["stress_CELL_DATA_val"] = []
        plotdata["stress_CELL_DATA_title"] = []

        # plotdata["strain_energy_CELL_DATA_val"] = []
        # plotdata["strain_energy_CELL_DATA_title"] = []

        plotdata["modes_POINT_DATA"] = []

        # ... any field in future
        # if "structural" in postprocset["COMPUTER"].keys():
        #     if "displ" in postprocset["COMPUTER"]["structural"].keys():
        if "structural" in postprocset["COMPUTER"].keys():
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
                ][:, 1:]
                plotdata["displ_POINT_DATA_title"] = postporc_result["solution"][st][
                    "title"
                ]

                if (
                    "stress" in postprocset["COMPUTER"]["structural"].keys()
                    and postprocset["COMPUTER"]["structural"]["stress"] == True
                ):
                    plotdata["stress_CELL_DATA_val"] = postporc_result["stress_elm"][
                        st
                    ]["val"]
                    plotdata["stress_CELL_DATA_title"] = postporc_result["stress_elm"][
                        st
                    ]["title"]

                if (
                    "average" in postprocset["COMPUTER"]["structural"].keys()
                    and postprocset["COMPUTER"]["structural"]["average"] == True
                ):
                    plotdata["stress_POINT_DATA_val"] = postporc_result["stress_avr"][
                        st
                    ]["val"]
                    plotdata["stress_POINT_DATA_title"] = postporc_result["stress_avr"][
                        st
                    ]["title"]

                # if "strain_energy" in postprocset["COMPUTER"]["structural"].keys():
                #     plotdata["strain_energy_CELL_DATA_val"] = postporc_result["strain_energy_elm"][st]["val"]
                #     plotdata["strain_energy_CELL_DATA_title"] = postporc_result["strain_energy_elm"][st]["title"]

                convert_to_vtk(plotdata)

        elif "dynamic" in postprocset["COMPUTER"].keys():
            # plotdata["filename"] = (path + "/" + postprocset["PLOTSET"]["filename"] + "_results_modes")
            plotdata["modes_POINT_DATA"] = postporc_result["solution"]
            # plotdata["modes_POINT_DATA_title"] = postporc_result["solution"][st]["title"]

        convert_to_vtk(plotdata)

    def __displ(self, U):
        """_summary_

        Arguments:
            U -- _description_

        Returns:
            _description_
        """
        return self.model.element.getElementDeformation(U, self.modelinfo)

    def __stress(self, U):
        """_summary_

        Arguments:
            U -- _description_

        Returns:
            _description_
        """
        stress_list = np.zeros(
            (self.modelinfo["nelem"], self.modelinfo["tensor"] + 1), dtype=float
        )
        strain_list = np.zeros(
            (self.modelinfo["nelem"], self.modelinfo["tensor"] + 1), dtype=float
        )
        compliance_list = np.zeros((self.modelinfo["nelem"], 1), dtype=float)
        factor_of_safety = np.zeros((self.modelinfo["nelem"], 1), dtype=float)

        pt, wt = gauss_points(self.modelinfo["type_shape"], 1)
        for ee in range(self.modelinfo["nelem"]):
            # tensor = get_stress(self.modelinfo, U, ee)
            epsilon, strain = self.model.material.getElementStrain(
                self.model, U, pt, ee
            )
            sigma, stress = self.model.material.getElementStress(
                self.model, epsilon, ee
            )
            strain_energy = self.model.material.getStrainEnergyDensity(
                sigma, epsilon, self.modelinfo["elemvol"][ee]
            )
            FoS = 0.0
            # FoS = self.model.material.getFailureCriteria(...)
            stress_list[ee, :] = stress
            strain_list[ee, :] = strain
            compliance_list[ee, :] = strain_energy
            factor_of_safety[ee, :] = FoS
        result = np.concatenate(
            (stress_list, stress_list, compliance_list, factor_of_safety), axis=1
        )
        tistrs = self.model.material.getTitleStress()
        tistrn = self.model.material.getTitleStrain()
        ticomp = self.model.material.getTitleCompliance()
        tifos = ["FoS_YIELD_VON_MISES"]  # self.model.material.getTitleFoS()
        title = np.concatenate((tistrs, tistrn, ticomp, tifos), axis=0)
        return result, title

    # def __strain_energy(self, sigma, epsilon, elemvol):
    #     result =  self.model.material.getStrainEnergyDensity(sigma, epsilon, elemvol)
    #     return result

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
            val_Y = postporc_result["solution"][plotset["step"]]["val"][
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
                        postporc_result["solution"][plotset["step"]]["val"][
                            hist_node[0] - 1, postprocset["TRACKER"]["frf"]["dof"]
                        ]
                    )
                )
                / 10e-12
            )
            val_X = postporc_result["solution"][plotset["step"]]["freq"]
            xlabel = "FREQUENCY RESPONSE [Hz]"
            ylabel = "FREQUENCY RESPONSE [dB]"

        else:
            val_X = 0
            val_Y = 0
            xlabel = "erro"
            ylabel = "erro"
        # plot(val_X, val_Y, xlabel, ylabel, plotset["fignumb"])
        return val_X, val_Y, xlabel, ylabel
