#!/usr/bin/env python
from myfempy.utils.utils import print_console
from myfempy.io.iovtk import convert_to_vtk
from myfempy.postprc.postset import get_stress, get_displ
from myfempy.felib.felemset import get_elemset
from myfempy.postprc.average import results_average
import scipy.sparse as sp
import os
import numpy as np

__doc__ = """
Post-Process Calculator
"""


class PostProcess:
    """_summary_"""

    def __init__(self, modelinfo):
        self.modelinfo = modelinfo
        self.dofe = modelinfo["nodecon"][0] * modelinfo["nodedof"][0]
        self.nodecon = modelinfo["nodecon"][0]
        self.fulldof = modelinfo["nodedof"][0] * len(modelinfo["coord"])
        self.nodedof = modelinfo["nodedof"][0]
        self.nelem = len(modelinfo["inci"])
        self.nnode = len(modelinfo["coord"])
        self.inci = modelinfo["inci"]
        self.coord = modelinfo["coord"]
        self.tabmat = modelinfo["tabmat"]
        self.tabgeo = modelinfo["tabgeo"]
        self.ntensor = modelinfo["ntensor"][0]

    def displ(self, U):
        """_summary_

        Arguments:
            U -- _description_

        Returns:
            _description_
        """
        result = get_displ(self.modelinfo, U)
        return result

    def stress(self, U):
        """_summary_

        Arguments:
            U -- _description_

        Returns:
            _description_
        """
        stress_list = np.zeros(
            (self.nelem, self.modelinfo["ntensor"][0] + 1), dtype=float
        )
        strain_list = np.zeros(
            (self.nelem, self.modelinfo["ntensor"][0] + 1), dtype=float
        )
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
        """_summary_

        Arguments:
            U -- _description_

        Returns:
            _description_
        """
        element = get_elemset(int(self.modelinfo["elemid"][0]))
        setelement = element(self.modelinfo)
        lines = [[None] * 2]
        if len(self.modelinfo["regions"]) == 0:
            for ll in range(len(self.modelinfo["inci"])):
                lines.append(
                    [
                        self.modelinfo["inci"][ll, 0],
                        self.modelinfo["inci"][ll, [4, 5]].tolist(),
                    ]
                )
            lines = lines[1::][::]
        else:
            lines = self.modelinfo["regions"][1][1]
        ifb, title = setelement.intforces(U, lines)
        return ifb, title

    def save_vtk(self, postporc_result, postprocset):
        """_summary_

        Arguments:
            postporc_result -- _description_
            postprocset -- _description_
        """
        path = os.getcwd()
        plotdata = dict()
        plotdata["inci"] = self.modelinfo["inci"]
        plotdata["nodecon"] = self.modelinfo["nodecon"]
        plotdata["elemid"] = self.modelinfo["elemid"]
        plotdata["filename"] = (
            path + "/" + postprocset["PLOTSET"]["filename"] + "_results_step-" + str(0)
        )
        plotdata["coord"] = self.modelinfo["coord"]
        plotdata["displ_POINT_DATA_val"] = []
        plotdata["displ_POINT_DATA_name"] = []
        plotdata["displ_POINT_DATA_title"] = []
        plotdata["stress_POINT_DATA_val"] = []
        plotdata["stress_POINT_DATA_name"] = []
        plotdata["stress_POINT_DATA_title"] = []
        plotdata["stress_CELL_DATA_val"] = []
        plotdata["stress_CELL_DATA_name"] = []
        plotdata["stress_CELL_DATA_title"] = []
        plotdata["modes_POINT_DATA"] = []
        # ... any field in future
        if "elasticity" in postprocset["COMPUTER"].keys():
            for st in range(len(postporc_result["displ"])):
                if postprocset["COMPUTER"]["elasticity"]["displ"]:
                    plotdata["filename"] = (
                        path
                        + "/"
                        + postprocset["PLOTSET"]["filename"]
                        + "_results_step-"
                        + str(st + 1)
                    )
                    plotdata["displ_POINT_DATA_val"] = postporc_result["displ"][st][
                        "val"
                    ][:, 1:]
                    plotdata["displ_POINT_DATA_title"] = postporc_result["displ"][st][
                        "title"
                    ]
                if postprocset["COMPUTER"]["elasticity"]["stress"]:
                    plotdata["stress_CELL_DATA_val"] = postporc_result["stress_elm"][
                        st
                    ]["val"]
                    plotdata["stress_CELL_DATA_title"] = postporc_result["stress_elm"][
                        st
                    ]["title"]
                if postprocset["COMPUTER"]["elasticity"]["average"]:
                    plotdata["stress_POINT_DATA_val"] = postporc_result["stress_avr"][
                        st
                    ]["val"]
                    plotdata["stress_POINT_DATA_title"] = postporc_result["stress_avr"][
                        st
                    ]["title"]
                convert_to_vtk(plotdata)
        if "eig" in postprocset["COMPUTER"].keys():
            if postprocset["COMPUTER"]["eig"]["modes"]:
                plotdata["filename"] = (
                    path + "/" + postprocset["PLOTSET"]["filename"] + "_results_modes"
                )
                plotdata["modes_POINT_DATA"] = postporc_result["modes"]
            convert_to_vtk(plotdata)

    # @profile
    def compute(self, postprocset):
        """_summary_

        Arguments:
            postprocset -- _description_

        Returns:
            _description_
        """
        # print_console("post")
        postporc_result = dict()
        postporc_result["displ"] = []
        postporc_result["stress_avr"] = []
        postporc_result["stress_elm"] = []
        postporc_result["intforces_plot"] = []
        postporc_result["modes"] = []
        postporc_result["frf"] = []

        if self.modelinfo["domain"] == "structural":
            SOLUTION = postprocset["SOLUTION"]

        result_disp = np.zeros(
            (SOLUTION["U"].shape[1], self.modelinfo["coord"].shape[0], 4)
        )
        for ns in range(SOLUTION["U"].shape[1]):
            result_disp[ns, :, :] = PostProcess.displ(self, SOLUTION["U"][:, ns])
        if "elasticity" in postprocset["COMPUTER"].keys():
            if postprocset["COMPUTER"]["elasticity"]["displ"]:
                for st in range(SOLUTION["U"].shape[1]):
                    postporc_result["displ"].append(
                        {
                            "val": result_disp[st][:][:],
                            "title": "DISPLACEMENT",
                            "avr": True,
                        }
                    )
            if postprocset["COMPUTER"]["elasticity"]["stress"]:
                result_stress = np.zeros(
                    (
                        SOLUTION["U"].shape[1],
                        self.modelinfo["inci"].shape[0],
                        2 * self.modelinfo["ntensor"][0] + 2,
                    )
                )
                for ns in range(SOLUTION["U"].shape[1]):
                    result_stress[ns, :, :], title = PostProcess.stress(
                        self, SOLUTION["U"][:, ns]
                    )
                results_avr = np.zeros(
                    (
                        SOLUTION["U"].shape[1],
                        self.modelinfo["coord"].shape[0],
                        2 * self.modelinfo["ntensor"][0] + 2,
                    )
                )
                for ns in range(SOLUTION["U"].shape[1]):
                    for nt in range(2 * self.modelinfo["ntensor"][0] + 2):
                        # results_avr[ns, :, nt] = PostComputer.results_average(self, result_stress[ns, :, nt])
                        # results_avr[ns, :, nt] = results_average(self, result_stress[ns, :, nt])
                        results_avr[ns, :, nt] = results_average(
                            result_stress[ns, :, nt],
                            self.nnode,
                            self.nelem,
                            self.dofe,
                            self.inci,
                        )

                result_stress_avr = results_avr
                for st in range(SOLUTION["U"].shape[1]):
                    postporc_result["stress_elm"].append(
                        {"val": result_stress[st][:][:], "title": title, "avr": False}
                    )
                    postporc_result["stress_avr"].append(
                        {
                            "val": result_stress_avr[st][:][:],
                            "title": title,
                            "avr": True,
                        }
                    )
            else:
                pass
        if "balance" in postprocset["COMPUTER"].keys():
            if postprocset["COMPUTER"]["balance"]["intforces_plot"]:
                for ns in range(SOLUTION["U"].shape[1]):
                    ifb, title = PostProcess.intforces(self, SOLUTION["U"][:, ns])
                    postporc_result["intforces_plot"].append(
                        {"val": [ifb["le"], ifb["val"]], "title": title}
                    )
            else:
                pass
        if "eig" in postprocset["COMPUTER"].keys():
            if postprocset["COMPUTER"]["eig"]["modes"]:
                for st in range(SOLUTION["U"].shape[1]):
                    postporc_result["modes"].append(
                        {
                            "val": result_disp[st][:][:],
                            "title": (
                                "MODE_"
                                + str(st + 1)
                                + "-"
                                + "FREQ_"
                                + str(round(SOLUTION["FREQ"][st][2], 2))
                                + "Hz"
                            ),
                            "avr": True,
                        }
                    )
            else:
                pass
        if "frf" in postprocset["COMPUTER"].keys():
            if postprocset["COMPUTER"]["frf"]["frf_plot"]:
                postporc_result["frf"].append(
                    {"val": SOLUTION["U"], "freqlog": SOLUTION["FREQ"]}
                )
            else:
                pass
        else:
            pass
        PostProcess.save_vtk(self, postporc_result, postprocset)  # save in vtk file
        return postporc_result
