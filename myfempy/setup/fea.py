from __future__ import annotations

import logging
import sys
from time import time

import numpy as np
import scipy.sparse as sp

from myfempy.core.utilities import setSteps
# from myfempy.core.solver import getSolver
from myfempy.io.controllers import (setElement, setGeometry, setMaterial,
                                    setMesh, setShape, setDomain, setCoupling,
                                    setPoints2NumericalIntegration)

from myfempy.plots.prevplot import preview_plot
from myfempy.setup.model import SetModel
from myfempy.setup.physics import SetPhysics
from myfempy.setup.results import setPostProcess
from myfempy.utils.utils import (clear_console, get_logo, get_version,
                                 loading_bar_v1, newDir, print_console)


class newAnalysis:
    """
    New analysis to FEA model simulation
    """

    def __init__(self, FEASolver) -> None:
        self.solver = FEASolver
        self.path = newDir("out")
        logging.basicConfig(
            filename=str(self.path) + "/" + "myfempy_log.log",
            encoding="utf-8",
            level=logging.DEBUG,
            filemode="w",
        )

    def Model(self, modeldata):
        """
        Model finite element model set

        Arguments:
            modeldata -- data information
        """
        clear_console()
        get_logo()
        print_console("mesh")
        try:
            modeldata["MESH"]["user_path"] = self.path
            Mesh = newAnalysis.__setMesh(modeldata)
            logging.info("TRY SET MESH -- SUCCESS")
        except:
            logging.warning("TRY SET MESH -- FAULT")
        try:
            Element = newAnalysis.__setElement(modeldata)
            logging.info("TRY SET ELEMENT -- SUCCESS")
        except:
            logging.warning("TRY SET ELEMENT -- FAULT")
        try:
            Shape = newAnalysis.__setShape(modeldata)
            logging.info("TRY SET SHAPE -- SUCCESS")
        except:
            logging.warning("TRY SET SHAPE -- FAULT")    
        try:
            Material = newAnalysis.__setMaterial(modeldata)
            logging.info("TRY SET MATERIAL -- SUCCESS")
        except:
            logging.warning("TRY SET MATERIAL -- FAULT")
        try:
            Geometry = newAnalysis.__setGeometry(modeldata)
            logging.info("TRY SET GEOMETRY -- SUCCESS")
        except:
            logging.warning("TRY SET GEOMETRY -- FAULT")
        try:
            GaussPoints = newAnalysis.__setIntGauss(modeldata)
            logging.info("TRY SET GAUSS POINTS -- SUCCESS")
        except:
            logging.warning("TRY SET GAUSS POINTS -- FAULT")            
        try:
            self.model = SetModel(Mesh, Element, Shape, Material, Geometry)
            self.model.modeldata = modeldata
            self.model.intgauss = GaussPoints
            logging.info("TRY SET FEMODEL -- SUCCESS")
        except:
            logging.warning("TRY SET FEMODEL -- FAULT")
            
        self.modelinfo = dict()
        self.modelinfo["inci"] = newAnalysis.getInci(self)
        self.modelinfo["coord"] = newAnalysis.getCoord(self)
        self.modelinfo["tabmat"] = newAnalysis.getTabmat(self)
        self.modelinfo["tabgeo"] = newAnalysis.getTabgeo(self)
        self.modelinfo["intgauss"] = GaussPoints
        
        try:
            self.modelinfo["regions"] = newAnalysis.getRegions(self)
        except:
            self.modelinfo["regions"] = []
        elem_set = self.model.element.getElementSet()
        self.modelinfo["tensor"] = len(elem_set["tensor"])
        self.modelinfo["dofs"] = elem_set["dofs"]
        self.modelinfo["nodedof"] = len(elem_set["dofs"]["d"])
        self.modelinfo["type_element"] = elem_set["key"]
        shape_set = self.model.shape.getShapeSet()
        self.modelinfo["shapeid"] = shape_set["id"]
        self.modelinfo["nodecon"] = len(shape_set["nodes"])
        self.modelinfo["elemdofs"] = len(shape_set["nodes"]) * self.modelinfo["nodedof"]
        self.modelinfo["type_shape"] = shape_set["key"]
        self.modelinfo["elemid"] = int(f'{elem_set["id"]}{shape_set["id"]}')
        self.modelinfo["nnode"] = len(self.model.coord)
        self.modelinfo["nelem"] = len(self.model.inci)
        self.modelinfo["fulldofs"] = len(elem_set["dofs"]["d"]) * len(self.model.coord)

    def Physic(self, physicdata):
        """
        Physic physics set load and boundary conditions

        Arguments:
            physicdata -- data information
        """
        print_console("pre")
        self.modelinfo["physic"] = physicdata["PHYSIC"]
        try:
            Loads, BoundCond = newAnalysis.__setDomain(physicdata)
            logging.info("TRY SET PHYSICS -- SUCCESS")
        except:
            logging.warning("TRY SET PHYSICS -- FAULT")

        self.physic = SetPhysics(self.model, Loads, BoundCond)
        # self.physic.physicdata = physicdata
        self.modelinfo["forces"] = newAnalysis.getLoadApply(self)

        constrains = newAnalysis.getBCApply(self)
        if any(constrains[:, 1] == 11):
            self.modelinfo["csleft"] = constrains[
                np.where(constrains[:, 1] == 11)[0], 0
            ]
            self.modelinfo["csright"] = constrains[
                np.where(constrains[:, 1] == 12)[0], 0
            ]
            self.modelinfo["constrains"] = constrains
        else:
            self.modelinfo["constrains"] = constrains

        if "COUPLING" in physicdata.keys():
            self.modelinfo["coupling"] = physicdata["COUPLING"]
            LoadCoup = newAnalysis.__setCoupling(physicdata)
            self.physic = SetPhysics(self.model, LoadCoup, BoundCond)
            self.modelinfo["forces"] = np.append(
                self.modelinfo["forces"], newAnalysis.getCouplingInterface(self), axis=0
            )

        # self.loadaply = FEANewAnalysis.getLoadApply(self)
        # self.constrains = FEANewAnalysis.getBCApply(self)

    def Assembly(self):
        """
        Assembly assembly of fe model algebric system

        Returns:
            matrix and vector from fe model
        """
        inci = self.modelinfo["inci"]
        coord = self.modelinfo["coord"]
        tabmat = self.modelinfo["tabmat"]
        tabgeo = self.modelinfo["tabgeo"]
        intgauss = self.modelinfo["intgauss"]
        # try:
        matrix = newAnalysis.getGlobalMatrix(
            self, inci, coord, tabmat, tabgeo, intgauss, self.symm, self.mp
        )
        #     logging.info("TRY RUN GLOBAL ASSEMBLY -- SUCCESS")
        # except:
        #     logging.warning("TRY RUN GLOBAL ASSEMBLY -- FAULT")
        loadaply = self.modelinfo["forces"]
        try:
            matrix = newAnalysis.getUpdateMatrix(self, matrix, loadaply)
            logging.info("TRY RUN UPDATE ASSEMBLY -- SUCCESS")
        except:
            logging.warning("TRY RUN UPDATE ASSEMBLY -- FAULT")
        try:
            forcelist = newAnalysis.getLoadArray(self, loadaply)
            logging.info("TRY RUN LOAD ASSEMBLY -- SUCCESS")
        except:
            logging.warning("TRY RUN LOAD ASSEMBLY -- FAULT")
        return matrix, forcelist

    def Solve(self, solverset):
        """
        runSolve run the solver set

        Arguments:
            solvedata -- data information

        Returns:
            solution
        """
        print_console("solver")
        # fulldofs = self.modelinfo['fulldofs']
        # self.modelinfo = dict()
        # self.modelinfo['coord'] = self.coord
        # self.modelinfo['regions'] = self.regions
        solverset["solverstatus"] = dict()
        # loading_bar_v1(5,"SOLVER")
        try:
            self.symm = solverset["SYMM"]
            if self.symm:
                solverset["solverstatus"]["typeasmb"] = "SYMMETRIC"
            else:
                solverset["solverstatus"]["typeasmb"] = "FULL"
        except:
            self.symm = False
            solverset["solverstatus"]["typeasmb"] = "FULL"
        try:
            self.mp = solverset["MP"]
            solverset["solverstatus"]["ncpu"] = (
                "PARALLEL_" + str(solverset["MP"]) + "_CORES"
            )
        except:
            self.mp = 0
            solverset["solverstatus"]["ncpu"] = "SERIAL_" + str(1) + "_CORE"
        # loading_bar_v1(10,"SOLVER")
        starttime = time()
        assembly, forcelist = newAnalysis.Assembly(self)
        endttime = time()
        solverset["solverstatus"]["timeasb"] = abs(endttime - starttime)
        solverset["solverstatus"]["memorysize"] = assembly["stiffness"].data.nbytes
        # loading_bar_v1(50,"SOLVER")
        constrains = self.modelinfo["constrains"]
        nsteps = setSteps(solverset["STEPSET"])
        constrainsdof = dict()
        try:
            freedof, fixedof, constdof = newAnalysis.getConstrains(self, constrains)
            logging.info("TRY RUN CONSTRAINS -- SUCCESS")
        except:
            logging.warning("TRY RUN CONSTRAINS -- FAULT")
        # loading_bar_v1(60,"SOLVER")
        constrainsdof["freedof"] = freedof
        constrainsdof["fixedof"] = fixedof
        constrainsdof["constdof"] = constdof
        try:
            Uc = newAnalysis.getDirichletNH(self, constrains)
            logging.info("TRY SET DNH CONSTRAINS -- SUCCESS")
        except:
            logging.warning("TRY SET DNH CONSTRAINS -- FAULT")
        if forcelist.shape[1] != nsteps:
            forcelist = np.repeat(forcelist, nsteps, axis=1)
        else:
            pass
        assembly["loads"] = forcelist
        if Uc.shape[1] != nsteps:
            Uc = np.repeat(Uc, nsteps, axis=1)
        else:
            pass
        assembly["bcdirnh"] = Uc
        # loading_bar_v1(80,"SOLVER")
        try:
            starttime = time()
            solverset["solution"] = self.solver.runSolve(
                assembly, constrainsdof, self.modelinfo, solverset
            )
            endttime = time()
            solverset["solverstatus"]["timesim"] = abs(endttime - starttime)
            logging.info("TRY RUN SOLVER -- SUCCESS")
        except:
            logging.warning("TRY RUN SOLVER -- FAULT")
        # loading_bar_v1(100,"SOLVER")
        return solverset

    def PreviewAnalysis(self, previewdata):
        """
        PreviewAnalysis preview the model+physic set

        Arguments:
            previewdata -- data information
        """
        # modelinfo = dict()
        # modelinfo["coord"] = self.coord   #FEANewAnalysis.getCoord(self)
        # modelinfo["inci"] = self.inci     #FEANewAnalysis.getInci(self)
        # modelinfo["tabgeo"] = self.tabgeo #FEANewAnalysis.getTabgeo(self)
        # shape_set = self.model.shape.getShapeSet()
        # elem_set = self.model.element.getElementSet()
        # self.modelinfo["nodecon"] = len(shape_set['nodes'])
        # self.modelinfo["elemid"] = elem_set['id']
        # try:
        #     self.regions = self.model.mesh.getRegionsList(self.model.mesh.getElementConection(self.model.modeldata['MESH']))
        # except:
        #     self.regions = []
        # modelinfo["regions"] = self.regions
        # try:
        #      modelinfo["forces"] = FEANewAnalysis.getLoadApply(self)
        # except:
        #     pass
        #     # modelinfo["forces"] = []
        # try:
        #     modelinfo["constrains"] = FEANewAnalysis.getBCApply(self)
        # except:
        #     pass
        #     modelinfo["constrains"] = []
        # modelinfo["forces"] = FEANewAnalysis.getLoadApply(self)
        # modelinfo["constrains"] = FEANewAnalysis.getBCApply(self)
        
        # try:
        preview_plot(previewdata, self.modelinfo, str(self.path))
        #     logging.info("TRY RUN PREVIEW PLOT -- SUCCESS")
        # except:
        #     logging.warning("TRY RUN PREVIEW PLOT -- FAULT")

    def PostProcess(self, postprocset):
        """
        PostProcess prost process the solution

        Arguments:
            postprocdata -- data information

        Returns:
            post process arrays
        """
        print_console("post")
        postprocdata = []

        self.modelinfo["elemvol"] = newAnalysis.getElementVolume(
            self,
            self.modelinfo["inci"],
            self.modelinfo["coord"],
            self.modelinfo["tabgeo"],
        )

        try:
            if "COMPUTER" in postprocset.keys():
                postprocdata = setPostProcess.getCompute(self, postprocset)
            if "TRACKER" in postprocset.keys():
                setPostProcess.getTracker(self, postprocset, postprocdata)
            if "OUTPUT" in postprocset.keys():
                postprocdata["log"] = []
                log_file = setPostProcess.getLog(self, postprocset, postprocdata)
                postprocdata["log"].append(log_file)
                logging.info("TRY GET POST PROCESS -- SUCCESS")
        except:
            logging.warning("TRY GET POST PROCESS -- FAULT")
        return postprocdata

    # gettings
    def getModelInfo(self):
        """
        getModelInfo _summary_

        Returns:
            _description_
        """
        return self.modelinfo

    def getInci(self):
        return self.model.getInci(self.model.modeldata)

    def getCoord(self):
        return self.model.getCoord(self.model.modeldata)

    def getTabmat(self):
        return self.model.getTabMat(self.model.modeldata)

    def getTabgeo(self):
        return self.model.getTabGeo(self.model.modeldata)

    def getIntGauss(self):
        return self.modelinfo["intgauss"]

    def getElementVolume(self, inci, coord, tabgeo):
        vol = np.zeros((inci.shape[0]))
        for ee in range(inci.shape[0]):
            vol[ee] = self.model.element.getElementVolume(
                self.model, inci, coord, tabgeo, ee
            )
        return vol

    def getElemStifLinearMat(
        self, inci, coord, tabmat, tabgeo, intgauss, element_number
    ):
        return self.model.element.getStifLinearMat(
            self.model, inci, coord, tabmat, tabgeo, intgauss, element_number
        )

    def getElemMassConsistentMat(
        self, inci, coord, tabmat, tabgeo, intgauss, element_number
    ):
        return self.model.element.getMassConsistentMat(
            self.model, inci, coord, tabmat, tabgeo, intgauss, element_number
        )

    def getGlobalMatrix(
        self, inci, coord, tabmat, tabgeo, intgauss, SYMM=None, MP=None
    ):
        return self.solver.getMatrixAssembler(
            self.model, inci, coord, tabmat, tabgeo, intgauss, SYMM=SYMM, MP=MP
        )

    def getForceList(self):
        return self.physic.getForceList(self.modelinfo["domain"])

    def getBoundCondList(self):
        return self.physic.getBoundCondList(self.modelinfo["domain"])

    def getLoadApply(self):
        return self.physic.getLoadApply(self.modelinfo)

    def getBCApply(self):
        return self.physic.getBoundCondApply(self.modelinfo)

    def getConstrains(self, constrains):
        nodetot = len(self.modelinfo["coord"])
        freedof, fixedof, constdof = self.solver.getConstrains(
            constrains, nodetot, self.modelinfo["nodedof"]
        )
        return freedof, fixedof, constdof

    def getDirichletNH(self, constrains):
        nodetot = len(self.modelinfo["coord"])
        return self.solver.getDirichletNH(
            constrains, nodetot, self.modelinfo["nodedof"]
        )

    def getLoadArray(self, loadaply):
        nodetot = len(self.modelinfo["coord"])
        return self.solver.getLoadAssembler(
            loadaply, nodetot, self.modelinfo["nodedof"]
        )

    def getCouplingInterface(self):
        return self.physic.getLoadCoup(self.modelinfo, self.modelinfo["coupling"])

    def getUpdateMatrix(self, matrix, addval):
        return self.physic.getUpdateMatrix(matrix, addval)

    def getRegions(self):
        return self.model.mesh.getRegionsList(
            self.model.mesh.getElementConection(self.model.modeldata["MESH"])
        )

    def getNodesFromRegions(self, set: int, type: str):
        if type == "point":
            domain_nodelist = newAnalysis.getRegions(self)[0][1][set - 1][1]
        elif type == "line":
            domain_nodelist = newAnalysis.getRegions(self)[1][1][set - 1][1]
        elif type == "plane":
            domain_nodelist = newAnalysis.getRegions(self)[2][1][set - 1][1]
        else:
            domain_nodelist = []
        return domain_nodelist

    def getElementFromNodesList(self, nodelist):
        return self.physic.getElementList(self.modelinfo["inci"], nodelist)
    

    # settings FEA ANALYSIS <privates>
    def __setMesh(modeldata):
        set_mesh = dict()
        set_mesh = modeldata["MESH"]
        set_mesh["SHAPE"] = modeldata["ELEMENT"]["SHAPE"]
        return setMesh(set_mesh)

    def __setShape(modeldata):
        return setShape(modeldata["ELEMENT"])

    def __setElement(modeldata):
        return setElement(modeldata["ELEMENT"])

    def __setMaterial(modeldata):
        return setMaterial(modeldata["MATERIAL"])

    def __setGeometry(modeldata):
        return setGeometry(modeldata["GEOMETRY"])

    def __setDomain(modeldata):
        return setDomain(modeldata["PHYSIC"])

    def __setCoupling(modeldata):
        return setCoupling(modeldata["COUPLING"])
    
    def __setIntGauss(modeldata):
        if "INTGAUSS" in modeldata["ELEMENT"].keys():
            intgauss = modeldata["ELEMENT"]["INTGAUSS"]
        else:
            intgauss = setPoints2NumericalIntegration(modeldata["ELEMENT"]["SHAPE"])
        return intgauss


    # def __setSolution(self, solvedata):
    #     # self.inci = FEANewAnalysis.getInci(self)
    #     # self.coord = FEANewAnalysis.getCoord(self)
    #     # self.tabmat = FEANewAnalysis.getTabmat(self)
    #     # self.tabgeo = FEANewAnalysis.getTabgeo(self)
    #     # self.intgauss = FEANewAnalysis.getIntGauss(self)
    #     # elem_set = self.model.element.getElementSet()
    #     # nodedof = len(elem_set["dofs"])
    #     self.solver.fulldofs = (self.modelinfo['dofs']) * len(self.modelinfo['nnode'])
    #     self.solver.solvedata = solvedata
