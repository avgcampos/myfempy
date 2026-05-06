from __future__ import annotations

import logging
import sys
from time import time

import numpy as np
import numpy.typing as npt

from myfempy.core.utilities import setSteps
# from myfempy.core.solver import getSolver
from myfempy.io.controllers import (setElement, setGeometry, setMaterial,
                                    setMesh, setShape, setDomain, setCoupling,
                                    setPoints2NumericalIntegration)

from myfempy.plots.prevplot import preview_plot
from myfempy.api.model import SetModel
from myfempy.api.physics import SetPhysics
from myfempy.api.results import setPostProcess
from myfempy.utils.utils import (clear_console, get_logo, get_version,
                                 loading_bar_v1, newDir, print_console, get_about)

__docformat__ = "google"

__doc__ = """
API module for performing finite element analysis with the myfempy package.

------------------------------------------------------------------------
                                        __                                
                     _ __ ___   _   _  / _|  ___  _ __ ___   _ __   _   _ 
                    | '_ ` _ \ | | | || |_  / _ \| '_ ` _ \ | '_ \ | | | |
                    | | | | | || |_| ||  _||  __/| | | | | || |_) || |_| |
                    |_| |_| |_| \__, ||_|   \___||_| |_| |_|| .__/  \__, |
                                |___/                       |_|     |___/ 

                    myfempy -- MultiphYsics Finite Element Module to PYthon    
                                COMPUTATIONAL ANALYSIS PROGRAM                   
                    Copyright (C) 2022-2026 Antonio Vinicius Garcia Campos

------------------------------------------------------------------------

This Python file is part of myfempy project.

myfempy is a python package based on finite element method to multiphysics
analysis. The code is open source and intended for educational and scientific
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

class newAnalysis:
    """
    Setup the New Analysis to FEA simulation
    """
    def __init__(self, FEASolver: object, path: str = None):
        """Initializes a new Finite Element Analysis (FEA) project environment.

        Sets up the solution directories, default logging configuration, and binds the
        numerical solver module to the simulation instance.

        Args:
            FEASolver: Class or module responsible for solving the physical state equations 
                (e.g., StaticLinear, SteadyStateLinear).
            path: Target directory path for exporting simulation logs and output files. 
                If None is passed, defaults to creating an "out" folder.

        Example:
            >>> from myfempy.api.api import newAnalysis
            >>> from myfempy.core.solver import SteadyStateLinear
            >>> FEA = newAnalysis(FEASolver=SteadyStateLinear, path="simulation_results")
        """
        try:
            self.solver = FEASolver
            self.path = newDir(path)
            text_init = "TRY SET NEW ANALYSIS AND SOLVER -- SUCCESS"
        except:
            self.path = newDir("out")
            text_init = "TRY SET NEW ANALYSIS AND SOLVER -- FAULT"
            
        logging.basicConfig(
            filename=str(self.path) + "/" + "myfempy_api-log.log",
            encoding="utf-8",
            level=logging.DEBUG,
            filemode="w",
        )
        logging.info(text_init)

    def Model(self, modeldata: dict) -> None:
        """Sets up the mesh, element configuration, materials, and domain geometry.

        Parses the unified configuration dictionary to build the mathematical model representation,
        computes element volumes, and compiles internal metadata arrays (such as coordinates
        and connectivity).

        Args:
            modeldata: A structured configuration dictionary containing modeling setup:
                - 'MESH': Coordinates and element incidence.
                - 'ELEMENT': Finite element formulation (e.g., 'structplane'), shape, and integration order.
                - 'MATERIAL': Constitutive model and properties.
                - 'GEOMETRY': Thickness, cross-section areas, or generic geometric descriptors.

        Example:
            >>> model_config = {
            ...     'MESH': {
            ...         'TYPE': 'manual',
            ...         'COORD': [[1, 0., 0.], [2, 1., 0.], [3, 1., 1.], [4, 0., 1.]],
            ...         'INCI': [[1, 1, 1, 1, 2, 3, 4]]
            ...     },
            ...     'ELEMENT': {
            ...         'TYPE': 'structplane',
            ...         'SHAPE': 'quad4',
            ...         'INTGAUSS': 2
            ...     },
            ...     'MATERIAL': {
            ...         'MAT': 'planestress',
            ...         'TYPE': 'isotropic',
            ...         'PROPMAT': [{'NAME': 'Steel', 'EXX': 210e9, 'VXY': 0.3}]
            ...     },
            ...     'GEOMETRY': {
            ...         'GEO': 'thickness',
            ...         'PROPGEO': [{'NAME': 'Pl_1', 'THICKN': 0.01}]
            ...     }
            ... }
            >>> FEA.Model(model_config)
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
            # self.model.intgauss = GaussPoints
            logging.info("TRY SET FEMODEL -- SUCCESS")
        except:
            logging.warning("TRY SET FEMODEL -- FAULT")
            
        self.model.inci = self.model.getInci(self.model.modeldata)
        self.model.coord = self.model.getCoord(self.model.modeldata)
        self.model.tabmat = self.model.getTabMat(self.model.modeldata)
        self.model.tabgeo = self.model.getTabGeo(self.model.modeldata)
        self.model.intgauss = GaussPoints

        self.model.modelinfo = dict()
        try:
            self.model.regions = self.model.mesh.getRegionsList(
            self.model.mesh.getElementConection(self.model.modeldata["MESH"])
        )
        except:
            self.model.regions = []
        elem_set = self.model.element.getElementSet()
        self.model.modelinfo["tensor"] = len(elem_set["tensor"])
        self.model.modelinfo["dofs"] = elem_set["dofs"]
        self.model.modelinfo["nodedof"] = len(elem_set["dofs"]["d"])
        self.model.modelinfo["type_element"] = elem_set["key"]
        shape_set = self.model.shape.getShapeSet()
        self.model.modelinfo["shapeid"] = shape_set["id"]
        self.model.modelinfo["nodecon"] = len(shape_set["nodes"])
        self.model.modelinfo["elemdofs"] = len(shape_set["nodes"]) * self.model.modelinfo["nodedof"]
        self.model.modelinfo["type_shape"] = shape_set["key"]
        self.model.modelinfo["elemid"] = int(f'{elem_set["id"]}{shape_set["id"]}')
        self.model.modelinfo["nnode"] = len(self.model.coord)
        self.model.modelinfo["nelem"] = len(self.model.inci)
        self.model.modelinfo["fulldofs"] = len(elem_set["dofs"]["d"]) * len(self.model.coord)

        self.model.elemvol = newAnalysis.getElementVolume(
            self,
            self.model.inci,
            self.model.coord,
            self.model.tabgeo,
        )

    def Physic(self, physicdata: dict) -> None:
        """Configures loads, coupling fields, and boundary constraints.

        Initializes force vectors, multiphysics couplings, and kinematic boundary
        constraints (such as essential Dirichlet conditions) on specified coordinates or nodal lists.

        Args:
            physicdata: Configuration dictionary specifying boundary conditions and loads.
                - 'PHYSIC': Holds physical 'DOMAIN' type, 'LOAD' array, and 'BOUNDCOND' array.
                - 'COUPLING' (optional): Holds interface multi-field configuration.

        Example:
            >>> physics_config = {
            ...     'PHYSIC': {
            ...         'DOMAIN': 'structural',
            ...         'LOAD': [
            ...             {'TYPE': 'force', 'DOF': 'uy', 'DIR': 'node', 'LOC': {'x': 1., 'y': 1.}, 'VAL': [-500.0]}
            ...         ],
            ...         'BOUNDCOND': [
            ...             {'TYPE': 'fixed', 'DOF': 'full', 'DIR': 'node', 'LOC': {'x': 0., 'y': 0.}}
            ...         ]
            ...     }
            ... }
            >>> FEA.Physic(physics_config)
        """
        print_console("phy")
        # self.model.modelinfo["physic"] = physicdata["PHYSIC"]
        try:
            Loads, BoundCond = newAnalysis.__setDomain(physicdata)
            self.physic = SetPhysics(self.model, Loads, BoundCond)
            self.physic.physicdata = physicdata
            logging.info("TRY SET PHYSICS -- SUCCESS")
        except:
            self.physic = []
            logging.warning("TRY SET PHYSICS -- FAULT")

        try:
            self.physic.forces = newAnalysis.getLoadApply(self)
            logging.info("TRY SET PHYSICS.FORCES -- SUCCESS")
        except:
            self.physic.forces = []
            logging.warning("TRY SET PHYSICS.FORCES -- FAULT")

        if "COUPLING" in physicdata.keys():
            coupling_load_zero = self.physic.forces
            LoadCoup, BoundCond = newAnalysis.__setCoupling(physicdata)
            self.physic = SetPhysics(self.model, LoadCoup, BoundCond)
            self.physic.physicdata = physicdata
            self.physic.forces = np.append(
                coupling_load_zero, newAnalysis.getCouplingInterface(self), axis=0
            )
            
        try:
            constrains = newAnalysis.getBCApply(self)
            if any(constrains[:, 1] == 11):
                self.physic.csleft = constrains[
                    np.where(constrains[:, 1] == 11)[0], 0
                ]
                self.physic.csright = constrains[
                    np.where(constrains[:, 1] == 12)[0], 0
                ]
                self.physic.constrains = constrains
            else:
                self.physic.constrains = constrains
            logging.info("TRY SET PHYSICS.CONSTRAINS -- SUCCESS")
        except:
            self.physic.constrains = []
            logging.warning("TRY SET PHYSICS.CONSTRAINS -- FAULT")

        # self.loadaply = FEANewAnalysis.getLoadApply(self)
        # self.constrains = FEANewAnalysis.getBCApply(self)

    def Assembly(self, Model: object) -> npt.NDArray[np.float64]:
        """Assembles the element equations into the global algebraic system of equations.

        Computes the global linearized matrix structure (e.g., global stiffness matrix) 
        and updates it incorporating force-term arrays under the current simulation model.

        Args:
            Model: Active Model instance containing geometry, material, and integration data.

        Returns:
            A tuple containing:
                - matrix: The assembled global matrix array (typically stiffness/governing behavior).
                - forcelist: The fully-assembled global force load vector list.

        Example:
            >>> K_global, F_global = FEA.Assembly(Model=FEA.model)
        """
        try:
            matrix = newAnalysis.getGlobalMatrix(self, Model, self.model.inci, self.model.coord, self.model.tabmat, self.model.tabgeo, self.model.intgauss, self.symm, self.mp)
            logging.info("TRY RUN GLOBAL ASSEMBLY -- SUCCESS")
        except:
            logging.warning("TRY RUN GLOBAL ASSEMBLY -- FAULT")
        loadaply = self.physic.forces
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
    # 
    def Solve(self, solverset=None) -> dict:
        """Executes the finite element equations solver over designated steps.

        Compiles matrix system transformations, processes Dirichlet non-homogeneous values,
        assembles loads across time/loading steps, and calls the bound numerical analysis core.

        Args:
            solverset: Configuration parameters dictionary for the numeric solver:
                - 'STEPSET': Steps mapping dictionary (type, start, end, step).
                - 'SYMM' (bool): Enable symmetric matrix optimizations.
                - 'MP' (bool/int): Enable multiprocessing CPU core count.

        Returns:
            A modified 'solverset' dictionary containing simulation status, calculation logs,
            times elapsed, and the final fields inside 'solverset["solution"]'.

        Example:
            >>> run_config = {
            ...     'STEPSET': {'type': 'table', 'start': 0.0, 'end': 1.0, 'step': 1.0},
            ...     'SYMM': False,
            ...     'MP': False
            ... }
            >>> results = FEA.Solve(run_config)
            >>> print(results["solution"])
        """
        print_console("solver")
        print(">>> RUNNING SOLVER:")
        print(self.solver.__doc__)
        try:
            solverset["solverstatus"] = dict()
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
        assembly, forcelist = newAnalysis.Assembly(self, Model=self.model)
        endttime = time()
        solverset["solverstatus"]["timeasb"] = abs(endttime - starttime)
        solverset["solverstatus"]["memorysize"] = (assembly["stiffness"].todense().nbytes)/1e6
        # loading_bar_v1(50,"SOLVER")
        try:
            constrains = self.physic.constrains
            freedof, fixedof, constdof = newAnalysis.getConstrains(self, constrains)
            logging.info("TRY RUN CONSTRAINS -- SUCCESS")
        except:
            logging.warning("TRY RUN CONSTRAINS -- FAULT")
        # loading_bar_v1(60,"SOLVER")
        constrainsdof = dict()
        constrainsdof["freedof"] = freedof
        constrainsdof["fixedof"] = fixedof
        constrainsdof["constdof"] = constdof
        try:
            Uc = newAnalysis.getDirichletNH(self, constrains)
            logging.info("TRY SET DNH CONSTRAINS -- SUCCESS")
        except:
            logging.warning("TRY SET DNH CONSTRAINS -- FAULT")
        nsteps = setSteps(solverset["STEPSET"])
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
            solverset["solution"] = self.solver.runSolve(self.model, self.physic, assembly, constrainsdof, solverset)
            endttime = time()
            solverset["solverstatus"]["timesim"] = abs(endttime - starttime)
            logging.info("TRY RUN SOLVER -- SUCCESS")
        except:
            logging.warning("TRY RUN SOLVER -- FAULT")
        # loading_bar_v1(100,"SOLVER")
        print_console("succ")
        return solverset

    def PreviewAnalysis(self, previewdata) -> None:
        """Renders pre-simulation plots for physical inspection of modeling items.

        Draws geometry shapes, elements, nodes, and applied load vectors before 
        running the solver.

        Args:
            previewdata: Graphical visualization configuration containing rendering attributes
                (e.g., 'RENDER' with keys like 'show', 'scale', 'lines', and 'savepng').

        Example:
            >>> preview_config = {
            ...     'RENDER': {
            ...         'filename': 'structural_preview',
            ...         'show': True,
            ...         'scale': 2.5,
            ...         'savepng': True,
            ...         'lines': True
            ...     }
            ... }
            >>> FEA.PreviewAnalysis(preview_config)
        """          
        try:
            preview_plot(self.model, previewdata, str(self.path), self.physic)
            logging.info("TRY RUN PREVIEW PLOT -- SUCCESS")
        except:
            preview_plot(self.model, previewdata, str(self.path))
            logging.warning("TRY RUN PREVIEW PLOT -- FAULT")

    def PostProcess(self, postprocset) -> dict:
        """Processes solutions, computes auxiliary fields, and builds output reports.

        Translates primary state parameters into derivative properties (such as strain, 
        stress tensor, and flux) and creates visualization plots and text records.

        Args:
            postprocset: Directives dictionary defining post-processing metrics:
                - "SOLVERDATA": Output simulation data containing solutions.
                - "COMPUTER": Physical fields to evaluate (e.g. displacement, stress).
                - "PLOTSET" or "PLOT": Map detailing export formats and CSV exports.
                - "REPORT": Config containing boolean flags for requested output tables.

        Returns:
            A post-processed analysis dataset dictionary mapping physical variables, plots, 
            and log file paths.

        Example:
            >>> config = {
            ...     "SOLVERDATA": results,
            ...     "COMPUTER": {'structural': {'displ': True, 'stress': True}},
            ...     "PLOTSET": {'show': True, 'filename': 'results_stress', 'savepng': True},
            ...     "REPORT": {'log': True, 'get': {'nelem': True, 'coord': True}}
            ... }
            >>> data = FEA.PostProcess(config)
        """
        print_console("post")
        postprocdata = []
        try:
            if "COMPUTER" in postprocset.keys():
                postprocdata = setPostProcess.getCompute(self, postprocset)
            if "PLOT" in postprocset.keys():
                setPostProcess.getPlotCSV(self, postprocset, postprocdata)
            if "REPORT" in postprocset.keys():
                postprocdata["log"] = []
                log_file = setPostProcess.getLog(self, postprocset, postprocdata)
                postprocdata["log"].append(log_file)
            logging.info("TRY GET POST PROCESS -- SUCCESS")
        except:
                logging.warning("TRY GET POST PROCESS -- FAULT")
        print_console("thank")
        return postprocdata

    # GET MODEL
    def getModel(self) -> object:
        """Retrieves the active finite element model container object.

        Returns:
            The SetModel instance describing mesh, connectivity, and formulations.

        Example:
            >>> model_obj = FEA.getModel()
        """
        return self.model

    def getModelInfo(self) -> dict:
        """Retrieves summary attributes and counts of degrees of freedom from the model.

        Returns:
            A dictionary containing structural details like elements count ('nelem'),
            nodes count ('nnode'), and total degrees of freedom ('fulldofs').

        Example:
            >>> info = FEA.getModelInfo()
            >>> print(info['nelem'], info['fulldofs'])
        """
        return self.model.modelinfo

    def getInci(self) -> npt.NDArray[np.float64]:
        """Retrieves the element incidence and connectivity matrix.

        Returns:
            An array where each row lists element tags, type identifiers, material IDs,
            and bounding node indices.

        Example:
            >>> connectivity = FEA.getInci()
        """
        return self.model.inci

    def getCoord(self) -> npt.NDArray[np.float64]:
        """Retrieves the global spatial coordinates table of all mesh nodes.

        Returns:
            An array of dimensions (N_nodes, 4) mapping node tags to [x, y, z] coordinates.

        Example:
            >>> nodes_coords = FEA.getCoord()
        """
        return self.model.coord

    def getTabmat(self) -> list:
        """Retrieves the material properties configuration table.

        Returns:
            A list of dictionary objects matching active materials.

        Example:
            >>> materials = FEA.getTabmat()
        """
        return self.model.tabmat
    def getTabgeo(self) -> list:
        """Retrieves the geometric cross-section/thickness attributes table.

        Returns:
            A list containing properties for geometric profiles.

        Example:
            >>> thickness_props = FEA.getTabgeo()
        """
        return self.model.tabgeo

    def getIntGauss(self) -> int:
        """Retrieves the number of Gauss integration points for numerical integration.

        Returns:
            The integration order/points parameter.

        Example:
            >>> order = FEA.getIntGauss()
        """
        return self.model.intgauss
    
    def getRegions(self) -> list:
        """Retrieves mesh entity group regions imported from GMSH files.

        Returns:
            A nested list of elements grouped by geometric entity tag definitions.

        Example:
            >>> physical_groups = FEA.getRegions()
        """
        return self.model.regions

    # MELHORAR ESTES COMANDOS
    def getElementVolume(self, inci:npt.NDArray[np.float64], coord:npt.NDArray[np.float64], tabgeo:list) -> npt.NDArray[np.float64]:
        """Calculates structural volumes (or areas/lengths) for all mesh elements.

        Calculates individual volumes/volumes-integrals for 1D, 2D, or 3D element profiles.

        Args:
            inci: Nodal incidence array of the elements.
            coord: Nodal spatial coordinate array.
            tabgeo: Geometry properties collection.

        Returns:
            A 1D numpy array containing the computed volume/area value for each element.

        Example:
            >>> volumes = FEA.getElementVolume(FEA.model.inci, FEA.model.coord, FEA.model.tabgeo)
        """
        vol = np.zeros((inci.shape[0]))
        for ee in range(inci.shape[0]):
            vol[ee] = self.model.element.getElementVolume(
                self.model, inci, coord, tabgeo, ee
            )
        return vol

    def getElemStifLinearMat(
        self, inci: npt.NDArray[np.float64], coord: npt.NDArray[np.float64], tabmat: list, tabgeo: list, intgauss: int, element_number: int
    ) -> npt.NDArray[np.float64]:
        """Computes the element linear stiffness matrix.

        Args:
            inci: Nodal incidence matrix.
            coord: Global coordinate coordinates.
            tabmat: Material configuration properties.
            tabgeo: Geometry properties profile.
            intgauss: Order count of integration.
            element_number: Index identifier of the element.

        Returns:
            A 2D array representing the local stiffness matrix.

        Example:
            >>> k_local = FEA.getElemStifLinearMat(
            ...     FEA.getInci(), FEA.getCoord(), FEA.getTabmat(), FEA.getTabgeo(), FEA.getIntGauss(), 0
            ... )
        """
        return self.model.element.getStifLinearMat(
            self.model, inci, coord, tabmat, tabgeo, intgauss, element_number
        )

    def getElemMassConsistentMat(
        self, inci: npt.NDArray[np.float64], coord: npt.NDArray[np.float64], tabmat: list, tabgeo: list, intgauss: int, element_number: int
    ) -> npt.NDArray[np.float64]:
        """Computes the element mass matrix using a consistent formulation.

        Args:
            inci: Nodal incidence matrix.
            coord: Global coordinate coordinates.
            tabmat: Material properties list.
            tabgeo: Geometry properties profile.
            intgauss: Gauss numerical integration order.
            element_number: Target element list index.

        Returns:
            A consistent local element mass matrix of shape (Dofs, Dofs).

        Example:
            >>> m_local = FEA.getElemMassConsistentMat(
            ...     FEA.getInci(), FEA.getCoord(), FEA.getTabmat(), FEA.getTabgeo(), FEA.getIntGauss(), 0
            ... )
        """
        return self.model.element.getMassConsistentMat(
            self.model, inci, coord, tabmat, tabgeo, intgauss, element_number
        )
    

    # GET SOLVER
    def getGlobalMatrix(self, Model, inci:npt.NDArray[np.float64] = None, coord:npt.NDArray[np.float64] = None, tabmat:list = None, tabgeo:list = None, intgauss:int = None, SYMM:bool=None, MP:bool=None) -> npt.NDArray[np.float64]:
        """Invokes the solver assembler to construct the global unconstrained system matrices.

        Args:
            Model: Active SetModel structure containing element/constitutive classes.
            inci: Optional connectivity matrix. Defaults to None.
            coord: Optional spatial coordinate array. Defaults to None.
            tabmat: Optional material lookup profile. Defaults to None.
            tabgeo: Optional geometry lookup profile. Defaults to None.
            intgauss: Optional Gaussian integration order. Defaults to None.
            SYMM: Enable storage/matrix symmetries mapping. Defaults to None.
            MP: Multiprocessing support option. Defaults to None.

        Returns:
            The raw assembled structural matrix system (unconstrained global equations).

        Example:
            >>> K_global = FEA.getGlobalMatrix(FEA.model, SYMM=False, MP=False)
        """
        return self.solver.getMatrixAssembler(Model, inci = inci, coord = coord, tabmat = tabmat, tabgeo = tabgeo, intgauss = intgauss, SYMM=SYMM, MP=MP)

    def getConstrains(self, constrains:list) -> npt.NDArray[np.float64]:
        """Maps boundary condition parameters to explicit indices classifications.

        Groups DOFs into free degrees of freedom, fixed nodes, and predefined Dirichlet
        condition values.

        Args:
            constrains: list of boundary conditions dict values.

        Returns:
            A tuple of three arrays containing index representations of:
                - freedof: Indices of unrestricted DOFs.
                - fixedof: Indices of constrained DOFs.
                - constdof: Prescribed displacement values array.

        Example:
            >>> free, fixed, vals = FEA.getConstrains(FEA.physic.constrains)
        """
        nodetot = len(self.model.coord)
        return self.solver.getConstrains(
            constrains, nodetot, self.model.modelinfo["nodedof"]
        )

    def getDirichletNH(self, constrains:list) -> npt.NDArray[np.float64]:
        """Builds Dirichlet non-homogeneous boundary value vectors.

        Args:
            constrains: Physical boundary constraints setup configuration list.

        Returns:
            An array vector defining prescribed non-zero values corresponding to restricted DOFs.

        Example:
            >>> U_dirichlet = FEA.getDirichletNH(FEA.physic.constrains)
        """
        nodetot = len(self.model.coord)
        return self.solver.getDirichletNH(
            constrains, nodetot, self.model.modelinfo["nodedof"]
        )

    def getLoadArray(self, loadaply:list) -> npt.NDArray[np.float64]:
        """Assembles local element/nodal actions into the global algebraic force vector.

        Args:
            loadaply: Array of nodal loads defined in the physics manager.

        Returns:
            An assembled global array containing load force definitions.

        Example:
            >>> F_load = FEA.getLoadArray(FEA.physic.forces)
        """
        nodetot = len(self.model.coord)
        return self.solver.getLoadAssembler(
            loadaply, nodetot, self.model.modelinfo["nodedof"]
        )

    # GET PHYSIC
    def getPhysic(self) -> object:
        """Retrieves the physical boundary conditions manager object.

        Returns:
            The SetPhysics instance managing loading definitions and constraints.

        Example:
            >>> physics_obj = FEA.getPhysic()
        """
        return self.physic

    # [bug]
    # def getForceList(self) -> list:
    #     """get forces list

    #     Returns:
    #         _description_
    #     """
    #     return self.physic.getForceList(self.modelinfo["domain"])

    # def getBoundCondList(self) -> list:
    #     """get boundary conditions list

    #     Returns:
    #         _description_
    #     """
    #     return self.physic.getBoundCondList(self.modelinfo["domain"])

    def getLoadApply(self) -> npt.NDArray[np.float64]:
        """Compiles active load definitions into mathematical vector profiles.

        Returns:
            An array detailing node indices, degrees of freedom indices, and scaling loads.

        Example:
            >>> active_loads = FEA.getLoadApply()
        """
        return self.physic.getLoadApply(self.physic.physicdata)

    def getBCApply(self) -> npt.NDArray[np.float64]:
        """Compiles boundary constraint rules into mathematical matrix profiles.

        Returns:
            An array of constraint properties specifying target indices, degree of freedom,
            and enforcement scales.

        Example:
            >>> active_constraints = FEA.getBCApply()
        """
        return self.physic.getBoundCondApply(self.physic.physicdata)
    
    def getCouplingInterface(self) -> list:
        """Retrieves interaction forces corresponding to physical domain couplings.

        Returns:
            A list detailing force components on multiphysics interface coupling nodes.

        Example:
            >>> coupling_forces = FEA.getCouplingInterface()
        """
        return self.physic.getLoadCoup(self.physic.physicdata)

    def getUpdateMatrix(self, matrix, addval) -> npt.NDArray[np.float64]:
        """Applies algebraic modifications to the global governing system matrix.

        Applies corrections for concentrated parameters (such as lumped point masses, elastic foundations,
        or penalty terms).

        Args:
            matrix: Assembled unconstrained system matrix array.
            addval: Additional modifications matrix array or force properties array.

        Returns:
            The corrected global system matrix.

        Example:
            >>> K_updated = FEA.getUpdateMatrix(K_global, active_loads)
        """
        return self.physic.getUpdateMatrix(matrix, addval)
    
    def getElementFromNodesList(self, nodelist) -> list:
        """Identifies elements associated with a specific list of nodal identifiers.

        Args:
            nodelist: Nodal tags list.

        Returns:
            A list containing element indices connected to those nodes.

        Example:
            >>> elements = FEA.getElementFromNodesList([5, 8])
        """
        return self.physic.getElementList(self.model.inci, nodelist)

    # OUTHERS
    def getNodesFromRegions(self, set: int, type: str) -> list:
        """Retrieves node tags belonging to a GMSH physical group region.

        Args:
            set: Region target integer ID index.
            type: Physical group geometric entity level type ('point', 'line', 'plane').

        Returns:
            A list containing nodes belonging to the requested domain region group.

        Example:
            >>> boundary_nodes = FEA.getNodesFromRegions(set=2, type="line")
        """
        if type == "point":
            domain_nodelist = newAnalysis.getRegions(self)[0][1][set - 1][1]
        elif type == "line":
            domain_nodelist = newAnalysis.getRegions(self)[1][1][set - 1][1]
        elif type == "plane":
            domain_nodelist = newAnalysis.getRegions(self)[2][1][set - 1][1]
        else:
            domain_nodelist = []
        return domain_nodelist
    
    def getAbout(self):
        """Prints program credits, author profile, and package version on terminal."""
        get_about()

    # settings FEA ANALYSIS <privates>
    def __setMesh(modeldata):
        set_mesh = dict()
        set_mesh = modeldata["MESH"]
        set_mesh["SHAPE"] = modeldata["ELEMENT"]["SHAPE"]
        set_mesh["user_path"] = modeldata["MESH"]["user_path"]
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
