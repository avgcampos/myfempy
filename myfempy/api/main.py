from __future__ import annotations

import logging
import sys
import os
import sysconfig
from time import time
from datetime import datetime

import numpy as np
import numpy.typing as npt

from myfempy.core.utilities import setSteps, gauss_points
# from myfempy.core.solver import getSolver
from myfempy.io.controllers import (setElement, setGeometry, setMaterial,
                                    setMesh, setShape, setDomain, setCoupling,
                                    setPoints2NumericalIntegration)

from myfempy.plots.prevplot import preview_plot
from myfempy.api.model import SetModel
from myfempy.api.physics import SetPhysics
from myfempy.api.results import setPostProcess
from myfempy.utils.utils import (clear_console, get_logo, get_version,
                                 loading_bar_v1, newDir, print_console)

__docformat__ = "google"
with open(os.getcwd()+'/myfempy/utils/about.txt', 'r', encoding='utf-8') as file:
    __doc__ = file.read()

class newAnalysis:
    """
    Setup the New Analysis to FEA simulation
    """
    def __init__(self, FEASolver: object, path: str = None) -> None:
        """Initializes a new Finite Element Analysis (FEA) project environment[cite: 1].

        Sets up the solution directories, default logging configuration, and binds the
        numerical solver module to the simulation instance[cite: 1].

        Args:
            FEASolver (object): Class or module responsible for solving the physical state equations 
                (e.g., StaticLinear, SteadyStateLinear).
            path (str, optional): Target directory path for exporting simulation logs and output files. 
                If None is passed, defaults to creating an "out" folder.

        Returns:
            None

        Example:
            >>> from myfempy.api.api import newAnalysis
            >>> from myfempy.core.solver import SteadyStateLinear
            >>> FEA = newAnalysis(FEASolver=SteadyStateLinear, path="simulation_results")
        """
        get_logo()
        try:
            now = datetime.now()
            # Format: Day/Month/Year Hour:Minute:Second
            self.timenow = now.strftime("%d/%m/%Y %H:%M:%S")
            self.solver = FEASolver
            text_init = "TRY SET NEW ANALYSIS AND SOLVER -- SUCCESS"
        except:
            text_init = "TRY SET NEW ANALYSIS AND SOLVER -- FAULT"
        try:
            self.path = newDir(path)
        except:
            print(">>> User save folder not found, creating 'out' folder")
            self.path = newDir("out")
            
        logging.basicConfig(
            filename=str(self.path) + "/api-log.log",
            encoding="utf-8",
            level=logging.DEBUG,
            filemode="w",
        )
        logging.info(text_init)

    def Model(self, modeldata: dict) -> None:
        """Sets up the mesh, element configuration, materials, and domain geometry[cite: 1].

        Parses the unified configuration dictionary to build the mathematical model representation,
        computes element volumes, and compiles internal metadata arrays (such as coordinates
        and connectivity)[cite: 1].

        Args:
            modeldata (dict): A structured configuration dictionary containing modeling setup:
                - 'MESH': Coordinates and element incidence.
                - 'ELEMENT': Finite element formulation, shape, and integration order.
                - 'MATERIAL': Constitutive model and properties.
                - 'GEOMETRY': Thickness, cross-section areas, or generic geometric descriptors.

        Returns:
            None

        Example:
            >>> model_config = {
            ...     'MESH': {'TYPE': 'manual', 'COORD': [...], 'INCI': [...]},
            ...     'ELEMENT': {'TYPE': 'structplane', 'SHAPE': 'quad4', 'INTGAUSS': 2},
            ...     'MATERIAL': {'MAT': 'planestress', 'TYPE': 'isotropic', 'PROPMAT': [...]},
            ...     'GEOMETRY': {'GEO': 'thickness', 'PROPGEO': [...]}
            ... }
            >>> FEA.Model(model_config)
        """
        print_console("mesh")
        try:
            modeldata["MESH"]["user_path"] = self.path
            Mesh = newAnalysis.__setMesh(modeldata)
            Element = newAnalysis.__setElement(modeldata)
            Shape = newAnalysis.__setShape(modeldata)
            Material = newAnalysis.__setMaterial(modeldata)
            Geometry = newAnalysis.__setGeometry(modeldata)
            GaussPoints = newAnalysis.__setIntGauss(modeldata)
            logging.info("TRY SET MODEL COMPONENTS -- SUCCESS")
        except KeyError as e:
            logging.error(f"Missing required key in modeldata configuration: {e}")
            raise
        except Exception as e:
            logging.error(f"TRY SET MODEL COMPONENTS -- FAULT: {e}")
            raise

        try:
            self.model = SetModel(Mesh, Element, Shape, Material, Geometry)
            self.model.modeldata = modeldata
            logging.info("TRY SET FEMODEL -- SUCCESS")
        except Exception as e:
            logging.error(f"TRY SET FEMODEL -- FAULT: {e}")
            raise
            
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
        except Exception:
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

        self.model.elemvol = newAnalysis.__setMeshElemVol(self)

    def Physic(self, physicdata: dict) -> None:
        """Configures loads, coupling fields, and boundary constraints[cite: 1].

        Initializes force vectors, multiphysics couplings, and kinematic boundary
        constraints (such as essential Dirichlet conditions) on specified coordinates or nodal lists[cite: 1].

        Args:
            physicdata (dict): Configuration dictionary specifying boundary conditions and loads:
                - 'PHYSIC': Holds physical 'DOMAIN' type, 'LOAD' array, and 'BOUNDCOND' array.
                - 'COUPLING' (optional): Holds interface multi-field configuration.

        Returns:
            None

        Example:
            >>> physics_config = {
            ...     'PHYSIC': {
            ...         'DOMAIN': 'structural',
            ...         'LOAD': [...],
            ...         'BOUNDCOND': [...]
            ...     }
            ... }
            >>> FEA.Physic(physics_config)
        """
        print_console("phy")
        try:
            Loads, BoundCond = newAnalysis.__setDomain(physicdata)
            self.physic = SetPhysics(self.model, Loads, BoundCond)
            self.physic.physicdata = physicdata
            logging.info("TRY SET PHYSICS -- SUCCESS")
        except Exception as e:
            self.physic = []
            logging.error(f"TRY SET PHYSICS -- FAULT: {e}")
            raise

        try:
            self.physic.forces = newAnalysis.getLoadApply(self)
            logging.info("TRY SET PHYSICS.FORCES -- SUCCESS")
        except Exception as e:
            self.physic.forces = []
            logging.warning(f"TRY SET PHYSICS.FORCES -- FAULT: {e}")

        if "COUPLING" in physicdata:
            coupling_load_zero = self.physic.forces
            LoadCoup, BoundCond = newAnalysis.__setCoupling(physicdata)
            self.physic = SetPhysics(self.model, LoadCoup, BoundCond)
            self.physic.physicdata = physicdata
            self.physic.forces = np.append(
                coupling_load_zero, newAnalysis.getCouplingInterface(self), axis=0
            )
            
        try:
            constrains = newAnalysis.getBCApply(self)
            if constrains.size > 0 and any(constrains[:, 1] == 11):
                self.physic.csleft = constrains[np.where(constrains[:, 1] == 11)[0], 0]
                self.physic.csright = constrains[np.where(constrains[:, 1] == 12)[0], 0]
            self.physic.constrains = constrains
            logging.info("TRY SET PHYSICS.CONSTRAINS -- SUCCESS")
        except Exception as e:
            self.physic.constrains = []
            logging.warning(f"TRY SET PHYSICS.CONSTRAINS -- FAULT: {e}")

    def Assembly(self, Model: object) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
        """Assembles the element equations into the global algebraic system of equations[cite: 1].

        Computes the global linearized matrix structure and updates it incorporating force-term 
        arrays under the current simulation model[cite: 1].

        Args:
            Model (object): Active Model instance containing geometry, material, and integration data[cite: 1].

        Returns:
            tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]: A tuple containing:
                - matrix: The assembled global matrix array.
                - forcelist: The fully-assembled global force load vector list.

        Example:
            >>> K_global, F_global = FEA.Assembly(Model=FEA.model)
        """
        is_free_threaded = sysconfig.get_config_var('Py_GIL_DISABLED') == 1
        self.MP = False
        if sys.version_info >= (3, 14) and is_free_threaded:
            self.MP = True
        try:
            matrix = newAnalysis.getGlobalMatrix(self, Model, self.model.inci, self.model.coord, self.model.tabmat, self.model.tabgeo, self.model.intgauss, self.MP)
            loadaply = self.physic.forces
            matrix = newAnalysis.getUpdateMatrix(self, matrix, loadaply)
            forcelist = newAnalysis.getLoadArray(self, loadaply)
            logging.info("TRY RUN GLOBAL ASSEMBLY -- SUCCESS")
            return matrix, forcelist
        except Exception as e:
            logging.error(f"TRY RUN GLOBAL ASSEMBLY -- FAULT: {e}")
            raise

    def Solve(self, solverset: dict = None) -> dict:
        """Executes the finite element equations solver over designated steps[cite: 1].

        Compiles matrix system transformations, processes Dirichlet non-homogeneous values,
        assembles loads across time/loading steps, and calls the bound numerical analysis core[cite: 1].

        Args:
            solverset (dict, optional): Configuration parameters dictionary for the numeric solver 
                including 'STEPSET'.

        Returns:
            dict: A modified 'solverset' dictionary containing simulation status, calculation logs,
            times elapsed, and simulation results under 'solverset["solution"]'.

        Example:
            >>> run_config = {'STEPSET': {'type': 'table', 'start': 0.0, 'end': 1.0, 'step': 1.0}}
            >>> results = FEA.Solve(run_config)
        """
        print_console("solver")
        print(">>> RUNNING SOLVER:")
        print(self.solver.__doc__)
        
        solverset["solverstatus"] = {
            "solvercore": self.solver.__doc__,
            "myfempyversion": get_version()
        }
        
        starttime = time()
        assembly, forcelist = self.Assembly(Model=self.model)
        solverset["solverstatus"]["timeasb"] = abs(time() - starttime)
        solverset["solverstatus"]["memorysize"] = (assembly["stiffness"].todense().nbytes) / 1e6
        if self.MP:
            solverset["solverstatus"]["typeasmb"] = (
                "PARALLEL_" + str(os.cpu_count()) + "_CORES"
            )
        else:
            solverset["solverstatus"]["typeasmb"] = "SERIAL_1_CORE"
        

        try:
            constrains = self.physic.constrains
            freedof, fixedof, constdof = newAnalysis.getConstrains(self, constrains)
            logging.info("TRY RUN CONSTRAINS -- SUCCESS")
        except Exception as e:
            logging.error(f"TRY RUN CONSTRAINS -- FAULT: {e}")
            raise

        constrainsdof = {
            "freedof": freedof,
            "fixedof": fixedof,
            "constdof": constdof
        }
        
        try:
            Uc = newAnalysis.getDirichletNH(self, constrains)
            logging.info("TRY SET DNH CONSTRAINS -- SUCCESS")
        except Exception as e:
            logging.error(f"TRY SET DNH CONSTRAINS -- FAULT: {e}")
            raise

        nsteps = setSteps(solverset["STEPSET"])
        if forcelist.shape[1] != nsteps:
            forcelist = np.repeat(forcelist, nsteps, axis=1)
        assembly["loads"] = forcelist

        if Uc.shape[1] != nsteps:
            Uc = np.repeat(Uc, nsteps, axis=1)
        assembly["bcdirnh"] = Uc

        try:
            starttime = time()
            solverset["solution"] = self.solver.runSolve(self.model, self.physic, assembly, constrainsdof, solverset)
            solverset["solverstatus"]["timesim"] = abs(time() - starttime)
            logging.info("TRY RUN SOLVER -- SUCCESS")
        except Exception as e:
            logging.error(f"TRY RUN SOLVER -- FAULT: {e}")
            raise

        print_console("succ")
        return solverset

    def PreviewAnalysis(self, previewdata: dict) -> None:
        """Renders pre-simulation plots for physical inspection of modeling items.

        Draws geometry shapes, elements, nodes, and applied load vectors before 
        running the solver.

        Args:
            previewdata (dict): Graphical visualization configuration containing rendering attributes.

        Returns:
            None

        Example:
            >>> preview_config = {'RENDER': {'show': True, 'scale': 2.5}}
            >>> FEA.PreviewAnalysis(preview_config)
        """        
        is_free_threaded = sysconfig.get_config_var('Py_GIL_DISABLED') == 1
        if sys.version_info >= (3, 14) and is_free_threaded:
            print("AVISO: O recurso de plotagem (preview_plot) é incompatível com a versão do Python instalada")
            logging.warning("PREVIEW PLOT SKIPPED -- Incompatible with Python Version")
            return
        try:
            preview_plot(self.model, previewdata, str(self.path), self.physic)
            logging.info("TRY RUN PREVIEW PLOT -- SUCCESS")
        except Exception:
            preview_plot(self.model, previewdata, str(self.path))
            logging.warning("TRY RUN PREVIEW PLOT (WITHOUT PHYSIC) -- FALLBACK")

    def PostProcess(self, postprocset: dict) -> list | dict:
        """Processes solutions, computes auxiliary fields, and builds output reports[cite: 1].

        Translates primary state parameters into derivative properties and creates 
        visualization plots and text records[cite: 1].

        Args:
            postprocset (dict): Directives dictionary defining post-processing metrics.

        Returns:
            list | dict: A post-processed analysis dataset mapping physical variables, plots, 
            and log file paths.

        Example:
            >>> config = {"COMPUTER": {'structural': {'displ': True}}}
            >>> data = FEA.PostProcess(config)
        """
        print_console("post")
        postprocdata = []
        try:
            if "COMPUTER" in postprocset:
                postprocdata = setPostProcess.getCompute(self, postprocset)
            if "PLOT" in postprocset:
                setPostProcess.getPlotCSV(self, postprocset, postprocdata)
            if "REPORT" in postprocset:
                postprocdata["log"] = []
                log_file = setPostProcess.getLog(self, postprocset, postprocdata)
                postprocdata["log"].append(log_file)
            logging.info("TRY GET POST PROCESS -- SUCCESS")
        except Exception as e:
            logging.error(f"TRY GET POST PROCESS -- FAULT: {e}")
        
        print_console("thank")
        return postprocdata

    def getModel(self) -> object:
        """Retrieves the active finite element model container object[cite: 1].

        Returns:
            object: The SetModel instance describing mesh, connectivity, and formulations.

        Example:
            >>> model_obj = FEA.getModel()
        """
        return self.model

    def getModelInfo(self) -> dict:
        """Retrieves summary attributes and counts of degrees of freedom from the model[cite: 1].

        Returns:
            dict: A dictionary containing structural details like elements count and DOFs.

        Example:
            >>> info = FEA.getModelInfo()
        """
        return self.model.modelinfo

    def getInci(self) -> npt.NDArray[np.float64]:
        """Retrieves the element incidence and connectivity matrix[cite: 1].

        Returns:
            npt.NDArray[np.float64]: An array listing element indices and node references.

        Example:
            >>> connectivity = FEA.getInci()
        """
        return self.model.inci

    def getCoord(self) -> npt.NDArray[np.float64]:
        """Retrieves the global spatial coordinates table of all mesh nodes[cite: 1].

        Returns:
            npt.NDArray[np.float64]: An array mapping node tags to spatial coordinates.

        Example:
            >>> nodes_coords = FEA.getCoord()
        """
        return self.model.coord

    def getTabmat(self) -> list:
        """Retrieves the material properties configuration table[cite: 1].

        Returns:
            list: A list of dictionary objects matching active materials.

        Example:
            >>> materials = FEA.getTabmat()
        """
        return self.model.tabmat

    def getTabgeo(self) -> list:
        """Retrieves the geometric cross-section/thickness attributes table[cite: 1].

        Returns:
            list: A list containing properties for geometric profiles.

        Example:
            >>> thickness_props = FEA.getTabgeo()
        """
        return self.model.tabgeo

    def getIntGauss(self) -> int:
        """Retrieves the number of Gauss integration points for numerical integration[cite: 1].

        Returns:
            int: The integration order/points parameter.

        Example:
            >>> order = FEA.getIntGauss()
        """
        return self.model.intgauss
    
    def getRegions(self) -> list:
        """Retrieves mesh entity group regions imported from files[cite: 1].

        Returns:
            list: A nested list of elements grouped by geometric entity tag definitions.

        Example:
            >>> physical_groups = FEA.getRegions()
        """
        return self.model.regions

    def getElementVolume(self) -> npt.NDArray[np.float64]:
        """Calculates structural volumes (or areas/lengths) for all mesh elements[cite: 1].

        Args:
            inci (npt.NDArray[np.float64]): Nodal incidence array of the elements[cite: 1].
            coord (npt.NDArray[np.float64]): Nodal spatial coordinate array[cite: 1].
            tabgeo (list): Geometry properties collection[cite: 1].

        Returns:
            npt.NDArray[np.float64]: A 1D numpy array containing the computed volume/area value for each element.

        Example:
            >>> volumes = FEA.getElementVolume(FEA.model.inci, FEA.model.coord, FEA.model.tabgeo)
        """
        return self.model.elemvol


    def getElemStifLinearMat(
        self, inci: npt.NDArray[np.float64], coord: npt.NDArray[np.float64], tabmat: list, tabgeo: list, intgauss: int, element_number: int
    ) -> npt.NDArray[np.float64]:
        """Computes the element linear stiffness matrix[cite: 1].

        Args:
            inci (npt.NDArray[np.float64]): Nodal incidence matrix[cite: 1].
            coord (npt.NDArray[np.float64]): Global coordinate coordinates[cite: 1].
            tabmat (list): Material configuration properties[cite: 1].
            tabgeo (list): Geometry properties profile[cite: 1].
            intgauss (int): Order count of integration[cite: 1].
            element_number (int): Index identifier of the element[cite: 1].

        Returns:
            npt.NDArray[np.float64]: A 2D array representing the local stiffness matrix.

        Example:
            >>> k_local = FEA.getElemStifLinearMat(FEA.getInci(), FEA.getCoord(), FEA.getTabmat(), FEA.getTabgeo(), FEA.getIntGauss(), 0)
        """
        C = self.model.material.getElasticTensor(tabmat, inci, element_number)
        shape_set = self.model.shape.getShapeSet()
        nodecon = len(shape_set["nodes"])
        elem_set = self.model.element.getElementSet()
        nodedof = len(elem_set["dofs"]["d"])
        elemdof = nodecon * nodedof
        nodelist = self.model.shape.getNodeList(self.model.inci, element_number)
        elementcoord = self.model.shape.getNodeCoord(self.model.coord, nodelist)
        getIntNum = self.model.shape.getIntNumK
        type_shape = shape_set["key"]
        point_gauss, weight_gauss = gauss_points(type_shape, intgauss)
        return self.model.element.getStifLinearMat(
            inci, coord, tabmat, tabgeo, elementcoord, C, elemdof, getIntNum, intgauss, point_gauss, weight_gauss, element_number
        )

    def getElemMassConsistentMat(
        self, inci: npt.NDArray[np.float64], coord: npt.NDArray[np.float64], tabmat: list, tabgeo: list, intgauss: int, element_number: int
    ) -> npt.NDArray[np.float64]:
        """Computes the element mass matrix using a consistent formulation[cite: 1].

        Args:
            inci (npt.NDArray[np.float64]): Nodal incidence matrix[cite: 1].
            coord (npt.NDArray[np.float64]): Global coordinate coordinates[cite: 1].
            tabmat (list): Material properties list[cite: 1].
            tabgeo (list): Geometry properties profile[cite: 1].
            intgauss (int): Gauss numerical integration order[cite: 1].
            element_number (int): Target element list index[cite: 1].

        Returns:
            npt.NDArray[np.float64]: A consistent local element mass matrix of shape (Dofs, Dofs).

        Example:
            >>> m_local = FEA.getElemMassConsistentMat(FEA.getInci(), FEA.getCoord(), FEA.getTabmat(), FEA.getTabgeo(), FEA.getIntGauss(), 0)
        """
        C = self.model.material.getElasticTensor(tabmat, inci, element_number)
        shape_set = self.model.shape.getShapeSet()
        nodecon = len(shape_set["nodes"])
        elem_set = self.model.element.getElementSet()
        nodedof = len(elem_set["dofs"]["d"])
        elemdof = nodecon * nodedof
        nodelist = self.model.shape.getNodeList(self.model.inci, element_number)
        elementcoord = self.model.shape.getNodeCoord(self.model.coord, nodelist)
        getIntNum = self.model.shape.getIntNumM
        type_shape = shape_set["key"]
        point_gauss, weight_gauss = gauss_points(type_shape, intgauss)
        return self.model.element.getMassConsistentMat(
            inci, coord, tabmat, tabgeo, elementcoord, C, elemdof, getIntNum, intgauss, point_gauss, weight_gauss, element_number
        )

    def getGlobalMatrix(self, Model: object, inci: npt.NDArray[np.float64] = None, coord: npt.NDArray[np.float64] = None, tabmat: list = None, tabgeo: list = None, intgauss: int = None, MP: bool = None, max_workers: int = None) -> npt.NDArray[np.float64]:
        """Invokes the solver assembler to construct the global unconstrained system matrices[cite: 1].

        Args:
            Model (object): Active SetModel structure containing element/constitutive classes[cite: 1].
            inci (npt.NDArray[np.float64], optional): Optional connectivity matrix. Defaults to None.
            coord (npt.NDArray[np.float64], optional): Optional spatial coordinate array. Defaults to None.
            tabmat (list, optional): Optional material lookup profile. Defaults to None.
            tabgeo (list, optional): Optional geometry lookup profile. Defaults to None.
            intgauss (int, optional): Optional Gaussian integration order. Defaults to None.
            MP (int, optional): Optional MultiProcessing Assembler. Defaults to None.
        Returns:
            npt.NDArray[np.float64]: The raw assembled structural matrix system.

        Example:
            >>> K_global = FEA.getGlobalMatrix(FEA.model)
        """
        return self.solver.getMatrixAssembler(Model, inci=inci, coord=coord, tabmat=tabmat, tabgeo=tabgeo, intgauss=intgauss, MP=MP, max_workers = max_workers)

    def getConstrains(self, constrains: list) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64], npt.NDArray[np.float64]]:
        """Maps boundary condition parameters to explicit indices classifications[cite: 1].

        Args:
            constrains (list): List of boundary conditions dict values[cite: 1].

        Returns:
            tuple[npt.NDArray[np.float64], npt.NDArray[np.float64], npt.NDArray[np.float64]]: Arrays for freedof, fixedof, and constdof.

        Example:
            >>> free, fixed, vals = FEA.getConstrains(FEA.physic.constrains)
        """
        nodetot = len(self.model.coord)
        return self.solver.getConstrains(constrains, nodetot, self.model.modelinfo["nodedof"])

    def getDirichletNH(self, constrains: list) -> npt.NDArray[np.float64]:
        """Builds Dirichlet non-homogeneous boundary value vectors[cite: 1].

        Args:
            constrains (list): Physical boundary constraints setup configuration list[cite: 1].

        Returns:
            npt.NDArray[np.float64]: An array vector defining prescribed non-zero values corresponding to restricted DOFs.

        Example:
            >>> U_dirichlet = FEA.getDirichletNH(FEA.physic.constrains)
        """
        nodetot = len(self.model.coord)
        return self.solver.getDirichletNH(constrains, nodetot, self.model.modelinfo["nodedof"])

    def getLoadArray(self, loadaply: list) -> npt.NDArray[np.float64]:
        """Assembles local element/nodal actions into the global algebraic force vector[cite: 1].

        Args:
            loadaply (list): Array of nodal loads defined in the physics manager[cite: 1].

        Returns:
            npt.NDArray[np.float64]: An assembled global array containing load force definitions.

        Example:
            >>> F_load = FEA.getLoadArray(FEA.physic.forces)
        """
        nodetot = len(self.model.coord)
        return self.solver.getLoadAssembler(loadaply, nodetot, self.model.modelinfo["nodedof"])

    def getPhysic(self) -> object:
        """Retrieves the physical boundary conditions manager object[cite: 1].

        Returns:
            object: The SetPhysics instance managing loading definitions and constraints.

        Example:
            >>> physics_obj = FEA.getPhysic()
        """
        return self.physic

    def getLoadApply(self) -> npt.NDArray[np.float64]:
        """Compiles active load definitions into mathematical vector profiles[cite: 1].

        Returns:
            npt.NDArray[np.float64]: An array detailing node indices, degrees of freedom indices, and scaling loads.

        Example:
            >>> active_loads = FEA.getLoadApply()
        """
        return self.physic.getLoadApply(self.physic.physicdata)

    def getBCApply(self) -> npt.NDArray[np.float64]:
        """Compiles boundary constraint rules into mathematical matrix profiles[cite: 1].

        Returns:
            npt.NDArray[np.float64]: An array of constraint properties.

        Example:
            >>> active_constraints = FEA.getBCApply()
        """
        return self.physic.getBoundCondApply(self.physic.physicdata)
    
    def getCouplingInterface(self) -> list:
        """Retrieves interaction forces corresponding to physical domain couplings[cite: 1].

        Returns:
            list: A list detailing force components on multiphysics interface coupling nodes.

        Example:
            >>> coupling_forces = FEA.getCouplingInterface()
        """
        return self.physic.getLoadCoup(self.physic.physicdata)

    def getUpdateMatrix(self, matrix: npt.NDArray[np.float64], addval: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        """Applies algebraic modifications to the global governing system matrix[cite: 1].

        Args:
            matrix (npt.NDArray[np.float64]): Assembled unconstrained system matrix array[cite: 1].
            addval (npt.NDArray[np.float64]): Additional modifications matrix array or force properties array[cite: 1].

        Returns:
            npt.NDArray[np.float64]: The corrected global system matrix.

        Example:
            >>> K_updated = FEA.getUpdateMatrix(K_global, active_loads)
        """
        return self.physic.getUpdateMatrix(matrix, addval)
    
    def getElementFromNodesList(self, nodelist: list) -> list:
        """Identifies elements associated with a specific list of nodal identifiers[cite: 1].

        Args:
            nodelist (list): Nodal tags list[cite: 1].

        Returns:
            list: A list containing element indices connected to those nodes.

        Example:
            >>> elements = FEA.getElementFromNodesList([5, 8])
        """
        return self.physic.getElementList(self.model.inci, nodelist)

    def getNodesFromRegions(self, set: int, type: str) -> list:
        """Retrieves node tags belonging to a physical group region[cite: 1].

        Args:
            set (int): Region target integer ID index[cite: 1].
            type (str): Physical group geometric entity level type ('point', 'line', 'plane')[cite: 1].

        Returns:
            list: A list containing nodes belonging to the requested domain region group.

        Example:
            >>> boundary_nodes = FEA.getNodesFromRegions(set=2, type="line")
        """
        regions = newAnalysis.getRegions(self)
        if type == "point":
            domain_nodelist = regions[0][1][set - 1][1]
        elif type == "line":
            domain_nodelist = regions[1][1][set - 1][1]
        elif type == "plane":
            domain_nodelist = regions[2][1][set - 1][1]
        else:
            domain_nodelist = []
        return domain_nodelist

    @staticmethod
    def __setMesh(modeldata: dict) -> object:
        set_mesh = dict(modeldata["MESH"])
        set_mesh["SHAPE"] = modeldata["ELEMENT"]["SHAPE"]
        set_mesh["user_path"] = modeldata["MESH"]["user_path"]
        return setMesh(set_mesh)

    @staticmethod
    def __setShape(modeldata: dict) -> object:
        return setShape(modeldata["ELEMENT"])

    @staticmethod
    def __setElement(modeldata: dict) -> object:
        return setElement(modeldata["ELEMENT"])

    @staticmethod
    def __setMaterial(modeldata: dict) -> object:
        return setMaterial(modeldata["MATERIAL"])

    @staticmethod
    def __setGeometry(modeldata: dict) -> object:
        return setGeometry(modeldata["GEOMETRY"])

    @staticmethod
    def __setDomain(modeldata: dict) -> tuple:
        return setDomain(modeldata["PHYSIC"])

    @staticmethod
    def __setCoupling(modeldata: dict) -> tuple:
        return setCoupling(modeldata["COUPLING"])
    
    @staticmethod
    def __setIntGauss(modeldata: dict) -> int:
        if "INTGAUSS" in modeldata["ELEMENT"]:
            intgauss = modeldata["ELEMENT"]["INTGAUSS"]
        else:
            intgauss = setPoints2NumericalIntegration(modeldata["ELEMENT"]["SHAPE"])
        return intgauss

    @staticmethod
    def __setMeshElemVol(self):
        vol = np.zeros((self.model.inci.shape[0]))
        for ee in range(self.model.inci.shape[0]):
            nodelist = self.model.shape.getNodeList(self.model.inci, ee)
            elementcoord = self.model.shape.getNodeCoord(self.model.coord, nodelist)
            vol[ee] = self.model.element.getElementVolume(self.model.inci, self.model.tabgeo,
                                                        self.model.shape.getVOL,
                                                        self.model.modelinfo["type_shape"],
                                                        elementcoord,
                                                        ee)
        return vol
        
