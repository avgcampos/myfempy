from __future__ import annotations

import numpy as np

from myfempy.core.mesh import getMesh
# from myfempy.core.solver import getSolver
from myfempy.felib.elements.element import getElement
from myfempy.felib.geometry.geometry import getGeometry
from myfempy.felib.material.material import getMaterial
from myfempy.felib.shapes.shape import getShape
from myfempy.plots.prevplot import preview_plot
from myfempy.setup.model import SetModel
from myfempy.setup.physics import SetPhysics
from myfempy.setup.results import setPostProcess


class newAnalysis():
    '''Finite Elements Analysis Class <ClassOrder>'''
    
    def __init__(self, FEASolver):
        self.solver = FEASolver
        
                       
    def Model(self, modeldata):
                
        Mesh = newAnalysis.setMesh(modeldata)
        Shape = newAnalysis.setShape(modeldata)
        Element = newAnalysis.setElement(modeldata)
        Material = newAnalysis.setMaterial(modeldata)
        Geometry = newAnalysis.setGeometry(modeldata)
               
        self.model = SetModel(Mesh, Element, Shape, Material, Geometry)
        self.model.modeldata = modeldata
        self.model.element = Element
        self.model.shape = Shape
        self.model.material = Material
        self.model.geoemtry = Geometry
        
        # self.inci = FEANewAnalysis.getInci(self)
        # self.coord = FEANewAnalysis.getCoord(self)
        # self.tabmat = FEANewAnalysis.getTabmat(self)
        # self.tabgeo = FEANewAnalysis.getTabgeo(self)
        # self.intgauss = FEANewAnalysis.getIntGauss(self)

        self.modelinfo = dict()
        self.modelinfo['inci'] = newAnalysis.getInci(self)
        self.modelinfo['coord'] = newAnalysis.getCoord(self)
        self.modelinfo['tabmat'] = newAnalysis.getTabmat(self)
        self.modelinfo['tabgeo'] = newAnalysis.getTabgeo(self)
        self.modelinfo['intgauss'] = newAnalysis.getIntGauss(self)
        try:
            self.modelinfo["regions"] = self.model.mesh.getRegionsList(self.model.mesh.getElementConection(self.model.modeldata['MESH']))
        except:
            self.modelinfo["regions"] = []
                    

        elem_set = self.model.element.getElementSet()
        
        self.modelinfo['tensor'] = len(elem_set['tensor'])
        self.modelinfo['dofs'] = elem_set["dofs"]
        self.modelinfo['nodedof'] = len(elem_set["dofs"]['d'])
        self.modelinfo['type_element'] = elem_set["key"]

        shape_set = self.model.shape.getShapeSet()
        self.modelinfo["shapeid"] = shape_set['id']
        self.modelinfo['nodecon'] =  len(shape_set['nodes'])
        self.modelinfo['type_shape'] = shape_set["key"]
        self.modelinfo["elemid"] = int(f'{elem_set["id"]}{shape_set["id"]}') #elem_set['id']

        self.modelinfo['nnode'] = len(self.model.coord)
        self.modelinfo['nelem'] = len(self.model.inci)
        self.modelinfo['fulldofs'] = len(elem_set["dofs"]['d']) * len(self.model.coord)
        self.modelinfo['elemvol'] = newAnalysis.getElementVolume(self, self.modelinfo['inci'], self.modelinfo['coord'], self.modelinfo['tabgeo'], self.modelinfo['intgauss'])

    def Physic(self, physicdata):

        Loads, BoundCond = newAnalysis.setDomain(physicdata)
                
        # Domain = NewAnalysis.setDomain(physicdata)
        
        # if "COUPLING" in physicdata.keys():
        #     coupling = NewAnalysis.setCoupling(physicdata)
        #     Domain.coupling = coupling
        
        self.physic = SetPhysics(Loads, BoundCond)
        self.physic.physicdata = physicdata
        self.modelinfo['domain'] = physicdata['DOMAIN']

        try:
              self.modelinfo["forces"] = newAnalysis.getLoadApply(self)
        except:
            pass
        
        try:
             self.modelinfo["constrains"] = newAnalysis.getBCApply(self)
        except:
            pass
        
        # self.loadaply = FEANewAnalysis.getLoadApply(self)
        # self.constrains = FEANewAnalysis.getBCApply(self)
        
    def Assembly(self):

        inci = self.modelinfo['inci']
        coord = self.modelinfo['coord']
        tabmat = self.modelinfo['tabmat']
        tabgeo = self.modelinfo['tabgeo']
        intgauss = self.modelinfo['intgauss']       
        
        loadaply = newAnalysis.getLoadApply(self)            
        matrix = newAnalysis.getGlobalMatrix(self, inci, coord, tabmat, tabgeo, intgauss) #self.solver.getMatrixAssembler(self.model, inci, coord, tabmat, tabgeo)
                            
        try:
            addSpring = np.where(loadaply[:,1]==16)
            addMass = np.where(loadaply[:,1]==15)
            if addSpring or addMass:
                addLoad = loadaply[addSpring, :][0]
                print(len(addLoad))
                for ii in range(len(addLoad)):
                    A_add = addLoad[ii, 2] * np.array([[1, -1], [-1, 1]])
                    loc = np.array([int(self.modelinfo['dofe']*addLoad[ii, 0]-(self.modelinfo['dofe'])),
                                    int(self.modelinfo['dofe']*addLoad[ii, 0]-(self.modelinfo['dofe'] - 1))])
                    matrix['stiffness'] = self.solver.addMatrix(matrix['stiffness'], A_add, loc)
                    
            if addMass:
                addLoad = loadaply[addMass, :][0]
                for ii in range(len(addLoad)):
                    A_add = addLoad[ii, 2] * np.array([[1, -1], [-1, 1]])
                    loc = np.array([int(self.modelinfo['dofe']*addLoad[ii, 0]-(self.modelinfo['dofe'])),
                                    int(self.modelinfo['dofe']*addLoad[ii, 0]-(self.modelinfo['dofe'] - 1))])
                    matrix['mass'] = self.solver.addMatrix(matrix['mass'], A_add, loc) 
        except:
            pass
                               
        forcelist = newAnalysis.getLoadArray(self, loadaply)
                                                            
        return matrix, forcelist
    
    def FEASolve(self, solvedata):

        # fulldofs = self.modelinfo['fulldofs']

        # self.modelinfo = dict()
        # self.modelinfo['coord'] = self.coord
        # self.modelinfo['regions'] = self.regions

        assembly, forcelist = newAnalysis.Assembly(self)
        
        constrains = newAnalysis.getBCApply(self)
        freedof, fixedof = newAnalysis.getConstrains(self, constrains)
        # forcelist = FEANewAnalysis.getLoadArray(self)
        return self.solver.Solve(self.modelinfo['fulldofs'], assembly, forcelist, freedof, solvedata)
                

    # seting
    def setMesh(modeldata):
        
        set_mesh = dict()
        
        set_mesh = modeldata['MESH']
        # set_mesh['type'] = modeldata['MESH']['TYPE'] 
        set_mesh['shape'] = modeldata['ELEMENT']['SHAPE']
    
        return getMesh(set_mesh)
    
    def setShape(modeldata):
        
        set_shape = dict()
        set_shape['type'] = modeldata['ELEMENT']['SHAPE']
        
        return getShape(set_shape)
    
    def setElement(modeldata):
        
        set_element = dict()
        set_element['type'] = modeldata['ELEMENT']['TYPE']
        
        return getElement(set_element)
    
    def setMaterial(modeldata):
        
        set_material = dict()
        
        set_material['mat'] = modeldata['MATERIAL']['MAT']
        set_material['type'] = modeldata['MATERIAL']['TYPE']
        
        return getMaterial(set_material)
    
    def setGeometry(modeldata):
        
        set_geometry = dict()
        set_geometry['geo'] = modeldata['GEOMETRY']['GEO']
        
        return getGeometry(set_geometry)
            
    
    def setDomain(physicdata):
        if physicdata['DOMAIN'] == 'structural':
            from myfempy.felib.physic.bcstruct import BoundCondStruct
            from myfempy.felib.physic.loadstruct import LoadStructural
            return LoadStructural, BoundCondStruct
        else:
            pass 
    

    def setSolution(self, solvedata):
        
        # self.inci = FEANewAnalysis.getInci(self)
        # self.coord = FEANewAnalysis.getCoord(self)
        # self.tabmat = FEANewAnalysis.getTabmat(self)
        # self.tabgeo = FEANewAnalysis.getTabgeo(self)
        # self.intgauss = FEANewAnalysis.getIntGauss(self)
        
        # elem_set = self.model.element.getElementSet()
        # nodedof = len(elem_set["dofs"])
        self.solver.fulldofs = (self.modelinfo['dofs']) * len(self.modelinfo['nnode'])
        self.solver.solvedata = solvedata


    # geting    
    def getInci(self):
        return self.model.getInci(self.model.modeldata)
    
    def getCoord(self):
        return self.model.getCoord(self.model.modeldata)
    
    def getTabmat(self):
        return self.model.getTabMat(self.model.modeldata)
        
    def getTabgeo(self):
        return self.model.getTabGeo(self.model.modeldata)
    
    def getIntGauss(self):
        return self.model.setIntGauss(self.model.modeldata)
    
    def getElementVolume(self, inci, coord, tabgeo, intgauss):
        vol = np.zeros((len(inci)))
        for ee in range(len(inci)):
            vol[ee] = self.model.element.getElementVolume(self.model, inci, coord, tabgeo, intgauss, ee)
        return vol

    def getElemStifLinearMat(self, inci, coord, tabmat, tabgeo, intgauss, element_number):        
        return self.model.element.getStifLinearMat(self.model, inci, coord, tabmat, tabgeo, intgauss, element_number)
         
    def getElemMassConsistentMat(self, inci, coord, tabmat, tabgeo, intgauss, element_number):        
        return self.model.element.getMassConsistentMat(self.model, inci, coord, tabmat, tabgeo, intgauss, element_number)
                     
    def getGlobalMatrix(self, inci, coord, tabmat, tabgeo, intgauss):
        return self.solver.getMatrixAssembler(self.model, inci, coord, tabmat, tabgeo, intgauss)
    
    def getForceList(self):
        return self.physic.getForceList(self.physic.physicdata)
    
    def getBoundCondList(self):
        return self.physic.getBoundCondList(self.physic.physicdata)
    
    def getLoadApply(self):
        # coord = FEANewAnalysis.getCoord(self)
        return self.physic.getLoadApply(self.physic.physicdata, self.modelinfo)
    
    def getBCApply(self):
        # coord = FEANewAnalysis.getCoord(self)
        return self.physic.getBoundCondApply(self.physic.physicdata, self.modelinfo)
        
    def getConstrains(self, constrains):
        # coord = FEANewAnalysis.getCoord(self)
        
        # constrains = FEANewAnalysis.getBCApply(self)
        
        nodetot = len(self.modelinfo['coord'])
        
        # elem_set = self.model.element.getElementSet()
        # nodedof = len(elem_set["dofs"])
        
        freedof, fixedof = self.solver.getConstrains(constrains, nodetot, self.modelinfo['nodedof'])
        
        return freedof, fixedof
    
    def getLoadArray(self, loadaply):
        # coord = FEANewAnalysis.getCoord(self)
                
        nodetot = len(self.modelinfo['coord'])
        
        # elem_set = self.model.element.getElementSet()
        # nodedof = len(elem_set["dofs"])
                
        # loadaply = FEANewAnalysis.getLoadApply(self)
                
        loadvec = self.solver.getLoadAssembler(loadaply, nodetot, self.modelinfo['nodedof'])
        
        return loadvec
  


    def PreviewAnalysis(self, previewdata):

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
        
        preview_plot(previewdata, self.modelinfo)
     

    def PostProcess(self, postprocdata):

        postporc_result = []
        # try:
        if "COMPUTER" in postprocdata.keys(): 
            postporc_result = setPostProcess.getCompute(self, postprocdata)
        
        if "TRACKER" in postprocdata.keys():
            setPostProcess.getTracker(self, postprocdata, postporc_result)
        # except:
        #     pass
        return postporc_result

