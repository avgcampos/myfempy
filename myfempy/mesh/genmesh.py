# -*- coding: utf-8 -*-
"""
========================================================================
~~~ MODULO DE SIMULACAO ESTRUTURAL PELO METODO DOS ELEMENTOS FINITOS ~~~
       	                    __                                
       	 _ __ ___   _   _  / _|  ___  _ __ ___   _ __   _   _ 
       	| '_ ` _ \ | | | || |_  / _ \| '_ ` _ \ | '_ \ | | | |
       	| | | | | || |_| ||  _||  __/| | | | | || |_) || |_| |
       	|_| |_| |_| \__, ||_|   \___||_| |_| |_|| .__/  \__, |
       	            |___/                       |_|     |___/ 

~~~      Mechanical studY with Finite Element Method in PYthon       ~~~
~~~                PROGRAMA DE ANÁLISE COMPUTACIONAL                 ~~~
~~~              copyright @ 2022, all rights reserved               ~~~
========================================================================
"""

#%% READ INPUTDATA FROM USER PATH
import numpy as np
import os
from myfempy.felib.felemset import get_elemset
from myfempy.felib.crossec import sec_def
from myfempy.felib.materset import mat_def, mat_beh
from myfempy.felib.crossec import sect_prop
from myfempy.felib.physicset import gen_force, gen_bound
from myfempy.felib.physics.loadsconstr import get_forces, get_constrain
from myfempy.felib.physics.getnode import nodes_from_regions, search_nodexyz
from myfempy.tools.tools import print_console
from myfempy.tools.path import create_user_path

class MeshSet:
    
    def mesh2elem_key(meshtype):
        
        l = {'line2':['truss21', 'beam21', 'frame21', 'frame22'],
             'tria3':['plane31'],
             'quad4':['plane41'],
             'hexa8':['solid81'],
             'tetr4':['solid41'],
            }
        
        return l[meshtype]
    
    
    # def mesh_dim(meshdim):
        
    #     m = {'line2':['edge'],
    #          'tria3':['surf'],
    #          'quad4':['surf'],
    #          'hexa8':['volu'],
    #          'tetr4':['volu'],
    #         }
        
    #     return m[meshdim]
    
    
    def get_coord(nodelist):
        
        nnod = len(nodelist)
        coord = np.zeros((nnod,4))
        for ii in range(0,nnod):
            coord[ii,:] = np.array(nodelist[ii][:])
        
        return coord
       
        
    def get_inci(elemlist, mat_lib, geo_lib, regions):
        
        MAXCONECELM = int(8)
        inci = [[None]*(1+3+MAXCONECELM)]
        # conec_elm = np.zeros((nelem,MAXCONECELM+1))
        # prop_elm = np.zeros((nelem,3))
        # mesh_type_list = dict()
        # contelm = int(0)
        # for dm in range(len(elemlist)):
        
        nelem = len(elemlist)
        conec_elm = np.zeros((nelem,MAXCONECELM+1))
        prop_elm = np.zeros((nelem,3))
        mesh_type_list = dict()        
        contelm = int(0)
        
        # if 
                       
        for kk in range(nelem):
            contelm += 1
            conec_elm[kk, 0] = contelm #elemlist[dm][kk][0]
            nodes = np.array(elemlist[kk][4])
            if len(nodes) != MAXCONECELM:
                nodes = np.append(nodes,np.zeros(int(MAXCONECELM-len(nodes))),axis=0)
            conec_elm[kk,1:] = nodes
            
            keyelem = elemlist[kk][1]
            elem = get_elemset(keyelem)
            elemset = elem.elemset()
            prop_elm[kk,0] = int(elemset['id'])
            mesh_type_list[keyelem] = [elemset['id'], len(elemset['dofs']), len(elemset['nnodes']), len(elemset['tensor'])]
            
            prop_elm[kk,1] = mat_lib[elemlist[kk][2]]
            prop_elm[kk,2] = geo_lib[elemlist[kk][3]]
    
        inci = np.concatenate((conec_elm[:,0][:, np.newaxis], prop_elm, conec_elm[:,1:]),axis=1)
        # inci.extend(inci_dm)
            
        return inci, mesh_type_list


    def get_tabmat(matlist):
        
        nmat = len(matlist)
        mat_lib = dict()
        mat_prop = dict()
        key_mat_list = {'EXX':'EXX','VXX':'VXX','GXX':'GXX',
                        'EYY':'EYY','VYY':'VYY','GYY':'GYY',
                        'RHO':'RHO',
                        'STIF':'STIF', 'DAMP':'DAMP'}   
        
        
        tabmat = np.zeros((nmat,len(key_mat_list)+2))
        for mm in range(nmat):
            mat_lib[matlist[mm]["NAME"]] = mm+1
            
           
            for pp in range(len(key_mat_list)):
                key = list(key_mat_list)[pp]
                if key in matlist[mm].keys():
                    mat_prop[key] = matlist[mm][key]
                else:
                    mat_prop[key] = 0.0
                    
            # mat_prop["EXX"] = matlist[mm]["EXX"]
            # mat_prop["VXX"] = matlist[mm]["VXX"]
            # mat_prop["GXX"] = matlist[mm]["GXX"]
            # mat_prop["EYY"] = matlist[mm]["EYY"]
            # mat_prop["VYY"] = matlist[mm]["VYY"]
            # mat_prop["GYY"] = matlist[mm]["GYY"]
            # mat_prop["RHO"] = matlist[mm]["RHO"]
            
            idmat = mat_def(matlist[mm]["MAT"])
            iddef = mat_beh(matlist[mm]["DEF"])        
            tabmat[mm,:] = [mat_prop["EXX"],
                            mat_prop["VXX"],
                            mat_prop["GXX"],
                            mat_prop["EYY"],
                            mat_prop["VYY"],
                            mat_prop["GYY"],
                            mat_prop["RHO"],
                            mat_prop["STIF"],
                            mat_prop["DAMP"],
                            idmat,
                            iddef,
                            ]
    
        return tabmat, mat_lib
    
    
    def get_tabgeo(geolist):
        
        ngeo = len(geolist)
        geo_lib = dict()
        geo_prop = dict()
        key_geo_list = {'AREACS':'AREACS','INERYY':'INERYY','INERZZ':'INERZZ',
                        'INERXX':'INERXX','THICKN':'THICKN',
                        'b':'b', 'h':'h', 't':'t', 'd':'d'
                        }
        
        tabgeo = np.zeros((ngeo,len(key_geo_list)+1))
        for gg in range(ngeo):
            geo_lib[geolist[gg]["NAME"]] = gg+1
            
            if 'SEC' in geolist[gg].keys():
                if "THICKN" in geolist[gg].keys():
                    Thk = geolist[gg]["THICKN"]
                else:
                    Thk = 0.0
                
                b = geolist[gg]["DIM"]['b']
                h = geolist[gg]["DIM"]['h']
                t = geolist[gg]["DIM"]['t']
                d = geolist[gg]["DIM"]['d']
                A, Izz, Iyy, Jxx = sect_prop(geolist[gg]["SEC"], geolist[gg]["DIM"])
                
                idsec = sec_def(geolist[gg]["SEC"])
                tabgeo[gg,:] = [A,
                                Izz,
                                Iyy,
                                Jxx,
                                Thk,
                                b,
                                h,
                                t,
                                d,
                                idsec,
                                ]
                
            else:
                for pp in range(len(key_geo_list)):
                    key = list(key_geo_list)[pp]
                    if key in geolist[gg].keys():
                        geo_prop[key] = geolist[gg][key]
                    else:
                        geo_prop[key] = 0.0
                
                idsecdef=int(0)
                tabgeo[gg,:] = [geo_prop["AREACS"],
                                geo_prop["INERZZ"],
                                geo_prop["INERYY"],
                                geo_prop["INERXX"],
                                geo_prop["THICKN"],
                                geo_prop["b"],
                                geo_prop["h"],
                                geo_prop["t"],
                                geo_prop["d"],
                                idsecdef,
                                ]
        
        return tabgeo, geo_lib


class MeshGen:
       
    def get_data_mesh(meshdata):
         
        mesh = str()
        regions = []
        elemlist = [[None]*5]
        nodelist = [[None]*4]
        
        if 'LEGACY' in meshdata.keys():
            mesh = 'legacyON'
            
            if meshdata['LEGACY']['mesh'] == 'line2':
                from myfempy.mesh.legacy import get_legacy_line2
                conec, coord = get_legacy_line2(meshdata['LEGACY'])
                l = MeshSet.mesh2elem_key(meshdata['LEGACY']['mesh'])
                index = l.index(meshdata['LEGACY']['elem'])
                elemtype = l[index]
            
            elif meshdata['LEGACY']['mesh'] == 'tria3':
                from myfempy.mesh.legacy import get_legacy_tria3
                conec, coord = get_legacy_tria3(meshdata['LEGACY'])
                l = MeshSet.mesh2elem_key(meshdata['LEGACY']['mesh'])
                index = l.index(meshdata['LEGACY']['elem'])
                elemtype = l[index]
            
            elif meshdata['LEGACY']['mesh'] == 'quad4':
                from myfempy.mesh.legacy import get_legacy_quad4
                conec, coord = get_legacy_quad4(meshdata['LEGACY'])
                l = MeshSet.mesh2elem_key(meshdata['LEGACY']['mesh'])
                index = l.index(meshdata['LEGACY']['elem'])
                elemtype = l[index]
            
            else:
                pass
           
            elem = [[None]*5]
            for ee in range(len(conec)):
                elem.append([int(conec[ee,0]), elemtype, meshdata["PROPMAT"][0]["NAME"], meshdata["PROPGEO"][0]["NAME"], conec[ee,1:].astype(int).tolist()])
            
            elem = elem[1::][::]
                
            nodes = [[None]*4]
            for nn in range(len(coord)):
                nodes.append([int(coord[nn,0]), coord[nn,1], coord[nn,2], coord[nn,3]])
            nodes = nodes[1::][::]
            
            elemlist.extend(elem)
            nodelist.extend(nodes)
            
        elif 'GMSH' in meshdata.keys():
            mesh = 'gmshON'
            
            from myfempy.mesh.gmsh import get_gmsh_geo,  get_gmsh_msh
            from myfempy.io.iomsh import convert_from_msh1
            from myfempy.io.iomsh import gmsh_elm_type
            from myfempy.felib.felemset import get_elemset
            
            path = os.getcwd()
            meshdata['GMSH']['filename'] = path+'/'+meshdata['GMSH']['filename']
            
            if 'meshimport' in meshdata['GMSH'].keys():
                conec, nodes = convert_from_msh1(meshdata['GMSH']['meshimport']['object'])
            else:
    
                get_gmsh_geo(meshdata)
                get_gmsh_msh(meshdata)
                conec, nodes = convert_from_msh1(meshdata['GMSH']['filename'])
             
            l = MeshSet.mesh2elem_key(meshdata['GMSH']['meshconfig']['mesh'])
            index = l.index(meshdata['GMSH']['meshconfig']['elem'])
            elemtype = l[index]
            
            # mesh_space = MeshSet.mesh_dim(meshdata['GMSH']['meshconfig']['mesh'])
            
            fe = get_elemset(elemtype)
            
            eset = fe.elemset()
            elemid = eset['id']
            
            meshid = gmsh_elm_type(str(elemid))
            nconec = int(len(eset['nnodes']))
                 
            elem = [[None]*5]   
            contelm = 0
            reglist = dict()
            reglist['reglist'] = []
            
            contptn = 0
            contedg = 0
            contsfr = 0
            contvlr = 0
            
            vlr0=0
            sfr0=0
            edg0=0
            ptn0=0
            
            for ee in range(len(conec)):    
                if int(conec[ee][1]) == meshid:
                    contelm += 1
                    elem.append([contelm, meshdata['GMSH']['meshconfig']['elem'], meshdata["PROPMAT"][0]["NAME"], meshdata["PROPGEO"][0]["NAME"], (np.array(conec[ee][5:])).astype(int).tolist()])
 
                # else:
                if int(conec[ee][1]) == 15:
                    ptn1 = int(conec[ee][3])
                    if ptn1 != ptn0:
                        contptn+=1
                    reglist['reglist'].append({'type':'point', 'list': str(contptn), 'nodes':(np.array(conec[ee][5:])).astype(int).tolist()})
                    ptn0 = ptn1
                    
                    
                elif int(conec[ee][1]) == 1:
                    edg1 = int(conec[ee][3])
                    if edg1 != edg0:
                        contedg+=1
                    reglist['reglist'].append({'type':'edge', 'list': str(contedg), 'nodes':(np.array(conec[ee][5:])).astype(int).tolist()})
                    edg0 = edg1
                    
                elif (int(conec[ee][1]) == 2) or (int(conec[ee][1]) == 3):
                    sfr1 = int(conec[ee][3])
                    if sfr1 != sfr0:
                        contsfr+=1
                    reglist['reglist'].append({'type':'surf', 'list': str(contsfr), 'nodes':(np.array(conec[ee][5:])).astype(int).tolist()})
                    sfr0 = sfr1
                        
                    
                # elif (int(conec[ee][1]) == 2) or (int(conec[ee][1]) == 3):
                #     sfr1 = int(conec[ee][3])
                #     if sfr1 != sfr0:
                #         contsfr+=1
                #     reglist['reglist'].append({'type':'volm', 'list': str(contsfr), 'nodes':(np.array(conec[ee][5:])).astype(int).tolist()})
                #     sfr0 = sfr1
                                
            elem = elem[1::][::]
            
            keys = ["nodes"]
            
            list_point = (np.arange(1,contptn+1).astype(str)).tolist()
            pointlist = { name.capitalize():{ key:[] for key in keys} for name in list_point}
            
            list_edge = (np.arange(1,contedg+1).astype(str)).tolist()
            edgelist = { name.capitalize():{ key:[] for key in keys} for name in list_edge}
            
            list_surf = (np.arange(1,contsfr+1).astype(str)).tolist()
            surflist = { name.capitalize():{ key:[] for key in keys} for name in list_surf}
            
            for rr in range(len(reglist['reglist'])):
                reg = reglist['reglist'][rr]
                if reg['type'] == 'point':
                    pointlist[reg['list']]['nodes'].extend(reg['nodes'])
                
                elif reg['type'] == 'edge':
                    edgelist[reg['list']]['nodes'].extend(reg['nodes'])         
                    
                elif reg['type'] == 'surf':
                    surflist[reg['list']]['nodes'].extend(reg['nodes'])
                    
                else:
                    pass
            
            regionlist = dict()
            regionlist['point'] = pointlist
            regionlist['edge'] = edgelist
            regionlist['surf'] = surflist
            
            regions = nodes_from_regions(regionlist)
                        
            elemlist.extend(elem)
            nodelist.extend(nodes)
            
        elemlist = elemlist[1::][::]
        nodelist = nodelist[1::][::]
        
        matlist = meshdata["PROPMAT"]
        geolist = meshdata["PROPGEO"]
            
        tabmat, mat_lib = MeshSet.get_tabmat(matlist)
        tabgeo, geo_lib = MeshSet.get_tabgeo(geolist)
        
        coord = MeshSet.get_coord(nodelist)
        
        if "ADD121" in meshdata.keys():
            if (mesh == 'legacyON') or (mesh == 'gmshON'):
                nodefind = [None]
                for nn in range(len(meshdata["NODELIST"])):
                    node = search_nodexyz(meshdata["NODELIST"][nn][1], meshdata["NODELIST"][nn][2], meshdata["NODELIST"][nn][3], coord, 2E-3)
                    if len(node) == 0:
                        nodefind.extend([len(coord)+1])
                        coord = np.append(coord,np.array([[len(coord)+1, meshdata["NODELIST"][nn][1], meshdata["NODELIST"][nn][2] ,meshdata["NODELIST"][nn][3]]]),axis=0)
                    else:
                        nodefind.extend(node)
                
                nodefind = nodefind[1::][::]
                elemlist.append([len(elemlist), meshdata["ELEMLIST"][0][1], meshdata["ELEMLIST"][0][2], meshdata["ELEMLIST"][0][3],
                                 nodefind])
                # inci, mesh_type_list = get_inci(elemlist, mat_lib, geo_lib)
            else:
                
                # regionlist = dict()
                # regionlist['edge'] = []
                
                elemlist = meshdata["ADD121"]
                nodelist = meshdata["NODELIST"]
                
                coord = MeshSet.get_coord(nodelist)
                # elemlist = elemlist[1::][::]
                # nodelist = nodelist[1::][::]
            
        else:
            pass
        
        # coord = get_coord(nodelist)
        inci, mesh_type_list = MeshSet.get_inci(elemlist, mat_lib, geo_lib, regions)
    
        return mesh_type_list, inci, coord, tabmat, tabgeo, regions



class ModelGen:
    
    def get_quadra(quadrature):
        
        if quadrature['meth'] == 'gaussian':
            quadra = [1, quadrature['npp']]
            return quadra
        
        elif quadrature['meth'] == 'no_interpol':
            quadra = [0, 1]
            return quadra

    def get_model(meshdata):
        '''meshdata running 
        
        
        Parameters
        ----------
        meshdata : dict()
            DESCRIPTION.
    
        Returns
        -------
        modelinfo : dict()
            DESCRIPTION.
    
        '''

        print_console('pre')
        print_console('mesh')

        mesh_type_list, inci, coord, tabmat, tabgeo, regions = MeshGen.get_data_mesh(meshdata)
        
        modelinfo = dict()
        modelinfo['elemid'] = []
        modelinfo["nodedof"] = []
        modelinfo["nodecon"] = []
        modelinfo["ntensor"] = []
        modelinfo['regions'] = regions
        # modelinfo["fulldof"] = []
        for dd in range(len(mesh_type_list)):
            keyelem = list(mesh_type_list)[dd]
            # elemid = mesh_type_list[keyelem][0]
            # ndof = mesh_type_list[keyelem][1]
            # nconec = mesh_type_list[keyelem][2]
            modelinfo['elemid'].append(mesh_type_list[keyelem][0])
            modelinfo["nodedof"].append(mesh_type_list[keyelem][1])
            modelinfo["nodecon"].append(mesh_type_list[keyelem][2])
            modelinfo["ntensor"].append(mesh_type_list[keyelem][3])
            # modelinfo["nnode"] = len(coord)
            # modelinfo["nelem"] = len(inci)
            # modelinfo["fulldof"].append(ndof*len(coord))
            modelinfo["inci"] = inci
            modelinfo["coord"] = coord
            modelinfo["tabmat"] = tabmat
            modelinfo["tabgeo"] = tabgeo
            
            
        if "FORCES" in meshdata.keys():
            flist = gen_force(meshdata["FORCES"])
            forces = get_forces(modelinfo, flist)
            modelinfo['forces'] = forces           
        
        if "BOUNDCOND" in meshdata.keys():
            blist = gen_bound(meshdata["BOUNDCOND"])
            constrains = get_constrain(modelinfo, blist)
            modelinfo['constrains'] = constrains
        
        if "QUADRATURE" in meshdata.keys():
            modelinfo['quadra'] = ModelGen.get_quadra(meshdata["QUADRATURE"])
        
        else:
            meshdata["QUADRATURE"] = {'meth':'no_interpol','npp':1}
            modelinfo['quadra'] = ModelGen.get_quadra(meshdata["QUADRATURE"])
        
              
        return modelinfo
    
    

    