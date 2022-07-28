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

import numpy as np
# import scipy.sparse as sp
# import matplotlib.pyplot as plt
# from myfempy.plots.plotter import postproc_show
# from myfempy.felib.modelset import get_element
# from myfempy.io.io import convert2vtk
# from myfempy.felib.loadsconstr import search_nodexyz
# from myfempy.plots.plotxy import tracker_plot


import os
# import numpy as np
import vtk
# import vedo as vd
from myfempy.io.iovtk_old import convert_to_vtk
from myfempy.plots.physics import view_listforce, view_bondcond_point, view_beam_crossSection, view_text_point
from myfempy.tools.tools import get_version
from myfempy.plots.meshquality import MeshProp


def preview_plot(previewset, modelinfo):
    
    path = os.getcwd()
    
    plotdata = dict()
    plotdata['coord'] = modelinfo['coord']
    plotdata['inci'] = modelinfo['inci']
    plotdata['nodecon'] = modelinfo['nodecon']
    plotdata['elemid']  = modelinfo['elemid']
    plotdata['filename'] = path+'/'+previewset['RENDER']['filename']
    plotdata['title'] = ['UNDEFORM_MESH']
    plotdata['solution'] = np.ones((len(modelinfo['inci']),1))
    plotdata['average'] = False
    
    convert_to_vtk(plotdata)

    if "scale" in previewset['RENDER'].keys():
        previewset['RENDER']["scale"] = (previewset['RENDER']["scale"]/100)*max([max(abs(modelinfo["coord"][:,1])), max(abs(modelinfo["coord"][:,2])), max(abs(modelinfo["coord"][:,3]))])
                   
    else:
        previewset['RENDER']["scale"] = 1
                     
        
    if 'lines' in previewset['RENDER'].keys():
        pass
    else:
        previewset['RENDER']['lines'] = True
        
        
    if 'plottags' in previewset['RENDER'].keys():
        
        if 'point' in previewset['RENDER']['plottags'].keys():
            if previewset['RENDER']['plottags']['point'] == True:
                previewset["regions"] = modelinfo["regions"][0]
            else:
                previewset["regions"] = [[],[]]
        
        
        if 'edge' in previewset['RENDER']['plottags'].keys():
            if previewset['RENDER']['plottags']['edge'] == True:
                previewset["regions"] = modelinfo["regions"][1]
            else:
                previewset["regions"] = [[],[]]
            
        elif 'surf' in previewset['RENDER']['plottags'].keys():
            if previewset['RENDER']['plottags']['surf'] == True:
                previewset["regions"] = modelinfo["regions"][2]
            else:
                previewset["regions"] = [[],[]]
    else:
        previewset["regions"] = [[],[]]
      
        
    if 'cs' in previewset['RENDER'].keys():
        pass
    else:
        previewset['RENDER']['cs'] = False
        
    
    previewset["coord"] = modelinfo["coord"]
    previewset["inci"] = modelinfo['inci']
    previewset["nnode"] = len(modelinfo["inci"])
    previewset["nodecon"] = modelinfo['nodecon'][0]
    # previewset["nnode"] = len(modelinfo["coord"])
    
    
    if "forces" in modelinfo.keys():
        previewset["forces"] = modelinfo["forces"]
        
    if "constrains" in modelinfo.keys():
        previewset["constrains"] = modelinfo["constrains"]
        
    else:
        pass
        # previewset["forces"] = []
        # previewset["constrains"] =[]

    previewset['tabcs'] = dict()
    previewset['tabcs']['typSection'] = []
    previewset['tabcs']['dimSection'] = []
    
    for gg in range(len(modelinfo["tabgeo"])):
        previewset['tabcs']['typSection'].append(int(modelinfo["tabgeo"][gg][-1]))
        previewset['tabcs']['dimSection'].append(modelinfo["tabgeo"][gg][5:9])

    build_preview(previewset)

                
    if previewset['QUALITY']['show'] == False:
        pass
    else:
        
        if 'lines' in previewset['QUALITY'].keys():
            pass
        else:
            previewset['QUALITY']['lines'] = False     
            
            
        if "scale" in previewset['QUALITY'].keys():
            previewset['QUALITY']["scale"] = (previewset['QUALITY']["scale"]/100)*max([max(abs(modelinfo["coord"][:,1])), max(abs(modelinfo["coord"][:,2])), max(abs(modelinfo["coord"][:,3]))])
                       
        else:
           previewset['QUALITY']["scale"] = 1    
            
                       
        # previewset['quality']['posx'] = max(abs(modelinfo['coord'][:,1]))+0.5*max(abs(modelinfo['coord'][:,1]))
        # previewset['quality']['posy'] = 0 #max(abs(modelinfo['coord'][:,2]))
        
        mesh = MeshProp(previewset)
        mesh.mesh_quality()

             
    if previewset['LABELS']['show'] == False:
        pass
    else:
        
        if 'lines' in previewset['LABELS'].keys():
            pass
        else:
            previewset['LABELS']['lines'] = False   
            
            
        if "scale" in previewset['LABELS'].keys():
            previewset['LABELS']["scale"] = (previewset['LABELS']["scale"]/100)*max([max(abs(modelinfo["coord"][:,1])), max(abs(modelinfo["coord"][:,2])), max(abs(modelinfo["coord"][:,3]))])
                       
        else:
           previewset['LABELS']["scale"] = 1     

        mesh = MeshProp(previewset)
        mesh.mesh_numbering()
       
           

#-----------------------------------------------------------------------------#
def build_preview(previewset):
    path = os.getcwd()
    # OS Linux
    # file_name = (previewset['file_name']+'.vtk')
    file_name = str(path+'/'+previewset['RENDER']['filename']+'.vtk')
    # file_savePNG = str(os.getcwd() + '/' +  previewset['file_savePNG'])
    # file_saveSTL = str(os.getcwd() + '/' +  previewset['file_saveSTL'])
    
    renderer = vtk.vtkRenderer()
    renderer_window = vtk.vtkRenderWindow()
    renderer_window.SetSize(600, 480)
    renderer.SetBackground(0.0, 0.0, 0.0) # Set background to gray -> 0.7,0.7,0.7
   
    # Read the source file.
    reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName(file_name)
    reader.Update()  # Needed because of GetScalarRange
    output = reader.GetOutput()
    output_port = reader.GetOutputPort()
    scalar_range = output.GetScalarRange()
        
    # Create the mapper that corresponds the objects of the vtk file
    # into graphics elements
    lut = vtk.vtkLookupTable()
    lut.SetNumberOfColors(8)
    lut.SetHueRange(0.6, 1) #jet color
    lut.Build()    
    
    mapper = vtk.vtkDataSetMapper()
    # mapper.vtkPolyDataMapper()
    mapper.SetLookupTable(lut)
    mapper.SetInputConnection(output_port)
    mapper.SetScalarRange(scalar_range)    
             
    
    #--------------------------------------------
    scala_view = previewset['RENDER']['scale']
    # view physics elements  
    
    
    dimfrlist = 0 
    if "forces" in previewset.keys():
        dimfrlist = previewset['forces'].shape[0]
        for num_lf in range(dimfrlist):
            exec(f'fr_point_actor_cone1_{num_lf},fr_point_actor_cone2_{num_lf} = view_listforce(previewset["coord"],previewset["forces"],scala_view,num_lf)')
                  
        
    dimbclist = 0
    if "constrains" in previewset.keys():
        dimbclist = previewset['constrains'].shape[0]
        for num_bc in range(dimbclist):
            exec(f'bc_point_actor_cone_{num_bc}, bc_point_actor_tdof_{num_bc} = view_bondcond_point(previewset["coord"],previewset["constrains"],scala_view,num_bc)')
          
    
    objs = 0 
    coordMax = [max(previewset["coord"][:,1]), max(previewset["coord"][:,2]), max(previewset["coord"][:,3])]
    # for num_rg in range(len(previewset["regions"])):
    for num_objs in range(len(previewset["regions"][1])):
        text = [previewset["regions"][0], str(previewset["regions"][1][num_objs][0])]
        
        coord = [sum(previewset["coord"][previewset["regions"][1][num_objs][1]-1,1])/len(previewset["coord"][previewset["regions"][1][num_objs][1]-1,1]),
                 sum(previewset["coord"][previewset["regions"][1][num_objs][1]-1,2])/len(previewset["coord"][previewset["regions"][1][num_objs][1]-1,2]),
                 sum(previewset["coord"][previewset["regions"][1][num_objs][1]-1,3])/len(previewset["coord"][previewset["regions"][1][num_objs][1]-1,2])]
        
        
        # coord =  [max(previewset["coord"][previewset["regions"][1][num_objs][1]-1,1]),
        #           max(previewset["coord"][previewset["regions"][1][num_objs][1]-1,2]),
        #           max(previewset["coord"][previewset["regions"][1][num_objs][1]-1,3])]
        
        # for num_nodes in range(len(previewset["regions"][1][num_objs][1])):
            
            # node = previewset["regions"][1][num_objs][1][num_nodes]
        exec(f'bc_text_actor_{objs} = view_text_point(coord, coordMax, scala_view, text)')
        objs += 1 
                
    
    i_node = 0
    if previewset['RENDER']['cs'] == True:
        dimcs = len(previewset['tabcs']['typSection'])
        for num_bcs in range(dimcs):
            inci_bcs = previewset["inci"][np.where(previewset["inci"][:,3]==num_bcs+1),:][0]
            typSection = int(previewset['tabcs']['typSection'][num_bcs])
            dimSection = previewset['tabcs']['dimSection'][num_bcs]
            for node_bcs in range(len(inci_bcs)):
                noi = int(inci_bcs[node_bcs,4])
                noj = int(inci_bcs[node_bcs,5])
                noix = previewset["coord"][noi-1,1]
                noiy = previewset["coord"][noi-1,2]
                noiz = previewset["coord"][noi-1,3]
                nojx = previewset["coord"][noj-1,1]
                nojy = previewset["coord"][noj-1,2]
                nojz = previewset["coord"][noj-1,3]
                coord_bcs = [noix,noiy,noiz,nojx,nojy,nojz]
                exec(f'beam_extrude_actor_{i_node} = view_beam_crossSection(dimSection, typSection, coord_bcs)')
                i_node +=1
    
        
    # Create the Actor
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    
    actor.GetProperty().EdgeVisibilityOn()
    
    if previewset['RENDER']['lines'] == False:
        actor.GetProperty().EdgeVisibilityOff()
        # actor.GetProperty().SetLineWidth(0.9)
    else:
        actor.GetProperty().SetLineWidth(3.0)
            
    text_logo =  vtk.vtkTextActor()
    text_logo.SetInput('MYFEMPY v'+get_version()+' < PreProc--Model >\nPress "q" to continue... ')
    txtprop=text_logo.GetTextProperty()
    txtprop.SetFontFamilyToArial()
    txtprop.SetFontSize(20)
    # txtprop.BoldOn()
    # txtprop.ShadowOn()
    txtprop.SetColor(1, 1, 1)
    text_logo.SetDisplayPosition(10,430)
    
    colors = vtk.vtkNamedColors()

    backgroundColor = colors.GetColor3d("DarkSlateGray")
    actorColor = colors.GetColor3d("Tomato")
    axis1Color = colors.GetColor3d("Salmon")
    axis2Color = colors.GetColor3d("PaleGreen")
    axis3Color = colors.GetColor3d("LightSkyBlue")
    
    cubeAxesActor = vtk.vtkCubeAxesActor()
    # cubeAxesActor.SetFlyModeToStaticTriad()
    # cubeAxesActor.SetUseTextActor3D(1)
    cubeAxesActor.SetBounds(actor.GetBounds())
    cubeAxesActor.SetCamera(renderer.GetActiveCamera())
    # cubeAxesActor.GetTitleTextProperty(0).SetColor(axis1Color)
    cubeAxesActor.GetTitleTextProperty(0).SetFontSize(25)
    cubeAxesActor.GetLabelTextProperty(0).SetColor(axis1Color)
    cubeAxesActor.GetLabelTextProperty(1).SetColor(axis2Color)
    cubeAxesActor.GetLabelTextProperty(2).SetColor(axis3Color)
    
    # cubeAxesActor.DrawXGridlinesOn()
    # cubeAxesActor.DrawYGridlinesOn()
    # cubeAxesActor.DrawZGridlinesOn()
    # cubeAxesActor.SetGridLineLocation(cubeAxesActor.VTK_GRID_LINES_FURTHEST)        
    # cubeAxesActor.SetFlyModeToStaticEdges() 
    # cubeAxesActor.SetFlyModeToStaticTriad() 
    # cubeAxesActor.SetFlyModeToOuterEdges()
          
    
    # Create the Renderer
    renderer.AddActor(actor)
    renderer.AddActor(text_logo)
    renderer.AddActor(cubeAxesActor)
    
            
    for ipt in range(dimbclist):
        exec(f'renderer.AddActor(bc_point_actor_cone_{ipt})')
        exec(f'renderer.AddActor(bc_point_actor_tdof_{ipt})')
        # exec(f'renderer.AddActor(bc_text_actor_{ipt})')
                
    for bcs in range(i_node):
        exec(f'renderer.AddActor(beam_extrude_actor_{bcs})')
        
        
    for ff in range(dimfrlist):
        exec(f'renderer.AddActor(fr_point_actor_cone1_{ff})')
        exec(f'renderer.AddActor(fr_point_actor_cone2_{ff})')
            
    if 'plottags' in previewset['RENDER'].keys():  
        for txt in range(objs):
            exec(f'renderer.AddActor(bc_text_actor_{txt})')
    else:
        pass
        
    renderer_window.AddRenderer(renderer)
    renderer.ResetCamera()
    # renderer.GetActiveCamera().Elevation(120.0)
    # renderer.GetActiveCamera().Azimuth(30.0)

    if previewset['RENDER']['savepng'] == False:
        pass
    else: 
        # screenshot code:
        im = vtk.vtkWindowToImageFilter()
        writer = vtk.vtkPNGWriter() #vtkSTLWriter()#vtkVRMLExporter()#vtkPolyDataWriter()#
        # writer.SetFileTypeToBinary()
        # writer.SetInputData()
        im.SetInput(renderer_window)
        im.Update()
        writer.SetInputConnection(im.GetOutputPort())
        writer.SetFileName((path+'/'+previewset['RENDER']['filename']+'.png'))
        writer.Write()
    
        # # Write the extrusion to stl ------------------------------------------------------
        
        # surface_filter = vtkDataSetSurfaceFilter()
        # surface_filter.SetInputConnection(renderer_window.GetOutputPort())
        
        # triangle_filter = vtkTriangleFilter()
        # triangle_filter.SetInputConnection(surface_filter.GetOutputPort())
        
        # stlWriter = vtkSTLWriter()
        # stlWriter.SetFileTypeToASCII()
        # # stlWriter.SetRenderWindow(window)
        # stlWriter.SetInputConnection(triangle_filter.GetOutputPort())
        # stlWriter.SetFileName(file_saveSTL)
        # stlWriter.Write()
    
    if previewset['RENDER']['show'] == False:
        pass
    
    else:
        # # # Create the RendererWindowInteractor and display the vtk_file
        interactor = vtk.vtkRenderWindowInteractor()
        interactor.SetRenderWindow(renderer_window)
        # interactor.RemoveObservers('LeftButtonPressEvent')
        interactor.RemoveObservers('RightButtonPressEvent')
        interactor.Initialize()
        renderer_window.Render() 
        interactor.Start()
        interactor.GetRenderWindow().Finalize()