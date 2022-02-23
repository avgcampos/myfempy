# -*- coding: utf-8 -*-
"""
Created on Sun Nov  1 16:00:31 2020
@author: ANTONIO VINICIUS GARCIA CAMPOS
@version: teste beta
_______________________________________________________________________________
 ~~~~~~~~~~ MODULO DE SIMULACAO PELO METODO DOS ELEMENTOS FINITOS ~~~~~~~~~~

ESTE MODULO CRIA UMA JANELA PARA VISUALIZACAO DE MALHA E DOS
RESULTADOS COM ARQUIVO DE ENTRADA EXT. vtk    
_______________________________________________________________________________
"""
# # upload de pacotes auxiliares
import os
import numpy as np
from vtk import *
from datetime import date

# visualizador de resultados de malha version teste - beta
def view_analysis(file_dir,file_savePNG,file_saveSTL,screen_on,save_screen,inci,coord,frcApy_vet,bondCond_vet,project_name,title_win,myfempy_version,scala_view,tabdimSection,tabtypSection,screen_beam_cs):
    
    # OS Linux
    file_name = str(os.getcwd() + '/' +  file_dir)
    file_savePNG = str(os.getcwd() + '/' +  file_savePNG)
    file_saveSTL = str(os.getcwd() + '/' +  file_saveSTL)
    
    renderer = vtkRenderer()
    renderer_window = vtkRenderWindow()
    renderer_window.SetSize(1024,768)
    renderer.SetBackground(0.0, 0.0, 0.0) # Set background to gray -> 0.7,0.7,0.7
       
    
    # interactor = vtkRenderWindowInteractor()
    # interactor.SetRenderWindow(renderer_window)
    
    # Read the source file.
    reader = vtkUnstructuredGridReader()
    reader.SetFileName(file_name)
    reader.Update()  # Needed because of GetScalarRange
    output = reader.GetOutput()
    output_port = reader.GetOutputPort()
    scalar_range = output.GetScalarRange()
        
    # Create the mapper that corresponds the objects of the vtk file
    # into graphics elements
    lut = vtkLookupTable()
    lut.SetNumberOfColors(8)
    lut.SetHueRange(0.5, 1) # jet color
    lut.Build()    
    
    mapper = vtkDataSetMapper()
    # mapper.vtkPolyDataMapper()
    mapper.SetLookupTable(lut)
    mapper.SetInputConnection(output_port)
    mapper.SetScalarRange(scalar_range)    
    
    #--------------------------------------------
    # view physics elements    
    dimfrlist = frcApy_vet.shape[0]
    for num_lf in range(dimfrlist):
        exec(f'fr_point_actor_cone1_{num_lf},fr_point_actor_cone2_{num_lf}, fr_text_actor_{num_lf} = view_listforce(coord,frcApy_vet,scala_view,num_lf)')
                
    # cont_pt = 0
    # cont_eg = 0
    # cont_sf = 0
    dimbclist = len(bondCond_vet)
    for num_bc in range(dimbclist):
        exec(f'bc_point_actor_cone_{num_bc},bc_point_actor_tdof_{num_bc},bc_text_actor_{num_bc} = view_bondcond_point(coord,bondCond_vet,scala_view,num_bc)')
        
    i_node = 0
    if screen_beam_cs == 'true':
        for num_bcs in range(len(tabtypSection)):
            inci_bcs = inci[np.where(inci[:,3]==num_bcs+1),:][0]
            typSection = tabtypSection[num_bcs]
            dimSection = tabdimSection[num_bcs]
            for node_bcs in range(len(inci_bcs)):
                noi = int(inci_bcs[node_bcs,4])
                noj = int(inci_bcs[node_bcs,5])
                noix = coord[noi-1,1]
                noiy = coord[noi-1,2]
                noiz = coord[noi-1,3]
                nojx = coord[noj-1,1]
                nojy = coord[noj-1,2]
                nojz = coord[noj-1,3]
                coord_bcs = [noix,noiy,noiz,nojx,nojy,nojz]
                exec(f'beam_extrude_actor_{i_node} = view_beam_crossSection(dimSection,typSection,coord_bcs)')
                i_node +=1

        
    #--------------------------------------------    
    # Create the Actor
    actor = vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetLineWidth(1.0)
    actor.GetProperty().EdgeVisibilityOn()
    # actor.GetProperty().SetColor(0,0,0)
    
    # # create a text actor
    # txt1 = vtkTextActor()
    # txt1.SetInput("{mouse3 - scroll} - zoom / {mouse3 - click} - move screen \n {w} - wireframe style / {s} - surface style \n {q} - close visualization")
    # txtprop1=txt1.GetTextProperty()
    # txtprop1.SetFontFamilyToArial()
    # txtprop1.SetFontSize(14)
    # txtprop1.BoldOn()
    # txtprop1.SetColor(1,1,1)
    # txt1.SetDisplayPosition(10,642)
    
    # txt2 =  vtkTextActor()
    # txt2.SetInput('close the visualization to continue the analysis')
    # txtprop2=txt2.GetTextProperty()
    # txtprop2.SetFontFamilyToArial()
    # txtprop2.SetFontSize(14)
    # txtprop2.BoldOn()
    # txtprop2.SetColor(1,1,1)
    # txt2.SetDisplayPosition(10,622)
    
    
    text_logo =  vtkTextActor()
    text_logo.SetInput(myfempy_version)
    txtprop=text_logo.GetTextProperty()
    txtprop.SetFontFamilyToArial()
    txtprop.SetFontSize(40)
    txtprop.BoldOn()
    txtprop.ShadowOn()
    txtprop.SetColor(1,1,1)
    text_logo.SetDisplayPosition(10,718)
    # text_logo.SpecifiedToDisplay()

    
    text_project =  vtkTextActor()
    text_project.SetInput('JOB NAME: '+project_name+" | DATA: {"+str(date.today())+"} ")
    txtprop=text_project.GetTextProperty()
    txtprop.SetFontFamilyToArial()
    txtprop.SetFontSize(20)
    txtprop.BoldOn()
    txtprop.ItalicOn()
    txtprop.SetColor(1,1,1)
    text_project.SetDisplayPosition(10,658)  
    
    # text_dir =  vtkTextActor()
    # text_dir.SetInput('PROJECT DIR: '+file_name)
    # txtprop=text_dir.GetTextProperty()
    # txtprop.SetFontFamilyToArial()
    # txtprop.SetFontSize(14)
    # txtprop.BoldOn()
    # txtprop.ItalicOn()
    # txtprop.SetColor(1,1,1)
    # text_dir.SetDisplayPosition(20,660)    
    
    text_title =  vtkTextActor()
    text_title.SetInput(title_win)
    txtprop=text_title.GetTextProperty()
    txtprop.SetFontFamilyToArial()
    txtprop.SetFontSize(20)
    txtprop.BoldOn()
    # txtprop.ItalicOn()
    txtprop.SetColor(1,1,1)
    text_title.SetDisplayPosition(10,688)
    
    # Axes Origin XYZ
    axes = vtkAxesActor()
    # axes.SetScale(1, 1, 1)
    axes.SetShaftTypeToCylinder()
    axes.SetXAxisLabelText('x')
    axes.SetYAxisLabelText('y')
    axes.SetZAxisLabelText('z')
    # axes.vtkTextProperty().SetFontSize(10)
        
    transform =vtkTransform()
    transform.Scale(2.6*scala_view,2.6*scala_view,2.6*scala_view)
      
    tprop = axes.GetXAxisCaptionActor2D().GetCaptionTextProperty()
    tprop.ItalicOn()
    tprop.ShadowOn()
    tprop.SetFontFamilyToArial()
    tprop.SetColor(1,1,1)
    # # Use the same text properties on the other two axes.
    axes.GetYAxisCaptionActor2D().GetCaptionTextProperty().ShallowCopy(tprop)
    axes.GetZAxisCaptionActor2D().GetCaptionTextProperty().ShallowCopy(tprop)
    axes.SetUserTransform(transform)
    
    
    colors = vtkNamedColors()

    backgroundColor = colors.GetColor3d("DarkSlateGray")
    actorColor = colors.GetColor3d("Tomato")
    axis1Color = colors.GetColor3d("Salmon")
    axis2Color = colors.GetColor3d("PaleGreen")
    axis3Color = colors.GetColor3d("LightSkyBlue")
    
    cubeAxesActor = vtkCubeAxesActor()
    # cubeAxesActor.SetFlyModeToStaticTriad()
    # cubeAxesActor.SetUseTextActor3D(1)
    cubeAxesActor.SetBounds(actor.GetBounds())
    cubeAxesActor.SetCamera(renderer.GetActiveCamera())
    # cubeAxesActor.GetTitleTextProperty(0).SetColor(axis1Color)
    # cubeAxesActor.GetTitleTextProperty(0).SetFontSize(48)
    # cubeAxesActor.GetLabelTextProperty(0).SetColor(axis1Color)
    
    cubeAxesActor.DrawXGridlinesOn()
    cubeAxesActor.DrawYGridlinesOn()
    cubeAxesActor.DrawZGridlinesOn()
    # cubeAxesActor.SetGridLineLocation(cubeAxesActor.VTK_GRID_LINES_FURTHEST)        
    # cubeAxesActor.SetFlyModeToStaticEdges() 
    # cubeAxesActor.SetFlyModeToStaticTriad() 
    cubeAxesActor.SetFlyModeToOuterEdges()
    
    # Create the Renderer
    renderer.AddActor(actor)
    # renderer.AddActor(txt1)
    # renderer.AddActor(txt2)
    renderer.AddActor(axes)
    renderer.AddActor(text_project)
    renderer.AddActor(text_logo)
    renderer.AddActor(text_title)
    renderer.AddActor(cubeAxesActor)
    
    for ff in range(dimfrlist):
        exec(f'renderer.AddActor(fr_point_actor_cone1_{ff})')
        exec(f'renderer.AddActor(fr_point_actor_cone2_{ff})')
        exec(f'renderer.AddActor(fr_text_actor_{ff})')
        
    
    for ipt in range(dimbclist):
        exec(f'renderer.AddActor(bc_point_actor_cone_{ipt})')
        exec(f'renderer.AddActor(bc_point_actor_tdof_{ipt})')
        exec(f'renderer.AddActor(bc_text_actor_{ipt})')
        
        
    for bcs in range(i_node):
        exec(f'renderer.AddActor(beam_extrude_actor_{bcs})')
        
    # for ieg in range(cont_eg):
    #     exec(f'renderer.AddActor(bc_edge_actor_{2}{ieg})')
    #     exec(f'renderer.AddActor(bc_center_edge_actor_{2}{ieg})')
        
        
    renderer_window.AddRenderer(renderer)
    renderer.ResetCamera()
    renderer.GetActiveCamera().Elevation(45.0)
    renderer.GetActiveCamera().Azimuth(45.0)
    
    if save_screen == 'false':
        print(' ')
    else: 
        # screenshot code:
        im = vtkWindowToImageFilter()
        writer = vtkPNGWriter() #vtkSTLWriter()#vtkVRMLExporter()#vtkPolyDataWriter()#
        # writer.SetFileTypeToBinary()
        # writer.SetInputData()
        im.SetInput(renderer_window)
        im.Update()
        writer.SetInputConnection(im.GetOutputPort())
        writer.SetFileName(file_savePNG)
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
    
    if screen_on == 'false':
        print(' ')
    else:
        # # # Create the RendererWindowInteractor and display the vtk_file
        interactor = vtkRenderWindowInteractor()
        interactor.SetRenderWindow(renderer_window)
        # interactor.RemoveObservers('LeftButtonPressEvent')
        interactor.RemoveObservers('RightButtonPressEvent')
        interactor.Initialize()
        renderer_window.Render() 
        interactor.Start()
 
 
def view_results(file_dir,file_save,scale_bar,text_title,screen_on,save_screen,project_name,title_win,myfempy_version,title_data,data_display,scala_view):
       

    # OS Linux
    file_name = str(os.getcwd() + '/' +  file_dir)
    file_save = str(os.getcwd() + '/' +  file_save)
    
    # Create the Renderer
    renderer = vtkRenderer()
    # Create the RendererWindow
    renderer_window = vtkRenderWindow()
    renderer_window.SetSize(1024,768)
    renderer.SetBackground(0.0, 0.0, 0.0) # Set background to gray
    
    interactor = vtkRenderWindowInteractor()
    interactor.SetRenderWindow(renderer_window)
    
    # Read the source file.
    reader = vtkUnstructuredGridReader()
    reader.SetFileName(file_name)
    reader.Update()  # Needed because of GetScalarRange
    output = reader.GetOutput()
    output_port = reader.GetOutputPort()
    scalar_range = output.GetScalarRange()
          
    
    # Create the mapper that corresponds the objects of the vtk file
    # into graphics elements    
    lut = vtkLookupTable()
    lut.SetNumberOfColors(16)
    lut.SetHueRange(0.6, 0.0) # jet color
    lut.Build()
    
    # Axes Origin XYZ
    axes = vtkAxesActor()
    # axes.SetScale(1, 1, 1)
    axes.SetShaftTypeToCylinder()
    axes.SetXAxisLabelText('z')
    axes.SetYAxisLabelText('y')
    axes.SetZAxisLabelText('x')
    # axes.vtkTextProperty().SetFontSize(10)
        
    transform =vtkTransform()
    transform.Scale(2.6*scala_view,2.6*scala_view,2.6*scala_view)
      
    tprop = axes.GetXAxisCaptionActor2D().GetCaptionTextProperty()
    tprop.ItalicOn()
    tprop.ShadowOn()
    tprop.SetFontFamilyToArial()
    tprop.SetColor(1,1,1)
    # # Use the same text properties on the other two axes.
    axes.GetYAxisCaptionActor2D().GetCaptionTextProperty().ShallowCopy(tprop)
    axes.GetZAxisCaptionActor2D().GetCaptionTextProperty().ShallowCopy(tprop)
    axes.SetUserTransform(transform)

    mapper = vtkDataSetMapper()
    # mapper.vtkPolyDataMapper()
    mapper.SetLookupTable(lut)
    mapper.SetInputConnection(output_port)
    mapper.SetScalarRange(scalar_range)
    # mapper.SetTitleTextProperty(propT)
                
    # Create the Actor
    actor = vtkActor()
    actor.SetMapper(mapper)
    if (data_display[-1][0] == 'spring20')or(data_display[-1][0] == 'truss22')or(data_display[-1][0] == 'beam21')or(data_display[-1][0] == 'frame22')or(data_display[-1][0] == 'frame23'):
        actor.GetProperty().SetLineWidth(12.0)
    else:
        actor.GetProperty().SetLineWidth(0.5)
    actor.GetProperty().EdgeVisibilityOn()
    
    if scale_bar == 'true':
        # # Add a scalar bar
        scalarBar = vtkScalarBarActor()
        # scalarBar.GetPositionCoordinate().SetCoordinateSystemToNormalizedViewport()
        scalarBar.GetPositionCoordinate().SetValue(0.01,0.01)
        scalarBar.SetOrientationToHorizontal()
        scalarBar.SetLookupTable(lut)
        propT = vtkTextProperty()
        propT.SetFontFamilyToArial()
        propT.SetFontSize(12)
        propT.ItalicOff()
        propT.BoldOff()
        propT.ShadowOn()
        propT.SetColor(1,1,1)
        scalarBar.SetTitleTextProperty(propT)
        scalarBar.SetLabelTextProperty(propT)
        scalarBar.SetTitle(text_title)
        scalarBar.SetWidth(1.0)
        scalarBar.SetHeight(0.12)
    
    
    # create a text actor
    # txt1 = vtkTextActor()
    # txt1.SetInput("{mouse3 - scroll} - zoom / {mouse3 - click} - move screen \n {w} - wireframe style / {s} - surface style \n {q} - close visualization")
    # txtprop1=txt1.GetTextProperty()
    # txtprop1.SetFontFamilyToArial()
    # txtprop1.SetFontSize(14)
    # txtprop1.BoldOn()
    # txtprop1.SetColor(1,1,1)
    # txt1.SetDisplayPosition(10,642)
    
    # txt2 =  vtkTextActor()
    # txt2.SetInput('close the visualization to continue the analysis')
    # txtprop2=txt2.GetTextProperty()
    # txtprop2.SetFontFamilyToArial()
    # txtprop2.SetFontSize(14)
    # txtprop2.BoldOn()
    # txtprop2.SetColor(1,1,1)
    # txt2.SetDisplayPosition(10,622)
    
    text_project =  vtkTextActor()
    text_project.SetInput('JOB NAME: '+project_name+" | DATA: {"+str(date.today())+"} ")
    txtprop=text_project.GetTextProperty()
    txtprop.SetFontFamilyToArial()
    txtprop.SetFontSize(20)
    txtprop.BoldOn()
    txtprop.ItalicOn()
    txtprop.SetColor(1,1,1)
    text_project.SetDisplayPosition(10,668)
    
    text_logo =  vtkTextActor()
    text_logo.SetInput(myfempy_version)
    txtprop=text_logo.GetTextProperty()
    txtprop.SetFontFamilyToArial()
    txtprop.SetFontSize(40)
    txtprop.BoldOn()
    txtprop.ShadowOn()
    txtprop.SetColor(1,1,1)
    text_logo.SetDisplayPosition(10,718)
    
    text_title2 =  vtkTextActor()
    text_title2.SetInput(title_win)
    txtprop=text_title2.GetTextProperty()
    txtprop.SetFontFamilyToArial()
    txtprop.SetFontSize(20)
    txtprop.BoldOn()
    # txtprop.ItalicOn()
    txtprop.SetColor(1,1,1)
    text_title2.SetDisplayPosition(10,648)
    
    text_disple =  vtkTextActor()
    text_disple.SetInput(title_data[0]+str(data_display[0,0]))
    txtprop=text_disple.GetTextProperty()
    txtprop.SetFontFamilyToArial()
    txtprop.SetFontSize(18)
    txtprop.BoldOn()
    # txtprop.ItalicOn()
    txtprop.SetColor(1,1,1)
    text_disple.SetDisplayPosition(10,628)
    
    text_stremx =  vtkTextActor()
    text_stremx.SetInput(title_data[1]+str(data_display[1,0]))
    txtprop=text_stremx.GetTextProperty()
    txtprop.SetFontFamilyToArial()
    txtprop.SetFontSize(18)
    txtprop.BoldOn()
    # txtprop.ItalicOn()
    txtprop.SetColor(1,1,1)
    text_stremx.SetDisplayPosition(10,608)
    
    text_stremd =  vtkTextActor()
    text_stremd.SetInput(title_data[2]+str(data_display[2,0]))
    txtprop=text_stremd.GetTextProperty()
    txtprop.SetFontFamilyToArial()
    txtprop.SetFontSize(18)
    txtprop.BoldOn()
    # txtprop.ItalicOn()
    txtprop.SetColor(1,1,1)
    text_stremd.SetDisplayPosition(10,588)
    
    text_stremn =  vtkTextActor()
    text_stremn.SetInput(title_data[3]+str(data_display[3,0]))
    txtprop=text_stremn.GetTextProperty()
    txtprop.SetFontFamilyToArial()
    txtprop.SetFontSize(18)
    txtprop.BoldOn()
    # txtprop.ItalicOn()
    txtprop.SetColor(1,1,1)
    text_stremn.SetDisplayPosition(10,568)
    
    # cubeAxesActor = vtkCubeAxesActor()
    # # cubeAxesActor.SetFlyModeToStaticTriad()
    # cubeAxesActor.SetUseTextActor3D(1)
    # cubeAxesActor.SetBounds(actor.GetBounds())
    # cubeAxesActor.SetCamera(renderer.GetActiveCamera())
    # # cubeAxesActor.GetTitleTextProperty(0).SetColor(axis1Color)
    # cubeAxesActor.GetTitleTextProperty(0).SetFontSize(48)
    # # cubeAxesActor.GetLabelTextProperty(0).SetColor(axis1Color)
    
    # cubeAxesActor.DrawXGridlinesOn()
    # cubeAxesActor.DrawYGridlinesOn()
    # cubeAxesActor.DrawZGridlinesOn()
    # # cubeAxesActor.SetGridLineLocation(cubeAxesActor.VTK_GRID_LINES_FURTHEST)        
    # # cubeAxesActor.SetFlyModeToStaticEdges()
    # # cubeAxesActor.SetFlyModeToStaticTriad()     
    # cubeAxesActor.SetFlyModeToOuterEdges()
    
    renderer.AddActor(actor)
    
    if scale_bar == 'true':
        renderer.AddActor(scalarBar)
        
    # renderer.AddActor(txt1)
    # renderer.AddActor(txt2)
    if scala_view!=0:
        renderer.AddActor(axes)
    
    renderer.AddActor(text_project)
    renderer.AddActor(text_logo)
    renderer.AddActor(text_title2)
    renderer.AddActor(text_disple)
    renderer.AddActor(text_stremx)
    renderer.AddActor(text_stremd)
    renderer.AddActor(text_stremn)
    # renderer.AddActor(cubeAxesActor)
    
    renderer_window.AddRenderer(renderer)
    renderer_window.Render()
    renderer.ResetCamera()
    # renderer.GetActiveCamera().Elevation(45.0)
    # renderer.GetActiveCamera().Azimuth(45.0)
    
       
    if save_screen == 'false':
        print(' ')
    else: 
        # screenshot code:
        im = vtkWindowToImageFilter()
        writer = vtkPNGWriter()
        im.SetInput(renderer_window)
        im.Update()
        writer.SetInputConnection(im.GetOutputPort())
        writer.SetFileName(file_save)
        writer.Write()
    
    if screen_on == 'false':
        print(' ')
    else:
        # # # Create the RendererWindowInteractor and display the vtk_file
        interactor = vtkRenderWindowInteractor()
        interactor.SetRenderWindow(renderer_window)
        # interactor.RemoveObservers('LeftButtonPressEvent')
        interactor.RemoveObservers('RightButtonPressEvent')
        interactor.Initialize()
        renderer_window.Render() 
        interactor.Start()
    
    
def view_listforce(coord,frcApy_vet,scala_view,num_lf):   
    coordX_force = coord[int(frcApy_vet[num_lf][0])-1,1]
    coordY_force = coord[int(frcApy_vet[num_lf][0])-1,2]
    coordZ_force = coord[int(frcApy_vet[num_lf][0])-1,3]
    
    fr_text = vtkVectorText()
    
    if frcApy_vet[num_lf][2] == 0.0:
        height_cone = 0.0
    else:
        height_cone = 2.3*scala_view*abs((frcApy_vet[num_lf][2])/(max(abs(frcApy_vet[:,2])))) #abs((scala_view*frcApy_vet[num_lf][2])/(0.5*max(abs(frcApy_vet[:,2]))))
    
    if frcApy_vet[num_lf][1] == 0:
        dir_cone = (0,0,0)
        center_cone1 = (0,0,0)
        center_cone2 = center_cone1
        color_fr = (0,0,0)
        height_cone = 0
        fr_text.SetText(' ')
        fr_text.Update()
    
    if frcApy_vet[num_lf][1] == 1:   #fx 
        dir_cone = (np.sign(frcApy_vet[num_lf][2]),0,0)
        center_cone1 = (coordX_force+height_cone/2,coordY_force,coordZ_force)
        center_cone2 = center_cone1
        color_fr = (1,0,0)

        
    elif frcApy_vet[num_lf][1] == 2:  #fy 
        dir_cone = (0,np.sign(frcApy_vet[num_lf][2]),0)
        center_cone1 = (coordX_force,coordY_force+height_cone/2,coordZ_force)
        center_cone2 = center_cone1
        color_fr = (1,0,0)

    elif frcApy_vet[num_lf][1] == 3:  #fz
        dir_cone = (0,0,np.sign(frcApy_vet[num_lf][2]))
        center_cone1 = (coordX_force,coordY_force,coordZ_force+height_cone/2)
        center_cone2 = center_cone1
        color_fr = (1,0,0)
        
    elif frcApy_vet[num_lf][1] == 4:  #tx
        height_cone = 1.2*height_cone
        dir_cone = (np.sign(frcApy_vet[num_lf][2]),0,0)
        center_cone1 = (coordX_force+height_cone/2,coordY_force,coordZ_force)
        center_cone2 = (coordX_force+height_cone,coordY_force,coordZ_force)
        color_fr = (0,1,0)
        
    elif frcApy_vet[num_lf][1] == 5:  #ty
        height_cone = 1.2*height_cone
        dir_cone = (0,np.sign(frcApy_vet[num_lf][2]),0)
        center_cone1 = (coordX_force,coordY_force+height_cone/2,coordZ_force)
        center_cone2 = (coordX_force,coordY_force+height_cone,coordZ_force)
        color_fr = (0,1,0)
     
    elif frcApy_vet[num_lf][1] == 6:  #tz
        height_cone = 1.2*height_cone
        dir_cone = (0,0,np.sign(frcApy_vet[num_lf][2]))
        center_cone1 = (coordX_force,coordY_force,coordZ_force+height_cone/2)
        center_cone2 = (coordX_force,coordY_force,coordZ_force+height_cone)
        color_fr = (0,1,0)
        
    elif frcApy_vet[num_lf][1] == 7:  # massadd
        height_cone = 3*scala_view
        dir_cone = (0,-1,0)
        center_cone1 = (coordX_force,coordY_force+height_cone/2,coordZ_force)
        center_cone2 = (coordX_force,coordY_force+height_cone/2,coordZ_force)
        color_fr = (0.4,0.1,0.8)
        
    elif frcApy_vet[num_lf][1] == 8:  # spring2gd
        height_cone = 1.8*scala_view
        dir_cone = (0,1,0)
        center_cone1 = (coordX_force,coordY_force-height_cone/2,coordZ_force)
        center_cone2 = (coordX_force,coordY_force-height_cone/2,coordZ_force)
        color_fr = (0,0.9,0.8)
        
    elif frcApy_vet[num_lf][1] == 9:  # damper2gd
        height_cone = 1.8*scala_view#*(0.5*max(abs(frcApy_vet[:,2])))/(0.5*max(abs(frcApy_vet[:,2])))
        dir_cone = (0,1,0)
        center_cone1 = (coordX_force,coordY_force,coordZ_force)
        center_cone2 = (coordX_force,coordY_force,coordZ_force)
        color_fr = (0.9,0.6,0.2)
    
    
    cone1 = vtkConeSource()
    cone1.SetResolution(1)
    cone1.SetHeight(height_cone)
    cone1.SetRadius(0.15*height_cone)
    cone1.SetCenter(center_cone1)
    cone1.SetDirection(dir_cone)
        
    forcemap1 = vtkPolyDataMapper()
    forcemap1.SetInputConnection(cone1.GetOutputPort() )
    fr_point_actor_cone1 = vtkActor()
    fr_point_actor_cone1.SetMapper(forcemap1)
    fr_point_actor_cone1.GetProperty().SetLineWidth(0.5)
    fr_point_actor_cone1.GetProperty().EdgeVisibilityOn()
    # color rotated cone blue
    fr_point_actor_cone1.GetProperty().SetColor(color_fr) # (R,G,B)  
    
    cone2 = vtkConeSource()
    cone2.SetResolution(1)
    cone2.SetHeight(height_cone)
    cone2.SetRadius(0.15*height_cone)
    cone2.SetCenter(center_cone2)
    cone2.SetDirection(dir_cone) 
        
    forcemap2 = vtkPolyDataMapper()
    forcemap2.SetInputConnection(cone2.GetOutputPort() )
    fr_point_actor_cone2 = vtkActor()
    fr_point_actor_cone2.SetMapper(forcemap2)
    fr_point_actor_cone2.GetProperty().SetLineWidth(0.5)
    fr_point_actor_cone2.GetProperty().EdgeVisibilityOn()
    # color rotated cone blue
    fr_point_actor_cone2.GetProperty().SetColor(color_fr) # (R,G,B)  
    
    fr_text.SetText(str(np.round(abs(frcApy_vet[num_lf][2]),2)))
    fr_text.Update()
    fr_text_map = vtkPolyDataMapper()
    fr_text_map.SetInputConnection(fr_text.GetOutputPort())
    fr_text_actor = vtkActor()
    fr_text_actor.SetMapper(fr_text_map)
    fr_text_actor.SetScale(3,3,3)
    fr_text_actor.SetPosition(center_cone1)
    # bc_text_actor.GetProperty().SetLineWidth(0.5)
    # bc_text_actor.GetProperty().EdgeVisibilityOn()
    fr_text_actor.GetProperty().SetColor(0,1,0)
    
    return fr_point_actor_cone1, fr_point_actor_cone2, fr_text_actor


def view_bondcond_point(coord,bondCond_vet,scala_view,num_bc):    
    coordX_bc = coord[int(bondCond_vet[num_bc,1])-1,1]
    coordY_bc = coord[int(bondCond_vet[num_bc,1])-1,2]
    coordZ_bc = coord[int(bondCond_vet[num_bc,1])-1,3]
    
    height_cone = 0.9*scala_view
    color_rgb = (1,1,0)
    max_coordX = max(coord[:,1])
    
    bc_text = vtkVectorText()
    
    if int(bondCond_vet[num_bc,0]) == 0: #fixed all dofs        
        if coordX_bc == max_coordX:
            dir_cone = (-1,0,0)
            center_cone = (coordX_bc+height_cone/2,coordY_bc,coordZ_bc)
            center_cube = (coordX_bc+2.25*height_cone/2,coordY_bc,coordZ_bc)
        else:
            dir_cone = (1,0,0)
            center_cone = (coordX_bc-height_cone/2,coordY_bc,coordZ_bc)
            center_cube = (coordX_bc-2.25*height_cone/2,coordY_bc,coordZ_bc)
        
        bc_text.SetText('F')
        bc_text.Update()
      
        
        cube = vtkCubeSource()
        cube.SetXLength(2*height_cone)
        cube.SetYLength(2*height_cone)
        cube.SetZLength(0*height_cone)
        cube.SetCenter(center_cube)
        
        bcmap = vtkPolyDataMapper()
        bcmap.SetInputConnection(cube.GetOutputPort())
        bc_point_actor_tdof = vtkActor()
        bc_point_actor_tdof.SetMapper(bcmap)
        bc_point_actor_tdof.GetProperty().SetLineWidth(0.5)
        bc_point_actor_tdof.GetProperty().EdgeVisibilityOn()
        # color rotated cone blue
        bc_point_actor_tdof.GetProperty().SetColor(1,1,0) # (R,G,B)
        # height_cone = 0.0
        
    elif int(bondCond_vet[num_bc,0]) == 1: #fixed ux dofs   
        dir_cone = (1,0,0)
        
        if coordX_bc == max_coordX:
            dir_cone = (-1,0,0)
            center_cone = (coordX_bc+height_cone/2,coordY_bc,coordZ_bc)
            # center_sphere = (coordX_bc+2.8*height_cone/2,coordY_bc,coordZ_bc)
            center_cube = (coordX_bc+3.5*height_cone/2,coordY_bc,coordZ_bc)
        else:
            dir_cone = (1,0,0)
            center_cone = (coordX_bc-height_cone/2,coordY_bc,coordZ_bc)
            # center_sphere = (coordX_bc-2.8*height_cone/2,coordY_bc,coordZ_bc)
            center_cube = (coordX_bc-3.5*height_cone/2,coordY_bc,coordZ_bc)
        
        # sphere = vtkSphereSource()
        # sphere.SetThetaResolution(20)
        # sphere.SetRadius(0.35*scala_view)
        # sphere.SetCenter(center_sphere)
        
        # # cylinder = vtkCylinderSource()
        # # cylinder.SetResolution(20)
        # # cylinder.SetHeight(0.4*height_cone)
        # # cylinder.SetRadius(0.4*height_cone)
        # # cylinder.SetCenter(coordX_bc,coordY_bc,0)
        
        
        # bcmap = vtkPolyDataMapper()
        # bcmap.SetInputConnection(sphere.GetOutputPort())
        # bc_point_actor_tdof = vtkActor()
        # bc_point_actor_tdof.SetMapper(bcmap)
        # # color rotated cone blue
        # bc_point_actor_tdof.GetProperty().SetColor(1,1,0) # (R,G,B)
        
        bc_text.SetText('UX')
        bc_text.Update()
        
        cube = vtkCubeSource()
        cube.SetXLength(0.5*height_cone)
        cube.SetYLength(2.5*height_cone)
        cube.SetZLength(0*height_cone)
        cube.SetCenter(center_cube)
        
        bcmap = vtkPolyDataMapper()
        bcmap.SetInputConnection(cube.GetOutputPort())
        bc_point_actor_tdof = vtkActor()
        bc_point_actor_tdof.SetMapper(bcmap)
        bc_point_actor_tdof.GetProperty().SetLineWidth(0.5)
        bc_point_actor_tdof.GetProperty().EdgeVisibilityOn()
        # color rotated cone blue
        bc_point_actor_tdof.GetProperty().SetColor(1,1,0) # (R,G,B)

        # translante = vtkTransform()
        # translante.Translate(coordX_bc,coordY_bc,0)
        # bc_point_actor_tdof
        
        # bc_point_actor_tdof.RotateX(90.0)
        # bc_point_actor_tdof.Translate(coordX_bc,coordY_bc,0)
        # bc_point_actor_tdof.RotateWXYZ(90,1,0,0)
        # bc_point_actor.RotateY(-45.0)
        
    elif int(bondCond_vet[num_bc,0]) == 2: #fixed uy dofs   
        dir_cone = (0,1,0)
        
        center_cone = (coordX_bc,coordY_bc-height_cone/2,coordZ_bc)
        # center_sphere = (coordX_bc+2.8*height_cone/2,coordY_bc,coordZ_bc)
        center_cube = (coordX_bc,coordY_bc-3.5*height_cone/2,coordZ_bc)
        
        # sphere = vtkSphereSource()
        # sphere.SetThetaResolution(20)
        # sphere.SetRadius(0.35*scala_view)
        # sphere.SetCenter(center_sphere)
        
        # # cylinder = vtkCylinderSource()
        # # cylinder.SetResolution(20)
        # # cylinder.SetHeight(0.4*height_cone)
        # # cylinder.SetRadius(0.4*height_cone)
        # # cylinder.SetCenter(coordX_bc,coordY_bc,0)
        
        
        # bcmap = vtkPolyDataMapper()
        # bcmap.SetInputConnection(sphere.GetOutputPort())
        # bc_point_actor_tdof = vtkActor()
        # bc_point_actor_tdof.SetMapper(bcmap)
        # # color rotated cone blue
        # bc_point_actor_tdof.GetProperty().SetColor(1,1,0) # (R,G,B)
        
        bc_text.SetText('UY')
        bc_text.Update()

        
        cube = vtkCubeSource()
        cube.SetXLength(2.5*height_cone)
        cube.SetYLength(0.5*height_cone)
        cube.SetZLength(0*height_cone)
        cube.SetCenter(center_cube)
        
        bcmap = vtkPolyDataMapper()
        bcmap.SetInputConnection(cube.GetOutputPort())
        bc_point_actor_tdof = vtkActor()
        bc_point_actor_tdof.SetMapper(bcmap)
        bc_point_actor_tdof.GetProperty().SetLineWidth(0.5)
        bc_point_actor_tdof.GetProperty().EdgeVisibilityOn()
        # color rotated cone blue
        bc_point_actor_tdof.GetProperty().SetColor(1,1,0) # (R,G,B)

    elif int(bondCond_vet[num_bc,0]) == 3: #fixed uz dofs   
        dir_cone = (0,0,1)

        center_cone = (coordX_bc,coordY_bc,coordZ_bc-height_cone/2)
        # center_sphere = (coordX_bc+2.8*height_cone/2,coordY_bc,coordZ_bc)
        center_cube = (coordX_bc,coordY_bc,coordZ_bc-3.5*height_cone/2)

        
        # sphere = vtkSphereSource()
        # sphere.SetThetaResolution(20)
        # sphere.SetRadius(0.35*scala_view)
        # sphere.SetCenter(center_sphere)
        
        # # cylinder = vtkCylinderSource()
        # # cylinder.SetResolution(20)
        # # cylinder.SetHeight(0.4*height_cone)
        # # cylinder.SetRadius(0.4*height_cone)
        # # cylinder.SetCenter(coordX_bc,coordY_bc,0)
        
        
        # bcmap = vtkPolyDataMapper()
        # bcmap.SetInputConnection(sphere.GetOutputPort())
        # bc_point_actor_tdof = vtkActor()
        # bc_point_actor_tdof.SetMapper(bcmap)
        # # color rotated cone blue
        # bc_point_actor_tdof.GetProperty().SetColor(1,1,0) # (R,G,B)
        
        bc_text.SetText('UZ')
        bc_text.Update()

        
        cube = vtkCubeSource()
        cube.SetXLength(0*height_cone)
        cube.SetYLength(2.5*height_cone)
        cube.SetZLength(0.5*height_cone)
        cube.SetCenter(center_cube)
        
        bcmap = vtkPolyDataMapper()
        bcmap.SetInputConnection(cube.GetOutputPort())
        bc_point_actor_tdof = vtkActor()
        bc_point_actor_tdof.SetMapper(bcmap)
        bc_point_actor_tdof.GetProperty().SetLineWidth(0.5)
        bc_point_actor_tdof.GetProperty().EdgeVisibilityOn()
        # color rotated cone blue
        bc_point_actor_tdof.GetProperty().SetColor(1,1,0) # (R,G,B)
        
    elif int(bondCond_vet[num_bc,0]) == 4: #fixed rx dofs   
        
        dir_cone = (0,0,1)

        center_cone = (coordX_bc,coordY_bc,coordZ_bc-height_cone/2)
        # center_sphere = (coordX_bc+2.8*height_cone/2,coordY_bc,coordZ_bc)
        center_cube = (coordX_bc,coordY_bc,coordZ_bc-3.5*height_cone/2)

        
        # sphere = vtkSphereSource()
        # sphere.SetThetaResolution(20)
        # sphere.SetRadius(0.35*scala_view)
        # sphere.SetCenter(center_sphere)
        
        # # cylinder = vtkCylinderSource()
        # # cylinder.SetResolution(20)
        # # cylinder.SetHeight(0.4*height_cone)
        # # cylinder.SetRadius(0.4*height_cone)
        # # cylinder.SetCenter(coordX_bc,coordY_bc,0)
        
        
        # bcmap = vtkPolyDataMapper()
        # bcmap.SetInputConnection(sphere.GetOutputPort())
        # bc_point_actor_tdof = vtkActor()
        # bc_point_actor_tdof.SetMapper(bcmap)
        # # color rotated cone blue
        # bc_point_actor_tdof.GetProperty().SetColor(1,1,0) # (R,G,B)
        
        bc_text.SetText('RX')
        bc_text.Update()
        
        cube = vtkCubeSource()
        cube.SetXLength(0*height_cone)
        cube.SetYLength(2.5*height_cone)
        cube.SetZLength(0.5*height_cone)
        cube.SetCenter(center_cube)
        
        bcmap = vtkPolyDataMapper()
        bcmap.SetInputConnection(cube.GetOutputPort())
        bc_point_actor_tdof = vtkActor()
        bc_point_actor_tdof.SetMapper(bcmap)
        bc_point_actor_tdof.GetProperty().SetLineWidth(0.5)
        bc_point_actor_tdof.GetProperty().EdgeVisibilityOn()
        # color rotated cone blue
        bc_point_actor_tdof.GetProperty().SetColor(1,1,0) # (R,G,B)

    elif int(bondCond_vet[num_bc,0]) == 5: #fixed ry dofs   
        
        dir_cone = (1,0,0)
        
        if coordX_bc == max_coordX:
            dir_cone = (-1,0,0)
            center_cone = (coordX_bc+height_cone/2,coordY_bc,coordZ_bc)
            # center_sphere = (coordX_bc+2.8*height_cone/2,coordY_bc,coordZ_bc)
            center_cube = (coordX_bc+3.5*height_cone/2,coordY_bc,coordZ_bc)
        else:
            dir_cone = (1,0,0)
            center_cone = (coordX_bc-height_cone/2,coordY_bc,coordZ_bc)
            # center_sphere = (coordX_bc-2.8*height_cone/2,coordY_bc,coordZ_bc)
            center_cube = (coordX_bc-3.5*height_cone/2,coordY_bc,coordZ_bc)
        
        # sphere = vtkSphereSource()
        # sphere.SetThetaResolution(20)
        # sphere.SetRadius(0.35*scala_view)
        # sphere.SetCenter(center_sphere)
        
        # # cylinder = vtkCylinderSource()
        # # cylinder.SetResolution(20)
        # # cylinder.SetHeight(0.4*height_cone)
        # # cylinder.SetRadius(0.4*height_cone)
        # # cylinder.SetCenter(coordX_bc,coordY_bc,0)
        
        
        # bcmap = vtkPolyDataMapper()
        # bcmap.SetInputConnection(sphere.GetOutputPort())
        # bc_point_actor_tdof = vtkActor()
        # bc_point_actor_tdof.SetMapper(bcmap)
        # # color rotated cone blue
        # bc_point_actor_tdof.GetProperty().SetColor(1,1,0) # (R,G,B)
        
        bc_text.SetText('RY')
        bc_text.Update()
        
        cube = vtkCubeSource()
        cube.SetXLength(0.5*height_cone)
        cube.SetYLength(2.5*height_cone)
        cube.SetZLength(0*height_cone)
        cube.SetCenter(center_cube)
        
        bcmap = vtkPolyDataMapper()
        bcmap.SetInputConnection(cube.GetOutputPort())
        bc_point_actor_tdof = vtkActor()
        bc_point_actor_tdof.SetMapper(bcmap)
        bc_point_actor_tdof.GetProperty().SetLineWidth(0.5)
        bc_point_actor_tdof.GetProperty().EdgeVisibilityOn()
        # color rotated cone blue
        bc_point_actor_tdof.GetProperty().SetColor(1,1,0) # (R,G,B)
        
        
        
        # translante = vtkTransform()
        # translante.Translate(coordX_bc,coordY_bc,0)
        # bc_point_actor_tdof
        
        # bc_point_actor_tdof.RotateX(90.0)
        # bc_point_actor_tdof.Translate(coordX_bc,coordY_bc,0)
        # bc_point_actor_tdof.RotateWXYZ(90,1,0,0)
        # bc_point_actor.RotateY(-45.0)
     
    elif int(bondCond_vet[num_bc,0]) == 6: #fixed rz dofs   
        dir_cone = (0,1,0)
        
        center_cone = (coordX_bc,coordY_bc-height_cone/2,coordZ_bc)
        # center_sphere = (coordX_bc+2.8*height_cone/2,coordY_bc,coordZ_bc)
        center_cube = (coordX_bc,coordY_bc-3.5*height_cone/2,coordZ_bc)
        
        # sphere = vtkSphereSource()
        # sphere.SetThetaResolution(20)
        # sphere.SetRadius(0.35*scala_view)
        # sphere.SetCenter(center_sphere)
        
        # # cylinder = vtkCylinderSource()
        # # cylinder.SetResolution(20)
        # # cylinder.SetHeight(0.4*height_cone)
        # # cylinder.SetRadius(0.4*height_cone)
        # # cylinder.SetCenter(coordX_bc,coordY_bc,0)
        
        
        # bcmap = vtkPolyDataMapper()
        # bcmap.SetInputConnection(sphere.GetOutputPort())
        # bc_point_actor_tdof = vtkActor()
        # bc_point_actor_tdof.SetMapper(bcmap)
        # # color rotated cone blue
        # bc_point_actor_tdof.GetProperty().SetColor(1,1,0) # (R,G,B)
    
        bc_text.SetText('RZ')
        bc_text.Update()
 
        
        cube = vtkCubeSource()
        cube.SetXLength(2.5*height_cone)
        cube.SetYLength(0.5*height_cone)
        cube.SetZLength(0*height_cone)
        cube.SetCenter(center_cube)
        
        bcmap = vtkPolyDataMapper()
        bcmap.SetInputConnection(cube.GetOutputPort())
        bc_point_actor_tdof = vtkActor()
        bc_point_actor_tdof.SetMapper(bcmap)
        bc_point_actor_tdof.GetProperty().SetLineWidth(0.5)
        bc_point_actor_tdof.GetProperty().EdgeVisibilityOn()
        # color rotated cone blue
        bc_point_actor_tdof.GetProperty().SetColor(1,1,0) # (R,G,B)

    
    # tdof = bondCond_vet[tbc,0]
    cone = vtkConeSource()
    cone.SetResolution(1)
    cone.SetHeight(height_cone)
    cone.SetRadius(0.5*height_cone)
    cone.SetCenter(center_cone)
    cone.SetDirection(dir_cone) 
    
    bccmap = vtkPolyDataMapper()
    bccmap.SetInputConnection(cone.GetOutputPort() )
    bc_point_actor_cone = vtkActor()
    bc_point_actor_cone.SetMapper(bccmap)
    bc_point_actor_cone.GetProperty().SetLineWidth(0.5)
    bc_point_actor_cone.GetProperty().EdgeVisibilityOn()
    # color rotated cone blue
    bc_point_actor_cone.GetProperty().SetColor(color_rgb) # (R,G,B)
    
    bc_text_map = vtkPolyDataMapper()
    bc_text_map.SetInputConnection(bc_text.GetOutputPort())
    bc_text_actor = vtkActor()
    bc_text_actor.SetMapper(bc_text_map)
    bc_text_actor.SetScale(5,5,5)
    bc_text_actor.SetPosition(center_cone)
    # bc_text_actor.GetProperty().SetLineWidth(0.5)
    # bc_text_actor.GetProperty().EdgeVisibilityOn()
    bc_text_actor.GetProperty().SetColor(0,0,1)
       
    return bc_point_actor_cone, bc_point_actor_tdof,bc_text_actor


# def view_bondcond_edge(coord,bondCond_vet,scala_view,tdof,node_list_eg):

#     coordX_bc = coord[node_list_eg-1,1]
#     coordY_bc = coord[node_list_eg-1,2]
    
#     minX_point_bc = min(coordX_bc)
#     minY_point_bc = min(coordY_bc)
#     maxX_point_bc = max(coordX_bc)
#     maxY_point_bc = max(coordY_bc)
        
#     Cx_bc = (minX_point_bc+maxX_point_bc)/2
#     Cy_bc = (minY_point_bc+maxY_point_bc)/2      
    
#     height_cone = 0.6*scala_view
            
#     if tdof == 0:
#         cube = vtkCubeSource()
#         cube.SetXLength(0.3*scala_view)
#         cube.SetYLength(0.3*scala_view)
#         cube.SetCenter(Cx_bc,Cy_bc,0)
        
#         bcmap = vtkPolyDataMapper()
#         bcmap.SetInputConnection(cube.GetOutputPort() )
#         bc_center_edge_actor = vtkActor()
#         bc_center_edge_actor.SetMapper(bcmap)
#         # color rotated cone blue
#         bc_center_edge_actor.GetProperty().SetColor(1,1,0) # (R,G,B)
        
#         line = vtkLineSource()
#         line.SetPoint1(minX_point_bc,minY_point_bc,0)
#         line.SetPoint2(maxX_point_bc,maxY_point_bc,0)
        
#         bcedgemap = vtkPolyDataMapper()
#         bcedgemap.SetInputConnection(line.GetOutputPort() )
#         bc_edge_actor = vtkActor()
#         bc_edge_actor.SetMapper(bcedgemap)
#         # color rotated cone blue
#         bc_edge_actor.GetProperty().SetColor(1,1,0) # (R,G,B)
#         bc_edge_actor.GetProperty().SetLineWidth(8.0)
        
        
#     elif tdof > 0:
#         # cylinder = vtkCylinderSource()
#         # cylinder.SetHeight(0.4*height_cone)
#         # cylinder.SetRadius(0.4*height_cone)
#         # cylinder.SetCenter(Cx_bc,Cy_bc,0)
#         # cylinder.SetDirection(1,0,0)
        
#         # bcmap = vtkPolyDataMapper()
#         # # bcmap.SetInputConnection(cylinder.GetOutputPort() )
        
#         sphere = vtkSphereSource()
#         sphere.SetThetaResolution(8)
#         sphere.SetRadius(0.2*scala_view)
#         sphere.SetCenter(Cx_bc,Cy_bc,0)
        
#         bcmap = vtkPolyDataMapper()
#         bcmap.SetInputConnection(sphere.GetOutputPort())
        
#         bc_center_edge_actor = vtkActor()
#         bc_center_edge_actor.SetMapper(bcmap)
#         # color rotated cone blue
#         bc_center_edge_actor.GetProperty().SetColor(0,1,0) # (R,G,B)
#         # bc_center_edge_actor.RotateX(90.0)
#         # bc_center_edge_actor.RotateY(0.0)
#         # bc_center_edge_actor.RotateZ(0.0)
        
#         line = vtkLineSource()
#         line.SetPoint1(minX_point_bc,minY_point_bc,0)
#         line.SetPoint2(maxX_point_bc,maxY_point_bc,0)
        
#         bcedgemap = vtkPolyDataMapper()
#         bcedgemap.SetInputConnection(line.GetOutputPort())
#         bc_edge_actor = vtkActor()
#         bc_edge_actor.SetMapper(bcedgemap)
#         # color rotated cone blue
#         bc_edge_actor.GetProperty().SetColor(0,1,0) # (R,G,B)
#         bc_edge_actor.GetProperty().SetLineWidth(8.0)
    
    
#     return bc_edge_actor, bc_center_edge_actor


def view_beam_crossSection(dimSection,typSection,coord_bcs):
    
    b = dimSection[0]
    h = dimSection[1]
    t = dimSection[2]
    d = dimSection[3]
    
    Lx = np.sqrt(((coord_bcs[3]-coord_bcs[0])**2))
    Ly = np.sqrt(((coord_bcs[4]-coord_bcs[1])**2))
    Lz = np.sqrt(((coord_bcs[5]-coord_bcs[2])**2))
    L = np.sqrt(Lx**2 + Ly**2 + Lz**2)
    
    if (Lx>0.0)and(Ly==0.0)and(Lz==0.0):
        translate = (coord_bcs[0], coord_bcs[1], coord_bcs[2])
        rotate1 = (0,1,0)
        rotate2 = (0,0,0)
        rotate3 = (0,0,0)
        ang1 = 90
        ang2 = 0
        ang3 = 0
    
        
    if (Lx==0.0)and(Ly>0.0)and(Lz==0.0):
        translate = (coord_bcs[0], coord_bcs[1], coord_bcs[2])
        rotate1 = (1,0,0)
        rotate2 = (0,0,0)
        rotate3 = (0,0,1)
        ang1 = 90
        ang2 = 0
        ang3 = 90

            
    if (Lx==0.0)and(Ly==0.0)and(Lz>0.0):
        translate = (coord_bcs[0], coord_bcs[1], coord_bcs[2])
        rotate1 = (0,0,0)
        rotate2 = (0,0,0)
        rotate3 = (0,1,0)
        ang1 = 0
        ang2 = 0
        ang3 = -90
    
    if (Lx>0.0)and(Ly>0.0)and(Lz==0.0):
        translate = (coord_bcs[0], coord_bcs[1], coord_bcs[2])
        
        
        l = (coord_bcs[3]-coord_bcs[0])/L
        m = (coord_bcs[1]-coord_bcs[4])/L
        n = (coord_bcs[5]-coord_bcs[2])/L
        d = np.sqrt(l**2 + m**2)
        
        ang1 = 90
        rotate1 = (0,m,0)
        if l<0:
            ang2 = (np.arctan(l/m))*(180/np.pi)
        else:
            ang2 = (np.arctan(-l/m))*(180/np.pi)
        
        rotate2 = (l,0,0)
        ang3 = 0
        rotate3 = (0,0,0)
        
    if (Lx>0.0)and(Ly==0.0)and(Lz>0.0):
        translate = (coord_bcs[0], coord_bcs[1], coord_bcs[2])
        
        l = (coord_bcs[3]-coord_bcs[0])/L
        m = (coord_bcs[1]-coord_bcs[4])/L
        n = (coord_bcs[5]-coord_bcs[2])/L
        d = np.sqrt(l**2 + m**2)
        
        ang1 = 90
        rotate1 = (0,m,0)
        ang2 = (np.arctan(l/n))*(180/np.pi)
        rotate2 = (0,m,0)
        ang3 = 0
        rotate3 = (0,0,0)
        
    if (Lx==0.0)and(Ly>0.0)and(Lz>0.0):
        translate = (coord_bcs[0], coord_bcs[1], coord_bcs[2])
        
        l = (coord_bcs[3]-coord_bcs[0])/L
        m = (coord_bcs[1]-coord_bcs[4])/L
        n = (coord_bcs[5]-coord_bcs[2])/L
        d = np.sqrt(l**2 + m**2)
        
        ang1 = 0
        rotate1 = (0,0,0)
        ang2 = (np.arctan(m/n))*(180/np.pi)
        rotate2 = (d,0,0)
        ang3 = 0
        rotate3 = (0,0,0)
        
    if (Lx>0.0)and(Ly>0.0)and(Lz>0.0):
        translate = (coord_bcs[0], coord_bcs[1], coord_bcs[2])
        
        l = (coord_bcs[3]-coord_bcs[0])/L
        m = (coord_bcs[1]-coord_bcs[4])/L
        n = (coord_bcs[5]-coord_bcs[2])/L
        d = np.sqrt(l**2 + m**2)
        
        ang1 = (np.arctan(l/m))*(180/np.pi)
        rotate1 = (0,m,0)
        ang2 = (np.arctan(l/m))*(180/np.pi)
        rotate2 = (l,0,0)
        ang3 = 0
        rotate3 = (0,0,0)
        
    if typSection == 'rectangle':
        
        # planeSource = vtkPlaneSource()
        # planeSource.SetOrigin(1000, 0, 0.0)
        # planeSource.SetNormal(1.0, 0.0, 0.0)
        # planeSource.SetPoint2(-10.0, 10.0, 0.0)
        # planeSource.SetPoint1(10.0, -10.0, 0.0)
        # planeSource.Update()        
        
        points = vtkPoints()
        points.SetNumberOfPoints(4)
        points.SetPoint(0, -b/2, -h/2, 0.0)
        points.SetPoint(1, b/2, -h/2, 0.0)
        points.SetPoint(2, b/2, h/2, 0.0)
        points.SetPoint(3, -b/2, h/2, 0.0)
        
        lines = vtkCellArray()
        lines.InsertNextCell(4)
        lines.InsertCellPoint(0)
        lines.InsertCellPoint(1)
        lines.InsertCellPoint(2)
        lines.InsertCellPoint(3)
        
        profile = vtkPolyData()
        profile.SetPoints(points)
        profile.SetPolys(lines)
        
        transform = vtkTransform()
        transform.Translate(translate)
        transform.RotateWXYZ(ang1,rotate1)
        transform.RotateWXYZ(ang2,rotate2)
        # transform.RotateWXYZ(45,(1,0,0))
        transform_filter = vtkTransformPolyDataFilter()
        transform_filter.SetTransform(transform)
        transform_filter.SetInputData(profile)
        transform_filter.Update()
        
        # extrude = vtkLinearExtrusionFilter()
        # extrude.SetInputConnection(transform_filter.GetOutputPort())
        # extrude.SetExtrusionTypeToNormalExtrusion()
        # extrude.SetVector(Lx, Ly,Lz)
        # # extrude.SetScaleFactor(1)
        # extrude.Update()
        
    if typSection == 'rectangle_tube':
        
        points = vtkPoints()
        points.SetNumberOfPoints(12)
        points.SetPoint(0, 0.0, h/2, 0.0)
        points.SetPoint(1, -b/2, h/2, 0.0)
        points.SetPoint(2, -b/2, -h/2, 0.0)
        points.SetPoint(3, b/2, -h/2, 0.0)
        points.SetPoint(4, b/2, h/2, 0.0)
        points.SetPoint(5, 0.0, h/2, 0.0)
        points.SetPoint(6, 0.0, h/2-t, 0.0)
        points.SetPoint(7, b/2-t, h/2-t, 0.0)
        points.SetPoint(8, b/2-t, -h/2+t, 0.0)
        points.SetPoint(9, -b/2+t, -h/2+t, 0.0)
        points.SetPoint(10, -b/2+t, h/2-t, 0.0)
        points.SetPoint(11, 0.0, h/2-t, 0.0)
        
        lines = vtkCellArray()
        lines.InsertNextCell(12)
        lines.InsertCellPoint(0)
        lines.InsertCellPoint(1)
        lines.InsertCellPoint(2)
        lines.InsertCellPoint(3)
        lines.InsertCellPoint(4)
        lines.InsertCellPoint(5)
        lines.InsertCellPoint(6)
        lines.InsertCellPoint(7)
        lines.InsertCellPoint(8)
        lines.InsertCellPoint(9)
        lines.InsertCellPoint(10)
        lines.InsertCellPoint(11)
        
        profile = vtkPolyData()
        profile.SetPoints(points)
        profile.SetPolys(lines)
        
        transform = vtkTransform()
        transform.Translate(translate)
        transform.RotateWXYZ(ang1,rotate1)
        transform.RotateWXYZ(ang2,rotate2)
        transform_filter = vtkTransformPolyDataFilter()
        transform_filter.SetTransform(transform)
        transform_filter.SetInputData(profile)
        transform_filter.Update()
        
        extrude = vtkLinearExtrusionFilter()
        extrude.SetInputConnection(transform_filter.GetOutputPort())
        extrude.SetExtrusionTypeToNormalExtrusion()
        extrude.SetVector(Lx, Ly,Lz)
        # extrude.SetScaleFactor(1)
        extrude.Update()

    if typSection == 'I_section':
                
        points = vtkPoints()
        points.SetNumberOfPoints(12)
        points.SetPoint(0, b/2,h/2,0)
        points.SetPoint(1, -b/2,h/2,0)
        points.SetPoint(2, -b/2,h/2-d,0)
        points.SetPoint(3, -t/2,h/2-d,0)
        points.SetPoint(4, -t/2,-h/2+d,0)
        points.SetPoint(5, -b/2,-h/2+d,0)
        points.SetPoint(6, -b/2,-h/2,0)
        points.SetPoint(7, b/2,-h/2,0)
        points.SetPoint(8, b/2,-h/2+d,0)
        points.SetPoint(9, t/2,-h/2+d,0)
        points.SetPoint(10, t/2,h/2-d,0)
        points.SetPoint(11, b/2,h/2-d,0)
        
        lines = vtkCellArray()
        lines.InsertNextCell(13)
        lines.InsertCellPoint(0)
        lines.InsertCellPoint(1)
        lines.InsertCellPoint(2)
        lines.InsertCellPoint(3)
        lines.InsertCellPoint(4)
        lines.InsertCellPoint(5)
        lines.InsertCellPoint(6)
        lines.InsertCellPoint(7)
        lines.InsertCellPoint(8)
        lines.InsertCellPoint(9)
        lines.InsertCellPoint(10)
        lines.InsertCellPoint(11)
        lines.InsertCellPoint(0)
        
        
        profile = vtkPolyData()
        profile.SetPoints(points)
        profile.SetPolys(lines)
        
        transform = vtkTransform()
        transform.Translate(translate)
        transform.RotateWXYZ(ang1,rotate1)
        transform.RotateWXYZ(ang2,rotate2)
        transform_filter = vtkTransformPolyDataFilter()
        transform_filter.SetTransform(transform)
        transform_filter.SetInputData(profile)
        transform_filter.Update()
        
        # extrude = vtkLinearExtrusionFilter()
        # extrude.SetInputConnection(transform_filter.GetOutputPort())
        # extrude.SetExtrusionTypeToNormalExtrusion()
        # extrude.SetVector(Lx, Ly,Lz)
        # # extrude.SetScaleFactor(1)
        # extrude.Update()


    if typSection == 'circle':
        
        profile = vtkDiskSource()
        profile.SetCircumferentialResolution(64)
        profile.SetRadialResolution(1)
        profile.SetOuterRadius(d/2)
        profile.SetInnerRadius(0)
        
        transform = vtkTransform()
        transform.Translate(translate)
        transform.RotateWXYZ(ang1,rotate1)
        transform.RotateWXYZ(ang2,rotate2)
        transform_filter = vtkTransformPolyDataFilter()
        transform_filter.SetTransform(transform)
        transform_filter.SetInputConnection(profile.GetOutputPort())
        transform_filter.Update()
        
        # extrude = vtkLinearExtrusionFilter()
        # extrude.SetInputConnection(transform_filter.GetOutputPort())
        # extrude.SetExtrusionTypeToNormalExtrusion()
        # extrude.SetVector(Lx, Ly,Lz)
        # # extrude.SetScaleFactor(1)
        # extrude.Update()
        
    if typSection == 'circle_tube':
        
        profile = vtkDiskSource()
        profile.SetCircumferentialResolution(64)
        profile.SetRadialResolution(1)
        profile.SetOuterRadius(d/2)
        profile.SetInnerRadius(d/2-t)
        
        transform = vtkTransform()
        transform.Translate(translate)
        transform.RotateWXYZ(ang1,rotate1)
        transform.RotateWXYZ(ang2,rotate2)
        
        transform_filter = vtkTransformPolyDataFilter()
        transform_filter.SetTransform(transform)
        transform_filter.SetInputConnection(profile.GetOutputPort())
        transform_filter.Update()
        
        # extrude = vtkLinearExtrusionFilter()
        # extrude.SetInputConnection(transform_filter.GetOutputPort())
        # extrude.SetExtrusionTypeToNormalExtrusion()
        # extrude.SetVector(Lx, Ly,Lz)
        # # extrude.SetScaleFactor(1)
        # extrude.Update()
        
    if typSection == 'spring':
        
        p0 = [0.0, 0.0, 0.0]
        p1 = [0.25*L, 0.0, 0.0]
        p2 = [0.5*L, 0.25*L, 0.0]
        p3 = [0.5*L, -0.25*L, 0.0]
        p4 = [0.75*L, 0.0, 0.0]
        p5 = [L, 0.0, 0.0]
    
        # Create a vtkPoints object and store the points in it
        points =vtkPoints()
        points.InsertNextPoint(p0)
        points.InsertNextPoint(p1)
        points.InsertNextPoint(p2)
        points.InsertNextPoint(p3)
        points.InsertNextPoint(p4)
        points.InsertNextPoint(p5)
    
        # Create a cell array to store the lines in and add the lines to it
        lines = vtkCellArray()
        lines.InsertNextCell(6)
        lines.InsertCellPoint(0)
        lines.InsertCellPoint(1)
        lines.InsertCellPoint(2)
        lines.InsertCellPoint(3)
        lines.InsertCellPoint(4)
        lines.InsertCellPoint(5)
    
        profile = vtkPolyData()
        profile.SetPoints(points)
        profile.SetLines(lines)
        
        transform = vtkTransform()
        transform.Translate(translate)
        transform.RotateWXYZ(ang3,rotate3)
        transform_filter = vtkTransformPolyDataFilter()
        transform_filter.SetTransform(transform)
        transform_filter.SetInputData(profile)
        transform_filter.Update()
        
        # sphere = vtkSphereSource()
        # sphere.SetThetaResolution(20)
        # sphere.SetRadius(0.35*L)
        # sphere.SetCenter(coord_bcs[3],coord_bcs[4],coord_bcs[5])    
        # sphere.Update()
              
        # transform_filter = vtkMultiBlockDataSet()
        # transform_filter.SetNumberOfBlocks(2)
        # transform_filter.SetBlock(0, transform_2.GetOutput())
        # transform_filter.SetBlock(1, sphere.GetOutput())
        
        # mapper = vtkCompositePolyDataMapper2()
        # mapper.SetInputDataObject(mbds)
        # cdsa = vtkCompositeDataDisplayAttributes()
        # mapper.SetCompositeDataDisplayAttributes(cdsa)
    


    
    beam_extrude = vtkPolyDataMapper()
    beam_extrude.SetInputConnection(transform_filter.GetOutputPort())
    beam_extrude_actor = vtkActor()        
    beam_extrude_actor.GetProperty().SetRepresentationToWireframe()  
    beam_extrude_actor.SetMapper(beam_extrude)
    beam_extrude_actor.GetProperty().SetLineWidth(5.0)
    beam_extrude_actor.GetProperty().EdgeVisibilityOn()
    # beam_extrude_actor.SetScale(5)
    beam_extrude_actor.GetProperty().SetColor(0.7,0.7,0.7) # (R,G,B)
    
    return beam_extrude_actor




def plotterXY_resultgraph(file_save,save_screen,screen_on,numPoints,X,Y,Xlabel,Ylabel):
    # OS Windows
    # file_name = str(os.getcwd() + '\\' +  file_dir)
    file_save = str(os.getcwd() + '/' +  file_save)
    
    chart = vtkChartXY()
    
    chart.GetAxis(0).SetTitle(Xlabel)
    chart.GetAxis(1).SetTitle(Ylabel)
    chart.GetAxis(0).GetTitleProperties().SetFontSize(20)
    chart.GetAxis(1).GetTitleProperties().SetFontSize(20)
    chart.GetAxis(0).GetLabelProperties().SetFontSize(20)
    chart.GetAxis(1).GetLabelProperties().SetFontSize(20)
    
    
    table = vtkTable()
    
    domX = vtkFloatArray()
    domX.SetName(Xlabel)

    domY = vtkFloatArray()
    domY.SetName(Ylabel)
    
    table.AddColumn(domX)
    table.AddColumn(domY)
        
    table.SetNumberOfRows(numPoints)    
    for i in range(numPoints):
        table.SetValue(i, 0, float(X[i]))
        table.SetValue(i, 1, float(Y[i]))
      
    line = chart.AddPlot(vtkChart.LINE)
    line.SetInputData(table, 0, 1)
    line.SetColor(0, 0, 255, 255)
    line.SetWidth(3.0)
    line.SetMarkerStyle(vtkPlotPoints.SQUARE)
    
    # Legend
    chart.SetShowLegend(True)
    legend = chart.GetLegend()
    legend.SetPoint(0, 20)
    
    view = vtkContextActor()
    view.GetScene().AddItem(chart)
    
    # Renderer
    renderer = vtkRenderer()
    renderer.SetBackground(1,1,1)
    renderer.AddViewProp(view)
    
    # Window
    renderer_window = vtkRenderWindow()
    renderer_window.AddRenderer(renderer)
    renderer_window.SetSize(1024,768)
    
    if save_screen == 'false':
        print(' ')
    else: 
        # screenshot code:
        im = vtkWindowToImageFilter()
        writer = vtkPNGWriter()
        im.SetInput(renderer_window)
        im.Update()
        writer.SetInputConnection(im.GetOutputPort())
        writer.SetFileName(file_save)
        writer.Write()
    
    if screen_on == 'false':
        print(' ')
    else:
        # # # Create the RendererWindowInteractor and display the vtk_file
        interactor = vtkRenderWindowInteractor()
        interactor.SetRenderWindow(renderer_window)
        interactor.RemoveObservers('LeftButtonPressEvent')
        interactor.RemoveObservers('RightButtonPressEvent')
        interactor.Initialize()
        renderer_window.Render() 
        interactor.Start()