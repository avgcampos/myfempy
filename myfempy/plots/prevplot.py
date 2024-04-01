#!/usr/bin/env python
__doc__ = """
Plotter Prev Process
"""
from os import environ
environ['OMP_NUM_THREADS'] = '3'
import numpy as np
import vtk

from myfempy.io.iovtk import convert_to_vtk
from myfempy.plots.meshquality import MeshProp
from myfempy.plots.physics import view_beam_crossSection, view_bondcond_point, view_listforce, view_text_point
# from myfempy.utils.utils import get_version


# @profile
def preview_plot(previewset: dict, modelinfo: dict, path: str):

    # path = os.getcwd()
    plotdata = dict()
    plotdata["coord"] = modelinfo["coord"]
    plotdata["inci"] = modelinfo["inci"]
    plotdata["nodecon"] = modelinfo["nodecon"]
    plotdata["elemid"] = modelinfo["elemid"]
    plotdata["filename"] = path + "/" + previewset["RENDER"]["filename"]
    plotdata["title"] = ["UNDEFORM_MESH"]
    plotdata["solution"] = np.ones((len(modelinfo["inci"]), 1))
    plotdata["average"] = False
    plotdata["displ_POINT_DATA_val"] = []
    plotdata["displ_POINT_DATA_name"] = []
    plotdata["displ_POINT_DATA_title"] = []
    plotdata["stress_CELL_DATA_val"] = ((np.array([(modelinfo["inci"][:,2])])).T).astype(int)
    plotdata["stress_CELL_DATA_name"] = ['Model']
    plotdata["stress_CELL_DATA_title"] = ['Model']
    plotdata["stress_POINT_DATA_val"] = []
    plotdata["stress_POINT_DATA_name"] = []
    plotdata["stress_POINT_DATA_title"] = []
    plotdata["modes_POINT_DATA"] = []
    plotdata["strain_energy_CELL_DATA_val"] = []
    plotdata["strain_energy_CELL_DATA_name"] = []
    plotdata["strain_energy_CELL_DATA_title"] = []
    
    convert_to_vtk(plotdata)
    
    if "scale" in previewset["RENDER"].keys():
        previewset["RENDER"]["scale"] = (previewset["RENDER"]["scale"] / 100) * max(
            [
                max(abs(modelinfo["coord"][:, 1])),
                max(abs(modelinfo["coord"][:, 2])),
                max(abs(modelinfo["coord"][:, 3])),
            ]
        )
    else:
        previewset["RENDER"]["scale"] = 1
    
    if "lines" in previewset["RENDER"].keys():
        pass
    else:
        previewset["RENDER"]["lines"] = True
    
    if "plottags" in previewset["RENDER"].keys():
        if "point" in previewset["RENDER"]["plottags"].keys() and previewset["RENDER"]["plottags"]["point"]==True:
            previewset["regions"] = modelinfo["regions"][0]
        else:
            pass
        
        if "line" in previewset["RENDER"]["plottags"].keys() and previewset["RENDER"]["plottags"]["line"] == True:
            previewset["regions"] = modelinfo["regions"][1]
        else:
            pass
        
        if "plane" in previewset["RENDER"]["plottags"].keys() and previewset["RENDER"]["plottags"]["plane"] == True:
            previewset["regions"] = modelinfo["regions"][2]
        else:
            pass
    else:
        previewset["regions"] = [[], []]
    
    if "cs" in previewset["RENDER"].keys():
        pass
    else:
        previewset["RENDER"]["cs"] = False
        
    previewset["tabcs"] = dict()
    previewset["tabcs"]["typSection"] = []
    previewset["tabcs"]["dimSection"] = []
    for gg in range(len(modelinfo["tabgeo"])):
        previewset["tabcs"]["typSection"].append(int(modelinfo["tabgeo"][gg][-1]))
        previewset["tabcs"]["dimSection"].append(modelinfo["tabgeo"][gg][5:9])
        
    if "forces" in modelinfo.keys():
        previewset["forces"] = modelinfo["forces"]
    else:
        pass
    
    if "constrains" in modelinfo.keys():
        previewset["constrains"] = modelinfo["constrains"]
    else:
        pass
    
    previewset["coord"] = modelinfo["coord"]
    previewset["inci"] = modelinfo["inci"]
    previewset["nnode"] = modelinfo["nnode"]
    previewset["nodecon"] = modelinfo["nodecon"]
    previewset["dofs"] = modelinfo["dofs"]
        
    build_preview(previewset, path)
    
    # if "QUALITY" in previewset.keys():
    #     if previewset["QUALITY"]["show"] == False:
    #         pass
    # else:
    #     # if "lines" in previewset["QUALITY"].keys() and previewset["QUALITY"]["lines"] == True:
    #     #     pass
    #     # else:
    #     #     previewset["QUALITY"]["lines"] = False
    #     # if "scale" in previewset["QUALITY"].keys():
    #     #     previewset["QUALITY"]["scale"] = (
    #     #         previewset["QUALITY"]["scale"] / 100
    #     #     ) * max(
    #     #         [
    #     #             max(abs(modelinfo["coord"][:, 1])),
    #     #             max(abs(modelinfo["coord"][:, 2])),
    #     #             max(abs(modelinfo["coord"][:, 3])),
    #     #         ]
    #     #     )

    #     # else:
    #     #     previewset["QUALITY"]["scale"] = 1
    #     mesh = MeshProp(previewset)
    #     mesh.mesh_quality()
    
    if "LABELS" in previewset.keys() and previewset["LABELS"]["show"] == True:
        if "lines" in previewset["LABELS"].keys():
            pass
        else:
            previewset["LABELS"]["lines"] = False
        
        if "scale" in previewset["LABELS"].keys():
            previewset["LABELS"]["scale"] = (previewset["LABELS"]["scale"] / 100) * max(
                [
                    max(abs(modelinfo["coord"][:, 1])),
                    max(abs(modelinfo["coord"][:, 2])),
                    max(abs(modelinfo["coord"][:, 3])),
                ]
            )
        else:
            previewset["LABELS"]["scale"] = 1
        mesh = MeshProp(previewset, path)
        mesh.mesh_numbering()
    else:
        pass

# @profile
def build_preview(previewset: dict, path):

    # path = os.getcwd()
    file_name = str(path + "/" + previewset["RENDER"]["filename"] + ".vtk")
    renderer = vtk.vtkRenderer()
    renderer_window = vtk.vtkRenderWindow()
    renderer_window.SetSize(600, 480)
    renderer.SetBackground(1.0, 1.0, 1.0)
    reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName(file_name)
    reader.Update()  # Needed because of GetScalarRange
    output = reader.GetOutput()
    output_port = reader.GetOutputPort()
    scalar_range = output.GetScalarRange()
    lut = vtk.vtkLookupTable()
    lut.SetNumberOfColors(8)
    lut.SetHueRange(0.6, 1)  # jet color
    lut.Build()
    mapper = vtk.vtkDataSetMapper()
    mapper.SetLookupTable(lut)
    mapper.SetInputConnection(output_port)
    mapper.SetScalarRange(scalar_range)
    scala_view = previewset["RENDER"]["scale"]
    dimfrlist = 0
    if "forces" in previewset.keys():
        
        key_list_fc = list(previewset["dofs"]['f'].keys())
        val_list_fc = list(previewset["dofs"]['f'].values())
        
        dimfrlist = previewset["forces"].shape[0]
        for num_lf in range(dimfrlist):
            frcApy_vet = previewset['forces'][[num_lf]][0]
            frcApy_vet[1] = __setLoadDof(key_list_fc[val_list_fc.index(frcApy_vet[1])])
            exec(
                f'fr_point_actor_cone1_{num_lf},fr_point_actor_cone2_{num_lf} = view_listforce(previewset["coord"],frcApy_vet,scala_view)'
            )
    dimbclist = 0
    if "constrains" in previewset.keys():
        
        key_list_bc = list(previewset["dofs"]['d'].keys())
        val_list_bc = list(previewset["dofs"]['d'].values())
        
        dimbclist = previewset["constrains"].shape[0]
        for num_bc in range(dimbclist):
            bondCond_vet = previewset['constrains'][[num_bc]][0]
            if int(bondCond_vet[1]) == 0:
                pass
            else:
                bondCond_vet[1] = __setBCDof(key_list_bc[val_list_bc.index(int(bondCond_vet[1]))])
            exec(
                f'bc_point_actor_cone_{num_bc}, bc_point_actor_tdof_{num_bc} = view_bondcond_point(previewset["coord"],bondCond_vet,scala_view)'
            )
    objs = 0
    coordMax = [
        max(previewset["coord"][:, 1]),
        max(previewset["coord"][:, 2]),
        max(previewset["coord"][:, 3]),
    ]
    # for num_objs in range(len(previewset["regions"])):
    for num_objs in range(len(previewset["regions"][1])):
        text = [previewset["regions"][0], str(previewset["regions"][1][num_objs][0])]
        coord = [
            sum(previewset["coord"][previewset["regions"][1][num_objs][1] - 1, 1])
            / len(previewset["coord"][previewset["regions"][1][num_objs][1] - 1, 1]),
            sum(previewset["coord"][previewset["regions"][1][num_objs][1] - 1, 2])
            / len(previewset["coord"][previewset["regions"][1][num_objs][1] - 1, 2]),
            sum(previewset["coord"][previewset["regions"][1][num_objs][1] - 1, 3])
            / len(previewset["coord"][previewset["regions"][1][num_objs][1] - 1, 2]),
        ]
        exec(
            f"bc_text_actor_{objs} = view_text_point(coord, coordMax, scala_view, text)"
        )
        objs += 1

    i_node = 0
    if previewset["RENDER"]["cs"] == True:
        dimcs = len(previewset["tabcs"]["typSection"])
        for num_bcs in range(dimcs):
            inci_bcs = previewset["inci"][
                np.where(previewset["inci"][:, 3] == num_bcs + 1), :
            ][0]
            typSection = int(previewset["tabcs"]["typSection"][num_bcs])
            dimSection = previewset["tabcs"]["dimSection"][num_bcs]
            for node_bcs in range(len(inci_bcs)):
                noi = int(inci_bcs[node_bcs, 4])
                noj = int(inci_bcs[node_bcs, 5])
                noix = previewset["coord"][noi - 1, 1]
                noiy = previewset["coord"][noi - 1, 2]
                noiz = previewset["coord"][noi - 1, 3]
                nojx = previewset["coord"][noj - 1, 1]
                nojy = previewset["coord"][noj - 1, 2]
                nojz = previewset["coord"][noj - 1, 3]
                coord_bcs = [noix, noiy, noiz, nojx, nojy, nojz]
                exec(
                    f"beam_extrude_actor_{i_node} = view_beam_crossSection(dimSection, typSection, coord_bcs)"
                )
                i_node += 1
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().EdgeVisibilityOn()
    if previewset["RENDER"]["lines"] == False:
        actor.GetProperty().EdgeVisibilityOff()
    else:
        actor.GetProperty().SetLineWidth(3.0)
    text_logo = vtk.vtkTextActor()
    text_logo.SetInput(
        "MYFEMPY " + ' < PreProc--Model >\nPress "w" to wireframe view \nPress "q" to exit\continue'
    )
    txtprop = text_logo.GetTextProperty()
    txtprop.SetFontFamilyToArial()
    txtprop.SetFontSize(20)
    txtprop.SetColor(0, 0, 0)
    text_logo.SetDisplayPosition(10, 400)
    colors = vtk.vtkNamedColors()
    backgroundColor = colors.GetColor3d("DarkSlateGray")
    actorColor = colors.GetColor3d("Tomato")
    axis1Color = colors.GetColor3d("Salmon")
    axis2Color = colors.GetColor3d("PaleGreen")
    axis3Color = colors.GetColor3d("LightSkyBlue")
    cubeAxesActor = vtk.vtkCubeAxesActor()
    cubeAxesActor.SetBounds(actor.GetBounds())
    cubeAxesActor.SetCamera(renderer.GetActiveCamera())
    cubeAxesActor.GetTitleTextProperty(0).SetFontSize(25)
    cubeAxesActor.GetLabelTextProperty(0).SetColor(axis1Color)
    cubeAxesActor.GetLabelTextProperty(1).SetColor(axis2Color)
    cubeAxesActor.GetLabelTextProperty(2).SetColor(axis3Color)
    renderer.AddActor(actor)
    renderer.AddActor(text_logo)
    renderer.AddActor(cubeAxesActor)
    for ipt in range(dimbclist):
        exec(f"renderer.AddActor(bc_point_actor_cone_{ipt})")
        exec(f"renderer.AddActor(bc_point_actor_tdof_{ipt})")
    for bcs in range(i_node):
        exec(f"renderer.AddActor(beam_extrude_actor_{bcs})")
    for ff in range(dimfrlist):
        exec(f"renderer.AddActor(fr_point_actor_cone1_{ff})")
        exec(f"renderer.AddActor(fr_point_actor_cone2_{ff})")
    if "plottags" in previewset["RENDER"].keys():
        for txt in range(objs):
            exec(f"renderer.AddActor(bc_text_actor_{txt})")
    else:
        pass
    renderer_window.AddRenderer(renderer)
    renderer.ResetCamera()
    if previewset["RENDER"]["savepng"] == False:
        pass
    else:
        im = vtk.vtkWindowToImageFilter()
        writer = (
            vtk.vtkPNGWriter()
        )  # vtkSTLWriter()#vtkVRMLExporter()#vtkPolyDataWriter()#
        im.SetInput(renderer_window)
        im.Update()
        writer.SetInputConnection(im.GetOutputPort())
        writer.SetFileName((path + "/" + previewset["RENDER"]["filename"] + ".png"))
        writer.Write()
    if previewset["RENDER"]["show"] == False:
        pass
    else:
        interactor = vtk.vtkRenderWindowInteractor()
        interactor.SetRenderWindow(renderer_window)
        interactor.RemoveObservers("RightButtonPressEvent")
        interactor.Initialize()
        renderer_window.Render()
        interactor.Start()
        interactor.GetRenderWindow().Finalize()

def __setLoadDof(forcedof):           
    fdoftype = {
        "fx": 1,
        "fy": 2,
        "fz": 3,
        "tx": 4,
        "ty": 5,
        "tz": 6,
        "masspoint": 15,
        "spring2ground": 16,
        "damper2ground": 17,
        "cg": 18,
    }
    return fdoftype[forcedof]

def __setBCDof(bcdof):           
    bcdoftype = {
        "ux": 1,
        "uy": 2,
        "uz": 3,
        "rx": 4,
        "ry": 5,
        "rz": 6,
        "all": 0,
    }
    return bcdoftype[bcdof]