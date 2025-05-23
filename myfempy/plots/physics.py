#!/usr/bin/env python
__doc__ = """
Physics Vtk Plot.
"""
from os import environ

environ["OMP_NUM_THREADS"] = "3"
import numpy as np
import vtk


def view_listforce(coord: np.ndarray, frcApy_vet: np.ndarray, scala_view: float):
    """vtk code"""

    coordX_force = coord[int(frcApy_vet[0]) - 1, 1]
    coordY_force = coord[int(frcApy_vet[0]) - 1, 2]
    coordZ_force = coord[int(frcApy_vet[0]) - 1, 3]
    fr_text = vtk.vtkVectorText()
    if frcApy_vet[2] == 0.0:
        height_cone = 0.0
    else:
        # 2.3*scala_view*abs((frcApy_vet[num_lf][2])/(max(abs(frcApy_vet[:,2])))) #abs((scala_view*frcApy_vet[num_lf][2])/(0.5*max(abs(frcApy_vet[:,2]))))
        height_cone = scala_view
    if frcApy_vet[1] == 0:
        dir_cone1 = (0, 0, 0)
        dir_cone2 = dir_cone1
        center_cone1 = (0, 0, 0)
        center_cone2 = center_cone1
        color_fr = (0, 0, 0)
        height_cone = 0
        fr_text.SetText(" ")
        fr_text.Update()
    if frcApy_vet[1] == 1:  # fx
        dir_cone1 = (np.sign(frcApy_vet[2]), 0, 0)
        dir_cone2 = dir_cone1
        center_cone1 = (coordX_force + height_cone / 2, coordY_force, coordZ_force)
        center_cone2 = center_cone1
        color_fr = (0, 1, 0)
    elif frcApy_vet[1] == 2:  # fy
        dir_cone1 = (0, np.sign(frcApy_vet[2]), 0)
        dir_cone2 = dir_cone1
        center_cone1 = (coordX_force, coordY_force + height_cone / 2, coordZ_force)
        center_cone2 = center_cone1
        color_fr = (0, 1, 0)
    elif frcApy_vet[1] == 3:  # fz
        dir_cone1 = (0, 0, np.sign(frcApy_vet[2]))
        dir_cone2 = dir_cone1
        center_cone1 = (coordX_force, coordY_force, coordZ_force + height_cone / 2)
        center_cone2 = center_cone1
        color_fr = (0, 1, 0)
    elif frcApy_vet[1] == 4:  # tx
        height_cone = 1.2 * height_cone
        dir_cone1 = (np.sign(frcApy_vet[2]), 0, 0)
        dir_cone2 = dir_cone1
        center_cone1 = (coordX_force + height_cone / 2, coordY_force, coordZ_force)
        center_cone2 = (coordX_force + height_cone, coordY_force, coordZ_force)
        color_fr = (1, 0, 1)
    elif frcApy_vet[1] == 5:  # ty
        height_cone = 1.2 * height_cone
        dir_cone1 = (0, np.sign(frcApy_vet[2]), 0)
        dir_cone2 = dir_cone1
        center_cone1 = (coordX_force, coordY_force + height_cone / 2, coordZ_force)
        center_cone2 = (coordX_force, coordY_force + height_cone, coordZ_force)
        color_fr = (1, 0, 1)
    elif frcApy_vet[1] == 6:  # tz
        height_cone = 1.2 * height_cone
        dir_cone1 = (0, 0, np.sign(frcApy_vet[2]))
        dir_cone2 = dir_cone1
        center_cone1 = (coordX_force, coordY_force, coordZ_force + height_cone / 2)
        center_cone2 = (coordX_force, coordY_force, coordZ_force + height_cone)
        color_fr = (1, 0, 1)
    elif frcApy_vet[1] == 15:  # massadd
        height_cone = 1.2 * scala_view
        dir_cone1 = (-1, 0, 0)
        dir_cone2 = (1, 0, 0)
        center_cone1 = (coordX_force, coordY_force, coordZ_force)
        center_cone2 = (coordX_force, coordY_force, coordZ_force)
        color_fr = (0.4, 0.1, 0.8)
    elif frcApy_vet[1] == 16:  # spring2gd
        height_cone = 1.2 * scala_view
        dir_cone1 = (0, 1, 0)
        dir_cone2 = (0, -1, 0)
        center_cone1 = (coordX_force, coordY_force - height_cone / 2, coordZ_force)
        center_cone2 = (coordX_force, coordY_force - 3 * height_cone / 2, coordZ_force)
        color_fr = (0, 0.9, 0.8)
    elif frcApy_vet[1] == 17:  # damper2gd
        # *(0.5*max(abs(frcApy_vet[:,2])))/(0.5*max(abs(frcApy_vet[:,2])))
        height_cone = 1.2 * scala_view
        dir_cone1 = (0, 1, 0)
        dir_cone2 = (0, -1, 0)
        center_cone1 = (coordX_force, coordY_force - height_cone / 2, coordZ_force)
        center_cone2 = (coordX_force, coordY_force - 3 * height_cone / 2, coordZ_force)
        color_fr = (0.9, 0.6, 0.2)
    cone1 = vtk.vtkConeSource()
    cone1.SetResolution(1)
    cone1.SetHeight(height_cone)
    cone1.SetRadius(0.15 * height_cone)
    cone1.SetCenter(center_cone1)
    cone1.SetDirection(dir_cone1)
    forcemap1 = vtk.vtkPolyDataMapper()
    forcemap1.SetInputConnection(cone1.GetOutputPort())
    fr_point_actor_cone1 = vtk.vtkActor()
    fr_point_actor_cone1.SetMapper(forcemap1)
    fr_point_actor_cone1.GetProperty().SetLineWidth(0.5)
    fr_point_actor_cone1.GetProperty().SetColor(color_fr)  # (R,G,B)
    cone2 = vtk.vtkConeSource()
    cone2.SetResolution(1)
    cone2.SetHeight(height_cone)
    cone2.SetRadius(0.15 * height_cone)
    cone2.SetCenter(center_cone2)
    cone2.SetDirection(dir_cone2)
    forcemap2 = vtk.vtkPolyDataMapper()
    forcemap2.SetInputConnection(cone2.GetOutputPort())
    fr_point_actor_cone2 = vtk.vtkActor()
    fr_point_actor_cone2.SetMapper(forcemap2)
    fr_point_actor_cone2.GetProperty().SetLineWidth(0.5)
    fr_point_actor_cone2.GetProperty().SetColor(color_fr)  # (R,G,B)
    return fr_point_actor_cone1, fr_point_actor_cone2  # , fr_text_actor


def view_bondcond_point(coord: np.ndarray, bondCond_vet: np.ndarray, scala_view: float):
    """vtk code"""

    coordX_bc = coord[int(bondCond_vet[0]) - 1, 1]
    coordY_bc = coord[int(bondCond_vet[0]) - 1, 2]
    coordZ_bc = coord[int(bondCond_vet[0]) - 1, 3]
    height_cone = 0.9 * scala_view
    color_rgb = (1, 1, 0)
    max_coordX = max(coord[:, 1])
    bc_text = vtk.vtkVectorText()
    if int(bondCond_vet[1]) == 0:  # fixed all dofs
        if coordX_bc == max_coordX:
            dir_cone = (-1, 0, 0)
            center_cone = (coordX_bc + height_cone / 2, coordY_bc, coordZ_bc)
            center_cube = (coordX_bc + 2 * height_cone / 2, coordY_bc, coordZ_bc)
        else:
            dir_cone = (1, 0, 0)
            center_cone = (coordX_bc - height_cone / 2, coordY_bc, coordZ_bc)
            center_cube = (coordX_bc - 2 * height_cone / 2, coordY_bc, coordZ_bc)
        bc_text.SetText("F")
        bc_text.Update()
        cube = vtk.vtkCubeSource()
        cube.SetXLength(1.2 * height_cone)
        cube.SetYLength(1.2 * height_cone)
        cube.SetZLength(1.2 * height_cone)
        cube.SetCenter(center_cube)
        bcmap = vtk.vtkPolyDataMapper()
        bcmap.SetInputConnection(cube.GetOutputPort())
        bc_point_actor_tdof = vtk.vtkActor()
        bc_point_actor_tdof.SetMapper(bcmap)
        bc_point_actor_tdof.GetProperty().SetLineWidth(0.5)
        bc_point_actor_tdof.GetProperty().SetColor(1, 1, 0)  # (R,G,B)

    elif int(bondCond_vet[1]) == 1:  # fixed ux dofs
        dir_cone = (1, 0, 0)
        if coordX_bc == max_coordX:
            dir_cone = (-1, 0, 0)
            center_cone = (coordX_bc + height_cone / 2, coordY_bc, coordZ_bc)
            center_cube = (coordX_bc + 3.5 * height_cone / 2, coordY_bc, coordZ_bc)
        else:
            dir_cone = (1, 0, 0)
            center_cone = (coordX_bc - height_cone / 2, coordY_bc, coordZ_bc)
            center_cube = (coordX_bc - 3.5 * height_cone / 2, coordY_bc, coordZ_bc)
        bc_text.SetText("UX")
        bc_text.Update()
        cube = vtk.vtkCubeSource()
        cube.SetXLength(0.5 * height_cone)
        cube.SetYLength(2.5 * height_cone)
        cube.SetZLength(0.5 * height_cone)
        cube.SetCenter(center_cube)
        bcmap = vtk.vtkPolyDataMapper()
        bcmap.SetInputConnection(cube.GetOutputPort())
        bc_point_actor_tdof = vtk.vtkActor()
        bc_point_actor_tdof.SetMapper(bcmap)
        bc_point_actor_tdof.GetProperty().SetLineWidth(0.5)
        bc_point_actor_tdof.GetProperty().SetColor(1, 1, 0)  # (R,G,B)

    elif int(bondCond_vet[1]) == 2:  # fixed uy dofs
        dir_cone = (0, 1, 0)
        center_cone = (coordX_bc, coordY_bc - height_cone / 2, coordZ_bc)
        center_cube = (coordX_bc, coordY_bc - 3.5 * height_cone / 2, coordZ_bc)
        bc_text.SetText("UY")
        bc_text.Update()
        cube = vtk.vtkCubeSource()
        cube.SetXLength(2.5 * height_cone)
        cube.SetYLength(0.5 * height_cone)
        cube.SetZLength(0.5 * height_cone)
        cube.SetCenter(center_cube)
        bcmap = vtk.vtkPolyDataMapper()
        bcmap.SetInputConnection(cube.GetOutputPort())
        bc_point_actor_tdof = vtk.vtkActor()
        bc_point_actor_tdof.SetMapper(bcmap)
        bc_point_actor_tdof.GetProperty().SetLineWidth(0.5)
        bc_point_actor_tdof.GetProperty().SetColor(1, 1, 0)  # (R,G,B)
    
    elif int(bondCond_vet[1]) == 3:  # fixed uz dofs
        dir_cone = (0, 0, 1)
        center_cone = (coordX_bc, coordY_bc, coordZ_bc - height_cone / 2)
        center_cube = (coordX_bc, coordY_bc, coordZ_bc - 3.5 * height_cone / 2)
        bc_text.SetText("UZ")
        bc_text.Update()
        cube = vtk.vtkCubeSource()
        cube.SetXLength(0.5 * height_cone)
        cube.SetYLength(2.5 * height_cone)
        cube.SetZLength(0.5 * height_cone)
        cube.SetCenter(center_cube)
        bcmap = vtk.vtkPolyDataMapper()
        bcmap.SetInputConnection(cube.GetOutputPort())
        bc_point_actor_tdof = vtk.vtkActor()
        bc_point_actor_tdof.SetMapper(bcmap)
        bc_point_actor_tdof.GetProperty().SetLineWidth(0.5)
        bc_point_actor_tdof.GetProperty().SetColor(1, 1, 0)  # (R,G,B)

    elif int(bondCond_vet[1]) == 4:  # fixed rx dofs
        dir_cone = (0, 0, 1)
        center_cone = (coordX_bc, coordY_bc, coordZ_bc - height_cone / 2)
        center_cube = (coordX_bc, coordY_bc, coordZ_bc - 3.5 * height_cone / 2)
        bc_text.SetText("RX")
        bc_text.Update()
        cube = vtk.vtkCubeSource()
        cube.SetXLength(0.5 * height_cone)
        cube.SetYLength(2.5 * height_cone)
        cube.SetZLength(0.5 * height_cone)
        cube.SetCenter(center_cube)
        bcmap = vtk.vtkPolyDataMapper()
        bcmap.SetInputConnection(cube.GetOutputPort())
        bc_point_actor_tdof = vtk.vtkActor()
        bc_point_actor_tdof.SetMapper(bcmap)
        bc_point_actor_tdof.GetProperty().SetLineWidth(0.5)
        bc_point_actor_tdof.GetProperty().SetColor(1, 1, 0)  # (R,G,B)

    elif int(bondCond_vet[1]) == 5:  # fixed ry dofs
        dir_cone = (1, 0, 0)
        if coordX_bc == max_coordX:
            dir_cone = (-1, 0, 0)
            center_cone = (coordX_bc + height_cone / 2, coordY_bc, coordZ_bc)
            center_cube = (coordX_bc + 3.5 * height_cone / 2, coordY_bc, coordZ_bc)
        else:
            dir_cone = (1, 0, 0)
            center_cone = (coordX_bc - height_cone / 2, coordY_bc, coordZ_bc)
            center_cube = (coordX_bc - 3.5 * height_cone / 2, coordY_bc, coordZ_bc)
        bc_text.SetText("RY")
        bc_text.Update()
        cube = vtk.vtkCubeSource()
        cube.SetXLength(0.5 * height_cone)
        cube.SetYLength(2.5 * height_cone)
        cube.SetZLength(0.5 * height_cone)
        cube.SetCenter(center_cube)
        bcmap = vtk.vtkPolyDataMapper()
        bcmap.SetInputConnection(cube.GetOutputPort())
        bc_point_actor_tdof = vtk.vtkActor()
        bc_point_actor_tdof.SetMapper(bcmap)
        bc_point_actor_tdof.GetProperty().SetLineWidth(0.5)
        bc_point_actor_tdof.GetProperty().SetColor(1, 1, 0)  # (R,G,B)

    elif int(bondCond_vet[1]) == 6:  # fixed rz dofs
        dir_cone = (0, 1, 0)
        center_cone = (coordX_bc, coordY_bc - height_cone / 2, coordZ_bc)
        center_cube = (coordX_bc, coordY_bc - 3.5 * height_cone / 2, coordZ_bc)
        bc_text.SetText("RZ")
        bc_text.Update()
        cube = vtk.vtkCubeSource()
        cube.SetXLength(2.5 * height_cone)
        cube.SetYLength(0.5 * height_cone)
        cube.SetZLength(0.5 * height_cone)
        cube.SetCenter(center_cube)
        bcmap = vtk.vtkPolyDataMapper()
        bcmap.SetInputConnection(cube.GetOutputPort())
        bc_point_actor_tdof = vtk.vtkActor()
        bc_point_actor_tdof.SetMapper(bcmap)
        bc_point_actor_tdof.GetProperty().SetLineWidth(0.5)
        bc_point_actor_tdof.GetProperty().SetColor(1, 1, 0)  # (R,G,B)

    cone = vtk.vtkConeSource()
    cone.SetResolution(1)
    cone.SetHeight(height_cone)
    cone.SetRadius(0.5 * height_cone)
    cone.SetCenter(center_cone)
    cone.SetDirection(dir_cone)
    bccmap = vtk.vtkPolyDataMapper()
    bccmap.SetInputConnection(cone.GetOutputPort())
    bc_point_actor_cone = vtk.vtkActor()
    bc_point_actor_cone.SetMapper(bccmap)
    bc_point_actor_cone.GetProperty().SetLineWidth(0.5)
    bc_point_actor_cone.GetProperty().SetColor(color_rgb)  # (R,G,B)
    bc_text_map = vtk.vtkPolyDataMapper()
    bc_text_map.SetInputConnection(bc_text.GetOutputPort())
    bc_text_actor = vtk.vtkActor()
    bc_text_actor.SetMapper(bc_text_map)
    bc_text_actor.SetScale(5, 5, 5)
    bc_text_actor.SetPosition(center_cone)
    bc_text_actor.GetProperty().SetColor(0, 0, 1)
    return bc_point_actor_cone, bc_point_actor_tdof


def view_text_point(coord: np.ndarray, coordMax: float, scala_view: float, text: str):
    """vtk code"""

    coordX = coord[0]
    coordY = coord[1]
    coordZ = coord[2]
    height_cone = 0.9 * scala_view
    scala = 0.5 * scala_view
    bc_text = vtk.vtkVectorText()
    bc_text.SetText(text[0] + "_" + text[1])
    bc_text.Update()
    center_cone = (coordX, coordY, coordZ)
    bc_text_map = vtk.vtkPolyDataMapper()
    bc_text_map.SetInputConnection(bc_text.GetOutputPort())
    bc_text_actor = vtk.vtkActor()
    bc_text_actor.SetMapper(bc_text_map)
    bc_text_actor.SetScale(scala, scala, scala)
    bc_text_actor.SetPosition(center_cone)
    bc_text_actor.GetProperty().SetColor(1.0, 1.0, 1.0)
    return bc_text_actor


def view_beam_crossSection(dimSection, typSection, coord_bcs, scala_view):
    """_summary_

    Arguments:
        dimSection -- _description_
        typSection -- _description_
        coord_bcs -- _description_

    Returns:
        _description_
    """
    scala = 0.008 * scala_view
    b = dimSection[0] * scala
    h = dimSection[1] * scala
    t = dimSection[2] * scala
    dia = dimSection[3] * scala
    Lx = np.sqrt(((coord_bcs[3] - coord_bcs[0]) ** 2))
    Ly = np.sqrt(((coord_bcs[4] - coord_bcs[1]) ** 2))
    Lz = np.sqrt(((coord_bcs[5] - coord_bcs[2]) ** 2))
    L = np.sqrt(Lx**2 + Ly**2 + Lz**2)
    if (Lx > 0.0) and (Ly == 0.0) and (Lz == 0.0):
        l = (coord_bcs[3] - coord_bcs[0]) / L
        m = (coord_bcs[1] - coord_bcs[4]) / L
        n = (coord_bcs[5] - coord_bcs[2]) / L
        rotate1 = (0, 1, 0)
        rotate2 = (0, 0, 0)
        rotate3 = (0, 0, 0)
        ang1 = 90
        ang2 = 0
        ang3 = 0
        translate = (
            coord_bcs[0] - np.cos(l * np.pi) * 0.5 * L,
            coord_bcs[1] - np.sin(m * np.pi) * 0.5 * L,
            coord_bcs[2] - np.sin(n * np.pi) * 0.5 * L,
        )
    elif (Lx == 0.0) and (Ly > 0.0) and (Lz == 0.0):
        l = (coord_bcs[3] - coord_bcs[0]) / L
        m = (coord_bcs[1] - coord_bcs[4]) / L
        n = (coord_bcs[5] - coord_bcs[2]) / L
        rotate1 = (1, 0, 0)
        rotate2 = (0, 0, 0)
        rotate3 = (0, 0, 1)
        ang1 = 90
        ang2 = 0
        ang3 = 90
        translate = (
            coord_bcs[0] - np.sin(l * np.pi) * 0.5 * L,
            coord_bcs[1] - np.cos(m * np.pi) * 0.5 * L,
            coord_bcs[2] - np.sin(n * np.pi) * 0.5 * L,
        )
    elif (Lx == 0.0) and (Ly == 0.0) and (Lz > 0.0):
        l = (coord_bcs[3] - coord_bcs[0]) / L
        m = (coord_bcs[1] - coord_bcs[4]) / L
        n = (coord_bcs[5] - coord_bcs[2]) / L
        rotate1 = (0, 0, 0)
        rotate2 = (0, 0, 0)
        rotate3 = (0, 1, 0)
        ang1 = 0
        ang2 = 0
        ang3 = -90
        translate = (
            coord_bcs[0] - np.sin(l * np.pi) * 0.5 * L,
            coord_bcs[1] - np.sin(m * np.pi) * 0.5 * L,
            coord_bcs[2] - np.cos(n * np.pi) * 0.5 * L,
        )
    elif (Lx > 0.0) and (Ly > 0.0) and (Lz == 0.0):
        l = (coord_bcs[3] - coord_bcs[0]) / L
        m = (coord_bcs[1] - coord_bcs[4]) / L
        n = (coord_bcs[5] - coord_bcs[2]) / L
        ang1 = 90
        rotate1 = (0, m, 0)
        if l < 0:
            ang2 = (np.arctan(l / m)) * (180 / np.pi)
        else:
            ang2 = (np.arctan(-l / m)) * (180 / np.pi)
        rotate2 = (l, 0, 0)
        ang3 = 0
        rotate3 = (0, 0, 0)
        translate = (
            coord_bcs[0] + l * 0.5 * L,
            coord_bcs[1] - m * 0.5 * L,
            coord_bcs[2] - n * 0.5 * L,
        )
    elif (Lx > 0.0) and (Ly == 0.0) and (Lz > 0.0):
        l = (coord_bcs[3] - coord_bcs[0]) / L
        m = (coord_bcs[1] - coord_bcs[4]) / L
        n = (coord_bcs[5] - coord_bcs[2]) / L
        ang1 = 90
        rotate1 = (0, 0, 0)
        rotate2 = (0, -l / n, 0)
        rotate3 = (0, 0, 0)
        if l < 0:
            ang2 = (np.arctan(l / n)) * (180 / np.pi)
        else:
            ang2 = (np.arctan(-l / n)) * (180 / np.pi)
        ang3 = 0
        translate = (
            coord_bcs[0] + l * 0.5 * L,
            coord_bcs[1] - m * 0.5 * L,
            coord_bcs[2] + n * 0.5 * L,
        )
    elif (Lx == 0.0) and (Ly > 0.0) and (Lz > 0.0):
        l = (coord_bcs[3] - coord_bcs[0]) / L
        m = (coord_bcs[1] - coord_bcs[4]) / L
        n = (coord_bcs[5] - coord_bcs[2]) / L
        d = np.sqrt(l**2 + m**2)
        ang1 = 0
        rotate1 = (0, 0, 0)
        ang2 = (np.arctan(m / n)) * (180 / np.pi)
        rotate2 = (d, 0, 0)
        ang3 = 0
        rotate3 = (0, 0, 0)
        translate = (
            coord_bcs[0] + l * 0.5 * L,
            coord_bcs[1] - m * 0.5 * L,
            coord_bcs[2] + n * 0.5 * L,
        )
    elif (Lx > 0.0) and (Ly > 0.0) and (Lz > 0.0):
        l = (coord_bcs[3] - coord_bcs[0]) / L
        m = (coord_bcs[1] - coord_bcs[4]) / L
        n = (coord_bcs[5] - coord_bcs[2]) / L
        ang1 = (np.arctan(l / m)) * (180 / np.pi)
        rotate1 = (0, m, 0)
        ang2 = (np.arctan(l / m)) * (180 / np.pi)
        rotate2 = (l, 0, 0)
        ang3 = 0
        rotate3 = (0, 0, 0)
        translate = (
            coord_bcs[0] + l * 0.5 * L,
            coord_bcs[1] + m * 0.5 * L,
            coord_bcs[2] + n * 0.5 * L,
        )
    if typSection == 10:  # Rectangle
        points = vtk.vtkPoints()
        points.SetNumberOfPoints(4)
        points.SetPoint(0, -b / 2, -h / 2, 0.0)
        points.SetPoint(1, b / 2, -h / 2, 0.0)
        points.SetPoint(2, b / 2, h / 2, 0.0)
        points.SetPoint(3, -b / 2, h / 2, 0.0)
        lines = vtk.vtkCellArray()
        lines.InsertNextCell(4)
        lines.InsertCellPoint(0)
        lines.InsertCellPoint(1)
        lines.InsertCellPoint(2)
        lines.InsertCellPoint(3)
        profile = vtk.vtkPolyData()
        profile.SetPoints(points)
        profile.SetPolys(lines)
        transform = vtk.vtkTransform()
        transform.Translate(translate)
        transform.RotateWXYZ(ang1, rotate1)
        transform.RotateWXYZ(ang2, rotate2)
        transform_filter = vtk.vtkTransformPolyDataFilter()
        transform_filter.SetTransform(transform)
        transform_filter.SetInputData(profile)
        transform_filter.Update()
    elif typSection == 11:  # Rectangle Tube
        points = vtk.vtkPoints()
        points.SetNumberOfPoints(12)
        points.SetPoint(0, 0.0, h / 2, 0.0)
        points.SetPoint(1, -b / 2, h / 2, 0.0)
        points.SetPoint(2, -b / 2, -h / 2, 0.0)
        points.SetPoint(3, b / 2, -h / 2, 0.0)
        points.SetPoint(4, b / 2, h / 2, 0.0)
        points.SetPoint(5, 0.0, h / 2, 0.0)
        points.SetPoint(6, 0.0, h / 2 - t, 0.0)
        points.SetPoint(7, b / 2 - t, h / 2 - t, 0.0)
        points.SetPoint(8, b / 2 - t, -h / 2 + t, 0.0)
        points.SetPoint(9, -b / 2 + t, -h / 2 + t, 0.0)
        points.SetPoint(10, -b / 2 + t, h / 2 - t, 0.0)
        points.SetPoint(11, 0.0, h / 2 - t, 0.0)
        lines = vtk.vtkCellArray()
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
        profile = vtk.vtkPolyData()
        profile.SetPoints(points)
        profile.SetPolys(lines)
        transform = vtk.vtkTransform()
        transform.Translate(translate)
        transform.RotateWXYZ(ang1, rotate1)
        transform.RotateWXYZ(ang2, rotate2)
        transform_filter = vtk.vtkTransformPolyDataFilter()
        transform_filter.SetTransform(transform)
        transform_filter.SetInputData(profile)
        transform_filter.Update()
        extrude = vtk.vtkLinearExtrusionFilter()
        extrude.SetInputConnection(transform_filter.GetOutputPort())
        extrude.SetExtrusionTypeToNormalExtrusion()
        extrude.SetVector(Lx, Ly, Lz)
        extrude.Update()
    elif typSection == 20:  # Circle
        profile = vtk.vtkDiskSource()
        profile.SetCircumferentialResolution(32)
        profile.SetOuterRadius(dia / 2)
        profile.SetInnerRadius(0)
        transform = vtk.vtkTransform()
        transform.Translate(translate)
        transform.RotateWXYZ(ang1, rotate1)
        transform.RotateWXYZ(ang2, rotate2)
        transform_filter = vtk.vtkTransformPolyDataFilter()
        transform_filter.SetTransform(transform)
        transform_filter.SetInputConnection(profile.GetOutputPort())
        transform_filter.Update()
    elif typSection == 21:  # Circle Tube
        profile = vtk.vtkDiskSource()
        profile.SetCircumferentialResolution(64)
        profile.SetOuterRadius(dia / 2)
        profile.SetInnerRadius(dia / 2 - t)
        transform = vtk.vtkTransform()
        transform.Translate(translate)
        transform.RotateWXYZ(ang1, rotate1)
        transform.RotateWXYZ(ang2, rotate2)
        transform_filter = vtk.vtkTransformPolyDataFilter()
        transform_filter.SetTransform(transform)
        transform_filter.SetInputConnection(profile.GetOutputPort())
        transform_filter.Update()
    elif typSection == 30:  # I Section
        points = vtk.vtkPoints()
        points.SetNumberOfPoints(12)
        points.SetPoint(0, b / 2, h / 2, 0)
        points.SetPoint(1, -b / 2, h / 2, 0)
        points.SetPoint(2, -b / 2, h / 2 - dia, 0)
        points.SetPoint(3, -t / 2, h / 2 - dia, 0)
        points.SetPoint(4, -t / 2, -h / 2 + dia, 0)
        points.SetPoint(5, -b / 2, -h / 2 + dia, 0)
        points.SetPoint(6, -b / 2, -h / 2, 0)
        points.SetPoint(7, b / 2, -h / 2, 0)
        points.SetPoint(8, b / 2, -h / 2 + dia, 0)
        points.SetPoint(9, t / 2, -h / 2 + dia, 0)
        points.SetPoint(10, t / 2, h / 2 - dia, 0)
        points.SetPoint(11, b / 2, h / 2 - dia, 0)
        lines = vtk.vtkCellArray()
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
        profile = vtk.vtkPolyData()
        profile.SetPoints(points)
        profile.SetPolys(lines)
        transform = vtk.vtkTransform()
        transform.Translate(translate)
        transform.RotateWXYZ(ang1, rotate1)
        transform.RotateWXYZ(ang2, rotate2)
        transform_filter = vtk.vtkTransformPolyDataFilter()
        transform_filter.SetTransform(transform)
        transform_filter.SetInputData(profile)
        transform_filter.Update()
    elif typSection == 2:  # Spring
        p0 = [0.0, 0.0, 0.0]
        p1 = [0.25 * L, 0.0, 0.0]
        p2 = [0.5 * L, 0.25 * L, 0.0]
        p3 = [0.5 * L, -0.25 * L, 0.0]
        p4 = [0.75 * L, 0.0, 0.0]
        p5 = [L, 0.0, 0.0]
        points = vtk.vtkPoints()
        points.InsertNextPoint(p0)
        points.InsertNextPoint(p1)
        points.InsertNextPoint(p2)
        points.InsertNextPoint(p3)
        points.InsertNextPoint(p4)
        points.InsertNextPoint(p5)
        lines = vtk.vtkCellArray()
        lines.InsertNextCell(6)
        lines.InsertCellPoint(0)
        lines.InsertCellPoint(1)
        lines.InsertCellPoint(2)
        lines.InsertCellPoint(3)
        lines.InsertCellPoint(4)
        lines.InsertCellPoint(5)
        profile = vtk.vtkPolyData()
        profile.SetPoints(points)
        profile.SetLines(lines)
        transform = vtk.vtkTransform()
        transform.Translate(translate)
        transform.RotateWXYZ(ang3, rotate3)
        transform_filter = vtk.vtkTransformPolyDataFilter()
        transform_filter.SetTransform(transform)
        transform_filter.SetInputData(profile)
        transform_filter.Update()
    elif typSection == 99:  # User Defined
        profile = vtk.vtkDiskSource()
        profile.SetCircumferentialResolution(32)
        profile.SetOuterRadius(scala)
        profile.SetInnerRadius(0)
        transform = vtk.vtkTransform()
        transform.Translate(translate)
        transform.RotateWXYZ(ang1, rotate1)
        transform.RotateWXYZ(ang2, rotate2)
        transform_filter = vtk.vtkTransformPolyDataFilter()
        transform_filter.SetTransform(transform)
        transform_filter.SetInputConnection(profile.GetOutputPort())
        transform_filter.Update()
    beam_extrude = vtk.vtkPolyDataMapper()
    beam_extrude.SetInputConnection(transform_filter.GetOutputPort())
    beam_extrude_actor = vtk.vtkActor()
    beam_extrude_actor.GetProperty().SetRepresentationToWireframe()
    beam_extrude_actor.SetMapper(beam_extrude)
    beam_extrude_actor.GetProperty().SetLineWidth(1.0)
    beam_extrude_actor.GetProperty().SetColor(1.0, 1.0, 1.0)  # (R,G,B)

    return beam_extrude_actor
