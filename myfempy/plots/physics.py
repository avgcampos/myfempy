#!/usr/bin/env python
__doc__ = """
Physics Vtk Plot.
"""

try:
    import vtk
    VTK_AVAILABLE = True
except ImportError:
    vtk = None
    VTK_AVAILABLE = False

import numpy as np

__docformat__ = "google"

__doc__ = """
This Python file is part of myfempy project.

myfempy is a python package based on finite element method to multiphysics
analysis. The code is open source and *intended for educational and scientific
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

if VTK_AVAILABLE:
    # versao nova
    def view_listforce(coord: np.ndarray, frcApy_vet: np.ndarray, scala_view: float):
        """Versão otimizada para plotagem de forças/momentos."""
        node_idx = int(frcApy_vet[0]) - 1
        coordX_force, coordY_force, coordZ_force = coord[node_idx, 1:4]
        
        val_type = int(frcApy_vet[1])
        val_magnitude = frcApy_vet[2]
        
        if val_magnitude == 0.0 or val_type == 0:
            empty_actor = vtk.vtkActor()
            return empty_actor, empty_actor

        height_cone = scala_view
        sign_mag = np.sign(val_magnitude)

        # Configurações padrão por tipo de força/momento
        # Formato: (h_mult, dir1, dir2, center1, center2, color)
        if 1 <= val_type <= 3:  # fx, fy, fz
            d = [0, 0, 0]
            d[val_type - 1] = sign_mag
            c = [coordX_force, coordY_force, coordZ_force]
            c[val_type - 1] += height_cone / 2
            config = (1.0, tuple(d), tuple(d), tuple(c), tuple(c), (0, 1, 0))
            
        elif 4 <= val_type <= 6:  # tx, ty, tz
            height_cone *= 1.2
            idx = val_type - 4
            d = [0, 0, 0]
            d[idx] = sign_mag
            c1, c2 = [coordX_force, coordY_force, coordZ_force], [coordX_force, coordY_force, coordZ_force]
            c1[idx] += height_cone / 2
            c2[idx] += height_cone
            config = (1.2, tuple(d), tuple(d), tuple(c1), tuple(c2), (1, 0, 1))
            
        elif val_type == 15:  # massadd
            height_cone *= 1.2
            c = (coordX_force, coordY_force, coordZ_force)
            config = (1.2, (-1, 0, 0), (1, 0, 0), c, c, (0.4, 0.1, 0.8))
            
        elif val_type == 16:  # spring2gd
            height_cone *= 1.2
            c1 = (coordX_force, coordY_force - height_cone / 2, coordZ_force)
            c2 = (coordX_force, coordY_force - 1.5 * height_cone, coordZ_force)
            config = (1.2, (0, 1, 0), (0, -1, 0), c1, c2, (0, 0.9, 0.8))
            
        elif val_type == 17:  # damper2gd
            height_cone *= 1.2
            c1 = (coordX_force, coordY_force - height_cone / 2, coordZ_force)
            c2 = (coordX_force, coordY_force - 1.5 * height_cone, coordZ_force)
            config = (1.2, (0, 1, 0), (0, -1, 0), c1, c2, (0.9, 0.6, 0.2))
        else:
            empty_actor = vtk.vtkActor()
            return empty_actor, empty_actor

        h_scale, dir1, dir2, center1, center2, color_fr = config

        def _create_cone(center, direction, h):
            cone = vtk.vtkConeSource()
            cone.SetResolution(3)
            cone.SetHeight(h)
            cone.SetRadius(0.15 * h)
            cone.SetCenter(center)
            cone.SetDirection(direction)
            
            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputConnection(cone.GetOutputPort())
            
            actor = vtk.vtkActor()
            actor.SetMapper(mapper)
            actor.GetProperty().SetLineWidth(0.5)
            actor.GetProperty().SetColor(color_fr)
            actor.GetProperty().SetOpacity(0.6)
            return actor

        fr_point_actor_cone1 = _create_cone(center1, dir1, height_cone)
        fr_point_actor_cone2 = _create_cone(center2, dir2, height_cone)

        return fr_point_actor_cone1, fr_point_actor_cone2


    def view_bondcond_point(coord: np.ndarray, bondCond_vet: np.ndarray, scala_view: float):
        """Versão otimizada para condições de contorno."""
        node_idx = int(bondCond_vet[0]) - 1
        coordX_bc, coordY_bc, coordZ_bc = coord[node_idx, 1:4]
        height_cone = 0.9 * scala_view
        color_rgb = (1, 1, 0)
        
        max_coords = np.max(coord[:, 1:4], axis=0)
        bc_type = int(bondCond_vet[1])
        
        bc_text = vtk.vtkVectorText()
        cube = vtk.vtkCubeSource()
        cube.SetLineWidth = 0.5  # Mantido compatibilidade
        
        if bc_type == 0:
            at_max = coordX_bc == max_coords[0]
            dir_cone = (-1, 0, 0) if at_max else (1, 0, 0)
            offset = 1 if at_max else -1
            center_cone = (coordX_bc + offset * height_cone / 2, coordY_bc, coordZ_bc)
            center_cube = (coordX_bc + offset * height_cone, coordY_bc, coordZ_bc)
            bc_text.SetText("F")
            cube.SetXLength(1.2 * height_cone)
            cube.SetYLength(1.2 * height_cone)
            cube.SetZLength(1.2 * height_cone)
            
        elif bc_type in (1, 5):  # UX ou RY
            at_max = coordX_bc == max_coords[0] if bc_type == 1 else True
            dir_cone = (-1, 0, 0) if at_max else (1, 0, 0)
            offset = 1 if at_max else -1
            center_cone = (coordX_bc + offset * height_cone / 2, coordY_bc, coordZ_bc)
            center_cube = (coordX_bc + offset * 1.75 * height_cone, coordY_bc, coordZ_bc)
            bc_text.SetText("UX" if bc_type == 1 else "RY")
            cube.SetXLength(0.5 * height_cone)
            cube.SetYLength(2.5 * height_cone)
            cube.SetZLength(0.5 * height_cone)
            
        elif bc_type in (2, 6):  # UY ou RZ
            at_max = coordY_bc == max_coords[1] if bc_type == 2 else False
            dir_cone = (0, -1, 0) if at_max else (0, 1, 0)
            offset = 1 if at_max else -1
            center_cone = (coordX_bc, coordY_bc + offset * height_cone / 2, coordZ_bc)
            center_cube = (coordX_bc, coordY_bc + offset * 1.75 * height_cone, coordZ_bc)
            bc_text.SetText("UY" if bc_type == 2 else "RZ")
            cube.SetXLength(2.5 * height_cone)
            cube.SetYLength(0.5 * height_cone)
            cube.SetZLength(0.5 * height_cone)
            
        elif bc_type in (3, 4):  # UZ ou RX
            at_max = coordZ_bc == max_coords[2] if bc_type == 3 else False
            dir_cone = (0, 0, -1) if at_max else (0, 0, 1)
            offset = 1 if at_max else -1
            center_cone = (coordX_bc, coordY_bc, coordZ_bc + offset * height_cone / 2)
            center_cube = (coordX_bc, coordY_bc, coordZ_bc + offset * 1.75 * height_cone)
            bc_text.SetText("UZ" if bc_type == 3 else "RX")
            cube.SetXLength(0.5 * height_cone)
            cube.SetYLength(2.5 * height_cone)
            cube.SetZLength(0.5 * height_cone)

        bc_text.Update()
        cube.SetCenter(center_cube)
        
        bcmap = vtk.vtkPolyDataMapper()
        bcmap.SetInputConnection(cube.GetOutputPort())
        bc_point_actor_tdof = vtk.vtkActor()
        bc_point_actor_tdof.SetMapper(bcmap)
        bc_point_actor_tdof.GetProperty().SetLineWidth(0.5)
        bc_point_actor_tdof.GetProperty().SetColor(1, 1, 0)
        bc_point_actor_tdof.GetProperty().SetOpacity(0.6)

        cone = vtk.vtkConeSource()
        cone.SetResolution(3)
        cone.SetHeight(height_cone)
        cone.SetRadius(0.5 * height_cone)
        cone.SetCenter(center_cone)
        cone.SetDirection(dir_cone)
        
        bccmap = vtk.vtkPolyDataMapper()
        bccmap.SetInputConnection(cone.GetOutputPort())
        bc_point_actor_cone = vtk.vtkActor()
        bc_point_actor_cone.SetMapper(bccmap)
        bc_point_actor_cone.GetProperty().SetLineWidth(0.5)
        bc_point_actor_cone.GetProperty().SetColor(color_rgb)
        bc_point_actor_cone.GetProperty().SetOpacity(0.6)

        return bc_point_actor_cone, bc_point_actor_tdof


    def view_text_point(coord: np.ndarray, coordMax: float, scala_view: float, text: list):
        """Versão limpa para tags de texto de regiões."""
        scala = 0.5 * scala_view
        bc_text = vtk.vtkVectorText()
        bc_text.SetText(f"{text[0]}_{text[1]}")
        bc_text.Update()
        
        bc_text_map = vtk.vtkPolyDataMapper()
        bc_text_map.SetInputConnection(bc_text.GetOutputPort())
        
        bc_text_actor = vtk.vtkActor()
        bc_text_actor.SetMapper(bc_text_map)
        bc_text_actor.SetScale(scala, scala, scala)
        bc_text_actor.SetPosition(tuple(coord))
        bc_text_actor.GetProperty().SetColor(1.0, 0.0, 0.0)
        
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
        b = dimSection[0] #* scala
        h = dimSection[1] #* scala
        t = dimSection[2] #* scala
        d = dimSection[3] #* scala
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
            profile.SetOuterRadius(d / 2)
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
            profile.SetOuterRadius(d / 2)
            profile.SetInnerRadius(d / 2 - t)
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
            points.SetPoint(2, -b / 2, h / 2 - d, 0)
            points.SetPoint(3, -t / 2, h / 2 - d, 0)
            points.SetPoint(4, -t / 2, -h / 2 + d, 0)
            points.SetPoint(5, -b / 2, -h / 2 + d, 0)
            points.SetPoint(6, -b / 2, -h / 2, 0)
            points.SetPoint(7, b / 2, -h / 2, 0)
            points.SetPoint(8, b / 2, -h / 2 + d, 0)
            points.SetPoint(9, t / 2, -h / 2 + d, 0)
            points.SetPoint(10, t / 2, h / 2 - d, 0)
            points.SetPoint(11, b / 2, h / 2 - d, 0)
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
        elif typSection == 40:  # T Section
                yc = (b*d**2 + t*(h-d)*(2*d+(h-d)))/(2*(d*b+t*(h-d)))
                points = vtk.vtkPoints()
                points.SetNumberOfPoints(8)
                points.SetPoint(0, t/2, -h + yc, 0)
                points.SetPoint(1, t/2, yc - d, 0)
                points.SetPoint(2, b/2, yc - d, 0)
                points.SetPoint(3, b/2, yc, 0)
                points.SetPoint(4, -b/2, yc, 0)
                points.SetPoint(5, -b/2, yc - d, 0)
                points.SetPoint(6, -t/2, yc - d, 0)
                points.SetPoint(7, -t/2, -h + yc, 0)
                lines = vtk.vtkCellArray()
                lines.InsertNextCell(8)
                lines.InsertCellPoint(0)
                lines.InsertCellPoint(1)
                lines.InsertCellPoint(2)
                lines.InsertCellPoint(3)
                lines.InsertCellPoint(4)
                lines.InsertCellPoint(5)
                lines.InsertCellPoint(6)
                lines.InsertCellPoint(7)
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
        elif typSection == 50:  # C Section
                xc = (1/(t*h + 2*d*(b-t))) * (((h - 2*d)*t**3)/2 + d*b**2)
                points = vtk.vtkPoints()
                points.SetNumberOfPoints(8)
                points.SetPoint(0, b - xc, - h/2, 0)
                points.SetPoint(1, b - xc, -h/2 + d, 0)
                points.SetPoint(2, -xc + t, -h/2 + d, 0)
                points.SetPoint(3, -xc + t, h/2 - d, 0)
                points.SetPoint(4, b - xc, h/2 - d, 0)
                points.SetPoint(5, b - xc, h/2, 0)
                points.SetPoint(6, -xc, h/2, 0)
                points.SetPoint(7, -xc, -h/2, 0)
                lines = vtk.vtkCellArray()
                lines.InsertNextCell(8)
                lines.InsertCellPoint(0)
                lines.InsertCellPoint(1)
                lines.InsertCellPoint(2)
                lines.InsertCellPoint(3)
                lines.InsertCellPoint(4)
                lines.InsertCellPoint(5)
                lines.InsertCellPoint(6)
                lines.InsertCellPoint(7)
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
        elif typSection == 60:  # L Section
                xc = t/(2*((h + b - t)*t)) * (b**2 + h*t -t**2)
                yc = t/(2*((h + b - t)*t)) * (h**2 + b*t -t**2)
                points = vtk.vtkPoints()
                points.SetNumberOfPoints(6)
                points.SetPoint(0, b - xc, - yc, 0)
                points.SetPoint(1, b - xc, -yc + t, 0)
                points.SetPoint(2, -xc + t, -yc + t, 0)
                points.SetPoint(3, -xc + t, h - yc, 0)
                points.SetPoint(4, -xc, h - yc, 0)
                points.SetPoint(5, -xc, -yc, 0)
                lines = vtk.vtkCellArray()
                lines.InsertNextCell(6)
                lines.InsertCellPoint(0)
                lines.InsertCellPoint(1)
                lines.InsertCellPoint(2)
                lines.InsertCellPoint(3)
                lines.InsertCellPoint(4)
                lines.InsertCellPoint(5)
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
            profile.SetOuterRadius(scala_view)
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