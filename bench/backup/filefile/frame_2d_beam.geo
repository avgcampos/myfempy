// Gmsh project created on Tue May 12 20:42:50 2020
SetFactory("OpenCASCADE");
Point(1) = {0, 0, 0};
Point(2) = {0, 4000, 0};
Line(1) = {1, 2};

Characteristic Length {1,2} = 5;
//Mesh 1;
//Save "mesh_from_GMSH_ex3.msh1";