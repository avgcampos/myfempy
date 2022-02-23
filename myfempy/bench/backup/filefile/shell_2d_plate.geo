// Gmsh project created on Tue May 12 20:42:50 2020
SetFactory("OpenCASCADE");

Point(1) = {0, 0, 0};
Point(2) = {200, 0, 0};
Point(3) = {200, 120,  0};
Point(4) = {0,  120, 0};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Curve Loop(1) = {3,2,1,4};
Plane Surface(1) = {1};
Characteristic Length {1,2,3,4} = 50;
Transfinite Surface {1};
Mesh 2;
//Save "mesh_from_GMSH_ex3.msh1";