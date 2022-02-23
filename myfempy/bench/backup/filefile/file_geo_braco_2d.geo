// Gmsh project created on Tue Jan 19 19:01:40 2021
SetFactory("OpenCASCADE");
Point(1) = {0, 6, 0, 1.0};
Point(2) = {0, 12, 0, 1.0};
Point(3) = {12, 6, 0, 1.0};
Point(4) = {12, 12, 0, 1.0};
Point(5) = {18, 0, 0, 1.0};
Point(6) = {100, 0, 0, 1.0};
Point(7) = {100, 6, 0, 1.0};
Point(8) = {18, 6, 0, 1.0};//+
Line(1) = {1, 3};
Line(2) = {3, 5};
Line(3) = {5, 6};
Line(4) = {6, 7};
Line(5) = {7, 8};
Line(6) = {8, 4};
Line(7) = {4, 2};
Line(8) = {1, 2};

Curve Loop(1) = {1,2,3,4,5,6,7,8};
Plane Surface(1) = {1};
Characteristic Length {1,2,3,4,5,6,7,8} = 0.5;
Mesh 2;
Save "mesh_from_GMSH_braco12.msh1";