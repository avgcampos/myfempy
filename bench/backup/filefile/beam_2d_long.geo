// Gmsh project created on Tue May 12 20:42:50 2020
SetFactor y("OpenCASCADE");
Point(1) = {0, 0, 0};
Point(2) = {4000, 0, 0};
Line(1) = {1, 2};

Transfinite Curve {1,2} = 161 Using Progression 1;
//Mesh 1;
//Save "mesh_from_GMSH_ex3.msh1";