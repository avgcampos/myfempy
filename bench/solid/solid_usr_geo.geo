// GMSH GEOMETRY FILE FROM MYFEMPY
Point(1) = {0,0,0};
Point(2) = {200,0,0};
Point(3) = {200,100,0};
Point(4) = {0,100,0};
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Curve Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};
Characteristic Length {1,2,3,4} = 50;
//FACE MAPPING AND EXTRUDE SURFACE
Transfinite Surface {1};
Recombine Surface {1};
Extrude {0, 0, 100} {Surface{1};Layers{2};Recombine;};
