// GMSH GEOMETRY FILE FROM MYFEMPY *SECTION RECTANGLE*
SetFactory("OpenCASCADE");
Point(1) = {-25.0,-50.0,0.0};
Point(2) = {25.0,-50.0,0.0};
Point(3) = {25.0,50.0,0.0};
Point(4) = {-25.0,50.0,0.0};
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Curve Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};
Transfinite Curve {1,2,3,4} = 10 Using Progression 1;
Transfinite Surface {1};
