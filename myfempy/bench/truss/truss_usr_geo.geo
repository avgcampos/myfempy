// GMSH GEOMETRY FILE FROM MYFEMPY
Point(1) = {0,0,0};
Point(2) = {1000,0,0};
Point(3) = {0,1000,0};
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,1};
Transfinite Curve {1,2,3} = 51 Using Progression 1;
