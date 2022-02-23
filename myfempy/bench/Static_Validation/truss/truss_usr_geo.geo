// GMSH GEOMETRY FILE FROM MYFEMPY
Point(1) = {0,0,0};
Point(2) = {20,0,0};
Point(3) = {0,20,0};
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,1};
Transfinite Curve {1,2,3} = 2 Using Progression 1;
