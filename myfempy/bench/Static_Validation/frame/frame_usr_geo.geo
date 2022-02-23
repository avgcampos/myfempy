// GMSH GEOMETRY FILE FROM MYFEMPY
Point(1) = {0,0,0};
Point(2) = {20,0,20};
Point(3) = {0,20,0};
Line(1) = {1,3};
Line(2) = {1,2};
Line(3) = {2,3};
Transfinite Curve {1,2,3} = 120 Using Progression 1;
