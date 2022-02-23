// GMSH GEOMETRY FILE FROM MYFEMPY
Point(1) = {0,0,0};
Point(2) = {0,20,0};
Point(3) = {20,20,0};
Line(1) = {1,2};
Line(2) = {2,3};
Transfinite Curve {1,2} = 21 Using Progression 1;
