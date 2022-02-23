// GMSH GEOMETRY FILE FROM MYFEMPY
Point(1) = {0,0,0};
Point(2) = {3,0,0};
Point(3) = {0,0,-3};
Point(4) = {0,-4,0};
Line(1) = {1,2};
Line(2) = {1,3};
Line(3) = {1,4};
Transfinite Curve {1,2,3} = 2 Using Progression 1;
