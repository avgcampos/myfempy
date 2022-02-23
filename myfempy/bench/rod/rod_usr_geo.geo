// GMSH GEOMETRY FILE FROM MYFEMPY
Point(1) = {0,0,0};
Point(2) = {0.03,0,0};
Point(3) = {0.07,0,0};
Point(4) = {0.13,0,0};
Point(5) = {0.17,0,0};
Point(6) = {0.23,0,0};
Point(7) = {0.27,0,0};
Point(8) = {0.30,0,0};
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Transfinite Curve {1,2,3,4,5,6,7} = 101 Using Progression 1;
