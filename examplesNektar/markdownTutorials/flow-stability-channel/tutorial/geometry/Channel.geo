Point(1) = {3.141592653589793, -1, 0, 1e+22};
Point(2) = {3.141592653589793, 1, 0, 1e+22};
Point(3) = {-3.141592653589793, 1, 0, 1e+22};
Point(4) = {-3.141592653589793, -1, 0, 1e+22};
Line(1) = {4, 1};
Line(2) = {1, 2};
Line(3) = {2, 3};
Line(4) = {3, 4};
Line Loop(5) = {2, 3, 4, 1};
Plane Surface(6) = {5};
Transfinite Line {1, 3} = 9 Using Progression 1;
Transfinite Line {4, 2} = 7 Using Bump 0.2;
Transfinite Surface{6};
Recombine Surface {6};

Physical Surface(0) = {6};
Physical Line(1) = {1, 3};
Physical Line(2) = {4};
Physical Line(3) = {2};

