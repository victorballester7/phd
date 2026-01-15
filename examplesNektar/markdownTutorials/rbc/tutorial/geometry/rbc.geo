Geometry.LineNumbers = 1;
Geometry.Color.Points = Red;
General.Color.Text = Blue;
Mesh.Color.Points = {255,0,0};

Lx = 2.02;
Ly = 1.0;
Ne = 11;
Point(1) = {0, 0, 0};
Point(2) = {Lx, 0, 0};
Point(3) = {Lx, Ly, 0};
Point(4) = {0, Ly, 0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};


Transfinite Line {1, 3} = Ne Using Progression 1;
Transfinite Line {4, 2} = Ne Using Bump 0.2;

Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Transfinite Surface{1};
Recombine Surface {1};

Physical Surface(0) = {1};
Physical Line(1) = {1};//Bottom Plate
Physical Line(2) = {3};//Top Plate
Physical Line(3) = {4};//Left surface
Physical Line(4) = {2};//Right surface
