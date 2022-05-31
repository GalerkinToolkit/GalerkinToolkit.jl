//+
Point(1) = {-1, -1, 0, 0.25};
//+
Point(2) = {0, -1, 0, 2.0};
//+
Point(3) = {0, -0, 0, 0.25};
//+
Point(4) = {-1, 0, 0, 2.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Curve Loop(1) = {4, 1, 2, 3};
//+
Plane Surface(1) = {1};

Periodic Line {3} = {1} Translate {0, 1, 0};
//+
Physical Point("PL") = {1};
//+
Physical Point("PR") = {2};
//+
Physical Point("PR") += {3};
//+
Physical Point("PL") += {4};
//+
Physical Curve("P") = {1};
//+
Physical Curve("R") = {2};
//+
Physical Curve("P") += {3};
//+
Physical Curve("L") = {4};
//+
Physical Surface("S") = {1};
