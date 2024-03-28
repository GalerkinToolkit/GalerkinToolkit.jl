// Puzzle piece cube with periodic surfaces
// TODO: Check physical groups!

//----------
// Geometry
//----------

SetFactory("OpenCASCADE");
Box(1) = {0, 0, 0, 1, 1, 1};
//+
Sphere(2) = {0.5, 0, 1, 0.25, -Pi/2, Pi/2, 2*Pi};
//+
Sphere(3) = {0.5, 1, 1, 0.25, -Pi/2, Pi/2, 2*Pi};
//+
Sphere(4) = {0.5, 1, 0, 0.25, -Pi/2, Pi/2, 2*Pi};
//+
Sphere(5) = {0.5, 0, 0, 0.25, -Pi/2, Pi/2, 2*Pi};
//+
Sphere(6) = {1, 0, 0.5, 0.25, -Pi/2, Pi/2, 2*Pi};
//+
Sphere(7) = {1, 1, 0.5, 0.25, -Pi/2, Pi/2, 2*Pi};
//+
Sphere(8) = {0, 1, 0.5, 0.25, -Pi/2, Pi/2, 2*Pi};
//+
Sphere(9) = {0, 0, 0.5, 0.25, -Pi/2, Pi/2, 2*Pi};
//+
Sphere(10) = {0, 0.5, 1, 0.25, -Pi/2, Pi/2, 2*Pi};
//+
Sphere(11) = {1, 0.5, 1, 0.25, -Pi/2, Pi/2, 2*Pi};
//+
Sphere(12) = {1, 0.5, 0, 0.25, -Pi/2, Pi/2, 2*Pi};
//+
Sphere(13) = {0, 0.5, 0, 0.25, -Pi/2, Pi/2, 2*Pi};
//+
BooleanDifference{ Volume{1}; Delete; }{ Volume{2}; Volume{11}; Volume{10}; Volume{9}; Volume{3}; Volume{13}; Volume{8}; Volume{4}; Volume{7}; Volume{12}; Volume{6}; Volume{5}; Delete; }

//----------
// Periodicity
//----------
// Back face copy of front face
Periodic Surface {9} = {4} Translate {0, 0, -1};
// Top face copy of bottom face
Periodic Surface {7} = {2} Translate {0, 1, 0};
// Right face copy of left face
Periodic Surface {12} = {1} Translate {1, 0, 0};
//+
Physical Point("0-face-1", 79) = {1};
//+
Physical Point("0-face-2", 80) = {16};
//+
Physical Point("0-face-3", 81) = {11};
//+
Physical Point("0-face-4", 82) = {32};
//+
Physical Point("0-face-5", 83) = {4};
//+
Physical Point("0-face-6", 84) = {19};
//+
Physical Point("0-face-7", 85) = {8};
//+
Physical Point("0-face-8", 86) = {28};
//+
Physical Curve("1-face-1", 87) = {24, 15};
//+
Physical Curve("1-face-2", 88) = {48, 41};
//+
Physical Curve("1-face-3", 89) = {11, 14};
//+
Physical Curve("1-face-4", 90) = {56, 58};
//+
Physical Curve("1-face-5", 91) = {22, 19};
//+
Physical Curve("1-face-6", 92) = {31, 33};
//+
Physical Curve("1-face-7", 93) = {7, 4};
//+
Physical Curve("1-face-8", 94) = {36, 34};
//+
Physical Curve("1-face-9", 95) = {3, 1};
//+
Physical Curve("1-face-10", 96) = {8, 10};
//+
Physical Curve("1-face-11", 97) = {18, 16};
//+
Physical Curve("1-face-12", 98) = {44, 42};
//+
Physical Surface("2-face-1", 99) = {9};
//+
Physical Surface("2-face-2", 100) = {4};
//+
Physical Surface("2-face-3", 101) = {2};
//+
Physical Surface("2-face-4", 102) = {7};
//+
Physical Surface("2-face-5", 103) = {1};
//+
Physical Surface("2-face-6", 104) = {12};
//+
Physical Volume("3-face-1", 105) = {1};
