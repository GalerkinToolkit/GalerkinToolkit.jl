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
// Physical Groups
//----------

Physical Point("0-face-1", 98) = {16};
//+
Physical Point("0-face-2", 99) = {15};
//+
Physical Point("0-face-3", 100) = {23};
//+
Physical Point("0-face-4", 101) = {1};
//+
Physical Point("0-face-5", 102) = {14};
//+
Physical Point("0-face-6", 103) = {12};
//+
Physical Point("0-face-7", 104) = {11};
//+
Physical Point("0-face-8", 105) = {36};
//+
Physical Point("0-face-9", 106) = {31};
//+
Physical Point("0-face-10", 107) = {32};
//+
Physical Point("0-face-11", 108) = {39};
//+
Physical Point("0-face-12", 109) = {40};
//+
Physical Point("0-face-13", 110) = {2};
//+
Physical Point("0-face-14", 111) = {3};
//+
Physical Point("0-face-15", 112) = {4};
//+
Physical Point("0-face-16", 113) = {5};
//+
Physical Point("0-face-17", 114) = {7};
//+
Physical Point("0-face-18", 115) = {8};
//+
Physical Point("0-face-19", 116) = {9};
//+
Physical Point("0-face-20", 117) = {10};
//+
Physical Point("0-face-21", 118) = {22};
//+
Physical Point("0-face-22", 119) = {20};
//+
Physical Point("0-face-23", 120) = {19};
//+
Physical Point("0-face-24", 121) = {30};
//+
Physical Point("0-face-25", 122) = {29};
//+
Physical Point("0-face-26", 123) = {28};
//+
Physical Point("0-face-27", 124) = {27};
//+
Physical Point("0-face-28", 125) = {26};
//+
Physical Point("0-face-29", 126) = {18};
//+
Physical Point("0-face-30", 127) = {17};
//+
Physical Point("0-face-31", 128) = {33};
//+
Physical Point("0-face-32", 129) = {34};
//+
Physical Curve("1-face-1", 79) = {22, 19};
//+
Physical Curve("1-face-2", 80) = {36, 34};
//+
Physical Curve("1-face-3", 81) = {33, 31};
//+
Physical Curve("1-face-4", 82) = {7, 4};
//+
Physical Curve("1-face-5", 83) = {18, 16};
//+
Physical Curve("1-face-6", 84) = {58, 56};
//+
Physical Curve("1-face-7", 85) = {42, 44};
//+
Physical Curve("1-face-8", 86) = {15, 24};
//+
Physical Curve("1-face-9", 87) = {14, 11};
//+
Physical Curve("1-face-10", 88) = {48, 41};
//+
Physical Curve("1-face-11", 89) = {1, 3};
//+
Physical Curve("1-face-12", 90) = {8, 10};
//+
Physical Surface("2-face-1", 92) = {4};
//+
Physical Surface("2-face-2", 93) = {12};
//+
Physical Surface("2-face-3", 94) = {9};
//+
Physical Surface("2-face-4", 95) = {1};
//+
Physical Surface("2-face-5", 96) = {7};
//+
Physical Surface("2-face-6", 97) = {2};
//+
Physical Volume("3-face-1", 91) = {1};

//----------
// Periodicity
//----------
// Back face copy of front face
Periodic Surface {9} = {4} Translate {0, 0, -1};
// Top face copy of bottom face
Periodic Surface {7} = {2} Translate {0, 1, 0};
// Right face copy of left face
Periodic Surface {12} = {1} Translate {1, 0, 0};
