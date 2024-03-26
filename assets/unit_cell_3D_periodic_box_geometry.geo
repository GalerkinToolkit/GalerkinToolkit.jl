// Physical group based on cartesian 3D glk mesh

//+
SetFactory("OpenCASCADE");
Box(1) = {0, 0, 0, 1, 1, 1};

//----------
// Periodicity
//----------
// Back face copy of front face
Periodic Surface {5} = {6} Translate {0, 0, -1};
// Top face copy of bottom face
Periodic Surface {4} = {3} Translate {0, 1, 0};
// Right face copy of left face
Periodic Surface {2} = {1} Translate {1, 0, 0};


//+
Physical Point("0-face-1", 13) = {2};
//+
Physical Point("0-face-2", 14) = {6};
//+
Physical Point("0-face-3", 15) = {4};
//+
Physical Point("0-face-4", 16) = {8};
//+
Physical Point("0-face-5", 17) = {1};
//+
Physical Point("0-face-6", 18) = {5};
//+
Physical Point("0-face-7", 19) = {3};
//+
Physical Point("0-face-8", 20) = {7};
//+
Physical Curve("1-face-1", 21) = {9};
//+
Physical Curve("1-face-2", 22) = {11};
//+
Physical Curve("1-face-3", 23) = {4};
//+
Physical Curve("1-face-4", 24) = {8};
//+
Physical Curve("1-face-5", 25) = {10};
//+
Physical Curve("1-face-6", 26) = {12};
//+
Physical Curve("1-face-7", 27) = {2};
//+
Physical Curve("1-face-8", 28) = {6};
//+
Physical Curve("1-face-9", 29) = {1};
//+
Physical Curve("1-face-10", 30) = {3};
//+
Physical Curve("1-face-11", 31) = {5};
//+
Physical Curve("1-face-12", 32) = {7};
//+
Physical Surface("2-face-1", 33) = {5};
//+
Physical Surface("2-face-2", 34) = {6};
//+
Physical Surface("2-face-3", 35) = {3};
//+
Physical Surface("2-face-4", 36) = {4};
//+
Physical Surface("2-face-5", 37) = {1};
//+
Physical Surface("2-face-6", 38) = {2};
//+
Physical Volume("3-face-1", 39) = {1};
