//+
SetFactory("OpenCASCADE");
Rectangle(1) = {0, 0, 0, 1, 1, 0};
// Right edge is periodic copy of left edge
Periodic Curve {2} = {4} Translate {1, 0, 0};
// Top edge is periodic copy of bottom edge
Periodic Curve {3} = {1} Translate {0, 1, 0};

//+
// Physical Point("0-face-1", 5) = {1};
// //+
// Physical Point("0-face-2", 6) = {2};
// //+
// Physical Point("0-face-3", 7) = {4};
// //+
// Physical Point("0-face-4", 8) = {3};
// //+
// Physical Curve("1-face-1", 9) = {1};
// //+
// Physical Curve("1-face-2", 11) = {2};
// //+
// Physical Curve("1-face-3", 10) = {3};
// //+
// Physical Curve("1-face-4", 12) = {4};
//+
Physical Surface("2-face-1", 13) = {1};
//+
Physical Point("0-face-1", 14) = {1};
//+
Physical Point("0-face-2", 15) = {2};
//+
Physical Point("0-face-3", 16) = {4};
//+
Physical Point("0-face-4", 17) = {3};
//+
Physical Curve("1-face-1", 18) = {1};
//+
Physical Curve("1-face-2", 19) = {3};
//+
Physical Curve("1-face-3", 20) = {4};
//+
Physical Curve("1-face-4", 21) = {2};
