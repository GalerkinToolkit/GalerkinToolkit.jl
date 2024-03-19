//+
SetFactory("OpenCASCADE");
Rectangle(1) = {0, 0, 0, 1, 1, 0};
// Right edge is periodic copy of left edge
Periodic Curve {2} = {4} Translate {1, 0, 0};
// Top edge is periodic copy of bottom edge
Periodic Curve {3} = {1} Translate {0, 1, 0};

//+
Physical Point("0-face-1", 5) = {1};
//+
Physical Point("0-face-2", 6) = {2};
//+
Physical Point("0-face-3", 7) = {4};
//+
Physical Point("0-face-4", 8) = {3};
//+
Physical Curve("1-face-1", 9) = {1};
//+
Physical Curve("1-face-2", 10) = {4};
//+
Physical Curve("1-face-3", 11) = {2};
//+
Physical Curve("1-face-4", 12) = {3};
//+
Physical Surface("2-face-1", 13) = {1};
