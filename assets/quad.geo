//+
SetFactory("OpenCASCADE");
Rectangle(1) = {0, 0, 0, 1, 1, 0};
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
Physical Curve("1-face-2", 10) = {2};
//+
Physical Curve("1-face-3", 11) = {3};
//+
Physical Curve("1-face-4", 12) = {4};
//+
Physical Surface("2-face-1", 13) = {1};
