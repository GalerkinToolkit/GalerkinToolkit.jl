SetFactory("OpenCASCADE");
Box(1) = {0, 0, 0, 1, 1, 1};
Cylinder(2) = {0, 0.5, 0.5, 1, 0, 0, 0.4, 2*Pi};
Cylinder(3) = {0.5, 0.0, 0.5, 0, 1, 0, 0.4, 2*Pi};
BooleanUnion{ Volume{3}; Delete; }{ Volume{2}; Delete; }
BooleanDifference{ Volume{1}; Delete; }{ Volume{2}; Delete; }
Physical Surface("top") = {14};
Physical Surface("bottom") = {15};
Physical Volume("volume") = {1};
