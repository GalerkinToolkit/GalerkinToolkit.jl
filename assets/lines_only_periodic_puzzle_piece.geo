// Only straight lines are periodics for this puzzle piece unit cell
SetFactory("OpenCASCADE");                                                                  
Rectangle(1) = {0, 0, 0, 1, 1, 0};                                                          
//+                                                                                         
Disk(2) = {1, 0.5, 0, 0.25, 0.25};                                                          
//+                                                                                         
Disk(3) = {-0, 0.5, 0, 0.25, 0.25};                                                         
//+                                                                                         
Disk(4) = {0.5, 1, 0, 0.25, 0.25};                                                          
//+                                                                                         
Disk(5) = {0.5, -0, 0, 0.25, 0.25};                                                         
//+                                                                                         
BooleanDifference{ Surface{1}; Delete; }{ Surface{3}; Surface{4}; Surface{2}; Surface{5}; Delete; }

// Right-side lines are periodic copy of left-side lines                                 
Periodic Curve {10} = {6} Translate {1, 0, 0};
Periodic Curve {12} = {3} Translate {1, 0, 0}; 

// Top lines are periodic copy of bottom line
Periodic Curve {9} = {13} Translate {0, 1, 0};
Periodic Curve {7} = {2} Translate {0, 1, 0}; 
                                   
//+
Physical Curve("1-face-1", 14) = {2, 13};
//+
Physical Curve("1-face-2", 15) = {12, 10};
//+
Physical Curve("1-face-3", 16) = {9, 7};
//+
Physical Curve("1-face-4", 17) = {6, 3};
//+
Physical Curve("1-face-4", 17) += {6, 3};
//+
Physical Point("0-face-1", 18) = {3};
//+
Physical Point("0-face-2", 19) = {2};
//+
Physical Point("0-face-3", 20) = {1};
//+
Physical Point("0-face-4", 21) = {13};
//+
Physical Point("0-face-5", 22) = {12};
//+
Physical Point("0-face-6", 23) = {11};
//+
Physical Point("0-face-7", 24) = {10};
//+
Physical Point("0-face-8", 25) = {9};
//+
Physical Point("0-face-9", 26) = {8};
//+
Physical Point("0-face-10", 27) = {7};
//+
Physical Point("0-face-11", 28) = {6};
//+
Physical Point("0-face-12", 29) = {4};
//+
Physical Surface("2-face-1", 30) = {1};