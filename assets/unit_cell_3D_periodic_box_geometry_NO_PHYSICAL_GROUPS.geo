// TODO: Physical group naming convention... check glk output for 3D?

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


