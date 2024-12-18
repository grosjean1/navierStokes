SetFactory("OpenCASCADE");

DefineConstant[
	L={5,Min 3, Max 10, Step 1,Name "00Longueur",Visible 1},
	l={1,Min 1, Max 2, Step 0.5,Name "00Largeur",Visible 1}
];
h=0.15;
Mesh. CharacteristicLengthMin = h;
Mesh. CharacteristicLengthMax = h;

Point(1) = {0, 1/2, 0,h};
Point(2) = {l, 1/2, 0,h};
Point(3) = {l, 0, 0,h};
Point(4) = {L, 0, 0,h};
Point(5) = {L, 1, 0,h};
Point(6) = {0, 1, 0,h};


Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,1};

Line Loop(1) = {1,2,3,4,5,6};
Plane Surface (1) = {1};
Physical Surface (0)={1};
Physical Line(10)={6};
Physical Line(20)={1,2,3};
Physical Line(30)={4};
Physical Line(40)={5};
