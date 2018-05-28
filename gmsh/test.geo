SetFactory("OpenCASCADE");

DefineConstant[
	L={5,Min 3, Max 10, Step 1,Name "00Longueur",Visible 0},
	l={1,Min 1, Max 2, Step 0.5,Name "00Largeur",Visible 0}
];
h=0.2;
Mesh. CharacteristicLengthMin = h;
Mesh. CharacteristicLengthMax = h;

Rectangle(1)= {0,0.5,0,l,0.5};
Rectangle(2)= {l,0,0,L,1};
BooleanUnion{Surface{1};Delete;}{Surface{2};Delete;}

Physical Surface (0)={1,2};
Physical Line(10)={2};
Physical Line(20)={1,6,5};
Physical Line(30)={7};
Physical Line(40)={8};


