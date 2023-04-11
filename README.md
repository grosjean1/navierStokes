Implémentation en C++ de la méthode des Elements Finis avec Fe space Taylor-Hood(P2-P1)
Navier-Stokes 
Backward-Facing step 2D, 



# dans /gmsh:
maillages .geo

# dans /maillages:
carre.msh (carre 10*10); test.msh (rectangle 1*2) ; test_triangle.msh (triangle de ref)

# carre.edp: test sur le carre 10*10
# projet.edp: comparaison avec Freefem

# C++ fichiers:
mainNS.cpp
Fonctions_Utiles.hpp
matNS.hpp
mesh.cpp
mesh.hpp
R2.hpp
plot.edp

Autres maillages
marche.msh ou projet.msh (differentes tailles de maillages)


Makefile 
