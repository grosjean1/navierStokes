# dans /gmsh:
maillages .geo

# dans /maillages:
carre.msh (carre 10*10); test.msh (rectangle 1*2) ; test_triangle.msh (triangle de ref)

# carre.edp: test sur le carre 10*10
# projet.edp: comparaison avec Freefem

# Autres fichiers:
mainNS.cpp
Fonctions_Utiles.hpp
matNS.hpp
mesh.cpp
mesh.hpp
R2.hpp
plot.edp

marche.msh ou projet.msh (differentes tailles de maillages)
Makefile 
