Pour génerer le Makefile :

cmake CMakeLists.txt 
------------------------------
Pour compiler le code et l'exécuter : 

make
mpirun -n  6 ./CHP 

Modifiez 6 par le nombre de processeurs que vous souhaitez activer
--------------------------------------

Pour regrouper les résultats :

g++ -o run read.cpp Concatrenation.cpp
./run

"Puis entrez le nombre de procs activé après la demande (par exemple 6)"

--------------------------------------
Pour afficher les résultats dans gnuplot :

load "plot.txt"

--------------------------------------
Pour modifier les paramètres de calculs :

Les paramètres de calcul tels que cas, xmin, xmax, ymin, ymax, Tf, Nx, Ny, CFL, space_scheme et time_scheme 
sont lus à partir du fichier param.txt dans cet ordre. Vous pouvez modifier ces paramètres dans ce fichier selon vos besoins. 


On a chosis 1 pour le shéma décentré et 2 pour le shéma centré .
On a choisis 1 pour le shéma implicite et 2 pour le shéma explicite .


--------------------------------------
// le cas de 1 proc est similaire au cas séquentiel , vous pouvez le tester pour le shéma explicite et implicite
mais il prend beaucoup plus de temps que le code séquentiel.
