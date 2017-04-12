/* Debut-Commentaires

Nom de la fonction: H6

Entrees: n entier, longueur du vecteur x ; x, un vecteur de n nombres de type double

Sorties: void

Description:
Calcul de H6(x), H6 etant le sixieme polynome de Legendre 
Modifie le vecteur x dans le programme principal, en le remplacant par H6(x)

Utilisation dans une fonction main:
-------------------------------------
#include <iostream.h>
#include "H6.cc"

int main()

{
int i,n;
double *x,*y,*xavant;

 n=5;   // longueur du vecteur x

x=new double[n]; // initialisation de x
y=new double [n];     // initialisation de y
xavant=new double [n];     // initialisation de xavant

x[0]=1;x[1]=2;x[2]=3;x[3]=4;x[4]=5;  // affectation de valeurs a x

for (i=0;i<n;i++){xavant[i]=x[i];} // On sauvegarde les valeurs de x dans xavant

 H6(n,x); // calcul de H6(x), remplace x par H6(x)

 for (i=0;i<n;i++){y[i]=x[i];}      // affecte a y le resultat
 for (i=0;i<n;i++){x[i]=xavant[i];}  // remet les bonnes valeurs dans x

 cout << "Affichage des valeurs de y\n";

for (i=0;i<n;i++)
cout << y[i] << " ";  // Affichage des valeurs de y=H6(x)

 cout << "\n";

 cout << "Affichage des valeurs de x\n";

for (i=0;i<n;i++)
cout << x[i] << " ";  // Affichage des valeurs de x

 cout << "\n";

delete [] xavant;  // On libere de la memoire

return(0);

}

-------------------------------------

Instructions de compilation:
g++ -Wall -O nom_du_fichier_contenant_la_fonction_main.cc

Fonctions exterieures appelees:


Auteur: Pierre Lafaye de Micheaux

Date: 04/01/2001

Fin-Commentaires */

#include<math.h>

void H6(int n, double *x) {
  int i;
  for (i = 0; i < n; i++) {
    x[i] = sqrt(13.0) * (231.0 * pow(x[i], 6.0) - 315.0 * pow(x[i], 4.0) + 105.0 * pow(x[i], 2.0) - 5.0) / 16.0;
  }
  return;
}

