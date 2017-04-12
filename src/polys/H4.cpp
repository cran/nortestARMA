/* Debut-Commentaires

Nom de la fonction: H4

Entrees: n entier, longueur du vecteur x ; x, un vecteur de n nombres de type double

Sorties: void

Description:
Calcul de H4(x), H4 etant le quatrieme polynome de Legendre 
Modifie le vecteur x dans le programme principal, en le remplacant par H4(x)

Utilisation dans une fonction main:
-------------------------------------
#include <iostream.h>
#include "H4.cc"

int main()

{
int i,n;
double *x,*y,*xavant;

 n=5;   // longueur du vecteur x

 x=new double[n]; // initialisation de x

 x[0]=1;x[1]=2;x[2]=3;x[3]=4;x[4]=5;  // affectation de valeurs a x

 xavant=x;  // On sauvegarde les valeurs de x

 H4(n,x); // calcul de H4(x)

 y=x;  // y contient le vecteur H4(x)

 x=xavant;  // On replace dans x ses valeurs

 for (i=0;i<n;i++){cout << y[i] << " ";}  // On affiche les valeurs du vecteur y=H4(x)

 delete [] xavant;  // On libere la memoire

return(0);

}

-------------------------------------

Auteur: Pierre Lafaye de Micheaux

Date: 04/01/2001

Fin-Commentaires */

#include<math.h>

void H4(int n, double *x) {
  int i;
  for (i = 0; i < n; i++) {
    x[i] = 3.0 * (35.0 * pow(x[i], 4.0) - 30.0 * pow(x[i], 2.0) + 3.0) / 8.0;
  }
  return;
}

