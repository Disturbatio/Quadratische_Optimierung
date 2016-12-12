#include <iostream>
#include <stdlib.h>
# include <cmath>

using namespace std;

double x[3];

/*

Funktion zum vereinen von 2 Matrizen. erwArr wird erstellt um WErte zu empfangen[wird so initialisiert weils ein 2D Pointer is]
Werte aus erstem Array eingefügt und Werte aus zweitem Array an letzte Stelle jeder Zeile gesetzt.

*/

double **joinMatrix(double a[3][3], double b[3]){
    double **erwArr = new double*[2];
    
    for(int i=0; i < 3; i++){
        erwArr[i] = new double[3];
    }
    
    for(int i = 0; i < 3;i++){
        for(int j =0; j < 4;j++){
            if(j<3){
                erwArr[i][j] = a[i][j];
            }
            else{
                erwArr[i][j] = b[i];
            }
        }
    }
    return erwArr;
}

/*

Funktion für Gauss Verfahren.
[int n nicht wirklich notwendig hier, wird benötigt wenn Array von unbestiummter Größe übergeben wird]

*/

double *gauss(double a[3][4]){
    int n = sizeof a[0] / sizeof(double);
    double **prepArr = new double*[2];
    
    for(int i=0; i < n-1; i++){
        prepArr[i] = new double[3];
    }

    for(int i = 0; i < n-1;i++){
        for(int j =0; j < n;j++){
                prepArr[i][j] = a[i][j];
        }
    }
    /*
    Dreiecksmatrix für letzen schritt wird erstellt
    Wird Reihe für Reihe abgearbeitet
    */
    cout << "Creating Triangle Matrix: " << endl;
    for(int i=0; i < n-1; i++){           
        double max = prepArr[i][i];     //[i][i]
        int row = i;
        
        /*
        Nach größtem Glied in Spalte suchen, Element und Spalte speichern
        */
        
        for (int k=i+1; k<n-1; k++) {    
            if (abs(prepArr[k][i]) > max) {
                max = abs(prepArr[k][i]);
                row = k;
            }
        }
        
        /*
        Schreibt Reihe die größtes Element beinhaltet in erste Reihe
        */

        for(int k = i; k<n; k++){               
            double swap = prepArr[row][k];
            prepArr[row][k] = prepArr[i][k];
            prepArr[i][k] = swap;
        }
        
        /*
        Löscht Elemente unter dem aktuellen größtem Element.
        Modifier wird erstellt durch aktuelles Element und aktuelles max, entspricht mathematischer Operation
        */
        
        for(int k=i+1; k < n-1; k++){
            double modifier = -prepArr[k][i]/prepArr[i][i];
            for(int j = i; j < n; j++){
                if (i==j){                      //wenn selber index in zeile = 0
                    prepArr[k][j] = 0;
                }
                else{                           //sonst modifizieren
                    prepArr[k][j] = prepArr[k][j] + modifier*prepArr[i][j];
                }
            }
        }
        
    
    }
    for(int i = 0; i<n-1; i++){
        cout << endl;
        for(int j = 0; j <4 ; j++){
            cout << " " << prepArr[i][j];
        }
    }
    cout << endl;
    
    /*
    löst die erstellte Dreiecksmatrix.
    Fängt unten an, modifiziert die werte darüber.
    */
                                               
    for(int i = 2; i>=0; i--){
        x[i] = prepArr[i][3]/prepArr[i][i];
        for(int j=i-1;j>=0; j--){
            prepArr[j][3] = prepArr[j][3] - prepArr[j][i] * x[i];
        }
    }
    return x;
    
}



int main(){
    double arrA[3][3]= {{9.0,25.5,92.2},{25.5,92.2,370.1},{92.2,370.1,1567.3}};
    double arrB[3]= {15.5,42.5,134.5};
    
    double arrFin[3][4];
    double finArr[3][4];
    
    
    cout << "Joined Matrix[finArr]: " << endl;
    double **joinArr = joinMatrix(arrA,arrB);
    for(int i = 0; i < 3;i++){
            for(int j =0; j < 4;j++){
                
                finArr[i][j] = joinArr[i][j];
                cout << finArr[i][j] << " ";
        }
        cout << endl;
    }

    delete[] joinArr;
    double *x = gauss(finArr);

    
    cout << endl <<"ergebnis: " << endl << "x0: " <<x[0]<<endl<< "x1: " << x[1]<<endl<< "x2: " <<x[2]<< endl;


    
}   