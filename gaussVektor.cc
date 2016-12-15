#include <iostream>
#include <stdlib.h>
#include <cmath>
#include<vector>

using namespace std;

double x[3];

/*

Funktion zum vereinen von 2 Matrizen. erwArr wird erstellt um WErte zu empfangen[wird so initialisiert weils ein 2D Pointer is]
Werte aus erstem Array eingefügt und Werte aus zweitem Array an letzte Stelle jeder Zeile gesetzt.

ACHTUNG: Muss noch für Vektoren modifiziert werden

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

void vectorPrint(vector<vector<double> >& matrix){
    cout << endl;
    for(int i = 0; i<matrix.size(); i++){
        cout << endl;
        for(int j = 0; j <matrix[0].size() ; j++){
            cout << " " << matrix[i][j];
        }
    }
    cout << endl;
}

/*

Funktion für Gauss Verfahren.

*/

double *gauss(vector<vector<double> >& matrix){
    
    cout << matrix.size() << "x" << matrix[0].size();
    /*
    Dreiecksmatrix für letzen schritt wird erstellt
    Wird Reihe für Reihe abgearbeitet
    */
    vectorPrint(matrix);
    cout << "Creating Triangle Matrix: " << endl;
    
    /*
        Nach größtem Glied in Spalte suchen, Element und Spalte speichern
    */
    for(int i=0; i < matrix.size(); i++){           
        double max = matrix[i][i];     //[i][i]
        int row = i;
        vectorPrint(matrix);
        
        
        
        for (int k=i+1; k<matrix.size(); k++) {    
            if (abs(matrix[k][i]) > max) {
                max = abs(matrix[k][i]);
                row = k;
            }
        }
        
        /*
        Schreibt Reihe die größtes Element beinhaltet in erste Reihe
        */

        for(int k = i; k<matrix[0].size(); k++){               
            double swap = matrix[row][k];
            matrix[row][k] = matrix[i][k];
            matrix[i][k] = swap;
        }
        vectorPrint(matrix);
        
        /*
        Löscht Elemente unter dem aktuellen größtem Element.
        Modifier wird erstellt durch aktuelles Element und aktuelles max, entspricht mathematischer Operation
        */
        
        for(int k=i+1; k < matrix.size(); k++){
            double modifier = -matrix[k][i]/matrix[i][i];
            for(int j = i; j <matrix[0].size(); j++){
                if (i==j){                      //wenn selber index in zeile = 0
                    matrix[k][j] = 0;
                }
                else{                           //sonst modifizieren
                    matrix[k][j] = matrix[k][j] + modifier*matrix[i][j];
                }
            }
        }
        
    
    }
    /*for(int i = 0; i<matrix.size(); i++){
        cout << endl;
        for(int j = 0; j <4 ; j++){
            cout << " " << matrix[i][j];
        }
    }
    cout << endl;*/
    vectorPrint(matrix);
    
    /*
    löst die erstellte Dreiecksmatrix.
    Fängt unten an, modifiziert die werte darüber.
    */
                                               
    for(int i = matrix.size() - 1; i>=0; i--){
        x[i] = matrix[i][matrix.size()]/matrix[i][i];
        for(int j=i-1;j>=0; j--){
            matrix[j][matrix.size()] = matrix[j][matrix.size()] - matrix[j][i] * x[i];
        }
    }
    return x;
    
}



int main(){
    double hardArrA[3][4]= {{9.0,25.5,92.2,15.5},{25.5,92.2,370.1,42.5},{92.2,370.1,1567.3,134.5}};
    double hardArrB[2][3]= {{1.0,2.0,3.0},{2.0,1.0,3.0}};

    
    int sizeA = 3;
    vector<double> vecColA(sizeA+1);
    vector<vector<double> > arrA(sizeA,vecColA);
    for(int i = 0; i < arrA.size(); i++){
        for(int j = 0; j < arrA[i].size();j++){ 
            (arrA[i])[j] = hardArrA[i][j];
            cout << arrA[i][j] << " ";
        }
        cout << endl;
    }
    double *x = gauss(arrA);
    cout << endl <<"ergebnis: " << endl << "x0: " <<x[0]<<endl<< "x1: " << x[1]<<endl<< "x2: " <<x[2]<< endl << endl;
    
    int sizeB = 2;
    vector<double> vecColB(sizeB+1);
    vector<vector<double> > arrB(sizeB,vecColB);
    for(int i = 0; i < arrB.size(); i++){
        for(int j = 0; j < arrB[i].size();j++){ 
            (arrB[i])[j] = hardArrB[i][j];
            cout << arrB[i][j] << " ";
        }
        cout << endl;
    }
    
    double *y = gauss(arrB);
    cout << endl <<"ergebnis: " << endl << "x0: " <<x[0]<<endl<< "x1: " << x[1]<<endl;//<< "x2: " <<x[2]<< endl << endl;
    
}   