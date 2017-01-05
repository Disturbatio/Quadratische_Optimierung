
//
//  quadopt.cpp
//  SOPT
//
//  Created by Stepan Kharin on 19.12.16.
//  Copyright © 2016 Stepan Kharin. All rights reserved.
//

#include <cmath>
#include <iostream>
#include <string>
#include "quad_op.h"



using namespace std;


int main(int argc, const char * argv[]) {


bool exit = false;
int selected = 1;
    if (argc==1){
        while(!exit){
            cout << "\nQuadratische Optimierung!\n=========================" << endl;
            cout << "\nChoose wisely: " << endl;
            cout << "\n(1)generate testfiles"; //ToDo eventuell gleich an optmierer uebergeben
            cout << "\n(2)show textfile";
            cout << "\n(3)load file and optimize!";//ToDo eventuell gleich an Gnuplot uebergeben und Ergebniss in File speicher?
            cout << "\n(4)GnuPlot";
            cout << "\n(0)Exit";

            cout << "\n\nYour turn: ";
            cin >> selected;

            switch(selected){
                case 0:
                    exit = true;
                    break;
                case 1:{
                    bool servantfailed = true;
                    while(servantfailed){

                        cout << "\nMaster, your humbly servant Polly is here to help you generating some numbers.";

                        string filename;
                        cout << "\nPolly need more informations - how do we call the file, master? ";
                        cin >> filename;
                        filename+=".txt";

                        int n=0, k=0;
                        cout << "\nMaster, Polly need to know the degree of the polynomial: ";
                        cin >> n;
                        cout << "\nWe are near the end, Master but I have to know the number of points which I should generate ... Master: ";
                        cin >> k;

                        if(test_gen(filename, k, n)){
                            cout << "\nMaster, we did it. I present you the glorious numbers:\n";
                            servantfailed = false;
                            if(test_gen_showfile(filename)){
                                cout << "\n\nMaster, it was a pleasure to serve you! Should I print your godlike numbers? (yes/no)\n";
/*
                                string servantanswer;
                                cin >> servantanswer;
                                if(servantanswer=="no")
                                    selected = 0;
                                else if(servantanswer=="yes")
                                    plot(filename);
*/
                            }
                            else
                                cout << "\nMaster, I have failed with the numbers, please forgive me.";
                        }
                        else{
                            cout << "\nMaster, your servant failed! Please, tell me how to do it right.";
                            servantfailed = true;
                        }
                    }

                    break;
                }
                case 2:{
                    string filename;
                    cout << "\nEnter Filename: ";
                    cin >> filename;
                    filename+=".txt";
                    if(test_gen_showfile(filename)){
                        cout << "Press Enter to continue"<< endl;
                        cin.ignore();
                        getchar();
                    }
                    else
                        cout << "\nOoops, something went wrong";
                    break;

                }
                case 3:{
                    string filename;
                    cout << "\nEnter filename: ";
                    cin >> filename;
                    double *result;
                    int resultsize=0;


                    result = optimize(filename, resultsize);
/*
//testdaten aus den Folien

                    double x[9]={0.5, 1, 1.5, 2.5, 3, 3.5, 4, 4.5, 5};
                    double y[9]={0.5, 1.5, 2, 3, 3, 2.5, 2, 1, 0};
                    int n=2;
                    int k=9;
                    resultsize = n;

                    result = polyfit(n,k,x,y);
*/
                    if(resultsize)
                    {
                        for (int i=0;i<resultsize;i++)
                            cout<<result[i]<<endl;
                    }

                    cout << "Press Enter to continue"<< endl;
                    cin.ignore();
                    break;
                }
                case 4:{
                //fuer Gnuplot muessen die Daten bearbeitet werden, vor allem die erste Zeile muss entfernt werden
                //mit system Funktion gnu script erstellen und dann aufrufen?
                    string filenamepoints;
                    string filenamecoeffs;
                    cout << "\nEnter filename for points: ";
                    cin >> filenamepoints;
//                    cout << "\nEnter filename for calculated coefficients: ";
//                    cin >> filenamecoeffs;


//                    plot(filenamepoints, filenamecoeffs);
                    plot(filenamepoints);

                    break;
                }
                default:
                    cout << "\nPlease enter a valid choice";
                    break;

            }
        }
    }
/*
    else if (argc==4){
        string filename = argv[1];      //name of the file to create
        int n = atoi(argv[2]);          //polynomial degree
        int k = atoi(argv[3]);          //number of points to generate

  //    call with this scheme double x[9]= {....}; double y[9]={....}; int n=..; int k=..;
  //    double* coeffs = polyfit(n,k,x,y);


        if(test_gen(filename, n, k))
           cout << "test_gen worked";
        else
            cout << "test_gen failed due to your incompentece";
    } else
        cout<<"Wrong number of arguments!"<<endl;



    //test data

    double x[9]={0.5, 1, 1.5, 2.5, 3, 3.5, 4, 4.5, 5};
    double y[9]={0.5, 1.5, 2, 3, 3, 2.5, 2, 1, 0};
    int n=2;
    int k=9;

    double* coeffs = polyfit(n, k, x, y);
    for (int i=0;i<n+1;i++)
        cout<<coeffs[i]<<endl;
*/

    return 0;
}
