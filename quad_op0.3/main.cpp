
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




    if (argc==1){
        add_log("Program started in menu");
        bool exit = false;
        int selected = 1;
        while(!exit){
            cout << "\nQuadratische Optimierung!\nB=========================D" << endl;
            cout << "\nChoose wisely: " << endl;
            cout << "\n(1)generate testfiles";
            cout << "\n(2)show textfile";
            cout << "\n(3)load file and optimize!";
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

                                string servantanswer;
                                cin >> servantanswer;
                                if(servantanswer=="no")
                                    selected = 0;
                                else if(servantanswer=="yes"){
                                    plot(filename);
                                }
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
                    optimize(filename);
                    cout << "Press Enter to continue"<< endl;
                    cin.ignore();
                    break;
                }
                case 4:{
                    string filenamepoints;
                    string filenamecoeffs;
                    cout << "\nEnter filename for points: ";
                    cin >> filenamepoints;

                    plot(filenamepoints);

                    break;
                }
                default:
                    cout << "\nPlease enter a valid choice";
                    break;

            }
        }
    }

    else if (argc==3){
    /*
        string filename = argv[1];      //name of the file to create
        int n = atoi(argv[2]);          //polynomial degree
        int k = atoi(argv[3]);          //number of points to generate

  //    call with this scheme double x[9]= {....}; double y[9]={....}; int n=..; int k=..;
  //    double* coeffs = polyfit(n,k,x,y);

*/
        string filename = argv[1]; //name of the file with points
        int n = atoi(argv[2]); //polynomial degree
        add_log("Program started via console with file "+filename);
        optimize(filename, n);
        //cout << "\nDone! Look for the file:"+filename+".gnu" << endl;
        cout << "\nDone!" << endl;
    }
    else{
        cout << "Wrong number of parameters - correct way quad_op <filename> <grade_of_polynomial>" << endl;
        add_log("Tried to start with wrong number of parameters via console");
    }
    return 0;
}
