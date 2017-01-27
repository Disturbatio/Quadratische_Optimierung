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
            int sel_number=1;
            cout << "\nQuadratische Optimierung!\n=========================" << endl;
            cout << "\nChoose wisely: " << endl;
            cout << "\n("<< sel_number++ <<")generate testfiles";
            cout << "\n("<< sel_number++ <<")load file and optimize!";
            cout << "\n("<< sel_number++ <<")GnuPlot";
            cout << "\n("<< sel_number++ <<")create 'btest' files with random Testgen points, result and GnuPlot Script";

            cout << "\n(0)Exit";

            cout << "\n\nYour turn: ";
            cin >> selected;

            switch(selected){
                case 0:
                    exit = true;
                    break;
                case 1:{
                    bool worked = false;
                    worked = randomGenerator();
                    if(worked){
                        std::cout << "\nGenerator worked" << endl;
                    }
                    break;
                }
                case 2:{
                    string filename;
                    cout << "\nEnter filename: ";
                    cin >> filename;
                    cout << "\nEnter degree: ";
                    int degree=0;
                    cin >> degree;
                    optimize(filename, degree);
                    break;
                }
                case 3:{
                    string filenamepoints;
                    string filenamecoeffs;
                    cout << "\nEnter filename for points: ";
                    cin >> filenamepoints;

                    plot(filenamepoints);

                    break;
                }
                case 4:{
                    bool test = false;
                    string btestfilename = "btest";
                    int grad = 0;
                    grad = 2;
                    test = test_gen(btestfilename, 54, grad);
                    if(test){}
                    optimize(btestfilename);
                    plot(btestfilename);

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
        cout << "\nDone!" << endl;
    }
    else{
        cout << "Wrong number of parameters - correct way quad_op <filename> <grade_of_polynomial>" << endl;
        add_log("Tried to start with wrong number of parameters via console");
    }
    return 0;
}
