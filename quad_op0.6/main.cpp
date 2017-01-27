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
                    bool servantfailed = true;
                    while(servantfailed){

                        cout << "\nMaster, your humbly servant Polly is here to help you generating some numbers.";

                        string filename;
                        cout << "\nPolly needs more informations - how do we call the file, master? ";
                        cin >> filename;

                        cout << "\nMaster, Polly needs to know the degree of the polynomial between 2 and 9: ";
                        int n=0, k=0;
                        bool qtest = true;
                        while(qtest){
                            cin >> n;
                            if(n>10||n<0){
                                cout << "\nMaster, I tried but I need a valid degree between 2 and 9" << endl;
                            } else{
                                qtest = false;
                            }
                        }

                        cout << "\nWe are near the end, Master but I have to know the number of points which I should generate ... Master: ";

                        while(qtest){
                            cin >> k;
                            if(k<n&&k<0){
                                cout << "\nMaster, I tried but I need a valid number of points. At least more than " <<n << endl;
                            } else{
                                qtest = false;
                            }
                        }

                        if(test_gen(filename, k, n)){
                            cout << "\nMaster, we did it. I present you the glorious numbers:\n";
                            optimize(filename, n);
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
                    cout << "\nEnter filename: ";
                    cin >> filename;
                    optimize(filename);
                    cout << "Press Enter to continue"<< endl;
                    cin.ignore();
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
