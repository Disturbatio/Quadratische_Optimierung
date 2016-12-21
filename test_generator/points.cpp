//
//  main.cpp
//  test_creator
//
//  Created by Stepan Kharin on 19.12.16.
//  Copyright © 2016 Stepan Kharin. All rights reserved.
//

#include <iomanip>
#include <cmath>
#include <fstream>
#include <iostream>
using namespace std;

//defining constants for coefficients values range
const int LOWC=-10;
const int HIGHC=10;
const int RANGEC=HIGHC-LOWC;

//defining constants for x values range
const int LOWX=-10;
const int HIGHX=10;
const int RANGEX=HIGHX-LOWX;

//value of maximal difference between real and generated y value
//DIST=0 for exact y values
const double DIST=0.2;

//number of digits after decimal in generated variables
const int PRECISION=2;

int main(int argc, const char * argv[]) {
    if (argc==4)
    {
        string filename = argv[1];      //name of the file to create
        int n = atoi(argv[2]);          //polynomial degree
        int k = atoi(argv[3]);          //number of points to generate
        
        if (n<1){
            cout<<"Polynonial degree must be greater then 0"<<endl;
            return 1;
        }
        
        if (k<n+1){
            cout<<"Number of points to generate must be greater than degree of the polynomial"<<endl;
            return 1;
        }
        
        ofstream outfile;
        outfile.open(filename);
        outfile<<k<<endl;
        srand(time(NULL));
        
        double coeffs[n+1];
        for (int i=0;i<n+1;i++){
            coeffs[i]=(LOWC+(rand()%(RANGEC*100))/100.)/pow(i+1,2);  //randomizing polynomial coefficients within the given range with 2 digits after comma precision; the greater the power the less the value
        }
        
        double step = (double)RANGEX/k; //calclulating step for randomizing x values, so they are more or less evenly DISTributed within the given range
        for (int i=0;i<k;i++){
            double x=LOWX+step*i+(rand()%(int)(step*100))/100.; //calculating x value somewhere within step range
            double y=0;
            double fakex=x-step*DIST+step*(rand()%(int)(DIST*200))/100.; //calculating fake x somewhere within step*distortion range
            for (int j=0;j<n+1;j++){
                y+=coeffs[j]*pow(fakex,j);   //calculating y value
            }
        outfile<<fixed<<setprecision(PRECISION)<<x<<" "<<y<<endl;
        }
        outfile.close();
    }
    else
        cout<<"Wrong number of arguments!"<<endl;
    return 0;
}
