#ifndef _quad_op_h_
#define _quad_op_h_

#include <string>
#include <fstream>
#include <iomanip> //setprecision()
#include "const_values.h"
#include <ctime>
#include <chrono>

void add_log(std::string, bool);
double** makeA(int,int,double*);
double* makeB(int,int,double*,double*);
double product(int, double*, double*);
double* product(int, double**, double**);
double* product(int, double, double*);
double* subtract(int, double*, double*);
double* add(int, double*, double*);
double* linsolve(int, double**, double*);
double* polyfit(int, int, double*, double*);
double* optimize(std::string, int);
void plot(std::string);
bool test_gen(std::string, int, int);
void processgnufile(std::string, bool);
bool test_gen_showfile(std::string);

void add_log(std::string message="", bool timestamp=true){

    std::chrono::system_clock::time_point time_point;
    time_point = std::chrono::system_clock::now();
    std::time_t now_c = std::chrono::system_clock::to_time_t(time_point);

    std::ifstream logfiletest("log.txt");
    std::ofstream log("log.txt", std::ios_base::app | std::ios_base::out);
    if(!logfiletest){
        if(log.is_open()){
            log << ctime(&now_c)<< "-logfile created\n" << message << std::endl;
        }
    }
    else{
        if(timestamp){
            log << ctime(&now_c) << " ";
        }
        log << message << std::endl;
    }
}


double** makeA(int n, int k, double* x){

    //calculating all required sums of powers of x
 //   double sum[2*n+1];
 //   add_log("called makeA");
    int nsize = 2*n+1;
    double* sum = new double[nsize];
    for (int i=0;i<2*n+1;i++){
        sum[i]=0.0;
        for (int j=0;j<k;j++){
            sum[i]+=pow(x[j],i);
        }
    }

    //filling up the matrix
    double** A = new double*[n+1];
    for(int i = 0; i < n+1; i++) {
        A[i] = new double[n+1];
        for(int j = 0; j < n+1; j++){
            A[i][j]=sum[i+j];
        }
    }
    return A;
}

double* makeB(int n, int k, double* x, double* y){
  //  add_log("called makeB");
    //calculating sums of y[i]x[i]^i
    double* B = new double[n+1];
    for (int i=0;i<n+1;i++){_quad_op_h_
        B[i]=0;
        for (int j=0;j<k;j++)
            B[i]+=y[j]*pow(x[j], i);
    }
    return B;
}

// product of two vectors of size n
double product(int n, double* a, double* b){
    double p=0;
    for (int i=0;i<n;i++)
        p+=a[i]*b[i];
    return p;
}

//product of n*n matrix and a vector of size n
double* product(int n, double** A, double* b){
    double* p = new double[n];
    for (int i=0;i<n;i++){
        p[i]=0;
        for (int j=0;j<n;j++)
            p[i]+=A[i][j]*b[j];
    }
    return p;
}

//product of scalar and a vector of size n
double* product(int n, double a, double* b){
    double* p = new double[n];
    for (int i=0;i<n;i++){
        p[i]=a*b[i];
    }
    return p;
}

//difference between two vectors of size n
double* subtract(int n, double* a, double* b){
    double* s=new double[n];
    for (int i=0;i<n;i++)
        s[i]=a[i]-b[i];
    return s;
}

//sum of two vectors of size n
double* add(int n, double* a, double* b){
    double* s=new double[n];
    for (int i=0;i<n;i++)
        s[i]=a[i]+b[i];
    return s;
}

//solving system of linear equations using Conjugate Gradient method
//n is the size of the n x n matrix
double* linsolve(int grad, double** A, double* b) {
    add_log("called linesolve");
    double* x = new double[grad];
    double* r;
    double* d;

    std::fill(x,x+grad,0.);
    r=b;
    d=r;

    double r0=product(grad, r, r);
    double r1;
    double* Ap;
    double alpha;

    for (int i=0;i<grad;i++){
        Ap=product(grad, A, d);
        alpha=r0/(product(grad, d, Ap));
        x=add(grad, x, product(grad, alpha, d));
        r=subtract(grad, r, product(grad, alpha, Ap));
        r1=product(grad, r, r);
        if (sqrt(r1)<1e-10)
            break;
        d=add(grad, r, product(grad, r1/r0, d));
        r0=r1;
    }
    return x;
}

//polyfit calculates the coefficients of the given points
//n is the grade of the approximated polynomial
//k is the lenght of x and y array
double* polyfit(int n, int k, double* x, double* y){
    add_log("called polyfit");
    double* coeffs = linsolve(n+1, makeA(n,k,x), makeB(n,k,x,y));

    //saving the polynom coefficients with given precision
    if (PRECISION){
        for (int i=0;i<n+1;i++){
            coeffs[i]=(int)(coeffs[i]*pow(10,PRECISION))/pow(10,PRECISION);
        }
    }
    return coeffs;
}

//optimize creates a resultfile from the result of polyfit
//n is the grade of the approximated polynomial

double* optimize(std::string filename, int grad=2){
    std::string add;
    add = "function optimize: filename: " + filename;
    add_log(add);
    if(grad<1||grad>10)
    {
        std::cout << "\nGrad of Polynomial is either lower 1 or higher 10. Correct it!";
        add_log("Grad was lower 1 or higher 10 Abort");
        return 0;
    }
    else{
        add= " grad=";
        add+= std::to_string(grad);
        add_log(add);

        std::string resultfilename;
        resultfilename = filename+"_result.txt";
        filename+=".txt";
        std::ifstream infile;

        infile.open(filename);
        int k=0;
        double* result;

        int resultsize = grad+1; //um die Laenge des coeffs Array fuer die Ausgabe zu haben
    //    std::cout << "n" << n << std::endl;

        if (infile.is_open()){
            add_log("optimize: successfully opended infile");
            infile >> k;

            double *x = new double[k];
            double *y = new double[k];
            for(int i=0; i<k; i++){
                infile >> x[i] >> y[i];
            }

            result = polyfit(grad, k, x, y);

            delete[] x;
            delete[] y;
            infile.close();

            std::ofstream resultfile;
            resultfile.open(resultfilename);

            if(resultfile.is_open()){
                add_log("optimize: successfully opened resultfile");
                resultfile << resultsize<<std::endl;
                for(int i=0; i<resultsize; i++){
                    resultfile << result[i]<<std::endl;
                    std::string logtemp;
                    logtemp = std::to_string(result[i]);
                    add_log(logtemp, false);
                }

                resultfile.close();
            }
            else{
                std::cout << "\nUnable to open resultfile" << std::endl;
                add_log("optimize: Unable to open resultfile "+resultfilename);
                return 0;
            }

            std::cout << "\nSuccesfully optmizied" << std::endl;
            add_log("optimize: Successfully optmizied");
            return result;
        }
        else{
            std::cout << "\nUnable to open datafile" << std::endl;
            add_log("optimize: Unable to open datafile named "+filename);
            return 0;
        }

    }


}

//plot() recieves the filename of the points, if a resultfile already exist, plot() creates a Gnuplot script which saves the points and the resulting polynomial in png file
//if plot() can't open the resultfile, it calls optimize() and create a resultfile on its own
void plot(std::string filenamepoints){

    std::ofstream gnufile;
    std::ifstream coeffs;
    int numberofcoeffs=0;
    double numberbridge=0.0;

    std::string gnufilename;
    std::string gnuprintname;
    std::string filenamecoeffs;
    std::string logtemp;

    filenamecoeffs = filenamepoints+"_result";
    gnuprintname = filenamepoints;
    gnufilename = filenamepoints+".gnu";

    coeffs.open(filenamecoeffs+".txt");
    gnufile.open(gnufilename);

//    std::cout <<"\nfilenamecoefss = " << filenamecoeffs << "\nfilenamepoints = " << filenamepoints << "\ngnufilename = " << gnufilename;

    if(gnufile.is_open()&&coeffs.is_open()){
        gnufile << "set autoscale x\nset autoscale y ";
        gnufile << "\nset terminal png";
        gnufile << "\nset output \""<<gnuprintname<<".png\"";
        gnufile << "\nf(x)=";

        logtemp="plot: coeffs for gnufile are: ";
        add_log(logtemp);
        logtemp="";

        coeffs >> numberofcoeffs;
        for(int i =0; i<numberofcoeffs; i++){
                gnufile << "+";
                coeffs >> numberbridge;
                logtemp+= std::to_string(numberbridge);
                logtemp+= " ";
                gnufile << numberbridge << "*x**" << i;
                }
        add_log(logtemp, false);
        gnufile <<"\nplot '"<< filenamepoints+".txt" << "' every::1, f(x)"; //every::1 um erste Zeile zu ueberspringen
        gnufile <<"\nreplot";
        add_log("plot: created GnuScriptfile named "+gnufilename);

        gnufile.close();
        coeffs.close();

        processgnufile(gnuprintname, true);
    }
    else if(gnufile.is_open()){
        std::cout << "\nCan't plot the given file!";
        logtemp="plot:Can't plot file with points named: ";
        logtemp+=filenamepoints;
        add_log(logtemp);

        const char* delfilename = gnufilename.c_str();
        remove(delfilename);
        gnufile.close();
        coeffs.close();
/*
//        std::cout << "\ndummyfile muss erstellt werden - filename" << filenamepoints;
        int dummyint=2;
        double* dummydouble;

        dummydouble = optimize(filenamepoints, dummyint);

        gnufile << "set autoscale x\nset autoscale y ";
        gnufile << "\nset terminal png";
        gnufile << "\nset output \""<<gnuprintname<<".png\"";
        gnufile << "\nf(x)=";

        numberofcoeffs = dummyint+1;
        logtemp="plot: coeffs for gnufile are: ";
        add_log(logtemp);
        for(int i =0; i<numberofcoeffs; i++){
                gnufile << "+";
                gnufile << dummydouble[i] << "*x**" << i;
                logtemp+= std::to_string(dummydouble[i]);
                logtemp+= " ";
                }
        add_log(logtemp, false);
        gnufile <<"\nplot '"<< filenamepoints+".txt" << "', f(x)";
        gnufile <<"\nreplot";
        add_log("plot: create GnuScriptfile with dummyfile named "+gnufilename);
*/
    }



}

//test_gen creates <filename>.txt with k random points following a polynomial with n terms
bool test_gen(std::string filename, int k, int grad){
        add_log("function: test_gen");

       if (grad<1){
            std::cout<<"Polymonial degree must be greater then 0"<<std::endl;
            add_log("Wrong Input - Polynomial degree <= 0");
            return false;
        }

        if (k<grad+1){
            std::cout<<"Number of points to generate must be greater than degree of the polynomial"<<std::endl;
            add_log("Wrong Input - Number of points were smaller than the polynomial degree");
            return false;
        }


        std::string logtext = "";

        logtext+= "test_gen Result von Gleichung vom Grad:";
        logtext+= std::to_string(grad);
        logtext+=" Anzahl Punkte:";
        logtext+= std::to_string(k);
        add_log(logtext);

        std::ofstream outfile;
        outfile.open(filename+".txt");
        outfile<<k<<std::endl;
        srand(time(NULL));

        int nsize = grad+1;
        double* coeffs = new double;
        double** pointvalues = new double*[k+1];

        logtext = "random coeffs: ";

        for (int i=0;i<nsize;i++){
            double temp = ((LOWC+(rand()%(RANGEC*100))/100.)/pow(i+1,2));
            coeffs[i]=temp;

            logtext+=std::to_string(coeffs[i]);
            logtext+=" ";
        }
        logtext+="\nVolitional error is between +- ";
        logtext+=std::to_string(DIST);
        logtext+="\nRandom Values: ";
        add_log(logtext);

        for (int i=0;i<k;i++){

            double x= (double)rand()/RAND_MAX;
            x = LOWX + x *(HIGHX-LOWX);
            double xvolitionalerror = (double)rand()/RAND_MAX;
            xvolitionalerror = -DIST + xvolitionalerror*(DIST+DIST);//volitional error for more random looking points
            x+=xvolitionalerror;

            pointvalues[i] = new double[2];
            pointvalues[i][0] = x;
            add_log("x: "+std::to_string(x), false);
            double y=0;
            double yvolitionalerror = (double)rand()/RAND_MAX;
            yvolitionalerror = -DIST + yvolitionalerror*(DIST+DIST);

            for (int j=0;j<nsize;j++){

                y+=coeffs[j]*pow(x,j);
            }
            y+=yvolitionalerror;
            pointvalues[i][1] = y;
            add_log("y: "+std::to_string(y), false);
            outfile<<std::fixed<<std::setprecision(PRECISION)<<x<<" "<<y<<std::endl;
        }

        outfile.close();
        return true;

}

//executes a gnu filescript created f.ex. with the plot function
void processgnufile(std::string filename, bool showpicture=true){
    std::string logstring;
    logstring ="processgnufile: ";
    logstring+=filename;
    logstring+=" showpic:";
    if(showpicture){
        logstring+="true";
    }else{
        logstring+"false";
    }
    add_log(logstring);

    std::string excommandtemp;
    excommandtemp = "gnuplot ";
    excommandtemp+=filename;
    excommandtemp+=".gnu";

    const char* excommand = excommandtemp.c_str();


    if(!system(excommand)){
            logstring = "processgnufile: successfully executed command was:";
            logstring+=excommand;
            add_log(logstring);
            std::cout << "\nGnuFile processed!";
            if(showpicture){
                excommandtemp = "xdg-open ";
                excommandtemp+= filename;
                excommandtemp+=".png";

                const char* excommandprint = excommandtemp.c_str();

                if(!system(excommandprint)){
                        logstring = "processgnufile: successfully executed command was:";
                        logstring+=excommandprint;
                        add_log(logstring);
                        std::cout << "\nPicture displayed!";
                }
                else{
                    logstring = "processgnufile: failed to execute command:";
                    logstring+=excommandprint;
                    add_log(logstring);
                }
            }
    }else{
        logstring = "processgnufile: failed to execute command:";
        logstring+=excommand;
        add_log(logstring);
    }


}
//shows the content of the files
bool test_gen_showfile(std::string filename){
    std::ifstream infile;
    infile.open(filename+".txt");
    std::string line;
    if (infile.is_open()){
        while (getline(infile, line)){
            std::cout << line << std::endl;
        }
        infile.close();
    return true;
    }
    else{
        std::cout << "\nUnable to open file" << std::endl;
        add_log("test_gen_showfile : Unable to open file");
        return false;
    }
}
#endif

