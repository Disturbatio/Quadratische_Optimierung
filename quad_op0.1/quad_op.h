#include <string>
#include <fstream>
#include <iomanip> //setprecision()
#include "const_values.h"


double** makeA(int n, int k, double* x){

    //calculating all required sums of powers of x
 //   double sum[2*n+1];
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

    //calculating sums of y[i]x[i]^i
    double* B = new double[n+1];
    for (int i=0;i<n+1;i++){
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
double* linsolve(int n, double** A, double* b) {

    double* x = new double[n];
    double* r;
    double* d;

    std::fill(x,x+n,0.);
    r=b;
    d=r;

    double r0=product(n, r, r);
    double r1;
    double* Ap;
    double alpha;

    for (int i=0;i<n;i++){
        Ap=product(n, A, d);
        alpha=r0/(product(n, d, Ap));
        x=add(n, x, product(n, alpha, d));
        r=subtract(n, r, product(n, alpha, Ap));
        r1=product(n, r, r);
        if (sqrt(r1)<1e-10)
            break;
        d=add(n, r, product(n, r1/r0, d));
        r0=r1;
    }

    return x;
}

double* polyfit(int n, int k, double* x, double* y){
    double* coeffs = linsolve(n+1, makeA(n,k,x), makeB(n,k,x,y));

    //saving the polynom coefficients with given precision
    if (PRECISION){
        for (int i=0;i<n+1;i++){
            coeffs[i]=(int)(coeffs[i]*pow(10,PRECISION))/pow(10,PRECISION);
        }
    }

    return coeffs;
}

double* optimize(std::string filename, int &resultsize, int n=2){
    std::cout << "\nin optimize filename: " << filename << std::endl;
    std::string resultfilename;
    resultfilename = filename+"_result.txt";
    filename+=".txt";
    std::ifstream infile;


    infile.open(filename);
//    std::string line;
//    int n = 2, k=0; //es soll ja eine quadratsiche Optimierung sein daher der Fixwert
    int k=0;
    double* result;

    resultsize = n+1; //um die Laenge des coeffs Array fuer die Ausgabe zu haben
//    std::cout << "n" << n << std::endl;

    if (infile.is_open()){
        std::cout << "\nin optimize: successfully opended infile";
        infile >> k;

        double *x = new double[k];
        double *y = new double[k];
        for(int i=0; i<k; i++){
            infile >> x[i] >> y[i];
        }

        result = polyfit(n, k, x, y);

        delete[] x;
        delete[] y;
        infile.close();

        std::ofstream resultfile;
        resultfile.open(resultfilename);

        if(resultfile.is_open()){
            std::cout << "\nin optimize: successfully opened resultfile";
            resultfile << resultsize<<std::endl;
            for(int i=0; i<resultsize; i++)
                resultfile << result[i]<<std::endl;

            resultfile.close();
        }
        else{
            std::cout << "\nUnable to open resultfile";
            return 0;
        }

        std::cout << "\noptmize success" << std::endl;
        return result;
    }
    else{
        std::cout << "\nUnable to open datafile" << std::endl;
        return 0;
    }
}

//void plot(std::string filenamepoints, std::string filenamecoeffs){
void plot(std::string filenamepoints){

    std::ofstream gnufile;
    std::ifstream coeffs;
    int numberofcoeffs=0;
    double numberbridge=0.0;

    std::string gnufilename;
    std::string gnuprintname;
    std::string filenamecoeffs;

    filenamecoeffs = filenamepoints+"_result";
    gnuprintname = filenamepoints;
    gnufilename = filenamepoints+".gnu";

 //   filenamepoints+=".txt";
    //filenamecoeffs+=".txt";
    coeffs.open(filenamecoeffs+".txt");
    gnufile.open(gnufilename);

//    std::cout <<"\nfilenamecoefss = " << filenamecoeffs << "\nfilenamepoints = " << filenamepoints << "\ngnufilename = " << gnufilename;

    if(gnufile.is_open()&&coeffs.is_open()){
        gnufile << "set autoscale x\nset autoscale y ";
        gnufile << "\nset terminal png";
        gnufile << "\nset output \""<<gnuprintname<<".png\"";
        gnufile << "\nf(x)=";


        coeffs >> numberofcoeffs;
        for(int i =0; i<numberofcoeffs; i++){
                gnufile << "+";
                coeffs >> numberbridge;
                gnufile << numberbridge << "*x**" << i;
                }
        gnufile <<"\nplot '"<< filenamepoints+".txt" << "', f(x)";
        gnufile <<"\nreplot";
    }
    else if(gnufile.is_open()){
        std::cout << "\ndummyfile muss erstellt werden - filename" << filenamepoints;
        int dummyint;
        double* dummydouble;

        dummydouble = optimize(filenamepoints, dummyint);

        gnufile << "set autoscale x\nset autoscale y ";
        gnufile << "\nset terminal png";
        gnufile << "\nset output \""<<gnuprintname<<".png\"";
        gnufile << "\nf(x)=";


        numberofcoeffs = dummyint;
        for(int i =0; i<numberofcoeffs; i++){
                gnufile << "+";
                gnufile << dummydouble[i] << "*x**" << i;
                }
        gnufile <<"\nplot '"<< filenamepoints+".txt" << "', f(x)";
        gnufile <<"\nreplot";
    }


    gnufile.close();
    coeffs.close();



}

//test_gen creates a file named "filename" with k random points following a polynomial with n terms

bool test_gen(std::string filename, int k, int n){
       if (n<1){
            std::cout<<"Polynonial degree must be greater then 0"<<std::endl;
            return false;
        }

        if (k<n+1){
            std::cout<<"Number of points to generate must be greater than degree of the polynomial"<<std::endl;
            return false;
        }
        std::ofstream outfile;
        outfile.open(filename);
        outfile<<k<<std::endl;
        srand(time(NULL));

        int nsize = n+1;
        double* coeffs = new double[nsize];
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
        outfile<<std::fixed<<std::setprecision(PRECISION)<<x<<" "<<y<<std::endl;
        }
        outfile.close();
        return true;

}
//shows the testfile
bool test_gen_showfile(std::string filename){
    std::ifstream infile;
    infile.open(filename);
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
        return false;
    }


}



