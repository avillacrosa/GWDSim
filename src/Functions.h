#ifndef FUNCTIONS
#define FUNCTIONS

//extern fftw_complex *(*Method)(int);
//extern double (*PSDC)(double);
//extern double (*Strain)(double);
//extern double CosF;

double RealFourier(double x, void * params);
double ComplexFourier(double x, void * params);
double ACVS (int param, int Part);
float Rand(int N, int j);
double *Polar (double Real, double Imag);
fftw_complex *fftw(int N, fftw_complex *in, int Direction);

#endif
