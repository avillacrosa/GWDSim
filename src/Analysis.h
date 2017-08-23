#ifndef ANALYSIS
#define ANALYSIS

fftw_complex *(*Method)(int);


double PSD();
double RefPSD();
double SNRInt (double *Params);
void FourierPhases ();
double ProtoMatchedFilter ();

#endif
