#ifndef SIGNAL
#define SIGNAL

extern fftw_complex *(*Method)(int);
extern double (*PSDC)(double);
extern double (*Strain)(double,double *);
extern double CosF;
extern int MT;
//extern int MaxTime;
extern double MaxTime;
extern double Thr[3];
extern char *GWPath;

double Lisa(double x);
double ILigo(double x);
double ALigo(double x);
double TAMA(double x);
double VIRGO(double x);
double *ReadFromData();
double Inspiral(double x,double *Params);
double Stochastic(double x,double *Params);
double AdvInspiral(double x,double *Params);

#endif
