//                      FUNCTIONS                       //
//------------------------------------------------------//
//Includes all kinds of functions that are used throug- //
//out the code. 										// 
//  1)RealFourier/ComplexFourier                        //
//		Prepares the PSD for the autocorrelation 		//
//		integration (fourier transform) by multiplying	//
//		for exp(-i*2*pi*f*tau)							//
//  2)ACVS    											//
//		Performs a GSL numerical integration to find    //
//		the autocorrelation sequence					//
// 	3)Rand												//
// 		Generates a random number using time as random	//
// 		seeding											//
// 	4)fftw												//  
// 		Calculates the fast fourier transform of the	//
// 		given array. The forward DFT is not normalized	//
//		wheres the backwards is (divided by MT)			//
//////////////////////////////////////////////////////////



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_integration.h>
#include <time.h>
#include "Signal.h"


double RealFourier(double x,void * params)
{
	int i = *(int *)params;
	double p = PSDC(x*CosF)*cos(2*M_PI*x*i);//Source of all trouble...
	return p;
}

double ComplexFourier(double x, void * params)
{	
	int i = *(int *)params;
	double p = CosF*PSDC(x)*sin(2*M_PI*x*i/CosF);
	return p;	
}

double ACVS(int param,int Part)
{
	gsl_integration_workspace *w=gsl_integration_workspace_alloc(10000);
	gsl_function INTEG;

	INTEG.params=&param;	
	double ACVSR=0,error;
	double relerror=Thr[2],abserror=0;//This can freeze program, make integrations not very precise, quite of a chaos
	if (Part==0)	
	{
		INTEG.function=&RealFourier;
		gsl_integration_qag(&INTEG,0,0.5,abserror,relerror,10000,6,w,&ACVSR,&error);
		gsl_integration_workspace_free (w);
		return 2*ACVSR;
	}
	if (Part==1)
	{
		gsl_integration_workspace_free (w);
		return 0;//We return 0 because the PSD is symmetric and the sine antisymmetric
		INTEG.function=&ComplexFourier;
		gsl_integration_qag(&INTEG,-0.5,0.5,abserror,relerror,10000,6,w,&ACVSR,&error);
		gsl_integration_workspace_free (w);
		//return 2*ACVSR;
	}
	return 0;
}

float Rand (int N, int j)
{
	const gsl_rng_type *T;
	gsl_rng *r;
	gsl_rng_env_setup();
	T=gsl_rng_default;
	r=gsl_rng_alloc(T);
	gsl_rng_set(r,time(NULL));//Random seeding
//	gsl_rng_set(r,0);//Non-Random seeding
	int i = 0;
	double f = 0;
	for (i=0;i<j+1;i++)
	{	
		f=gsl_ran_gaussian(r,N);
	}
	gsl_rng_free(r);
	return f;
}

double *Polar(double Real,double Imag)
{
	double *P=malloc(sizeof(double)*2);
	P[0]=0;
	P[1]=0;
	P[0] = sqrt(pow(Real,2)+pow(Imag,2));
	P[1] = atan(Imag/Real);
	return P;
}

fftw_complex *fftw(int N,fftw_complex *in,int Direction)
{
	fftw_complex *out;
	fftw_plan p;
	out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*(N));
	p = fftw_plan_dft_1d(N,in,out,Direction,FFTW_ESTIMATE);
	fftw_execute(p);
	fftw_destroy_plan(p);
	return out;
}

