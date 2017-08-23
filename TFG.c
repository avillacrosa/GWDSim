/*Structure of run, rest is found inside src folder*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <fftw3.h> //For fast fourier transforms
#include <gsl/gsl_math.h> //For mathematical functions
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_integration.h>
#include "./src/Process.h"
#include "./src/Analysis.h"


int main (int argc, char *argv[])
{
	double SNR,SNRFilt;
	int i = 0;
	ChooseDetector();
	ChooseMethod();
	ChooseGW();
	printf("Time of measurement in seconds (integer) : ");
	scanf(" %d",&MaxTime);//Maximum time of measurement
	printf("Sampling frequency in Hz (integer) : ");
	scanf(" %lf",&CosF);//Sampling frequency
	MT=CosF*MaxTime;//Number of points in the final signal series
	fftw_complex *Sig;
	Sig = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*MT);
	for (i = 0;i<MT;i++)
	{
		Sig[i][0]=0;
		Sig[i][1]=0;
	} 
	printf("---GENERATING TIME SERIES--- \n");
	Sig = GetSignal(MT);
	FourierPhases();
	printf("---CALCULATING SNR--- \n");
	SNR=SNRInt(ReadFromData());
	printf("Approximate signal to noise ratio (SNR) : %f \n",SNR);
	PSD(MT);
	RefPSD(MT);
	SNRFilt=ProtoMatchedFilter();
	printf("Matched (Wiener) filter SNR : %f \n",SNRFilt);
	if (SNR==SNRFilt)
		printf("Parameters are perfectly reconstructed (Wiener SNR is equal to the SNR)  \n");
	printf("END \n");
	return 0;
}
