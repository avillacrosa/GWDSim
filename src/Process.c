//                      PROCESS                         //
//------------------------------------------------------//
//Manages all inputs for the program, sets the global	//
//variables and returns in the end the desired signal	// 
//It does NOT perform any kind of analysis				// 
//  1) ChooseDetector()                                 //
// 		Prints all the available detectors and let's you//
// 		choose one of them, stored in the global varia-	//  
// 		ble	PSDC										//
// 	2) ChooseMethod()									//
// 		Same as before but with the methods (exact time,//
// 		aproximate frequency and exat frequency). Method//
// 		is stored in the global variable Method			//
//  3) ChooseGW                                         //
// 		Same with gravitational wave functions, this 	//  
// 		time stored in global GWStrain					//
// 	4) GetSignal()										//
//		Generates the time series that includes the 	//
// 		noise and the gravitational wave, and stores	//  
// 		each of them in Noise.dat, Wave.dat and 		//
// 		Signal.dat										//
//////////////////////////////////////////////////////////


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>
#include "Methods.h"
#include "Signal.h"
#include "Functions.h"

//Typedefs of PSD (SFunc), method (MFunc), strain (StrFunc)
typedef double (*SFunc)(double);
typedef fftw_complex *(*MFunc)(int);
typedef double (*StrFunc)(double,double*);

void ChooseDetector()
{
	int i = 1;
	int Sig = 0;
	char StrSig[3];
	char *Signals[50];
	SFunc SignalsAddress[5];
	SignalsAddress[1]=&Lisa;
	SignalsAddress[2]=&ILigo;	
	SignalsAddress[3]=&ALigo;	
	SignalsAddress[4]=&TAMA;	
	SignalsAddress[5]=&VIRGO;	

	Signals[1]="Lisa                (10^-5 to 10 Hz)";
	Signals[2]="Initial Ligo        (40 to 10^4 Hz)";
	Signals[3]="Advanced Ligo       (20 to 10^4 Hz)";
	Signals[4]="TAMA                (75 to 10^4 Hz)";
	Signals[5]="VIRGO               (20 to 10^4 Hz)";
	
	printf("Choose a PSD to generate the noise \n");
	for (i=1;i<6;i++)
	{
		printf(" %d - %s \n",i,Signals[i]);
	}
	printf("Write the number for the PSD that you want : ");
	fgets(StrSig,3,stdin);
	sscanf(StrSig, "%d",&Sig);
	PSDC=SignalsAddress[Sig];	

	switch (Sig)//Threshholds of every detector. Thr[0] = Lowest frequency, Thr[1] = Maximum frequency, Thr[2] = Precision necessary for integration (ExactTime method)
	{
		case 1: //Lisa
			Thr[0]=pow(10.0,-5.0);
			Thr[1]=10;
			Thr[2]=1e-9;
			break;
		case 2: //ILigo
			Thr[0]=40.0;
			Thr[1]=pow(10,4);
			Thr[2]=1e-9;
			break;
		case 3: //ALigo
			Thr[0]=20.0;
			Thr[1]=pow(10,4);
			Thr[2]=1e-9;
			break;
		case 4: //TAMA
			Thr[0]=75.0;
			Thr[1]=pow(10,4);
			Thr[2]=1e-9;
			break;
		case 5: //VIRGO
			Thr[0]=20;
			Thr[1]=pow(10,4);
			Thr[2]=1e-9;
			break;
		default :	
			break;
	}	
}

void ChooseMethod()
{
	int i=1,Met=0;
	char StrMet[3];
	char *Methods[4];
	MFunc MethodsAddress[4];
	Methods[1]="Exact Frequency";
	Methods[2]="Aproximate Frequency (Most efficient and reliable)";
	Methods[3]="Exact Time";
	
	MethodsAddress[1]=&ExactFreq;
	MethodsAddress[2]=&AproxFreq;
	MethodsAddress[3]=&ExactTime;
	printf("Choose a Method to generate the time series \n");
	for (i=1;i<4;i++)
	{
		printf(" %d - %s \n",i,Methods[i]);
	}
	printf("Write the number for the method that you want : ");	
	fgets(StrMet,3,stdin);
	sscanf(StrMet, "%d",&Met);
	Method=MethodsAddress[Met];
}

void ChooseGW()
{
	int i=1,GWi=0;
	char* GW[20];
	char* GWData[50];
	StrFunc GWStrain[10];
	
	GW[1]="Advanced Inspiral (2.5 PostNewtonian (PN))";
	GW[2]="GW150914 (First observation of GW) (2.5 PN)";
	GW[3]="Stochastic background";
	GW[4]="Basic Inspiral";

	GWData[1]="./data/Inspiral.dat";
	GWData[2]="./data/GW150914.dat";
	GWData[3]="./data/StochasticBG.dat";
	GWData[4]="./data/BasicInspiral.dat";
	
	GWStrain[1]=&AdvInspiral;
	GWStrain[2]=&AdvInspiral;
	GWStrain[3]=&Stochastic;	
	GWStrain[4]=&Inspiral;	

	printf("Choose the source of the GW : \n");
	for(i=1;i<5;i++)
	{
		printf(" %d - %s \n",i,GW[i]);
	}
	printf("Write the number for the selected GW : ");
	scanf(" %d",&GWi);
	GWPath=GWData[GWi];
	Strain=GWStrain[GWi];
}


fftw_complex *GetSignal()
{
	FILE *fout;
	FILE *foutS;
	FILE *foutN;
	fout=fopen("./output/Signal.dat","w");
	foutN=fopen("./output/Noise.dat","w");
	foutS=fopen("./output/Wave.dat","w");
	double x=MT;
	int N=0;
	int Len = x+1;
	fftw_complex *Noise;
	fftw_complex *Wave;
	fftw_complex *Signal;
	
	Noise=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*(Len));
	Wave=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*(Len));
	Signal=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*(Len));
	Noise=Method(MT);
	for(N=0;N<MT;N++)
	{	
		Wave[N][0]=Strain(N/CosF,ReadFromData());
		Wave[N][1]=0;
		Signal[N][0]=Wave[N][0]+Noise[N][0];
		Signal[N][1]=0;
		fprintf(fout," %.2e %.5e \n",N/CosF,Signal[N][0]);
		fprintf(foutN," %.2e %.5e \n",N/CosF,Noise[N][0]);
		fprintf(foutS," %.2e %.5e \n",N/CosF,Wave[N][0]);
	}
	fclose(fout);
	fclose(foutS);
	fclose(foutN);
	fftw_free(Noise);
	fftw_free(Wave);
	return Signal;
}
