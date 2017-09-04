//                      ANALYSIS                        //
//------------------------------------------------------//
//Calculates desired parameters for the signal output,  //
//such as signal to noise ratio, fourier phases, and 	// 
//a very basic matched (Wiener) filter					// 
//  1) PSD                                   			//
// 		Calculates the PSD of the generated noise by the//
// 		method of the periodogram. 						//  
//  2) RefPSD                                   		//
//		Calculates the "reference" PSD from the function//
// 		in Signal.c (PSDC).								//
//  3) SNRInt                                           //
// 		Calculates the Signal to Noise Ratio through 	//  
// 		basic rectangular numerical integration			//
//  4) FourierPhases                                    //
// 		Calculates the phases between complex and real 	//
// 		part of the fourier transform of the noise to	//  
// 		search for correlations							//
// 	5) ProtoMatchedFilter								//							
//		Calculates the mass of one of the objects of the//
//		inspiral that maximizes the overlap. This is 	//
//		more than anything an exercise since all of the //
//		other parameters are fixed.						//
//////////////////////////////////////////////////////////


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_integration.h>
#include <time.h>
#include "Process.h"

#define SolM 1.989*pow(10.0,30.0);
#define G 6.674*pow(10.0,-11.0);
#define co 3.0*pow(10.0,8.0);
#define MPc 3.08*pow(10.0,22.0);

double PSD(int M)
{
	printf("---CALCULATING PSD--- \n");
	FILE *fout;
	fout=fopen("./output/PSD.dat","w");
	int N,L;
	fftw_complex *Result;
	Result = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*MT);
	double Modulus[M];
	double k = 1000;
	double MeanPSD,x;
		for (N=0;N<MT;N++)//Arrays must be initialized
		{
			Modulus[N]=0;
		}
		for (N=0;N<k;N++)//Realizations (average of the ensemble). Modifiyable
		{
			Result = fftw(MT,Method(MT),-1);
			for (L=0;L<MT/2;L++)//Only runs up to MT/2 since it is symmetric after that point
			{
				Result[L][0]=Result[L][0]/sqrt(CosF*MT);
				Result[L][1]=Result[L][1]/sqrt(CosF*MT);
				Modulus[L]+=pow(Result[L][0],2)+pow(Result[L][1],2);
			}
		}
		for (N=0;N<MT/2;N++)
		{
			x=N;
			x=x/(MaxTime);
			Modulus[N]=sqrt(Modulus[N]*CosF)/sqrt(k);//TAKING SQRT OF PSD!!!
			fprintf(fout," %f %e \n",x,Modulus[N]);
			MeanPSD+=Modulus[N];	
		}
		MeanPSD=MeanPSD/(MT/2);
	fclose(fout);
	fftw_free(Result);
	return 0;
}

double RefPSD ()
{
	printf("---CALCULATING REFERENCE PSD--- \n");
	FILE *fout;
	fout=fopen("./output/RefPSD.dat","w");
	double x,step=0.0001;
	for (x=Thr[0]+step;x<Thr[1];x+=step)//Recall that Thr[0] is the minimun frequency and Thr[1] the maximum one
	{
		if (x<=(100*Thr[0]))
		{
			step=Thr[0];
		}
		else
			step = 100*Thr[0];
		fprintf(fout," %f %e \n",x,sqrt(PSDC(x)));
	}	
	return 0;
}

double SNRInt (double *Params)
{
    int j=0,N=0;
    double SNR=0;
	double StrainMod[MT];
	fftw_complex *Wave;
	Wave = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*MT);
	N=0;
	for(N=0;N<MT;N++)
	{
		Wave[N][0]=Strain(N/CosF,Params);
		Wave[N][1]=0;
	}

    fftw_complex *StrainSNR;
    StrainSNR = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*MT);
    StrainSNR=fftw(MT,Wave,-1);//Fourier transform of the wave

    for (j=0;j<MT;j++)
    {	

		StrainSNR[j][0]=StrainSNR[j][0]/pow(MT,1);//Normalization of the DFT
		StrainSNR[j][1]=StrainSNR[j][1]/pow(MT,1);

		if(PSDC(CosF*j/MT)!=0 && PSDC(CosF*j/MT)==PSDC(CosF*j/MT))
		{
        	StrainMod[j]=(pow(StrainSNR[j][0],2)+pow(StrainSNR[j][1],2));//Modulus
        	SNR+=StrainMod[j]/(PSDC(j*CosF/MT))*CosF/MT;//Definition of SNR
		}
    }
	SNR=sqrt(4.0*SNR);
	fftw_free(Wave);
	fftw_free(StrainSNR);
    return SNR;
}

void FourierPhases ()
{
	printf("---CALCULATING FOURIER PHASES--- \n");
	FILE *fout;
	fout=fopen("./output/Phases.dat","w");
	int j=0;
	double x=0,Phase=0;
	fftw_complex *NoiseFFT;
	NoiseFFT=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*MT);
	NoiseFFT=fftw(MT,Method(MT),-1);
	for(x=0;x<CosF;x+=CosF/MT)
	{
		Phase=atan(NoiseFFT[j][1]/NoiseFFT[j][0]);
		if(NoiseFFT[j][0] > 0 && NoiseFFT[j][1] > 0)
				Phase=Phase;
		if(NoiseFFT[j][0] < 0 && NoiseFFT[j][1] > 0)
				Phase=Phase+M_PI;	
		if(NoiseFFT[j][0] > 0 && NoiseFFT[j][1] < 0)
				Phase=Phase;
		if(NoiseFFT[j][0] < 0 && NoiseFFT[j][1] < 0)
				Phase=Phase-M_PI;
		fprintf(fout," %f %f \n",x,Phase);//Phas is just the arctangent of the complex number by the real number
		j++;	
	}
	fclose(fout);
	fftw_free(NoiseFFT);
}

double ProtoMatchedFilter ()
{
	printf("---CALCULATING WIENER FILTER--- \n");
	//Recall that order is important and not given in data. For AdvInspiral : 0=M1,1=M2,2=r,3=i,4=PhiC,5=Omega0
	FILE *fout;
	fout=fopen("./output/Filter.dat","w");
	int i,j,lag=0;
	double SNR=0,spre=0,Mass,smod;
	double Params[6];
	double s[2];
	s[0]=0;
	s[1]=0;
	fftw_complex *DetOut;//Detector output
	fftw_complex *WaveM;//Wave
	fftw_complex *FWave;//Fourier transform of the wave
	fftw_complex *Wiener;//Wiener Filter
	fftw_complex *FWiener;//Fourier transform of the wiener filter
	WaveM=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*(MT));
	DetOut=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*(MT));
	FWave=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*(MT));
	Wiener=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*(MT));
	FWiener=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*(MT));
	for (j=0;j<100;j++)//Every j adds a mass to the Wiener filter (100 in this case). This can be modified with no problem
	{
		Params[0]=36.0;
		Params[1]++;
		Params[2]=410.0;
		Params[3]=1.42;
		Params[4]=0.0;
		Params[5]=1.0;
		DetOut=GetSignal();
		i=0;
		for(i=0;i<MT;i++)
		{
			WaveM[i][0]=Strain(i/CosF,Params);
			WaveM[i][1]=0;
		}
		FWave=fftw(MT,WaveM,-1);
		for(i=0;i<MT;i++)
		{
			FWave[i][0]=FWave[i][0]/MT;//Not necessary since we're maximizing but still...
			FWave[i][1]=FWave[i][1]/MT;//
			if(PSDC(CosF*i/MT)!=0.0)
			{
				FWiener[i][0]=FWave[i][0]/PSDC(CosF*i/MT);
				FWiener[i][1]=FWave[i][1]/PSDC(CosF*i/MT);
			}
			if(PSDC(CosF*i/MT)==0.0)	
			{
				FWiener[i][0]=0;
				FWiener[i][1]=0;
			}	
		}
		Wiener=fftw(MT,FWiener,1);
		for(i=0;i<MT;i++)
		{
			if((i-lag)>0)//lag is generally zero. In the future I may add the option to make the wave of a different time length than the noise for which case this will be important.
			{
				s[0]+=DetOut[i][0]*Wiener[i-lag][0];	
				s[1]+=DetOut[i][0]*Wiener[i-lag][1];
			}
		}
		smod=sqrt(pow(s[0],2)+pow(s[1],2));
		if(spre<s[0])//If the overlap is higher than before
		{	
			fprintf(fout," Mass = %d \n",j);
			for(i=0;i<MT;i++)
			{
				fprintf(fout," %.2e %.10e %.2e %.2e \n",Wiener[i][0],DetOut[i][0],Wiener[i][1],DetOut[i][0]*Wiener[i][0]);
			}
			fprintf(fout," s0 %.2e s1 %.2e smod %.2e\n \n",s[0],s[1],smod);
			Mass=j+1;
			spre=s[0];
		}
		s[0]=0;
		s[1]=0;
	}
	Params[1]=Mass;
	SNR=SNRInt(Params);
	printf("Maximum happened for M2= %f \n",Mass);
	fftw_free(WaveM);
	fftw_free(DetOut);
	fftw_free(FWave);
	fftw_free(Wiener);
	fftw_free(FWiener);
	return SNR;
}
