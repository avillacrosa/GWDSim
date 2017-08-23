//                      METHODS                         //
//------------------------------------------------------//
//Generates the noise of the detector given its PSD 	//
//using the methods described in Percival's paper		//
//	1) ExactFreq										//
//		Generates a time series with the desired PSD	//
//      by finding first its autocorrelation sequence	//
//		but it generally fails because of a non-negative//
//		condition. Therefore, almost never is useful	//
//  2) AproxFreq                                        //
//		Generates a time series with the desired PSD	//
//		by evaluating directly on the PSD and using some// 
//		Fourier properties. It is the most efficient and//
//		the most robust of the three and therefore its	//
//		use is recommended								//
//	3) ExactTime										//
//		Generates a time series with the desired PSD	//  
//		using again the autocorrelation function but 	//
//		not using any kind of DFT. This should work for //
//		PSD theoretically but practically it does not.	//
//		This might be because of the non-continuity of	//  
//		the PSD's but it has to be tested. As far as I  //
//		know it works just fine for LISA's PSD with freq//                          
//		encies close to unity							//
//////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
#include "Functions.h"
#include "Signal.h"

fftw_complex *ExactFreq ()
{
	FILE *fout;
	int M = 2*MT;//Since N<=M/2 for the method to converge we double the array length to then just talk the half
	fout = fopen("./output/ExactFreq.dat","w");
	int  k = 0;
	double W[M];
	fftw_complex *s;
	s = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*M);
	fftw_complex *V;
	V = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*M);
	fftw_complex *Vend;
	Vend = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*M);
	for (k = 0 ; k <= M/2 ; k++)
	{
		s[k][0]=ACVS(k,0);
		s[k][1]=ACVS(k,1);
        if (k != 0)
	    {
			s[M-k][0]=s[k][0];
			s[M-k][1]=s[k][1];
		}
	}
	fftw_complex *S;//From Percival's paper, it is somewhat implicit that it is purely real but i'm not sure
	S = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*M);
	S = fftw(M,s,-1);

	for (k = 0; k<M ; k++)
	{
		W[k]=Rand(1,k);//Array wind random numbers
	}
	
	for (k = 0 ; k < M ; k++)
	{
		if(S[k][0]<0)
		{
			printf(" Method failed at step %d! (PSD does not meet requierements) \n",k);
		}
		if(k==0)
		{
			V[k][0]=sqrt(S[k][0])*W[k];
			V[k][1]=0; 
		}
		if(1<=k && k<M/2)
		{
	  		V[k][0]=sqrt(0.5*S[k][0])*W[2*k-1];
			V[k][1]=sqrt(0.5*S[k][0])*W[2*k]; 
		}
		if(k==M/2)
		{
	  		V[k][0]=sqrt(S[k][0])*W[M-1];
			V[k][1]=0; 
		}
		if(M/2<k && k<M)//Complex conjugate
		{
			V[k][0]=V[M-k][0];
			V[k][1]=-V[M-k][1];
		}
	}
		
	fftw_complex *Vf;
	Vf = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*M);
	Vf = fftw(M,V,-1);
	for (k=0;k<M/2;k++)
	{
		Vend[k][0]=Vf[k][0];//Only take the real part
		Vend[k][1]=0;		//Vf[k][0] should be smaller that Vf[k][0]
		fprintf(fout," %d %.5e %.5e \n",k,Vend[k][0],Vend[k][1]);
	}
	fftw_free(V);
	fftw_free(s);
	fftw_free(S);
	fftw_free(Vf);
	fclose(fout);
	return Vend;
}

fftw_complex *AproxFreq()
{
	FILE *fout;
	fout=fopen("./output/AproxFreq.dat","w");
	int k=0, j=0;
	double W[MT],Auxi,Auxi2,x=0;
	fftw_complex *Ui;
	Ui = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*MT);
	fftw_complex *U;
	U = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*MT);
	for (k=0;k<MT;k++)
	{
		W[k]=Rand(1,k);//Random numbers
	}
	for (j=0;j<MT;j++)

	{
		if(j==0)
		{
			Ui[j][0]=sqrt(PSDC(0))*W[j];
			Ui[j][1]=0;
		}
		if(1<=j && j< MT/2)
		{
			Auxi = j;
			Auxi2 = MT;
			Ui[j][0]=sqrt(0.5*PSDC(CosF*Auxi/Auxi2))*W[2*j-1];
			Ui[j][1]=sqrt(0.5*PSDC(CosF*Auxi/Auxi2))*W[2*j];
		}
		if(j==MT/2)
		{
			Ui[j][0]=sqrt(PSDC(CosF*0.5))*W[MT-1];
			Ui[j][1]=0;
		}
		if(MT/2<j && j<=MT-1)//Complex conjugate
		{
			Ui[j][0]=Ui[MT-j][0];
			Ui[j][1]=-Ui[MT-j][1];
		}
	}
	U = fftw(MT,Ui,-1);
	k=0;
	for (x=0;x<MaxTime;x+=1/CosF)
	{
		U[k][0]=U[k][0]/sqrt(MT);
		fprintf(fout," %f %.5e %.5e \n",x,U[k][0],U[k][1]); 
		k++;
	}
	fclose (fout);	
	fftw_free(Ui);
	return U;
}

fftw_complex *ExactTime()//This is definetely hard to read from the code. I suggest looking at Percival's paper.
{
	FILE *fout;
	fout=fopen("./output/ExactTime.dat","w");	
	int N=MT;
	double sigma[N], Phi[N][N], W[N], s[N][2], SumPhi=0, SumPhiG=0;
	double x;
	fftw_complex *Y;
	Y = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N);
	int k = 0, t = 0, j = 0, g = 0;
	Phi[0][0]=0;
	for (k = 0 ; k < N ; k++)
	{
	    s[k][0]=ACVS(k,0);
		s[k][1]=ACVS(k,1);
		W[k]=Rand(1,k);
			for(j=0 ; j < N ; j++)
			{
				Phi[k][j]=0;
			}
		Y[k][0]=0;
		Y[k][1]=0;
		sigma[k]=0;
	}
	sigma[0]=sqrt(s[0][0]);
	Y[0][0]=W[0]*s[0][0];
	Y[0][1]=0;
	for (t = 1; t < N ; t++)
	{	
		if (t==1)
		{
			Phi[1][1]=s[1][0]/pow(sigma[0],2);
			sigma[1]=sigma[0]*sqrt(1-pow(Phi[1][1],2));
			SumPhi=0;
		}
		else
		{
			for (j = 1; j < t ; j++)
			{
				for (g = 1; g <= t-1;g++)
				{
					SumPhiG+=Phi[g][t-1]*s[t-g][0];
				}
				Phi[t][t]=(s[t][0]-SumPhiG)/pow(sigma[t-1],2);

				sigma[t]=sigma[t-1]*sqrt(1-pow(Phi[t][t],2));
				Phi[j][t]=Phi[j][t-1]-Phi[t][t]*Phi[t-j][t-1];
				SumPhi+=Phi[j][t]*Y[t-j][0];
				SumPhiG=0;
			}
		}
		SumPhi+=Phi[t][t]*Y[0][0];
		Y[t][0]=SumPhi+W[t]*sigma[t];
		SumPhi=0;
	}
	t=0;
	for (x=0;x<MaxTime;x+=1/CosF)
	{
		fprintf(fout," %f %e %e \n",x,Y[t][0],Y[t][1]);
		t++;
	}
	fclose(fout);
	return Y;
}
