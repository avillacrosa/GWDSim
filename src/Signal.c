//		            	SIGNAL							//
//------------------------------------------------------//
//Includes all functions for PSD's of multiple detectors//
// as well as some basic gravitational wave functions of// 
//basic and a 2.5 PN inspiral (Blanchet).				// 
//	1) PSD's											//
//		Under a certain value (Thr[2]), 				//
//		exclusive of every detector, each PSD			//	
//		function returns 0. Otherwise, the 				//
//		value taken from Schutz and Sathyaprakash		//
//		If you want to add a new PSD, you must also		//
//		include it in the ChooseSignal() function		//
//		from Process.c.									//
//	2) GW Strains										//
//		GW functions take its parameters from the		//	
//		data files in the data/ folder through the		//
//		ReadData function.								//
//		If you want to add a new GWStrain, you must		//
//		add a new one here and also in the ChooseGW()	//	
//		function from Process.c . If may or may not		//
//		generate a data file in the same way done for	//
//		my functions. To do that, keep in mind that		//
//		lines starting with '#' are read as comments	//	
//////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <fftw3.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_integration.h>

//Constants used (co=speed of light, G=grav. constant
//MPc=megaparsec, SolM=Solar mass)
#define WN 0
#define co 3.0*pow(10.0,8.0)
#define G 6.674*pow(10.0,-11.0)
#define MPc 3.08*pow(10.0,22.0)
#define SolM 1.989*pow(10.0,30.0)

//Global (accesable from every part of the code)
fftw_complex *(*Method)(int);
double (*PSDC)(double);
double *(*Strain)(double,double *);
double CosF;
int MT;
//int MaxTime;
double MaxTime;
double Thr[3];
char *GWPath;

//DETECTOR PSD's
double Lisa(double x)//Threshold=10^-5Hz; 
{
	if (x<=Thr[0])
	{
		return 0;
	}
	else
	{
		double q = 1.59*pow(10,-41)+9.18*pow(10,-52)/pow(x,4)+9.18*pow(10,-38)*pow(x,2);    
		return q;
	}
}

double ILigo(double x)//Thresh=40Hz
{
	if(x<=Thr[0])
		return 0;
	else
	{
		double xf=x/150.0;
	//	double p = 9.0*pow(10,-46)*(4.49*pow(xf,-56.0)+0.16*pow(xf,-4.52)+0.52+0.32*pow(xf,2.0));//First term introduces high values for low freqs.	
		double p = 9.0*pow(10,-46)*(0.16*pow(xf,-4.52)+0.52+0.32*pow(xf,2.0));//Same as commented p of the last line but without first term	
		return p;
	}
}

double ALigo(double x)//Thresh=20Hz
{
	if(x<=Thr[0])
		return 0;
	else
	{
		double xf=x/215.0;
		double p = 1.0*pow(10.0,-49.0)*(pow(xf,-4.14)-5.0*pow(xf,-2.0)+111.0*(1.0-pow(xf,2.0)+0.5*pow(xf,4.0))/(1.0+0.5*pow(xf,2.0)));
		return p;
	}
}

double TAMA(double x)//Thresh=75Hz
{
	if(x<=Thr[0])
	{
		return 0;
	}
	else
	{
		double xf=x/400;
		double p = 7.5*pow(10,-46)*(pow(xf,-5)+13.0/xf+9.0*(1.0+pow(xf,2)));
		return p;
	}
}

double VIRGO(double x)//Thresh=20Hz
{
	if(x<=Thr[0])
	{
		return 0;
	}
	else
	{
		double xf=x/500;
		double p=3.2*pow(10,-46)*(pow((7.8*xf),-5)+2.0/x+0.63+pow(xf,2));
		return p;
	}
}


double *ReadFromData()//Functions used to read GW parameters. Order is important since there is no given information anywhere on what is each parameter
{
	char line[100];
	char s1[10];
	FILE *dataf;
	int k = 0;
	dataf=fopen(GWPath,"r");
	double *DataParams = malloc(sizeof(double) * 10);//M1,M2,r,i,phi_o,omega_o
	while (fgets(line, sizeof line, dataf))
	{
		if (*line=='#') continue;
		if(sscanf(line,"%9s",s1)!=1)
		{
			printf("Error while reading wave file");
		}
		else
		{
			sscanf(s1, "%lf",&DataParams[k]);
		}
		k++;
	}
	fclose(dataf);
	return DataParams;
}

//GRAVITATIONAL WAVE FUNCTIONS

double Inspiral(double x,double *Params)//From C.J. Moore, hasn't been tested
{	
	double fdot = 3,ho=5;
	double result = sqrt(2/fdot)*x*ho;
	return result;
}

double Stochastic(double x, double *Params)//From C.J. Moore, hasn't been tested
{
	double fo=12,alpha=-0.5,A=20;
	double result = A*pow((x/fo),alpha);
	return result;
}

double AdvInspiral(double t,double *Params)//From Blanchet, has been tested. It might blow up for t close to MaxTime and high sampling frequencies
{	
	double M1=Params[0]*SolM,M2=Params[1]*SolM,r=Params[2]*MPc,i=Params[3],PhiC=Params[4],Omega0=Params[5];
	double Tc=MaxTime;
	double Theta,Phi,Omega,Psi,x,M,eta,DeltaM;//Derived parameters
	double Hplus[5],Hcross[5];//PN Terms
	double hplus,hcross;//Strains
	double Azimuth = M_PI;//Angles of observation for antenna patern functions
	double OrbAngle = M_PI;//Angle of observation for antenna patern functions
	DeltaM=M1-M2;
	M=M1+M2;
	eta = (M1*M2)/pow(M,2.0);
	double s = sin(i);
	double c = cos(i);
	double h = 0,Fplus,Fcross;
	
	Fplus=0.5*(1+pow(c,2.0))*cos(2*Azimuth)*cos(2*OrbAngle)-c*sin(2*Azimuth)*sin(2*OrbAngle);
	Fcross=0.5*(1+pow(c,2.0))*cos(2*Azimuth)*cos(2*OrbAngle)+c*sin(2*Azimuth)*sin(2*OrbAngle);

	Theta=pow(co,3.0)*eta*(Tc-t)/(5.0*G*M);
	
	Phi=PhiC-(1.0/eta)*(pow(Theta,5.0/8.0)+(3715.0/8064.0+55.0*eta/96.0)*pow(Theta,3.0/8.0)-3.0*M_PI*pow(Theta,1.0/4.0)/4.0+
	(9275495.0/14450688.0+284875.0*eta/258048.0+1855.0*pow(eta,2.0)/2048.0)*pow(Theta,1.0/8.0));

	Omega=(pow(co,3.0)/(8.0*G*M))*((pow(Theta,-3.0/8.0))+(743.0/2688.0+11.0*eta/32.0)*pow(Theta,-5.0/8.0)-3.0*M_PI*pow(Theta,-3.0/4.0)/10.0
	+(1855099.0/14450688.0+56975.0*eta/258048.0+371.0*pow(eta,2.0)/2048.0)*pow(Theta,-7.0/8.0));//OMEGA IS BECOMING NEGATIVE BUt STILL A NUMBER
	
	x=pow(G*M*Omega/pow(co,3.0),2.0/3.0);
	Psi=Phi-(2.0*G*M*Omega/pow(co,3.0))*log(Omega/Omega0);	
	
	Hplus[0]=-(1.0+pow(c,2.0))*cos(2.0*Psi);
	Hplus[1]=-(s*DeltaM/(8.0*M))*((5.0+pow(c,2.0))*cos(Psi)-9.0*(1.0+pow(c,2.0))*cos(3.0*Psi));
	Hplus[2]=1.0/6.0*(19.0+9.0*pow(c,2.0)-2*pow(c,4.0)-eta*(19.0-11.0*pow(c,2.0)-6.0*pow(c,4.0)))*cos(2.0*Psi)
	-4.0/3.0*pow(s,2.0)*(1.0+pow(c,2))*(1.0-3.0*eta)*cos(4.0*Psi);		
	Hplus[3]=s/192.0*DeltaM/M*((57.0+60.0*pow(c,2.0)-pow(c,4.0)-2.0*eta*(49.0-12.0*pow(c,2.0)-pow(c,4.0)))*cos(Psi)
	-(27.0/2.0)*(73.0+40.0*pow(c,2)-9.0*pow(c,4.0)-2.0*eta*(25.0-8.0*pow(c,2.0)-9.0*pow(c,4.0)))*cos(3.0*Psi)
	+625.0/2.0*(1.0-2.0*eta)*pow(s,2)*(1.0+pow(c,2.0))*cos(5.0*Psi))-2.0*M_PI*(1.0+pow(c,2.0))*cos(2.0*Psi);
	Hplus[4]=1.0/120.0*(22.0+396.0*pow(c,2)+145.0*pow(c,4.0)-5.0*pow(c,6)+5.0/3.0*eta*(706.0-216.0*pow(c,2.0)
	-251.0*pow(c,4)+15.0*pow(c,6))-5.0*pow(eta,2.0)*(98.0-108.0*pow(c,2)+7.0*pow(c,4.0)+5.0*pow(c,6.0)))*cos(2.0*Psi)
	+2.0/15.0*pow(s,2.0)*(59.0+35.0*pow(c,2.0)-8*pow(c,4.0)-5.0/3.0*eta*(131.0+59.0*pow(c,2)-24.0*pow(c,4.0))
	+5.0*pow(eta,2.0)*(21.0-3*pow(c,2.0)-8.0*pow(c,4.0)))*cos(4.0*Psi)-81.0/40.0*(1.0-5.0*eta+5.0*pow(eta,2))*pow(s,4.0)
	*(1.0+pow(c,2.0))*cos(6.0*Psi)+s/40.0*DeltaM/M*((11.0+7.0*pow(c,2)+10*(5.0+pow(c,2.0))*log(2.0))*sin(Psi)
	-5.0*M_PI*(5.0+pow(c,2.0))*cos(Psi)-27.0*(7.0-10.0*log(3.0/2.0))*(1.0+pow(c,2.0))*sin(3.0*Psi)+135.0*M_PI*(1+pow(c,2.0))*cos(3.0*Psi));
	
	Hcross[0]=-2.0*c*sin(2.0*Psi);
	Hcross[1]=-3.0/4.0*s*c*DeltaM/M*(sin(Psi)-3.0*sin(3.0*Psi));
	Hcross[2]=c/3.0*(17.0-4.0*pow(c,2)-eta*(13.0-12.0*pow(c,2.0)))*sin(2.0*Psi)-8.0/3.0*(1.0-3.0*eta)*c*pow(s,2.0)*sin(4.0*Psi);
	Hcross[3]=s*c*DeltaM/(96.0*M)*((63.0-5.0*pow(c,2.0)-2.0*eta*(23.0-5.0*pow(c,2.0)))*sin(Psi)
	-27.0/2.0*(67.0-15.0*pow(c,2.0)-2.0*eta*(19.0-15.0*pow(c,2.0)))*sin(3.0*Psi)+625.0/2.0*(1.0-2.0*eta)*pow(s,2.0)*sin(5.0*Psi))
	-4.0*M_PI*c*sin(2.0*Psi);
	Hcross[4]=c/60.0*(68.0+226.0*pow(c,2)-15.0*pow(c,4.0)+5.0/3.0*eta*(572.0-490.0*pow(c,2.0)+45.0*pow(c,4.0))
	-5.0*pow(eta,2.0)*(56.0-70.0*pow(c,2.0)+15.0*pow(c,4.0)))*sin(2.0*Psi)+4.0/15.0*c*pow(s,2.0)*(55.0-12.0*pow(c,2.0)
	-5.0/3.0*eta*(119.0-36.0*pow(c,2.0))+5.0*pow(eta,2.0)*(17-0-12.0*pow(c,2.0)))*sin(4.0*Psi)-81.0/20.0*(1.0-5.0*eta
	+5.0*pow(eta,2.0))*c*pow(s,4.0)*sin(6.0*Psi)-3.0/20.0*s*c*DeltaM/M*((3.0+10.0*log(2.0))*cos(Psi)+5.0*M_PI*sin(Psi)
	-9.0*(7.0-10.0*log(3.0/2.0))*cos(3*Psi)-45.0*M_PI*sin(3.0*Psi));    
	
	hplus = (2.0*G*M*eta/(pow(co,2.0)*r))*pow(G*M*Omega/pow(co,3.0),2.0/3.0)*(Hplus[0]+pow(x,0.5)*Hplus[1]+x*Hplus[2]+pow(x,1.5)*Hplus[3]+pow(x,2.0)*Hplus[4]);
	hcross = (2.0*G*M*eta/(pow(co,2.0)*r))*pow(G*M*Omega/pow(co,3.0),2.0/3.0)*(Hcross[0]+pow(x,0.5)*Hcross[1]+x*Hcross[2]+pow(x,1.5)*Hcross[3]+pow(x,2.0)*Hcross[4]);
	h=hplus*Fplus+Fcross*hcross;
	if (h!=h)//In the case where h is -nan, we return 0 instead.
		return 0;
	return h;
}
