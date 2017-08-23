
#include "Methods.h"
#include "Signal.h"
#include "Functions.h"

typedef double (*SFunc)(double);
typedef fftw_complex *(*MFunc)(int);
typedef double *(StrFunc)(double,int);

void ChooseDetector ();
void ChooseMethod ();
void ChooseGW();
fftw_complex *GetSignal ();
