#!/bin/sh
# sh scriptname 1

gcc TFG.c ./src/Process.c ./src/Signal.c ./src/Analysis.c ./src/Functions.c ./src/Methods.c -lfftw3 -lm -lpthread -lgsl -lgslcblas -g -Wall -o TFG.exe 
./TFG.exe "$1" "$2"

if [ "$1" = "af" ]
then
	vim ./output/AproxFreq.dat
fi

if [ "$1" = "ef" ]
then
	vim ./output/ExactFreq.dat
fi

if [ "$1" = "et" ]
then
	vim ./output/ExactTime.dat
fi

if [ "$1" = "st" ]
then
	vim -o ./output/Noise.dat
fi

<<COMMENT1
gnuplot -e "set term qt 0;
			set format y \"%.2e\";
			set ylabel \"Senyal\";
			set xlabel \"Temps (s)\";
			set for[i=1:5] linetype 3 dt 3;
			plot \"./output/Signal.dat\" with lines title \"h(t)+n(t)\", \"./output/Wave.dat\" with lines lt 3 title \"h(t) \";pause 1;
			set term qt 1;
			set format x \"%.2e\"; 
 			set logscale x;
			set xlabel \"log f (Hz)\";
			set format y \"%.2e\";
			set logscale y;
			set ylabel \"log(sqrt(PSD))\";
			plot \"./output/RefPSD.dat\" with lines, \"./output/PSD.dat\" with lines; pause -1;"
COMMENT1
