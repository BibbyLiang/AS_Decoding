#include <stdio.h>
#include <stdlib.h>
#include "debug_info.h"
#include "gf_cal.h"
#include "rnd.h"
#include "math.h"

float awgn_gen(float snr)
{
	float val = 0;

	#define E_B		1
	float w = 0, r = 0;

	w = (float)rand() / (float)RAND_MAX;
	if (w == 1.0)
	{
		w = 0.999999;
	}

	r = gaussrand();

	/*convert it to modulation constellation*/
	val = (float)r * cos(2 * PI * w);

	val = val * (sqrt(E_B / pow(10, snr / 10) / 2));

	//printf("gauss_val: %f\n", val);

	return val;
}
