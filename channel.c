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

#if 1
	r = gaussrand();
	val = r * (sqrt(E_B / pow(10, snr / 10) / 2));
#else
	r = sqrt(2.0 * log(1.0 / (1.0 - w)));

	/*convert it to modulation constellation*/
	val = r * (float)cos(2 * PI * w);

	val = val * (sqrt(E_B / pow(10, snr / 10) / 2));
#endif
	DEBUG_NOTICE("val: %f %f %f\n", r, w, val);
	DEBUG_NOTICE("gauss_val: %f\n", val);

	return val;
}
