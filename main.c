#include "gf_cal.h"
#include "as_decoding.h"
#include "encoding.h"

void main()
{
	unsigned char i = 0;

	systematic_encoding();

	/*transmission through channel*/
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		received_polynomial[i] = gf_add(encoded_polynomial[i], error_polynomial[i]);
	}

	mul_assign();
	re_encoding();
	
	as_decoding();

	return;
}
