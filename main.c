#include "debug_info.h"
#include "gf_cal.h"
#include "as_decoding.h"
#include "encoding.h"

void main()
{
	unsigned char i = 0;

#if (1 == SYS_ENC)
	systematic_encoding();
#else
	evaluation_encoding();
#endif

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
