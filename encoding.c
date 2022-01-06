#include <stdio.h>
#include <string.h>
#include "debug_info.h"
#include "gf_cal.h"
#include "encoding.h"

#if (3 == GF_Q)
/*n = 7, k = 3*/
unsigned char generator_polynomial[CODEWORD_LEN - MESSAGE_LEN + 1] =
{
	/*generated by MATLAB: genpoly = rsgenpolycoeffs(n, k), in add-style*/
#if 0	
	0x3,
	0x1,
	0x0,
	0x3,
	0x0
#else
	1,
	3,
	1,
	2,
	3
#endif
};
#endif

#if (4 == GF_Q)
/*n = 15, k = 7*/
unsigned char generator_polynomial[CODEWORD_LEN - MESSAGE_LEN + 1] =
{
	1,
	9,
	4,
	3,
	4,
	13,
	6,
	14,
	12
};
#endif

#if (6 == GF_Q)
/*n = 63, k = 21*/
unsigned char generator_polynomial[CODEWORD_LEN - MESSAGE_LEN + 1] =
{
	/*generated by MATLAB: genpoly = rsgenpolycoeffs(n, k), in add-style*/
	1,
	19,
	36,
	57,
	60,
	5,
	49,
	18,
	3,
	2,
	56,
	15,
	33,
	41,
	25,
	40,
	19,
	31,
	63,
	34,
	2,
	51,
	45,
	13,
	49,
	54,
	58,
	59,
	37,
	21,
	12,
	6,
	49,
	10,
	49,
	32,
	10,
	6,
	63,
	40,
	9,
	62,
	59
};
#endif

unsigned char message_polynomial[MESSAGE_LEN] = 
{
	/*power representation*/
	0x1,
	0x2,
	0x0
};

unsigned char encoded_polynomial[CODEWORD_LEN] =
{
	0xFF,
	0xFF,
	0xFF,
	0xFF,
	0xFF,
	0xFF,
	0xFF
};

unsigned char error_polynomial[CODEWORD_LEN] =
{
	0xFF,
	0xFF,
	0xFF,
	0xFF,
	0xFF,
	0xFF,
	0xFF
};

unsigned char evaluation_encoding()
{
	unsigned char i = 0, j = 0;
	unsigned char tmp_power_result[MESSAGE_LEN], tmp_add_result = 0;

	unsigned char codeword[CODEWORD_LEN];

	for(i = 0; i < CODEWORD_LEN; i++)
	{
		for(j = 0; j < MESSAGE_LEN; j++)
		{
			tmp_power_result[j] = gf_multp(message_polynomial[j], (i * gf_location(j)) % (GF_FIELD - 1));
			//if(0xd == i)printf("tmp_power_result: %d = %d * %d\n", tmp_power_result[j], message_polynomial[j], i * gf_location(j));
		}

		tmp_add_result = 0xFF;
		for(j = 0; j < MESSAGE_LEN; j++)
		{
			//tmp_add_result = gf_pow2poly(tmp_power_result[j]) ^ tmp_add_result;
			//if(0xd == i)printf("tmp_add_result: %d = %d + %d\n", gf_add(tmp_power_result[j], tmp_add_result), tmp_power_result[j], tmp_add_result);
			tmp_add_result = gf_add(tmp_power_result[j], tmp_add_result);
		}

		//codeword[i] = gf_poly2pow(tmp_add_result);
		codeword[i] = tmp_add_result;
	}

#if 1
	memcpy(encoded_polynomial, codeword, sizeof(unsigned char) * CODEWORD_LEN);
#else
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		encoded_polynomial[i] = codeword[CODEWORD_LEN - 1 - i];
	}
#endif

#if (1 == CFG_DEBUG_IMPOTANT)
	DEBUG_IMPOTANT("Evaluation Encoding Codeword:\n");
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		DEBUG_IMPOTANT("%x ", encoded_polynomial[i]);
	}
	DEBUG_IMPOTANT("\n");
#endif

	return 0;
}

unsigned char systematic_encoding()
{
	unsigned long long i = 0;
	
	unsigned long long clock_cycle = 0;
	unsigned char input_message = 0;
	unsigned char feed_back_reg = 0, reg_delay[CODEWORD_LEN - MESSAGE_LEN];
	unsigned char codeword[CODEWORD_LEN];

	/*for iterating calculation, tmp for C, not in HDL*/
	unsigned char feed_back_reg_prev = feed_back_reg, reg_delay_prev[CODEWORD_LEN - MESSAGE_LEN];

	memset(reg_delay, 0xFF, sizeof(unsigned char) * (CODEWORD_LEN - MESSAGE_LEN));
	memset(codeword, 0xFF, sizeof(unsigned char) * CODEWORD_LEN);
	memcpy(reg_delay_prev, reg_delay, sizeof(unsigned char) * (CODEWORD_LEN - MESSAGE_LEN));

	if(0 == clock_cycle)
	{
		input_message = message_polynomial[MESSAGE_LEN - 1 - clock_cycle];
		feed_back_reg = gf_add(input_message, reg_delay[CODEWORD_LEN - MESSAGE_LEN - 1]);
		feed_back_reg_prev = feed_back_reg;
		memcpy(reg_delay_prev, reg_delay, sizeof(unsigned char) * (CODEWORD_LEN - MESSAGE_LEN));
	}
	for(clock_cycle = 1; clock_cycle < (CODEWORD_LEN - MESSAGE_LEN); clock_cycle++)
	{
		for(i = CODEWORD_LEN - MESSAGE_LEN - 1; i > 0; i--)
		{
			reg_delay[i] = gf_add(gf_multp(feed_back_reg_prev, generator_polynomial[i]), reg_delay_prev[i - 1]);
		}
		reg_delay[0] = gf_multp(feed_back_reg_prev, generator_polynomial[0]);
		input_message = message_polynomial[MESSAGE_LEN - 1 - clock_cycle];		
		feed_back_reg = gf_add(input_message, reg_delay[CODEWORD_LEN - MESSAGE_LEN - 1]);

		feed_back_reg_prev = feed_back_reg;
		memcpy(reg_delay_prev, reg_delay, sizeof(unsigned char) * (CODEWORD_LEN - MESSAGE_LEN));
	}

	for(i = CODEWORD_LEN - MESSAGE_LEN; i < CODEWORD_LEN; i++)
	{
		codeword[i] = message_polynomial[i - (CODEWORD_LEN - MESSAGE_LEN)];
	}
	for(i = 0; i < CODEWORD_LEN - MESSAGE_LEN; i++)
	{
		codeword[i] = reg_delay[i];
	}
	memcpy(encoded_polynomial, codeword, sizeof(unsigned char) * CODEWORD_LEN);

#if (1 == CFG_DEBUG_IMPOTANT)
	DEBUG_IMPOTANT("Systematic Encoding Codeword:\n");
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		DEBUG_IMPOTANT("%x ", encoded_polynomial[i]);
	}
	DEBUG_IMPOTANT("\n");
#endif

	return 0;
}

unsigned char evaluation_encoding_v2(unsigned char *message,
											unsigned char *codeword_output)
{
	unsigned long long i = 0, j = 0;
	unsigned char tmp_power_result[MESSAGE_LEN], tmp_add_result = 0;

	unsigned char codeword[CODEWORD_LEN];

	for(i = 0; i < CODEWORD_LEN; i++)
	{
		for(j = 0; j < MESSAGE_LEN; j++)
		{
			tmp_power_result[j] = gf_multp(message[j], i * gf_location(j));
		}
		tmp_add_result = 0;
		for(j = 0; j < MESSAGE_LEN; j++)
		{
			tmp_add_result = gf_pow2poly(tmp_power_result[j]) ^ tmp_add_result;
		}

		for(j = 0; j < GF_FIELD; j++)
		{
			codeword[i] = gf_poly2pow(tmp_add_result);
		}
	}

#if 1
	memcpy(codeword_output, codeword, sizeof(unsigned char) * CODEWORD_LEN);
#else
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		encoded_polynomial[i] = codeword[CODEWORD_LEN - 1 - i];
	}
#endif
#if 0
	DEBUG_INFO("Evaluation Encoding Codeword:\n");
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		DEBUG_INFO("%x ", codeword_output[i]);
	}
	DEBUG_INFO("\n");
#endif
	return 0;
}

unsigned char systematic_encoding_v2(unsigned char *message,
											unsigned char *codeword_output)
{
	unsigned long long i = 0;
	
	unsigned long long clock_cycle = 0;
	unsigned char input_message = 0;
	unsigned char feed_back_reg = 0, reg_delay[CODEWORD_LEN - MESSAGE_LEN];
	unsigned char codeword[CODEWORD_LEN];

	/*for iterating calculation, tmp for C, not in HDL*/
	unsigned char feed_back_reg_prev = feed_back_reg, reg_delay_prev[CODEWORD_LEN - MESSAGE_LEN];

	memset(reg_delay, 0xFF, sizeof(unsigned char) * (CODEWORD_LEN - MESSAGE_LEN));
	memset(codeword, 0xFF, sizeof(unsigned char) * CODEWORD_LEN);
	memcpy(reg_delay_prev, reg_delay, sizeof(unsigned char) * (CODEWORD_LEN - MESSAGE_LEN));

	if(0 == clock_cycle)
	{
		input_message = message[MESSAGE_LEN - 1 - clock_cycle];
		feed_back_reg = gf_add(input_message, reg_delay[CODEWORD_LEN - MESSAGE_LEN - 1]);
		feed_back_reg_prev = feed_back_reg;
		memcpy(reg_delay_prev, reg_delay, sizeof(unsigned char) * (CODEWORD_LEN - MESSAGE_LEN));
	}
	for(clock_cycle = 1; clock_cycle < (CODEWORD_LEN - MESSAGE_LEN); clock_cycle++)
	{
		for(i = CODEWORD_LEN - MESSAGE_LEN - 1; i >0; i--)
		{
			reg_delay[i] = gf_add(gf_multp(feed_back_reg_prev, generator_polynomial[i]), reg_delay_prev[i - 1]);
		}
		reg_delay[0] = gf_multp(feed_back_reg_prev, generator_polynomial[0]);
		input_message = message[MESSAGE_LEN - 1 - clock_cycle];		
		feed_back_reg = gf_add(input_message, reg_delay[CODEWORD_LEN - MESSAGE_LEN - 1]);

		feed_back_reg_prev = feed_back_reg;
		memcpy(reg_delay_prev, reg_delay, sizeof(unsigned char) * (CODEWORD_LEN - MESSAGE_LEN));
	}

	for(i = CODEWORD_LEN - MESSAGE_LEN; i < CODEWORD_LEN; i++)
	{
		codeword[i] = message[i - (CODEWORD_LEN - MESSAGE_LEN)];
	}
	for(i = 0; i < CODEWORD_LEN - MESSAGE_LEN; i++)
	{
		codeword[i] = reg_delay[i];
	}
	memcpy(codeword_output, codeword, sizeof(unsigned char) * CODEWORD_LEN);

#if 0
	DEBUG_INFO("Systematic Encoding Codeword:\n");
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		DEBUG_INFO("%x ", codeword_output[i]);
	}
	DEBUG_INFO("\n");
#endif	

	return 0;
}

int gen_poly_trans()
{
	unsigned long long i = 0, j = 0;
	
	unsigned char tmp_generator_polynomial[CODEWORD_LEN - MESSAGE_LEN + 1];
	memset(tmp_generator_polynomial, 0xFF, sizeof(unsigned char) * (CODEWORD_LEN - MESSAGE_LEN + 1));

	for(i = 0; i < (CODEWORD_LEN - MESSAGE_LEN + 1); i++)
	{
		for(j = 0; j < GF_FIELD; j++)
		{
			if(generator_polynomial[i] == power_polynomial_table[j][1])
			{
				tmp_generator_polynomial[(CODEWORD_LEN - MESSAGE_LEN + 1) - 1 - i] = power_polynomial_table[j][0];
				DEBUG_NOTICE("gen-%d: %x\n", i, tmp_generator_polynomial[(CODEWORD_LEN - MESSAGE_LEN + 1) - 1 - i]);
				break;
			}
		}
	}

	memcpy(generator_polynomial, tmp_generator_polynomial, sizeof(unsigned char) * (CODEWORD_LEN - MESSAGE_LEN + 1));

	for(i = 0; i < (CODEWORD_LEN - MESSAGE_LEN + 1); i++)
	{
		printf("gen-%d: %x\n", i, generator_polynomial[i]);
	}

	return 0;
}
