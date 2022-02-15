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
#if (2 == MESSAGE_LEN)
	1,
	4,
	3,
	5,
	6,
	2
#endif
	
#if (3 == MESSAGE_LEN)	
	1,
	3,
	1,
	2,
	3
#endif

#if (5 == MESSAGE_LEN)
	1,
	6,
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
#if (15 == MESSAGE_LEN)
	1,
	8,
	48,
	35,
	1,
	55,
	17,
	7,
	20,
	18,
	34,
	57,
	43,
	60,
	63,
	15,
	38,
	42,
	40,
	17,
	25,
	4,
	46,
	6,
	27,
	31,
	14,
	41,
	21,
	38,
	62,
	15,
	37,
	15,
	45,
	46,
	24,
	16,
	15,
	5,
	58,
	47,
	17,
	56,
	11,
	9,
	38,
	61,
	58
#endif

#if (31 == MESSAGE_LEN)
	1,
	49,
	53,
	54,
	19,
	1,
	41,
	35,
	25,
	10,
	47,
	10,
	4,
	54,
	62,
	10,
	55,
	50,
	57,
	5,
	12,
	21,
	26,
	46,
	62,
	31,
	27,
	13,
	36,
	14,
	55,
	21,
	17
#endif

#if (21 == MESSAGE_LEN)//k = 21
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
#endif	
};
#endif

#if (8 == GF_Q)
/*n = 255, k = 239*/
unsigned char generator_polynomial[CODEWORD_LEN - MESSAGE_LEN + 1] =
{
	/*generated by MATLAB: genpoly = rsgenpolycoeffs(n, k), in add-style*/
	1,
	118,
	52,
	103,
	31,
	104,
	126,
	187,
	232,
	17,
	56,
	183,
	49,
	100,
	81,
	44,
	79
};
#endif

unsigned char message_polynomial[MESSAGE_LEN] = 
{
	/*power representation*/
	//0x1,
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
	0x0,
	0x0,
	0xFF,
	0xFF,
	0xFF,
	0xFF,
	0xFF
};

unsigned char evaluation_encoding()
{
	long long i = 0, j = 0;
	unsigned char tmp_power_result[MESSAGE_LEN], tmp_add_result = 0;

	unsigned char codeword[CODEWORD_LEN];

	for(i = 0; i < CODEWORD_LEN; i++)
	{
		for(j = 0; j < MESSAGE_LEN; j++)
		{
			tmp_power_result[j] = gf_multp(message_polynomial[j], (i * gf_location(j)) % (GF_FIELD - 1));
			//DEBUG_NOTICE("tmp_power_result: %d = %d * %d\n", tmp_power_result[j], message_polynomial[j], i * gf_location(j));
		}

		tmp_add_result = 0xFF;
		for(j = 0; j < MESSAGE_LEN; j++)
		{
			//tmp_add_result = gf_pow2poly(tmp_power_result[j]) ^ tmp_add_result;
			//DEBUG_NOTICE("tmp_add_result: %d = %d + %d\n", gf_add(tmp_power_result[j], tmp_add_result), tmp_power_result[j], tmp_add_result);
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
	long long i = 0;
	
	long long clock_cycle = 0;
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
			reg_delay[i] = gf_add(gf_multp(feed_back_reg_prev, (generator_polynomial[i])), reg_delay_prev[i - 1]);
		}
		reg_delay[0] = gf_multp(feed_back_reg_prev, (generator_polynomial[0]));
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
	DEBUG_IMPOTANT("Systematic Encoding Codeword = {\n");
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		//DEBUG_IMPOTANT("%x ", encoded_polynomial[i]);
		if(0xFF == encoded_polynomial[i])
		{
			DEBUG_IMPOTANT("%d, ", power_polynomial_table[0][1]);
		}
		else
		{
			DEBUG_IMPOTANT("%d, ", power_polynomial_table[encoded_polynomial[i] + 1][1]);
		}
	}
	DEBUG_IMPOTANT("};\n");
#endif

	return 0;
}

unsigned char evaluation_encoding_v2(unsigned char *message,
											unsigned char *codeword_output)
{
	long long i = 0, j = 0;
	unsigned char tmp_power_result[MESSAGE_LEN], tmp_add_result = 0;

	unsigned char codeword[CODEWORD_LEN];

	for(i = 0; i < CODEWORD_LEN; i++)
	{
		for(j = 0; j < MESSAGE_LEN; j++)
		{
			tmp_power_result[j] = gf_multp(message[j], (i * gf_location(j)) % (GF_FIELD - 1));
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
#if 1
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
	long long i = 0;
	
	long long clock_cycle = 0;
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
			reg_delay[i] = gf_add(gf_multp(feed_back_reg_prev, (generator_polynomial[i])), reg_delay_prev[i - 1]);
		}
		reg_delay[0] = gf_multp(feed_back_reg_prev, (generator_polynomial[0]));
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
	long long i = 0, j = 0;
	
	unsigned char tmp_generator_polynomial[CODEWORD_LEN - MESSAGE_LEN + 1];
	memset(tmp_generator_polynomial, 0xFF, sizeof(unsigned char) * (CODEWORD_LEN - MESSAGE_LEN + 1));

	for(i = 0; i < (CODEWORD_LEN - MESSAGE_LEN + 1); i++)
	{
		for(j = 0; j < GF_FIELD; j++)
		{
			if(generator_polynomial[i] == power_polynomial_table[j][1])
			{
				tmp_generator_polynomial[(CODEWORD_LEN - MESSAGE_LEN + 1) - 1 - i] = power_polynomial_table[j][0];
				DEBUG_NOTICE("gen-%d: %d\n", i, tmp_generator_polynomial[(CODEWORD_LEN - MESSAGE_LEN + 1) - 1 - i]);
				break;
			}
		}
	}

	memcpy(generator_polynomial, tmp_generator_polynomial, sizeof(unsigned char) * (CODEWORD_LEN - MESSAGE_LEN + 1));

#if 0
	for(i = 0; i < (CODEWORD_LEN - MESSAGE_LEN + 1); i++)
	{
		DEBUG_NOTICE("gen-%d: %x\n", i, generator_polynomial[i]);
	}
#endif

	DEBUG_SYS("gen_poly_trans OK\n");

	return 0;
}

void encoding2()
{
    int i, j, k;

	unsigned char msg[MESSAGE_LEN], tmp = 0xFF;
	memcpy(msg, message_polynomial, sizeof(unsigned char) * MESSAGE_LEN);
	unsigned char cwd[CODEWORD_LEN];
	memset(cwd, 0xFF, sizeof(unsigned char) * CODEWORD_LEN);
	
	unsigned char Quotient[MESSAGE_LEN] = { 0 };
    unsigned char gx1[CODEWORD_LEN - MESSAGE_LEN + 1];
    for (i = 0; i < (CODEWORD_LEN - MESSAGE_LEN + 1); i++)//n-k+1 denote the numbers of term (generator polynomial)
    {
        gx1[i] = gf_div(0x0, generator_polynomial[i]);
		//printf("gx1_tmp: %d\n", gx1[i]);
		//gx1[i] = gf_pow2poly(gx1[i]);
		//printf("gx1: %d %d %d %d\n", generator_polynomial[i], tmp, i, gx1[i]);
    }
    for (i = 0; i < MESSAGE_LEN; i++)
        cwd[i] = msg[i];
    for (i = MESSAGE_LEN; i < CODEWORD_LEN; i++)
        cwd[i] = 0xFF;
    for (j = 0; j < MESSAGE_LEN; j++)//find the remainder (long division)
    {
        Quotient[j] = gf_multp(cwd[j], gx1[0]);
		//DEBUG_NOTICE("Quotient: %d = %d * %d\n", gf_pow2poly(Quotient[j]), gf_pow2poly(cwd[j]), gf_pow2poly(gx1[0]));
        for (i = j; i < (CODEWORD_LEN - MESSAGE_LEN + 1 + j); i++)
        {
			tmp = cwd[i];
            cwd[i] = gf_add(cwd[i], gf_multp(Quotient[j], (generator_polynomial[i - j])));
			//DEBUG_NOTICE("cwd[i]: %d %d %d %d\n", gf_pow2poly(cwd[i]), gf_pow2poly(tmp), gf_pow2poly(Quotient[j]), gf_pow2poly(generator_polynomial[i - j]));
        }
    }

    for (i = 0; i < (CODEWORD_LEN - MESSAGE_LEN); i++)
    {
        cwd[i] = cwd[MESSAGE_LEN + i];
    }
    for (i = (CODEWORD_LEN - MESSAGE_LEN); i < CODEWORD_LEN; i++)
        cwd[i] = msg[i - CODEWORD_LEN + MESSAGE_LEN];

	for(i = 0; i < MESSAGE_LEN; i++)
	{
		DEBUG_NOTICE("msg: %d | %d %d\n", i, message_polynomial[i], gf_pow2poly(message_polynomial[i]));
	}
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		DEBUG_NOTICE("cwd: %d | %d %d\n", i, cwd[i], gf_pow2poly(cwd[i]));
	}

	memcpy(encoded_polynomial, cwd, sizeof(unsigned char) * CODEWORD_LEN);
}
