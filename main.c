#include <stdlib.h>
#include "debug_info.h"
#include "gf_cal.h"
#include "as_decoding.h"
#include "encoding.h"
#include "mod.h"
#include "rnd.h"
#include "channel.h"
#include "time.h"

#define SIMULATION_TIMES	1000

void init_simulation()
{
	srand(time(NULL));
	init_genrand((long)(time(NULL)));

	return;
}

void main()
{
	unsigned char i = 0, j = 0;
	unsigned int symbol_num = CODEWORD_LEN * GF_Q * BITS_PER_SYMBOL_BPSK;
	unsigned long iter = 0;
	unsigned long bit_err = 0, symbol_err = 0, frame_err = 0;
	unsigned long uncoded_bit_err = 0, uncoded_symbol_err = 0, uncoded_frame_err = 0;
	unsigned char frame_err_flag = 0, uncoded_frame_err_flag = 0;
	unsigned char tmp_mes = 0, tmp_dec = 0;

	float snr = 7;

	float **mod_seq;
	mod_seq = (float**)malloc(sizeof(float*) * symbol_num);
	for (i = 0; i < symbol_num; i++)
	{
		mod_seq[i] = (float*)malloc(sizeof(float) * 2);
	}

	init_simulation();

	for(iter = 0; iter < SIMULATION_TIMES; iter++)
	{

#if (1 == SYS_ENC)
		systematic_encoding();
#else
		evaluation_encoding();
#endif

		bpsk_mod(encoded_polynomial,
				 CODEWORD_LEN,
				 (float **)mod_seq,
				 symbol_num);

#if 1
		DEBUG_IMPOTANT("Transmission over Channel:\n");
		for(i = 0; i < symbol_num; i++)
		{
			mod_seq[i][0] = mod_seq[i][0] + awgn_gen(snr);
			//mod_seq[i][1] = mod_seq[i][1] + awgn_gen(snr);
			DEBUG_IMPOTANT("%f %f\n", mod_seq[i][0], mod_seq[i][1]);
		}
		DEBUG_IMPOTANT("\n");
#endif

		bpsk_demod((float **)mod_seq,
				 symbol_num,
				 received_polynomial,
				 CODEWORD_LEN);

		for(i = 0; i < CODEWORD_LEN; i++)
		{
			DEBUG_NOTICE("%x %x\n", received_polynomial[i], encoded_polynomial[i]);
			if(received_polynomial[i] != encoded_polynomial[i])
			{
				uncoded_symbol_err = uncoded_symbol_err + 1;
				if(0 == uncoded_frame_err_flag)
				{
					uncoded_frame_err_flag = 1;
					uncoded_frame_err = uncoded_frame_err + 1;
				}
				for(j = 0; j < GF_Q; j++)
				{
					tmp_mes = ((received_polynomial[i] >> j) & 0x1);
					tmp_dec = ((encoded_polynomial[i] >> j) & 0x1);
					if(tmp_mes != tmp_dec)
					{
						uncoded_bit_err = uncoded_bit_err + 1;
					}
				}
				
			}
		}
		uncoded_frame_err_flag = 0;

#if 0
		/*transmission through channel*/
		for(i = 0; i < CODEWORD_LEN; i++)
		{
			received_polynomial[i] = gf_add(encoded_polynomial[i], error_polynomial[i]);
		}
#endif

		mul_assign();
		
		re_encoding();
		
		as_decoding();

		for(i = 0; i < MESSAGE_LEN; i++)
		{
			DEBUG_NOTICE("%x %x\n", decoded_message[i], message_polynomial[i]);
			if(decoded_message[i] != message_polynomial[i])
			{
				symbol_err = symbol_err + 1;
				if(0 == frame_err_flag)
				{
					frame_err_flag = 1;
					frame_err = frame_err + 1;
				}
				for(j = 0; j < GF_Q; j++)
				{
					tmp_mes = ((message_polynomial[i] >> j) & 0x1);
					tmp_dec = ((decoded_message[i] >> j) & 0x1);
					if(tmp_mes != tmp_dec)
					{
						bit_err = bit_err + 1;
					}
				}
				
			}
		}
		frame_err_flag = 0;

	}

	DEBUG_SYS("Frame: %ld\n", SIMULATION_TIMES);
	DEBUG_SYS("Uncoded Frame Error: %ld\n", uncoded_frame_err);
	DEBUG_SYS("Uncoded Symbol Error: %ld\n", uncoded_symbol_err);
	DEBUG_SYS("Uncoded Bit Error: %ld\n", uncoded_bit_err);
	DEBUG_SYS("Frame Error: %ld\n", frame_err);
	DEBUG_SYS("Symbol Error: %ld\n", symbol_err);
	DEBUG_SYS("Bit Error: %ld\n", bit_err);

	for (i = 0; i < symbol_num; i++)
	{
  		free(mod_seq[i]);
		mod_seq[i] = NULL;
  	}
	free(mod_seq);
	mod_seq = NULL;

	return;
}
