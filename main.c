#include <stdlib.h>
#include "debug_info.h"
#include "gf_cal.h"
#include "as_decoding.h"
#include "encoding.h"
#include "mod.h"
#include "rnd.h"
#include "channel.h"
#include "time.h"

void init_simulation()
{
	srand(time(NULL));
	init_genrand((long)(time(NULL)));

	return;
}

void main()
{
	unsigned long long i = 0, j = 0;
	unsigned long long symbol_num = CODEWORD_LEN * GF_Q * BITS_PER_SYMBOL_BPSK;
	unsigned long long iter = 0;
	unsigned long long bit_err = 0, symbol_err = 0, frame_err = 0;
	unsigned long long uncoded_bit_err = 0, uncoded_symbol_err = 0, uncoded_frame_err = 0;
	unsigned char frame_err_flag = 0, uncoded_frame_err_flag = 0;
	unsigned char tmp_mes = 0, tmp_dec = 0;

	clock_t start, stop;
	float runtime;

	float snr_start = 15, snr_stop = 15, snr_step = 1, snr = 15;
	unsigned long iter_cnt = 0, monitor_cnt = 1;
#if (0 == TEST_MODE)
	printf("Please Input SNR Start: ");
	scanf("%f", &snr_start);
	printf("Please Input SNR Stop: ");
	scanf("%f", &snr_stop);
	printf("Please Input SNR Step: ");
	scanf("%f", &snr_step);
	printf("Please Input Simulation Times: ");
	scanf("%ld", &iter_cnt);
	printf("Please Monitor Times: ");
	scanf("%ld", &monitor_cnt);
#endif
	FILE *frc;

	float **mod_seq;
	mod_seq = (float**)malloc(sizeof(float*) * symbol_num);
	for (i = 0; i < symbol_num; i++)
	{
		mod_seq[i] = (float*)malloc(sizeof(float) * 2);
		//DEBUG_NOTICE("%d %d\n", symbol_num, i);
	}

	init_simulation();

#if (1 == SYS_ENC)
	gen_poly_trans();
#endif

	start = clock();

	for(snr = snr_start; snr <= snr_stop; snr = snr + snr_step)
	{
		bit_err = 0;
		symbol_err = 0;
		frame_err = 0;
		uncoded_bit_err = 0;
		uncoded_symbol_err = 0;
		uncoded_frame_err = 0;

		for(iter = 0; iter <= iter_cnt; iter++)
		{

			for(i = 0; i < MESSAGE_LEN; i++)
			{
				j = genrand_int32() % GF_FIELD;
#if (1 == TEST_MODE)			
				message_polynomial[i] = 0x0;
#else
				message_polynomial[i] = power_polynomial_table[j][0];
#endif			
			}

#if (1 == SYS_ENC)
			systematic_encoding();
#else
			evaluation_encoding();
#endif

			bpsk_mod(encoded_polynomial,
					 CODEWORD_LEN,
					 (float **)mod_seq,
					 symbol_num);

			DEBUG_IMPOTANT("Transmission over Channel:\n");
			for(i = 0; i < symbol_num; i++)
			{
				mod_seq[i][0] = mod_seq[i][0] + awgn_gen(snr);
				mod_seq[i][1] = mod_seq[i][1] + awgn_gen(snr);
				DEBUG_IMPOTANT("%f %f\n", mod_seq[i][0], mod_seq[i][1]);
			}
			DEBUG_IMPOTANT("\n");

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

#if 1
			mul_assign();
			
			re_encoding();
			
			as_decoding();
#endif

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

			if(0 == (iter % monitor_cnt))
			{
				stop = clock();
				runtime = (stop - start) / 1000.0000;
				
				DEBUG_SYS("---------------------\n");
				DEBUG_SYS("Time: %fs\n", runtime);
				DEBUG_SYS("SNR: %f dB\n", snr);
				DEBUG_SYS("Frame: %ld\n", iter);
				DEBUG_SYS("Uncoded Frame Error: %ld\n", uncoded_frame_err);
				DEBUG_SYS("Uncoded Symbol Error: %ld\n", uncoded_symbol_err);
				DEBUG_SYS("Uncoded Bit Error: %ld\n", uncoded_bit_err);
				DEBUG_SYS("Frame Error: %ld\n", frame_err);
				DEBUG_SYS("Symbol Error: %ld\n", symbol_err);
				DEBUG_SYS("Bit Error: %ld\n", bit_err);

				frc = fopen("runing_log.txt", "a+");
				fprintf(frc, "---------------------\n");
				fprintf(frc, "Time: %fs\n", runtime);
				fprintf(frc, "SNR: %f dB\n", snr);
				fprintf(frc, "Frame: %ld\n", iter);
				fprintf(frc, "Uncoded Frame Error: %ld\n", uncoded_frame_err);
				fprintf(frc, "Uncoded Symbol Error: %ld\n", uncoded_symbol_err);
				fprintf(frc, "Uncoded Bit Error: %ld\n", uncoded_bit_err);
				fprintf(frc, "Frame Error: %ld\n", frame_err);
				fprintf(frc, "Symbol Error: %ld\n", symbol_err);
				fprintf(frc, "Bit Error: %ld\n", bit_err);
			    fclose(frc);
				frc = NULL;
			}

		}

		stop = clock();
		runtime = (stop - start) / 1000.0000;

		DEBUG_SYS("*********************************\n");
		DEBUG_SYS("Time: %fs\n", runtime);
		DEBUG_SYS("SNR: %f dB\n", snr);
		DEBUG_SYS("Frame: %ld\n", iter_cnt);
		DEBUG_SYS("Uncoded Frame Error: %ld\n", uncoded_frame_err);
		DEBUG_SYS("Uncoded Symbol Error: %ld\n", uncoded_symbol_err);
		DEBUG_SYS("Uncoded Bit Error: %ld\n", uncoded_bit_err);
		DEBUG_SYS("Frame Error: %ld\n", frame_err);
		DEBUG_SYS("Symbol Error: %ld\n", symbol_err);
		DEBUG_SYS("Bit Error: %ld\n", bit_err);
		frc = fopen("runing_log.txt", "a+");
		fprintf(frc, "*********************************\n");
		fprintf(frc, "Time: %fs\n", runtime);
		fprintf(frc, "SNR: %f dB\n", snr);
		fprintf(frc, "Frame: %ld\n", iter);
		fprintf(frc, "Uncoded Frame Error: %ld\n", uncoded_frame_err);
		fprintf(frc, "Uncoded Symbol Error: %ld\n", uncoded_symbol_err);
		fprintf(frc, "Uncoded Bit Error: %ld\n", uncoded_bit_err);
		fprintf(frc, "Frame Error: %ld\n", frame_err);
		fprintf(frc, "Symbol Error: %ld\n", symbol_err);
		fprintf(frc, "Bit Error: %ld\n", bit_err);
	    fclose(frc);
		frc = NULL;

	}

	for (i = 0; i < symbol_num; i++)
	{
  		free(mod_seq[i]);
		mod_seq[i] = NULL;
  	}
	free(mod_seq);
	mod_seq = NULL;

	return;
}
