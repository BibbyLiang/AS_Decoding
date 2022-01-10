#include <stdlib.h>
#include <string.h>
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
	unsigned long long hamm_err = 0;

	clock_t start, stop;
	float runtime;

	float eb2n0_start = 15, eb2n0_stop = 15, eb2n0_step = 1, eb2n0 = 15;
	unsigned long iter_cnt = 0, monitor_cnt = 1;
#if (0 == TEST_MODE)
	printf("Please Input Eb/N0 Start: ");
	scanf("%f", &eb2n0_start);
	printf("Please Input Eb/N0 Stop: ");
	scanf("%f", &eb2n0_stop);
	printf("Please Input Eb/N0 Step: ");
	scanf("%f", &eb2n0_step);
	printf("Please Input Simulation Times: ");
	scanf("%ld", &iter_cnt);
	printf("Please Monitor Times: ");
	scanf("%ld", &monitor_cnt);
#endif

	char log_name[255];
	sprintf(log_name, "n_%d-k_%d-m_%d-sasd_%d-snr_%f_%f_%f-cnt_%ld_%ld.txt",
					  CODEWORD_LEN,
					  MESSAGE_LEN,
					  S_MUL,
					  SIMPLE_ASD,
					  eb2n0_start,
					  eb2n0_step,
					  eb2n0_stop,
					  iter_cnt,
					  monitor_cnt);
	FILE *frc;

#if 0
	float **mod_seq;
	mod_seq = (float**)malloc(sizeof(float*) * symbol_num);
	for (i = 0; i < symbol_num; i++)
	{
		mod_seq[i] = (float*)malloc(sizeof(float) * 2);
		//DEBUG_NOTICE("%d %d\n", symbol_num, i);
	}
#endif
	g_term_malloc();

	mod_init();

	init_simulation();

#if (1 == SYS_ENC)
	gen_poly_trans();
#endif

	start = clock();

	for(eb2n0 = eb2n0_start; eb2n0 <= eb2n0_stop; eb2n0 = eb2n0 + eb2n0_step)
	{
		bit_err = 0;
		symbol_err = 0;
		frame_err = 0;
		uncoded_bit_err = 0;
		uncoded_symbol_err = 0;
		uncoded_frame_err = 0;
		hamm_err = 0;

		for(iter = 0; iter <= iter_cnt; iter++)
		{

			decoding_ok_flag = 0;
			err_num = 0;
			memset(recv_rel, 0.0, sizeof(float) * CODEWORD_LEN);

			for(i = 0; i < MESSAGE_LEN; i++)
			{
				j = genrand_int32() % GF_FIELD;
#if (1 == TEST_MODE)
				message_polynomial[i] = 0x0;
#else
				message_polynomial[i] = power_polynomial_table[j][0];
#endif			
			}

#if (1 == TEST_MODE)//test
			message_polynomial[0] = 0x4;
			message_polynomial[1] = 0x3;
			message_polynomial[2] = 0x5;
#endif

#if (1 == SYS_ENC)
			systematic_encoding();
#else
			evaluation_encoding();
#endif

#if 0
			bpsk_mod(encoded_polynomial,
					 CODEWORD_LEN,
					 (float **)mod_seq,
					 symbol_num);

			DEBUG_IMPOTANT("Transmission over Channel:\n");
			for(i = 0; i < symbol_num; i++)
			{
				mod_seq[i][0] = mod_seq[i][0] + awgn_gen(eb2n0);
				mod_seq[i][1] = mod_seq[i][1] + awgn_gen(eb2n0);
				DEBUG_IMPOTANT("%f %f\n", mod_seq[i][0], mod_seq[i][1]);

				recv_seq[i][0] = mod_seq[i][0];
				recv_seq[i][1] = mod_seq[i][1];
			}
			DEBUG_IMPOTANT("\n");

			bpsk_demod((float **)mod_seq,
					 symbol_num,
					 received_polynomial,
					 CODEWORD_LEN);
#else
			bpsk_mod(encoded_polynomial,
					 CODEWORD_LEN,
					 recv_seq,
					 symbol_num);

			DEBUG_IMPOTANT("Transmission over Channel:\n");
			for(i = 0; i < symbol_num; i++)
			{
				recv_seq[i][0] = recv_seq[i][0] + awgn_gen(eb2n0);
				recv_seq[i][1] = recv_seq[i][1] + awgn_gen(eb2n0);
				DEBUG_NOTICE("%f %f\n", recv_seq[i][0], recv_seq[i][1]);
			}
			DEBUG_IMPOTANT("\n");

			bpsk_demod((float **)recv_seq,
					 symbol_num,
					 received_polynomial,
					 CODEWORD_LEN);
#endif

			for(i = 0; i < CODEWORD_LEN; i++)
			{
				DEBUG_NOTICE("%x %x\n", received_polynomial[i], encoded_polynomial[i]);
				if(received_polynomial[i] != encoded_polynomial[i])
				{
					uncoded_symbol_err = uncoded_symbol_err + 1;
					err_num = err_num+ 1;
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

#if (1 == TEST_MODE)//test
			/*transmission through channel*/
			for(i = 0; i < CODEWORD_LEN; i++)
			{
				received_polynomial[i] = gf_add(encoded_polynomial[i], error_polynomial[i]);
			}
#endif

#if 1
			chnl_rel_cal();

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
			if((2 == decoding_ok_flag)
				&& (1 == frame_err_flag))
			{
				if(0 == hamm_distance_debug)
				{
#if 0					
					DEBUG_SYS("hamm_distance_cal_err\n");
					frc = fopen(log_name, "a+");
					fprintf(frc, "hamm_distance_cal_err\n");
					fclose(frc);
					frc = NULL;
#endif
					hamm_err = hamm_err + 1;
				}
				else
				{
					DEBUG_SYS("Prog. Err. for Decoding\n");
					frc = fopen(log_name, "a+");
					fprintf(frc, "Prog. Err. for Decoding\n");

					DEBUG_SYS("Para.: %ld %ld %ld\n",
							  (S_MUL * (CODEWORD_LEN - err_num)),
							  weight_stored,
							  hamm_distance_debug);
					fprintf(frc, "Para.: %ld %ld %ld\n",
							(S_MUL * (CODEWORD_LEN - err_num)),
							weight_stored,
							hamm_distance_debug);
					
					for(i = 0; i < MESSAGE_LEN; i++)
					{
						DEBUG_SYS("message: %x\n", message_polynomial[i]);
						fprintf(frc, "message: %x\n", message_polynomial[i]);
					}
					for(i = 0; i < CODEWORD_LEN; i++)
					{
						DEBUG_SYS("encoded: %x\n", encoded_polynomial[i]);
						fprintf(frc, "encoded: %x\n", encoded_polynomial[i]);
					}
					for(i = 0; i < CODEWORD_LEN; i++)
					{
						DEBUG_SYS("recv: %x\n", received_polynomial[i]);
						fprintf(frc, "recv: %x\n", received_polynomial[i]);
					}
					for(i = 0; i < CODEWORD_LEN; i++)
					{
						DEBUG_SYS("err: %x\n", gf_add(received_polynomial[i], encoded_polynomial[i]));
						fprintf(frc, "err: %x\n", gf_add(received_polynomial[i], encoded_polynomial[i]));
					}

					fclose(frc);
					frc = NULL;
				}
			}
			frame_err_flag = 0;

			if(0 == (iter % monitor_cnt))
			{
				stop = clock();
				runtime = (stop - start) / 1000.0000;
				
				DEBUG_SYS("---------------------\n");
				DEBUG_SYS("Time: %fs\n", runtime);
				DEBUG_SYS("Eb/N0: %f dB\n", eb2n0);
				DEBUG_SYS("Frame: %ld\n", iter);
				DEBUG_SYS("Uncoded Frame Error: %ld\n", uncoded_frame_err);
				DEBUG_SYS("Uncoded Symbol Error: %ld\n", uncoded_symbol_err);
				DEBUG_SYS("Uncoded Bit Error: %ld\n", uncoded_bit_err);
				DEBUG_SYS("Frame Error: %ld\n", frame_err);
				DEBUG_SYS("Symbol Error: %ld\n", symbol_err);
				DEBUG_SYS("Bit Error: %ld\n", bit_err);
				DEBUG_SYS("Hamming Error: %ld\n", hamm_err);

				frc = fopen(log_name, "a+");
				fprintf(frc, "---------------------\n");
				fprintf(frc, "Time: %fs\n", runtime);
				fprintf(frc, "Eb/N0: %f dB\n", eb2n0);
				fprintf(frc, "Frame: %ld\n", iter);
				fprintf(frc, "Uncoded Frame Error: %ld\n", uncoded_frame_err);
				fprintf(frc, "Uncoded Symbol Error: %ld\n", uncoded_symbol_err);
				fprintf(frc, "Uncoded Bit Error: %ld\n", uncoded_bit_err);
				fprintf(frc, "Frame Error: %ld\n", frame_err);
				fprintf(frc, "Symbol Error: %ld\n", symbol_err);
				fprintf(frc, "Bit Error: %ld\n", bit_err);
				fprintf(frc, "Hamming Error: %ld\n", hamm_err);
			    fclose(frc);
				frc = NULL;
			}

		}

		stop = clock();
		runtime = (stop - start) / 1000.0000;
		iter = iter - 1;

		DEBUG_SYS("*********************************\n");
		DEBUG_SYS("Time: %fs\n", runtime);
		DEBUG_SYS("Eb/N0: %f dB\n", eb2n0);
		DEBUG_SYS("Frame: %ld\n", iter_cnt);
		DEBUG_SYS("Uncoded Frame Error: %ld\n", uncoded_frame_err);
		DEBUG_SYS("Uncoded Symbol Error: %ld\n", uncoded_symbol_err);
		DEBUG_SYS("Uncoded Bit Error: %ld\n", uncoded_bit_err);
		DEBUG_SYS("Frame Error: %ld\n", frame_err);
		DEBUG_SYS("Symbol Error: %ld\n", symbol_err);
		DEBUG_SYS("Bit Error: %ld\n", bit_err);
		DEBUG_SYS("Hamming Error: %ld\n", hamm_err);
		frc = fopen(log_name, "a+");
		fprintf(frc, "*********************************\n");
		fprintf(frc, "Time: %fs\n", runtime);
		fprintf(frc, "Eb/N0: %f dB\n", eb2n0);
		fprintf(frc, "Frame: %ld\n", iter);
		fprintf(frc, "Uncoded Frame Error: %ld\n", uncoded_frame_err);
		fprintf(frc, "Uncoded Symbol Error: %ld\n", uncoded_symbol_err);
		fprintf(frc, "Uncoded Bit Error: %ld\n", uncoded_bit_err);
		fprintf(frc, "Frame Error: %ld\n", frame_err);
		fprintf(frc, "Symbol Error: %ld\n", symbol_err);
		fprintf(frc, "Bit Error: %ld\n", bit_err);
		fprintf(frc, "Hamming Error: %ld\n", hamm_err);
		fprintf(frc, "Uncoded Results: %.10lf %.10lf %.10lf\n", 
			    (double)uncoded_frame_err / (double)iter,
			    (double)uncoded_symbol_err / (double)iter / CODEWORD_LEN * BITS_PER_SYMBOL_BPSK,
			    (double)uncoded_bit_err / (double)iter / CODEWORD_LEN * BITS_PER_SYMBOL_BPSK / GF_Q);
		fprintf(frc, "Decoded Results: %.10lf %.10lf %.10lf %.10lf\n", 
				(double)frame_err / (double)iter,
				(double)symbol_err / (double)iter / CODEWORD_LEN * BITS_PER_SYMBOL_BPSK,
				(double)bit_err / (double)iter / CODEWORD_LEN * BITS_PER_SYMBOL_BPSK / GF_Q,
				(double)(frame_err - hamm_err) / (double)iter);
	    fclose(frc);
		frc = NULL;

	}

#if 0
	for (i = 0; i < symbol_num; i++)
	{
  		free(mod_seq[i]);
		mod_seq[i] = NULL;
  	}
	free(mod_seq);
	mod_seq = NULL;
#endif
	mod_exit();

	g_term_destroy();

	return;
}
