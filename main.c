#include <stdlib.h>
#include <string.h>
#include "cfg_decoding.h"
#include "debug_info.h"
#include "gf_cal.h"
#include "as_decoding.h"
#include "re_encoding.h"
#include "encoding.h"
#include "mod.h"
#include "rnd.h"
#include "channel.h"
#include "time.h"

void init_simulation()
{
	srand(time(NULL));
	init_genrand((long)(time(NULL)));

#if (1 == GF_CAL_COUNT)
	memset(err_hist, 0, sizeof(long long) * (CODEWORD_LEN - MESSAGE_LEN - 1));
	memset(add_cnt_hist, 0, sizeof(long long) * (CODEWORD_LEN - MESSAGE_LEN - 1));
	memset(mul_cnt_hist, 0, sizeof(long long) * (CODEWORD_LEN - MESSAGE_LEN - 1));
	memset(div_cnt_hist, 0, sizeof(long long) * (CODEWORD_LEN - MESSAGE_LEN - 1));
	memset(real_cbm_cnt_hist, 0, sizeof(long long) * (CODEWORD_LEN - MESSAGE_LEN - 1));
	memset(real_mul_ff_cnt_hist, 0, sizeof(long long) * (CODEWORD_LEN - MESSAGE_LEN - 1));
	memset(pow_cnt_hist, 0, sizeof(long long) * (CODEWORD_LEN - MESSAGE_LEN - 1));
#endif
	DEBUG_SYS("init_simulation OK\n");

	return;
}

void main()
{
	long long i = 0, j = 0, k = 0;
	long long symbol_num = CODEWORD_LEN * GF_Q * BITS_PER_SYMBOL_BPSK;
	long long iter = 0;
	long long bit_err = 0, symbol_err = 0, symbol_err_prev = 0, frame_err = 0;
	long long uncoded_bit_err = 0, uncoded_symbol_err = 0, uncoded_frame_err = 0;
	unsigned char frame_err_flag = 0, uncoded_frame_err_flag = 0;
	unsigned char tmp_mes = 0, tmp_dec = 0;
	long long hamm_err = 0;
#if (0 == RECUR_RR)	
	long long rr_err_cnt = 0;
#endif
#if (1 == GF_CAL_COUNT)	
	long long symbol_err_this_frame = 0;
#endif

	clock_t start, stop;
	float runtime;

	float eb2n0_start = 10, eb2n0_stop = 10, eb2n0_step = 1, eb2n0 = 10;
	unsigned long iter_cnt = 1, monitor_cnt = 1;
#if (0 == TEST_MODE)
#if 1
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
#endif

#if (1 == OUTPUT_LOG)
	char log_name[255];
	sprintf(log_name, "n_%d-k_%d-m_%d-ret_%d-snr_%f_%f_%f-cnt_%ld_%ld.txt",
					  CODEWORD_LEN,
					  MESSAGE_LEN,
					  S_MUL,
					  RE_ENCODING,
					  eb2n0_start,
					  eb2n0_step,
					  eb2n0_stop,
					  iter_cnt,
					  monitor_cnt);
	FILE *frc;
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
		symbol_err_prev = 0;
		frame_err = 0;
		uncoded_bit_err = 0;
		uncoded_symbol_err = 0;
		uncoded_frame_err = 0;
		hamm_err = 0;

		for(iter = 0; iter < iter_cnt; iter++)
		{
			decoding_ok_flag = 0;
			err_num = 0;
			memset(recv_rel, 0.0, sizeof(float) * CODEWORD_LEN);
			symbol_err_prev = symbol_err;

			for(i = 0; i < MESSAGE_LEN; i++)
			{
				j = genrand_int32() % GF_FIELD;
#if (1 == TEST_MODE)
				message_polynomial[i] = 0x0;
#else
				message_polynomial[i] = power_polynomial_table[j][0];
#endif
			}
			//memset(message_polynomial, 0x0, sizeof(unsigned char) * MESSAGE_LEN);

#if (1 == TEST_MODE)//test
			message_polynomial[0] = 0x0;
			message_polynomial[1] = 0x1;
			message_polynomial[2] = 0x2;
			message_polynomial[3] = 0x3;
			message_polynomial[4] = 0x4;
			//message_polynomial[5] = 0xC;
			//message_polynomial[6] = 0xA;
			//memset(message_polynomial, 0x0, sizeof(unsigned char) * MESSAGE_LEN);
#endif

#if (1 == SYS_ENC)
			//systematic_encoding();
			encoding2();
#else
			evaluation_encoding();
#endif

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

#if (1 == TEST_MODE)//test
			/*transmission through channel*/
#if 0
			for(i = 0; i < CODEWORD_LEN; i++)
			{
				received_polynomial[i] = gf_add(encoded_polynomial[i], error_polynomial[i]);
			}
#else
			memcpy(received_polynomial, encoded_polynomial, sizeof(unsigned char) * CODEWORD_LEN);
			received_polynomial[1] = gf_add(encoded_polynomial[1], 0x2);
			//received_polynomial[3] = gf_add(encoded_polynomial[3], 0x0);
			//received_polynomial[4] = gf_add(encoded_polynomial[4], 0x5);
			//received_polynomial[5] = gf_add(encoded_polynomial[5], 0x0);
			//received_polynomial[9] = gf_add(encoded_polynomial[9], 0x2);
			//received_polynomial[0] = 0x5;
			//received_polynomial[1] = 0x4;
			//received_polynomial[2] = 0x4;
			//received_polynomial[3] = 0x0;
			//received_polynomial[4] = 0x1;
			//received_polynomial[5] = 0x1;
			//received_polynomial[6] = 0x3;
#endif
#endif

			for(i = 0; i < CODEWORD_LEN; i++)
			{
				DEBUG_NOTICE("%x %x\n", received_polynomial[i], encoded_polynomial[i]);
				if(received_polynomial[i] != encoded_polynomial[i])
				{
					uncoded_symbol_err = uncoded_symbol_err + 1;
					err_num = err_num + 1;
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
#if (1 == GF_CAL_COUNT)
					symbol_err_this_frame++;
#endif
				}
			}

#if (1 == RE_ENCODING)
			chnl_rel_cal(recv_seq, symbol_num);
#endif

			mul_assign();
			
			re_encoding();
			
			as_decoding();

#if (1 == GF_CAL_COUNT)			
			gf_count_hist(symbol_err_this_frame);
			symbol_err_this_frame = 0;
#endif
			uncoded_frame_err_flag = 0;

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

			if((1 == decoding_ok_flag)
				&& (1 == frame_err_flag))
			{
				if(0 == hamm_distance_debug)
				{
					hamm_err = hamm_err + 1;

					for(i = 0; i < MESSAGE_LEN; i++)
					{
						DEBUG_SYS("Strage Err: %x %x\n",
							      decoded_message[i],
							      message_polynomial[i]);
					}
				}
				else
				{
					DEBUG_SYS("Prog. Err. for Decoding\n");
#if (1 == OUTPUT_LOG)					
					frc = fopen(log_name, "a+");
					fprintf(frc, "Prog. Err. for Decoding\n");
#endif					

					DEBUG_SYS("Para.: %ld %ld %ld %ld\n",
							  (S_MUL * (CODEWORD_LEN - err_num)),
							  weight_stored,
							  hamm_distance_debug,
							  symbol_err - symbol_err_prev);
#if (1 == OUTPUT_LOG)
					fprintf(frc, "Para.: %ld %ld %ld %ld\n",
							(S_MUL * (CODEWORD_LEN - err_num)),
							weight_stored,
							hamm_distance_debug,
							symbol_err - symbol_err_prev);
#endif

#if (0 == RE_ENCODING)

					for(i = 0; i < MESSAGE_LEN; i++)
					{
						DEBUG_SYS("message: %x\n", message_polynomial[i]);
#if (1 == OUTPUT_LOG)						
						fprintf(frc, "message: %x\n", message_polynomial[i]);
#endif
					}
					for(i = 0; i < CODEWORD_LEN; i++)
					{
						DEBUG_SYS("encoded: %x\n", encoded_polynomial[i]);
#if (1 == OUTPUT_LOG)						
						fprintf(frc, "encoded: %x\n", encoded_polynomial[i]);
#endif
					}
					for(i = 0; i < CODEWORD_LEN; i++)
					{
						DEBUG_SYS("recv: %x\n", received_polynomial[i]);
#if (1 == OUTPUT_LOG)						
						fprintf(frc, "recv: %x\n", received_polynomial[i]);
#endif
					}
					for(i = 0; i < CODEWORD_LEN; i++)
					{
						DEBUG_SYS("err: %x\n", gf_add(received_polynomial[i], encoded_polynomial[i]));
#if (1 == OUTPUT_LOG)						
						fprintf(frc, "err: %x\n", gf_add(received_polynomial[i], encoded_polynomial[i]));
#endif
					}

#endif
#if (1 == OUTPUT_LOG)
					fclose(frc);
					frc = NULL;
#endif
				}
			}
			frame_err_flag = 0;
#if (0 == RECUR_RR)
			if(1 == rr_err)
			{
				rr_err_cnt = rr_err_cnt + 1;
			}
#endif
			if(0 == (iter % monitor_cnt))
			{
				stop = clock();
				runtime = (stop - start) / 1000.0000;
				
				DEBUG_SYS("---------------------\n");
				DEBUG_SYS("Time: %fs\n", runtime);
				DEBUG_SYS("Eb/N0: %f dB\n", eb2n0);
				DEBUG_SYS("Frame: %ld\n", iter + 1);
				DEBUG_SYS("Uncoded Frame Error: %ld\n", uncoded_frame_err);
				DEBUG_SYS("Uncoded Symbol Error: %ld\n", uncoded_symbol_err);
				DEBUG_SYS("Uncoded Bit Error: %ld\n", uncoded_bit_err);
				DEBUG_SYS("Frame Error: %ld\n", frame_err);
				DEBUG_SYS("Symbol Error: %ld\n", symbol_err);
				DEBUG_SYS("Bit Error: %ld\n", bit_err);
				DEBUG_SYS("Hamming Error: %ld\n", hamm_err);
#if (0 == RECUR_RR)				
				DEBUG_SYS("RR Error: %ld\n", rr_err_cnt);
#endif
#if (1 == GF_CAL_COUNT)				
				DEBUG_SYS("Add Cnt: %ld\n", add_cnt);
				DEBUG_SYS("Mul Cnt: %ld\n", mul_cnt);
				DEBUG_SYS("Div Cnt: %ld\n", div_cnt);
				DEBUG_SYS("RCB Cnt: %ld\n", real_cbm_cnt);
				DEBUG_SYS("RMF Cnt: %ld\n", real_mul_ff_cnt);
				DEBUG_SYS("Pow Cnt: %ld\n", pow_cnt);
#endif				
				DEBUG_SYS("Max DX: %ld\n", max_dx);
				DEBUG_SYS("Max DY: %ld\n", max_dy);
				DEBUG_SYS("TERM_SIZE: %ld %ld\n", term_size_x, term_size_y);

#if (1 == OUTPUT_LOG)
				frc = fopen(log_name, "a+");
				fprintf(frc, "---------------------\n");
				fprintf(frc, "Time: %fs\n", runtime);
				fprintf(frc, "Eb/N0: %f dB\n", eb2n0);
				fprintf(frc, "Frame: %ld\n", iter + 1);
				fprintf(frc, "Uncoded Frame Error: %ld\n", uncoded_frame_err);
				fprintf(frc, "Uncoded Symbol Error: %ld\n", uncoded_symbol_err);
				fprintf(frc, "Uncoded Bit Error: %ld\n", uncoded_bit_err);
				fprintf(frc, "Frame Error: %ld\n", frame_err);
				fprintf(frc, "Symbol Error: %ld\n", symbol_err);
				fprintf(frc, "Bit Error: %ld\n", bit_err);
				fprintf(frc, "Hamming Error: %ld\n", hamm_err);
#if (0 == RECUR_RR)				
				fprintf(frc, "RR Error: %ld\n", rr_err_cnt);
#endif
#if (1 == GF_CAL_COUNT)				
				fprintf(frc, "Add Cnt: %ld\n", add_cnt);
				fprintf(frc, "Mul Cnt: %ld\n", mul_cnt);
				fprintf(frc, "Div Cnt: %ld\n", div_cnt);
				fprintf(frc, "RCB Cnt: %ld\n", real_cbm_cnt);
				fprintf(frc, "RMF Cnt: %ld\n", real_mul_ff_cnt);
				fprintf(frc, "Pow Cnt: %ld\n", pow_cnt);
#endif				
			    fclose(frc);
				frc = NULL;
#endif

#if (1 == GF_CAL_COUNT)
				for(k = 0; k < (CODEWORD_LEN - MESSAGE_LEN - 1); k++)
				{
					if(0 != err_hist[k])
					{
						//DEBUG_SYS("-------------------\n");
						//DEBUG_SYS("Err Hist: %ld %ld\n", k, err_hist[k]);
						DEBUG_SYS("Add/Mul Hist: %ld %f %f\n", k, (float)(add_cnt_hist[k] / err_hist[k]), (float)(mul_cnt_hist[k] / err_hist[k]));
						//DEBUG_SYS("Mul Hist: %ld %f\n", k, (float)(mul_cnt_hist[k] / err_hist[k]));
						//DEBUG_SYS("Div Hist: %ld %f\n", k, (float)(div_cnt_hist[k] / err_hist[k]));
						//DEBUG_SYS("RCB Hist: %ld %f\n", k, (float)(real_cbm_cnt_hist[k] / err_hist[k]));
						//DEBUG_SYS("RMF Hist: %ld %f\n", k, (float)(real_mul_ff_cnt_hist[k] / err_hist[k]));
						//DEBUG_SYS("Pow Hist: %ld %f\n", k, (float)(pow_cnt_hist[k] / err_hist[k]));

#if (1 == OUTPUT_LOG)
						frc = fopen(log_name, "a+");
						//fprintf(frc, "-------------------\n");
						//fprintf(frc, "Err Hist: %ld %ld\n", k, err_hist[k]);
						fprintf(frc, "Add/Mul Hist: %ld %f %f\n", k, (float)(add_cnt_hist[k] / err_hist[k]), (float)(mul_cnt_hist[k] / err_hist[k]));
						//fprintf(frc, "Mul Hist: %ld %f\n", k, (float)(mul_cnt_hist[k] / err_hist[k]));
						//fprintf(frc, "Div Hist: %ld %f\n", k, (float)(div_cnt_hist[k] / err_hist[k]));
						//fprintf(frc, "RCB Hist: %ld %f\n", k, (float)(real_cbm_cnt_hist[k] / err_hist[k]));
						//fprintf(frc, "RMF Hist: %ld %f\n", k, (float)(real_mul_ff_cnt_hist[k] / err_hist[k]));
						//fprintf(frc, "Pow Hist: %ld %f\n", k, (float)(pow_cnt_hist[k] / err_hist[k]));
						fclose(frc);
						frc = NULL;
#endif						
					}
				}
#endif				
			}

/*early termination for simulation, reduce iteration times*/
/*more than 10 errors are found, and 10% simulation times have been excuted*/
#if (1 == EARLY_TERMINATION)
			if((EARLY_TERMINATION_NUM <= (frame_err - hamm_err))
				&& ((iter_cnt / 10) < iter))
			{
				/*simulation times are enough, go to next Eb/N0 point*/
				break;
			}
#endif

		}

		stop = clock();
		runtime = (stop - start) / 1000.0000;

		DEBUG_SYS("*********************************\n");
		DEBUG_SYS("Time: %fs\n", runtime);
		DEBUG_SYS("Eb/N0: %f dB\n", eb2n0);
		DEBUG_SYS("Frame: %ld\n", iter_cnt + 1);
		DEBUG_SYS("Uncoded Frame Error: %ld\n", uncoded_frame_err);
		DEBUG_SYS("Uncoded Symbol Error: %ld\n", uncoded_symbol_err);
		DEBUG_SYS("Uncoded Bit Error: %ld\n", uncoded_bit_err);
		DEBUG_SYS("Frame Error: %ld\n", frame_err);
		DEBUG_SYS("Symbol Error: %ld\n", symbol_err);
		DEBUG_SYS("Bit Error: %ld\n", bit_err);
		DEBUG_SYS("Hamming Error: %ld\n", hamm_err);
#if (0 == RECUR_RR)		
		DEBUG_SYS("RR Error: %ld\n", rr_err_cnt);
#endif
#if (1 == GF_CAL_COUNT)		
		DEBUG_SYS("Add/Mul Cnt: %ld %ld\n", add_cnt, mul_cnt);
		DEBUG_SYS("Mul Cnt: %ld\n", mul_cnt);
		DEBUG_SYS("Div Cnt: %ld\n", div_cnt);
		DEBUG_SYS("RCB Cnt: %ld\n", real_cbm_cnt);
		DEBUG_SYS("RMF Cnt: %ld\n", real_mul_ff_cnt);
		DEBUG_SYS("Pow Cnt: %ld\n", pow_cnt);
#endif
		DEBUG_SYS("Uncoded Results: %.10lf %.10lf %.10lf\n", 
			    (double)uncoded_frame_err / (double)(iter + 1),
			    (double)uncoded_symbol_err / (double)(iter + 1) / CODEWORD_LEN * BITS_PER_SYMBOL_BPSK,
			    (double)uncoded_bit_err / (double)(iter + 1) / CODEWORD_LEN * BITS_PER_SYMBOL_BPSK / GF_Q);
		DEBUG_SYS("Decoded Results: %.10lf %.10lf %.10lf %.10lf\n", 
				(double)frame_err / (double)(iter + 1),
				(double)symbol_err / (double)(iter + 1) / CODEWORD_LEN * BITS_PER_SYMBOL_BPSK,
				(double)bit_err / (double)(iter + 1) / CODEWORD_LEN * BITS_PER_SYMBOL_BPSK / GF_Q,
				(double)(frame_err - hamm_err) / (double)(iter + 1));

#if (1 == OUTPUT_LOG)
		frc = fopen(log_name, "a+");
		fprintf(frc, "*********************************\n");
		fprintf(frc, "Time: %fs\n", runtime);
		fprintf(frc, "Eb/N0: %f dB\n", eb2n0);
		fprintf(frc, "Frame: %ld\n", iter_cnt + 1);
		fprintf(frc, "Uncoded Frame Error: %ld\n", uncoded_frame_err);
		fprintf(frc, "Uncoded Symbol Error: %ld\n", uncoded_symbol_err);
		fprintf(frc, "Uncoded Bit Error: %ld\n", uncoded_bit_err);
		fprintf(frc, "Frame Error: %ld\n", frame_err);
		fprintf(frc, "Symbol Error: %ld\n", symbol_err);
		fprintf(frc, "Bit Error: %ld\n", bit_err);
		fprintf(frc, "Hamming Error: %ld\n", hamm_err);
#if (0 == RECUR_RR)		
		fprintf(frc, "RR Error: %ld\n", rr_err_cnt);
#endif
#if (1 == GF_CAL_COUNT)		
		fprintf(frc, "Add Cnt: %ld\n", add_cnt);
		fprintf(frc, "Mul Cnt: %ld\n", mul_cnt);
		fprintf(frc, "Div Cnt: %ld\n", div_cnt);
		fprintf(frc, "RCB Cnt: %ld\n", real_cbm_cnt);
		fprintf(frc, "RMF Cnt: %ld\n", real_mul_ff_cnt);
		fprintf(frc, "Pow Cnt: %ld\n", pow_cnt);
#endif		
		fprintf(frc, "Uncoded Results: %.10lf %.10lf %.10lf\n", 
			    (double)uncoded_frame_err / (double)(iter + 1),
			    (double)uncoded_symbol_err / (double)(iter + 1) / CODEWORD_LEN * BITS_PER_SYMBOL_BPSK,
			    (double)uncoded_bit_err / (double)(iter + 1) / CODEWORD_LEN * BITS_PER_SYMBOL_BPSK / GF_Q);
		fprintf(frc, "Decoded Results: %.10lf %.10lf %.10lf %.10lf\n", 
				(double)frame_err / (double)(iter + 1),
				(double)symbol_err / (double)(iter + 1) / CODEWORD_LEN * BITS_PER_SYMBOL_BPSK,
				(double)bit_err / (double)(iter + 1) / CODEWORD_LEN * BITS_PER_SYMBOL_BPSK / GF_Q,
				(double)(frame_err - hamm_err) / (double)(iter + 1));
	    fclose(frc);
		frc = NULL;
#endif

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

#if (0 == TEST_MODE)
	g_term_destroy();
#endif

	return;
}
