#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "debug_info.h"
#include "gf_cal.h"
#include "encoding.h"
#include "mod.h"
#include "as_decoding.h"

unsigned char received_polynomial[CODEWORD_LEN] =
{
	0xFF,
	0xFF,
	0xFF,
	0xFF,
	0xFF,
	0xFF,
	0xFF
};

unsigned char output_polynomial[CODEWORD_LEN] =
{
	0xFF,
	0xFF,
	0xFF,
	0xFF,
	0xFF,
	0xFF,
	0xFF
};

/*col(locator, 0xff~GF_FIELD-2)-row(mesg, CODEWORD_LEN), same as matlab, contrary to most papers*/
float chnl_rel_matrix[CODEWORD_LEN + 1][CODEWORD_LEN];
float chnl_rel_matrix_tmp[CODEWORD_LEN + 1][CODEWORD_LEN];
unsigned char mul_matrix[CODEWORD_LEN + 1][CODEWORD_LEN];
unsigned char beta_matrix[CODEWORD_LEN + 1][CODEWORD_LEN];
unsigned char rel_group_seq[MESSAGE_LEN];
unsigned char unrel_group_seq[CODEWORD_LEN - MESSAGE_LEN];
unsigned char syndrome[CODEWORD_LEN - MESSAGE_LEN];
unsigned char tao[CODEWORD_LEN - MESSAGE_LEN + 1];
unsigned char omega[CODEWORD_LEN - MESSAGE_LEN];
unsigned char sigma[((CODEWORD_LEN - MESSAGE_LEN) + (CODEWORD_LEN - MESSAGE_LEN + 1) - 1) - (CODEWORD_LEN - MESSAGE_LEN)];
unsigned char erasure_polynomial[CODEWORD_LEN];
unsigned char phi[CODEWORD_LEN];
//unsigned char g_term_c[LAYER_NUM][POLY_NUM][TERM_SIZE * TERM_SIZE];
unsigned long long g_term_x[TERM_SIZE * TERM_SIZE];
unsigned long long g_term_y[TERM_SIZE * TERM_SIZE];
unsigned char g_term_0_y_c[LAYER_NUM][TERM_SIZE * TERM_SIZE];
unsigned char g_term_x_0_c[LAYER_NUM][TERM_SIZE * TERM_SIZE];
unsigned char f_root_val[ROOT_SIZE][ROOT_SIZE];
unsigned char f_root_prev[ROOT_SIZE][ROOT_SIZE];
unsigned long long f_root_cnt[ROOT_SIZE + 1];//used for next layer
unsigned long long sml_poly = 0xFF;
unsigned char decoded_codeword[CODEWORD_LEN];
unsigned char decoded_message[MESSAGE_LEN];

unsigned char ***g_term_c_p;
unsigned char g_term_phase = 0;

unsigned long long err_num = 0;
unsigned char decoding_ok_flag = 0;
unsigned long long weight_stored = 0;
unsigned long long hamm_distance_debug = 0xFFFF;

void find_max_val(float matrix[][CODEWORD_LEN], unsigned long long col,
					 unsigned char* m_ptr, unsigned char* n_ptr)
{
	unsigned long long i = 0, j = 0;
	float max_val = 0;

	for(i = 0; i < col; i++)
	{
		for(j = 0; j < CODEWORD_LEN; j++)
		{
			if(matrix[i][j] >= max_val)
			{
				max_val = matrix[i][j];
				*m_ptr = i;
				*n_ptr = j;
			}
		}
	}

	//DEBUG_INFO("i: %d, j: %d, val: %f\n", *m_ptr, *n_ptr, max_val);

	return;
}

int chnl_rel_init()
{
	unsigned long long i = 0, j = 0;

	for(i = 0; i < CODEWORD_LEN + 1; i++)
	{
		for(j = 0; j < CODEWORD_LEN; j++)
		{
			chnl_rel_matrix[i][j] = 0;
		}
	}

	chnl_rel_matrix[0][1] = 0.12;
	chnl_rel_matrix[0][2] = 0.01;
	chnl_rel_matrix[1][0] = 0.02;
	chnl_rel_matrix[1][1] = 0.01;
	chnl_rel_matrix[1][5] = 0.14;
	chnl_rel_matrix[1][6] = 0.08;
	chnl_rel_matrix[2][0] = 0.01;
	chnl_rel_matrix[2][1] = 0.79;
	chnl_rel_matrix[2][2] = 0.19;
	chnl_rel_matrix[2][3] = 0.13;
	chnl_rel_matrix[2][4] = 0.91;
	chnl_rel_matrix[2][5] = 0.02;
	chnl_rel_matrix[3][0] = 0.64;
	chnl_rel_matrix[3][1] = 0.06;
	chnl_rel_matrix[3][2] = 0.08;
	chnl_rel_matrix[3][4] = 0.01;
	chnl_rel_matrix[3][5] = 0.83;
	chnl_rel_matrix[3][6] = 0.19;
	chnl_rel_matrix[4][2] = 0.01;
	chnl_rel_matrix[4][3] = 0.01;
	chnl_rel_matrix[5][0] = 0.01;
	chnl_rel_matrix[5][2] = 0.01;
	chnl_rel_matrix[5][6] = 0.21;
	chnl_rel_matrix[6][1] = 0.02;
	chnl_rel_matrix[6][2] = 0.49;
	chnl_rel_matrix[6][3] = 0.85;
	chnl_rel_matrix[6][4] = 0.08;
	chnl_rel_matrix[6][6] = 0.01;
	chnl_rel_matrix[7][0] = 0.32;
	chnl_rel_matrix[7][2] = 0.21;
	chnl_rel_matrix[7][6] = 0.51;

	for(i = 0; i < CODEWORD_LEN + 1; i++)
	{
		for(j = 0; j < CODEWORD_LEN; j++)
		{
			chnl_rel_matrix_tmp[i][j] = chnl_rel_matrix[i][j];
		}
	}

	return 0;
}

int chnl_rel_cal()
{
	unsigned long long i = 0, j = 0;

	for(i = 0; i < CODEWORD_LEN; i++)
	{
		for(j = 0; j < (GF_Q / BITS_PER_SYMBOL_BPSK); j++)
		{
			DEBUG_NOTICE("recv_seq: %f\n", fabs(recv_seq[i * (GF_Q / BITS_PER_SYMBOL_BPSK) + j][0]));
			recv_rel[i] = recv_rel[i] + fabs(recv_seq[i * (GF_Q / BITS_PER_SYMBOL_BPSK) + j][0]);
		}
		DEBUG_NOTICE("recv_rel: %f\n", recv_rel[i]);
	}

	return 0;
}

int mul_assign()
{
	unsigned long long i = 0, j = 0;
#if 0	
	unsigned long long s = 0;
	unsigned char *m_ptr = (unsigned char*)malloc(sizeof(unsigned char));
	unsigned char *n_ptr = (unsigned char*)malloc(sizeof(unsigned char));

	chnl_rel_init();

	for(i = 0; i < CODEWORD_LEN + 1; i++)
	{
		for(j = 0; j < CODEWORD_LEN; j++)
		{
			mul_matrix[i][j] = 0;
		}
	}

	for(s = 0; s < S_MUL; s++)
	{
		find_max_val(chnl_rel_matrix_tmp, CODEWORD_LEN + 1,
					  m_ptr, n_ptr);
		i = *m_ptr;
		j = *n_ptr;
		chnl_rel_matrix_tmp[i][j] = chnl_rel_matrix[i][j] / (mul_matrix[i][j] + 2);
		mul_matrix[i][j] = mul_matrix[i][j] + 1;
	}

	DEBUG_IMPOTANT("Multiplicity Assignment:\n");
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		for(j = 0; j < CODEWORD_LEN + 1; j++)
		{
			DEBUG_IMPOTANT("%d ", mul_matrix[j][i]);
		}
		DEBUG_IMPOTANT("\n");
	}

	free(m_ptr);
	free(n_ptr);
#else//GS algorithm
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		for(j = 0; j < CODEWORD_LEN + 1; j++)
		{
			mul_matrix[j][i] = 0;

			if(power_polynomial_table[j][0] == received_polynomial[i])
			{
				mul_matrix[j][i] = S_MUL;
#if 0//for ASD test
				if(encoded_polynomial[i] != received_polynomial[i])
				{
					mul_matrix[j][i] = 1;
				}
#endif
#if (1 == SIMPLE_ASD)//for ASD test, a simple multiplicity assignment strategy
				if(GF_Q > recv_rel[i])
				{
					mul_matrix[j][i] = S_MUL / 2;
				}
				if((GF_Q / 2) > recv_rel[i])
				{
					mul_matrix[j][i] = S_MUL / 4;
				}
#endif

			}
		}
	}

#if (1 == CFG_DEBUG_IMPOTANT)
	DEBUG_IMPOTANT("Multiplicity Assignment:\n");
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		for(j = 0; j < CODEWORD_LEN + 1; j++)
		{
			DEBUG_IMPOTANT("%d ", mul_matrix[j][i]);
		}
		DEBUG_IMPOTANT("\n");
	}
#endif	
#endif
	
	return 0;
}

unsigned char syndrome_cal(unsigned char *recv, unsigned char *synd,
								unsigned long long cw_len, unsigned long long msg_len)
{
	unsigned long long i = 0, j = 0;
	unsigned char tmp = 0xFF, tmp_sum = 0xFF;

	for(i = 0; i < cw_len - msg_len; i++)
	{
		tmp = 0xFF;
		tmp_sum = 0xFF;
		for(j = 0; j < cw_len; j++)
		{
			tmp = gf_multp(recv[j], (i + 1) * j);
			tmp_sum = gf_add(tmp, tmp_sum);
			//DEBUG_INFO("%x %x\n", tmp, tmp_sum);
		}
		synd[i] = tmp_sum;
	}
#if (1 == CFG_DEBUG_INFO)	
	DEBUG_INFO("Syndrome Polynomial:\n");
	for(i = 0; i < cw_len - msg_len; i++)
	{
		DEBUG_INFO("%x ", synd[i]);
	}
	DEBUG_INFO("\n");
#endif	
}

int rel_group()
{
	unsigned long long i = 0, j = 0;

	float rel_thrd = 0.7;
	unsigned long long rel_cnt = 0, rel_flag = 0;

	while(MESSAGE_LEN != rel_cnt)
	{
		rel_cnt = 0;
		rel_flag = 0;
		memset(rel_group_seq, 0, sizeof(unsigned char) * MESSAGE_LEN);
		memset(unrel_group_seq, 0, sizeof(unsigned char) * (CODEWORD_LEN - MESSAGE_LEN));
	
		for(j = 0; j < CODEWORD_LEN; j++)
		{
			for(i = 0; i < CODEWORD_LEN + 1; i++)
			{
				if(rel_thrd < chnl_rel_matrix[i][j])
				{
					rel_flag = 1;
					//DEBUG_NOTICE("rel_val: %d %d %f %f\n", i, j, rel_thrd, chnl_rel_matrix[i][j]);
					break;
				}
			}
			if(1 == rel_flag)
			{
				if(MESSAGE_LEN > rel_cnt)
				{
					rel_group_seq[rel_cnt] = j;
				}
				rel_cnt = rel_cnt + 1;
			}
			else
			{
				if(MESSAGE_LEN >= rel_cnt)
				{
					unrel_group_seq[j - rel_cnt] = j;
				}
			}
			rel_flag = 0;
		}

		//DEBUG_NOTICE("rel_group: %f %d\n", rel_thrd, rel_cnt);

		if(MESSAGE_LEN < rel_cnt)
		{
			rel_thrd = rel_thrd + 0.05;
		}
		else
		{
			rel_thrd = rel_thrd - 0.05;
		}
	}

#if (1 == CFG_DEBUG_INFO)
	DEBUG_INFO("rel_seq: ");
	for(i = 0; i < MESSAGE_LEN; i++)
	{
		DEBUG_INFO("%d ", rel_group_seq[i]);
	}
	DEBUG_INFO("\n");
	DEBUG_INFO("unrel_seq: ");
	for(i = 0; i < CODEWORD_LEN - MESSAGE_LEN; i++)
	{
		DEBUG_INFO("%d ", unrel_group_seq[i]);
	}
	DEBUG_INFO("\n");
#endif
	
	return 0;
}

int tao_cal()
{
	unsigned long long i = 0, j = 0;
	memset(tao, 0xFF, sizeof(unsigned char) * (CODEWORD_LEN - MESSAGE_LEN + 1));

	unsigned char a[2] = {0, 0};
	unsigned char reg[CODEWORD_LEN - MESSAGE_LEN + 1];/*n-k+1, R' has (n-k) terms*/
	memset(reg, 0xFF, sizeof(unsigned char) * (CODEWORD_LEN - MESSAGE_LEN + 1));
#if 0
	reg[0] = unrel_group_seq[0];
	reg[1] = 0;
#endif
	unsigned char b[CODEWORD_LEN - MESSAGE_LEN];
	memset(b, 0xFF, sizeof(unsigned char) * (CODEWORD_LEN - MESSAGE_LEN));
#if 0
	for(i = 1; i < (CODEWORD_LEN - MESSAGE_LEN); i++)
	{
		a[0] = unrel_group_seq[i];
		a[1] = 0;
		memcpy(b, reg, sizeof(unsigned char) * (1 + i));
		gf_multp_poly_hw(a, 2,
						  b, (1 + i),
						  reg, (2 + i));
	}
#else
	reg[0] = 0;
	reg[1] = unrel_group_seq[0];
	for(i = 1; i < (CODEWORD_LEN - MESSAGE_LEN); i++)
	{
		a[0] = 0;
		a[1] = unrel_group_seq[i];
		memcpy(b, reg, sizeof(unsigned char) * (1 + i));
		gf_multp_poly_hw(a, 2,
						  b, (1 + i),
						  reg, (2 + i));
	}
#endif
	memcpy(tao, reg, sizeof(unsigned char) * (CODEWORD_LEN - MESSAGE_LEN + 1));

#if (1 == CFG_DEBUG_INFO)
	DEBUG_INFO("tao: ");
	for(j = 0; j < (CODEWORD_LEN - MESSAGE_LEN + 1); j++)
	{
		DEBUG_INFO("%x ", reg[j]);
	}
	DEBUG_INFO("\n");
#endif

	return 0;
}

int sigma_cal()
{
	unsigned long long i = 0, j = 0;
	memset(sigma, 0xFF, sizeof(unsigned char) * (((CODEWORD_LEN - MESSAGE_LEN) + (CODEWORD_LEN - MESSAGE_LEN + 1) - 1) - (CODEWORD_LEN - MESSAGE_LEN)));

	unsigned char tmp[(CODEWORD_LEN - MESSAGE_LEN) + (CODEWORD_LEN - MESSAGE_LEN + 1) - 1];
	memset(tmp, 0xFF, sizeof(unsigned char) * ((CODEWORD_LEN - MESSAGE_LEN) + (CODEWORD_LEN - MESSAGE_LEN + 1) - 1));

	gf_multp_poly_hw(syndrome, (CODEWORD_LEN - MESSAGE_LEN),
					  tao, 		  (CODEWORD_LEN - MESSAGE_LEN + 1),
					  tmp, 		  ((CODEWORD_LEN - MESSAGE_LEN) + (CODEWORD_LEN - MESSAGE_LEN + 1) - 1));

	memcpy(sigma, tmp + (CODEWORD_LEN - MESSAGE_LEN), sizeof(unsigned char) * (((CODEWORD_LEN - MESSAGE_LEN) + (CODEWORD_LEN - MESSAGE_LEN + 1) - 1) - (CODEWORD_LEN - MESSAGE_LEN)));

#if (1 == CFG_DEBUG_NOTICE)
	DEBUG_NOTICE("sigma: ");
	for(i = 0; i < (((CODEWORD_LEN - MESSAGE_LEN) + (CODEWORD_LEN - MESSAGE_LEN + 1) - 1) - (CODEWORD_LEN - MESSAGE_LEN)); i++)
	{
		DEBUG_NOTICE("%x ", sigma[i]);
	}
	DEBUG_NOTICE("\n");
#endif

	memcpy(omega, tmp, sizeof(unsigned char) * (CODEWORD_LEN - MESSAGE_LEN));

#if (1 == CFG_DEBUG_NOTICE)
	DEBUG_NOTICE("omega: ");
	for(i = 0; i < (CODEWORD_LEN - MESSAGE_LEN); i++)
	{
		DEBUG_NOTICE("%x ", omega[i]);
	}
	DEBUG_NOTICE("\n");
#endif

	return 0;
}

unsigned char poly_eva(unsigned char *poly, unsigned char poly_len, unsigned char input_val)
{
	unsigned long long i = 0;
	unsigned char poly_val = 0xFF, tmp_product = 0;

	for(i = 0; i < poly_len; i++)
	{
		tmp_product = gf_multp(*(poly + i), (input_val * i) % (GF_FIELD - 1));
		//DEBUG_NOTICE("poly_eva: %d %x %x %x %x %x\n", i, *(poly + i), input_val, (input_val * i) % (GF_FIELD - 1), tmp_product, poly_val);
		poly_val = gf_add(tmp_product, poly_val);
	}
	
	return poly_val;
}

int phi_cal()
{
	unsigned long long i = 0, k = 0;
	unsigned char find_flag = 0;
	unsigned char tao_dev[(CODEWORD_LEN - MESSAGE_LEN + 1) - 1]; 
	unsigned char tmp = 0xFF, locator = 0xFF;
#if 0
	unrel_group_seq[0] = 0;
	unrel_group_seq[1] = 3;
	unrel_group_seq[2] = 5;
	unrel_group_seq[3] = 1;
	//tao[0] = 0;
	//tao[1] = 5;
	//tao[2] = 4;
	//tao[3] = 0;
	//tao[4] = 4;
	received_polynomial[0] = 0xFF;
	received_polynomial[1] = 0xFF;
	received_polynomial[2] = 0;
	received_polynomial[3] = 0xFF;
	received_polynomial[4] = 0;
	received_polynomial[5] = 0xFF;
	received_polynomial[6] = 0xFF;
	//omega[0] = 3;
	//omega[1] = 2;
	//omega[2] = 6;
	//omega[3] = 5;
#endif
	memset(phi, 0xFF, sizeof(unsigned char) * CODEWORD_LEN);
	memcpy(erasure_polynomial, received_polynomial, sizeof(unsigned char) * CODEWORD_LEN);
	for(i = 0; i < (CODEWORD_LEN - MESSAGE_LEN); i++)
	{
		for(k = 0; k < (CODEWORD_LEN - MESSAGE_LEN); k++)
		{
			if(unrel_group_seq[k] == i)
			{
				find_flag = 1;
				break;
			}
			else
			{
				find_flag = 0;
			}
		}

		if(1 == find_flag)
		{
			erasure_polynomial[i] = 0xFF;
		}
	}
	find_flag = 0;
	
	syndrome_cal(erasure_polynomial, syndrome,
				  CODEWORD_LEN, MESSAGE_LEN);
	tao_cal();
	sigma_cal();

	for(i = 0; i < (CODEWORD_LEN - MESSAGE_LEN); i++)
	{
		if(0 != (i % 2))
		{
			tao_dev[i] = 0xFF;
		}
		else
		{
			tao_dev[i] = tao[i + 1];
		}
	}

#if (1 == CFG_DEBUG_NOTICE)	
	DEBUG_NOTICE("tao_dev: ");
	for(i = 0; i < (CODEWORD_LEN - MESSAGE_LEN); i++)
	{
		DEBUG_NOTICE("%x ", tao_dev[i]);
	}
	DEBUG_NOTICE("\n");
#endif

	for(i = 0; i < CODEWORD_LEN; i++)
	{
		for(k = 0; k < (CODEWORD_LEN - MESSAGE_LEN); k++)
		{
			if(unrel_group_seq[k] == i)
			{
				find_flag = 1;
				break;
			}
			else
			{
				find_flag = 0;
			}
		}

		if(1 == find_flag)
		{
			locator = power_polynomial_table[i + 1][0];
			locator = ((GF_FIELD - 1) - locator) % (GF_FIELD - 1);
			//DEBUG_NOTICE("%x: %x\n", i, locator);
			tmp = poly_eva(omega, CODEWORD_LEN - MESSAGE_LEN, locator);
			//DEBUG_NOTICE("%x ", tmp);
			//tmp = gf_div(tmp, (locator * (CODEWORD_LEN - MESSAGE_LEN + 1)) % (GF_FIELD - 1));
			//DEBUG_NOTICE("%x ", tmp);
			tmp = gf_div(tmp, poly_eva(tao_dev, CODEWORD_LEN - MESSAGE_LEN, locator));
			//DEBUG_NOTICE("%x ", tmp);
			//phi[i] = gf_add(tmp, received_polynomial[i]);
			phi[i] = tmp;
			//DEBUG_NOTICE("\n");
		}
		else
		{
			phi[i] = received_polynomial[i];
		}
	}

#if (1 == CFG_DEBUG_INFO)
	DEBUG_INFO("phi: ");
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		DEBUG_INFO("%x ", phi[i]);
	}
	DEBUG_INFO("\n");
#if 0
	DEBUG_INFO("code_val: ");
	for(i = 0; i < (GF_FIELD - 1); i++)
	{
		DEBUG_INFO("%x ", poly_eva(phi, (CODEWORD_LEN), power_polynomial_table[i + 1][0]));
	}
	DEBUG_INFO("\n");
#endif	
#endif

	return 0;
}

unsigned char coordinate_trans(unsigned char locator, unsigned char r, unsigned char r_hd)
{
	unsigned long long i = 0, j = 0;
	unsigned char coordinate = 0xFF;
	unsigned char locator_product = 0, tmp_sum = 0xFF;
#if 0
	for(i = 0; i < CODEWORD_LEN - MESSAGE_LEN; i++)
	{
		if(locator == unrel_group_seq[i])
		{
			continue;
		}
		tmp_sum = gf_add(locator, unrel_group_seq[i]);
		locator_product = gf_multp(locator_product, tmp_sum);
	}
	//locator_product = gf_multp(0, locator);
	locator_product = gf_multp(locator_product, locator);
	tmp_sum = gf_add(r_hd, r);
	locator_product = gf_multp(tmp_sum, locator_product);

	tmp_sum = poly_eva(sigma, ((CODEWORD_LEN - MESSAGE_LEN) + (CODEWORD_LEN - MESSAGE_LEN + 1) - 1) - (CODEWORD_LEN - MESSAGE_LEN), locator);

	coordinate = gf_add(locator_product, tmp_sum);
#else
	if(0xFF == gf_add(r_hd, r))
	{
		return 0xFF;
	}

	for(i = 0; i < MESSAGE_LEN; i++)
	{
		tmp_sum = gf_add(locator, rel_group_seq[i]);
		locator_product = gf_multp(locator_product, tmp_sum);
	}
	coordinate = gf_div(gf_add(r_hd, r), locator_product);
#endif
	return coordinate;
}

int re_encoding()
{
	unsigned long long i = 0, j = 0, k = 0, l = 0;
	unsigned long long tmp = 0;
	unsigned char find_flag = 0;
	memset(syndrome, 0xFF, sizeof(unsigned char) * (CODEWORD_LEN - MESSAGE_LEN));

	syndrome_cal(received_polynomial, syndrome,
				  CODEWORD_LEN, MESSAGE_LEN);
#if 0//use re-encoding
#if 0
	DEBUG_INFO("syn_val: ");
	for(i = 0; i < (GF_FIELD - 1); i++)
	{
		DEBUG_INFO("%x ", poly_eva(syndrome, (CODEWORD_LEN - MESSAGE_LEN), power_polynomial_table[i + 1][0]));
	}
	DEBUG_INFO("\n");
	DEBUG_INFO("code_val: ");
	for(i = 0; i < (GF_FIELD - 1); i++)
	{
		DEBUG_INFO("%x ", poly_eva(encoded_polynomial, (CODEWORD_LEN), power_polynomial_table[i + 1][0]));
	}
	DEBUG_INFO("\n");
#endif

	rel_group();
#if 0
	tao_cal();
	sigma_cal();
#endif

	phi_cal();
#if 0
	for(i = 0; i < (CODEWORD_LEN + 1); i++)
	{
		for(j = 0; j < CODEWORD_LEN; j++)
		{
			find_flag = 0;
			for(k = 0; k < MESSAGE_LEN; k++)/*it is inefficient in C, but it is easy in verilog*/
			{
				if(rel_group_seq[k] == j)
				{
					find_flag = 1;
					break;
				}
			}
#if 1			
			//if((0 == mul_matrix[i][j]) || (1 == find_flag))
			if(1 == find_flag)
			{
				beta_matrix[i][j] = 0xFF;
			}
			else
			{
				//DEBUG_INFO("%d %d %x\n", i, j, mul_matrix[i][j]);
				beta_matrix[i][j] = coordinate_trans(power_polynomial_table[j + 1][0], power_polynomial_table[i][0], received_polynomial[j]);
			}
#else
			beta_matrix[i][j] = coordinate_trans(power_polynomial_table[j + 1][0], power_polynomial_table[i][0], received_polynomial[j]);
#endif
		}
	}
#else
	for(i = 0; i < (CODEWORD_LEN + 1); i++)
	{
		for(j = 0; j < CODEWORD_LEN; j++)
		{
			find_flag = 0;
			for(k = 0; k < MESSAGE_LEN; k++)/*it is inefficient in C, but it is easy in verilog*/
			{
				if(rel_group_seq[k] == j)
				{
					find_flag = 1;
					break;
				}
			}

			if(1 == find_flag)
			{
				beta_matrix[i][j] = 0xFF;
			}
			else
			{
				beta_matrix[i][j] = coordinate_trans(power_polynomial_table[j + 1][0], power_polynomial_table[i][0], phi[j]);
			}
		}
	}
#endif
#else//do not use re-encoding
	for(j = 0; j < CODEWORD_LEN; j++)
	{
		for(i = 0; i < (CODEWORD_LEN + 1); i++)
		{
			if(0 != mul_matrix[i][j])
			{
				beta_matrix[i][j] = power_polynomial_table[i][0];
			}
			else
			{
				beta_matrix[i][j] = 0xFF;
			}
		}
	}
#endif

#if (1 == CFG_DEBUG_IMPOTANT)
	DEBUG_IMPOTANT("beta:\n");
	for(j = 0; j < CODEWORD_LEN; j++)
	{
		for(i = 0; i < (CODEWORD_LEN + 1); i++)
		{
			DEBUG_IMPOTANT("%x ", beta_matrix[i][j]);
		}
		DEBUG_IMPOTANT("\n");
	}
	DEBUG_IMPOTANT("\n");
#endif	

	return 0;
}

int lex_order(unsigned long long **lex_table, unsigned long long d_x, unsigned long long d_y)
{
	unsigned long long i = 0, j = 0, k = 0;
	unsigned long long max_degree = 0, degree_index = 0;
	unsigned long long tmp_lex_table[d_x][d_y];

	DEBUG_INFO("lex_order: %ld %ld\n", d_x, d_y);

	DEBUG_INFO("tmp_lex_table:\n");
	for(i = 0; i < d_x; i++)
	{
		for(j = 0; j < d_y; j++)
		{
#if 1	
			tmp_lex_table[i][j] = 1 * i + (MESSAGE_LEN - 1) * j;
#else//for test
			*((unsigned long long *)lex_table + i * d_y + j) = 1 * i + 1 * j;
#endif
			if(max_degree < tmp_lex_table[i][j])
			{
				max_degree = tmp_lex_table[i][j];
			}
			DEBUG_INFO("%ld ", tmp_lex_table[i][j]);
		}
		DEBUG_INFO("\n");
	}
	DEBUG_INFO("\n");

	DEBUG_INFO("max_degree: %ld\n", max_degree);

	for(k = 0; k <= max_degree; k++)
	{
#if(1 == RELEX_ORDER)
		for(j = 0; j < d_y; j++)
		{
			for(i = 0; i < d_x; i++)
#else
		for(i = 0; i < d_x; i++)
		{
			for(j = 0; j < d_y; j++)
#endif
			{
				if(k == tmp_lex_table[i][j])
				{
					*((unsigned long long *)lex_table + i * d_y + j) = degree_index;
					degree_index = degree_index + 1;
					//DEBUG_INFO("lex_table: %ld %ld\n", degree_index, *((unsigned long long *)lex_table + i * d_y + j));
				}
			}
		}
	}
}

int koetter_interpolation()
{
	unsigned long long i = 0, j = 0, k = 0, m = 0, n = 0;
	unsigned long long a = 0, b = 0, v = 0;
	unsigned long long tmp_sum = 0, tmp_real = 0;
	unsigned char tmp_ff = 0xFF;
	unsigned long long l_s = 0xFFFF, l_w = 0, l_o = 0xFFFF;;

	unsigned long long d_x = 0, d_y = 0, c = 0, term_num = 0, term_num_real = 0;
	unsigned long long d_x_max = 0, d_y_max = 0;
	for(i = 0; i < (CODEWORD_LEN + 1); i++)
	{
		for(j = 0; j < CODEWORD_LEN; j++)
		{
			tmp_sum = tmp_sum + mul_matrix[i][j] * (mul_matrix[i][j] + 1);
		}
	}
	c = tmp_sum / 2;
	tmp_sum = pow((1 + 8 * (float)c / (MESSAGE_LEN - 1)), 0.5);
	tmp_sum = (1 + ((unsigned long long)tmp_sum)) / 2;
	d_y = (unsigned long long)(floor(tmp_sum));
	//d_y = floor((1 + (unsigned long long)pow((1 + 8 * (float)c / (MESSAGE_LEN - 1)), 0.5)) / 2);
	d_x = floor(c / (d_y + 1) + d_y * (CODEWORD_LEN - MESSAGE_LEN - 1) / 2);

	//d_x = d_x << LEX_TABLE_EXPAND_SIZE;
	//d_y = d_y << LEX_TABLE_EXPAND_SIZE;
	d_x = d_x * LEX_TABLE_EXPAND_SIZE;
	d_y = d_y * LEX_TABLE_EXPAND_SIZE;

	unsigned long long lex_order_table[d_x][d_y];
	lex_order((unsigned long long **)lex_order_table, d_x, d_y);
	DEBUG_INFO("lex_order_table:\n");
	for(i = 0; i < d_x; i++)
	{
		for(j = 0; j < d_y; j++)
		{
			DEBUG_INFO("%ld ", lex_order_table[i][j]);
			if(c >= lex_order_table[i][j])
			{
				//term_num = term_num + 1;
#if 0				
				if(d_x_max < i)
				{
					d_x_max = i;
				}
				if(d_y_max < j)
				{
					d_y_max = j;
				}
#endif
				if(d_y_max < j)
				{
					d_y_max = j;
				}
			}
		}
		DEBUG_INFO("\n");
	}

	d_x_max = d_x;

	for(i = 0; i < (d_x + 1); i++)
	{
		//for(j = i; j < (d_y_max + 1); j++)
		for(j = 0; j < (d_y + 1); j++)
		{
			term_num_real = term_num_real + 1;
		}
	}
	
	DEBUG_INFO("constraint: %d %d %d %d %d %d\n", c, d_x, d_y, term_num_real, d_x_max, d_y_max);
#if 0
	unsigned char g_table_c[d_y_max + 1][term_num];
	unsigned char g_table_x[term_num];
	unsigned char g_table_y[term_num];
	unsigned char discrepancy[d_y_max + 1];
	unsigned char weight_pol[d_y_max + 1];
	unsigned char tmp_table_c[term_num];
#endif
	unsigned char **g_table_c;
	unsigned char **g_table_c_prev;
	unsigned long long *g_table_x;
	unsigned long long *g_table_y;
	unsigned char *discrepancy;
	unsigned long long *weight_pol;
	unsigned long long *lexorder_pol;
	unsigned char *tmp_table_c;
	unsigned long long *term_use_index;//malloc later

	g_table_c = (unsigned char**)malloc(sizeof(unsigned char*) * (d_y_max + 1));
	g_table_c_prev = (unsigned char**)malloc(sizeof(unsigned char*) * (d_y_max + 1));
	for (i = 0; i < (d_y_max + 1); i++)
	{
  		g_table_c[i] = (unsigned char*)malloc(sizeof(unsigned char) * term_num_real);
		g_table_c_prev[i] = (unsigned char*)malloc(sizeof(unsigned char) * term_num_real);
  	}
	g_table_x = (unsigned long long*)malloc(sizeof(unsigned long long) * term_num_real);
	g_table_y = (unsigned long long*)malloc(sizeof(unsigned long long) * term_num_real);
	discrepancy = (unsigned char*)malloc(sizeof(unsigned char) * (d_y_max + 1));
	weight_pol = (unsigned long long*)malloc(sizeof(unsigned long long) * (d_y_max + 1));
	lexorder_pol = (unsigned long long*)malloc(sizeof(unsigned long long) * (d_y_max + 1));
	tmp_table_c = (unsigned char*)malloc(sizeof(unsigned char) * term_num_real);

	for(i = 0; i < (d_y_max + 1); i++)
	{
		for(j = 0; j < term_num_real; j++)
		{
			g_table_c[i][j] = 0xFF;
			g_table_c_prev[i][j] = 0xFF;
		}
	}
	memset(g_table_x, 0, sizeof(unsigned long long) * term_num_real);
	memset(g_table_y, 0, sizeof(unsigned long long) * term_num_real);
	memset(discrepancy, 0xFF, sizeof(unsigned char) * (d_y_max + 1));
	memset(tmp_table_c, 0xFF, sizeof(unsigned char) * term_num_real);
	memset(weight_pol, 0, sizeof(unsigned long long) * (d_y_max + 1));
	memset(lexorder_pol, 0xFFFF, sizeof(unsigned long long) * (d_y_max + 1));

	/*init (a, b) pairs*/
	k = 0;
	for(i = 0; i < (d_x + 1); i++)
	{
		for(j = 0; j < (d_y + 1); j++)
		{
			//if(S_MUL > (i + j))
			{
				g_table_x[k] = i;
				g_table_y[k] = j;
				DEBUG_NOTICE("g_table: %d %d %d\n", k, g_table_x[k], g_table_y[k]);
				k = k + 1;
			}
		}
	}
	
	//unsigned char g_term[d_y_max + 1];//j, degree of x, degree of y, pick indexes form g_table
	for(i = 0; i <= d_y_max; i++)
	{
		for(j = 0; j < term_num_real; j++)
		{
			if((0 == g_table_x[j]) && (i == g_table_y[j]))
			{
				g_table_c[i][j] = 0; //set coefficient of initial g to 1
				DEBUG_NOTICE("g_table_c: %d %d %x\n", i, j, g_table_c[i][j]);
				g_table_c_prev[i][j] = g_table_c[i][j];
				break;
			}
		}
	}
	for(i = 0; i < (d_y_max + 1); i++)
	{
		for(j = 0; j < term_num_real; j++)
		{
			if(0xFF != g_table_c[i][j])
			{
				//DEBUG_NOTICE("weight_pol_searching: %d %d\n", weight_pol[i], lex_order_table[g_table_x[j]][g_table_y[j]]);
#if 0				
				if(weight_pol[i] < lex_order_table[g_table_x[j]][g_table_y[j]])
				{
					weight_pol[i] = lex_order_table[g_table_x[j]][g_table_y[j]];
				}
#else
				if((weight_pol[i] < (g_table_x[j] + (MESSAGE_LEN - 1) * g_table_y[j]))
					|| ((weight_pol[i] == (g_table_x[j] + (MESSAGE_LEN - 1) * g_table_y[j]))
						&& (lexorder_pol[i] > lex_order_table[g_table_x[j]][g_table_y[j]])))
				{
					weight_pol[i] = (g_table_x[j] + (MESSAGE_LEN - 1) * g_table_y[j]);
					lexorder_pol[i] = lex_order_table[g_table_x[j]][g_table_y[j]];
				}
#endif
			}
		}
		//weight_pol[i] = (MESSAGE_LEN - 1) * i;
		DEBUG_INFO("pol: %d %d\n", weight_pol[i], lexorder_pol[i]);
	}

	for(j = 0; j < CODEWORD_LEN; j++) //for I, as 2^q - 1
	{
		for(i = 0; i < (CODEWORD_LEN + 1); i++)//for I, as n
		{
			
			//DEBUG_INFO("point: %d %d | %d %d | (%x %x)\n", CODEWORD_LEN, (CODEWORD_LEN + 1), i, j, power_polynomial_table[j + 1][0], beta_matrix[i][j]);
			if(0 == mul_matrix[i][j])
			{
				continue;
			}
			DEBUG_INFO("-------------------------------------------\n");
			DEBUG_INFO("point: %d %d | %d %d | (%x %x)\n", CODEWORD_LEN, (CODEWORD_LEN + 1), i, j, power_polynomial_table[j + 1][0], beta_matrix[i][j]);

			/*the (a, b) pairs should be init here*/
			term_num = 0;//clear here
			for(m = 0; m < d_x; m++)
			{
				for(n = 0; n < d_y; n++)
				{
					if((mul_matrix[i][j] - 1) >= (m + n))
					{
						term_num = term_num + 1;
					}
				}
			}

			term_use_index = (unsigned long long*)malloc(sizeof(unsigned long long) * term_num);
			m = 0;
			for(k = 0; k < term_num_real; k++)
			{
				if((mul_matrix[i][j] - 1) >= (g_table_x[k] + g_table_y[k]))
				{
					term_use_index[m] = k;
					DEBUG_NOTICE("term_use_index: %d | %d | %d %d\n", m, term_use_index[m], g_table_x[k], g_table_y[k]);
					m = m + 1;
				}
			}

			/*notice that only (a, b) pairs is constained by (a + b) < m directly.*/
			/*(r, s) pairs is related to the point (a, b), but not the constrain.*/
			for(k = 0; k < term_num; k++) //for II, as (a, b) pairs
			{
				DEBUG_NOTICE("*********************\n");
				DEBUG_NOTICE("k ready: %d %d | %d | (%d %d)\n", term_num, k, term_use_index[k], g_table_x[term_use_index[k]], g_table_y[term_use_index[k]]);
				memset(discrepancy, 0xFF, sizeof(unsigned char) * (d_y_max + 1));

				//l_s = 0xFF;
				//l_w = 0;

				for(m = 0; m < (d_y_max + 1); m++) //for III, as Q_d_y_max
				{
					tmp_sum = 0xFF;

					for(n = 0; n < term_num_real; n++) //for IV, as (r, s) pairs
					{
#if 0					
						if(4 == k)
						{
							DEBUG_NOTICE("m, n: %d %d | %d %d | %d %d | %x %x\n",
								    m, n,
									g_table_x[n], g_table_y[n],
									g_table_x[term_use_index[k]], g_table_y[term_use_index[k]],
									g_table_c[m][n], beta_matrix[i][j]);
						}
#endif						
						if((g_table_x[n] < g_table_x[term_use_index[k]])
							|| (g_table_y[n] < g_table_y[term_use_index[k]])
							|| (0xFF == g_table_c_prev[m][n])
							|| ((0xFF == beta_matrix[i][j]) && (g_table_y[n] != g_table_y[term_use_index[k]])))
						{
							continue;
						}
						
#if 0
						tmp_sum = tmp_sum
								+ real_combine(g_table_x[n], g_table_x[k])
								* real_combine(g_table_y[n], g_table_y[k])
								* g_table_c[n]
								* (power_polynomial_table[i + 1][0])^(g_table_x[n] - g_table_x[k])
								* (beta_matrix[i][j])^(g_table_y[n] - g_table_y[k]);
#endif
						DEBUG_NOTICE("++++++++++\n");
						DEBUG_NOTICE("(m, n) ready: %d (%d %d) %x\n", m, g_table_x[n], g_table_y[n], g_table_c_prev[m][n]);
						tmp_real = real_combine(g_table_x[n], g_table_x[term_use_index[k]]) * real_combine(g_table_y[n], g_table_y[term_use_index[k]]);
						//DEBUG_NOTICE("tmp_real: %d | %d %d %d | %d %d %d\n", tmp_real, g_table_x[n], g_table_x[term_use_index[k]], real_combine(g_table_x[n], g_table_x[term_use_index[k]]), g_table_y[n], g_table_y[term_use_index[k]], real_combine(g_table_y[n], g_table_y[term_use_index[k]]));
						tmp_ff = gf_real_mutp_ff(tmp_real, g_table_c_prev[m][n]);
						//DEBUG_NOTICE("tmp_ff: %x = %d mul %x\n", tmp_ff, tmp_real, g_table_c_prev[m][n]);
						//DEBUG_NOTICE("tmp_ff: %x * (%x = %x^%x)\n", tmp_ff, gf_pow_cal(power_polynomial_table[j + 1][0], (g_table_x[n] - g_table_x[term_use_index[k]])), power_polynomial_table[j + 1][0], (g_table_x[n] - g_table_x[term_use_index[k]]));
						tmp_ff = gf_multp(tmp_ff, gf_pow_cal(power_polynomial_table[j + 1][0], (g_table_x[n] - g_table_x[term_use_index[k]])));
						//DEBUG_NOTICE("tmp_ff: %x * (%x = %x^%x)\n", tmp_ff, gf_pow_cal(beta_matrix[i][j], (g_table_y[n] - g_table_y[term_use_index[k]])), beta_matrix[i][j], (g_table_y[n] - g_table_y[term_use_index[k]]));
						tmp_ff = gf_multp(tmp_ff, gf_pow_cal(beta_matrix[i][j], (g_table_y[n] - g_table_y[term_use_index[k]])));
						//DEBUG_NOTICE("gf_add: %x + %x\n", tmp_sum, tmp_ff);
						tmp_sum = gf_add(tmp_sum, tmp_ff);
#if (1 == CFG_DEBUG_NOTICE)						
						if(0xFF != tmp_sum)
						{
							DEBUG_NOTICE("tmp_sum: %x | %x\n", tmp_sum, tmp_ff);
						}
#endif						
					}
					
					if(0xFF != tmp_sum)
					{
						discrepancy[m] = tmp_sum;
						DEBUG_NOTICE("d: %d %d | %d %d | %d | %d\n", i, j, g_table_x[term_use_index[k]], g_table_y[term_use_index[k]], m, discrepancy[m]);
					}

					if(0xFF != discrepancy[m])
					{
						DEBUG_NOTICE("updating place center: %d | %d vs %d | %d vs %d\n", l_s, l_w, weight_pol[m], l_o, lexorder_pol[m]);
#if 1
						if(((l_w > weight_pol[m])
								|| ((l_w == weight_pol[m])
									&& (l_o > lexorder_pol[m])))
#else//there may be some err
						if((l_w >= weight_pol[m])
#endif
							|| (0xFF == discrepancy[l_s]))
						{
							l_s = m;
							l_w = weight_pol[m];
							l_o = lexorder_pol[m];
							DEBUG_NOTICE("updated place center: %d | %d | %d\n", l_s, l_w, l_o);
						}
					}
#if 0					
					if((0xFF != discrepancy[m])//find l_s
						&& (l_w >= weight_pol[m]))
					{
						l_s = m;
						l_w = weight_pol[m];
						DEBUG_NOTICE("updated place center: %d %d\n", l_s, l_w);
					}
#endif
					//DEBUG_NOTICE("g_table_c[0][0]: %x\n", g_table_c[0][0]);
				}

				DEBUG_NOTICE("l_s: %d\n", l_s);

				//DEBUG_NOTICE("g_table_c[0][0]: %x\n", g_table_c[0][0]);
				//DEBUG_NOTICE("update polynomial ready: %d %d\n", l_s, l_w);
				for(m = 0; m < (d_y_max + 1); m++)//update normal Q_l
				{
					if((m == l_s)
						|| (0xFF == discrepancy[m]))
					{
						continue;
					}

					for(n = 0; n < term_num_real; n++)
					{
#if 0
						if((0xFF != g_table_c[m][n])
							&& (0xFF != discrepancy[l_s]))
						{
							DEBUG_NOTICE("delta_l_s * q_l: %d %x %x %x\n", n,
								                                  gf_multp(g_table_c[m][n], discrepancy[l_s]),
														          g_table_c[m][n],
														          discrepancy[l_s]);
						}
#endif
						g_table_c[m][n] = gf_multp(g_table_c_prev[m][n], discrepancy[l_s]);
					}

					for(n = 0; n < term_num_real; n++)
					{
#if 0
						if((0xFF != g_table_c[l_s][n])
							&& (0xFF != discrepancy[m]))
						{
							DEBUG_NOTICE("delta_l * q_l_s: %d %x %x %x\n", n,
																  gf_multp(g_table_c[l_s][n], discrepancy[m]),
														  		  g_table_c[l_s][n],
														  		  discrepancy[m]);
						}
#endif
						tmp_table_c[n] = gf_multp(g_table_c_prev[l_s][n], discrepancy[m]);
					}

					for(n = 0; n < term_num_real; n++)
					{
#if 0
						if((0xFF != g_table_c[m][n])
							&& (0xFF != tmp_table_c[n]))
						{
							DEBUG_NOTICE("q_l add: %d %x %x %x\n", n,
														  gf_add(g_table_c[m][n], tmp_table_c[n]),
														  g_table_c[m][n],
														  tmp_table_c[n]);
						}
#endif
						g_table_c[m][n] = gf_add(g_table_c[m][n], tmp_table_c[n]);
					}
#if (1 == CFG_DEBUG_NOTICE)
					for(n = 0; n < term_num_real; n++)
					{
						if(0xFF != g_table_c[m][n])
						{
							DEBUG_NOTICE("Q_l: %d | %d %d | %x\n", m, g_table_x[n], g_table_y[n], g_table_c[m][n]);
						}
					}
#endif					
				}
				//DEBUG_NOTICE("update normal Q_l OK\n");

				//DEBUG_NOTICE("g_table_c[0][0]: %x\n", g_table_c[0][0]);
				memset(tmp_table_c, 0xFF, sizeof(unsigned char) * term_num_real);
				for(n = 0; n < term_num_real; n++)//update Q_l_s
				{
					if(0xFF == discrepancy[l_s])
					{
						DEBUG_NOTICE("no update for this (a, b) pair: %d\n", l_s);
						break;
					}
					
					if(0xFF != g_table_c_prev[l_s][n])
					{
						//DEBUG_NOTICE("update Q_l_s: %d %d | %d %d | %x\n", l_s, n, g_table_x[n], g_table_y[n], g_table_c[l_s][n]);
						for(m = 0; m < term_num_real; m++)
						{
							//DEBUG_NOTICE("checking m: %d | %d %d | %x\n", m, g_table_x[m], g_table_y[m], g_table_c[l_s][m]);
							if((g_table_x[m] == (g_table_x[n] + 1))
								&& (g_table_y[m] == g_table_y[n]))
							{
								tmp_table_c[m] = gf_add(g_table_c_prev[l_s][n], tmp_table_c[m]);
								//DEBUG_NOTICE("update tmp_table_c[m]: %d %d | %d %d | %x\n", m, n, g_table_x[m], g_table_y[m], tmp_table_c[m]);
								break;
							}
						}

						//DEBUG_NOTICE("updating g_table_c[l_s][term_use_index[n]]: %x %x\n", g_table_c[l_s][term_use_index[n]], power_polynomial_table[j + 1][0]);
						//g_table_c[l_s][n] = gf_multp(g_table_c[l_s][n], power_polynomial_table[j + 1][0]);
						tmp_table_c[n] = gf_add(gf_multp(g_table_c_prev[l_s][n], power_polynomial_table[j + 1][0]), tmp_table_c[n]);
						//DEBUG_NOTICE("update tmp_table_c[n]: %d %d | %x\n", g_table_x[m], g_table_y[m], tmp_table_c[n]);
						//g_table_c[l_s][n] = gf_add(g_table_c[l_s][n], tmp_table_c[n]);
						//g_table_c[l_s][m] = gf_add(g_table_c[l_s][m], tmp_table_c[m]);
						//DEBUG_NOTICE("update_2 g_table_c[l_s][index]: %d %d %x %x\n", n, m, g_table_c[l_s][n], g_table_c[l_s][m]);
						//break;
					}
				}
				for(n = 0; n < term_num_real; n++)
				{
			#if 0		
					if(0xFF != g_table_c[l_s][n])
					{
						DEBUG_NOTICE("Q_l_s_c: %d | %d %d | %x\n", l_s, g_table_x[n], g_table_y[n], g_table_c[l_s][n]);
					}
			#endif
					if(0xFF == discrepancy[l_s])
					{
						DEBUG_NOTICE("no update for this (a, b) pair: %d\n", l_s);
						break;
					}

					g_table_c[l_s][n] = 0xFF;

					if(0xFF != tmp_table_c[n])
					{
						g_table_c[l_s][n] = gf_multp(tmp_table_c[n], discrepancy[l_s]);
						DEBUG_NOTICE("Q_l_s_c: %d | %d %d | %x\n", l_s, g_table_x[n], g_table_y[n], g_table_c[l_s][n]);
					}
				}
				//DEBUG_NOTICE("update Q_l_s OK\n");

				for(m = 0; m < (d_y_max + 1); m++)
				{
					weight_pol[m] = 0;
					lexorder_pol[m] = 0xFFFF;
					
					for(n = 0; n < term_num_real; n++)
					{
						if(0xFF != g_table_c[m][n])
						{
#if 0							
							if(weight_pol[m] < lex_order_table[g_table_x[n]][g_table_y[n]])
							{
								weight_pol[m] = lex_order_table[g_table_x[n]][g_table_y[n]];
							}
#else
							if((weight_pol[m] < (g_table_x[n] + (MESSAGE_LEN - 1) * g_table_y[n]))
								|| ((weight_pol[m] == (g_table_x[n] + (MESSAGE_LEN - 1) * g_table_y[n])) 
									&& (lexorder_pol[m] > lex_order_table[g_table_x[n]][g_table_y[n]])))
							{
								weight_pol[m] = (g_table_x[n] + (MESSAGE_LEN - 1) * g_table_y[n]);
								lexorder_pol[m] = lex_order_table[g_table_x[n]][g_table_y[n]];
								DEBUG_NOTICE("pol_updating: %d %d\n", weight_pol[m], lexorder_pol[m]);
							}
#endif
						}
					}
					//weight_pol[i] = (MESSAGE_LEN - 1) * i;
					DEBUG_NOTICE("pol: %d %d\n", weight_pol[m], lexorder_pol[m]);
				}
				//weight_pol[l_s] = weight_pol[l_s] + 1;//update w_l_s
				l_w = weight_pol[l_s];
				l_o = lexorder_pol[l_s];
				//DEBUG_NOTICE("update w_l_s OK\n");

				/*store g as prev_poly*/
				for(m = 0; m < (d_y_max + 1); m++)
				{
					for(n = 0; n < term_num_real; n++)
					{
						g_table_c_prev[m][n] = g_table_c[m][n];
#if (1 == CFG_DEBUG_NOTICE)						
						if(0xFF != g_table_c_prev[m][n])
						{
							DEBUG_NOTICE("g_table_c_prev: %d | %d %d | %d %d | %x\n", m, lex_order_table[g_table_x[n]][g_table_y[n]], (g_table_x[n] + (MESSAGE_LEN - 1) * g_table_y[n]), g_table_x[n], g_table_y[n], g_table_c_prev[m][n]);
						}
#endif						
					}
				}
			}

			free(term_use_index);
			term_use_index = NULL;
		}
	}

	/*find the place of smallest polynomial*/
	tmp_real = 0xFF;
	for(m = 0; m < (d_y_max + 1); m++)
	{
		if(tmp_real > weight_pol[m])
		{
			DEBUG_NOTICE("sml_updated: %d | %d %d | %d\n", m, tmp_real, weight_pol[m], lexorder_pol[m]);
			sml_poly = m;
			tmp_real = weight_pol[m];
		}
		if((tmp_real == weight_pol[m])
			&& (lexorder_pol[sml_poly] > lexorder_pol[m]))
		{
			DEBUG_NOTICE("sml_updated: %d | %d %d | %d\n", m, tmp_real, weight_pol[m], lexorder_pol[m]);
			sml_poly = m;
		}
	}

	/*copy the interpolated polynomials to global memory*/
#if 0	
	for(m = 0; m < (d_y_max + 1); m++)
	{
		for(n = 0; n < term_num_real; n++)
		{
			if(0xFF != g_table_c[m][n])
			{
				for(k = 0; k < (TERM_SIZE * TERM_SIZE); k++)
				{
					if((g_table_x[n] == g_term_x[k])
						&& (g_table_y[n] == g_term_y[k]))
					{
						g_term_c[m][k] = g_table_c[m][n];
						DEBUG_NOTICE("g_term: %d | %d %d | %x\n", m,
															g_term_x[k],
															g_term_y[k],
															g_term_c[m][k]);
						break;
					}
				}
			}
		}
	}
#else
	//sml_poly = 5;//test

	for(n = 0; n < term_num_real; n++)
	{
		if(0xFF != g_table_c[sml_poly][n])
		{
			for(k = 0; k < (TERM_SIZE * TERM_SIZE); k++)
			{
				if((g_table_x[n] == g_term_x[k])
					&& (g_table_y[n] == g_term_y[k]))
				{
#if 0
					g_term_c[0][0][k] = g_table_c[sml_poly][n];
					DEBUG_IMPOTANT("g_term: %d | %d %d | %x\n", sml_poly,
														g_term_x[k],
														g_term_y[k],
														g_term_c[0][0][k]);
#else
					g_term_c_p[g_term_phase][0][k] = g_table_c[sml_poly][n];
					DEBUG_IMPOTANT("g_term: %d | %d %d | %x\n", sml_poly,
														g_term_x[k],
														g_term_y[k],
														g_term_c_p[g_term_phase][0][k]);
#endif
					break;
				}
			}
		}
	}
#endif

	if((S_MUL * (CODEWORD_LEN - err_num)) > weight_pol[sml_poly])
	{
		if(0 == decoding_ok_flag)
		{
			decoding_ok_flag = 1;
		}
		weight_stored = weight_pol[sml_poly];
	}

	/*free resources*/
	for (i = 0; i < (d_y_max + 1); i++)
	{
  		free(g_table_c[i]);
		g_table_c[i] = NULL;

		free(g_table_c_prev[i]);
		g_table_c_prev[i] = NULL;
  	}
	free(g_table_c);
	g_table_c = NULL;
	free(g_table_c_prev);
	g_table_c_prev = NULL;
	free(g_table_x);
	g_table_x = NULL;
	free(g_table_y);
	g_table_y = NULL;
	free(discrepancy);
	discrepancy = NULL;
	free(weight_pol);
	weight_pol = NULL;
	free(lexorder_pol);
	lexorder_pol = NULL;
	free(tmp_table_c);
	tmp_table_c = NULL;

	return 0;
}

int g_term_malloc()
{
	unsigned long long i = 0, j = 0;

	g_term_c_p = (unsigned char***)malloc(sizeof(unsigned char**) * LAYER_NUM);
	for (i = 0; i < LAYER_NUM; i++)
	{
		DEBUG_SYS("malloc g_term: %d\n", i);
	
  		g_term_c_p[i] = (unsigned char**)malloc(sizeof(unsigned char*) * POLY_NUM);

		for(j = 0; j < POLY_NUM; j++)
		{
			g_term_c_p[i][j] = (unsigned char*)malloc(sizeof(unsigned char) * (TERM_SIZE * TERM_SIZE));
		}
  	}
	DEBUG_SYS("malloc g_term OK\n", i);
	
	return 0;
}

int g_term_destroy()
{
	unsigned long long i = 0, j = 0;

	for (i = 0; i < LAYER_NUM; i++)
	{
		for(j = 0; j < POLY_NUM; j++)
		{
			free(g_term_c_p[i][j]);
			g_term_c_p[i][j] = NULL;;
		}

		free(g_term_c_p[i]);
		g_term_c_p[i] = NULL;
	}

	free(g_term_c_p);
	g_term_c_p = NULL;
	
	return 0;
}

int g_term_init()
{
	unsigned long long i = 0, j = 0, k = 0;

	for(k = 0; k < LAYER_NUM; k++)
	{
		//for(i = 0; i < POLY_NUM; i++)
		for(j = 0; j < (TERM_SIZE * TERM_SIZE); j++)
		{
#if 0			
			//for(j = 0; j < (TERM_SIZE * TERM_SIZE); j++)
			for(i = 0; i < POLY_NUM; i++)
			{
				g_term_c[k][i][j] = 0xFF;
			}
#endif			
			g_term_0_y_c[k][j] = 0xFF;
			g_term_x_0_c[k][j] = 0xFF;

			for(i = 0; i < POLY_NUM; i++)
			{
				g_term_c_p[k][i][j] = 0xFF;
			}
		}
	}

	k = 0;
	for(i = 0; i < TERM_SIZE; i++)
	{
		for(j = 0; j < TERM_SIZE; j++)
		{
			g_term_x[k] = i;
			g_term_y[k] = j;
			k++;
		}
	}

	return 0;
}

int f_root_init()
{
	unsigned long long i = 0, j = 0;

	for(i = 0; i < ROOT_SIZE; i++)
	{
		for(j = 0; j < ROOT_SIZE; j++)
		{
			f_root_val[i][j] = 0xFF;
			f_root_prev[i][j] = 0xFF;
		}
		f_root_cnt[i + 1] = 0;
	}
	f_root_cnt[0] = 1;
	
	return 0;
}

int g_term_new_gen(unsigned long long layer_idx, unsigned long long tern_idx, unsigned char root_insert)
{
	unsigned long long i = 0, j = 0, k = 0, l = 0, r = 0, s = 0;

	unsigned char tmp = 0xFF;
	unsigned char g_term_c_expand[TERM_SIZE * TERM_SIZE];
	unsigned char tmp_g_term_c_expand[TERM_SIZE * TERM_SIZE];
	unsigned char mul_g_term_c_expand[TERM_SIZE * TERM_SIZE];
	unsigned char g_term_c_expand_store[TERM_SIZE * TERM_SIZE];
	for(k = 0; k < (TERM_SIZE * TERM_SIZE); k++)
	{
#if 0		
		if((1 == g_term_x[k])
			&& (1 == g_term_y[k]))
		{
			g_term_c_expand[k] = 0x0;
		}
		else if((0 == g_term_x[k])
			&& (0 == g_term_y[k]))
		{
			g_term_c_expand[k] = root_insert;
		}
		else
		{
			g_term_c_expand[k] = 0xFF;
		}

		mul_g_term_c_expand[k] = g_term_c_expand[k];
#else
		if((1 == g_term_x[k])
			&& (1 == g_term_y[k]))
		{
			mul_g_term_c_expand[k] = 0x0;
			g_term_c_expand[k] = 0xFF;
			DEBUG_NOTICE("g_term_init: %d %d | %x %x\n",
					g_term_x[k],
					g_term_y[k],
					g_term_c_expand[k],
					mul_g_term_c_expand[k]);
		}
		else if((0 == g_term_x[k])
			&& (0 == g_term_y[k]))
		{
			mul_g_term_c_expand[k] = root_insert;
			g_term_c_expand[k] = 0x0;
			DEBUG_NOTICE("g_term_init: %d %d | %x %x\n",
					g_term_x[k],
					g_term_y[k],
					g_term_c_expand[k],
					mul_g_term_c_expand[k]);
		}
		else
		{
			mul_g_term_c_expand[k] = 0xFF;
			g_term_c_expand[k] = 0xFF;
		}
#endif
		g_term_c_expand_store[k] = 0xFF;
		tmp_g_term_c_expand[k] = 0xFF;
	}

	/*init another g*/
	if((0 == layer_idx)
		&& (0 != tern_idx))
	{
#if 0		
		memcpy(g_term_c[layer_idx][tern_idx],
			   g_term_c[layer_idx][tern_idx - 1],
			   sizeof(unsigned char) * k);
#else
		memcpy(g_term_c_p[g_term_phase][tern_idx],
			   g_term_c_p[g_term_phase][tern_idx - 1],
			   sizeof(unsigned char) * k);
#endif
	}

	for(k = 0; k < (TERM_SIZE * TERM_SIZE); k++)//for every term contain "y"
	{
		if(0 == k)
		{
			DEBUG_NOTICE("g_term_check_err_k_1: %d | %d %d\n",
			k,
			g_term_x[k],
			g_term_y[k]);
			if((0 != g_term_x[k])
				|| (0 != g_term_y[k]))
			{
				while(1);
			}
		}

#if 0
		g_term_c[layer_idx + 1][tern_idx][k] = g_term_c[layer_idx][tern_idx][k];
#else
		g_term_c_p[phase_trans(g_term_phase)][tern_idx][k] = g_term_c_p[g_term_phase][tern_idx][k];
#endif

		if(0 == k)
		{
			DEBUG_NOTICE("g_term_check_err_k_2: %d %d %d | %d | %d %d\n",
			layer_idx,
			tern_idx,
			POLY_NUM,
			k,
			g_term_x[k],
			g_term_y[k]);
			if((0 != g_term_x[k])
				|| (0 != g_term_y[k]))
			{
				while(1);
			}
		}
		
		if((0 != g_term_y[k])
#if 0			
			&& (0xFF != g_term_c[layer_idx][tern_idx][k]))
#else			
			&& (0xFF != g_term_c_p[g_term_phase][tern_idx][k]))
#endif			
		{
			DEBUG_NOTICE("*********************\n");
			DEBUG_NOTICE("g_term_have_y: %d | %d %d | %x\n",
		    tern_idx,
		    g_term_x[k],
		    g_term_y[k],
#if 0
		    g_term_c[layer_idx][tern_idx][k]);
#else
			g_term_c_p[g_term_phase][tern_idx][k]);
#endif
			
			for(l = 0; l < (g_term_y[k] - 0); l++)//for pow_cal, g*m => tmp
			{
				//DEBUG_NOTICE("*********************\n");
				for(i = 0; i < (TERM_SIZE * TERM_SIZE); i++)//for every a
				{
					if(0 == i)
					{
						DEBUG_NOTICE("g_term_check_err_i: %ld | %ld %ld | %x\n",
						i,
						g_term_x[i],
						g_term_y[i],
						g_term_c_expand[i]);
						if((0 != g_term_x[i])
							|| (0 != g_term_y[i]))
						{
							while(1);
						}
					}
					
					if(0xFF != g_term_c_expand[i])
					{
						DEBUG_NOTICE("g_term_c_expand: %ld | %ld %ld | %x\n",
						i,
					    g_term_x[i],
					    g_term_y[i],
					    g_term_c_expand[i]);

						for(j = 0; j < (TERM_SIZE * TERM_SIZE); j++)//for every b
						{
							if(0xFF != mul_g_term_c_expand[j])
							{
								DEBUG_NOTICE("mul_g_term_c_expand: %d %d | %x\n",
							    g_term_x[j],
							    g_term_y[j],
							    mul_g_term_c_expand[j]);

								for(r = 0; r < (TERM_SIZE * TERM_SIZE); r++)//find the place
								{
									if(((g_term_x[i] + g_term_x[j]) == g_term_x[r])
										&& ((g_term_y[i] + g_term_y[j]) == g_term_y[r]))//r <= place(i, j)
									{
										tmp_g_term_c_expand[r] = gf_add(tmp_g_term_c_expand[r], gf_multp(g_term_c_expand[i], mul_g_term_c_expand[j]));
										DEBUG_NOTICE("tmp_g_term_c_expand: %d %d | %x\n",
											     g_term_x[r],
											     g_term_y[r],
											     tmp_g_term_c_expand[r]);
										break;
									}
								}
							}
						}
					}
				}

				/*after every mulp cal, copy and clear*/
				for(r = 0; r < (TERM_SIZE * TERM_SIZE); r++)
				{
					g_term_c_expand[r] = tmp_g_term_c_expand[r];
					tmp_g_term_c_expand[r] = 0xFF;
#if (1 == CFG_DEBUG_NOTICE)
					if(0xFF != g_term_c_expand[r])
					{
						DEBUG_NOTICE("g_term_expand_update: %d %d | %x\n",
							    g_term_x[r],
							    g_term_y[r],
							    g_term_c_expand[r]);
					}
#endif					
				}
			}

			for(r = 0; r < (TERM_SIZE * TERM_SIZE); r++)
			{
				if(0xFF != g_term_c_expand[r])
				{
					for(s = 0; s < (TERM_SIZE * TERM_SIZE); s++)
					{
						if(((g_term_x[r] + g_term_x[k]) == g_term_x[s])
						&& (g_term_y[r] == g_term_y[s]))
						{
							DEBUG_NOTICE("g_term_expand_store updating: %d + %d = %d, %d | %x %x %x\n",
								    g_term_x[r],
								    g_term_x[k],
							    	g_term_x[s],
							    	g_term_y[s],
							    	g_term_c_expand[r],
#if 0
							    	g_term_c[layer_idx][tern_idx][k],
#else
							    	g_term_c_p[g_term_phase][tern_idx][k],
#endif
							    	g_term_c_expand_store[s]);
							//g_term_c_expand_store[s] = gf_add(g_term_c_expand_store[s], g_term_c_expand[r]);
							//g_term_c_expand_store[s] = gf_multp(g_term_c_expand_store[s], g_term_c[layer_idx][tern_idx][k]);
#if 0
							tmp = gf_multp(g_term_c_expand[r], g_term_c[layer_idx][tern_idx][k]);
#else
							tmp = gf_multp(g_term_c_expand[r], g_term_c_p[g_term_phase][tern_idx][k]);
#endif
							g_term_c_expand_store[s] = gf_add(g_term_c_expand_store[s], tmp);
							DEBUG_NOTICE("g_term_expand_store: %d %d | %x\n",
							    	g_term_x[s],
							    	g_term_y[s],
							    	g_term_c_expand_store[s]);
						}
					}
					g_term_c_expand[r] = 0xFF;
				}
			}

			for(r = 0; r < (TERM_SIZE * TERM_SIZE); r++)
			{
#if 0		
				if((1 == g_term_x[r])
					&& (1 == g_term_y[r]))
				{
					g_term_c_expand[r] = 0x0;
				}
				else if((0 == g_term_x[r])
					&& (0 == g_term_y[r]))
				{
					g_term_c_expand[r] = root_insert;
				}
				else
				{
					g_term_c_expand[r] = 0xFF;
				}

				mul_g_term_c_expand[r] = g_term_c_expand[r];
#else
				if((1 == g_term_x[r])
					&& (1 == g_term_y[r]))
				{
					mul_g_term_c_expand[r] = 0x0;
					g_term_c_expand[r] = 0xFF;
				}
				else if((0 == g_term_x[r])
					&& (0 == g_term_y[r]))
				{
					mul_g_term_c_expand[r] = root_insert;
					g_term_c_expand[r] = 0x0;
				}
				else
				{
					mul_g_term_c_expand[r] = 0xFF;
					g_term_c_expand[r] = 0xFF;
				}
#endif
				//g_term_c_expand_store[r] = 0xFF;
				tmp_g_term_c_expand[r] = 0xFF;
			}

			/*clear terms have y*/
#if 0			
			g_term_c[layer_idx + 1][tern_idx][k] = 0xFF;
#else
			g_term_c_p[phase_trans(g_term_phase)][tern_idx][k] = 0xFF;
#endif
		}
	}

	DEBUG_NOTICE("*********************\n");
	/*update g_term*/
	for(r = 0; r < (TERM_SIZE * TERM_SIZE); r++)
	{
#if 0
		if((0 != g_term_y[r])
			&& (0xFF != g_term_c[tern_idx][r]))
		{
			g_term_c[tern_idx][r] = 0xFF;
		}
#endif
#if 0
		g_term_c[layer_idx + 1][tern_idx][r] = gf_add(g_term_c_expand_store[r], g_term_c[layer_idx + 1][tern_idx][r]);
#else
		g_term_c_p[phase_trans(g_term_phase)][tern_idx][r] = gf_add(g_term_c_expand_store[r], g_term_c_p[phase_trans(g_term_phase)][tern_idx][r]);
#endif
#if (1 == CFG_DEBUG_NOTICE)
#if 0
		if(0xFF != g_term_c[layer_idx + 1][tern_idx][r])
#else
		if(0xFF != g_term_c_p[phase_trans(g_term_phase)][tern_idx][r])
#endif
		{
			DEBUG_NOTICE("g_term_update: %d | %d %d | %x\n",
				    tern_idx,
				    g_term_x[r],
				    g_term_y[r],
#if 0
				    g_term_c[layer_idx + 1][tern_idx][r]);
#else
					g_term_c_p[phase_trans(g_term_phase)][tern_idx][r]);
#endif
		}
#endif		
	}

	DEBUG_NOTICE("*********************\n");
	/*eliminate common factor*/
	unsigned char cmn_factor = 0xFF;//deal with x
	for(r = 0; r < (TERM_SIZE * TERM_SIZE); r++)
	{
#if 0		
		if(0xFF != g_term_c[layer_idx + 1][tern_idx][r])
#else
		if(0xFF != g_term_c_p[phase_trans(g_term_phase)][tern_idx][r])
#endif
		{
			if(cmn_factor > g_term_x[r])
			{
				cmn_factor = g_term_x[r];
			}
		}
	}
	DEBUG_NOTICE("cmn_factor for x: %d\n", cmn_factor);
	if(0 != cmn_factor)
	{
		for(r = 0; r < (TERM_SIZE * TERM_SIZE); r++)
		{
#if 0
			if(0xFF != g_term_c[layer_idx + 1][tern_idx][r])
#else
			if(0xFF != g_term_c_p[phase_trans(g_term_phase)][tern_idx][r])
#endif
			{
				for(s = 0; s < (TERM_SIZE * TERM_SIZE); s++)
				{
					if((g_term_x[s] == (g_term_x[r] - cmn_factor))
						&& (g_term_y[s] == g_term_y[r]))
					{
#if 0
						g_term_c[layer_idx + 1][tern_idx][s] = g_term_c[layer_idx + 1][tern_idx][r];
						g_term_c[layer_idx + 1][tern_idx][r] = 0xFF;
						DEBUG_NOTICE("g_term_update_x: %d | %d %d | %x\n",
							    tern_idx,
							    g_term_x[s],
							    g_term_y[s],
							    g_term_c[layer_idx + 1][tern_idx][s]);
#else
						g_term_c_p[phase_trans(g_term_phase)][tern_idx][s] = g_term_c_p[phase_trans(g_term_phase)][tern_idx][r];
						g_term_c_p[phase_trans(g_term_phase)][tern_idx][r] = 0xFF;
						DEBUG_NOTICE("g_term_update_x: %d | %d %d | %x\n",
							    tern_idx,
							    g_term_x[s],
							    g_term_y[s],
							    g_term_c_p[phase_trans(g_term_phase)][tern_idx][s]);
#endif
						break;
					}
				}
			}
		}
	}
#if 0//you cannot deal with y, since y may be 0
	cmn_factor = 0xFF;//deal with y
	for(r = 0; r < (TERM_SIZE * TERM_SIZE); r++)
	{
		if(0xFF != g_term_c[tern_idx][r])
		{
			if(cmn_factor > g_term_y[r])
			{
				cmn_factor = g_term_y[r];
			}
		}
	}
	DEBUG_NOTICE("cmn_factor for y: %d\n", cmn_factor);
	if(0 != cmn_factor)
	{
		for(r = 0; r < (TERM_SIZE * TERM_SIZE); r++)
		{
			if(0xFF != g_term_c[tern_idx][r])
			{
				for(s = 0; s < (TERM_SIZE * TERM_SIZE); s++)
				{
					if((g_term_y[s] == (g_term_y[r] - cmn_factor))
						&& (g_term_x[s] == g_term_x[r]))
					{
						g_term_c[tern_idx][s] = g_term_c[tern_idx][r];
						g_term_c[tern_idx][r] = 0xFF;
						DEBUG_NOTICE("g_term_update_y: %d | %d %d | %x\n",
							    tern_idx,
							    g_term_x[s],
							    g_term_y[s],
							    g_term_c[tern_idx][s]);
						break;
					}
				}
			}
		}
	}
#endif

#if (1 == CFG_DEBUG_NOTICE)
	DEBUG_NOTICE("*********************\n");
	for(r = 0; r < (TERM_SIZE * TERM_SIZE); r++)
	{
#if 0
		if(0xFF != g_term_c[layer_idx + 1][tern_idx][r])
		{
			DEBUG_INFO("g_term_update_OK: %d | %d %d | %x\n",
				    tern_idx,
				    g_term_x[r],
				    g_term_y[r],
				    g_term_c[layer_idx + 1][tern_idx][r]);
		}
#else
		if(0xFF != g_term_c_p[phase_trans(g_term_phase)][tern_idx][r])
		{
			DEBUG_INFO("g_term_update_OK: %d | %d %d | %x\n",
				    tern_idx,
				    g_term_x[r],
				    g_term_y[r],
				    g_term_c_p[phase_trans(g_term_phase)][tern_idx][r]);
		}
#endif
	}
#endif
	
	return 0;
}

int g_term_0_y_cal(unsigned long long layer_idx, unsigned long long tern_idx)
{
	unsigned long long i = 0, j = 0;

#if 0
	/*clear*/
	for(i = 0; i < POLY_NUM; i++)
	{
		for(j = 0; j < (TERM_SIZE * TERM_SIZE); j++)
		{
			g_term_0_y_c[i][j] = 0xFF;
		}
	}

	for(i = 0; i < POLY_NUM; i++)
	{
		for(j = 0; j < (TERM_SIZE * TERM_SIZE); j++)
		{
			if((0 != g_term_x[j])
				&& (0xFF != g_term_c[i][j]))
			{
				g_term_0_y_c[i][j] = 0xFF;
				DEBUG_NOTICE("g_term_set_to_zero: %d | %d %d | %x %x\n",
					    i,
					    g_term_x[j],
					    g_term_y[j],
					    g_term_c[i][j],
					    g_term_0_y_c[i][j]);
			}
			else
			{
				g_term_0_y_c[i][j] = g_term_c[i][j];
				if(0xFF != g_term_0_y_c[i][j])
				{
					DEBUG_NOTICE("g_term_0_y: %d | %d %d | %x\n",
							i,
							g_term_x[j],
							g_term_y[j],
							g_term_0_y_c[i][j]);
				}
			}
		}
	}
#else
	/*clear*/
	for(j = 0; j < (TERM_SIZE * TERM_SIZE); j++)
	{
		//DEBUG_INFO("layer_idx: %ld, j: %ld\n", layer_idx, j);
		g_term_0_y_c[layer_idx][j] = 0xFF;
	}

	for(j = 0; j < (TERM_SIZE * TERM_SIZE); j++)
	{
		if(0 == j)
		{
			DEBUG_NOTICE("g_term_check_err_j: %d | %d %d\n",
			j,
			g_term_x[j],
			g_term_y[j]);
			if((0 != g_term_x[j])
				|| (0 != g_term_y[j]))
			{
				while(1);
			}
		}
		
		if((0 != g_term_x[j])
#if 0
			&& (0xFF != g_term_c[layer_idx][tern_idx][j]))
#else
			&& (0xFF != g_term_c_p[g_term_phase][tern_idx][j]))
#endif
		{
			g_term_0_y_c[layer_idx][j] = 0xFF;
			DEBUG_NOTICE("g_term_set_to_zero: %d | %d %d | %x %x\n",
				    i,
				    g_term_x[j],
				    g_term_y[j],
#if 0
				    g_term_c[layer_idx][tern_idx][j],
#else
					g_term_c_p[g_term_phase][tern_idx][j],
#endif
				    g_term_0_y_c[layer_idx][j]);
		}
		else
		{
#if 0
			g_term_0_y_c[layer_idx][j] = g_term_c[layer_idx][tern_idx][j];
#else
			g_term_0_y_c[layer_idx][j] = g_term_c_p[g_term_phase][tern_idx][j];
#endif
#if (1 == CFG_DEBUG_NOTICE)
			if(0xFF != g_term_0_y_c[layer_idx][j])
			{
				DEBUG_NOTICE("g_term_0_y: %d | %d %d | %x\n",
						i,
						g_term_x[j],
						g_term_y[j],
						g_term_0_y_c[layer_idx][j]);
			}
#endif
		}
	}
#endif

	return 0;
}

unsigned char g_term_x_0_cal(unsigned long long layer_idx, unsigned long long tern_idx)
{
	unsigned char val = 0xFF;

	unsigned long long i = 0, j = 0;

	/*clear*/
	for(j = 0; j < (TERM_SIZE * TERM_SIZE); j++)
	{
		g_term_x_0_c[layer_idx][j] = 0xFF;
	}

	for(j = 0; j < (TERM_SIZE * TERM_SIZE); j++)
	{
		if((0 != g_term_y[j])
#if 0
			&& (0xFF != g_term_c[layer_idx][tern_idx][j]))
#else
			&& (0xFF != g_term_c_p[g_term_phase][tern_idx][j]))
#endif
		{
			g_term_x_0_c[layer_idx][j] = 0xFF;
#if 0			
			DEBUG_NOTICE("g_term_set_to_zero: %d | %d %d | %x %x\n",
				    i,
				    g_term_x[j],
				    g_term_y[j],
				    g_term_c[layer_idx][tern_idx][j],
				    g_term_x_0_c[layer_idx][tern_idx][j]);
#endif
		}
		else
		{
#if 0
			g_term_x_0_c[layer_idx][j] = g_term_c[layer_idx][tern_idx][j];
#else
			g_term_x_0_c[layer_idx][j] = g_term_c_p[g_term_phase][tern_idx][j];
#endif
#if 0
			if(0xFF != g_term_x_0_c[layer_idx][tern_idx][j])
			{
				DEBUG_NOTICE("g_term_0_y: %d | %d %d | %x\n",
						i,
						g_term_x[j],
						g_term_y[j],
						g_term_x_0_c[layer_idx][tern_idx][j]);
			}
#endif			
		}

		val = gf_add(val, g_term_x_0_c[layer_idx][j]);
	}
	
	return val;
}

int chien_searching_for_g_0_y(unsigned long long layer_idx, unsigned long long tern_idx, unsigned char root_test)
{
	long long is_root = 0;
	unsigned long long k = 0;
	unsigned char tmp = 0xFF, tmp_sum = 0xFF;
	
	for(k = 0; k < (TERM_SIZE * TERM_SIZE); k++)
	{
		if(0xFF != g_term_0_y_c[layer_idx][k])
		{
			tmp = 0xFF;
			tmp = gf_pow_cal(root_test, g_term_y[k]);
			//DEBUG_NOTICE("tmp: %x | %x^%d\n", tmp, root_test, g_term_y[k]);
			//DEBUG_NOTICE("tmp: %x | %x*%x\n", gf_multp(tmp, g_term_0_y_c[tern_idx][k]), tmp, g_term_0_y_c[tern_idx][k]);
			tmp = gf_multp(tmp, g_term_0_y_c[layer_idx][k]);
			//DEBUG_NOTICE("tmp_sum: %x | %x+%x\n", gf_add(tmp_sum, tmp), tmp_sum, tmp);
			tmp_sum = gf_add(tmp_sum, tmp);
		}
	}

	if(0xFF == tmp_sum)
	{
		is_root = 1;
		DEBUG_INFO("is_root: %d %x\n", tern_idx, root_test);
	}
	else
	{
		is_root = 0;
	}

	return is_root;
}

#if 0
int test_factorization_init()
{
	int k = 0;

	for(k = 0; k < (TERM_SIZE * TERM_SIZE); k++)
	{
		g_term_c[0][0][k] = 0xFF;

		if((0 == g_term_x[k])
			&& (7 == g_term_y[k]))
		{
			g_term_c[0][0][k] = 2;
		}

		if((3 == g_term_x[k])
			&& (6 == g_term_y[k]))
		{
			g_term_c[0][0][k] = 5;
		}
		if((2 == g_term_x[k])
			&& (6 == g_term_y[k]))
		{
			g_term_c[0][0][k] = 1;
		}
		if((1 == g_term_x[k])
			&& (6 == g_term_y[k]))
		{
			g_term_c[0][0][k] = 6;
		}
		if((0 == g_term_x[k])
			&& (6 == g_term_y[k]))
		{
			g_term_c[0][0][k] = 6;
		}

		if((5 == g_term_x[k])
			&& (5 == g_term_y[k]))
		{
			g_term_c[0][0][k] = 1;
		}
		if((4 == g_term_x[k])
			&& (5 == g_term_y[k]))
		{
			g_term_c[0][0][k] = 6;
		}
		if((3 == g_term_x[k])
			&& (5 == g_term_y[k]))
		{
			g_term_c[0][0][k] = 4;
		}
		if((2 == g_term_x[k])
			&& (5 == g_term_y[k]))
		{
			g_term_c[0][0][k] = 4;
		}
		if((0 == g_term_x[k])
			&& (5 == g_term_y[k]))
		{
			g_term_c[0][0][k] = 2;
		}

		if((7 == g_term_x[k])
			&& (4 == g_term_y[k]))
		{
			g_term_c[0][0][k] = 0;
		}
		if((6 == g_term_x[k])
			&& (4 == g_term_y[k]))
		{
			g_term_c[0][0][k] = 3;
		}
		if((5 == g_term_x[k])
			&& (4 == g_term_y[k]))
		{
			g_term_c[0][0][k] = 1;
		}
		if((4 == g_term_x[k])
			&& (4 == g_term_y[k]))
		{
			g_term_c[0][0][k] = 0;
		}
		if((3 == g_term_x[k])
			&& (4 == g_term_y[k]))
		{
			g_term_c[0][0][k] = 2;
		}
		if((1 == g_term_x[k])
			&& (4 == g_term_y[k]))
		{
			g_term_c[0][0][k] = 5;
		}
		if((0 == g_term_x[k])
			&& (4 == g_term_y[k]))
		{
			g_term_c[0][0][k] = 2;
		}

		if((9 == g_term_x[k])
			&& (3 == g_term_y[k]))
		{
			g_term_c[0][0][k] = 1;
		}
		if((8 == g_term_x[k])
			&& (3 == g_term_y[k]))
		{
			g_term_c[0][0][k] = 6;
		}
		if((7 == g_term_x[k])
			&& (3 == g_term_y[k]))
		{
			g_term_c[0][0][k] = 6;
		}
		if((6 == g_term_x[k])
			&& (3 == g_term_y[k]))
		{
			g_term_c[0][0][k] = 4;
		}
		if((5 == g_term_x[k])
			&& (3 == g_term_y[k]))
		{
			g_term_c[0][0][k] = 6;
		}
		if((4 == g_term_x[k])
			&& (3 == g_term_y[k]))
		{
			g_term_c[0][0][k] = 2;
		}
		if((3 == g_term_x[k])
			&& (3 == g_term_y[k]))
		{
			g_term_c[0][0][k] = 3;
		}
		if((2 == g_term_x[k])
			&& (3 == g_term_y[k]))
		{
			g_term_c[0][0][k] = 3;
		}
		if((1 == g_term_x[k])
			&& (3 == g_term_y[k]))
		{
			g_term_c[0][0][k] = 3;
		}
		if((0 == g_term_x[k])
			&& (3 == g_term_y[k]))
		{
			g_term_c[0][0][k] = 2;
		}

		if((11 == g_term_x[k])
			&& (2 == g_term_y[k]))
		{
			g_term_c[0][0][k] = 2;
		}
		if((10 == g_term_x[k])
			&& (2 == g_term_y[k]))
		{
			g_term_c[0][0][k] = 4;
		}
		if((9 == g_term_x[k])
			&& (2 == g_term_y[k]))
		{
			g_term_c[0][0][k] = 4;
		}
		if((8 == g_term_x[k])
			&& (2 == g_term_y[k]))
		{
			g_term_c[0][0][k] = 1;
		}
		if((7 == g_term_x[k])
			&& (2 == g_term_y[k]))
		{
			g_term_c[0][0][k] = 5;
		}
		if((6 == g_term_x[k])
			&& (2 == g_term_y[k]))
		{
			g_term_c[0][0][k] = 4;
		}
		if((5 == g_term_x[k])
			&& (2 == g_term_y[k]))
		{
			g_term_c[0][0][k] = 0;
		}
		if((4 == g_term_x[k])
			&& (2 == g_term_y[k]))
		{
			g_term_c[0][0][k] = 2;
		}
		if((3 == g_term_x[k])
			&& (2 == g_term_y[k]))
		{
			g_term_c[0][0][k] = 0;
		}
		if((2 == g_term_x[k])
			&& (2 == g_term_y[k]))
		{
			g_term_c[0][0][k] = 4;
		}
		if((1 == g_term_x[k])
			&& (2 == g_term_y[k]))
		{
			g_term_c[0][0][k] = 0;
		}
		if((0 == g_term_x[k])
			&& (2 == g_term_y[k]))
		{
			g_term_c[0][0][k] = 3;
		}

		if((13 == g_term_x[k])
			&& (1 == g_term_y[k]))
		{
			g_term_c[0][0][k] = 1;
		}
		if((12 == g_term_x[k])
			&& (1 == g_term_y[k]))
		{
			g_term_c[0][0][k] = 6;
		}
		if((11 == g_term_x[k])
			&& (1 == g_term_y[k]))
		{
			g_term_c[0][0][k] = 6;
		}
		if((9 == g_term_x[k])
			&& (1 == g_term_y[k]))
		{
			g_term_c[0][0][k] = 0;
		}
		if((8 == g_term_x[k])
			&& (1 == g_term_y[k]))
		{
			g_term_c[0][0][k] = 4;
		}
		if((7 == g_term_x[k])
			&& (1 == g_term_y[k]))
		{
			g_term_c[0][0][k] = 5;
		}
		if((5 == g_term_x[k])
			&& (1 == g_term_y[k]))
		{
			g_term_c[0][0][k] = 2;
		}
		if((4 == g_term_x[k])
			&& (1 == g_term_y[k]))
		{
			g_term_c[0][0][k] = 6;
		}
		if((3 == g_term_x[k])
			&& (1 == g_term_y[k]))
		{
			g_term_c[0][0][k] = 5;
		}
		if((1 == g_term_x[k])
			&& (1 == g_term_y[k]))
		{
			g_term_c[0][0][k] = 4;
		}
		if((0 == g_term_x[k])
			&& (1 == g_term_y[k]))
		{
			g_term_c[0][0][k] = 4;
		}

		if((14 == g_term_x[k])
			&& (0 == g_term_y[k]))
		{
			g_term_c[0][0][k] = 6;
		}
		if((13 == g_term_x[k])
			&& (0 == g_term_y[k]))
		{
			g_term_c[0][0][k] = 0;
		}
		if((12 == g_term_x[k])
			&& (0 == g_term_y[k]))
		{
			g_term_c[0][0][k] = 3;
		}
		if((11 == g_term_x[k])
			&& (0 == g_term_y[k]))
		{
			g_term_c[0][0][k] = 6;
		}
		if((10 == g_term_x[k])
			&& (0 == g_term_y[k]))
		{
			g_term_c[0][0][k] = 3;
		}
		if((9 == g_term_x[k])
			&& (0 == g_term_y[k]))
		{
			g_term_c[0][0][k] = 1;
		}
		if((8 == g_term_x[k])
			&& (0 == g_term_y[k]))
		{
			g_term_c[0][0][k] = 5;
		}
		if((7 == g_term_x[k])
			&& (0 == g_term_y[k]))
		{
			g_term_c[0][0][k] = 4;
		}
		if((6 == g_term_x[k])
			&& (0 == g_term_y[k]))
		{
			g_term_c[0][0][k] = 1;
		}
		if((4 == g_term_x[k])
			&& (0 == g_term_y[k]))
		{
			g_term_c[0][0][k] = 4;
		}
		if((2 == g_term_x[k])
			&& (0 == g_term_y[k]))
		{
			g_term_c[0][0][k] = 0;
		}
		if((1 == g_term_x[k])
			&& (0 == g_term_y[k]))
		{
			g_term_c[0][0][k] = 2;
		}
	}
	
	return 0;
}
#endif

int rr_factorization()
{
	unsigned long long i = 0, j = 0, k = 0, l = 0, r = 0, s = 0;
	long long is_root = 0;

	//test_factorization_init();
	f_root_init();

	//chien_searching_for_g_0_y(1, power_polynomial_table[0][0]);//test
	//g_term_new_gen(1, power_polynomial_table[0][0]);//test

	for(s = 0; s < MESSAGE_LEN; s++)//layer
	{
		DEBUG_INFO("-------------------------------------------\n");
		//DEBUG_INFO("Layer: %ld, g_term_phase: %ld\n", s, g_term_phase);
		l = 0;
		for(r = 0; r < f_root_cnt[s]; r++)//roots per layer
		{
			//DEBUG_INFO("Root: %ld, f_root_cnt: %ld\n", r, f_root_cnt[s]);
			g_term_0_y_cal(s, r);

			for(k = 0; k < GF_FIELD; k++)//test roots
			{
				is_root = chien_searching_for_g_0_y(s, r, power_polynomial_table[k][0]);
				if(1 == is_root)
				{
					f_root_val[s][l] = power_polynomial_table[k][0];
					f_root_prev[s][l] = r;
					f_root_cnt[s + 1] = f_root_cnt[s + 1] + 1;

					DEBUG_INFO("f_root: %d %d | %x | %x | %d\n",
								s,
								l,
								f_root_val[s][l],
								f_root_prev[s][l],
								f_root_cnt[s + 1]);

					g_term_new_gen(s, l, f_root_val[s][l]);

					l = l + 1;
					if(POLY_NUM <= l)
					{
						/*just exit*/
						k = GF_FIELD;
						r = f_root_cnt[s];
					}
				}
			}
		}

		/*phase trans.*/
		if(0 == g_term_phase)
		{
			g_term_phase = 1;
		}
		else
		{
			g_term_phase = 0;
		}
	}

	return 0;
}

float euc_distance_code_cal(unsigned char *a,
								   float **b,
								   unsigned long long len)
{
	float euc_distance = 0;

	unsigned long long i = 0;
	unsigned long long symbol_num = CODEWORD_LEN * GF_Q * BITS_PER_SYMBOL_BPSK;
	float **a_mod;
	a_mod = (float**)malloc(sizeof(float*) * symbol_num);
	for (i = 0; i < symbol_num; i++)
	{
		a_mod[i] = (float*)malloc(sizeof(float) * 2);
	}

	bpsk_mod(a,
			 len,
			 a_mod,
			 symbol_num);

	for(i = 0; i < symbol_num; i++)
	{
		euc_distance = euc_distance + abs(a_mod[i][0] - b[i][0])
									+ abs(a_mod[i][1] - b[i][1]);
	}

	for (i = 0; i < symbol_num; i++)
	{
  		free(a_mod[i]);
		a_mod[i] = NULL;
  	}
	free(a_mod);
	a_mod = NULL;

	//DEBUG_IMPOTANT("Euc. Distance Test: %d\n", euc_distance);

	return euc_distance;
}

unsigned long long hamm_distance_code_cal(unsigned char *a,
									  unsigned char *b,
									  unsigned long long len)
{
	unsigned long long hamm_distance = 0;

	long long i = 0;

	for(i = 0; i < len; i++)
	{
		if(a[i] != b[i])
		{
			hamm_distance = hamm_distance + 1;
		}
	}

	//DEBUG_IMPOTANT("Hamming Distance Test: %d\n", hamm_distance);

	return hamm_distance;
}

unsigned long long hamm_distance_bit_cal(unsigned char *a,
									  unsigned char *b,
									  unsigned long long len)
{
	unsigned long long hamm_distance_bit = 0;

	long long i = 0, j = 0;
	unsigned char tmp_a = 0, tmp_b = 0;

	for(i = 0; i < len; i++)
	{
		if(a[i] != b[i])
		{
			for(j = 0; j < GF_Q; j++)
			{
				/*notice these indexes, 0xFF + 0x1 = 0, ..., 0x6 + 0x1 = 7*/
				tmp_a = (power_polynomial_table[a[i] + 0x1][1] >> j) & 0x1;
				tmp_b = (power_polynomial_table[b[i] + 0x1][1] >> j) & 0x1;

				if(tmp_a != tmp_b)
				{
					hamm_distance_bit = hamm_distance_bit + 1;
				}
			}
		}
	}

	//DEBUG_IMPOTANT("Hamming Distance Test: %d\n", hamm_distance);

	return hamm_distance_bit;
}

int check_rr_decoded_result()
{
	long long r = 0, s = 0, prev_r = 0;
	unsigned long long hamm_distance_code = 0xFFFF, hamm_distance_bit = 0xFFFF, tmp_code = 0, tmp_bit = 0;
	unsigned char tmp_decoded_message[MESSAGE_LEN];
	unsigned char tmp_decoded_codeword[CODEWORD_LEN];
	//unsigned long long hamm_distance_debug = 0xFFFF;

	FILE *frc;

	memset(tmp_decoded_message, 0xFF, sizeof(unsigned char) * MESSAGE_LEN);
	memset(tmp_decoded_codeword, 0xFF, sizeof(unsigned char) * CODEWORD_LEN);

	for(r = 0; r < f_root_cnt[MESSAGE_LEN]; r++)
	{
		if(0xFF == g_term_x_0_cal(MESSAGE_LEN, r))
		{
			DEBUG_IMPOTANT("Decoding OK!\n");
			if(1 == decoding_ok_flag)
			{
				decoding_ok_flag = 2;
			}

			DEBUG_IMPOTANT("Message: %d %d | %x\n", MESSAGE_LEN - 1, r, f_root_val[MESSAGE_LEN - 1][r]);

			memset(tmp_decoded_message, 0xFF, sizeof(unsigned char) * MESSAGE_LEN);
			tmp_decoded_message[MESSAGE_LEN - 1] = f_root_val[MESSAGE_LEN - 1][r];

			prev_r = r;
			for(s = MESSAGE_LEN - 2; s >= 0; s--)
			{
				DEBUG_IMPOTANT("Message: %d %d | %x\n", s, prev_r, f_root_val[s][f_root_prev[s + 1][prev_r]]);

				tmp_decoded_message[s] = f_root_val[s][f_root_prev[s + 1][prev_r]];

				prev_r = f_root_prev[s + 1][prev_r];
			}

			evaluation_encoding_v2(tmp_decoded_message, tmp_decoded_codeword);

			tmp_code = hamm_distance_code_cal(tmp_decoded_codeword, received_polynomial, CODEWORD_LEN);
			//tmp_code = (unsigned long long)euc_distance_code_cal(tmp_decoded_codeword, recv_seq, CODEWORD_LEN);
			tmp_bit = hamm_distance_bit_cal(tmp_decoded_codeword, received_polynomial, CODEWORD_LEN);
			if(tmp_code < hamm_distance_code)
			{
				hamm_distance_code = tmp_code;
				hamm_distance_bit = tmp_bit;
				memcpy(decoded_codeword, tmp_decoded_codeword, sizeof(unsigned char) * CODEWORD_LEN);
#if (1 == SYS_ENC)				
				memcpy(decoded_message, (tmp_decoded_codeword + (CODEWORD_LEN - MESSAGE_LEN)), sizeof(unsigned char) * MESSAGE_LEN);
#else
				memcpy(decoded_message, tmp_decoded_message, sizeof(unsigned char) * MESSAGE_LEN);
#endif
			}
			if(tmp_code == hamm_distance_code)
			{
				if(tmp_bit < hamm_distance_bit)
				{
					hamm_distance_bit = tmp_bit;
					memcpy(decoded_codeword, tmp_decoded_codeword, sizeof(unsigned char) * CODEWORD_LEN);
#if (1 == SYS_ENC)				
					memcpy(decoded_message, (tmp_decoded_codeword + (CODEWORD_LEN - MESSAGE_LEN)), sizeof(unsigned char) * MESSAGE_LEN);
#else
					memcpy(decoded_message, tmp_decoded_message, sizeof(unsigned char) * MESSAGE_LEN);
#endif
				}
			}

			if(0 != hamm_distance_debug)
			{
				hamm_distance_debug = hamm_distance_code_cal(tmp_decoded_message, message_polynomial, MESSAGE_LEN);		
				DEBUG_IMPOTANT("Hamming Distance Debug: %ld\n", hamm_distance_debug);
			}
		}
		else
		{
			DEBUG_IMPOTANT("Decoding Fail for Root-%d\n", r);
		}
	}

#if (1 == CFG_DEBUG_IMPOTANT)
	if(0xFFFF != hamm_distance_code)
	{
		DEBUG_IMPOTANT("Received Codeword:\n");
		for(r = 0; r < CODEWORD_LEN; r++)
		{
			DEBUG_IMPOTANT("%x ", received_polynomial[r]);
		}
		DEBUG_IMPOTANT("\n");

		DEBUG_IMPOTANT("Decoding Result:\n");
		for(r = 0; r < CODEWORD_LEN; r++)
		{
			DEBUG_IMPOTANT("%x ", decoded_codeword[r]);
		}
		DEBUG_IMPOTANT("\n");

		DEBUG_IMPOTANT("Hamming Distance: %ld %ld\n", hamm_distance_code, hamm_distance_bit);

		DEBUG_IMPOTANT("Decoded Message:\n");
		for(r = 0; r < MESSAGE_LEN; r++)
		{
			DEBUG_IMPOTANT("%x ", decoded_message[r]);
		}
		DEBUG_IMPOTANT("\n");
	}
	else
	{
		DEBUG_IMPOTANT("Decoding Fail!\n");
	}
#endif

	return 0;
}

int as_decoding()
{
#if (1 == CFG_DEBUG_IMPOTANT)	
	unsigned long long i = 0;

	DEBUG_IMPOTANT("Received:\n");
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		DEBUG_IMPOTANT("%x ", received_polynomial[i]);
	}
	DEBUG_IMPOTANT("\n");
#endif

	memset(decoded_codeword, 0xFF, sizeof(unsigned char) * CODEWORD_LEN);
	memset(decoded_message, 0xFF, sizeof(unsigned char) * MESSAGE_LEN);

	g_term_init();

	koetter_interpolation();

	rr_factorization();

	check_rr_decoded_result();

	return 0;
}
