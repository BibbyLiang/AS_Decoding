#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "gf_cal.h"
#include "encoding.h"
#include "as_decoding.h"

#define S_MUL		4
#define K_M			S_MUL * (S_MUL - 1) * MESSAGE_LEN
#define LAYER_NUM	5
#define POLY_NUM	100
#define TERM_SIZE	100//REAL_SIZE = TERM_SIZE^2
#define ROOT_SIZE	100//REAL_SIZE = TERM_SIZE^2

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
unsigned char g_term_c[LAYER_NUM][POLY_NUM][TERM_SIZE * TERM_SIZE];
unsigned char g_term_x[TERM_SIZE * TERM_SIZE];
unsigned char g_term_y[TERM_SIZE * TERM_SIZE];
unsigned char g_term_0_y_c[LAYER_NUM][POLY_NUM][TERM_SIZE * TERM_SIZE];
unsigned char g_term_x_0_c[LAYER_NUM][POLY_NUM][TERM_SIZE * TERM_SIZE];
unsigned char f_root_val[ROOT_SIZE][ROOT_SIZE];
unsigned char f_root_prev[ROOT_SIZE][ROOT_SIZE];
unsigned char f_root_cnt[ROOT_SIZE + 1];//used for next layer
unsigned char sml_poly = 0xFF;

void find_max_val(float matrix[][CODEWORD_LEN], unsigned char col,
					 unsigned char* m_ptr, unsigned char* n_ptr)
{
	unsigned char i = 0, j = 0;
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

	//printf("i: %d, j: %d, val: %f\n", *m_ptr, *n_ptr, max_val);

	return;
}

int chnl_rel_init()
{
	unsigned char i = 0, j = 0;

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

int mul_assign()
{
	unsigned char i = 0, j = 0;
#if 0	
	unsigned char s = 0;
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

	printf("Multiplicity Assignment:\n");
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		for(j = 0; j < CODEWORD_LEN + 1; j++)
		{
			printf("%d ", mul_matrix[j][i]);
		}
		printf("\n");
	}

	free(m_ptr);
	free(n_ptr);
#else//GS algorithm
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		for(j = 0; j < CODEWORD_LEN + 1; j++)
		{
			if(power_polynomial_table[j][0] == received_polynomial[i])
			{
				mul_matrix[j][i] = S_MUL;
			}
		}
	}

	printf("Multiplicity Assignment:\n");
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		for(j = 0; j < CODEWORD_LEN + 1; j++)
		{
			printf("%d ", mul_matrix[j][i]);
		}
		printf("\n");
	}
#endif
	
	return 0;
}

unsigned char syndrome_cal(unsigned char *recv, unsigned char *synd,
								unsigned int cw_len, unsigned int msg_len)
{
	unsigned int i = 0, j = 0;
	unsigned char tmp = 0xFF, tmp_sum = 0xFF;

	for(i = 0; i < cw_len - msg_len; i++)
	{
		tmp = 0xFF;
		tmp_sum = 0xFF;
		for(j = 0; j < cw_len; j++)
		{
			tmp = gf_multp(recv[j], (i + 1) * j);
			tmp_sum = gf_add(tmp, tmp_sum);
			//printf("%x %x\n", tmp, tmp_sum);
		}
		synd[i] = tmp_sum;
	}
	printf("Syndrome Polynomial:\n");
	for(i = 0; i < cw_len - msg_len; i++)
	{
		printf("%x ", synd[i]);
	}
	printf("\n");
}

int rel_group()
{
	unsigned int i = 0, j = 0;

	float rel_thrd = 0.7;
	unsigned int rel_cnt = 0, rel_flag = 0;

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
					//printf("rel_val: %d %d %f %f\n", i, j, rel_thrd, chnl_rel_matrix[i][j]);
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

		//printf("rel_group: %f %d\n", rel_thrd, rel_cnt);

		if(MESSAGE_LEN < rel_cnt)
		{
			rel_thrd = rel_thrd + 0.05;
		}
		else
		{
			rel_thrd = rel_thrd - 0.05;
		}
	}

	printf("rel_seq: ");
	for(i = 0; i < MESSAGE_LEN; i++)
	{
		printf("%d ", rel_group_seq[i]);
	}
	printf("\n");
	printf("unrel_seq: ");
	for(i = 0; i < CODEWORD_LEN - MESSAGE_LEN; i++)
	{
		printf("%d ", unrel_group_seq[i]);
	}
	printf("\n");
	
	return 0;
}

int tao_cal()
{
	unsigned int i = 0, j = 0;
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
	printf("tao: ");
	for(j = 0; j < (CODEWORD_LEN - MESSAGE_LEN + 1); j++)
	{
		printf("%x ", reg[j]);
	}
	printf("\n");
	
	return 0;
}

int sigma_cal()
{
	unsigned int i = 0, j = 0;
	memset(sigma, 0xFF, sizeof(unsigned char) * (((CODEWORD_LEN - MESSAGE_LEN) + (CODEWORD_LEN - MESSAGE_LEN + 1) - 1) - (CODEWORD_LEN - MESSAGE_LEN)));

	unsigned char tmp[(CODEWORD_LEN - MESSAGE_LEN) + (CODEWORD_LEN - MESSAGE_LEN + 1) - 1];
	memset(tmp, 0xFF, sizeof(unsigned char) * ((CODEWORD_LEN - MESSAGE_LEN) + (CODEWORD_LEN - MESSAGE_LEN + 1) - 1));

	gf_multp_poly_hw(syndrome, (CODEWORD_LEN - MESSAGE_LEN),
					  tao, 		  (CODEWORD_LEN - MESSAGE_LEN + 1),
					  tmp, 		  ((CODEWORD_LEN - MESSAGE_LEN) + (CODEWORD_LEN - MESSAGE_LEN + 1) - 1));

	memcpy(sigma, tmp + (CODEWORD_LEN - MESSAGE_LEN), sizeof(unsigned char) * (((CODEWORD_LEN - MESSAGE_LEN) + (CODEWORD_LEN - MESSAGE_LEN + 1) - 1) - (CODEWORD_LEN - MESSAGE_LEN)));

	printf("sigma: ");
	for(i = 0; i < (((CODEWORD_LEN - MESSAGE_LEN) + (CODEWORD_LEN - MESSAGE_LEN + 1) - 1) - (CODEWORD_LEN - MESSAGE_LEN)); i++)
	{
		printf("%x ", sigma[i]);
	}
	printf("\n");

	memcpy(omega, tmp, sizeof(unsigned char) * (CODEWORD_LEN - MESSAGE_LEN));
	printf("omega: ");
	for(i = 0; i < (CODEWORD_LEN - MESSAGE_LEN); i++)
	{
		printf("%x ", omega[i]);
	}
	printf("\n");

	return 0;
}

unsigned char poly_eva(unsigned char *poly, unsigned char poly_len, unsigned char input_val)
{
	unsigned int i = 0;
	unsigned char poly_val = 0xFF, tmp_product = 0;

	for(i = 0; i < poly_len; i++)
	{
		tmp_product = gf_multp(*(poly + i), (input_val * i) % (GF_FIELD - 1));
		//printf("poly_eva: %d %x %x %x %x %x\n", i, *(poly + i), input_val, (input_val * i) % (GF_FIELD - 1), tmp_product, poly_val);
		poly_val = gf_add(tmp_product, poly_val);
	}
	
	return poly_val;
}

int phi_cal()
{
	unsigned int i = 0, k = 0;
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
	printf("tao_dev: ");
	for(i = 0; i < (CODEWORD_LEN - MESSAGE_LEN); i++)
	{
		printf("%x ", tao_dev[i]);
	}
	printf("\n");

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
			//printf("%x: %x\n", i, locator);
			tmp = poly_eva(omega, CODEWORD_LEN - MESSAGE_LEN, locator);
			//printf("%x ", tmp);
			//tmp = gf_div(tmp, (locator * (CODEWORD_LEN - MESSAGE_LEN + 1)) % (GF_FIELD - 1));
			//printf("%x ", tmp);
			tmp = gf_div(tmp, poly_eva(tao_dev, CODEWORD_LEN - MESSAGE_LEN, locator));
			//printf("%x ", tmp);
			//phi[i] = gf_add(tmp, received_polynomial[i]);
			phi[i] = tmp;
			//printf("\n");
		}
		else
		{
			phi[i] = received_polynomial[i];
		}
	}

	printf("phi: ");
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		printf("%x ", phi[i]);
	}
	printf("\n");
#if 0
	printf("code_val: ");
	for(i = 0; i < (GF_FIELD - 1); i++)
	{
		printf("%x ", poly_eva(phi, (CODEWORD_LEN), power_polynomial_table[i + 1][0]));
	}
	printf("\n");
#endif	
	return 0;
}

unsigned char coordinate_trans(unsigned char locator, unsigned char r, unsigned char r_hd)
{
	unsigned int i = 0, j = 0;
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
	unsigned int i = 0, j = 0, k = 0, l = 0;
	unsigned int tmp = 0;
	unsigned char find_flag = 0;
	memset(syndrome, 0xFF, sizeof(unsigned char) * (CODEWORD_LEN - MESSAGE_LEN));

	syndrome_cal(received_polynomial, syndrome,
				  CODEWORD_LEN, MESSAGE_LEN);
#if 0//use re-encoding
#if 0
	printf("syn_val: ");
	for(i = 0; i < (GF_FIELD - 1); i++)
	{
		printf("%x ", poly_eva(syndrome, (CODEWORD_LEN - MESSAGE_LEN), power_polynomial_table[i + 1][0]));
	}
	printf("\n");
	printf("code_val: ");
	for(i = 0; i < (GF_FIELD - 1); i++)
	{
		printf("%x ", poly_eva(encoded_polynomial, (CODEWORD_LEN), power_polynomial_table[i + 1][0]));
	}
	printf("\n");
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
				//printf("%d %d %x\n", i, j, mul_matrix[i][j]);
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

	printf("beta:\n");
	for(j = 0; j < CODEWORD_LEN; j++)
	{
		for(i = 0; i < (CODEWORD_LEN + 1); i++)
		{
			printf("%x ", beta_matrix[i][j]);
		}
		printf("\n");
	}
	printf("\n");

	return 0;
}

int lex_order(unsigned int **lex_table, unsigned int d_x, unsigned int d_y)
{
	unsigned int i = 0, j = 0, k = 0;
	unsigned int max_degree = 0, degree_index = 0;
	unsigned int tmp_lex_table[d_x][d_y];

	printf("tmp_lex_table:\n");
	for(i = 0; i < d_x; i++)
	{
		for(j = 0; j < d_y; j++)
		{
#if 1	
			tmp_lex_table[i][j] = 1 * i + (MESSAGE_LEN - 1) * j;
#else//for test
			*((unsigned int *)lex_table + i * d_y + j) = 1 * i + 1 * j;
#endif
			if(max_degree < tmp_lex_table[i][j])
			{
				max_degree = tmp_lex_table[i][j];
			}
			printf("%d ", tmp_lex_table[i][j]);
		}
		printf("\n");
	}
	printf("\n");

	for(k = 0; k <= max_degree; k++)
	{
#if 1
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
					*((unsigned int *)lex_table + i * d_y + j) = degree_index;
					degree_index = degree_index + 1;
				}
			}
		}
	}
}

int koetter_interpolation()
{
	int i = 0, j = 0, k = 0, m = 0, n = 0;
	unsigned int a = 0, b = 0, v = 0;
	unsigned int tmp_sum = 0, tmp_real = 0;
	unsigned char tmp_ff = 0xFF;
	unsigned int l_s = 0xFFFF, l_w = 0, l_o = 0xFFFF;;

	unsigned int d_x = 0, d_y = 0, c = 0, term_num = 0, term_num_real = 0;
	unsigned int d_x_max = 0, d_y_max = 0;
	for(i = 0; i < (CODEWORD_LEN + 1); i++)
	{
		for(j = 0; j < CODEWORD_LEN; j++)
		{
			tmp_sum = tmp_sum + mul_matrix[i][j] * (mul_matrix[i][j] + 1);
		}
	}
	c = tmp_sum / 2;
	tmp_sum = pow((1 + 8 * (float)c / (MESSAGE_LEN - 1)), 0.5);
	tmp_sum = (1 + ((unsigned int)tmp_sum)) / 2;
	d_y = (unsigned int)(floor(tmp_sum));
	//d_y = floor((1 + (unsigned int)pow((1 + 8 * (float)c / (MESSAGE_LEN - 1)), 0.5)) / 2);
	d_x = floor(c / (d_y + 1) + d_y * (CODEWORD_LEN - MESSAGE_LEN - 1) / 2);

	d_x = d_x << 1;
	d_y = d_y << 1;

	unsigned int lex_order_table[d_x][d_y];
	lex_order((unsigned int **)lex_order_table, d_x, d_y);
	printf("lex_order_table:\n");
	for(i = 0; i < d_x; i++)
	{
		for(j = 0; j < d_y; j++)
		{
			printf("%d ", lex_order_table[i][j]);
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
		printf("\n");
	}
	//printf("constraint: %d %d %d %d\n", c, d_x, d_y, term_num);
	d_x_max = d_x;

	for(i = 0; i < (d_x + 1); i++)
	{
		//for(j = i; j < (d_y_max + 1); j++)
		for(j = 0; j < (d_y + 1); j++)
		{
			term_num_real = term_num_real + 1;
		}
	}
	
	printf("constraint: %d %d %d %d %d %d\n", c, d_x, d_y, term_num_real, d_x_max, d_y_max);
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
	unsigned char *g_table_x;
	unsigned char *g_table_y;
	unsigned char *discrepancy;
	unsigned char *weight_pol;
	unsigned char *lexorder_pol;
	unsigned char *tmp_table_c;
	unsigned char *term_use_index;//malloc later

	g_table_c = (unsigned char**)malloc(sizeof(unsigned char*) * (d_y_max + 1));
	g_table_c_prev = (unsigned char**)malloc(sizeof(unsigned char*) * (d_y_max + 1));
	for (i = 0; i < (d_y_max + 1); i++)
	{
  		g_table_c[i] = (unsigned char*)malloc(sizeof(unsigned char) * term_num_real);
		g_table_c_prev[i] = (unsigned char*)malloc(sizeof(unsigned char) * term_num_real);
  	}
	g_table_x = (unsigned char*)malloc(sizeof(unsigned char) * term_num_real);
	g_table_y = (unsigned char*)malloc(sizeof(unsigned char) * term_num_real);
	discrepancy = (unsigned char*)malloc(sizeof(unsigned char) * (d_y_max + 1));
	weight_pol = (unsigned char*)malloc(sizeof(unsigned char) * (d_y_max + 1));
	lexorder_pol = (unsigned char*)malloc(sizeof(unsigned char) * (d_y_max + 1));
	tmp_table_c = (unsigned char*)malloc(sizeof(unsigned char) * term_num_real);

	for(i = 0; i < (d_y_max + 1); i++)
	{
		for(j = 0; j < term_num_real; j++)
		{
			g_table_c[i][j] = 0xFF;
			g_table_c_prev[i][j] = 0xFF;
		}
	}
	memset(g_table_x, 0, sizeof(unsigned char) * term_num_real);
	memset(g_table_y, 0, sizeof(unsigned char) * term_num_real);
	memset(discrepancy, 0xFF, sizeof(unsigned char) * (d_y_max + 1));
	memset(tmp_table_c, 0xFF, sizeof(unsigned char) * term_num_real);
	memset(weight_pol, 0, sizeof(unsigned char) * (d_y_max + 1));
	memset(lexorder_pol, 0, sizeof(unsigned char) * (d_y_max + 1));

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
				printf("g_table: %d %d %d\n", k, g_table_x[k], g_table_y[k]);
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
				printf("g_table_c: %d %d %x\n", i, j, g_table_c[i][j]);
				g_table_c_prev[i][j] = g_table_c[i][j];
				break;
			}
		}
		//g_term[i] = j;
		//printf("g_term: %d %d\n", i, g_term[i]);
	}
	for(i = 0; i < (d_y_max + 1); i++)
	{
		for(j = 0; j < term_num_real; j++)
		{
			if(0xFF != g_table_c[i][j])
			{
				//printf("weight_pol_searching: %d %d\n", weight_pol[i], lex_order_table[g_table_x[j]][g_table_y[j]]);
#if 0				
				if(weight_pol[i] < lex_order_table[g_table_x[j]][g_table_y[j]])
				{
					weight_pol[i] = lex_order_table[g_table_x[j]][g_table_y[j]];
				}
#else
				if(weight_pol[i] < (g_table_x[j] + (MESSAGE_LEN - 1) * g_table_y[j]))
				{
					weight_pol[i] = (g_table_x[j] + (MESSAGE_LEN - 1) * g_table_y[j]);
				}
				if(lexorder_pol[i] < lex_order_table[g_table_x[j]][g_table_y[j]])
				{
					lexorder_pol[i] = lex_order_table[g_table_x[j]][g_table_y[j]];
				}
#endif
			}
		}
		//weight_pol[i] = (MESSAGE_LEN - 1) * i;
		printf("weight_pol: %d\n", weight_pol[i]);
	}
	//printf("g_table_c[0][0]: %x\n", g_table_c[0][0]);

	for(j = 0; j < CODEWORD_LEN; j++) //for I, as 2^q - 1
	{
		for(i = 0; i < (CODEWORD_LEN + 1); i++)//for I, as n
		{
			if(0 == mul_matrix[i][j])
			{
				continue;
			}
			printf("-------------------------------------------\n");
			printf("point: (%x %x)\n", power_polynomial_table[j + 1][0], beta_matrix[i][j]);

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
			term_use_index = (unsigned char*)malloc(sizeof(unsigned char) * term_num);
			m = 0;
			for(k = 0; k < term_num_real; k++)
			{
				if((mul_matrix[i][j] - 1) >= (g_table_x[k] + g_table_y[k]))
				{
					term_use_index[m] = k;
					printf("term_use_index: %d | %d | %d %d\n", m, term_use_index[m], g_table_x[k], g_table_y[k]);
					m = m + 1;
				}
			}

			/*notice that only (a, b) pairs is constained by (a + b) < m directly.*/
			/*(r, s) pairs is related to the point (a, b), but not the constrain.*/
			for(k = 0; k < term_num; k++) //for II, as (a, b) pairs
			{
				//printf("*********************\n");
				//printf("k ready: %d %d | %d | (%d %d)\n", term_num, k, term_use_index[k], g_table_x[term_use_index[k]], g_table_y[term_use_index[k]]);
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
							printf("m, n: %d %d | %d %d | %d %d | %x %x\n",
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
						//printf("++++++++++\n");
						//printf("(m, n) ready: %d (%d %d) %x\n", m, g_table_x[n], g_table_y[n], g_table_c_prev[m][n]);
						tmp_real = real_combine(g_table_x[n], g_table_x[term_use_index[k]]) * real_combine(g_table_y[n], g_table_y[term_use_index[k]]);
						//printf("tmp_real: %d | %d %d %d | %d %d %d\n", tmp_real, g_table_x[n], g_table_x[term_use_index[k]], real_combine(g_table_x[n], g_table_x[term_use_index[k]]), g_table_y[n], g_table_y[term_use_index[k]], real_combine(g_table_y[n], g_table_y[term_use_index[k]]));
						tmp_ff = gf_real_mutp_ff(tmp_real, g_table_c_prev[m][n]);
						//printf("tmp_ff: %x = %d mul %x\n", tmp_ff, tmp_real, g_table_c_prev[m][n]);
						//printf("tmp_ff: %x * (%x = %x^%x)\n", tmp_ff, gf_pow_cal(power_polynomial_table[j + 1][0], (g_table_x[n] - g_table_x[term_use_index[k]])), power_polynomial_table[j + 1][0], (g_table_x[n] - g_table_x[term_use_index[k]]));
						tmp_ff = gf_multp(tmp_ff, gf_pow_cal(power_polynomial_table[j + 1][0], (g_table_x[n] - g_table_x[term_use_index[k]])));
						//printf("tmp_ff: %x * (%x = %x^%x)\n", tmp_ff, gf_pow_cal(beta_matrix[i][j], (g_table_y[n] - g_table_y[term_use_index[k]])), beta_matrix[i][j], (g_table_y[n] - g_table_y[term_use_index[k]]));
						tmp_ff = gf_multp(tmp_ff, gf_pow_cal(beta_matrix[i][j], (g_table_y[n] - g_table_y[term_use_index[k]])));
						tmp_sum = gf_add(tmp_sum, tmp_ff);
						if(0xFF != tmp_sum)
						{
							//printf("tmp_sum: %x | %x\n", tmp_sum, tmp_ff);
						}
					}
					
					if(0xFF != tmp_sum)
					{
						discrepancy[m] = tmp_sum;
						printf("d: %d %d | %d %d | %d | %x\n", i, j, g_table_x[term_use_index[k]], g_table_y[term_use_index[k]], m, discrepancy[m]);
					}

					if(0xFF != discrepancy[m])
					{
						printf("updating place center: %d | %d vs %d | %d vs %d\n", l_s, l_w, weight_pol[m], l_o, lexorder_pol[m]);
						if(((l_w >= weight_pol[m]) && (l_o > lexorder_pol[m]))
							|| (0xFF == discrepancy[l_s]))
						{
							l_s = m;
							l_w = weight_pol[m];
							l_o = lexorder_pol[m];
							printf("updated place center: %d | %d | %d\n", l_s, l_w, l_o);
						}
					}
#if 0					
					if((0xFF != discrepancy[m])//find l_s
						&& (l_w >= weight_pol[m]))
					{
						l_s = m;
						l_w = weight_pol[m];
						printf("updated place center: %d %d\n", l_s, l_w);
					}
#endif
					//printf("g_table_c[0][0]: %x\n", g_table_c[0][0]);
				}

				printf("l_s: %d\n", l_s);

				//printf("g_table_c[0][0]: %x\n", g_table_c[0][0]);
				//printf("update polynomial ready: %d %d\n", l_s, l_w);
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
							printf("delta_l_s * q_l: %d %x %x %x\n", n,
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
							printf("delta_l * q_l_s: %d %x %x %x\n", n,
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
							printf("q_l add: %d %x %x %x\n", n,
														  gf_add(g_table_c[m][n], tmp_table_c[n]),
														  g_table_c[m][n],
														  tmp_table_c[n]);
						}
#endif
						g_table_c[m][n] = gf_add(g_table_c[m][n], tmp_table_c[n]);
					}

					for(n = 0; n < term_num_real; n++)
					{
						if(0xFF != g_table_c[m][n])
						{
							printf("Q_l: %d | %d %d | %x\n", m, g_table_x[n], g_table_y[n], g_table_c[m][n]);
						}
					}
				}
				//printf("update normal Q_l OK\n");

				//printf("g_table_c[0][0]: %x\n", g_table_c[0][0]);
				memset(tmp_table_c, 0xFF, sizeof(unsigned char) * term_num_real);
				for(n = 0; n < term_num_real; n++)//update Q_l_s
				{
					if(0xFF == discrepancy[l_s])
					{
						printf("no update for this (a, b) pair: %d\n", l_s);
						break;
					}
					
					if(0xFF != g_table_c_prev[l_s][n])
					{
						//printf("update Q_l_s: %d %d | %d %d | %x\n", l_s, n, g_table_x[n], g_table_y[n], g_table_c[l_s][n]);
						for(m = 0; m < term_num_real; m++)
						{
							//printf("checking m: %d | %d %d | %x\n", m, g_table_x[m], g_table_y[m], g_table_c[l_s][m]);
							if((g_table_x[m] == (g_table_x[n] + 1))
								&& (g_table_y[m] == g_table_y[n]))
							{
								tmp_table_c[m] = gf_add(g_table_c_prev[l_s][n], tmp_table_c[m]);
								//printf("update tmp_table_c[m]: %d %d | %d %d | %x\n", m, n, g_table_x[m], g_table_y[m], tmp_table_c[m]);
								break;
							}
						}

						//printf("updating g_table_c[l_s][term_use_index[n]]: %x %x\n", g_table_c[l_s][term_use_index[n]], power_polynomial_table[j + 1][0]);
						//g_table_c[l_s][n] = gf_multp(g_table_c[l_s][n], power_polynomial_table[j + 1][0]);
						tmp_table_c[n] = gf_add(gf_multp(g_table_c_prev[l_s][n], power_polynomial_table[j + 1][0]), tmp_table_c[n]);
						//printf("update tmp_table_c[n]: %d %d | %x\n", g_table_x[m], g_table_y[m], tmp_table_c[n]);
						//g_table_c[l_s][n] = gf_add(g_table_c[l_s][n], tmp_table_c[n]);
						//g_table_c[l_s][m] = gf_add(g_table_c[l_s][m], tmp_table_c[m]);
						//printf("update_2 g_table_c[l_s][index]: %d %d %x %x\n", n, m, g_table_c[l_s][n], g_table_c[l_s][m]);
						//break;
					}
				}
				for(n = 0; n < term_num_real; n++)
				{
			#if 0		
					if(0xFF != g_table_c[l_s][n])
					{
						printf("Q_l_s_c: %d | %d %d | %x\n", l_s, g_table_x[n], g_table_y[n], g_table_c[l_s][n]);
					}
			#endif
					if(0xFF == discrepancy[l_s])
					{
						printf("no update for this (a, b) pair: %d\n", l_s);
						break;
					}

					g_table_c[l_s][n] = 0xFF;

					if(0xFF != tmp_table_c[n])
					{
						g_table_c[l_s][n] = gf_multp(tmp_table_c[n], discrepancy[l_s]);
						printf("Q_l_s_c: %d | %d %d | %x\n", l_s, g_table_x[n], g_table_y[n], g_table_c[l_s][n]);
					}
				}
				//printf("update Q_l_s OK\n");

				for(m = 0; m < (d_y_max + 1); m++)
				{
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
							if(weight_pol[m] < (g_table_x[n] + (MESSAGE_LEN - 1) * g_table_y[n]))
							{
								weight_pol[m] = (g_table_x[n] + (MESSAGE_LEN - 1) * g_table_y[n]);
							}
							if(lexorder_pol[m] < lex_order_table[g_table_x[n]][g_table_y[n]])
							{
								lexorder_pol[m] = lex_order_table[g_table_x[n]][g_table_y[n]];
							}
#endif
						}
					}
					//weight_pol[i] = (MESSAGE_LEN - 1) * i;
					printf("weight_pol: %d\n", weight_pol[m]);
				}
				//weight_pol[l_s] = weight_pol[l_s] + 1;//update w_l_s
				l_w = weight_pol[l_s];
				l_o = lexorder_pol[l_s];
				//printf("update w_l_s OK\n");

				/*store g as prev_poly*/
				for(m = 0; m < (d_y_max + 1); m++)
				{
					for(n = 0; n < term_num_real; n++)
					{
						g_table_c_prev[m][n] = g_table_c[m][n];
						if(0xFF != g_table_c_prev[m][n])
						{
							printf("g_table_c_prev: %d | %d %d | %d %d | %x\n", m, lex_order_table[g_table_x[n]][g_table_y[n]], (g_table_x[n] + (MESSAGE_LEN - 1) * g_table_y[n]), g_table_x[n], g_table_y[n], g_table_c_prev[m][n]);
						}
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
			printf("sml_updated: %d | %d %d | %d\n", m, tmp_real, weight_pol[m], lexorder_pol[m]);
			sml_poly = m;
			tmp_real = weight_pol[m];
		}
		if((tmp_real == weight_pol[m])
			&& (lexorder_pol[sml_poly] > lexorder_pol[m]))
		{
			printf("sml_updated: %d | %d %d | %d\n", m, tmp_real, weight_pol[m], lexorder_pol[m]);
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
						printf("g_term: %d | %d %d | %x\n", m,
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
					g_term_c[0][0][k] = g_table_c[sml_poly][n];
					printf("g_term: %d | %d %d | %x\n", sml_poly,
														g_term_x[k],
														g_term_y[k],
														g_term_c[0][0][k]);
					break;
				}
			}
		}
	}
#endif

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

int g_term_init()
{
	unsigned int i = 0, j = 0, k = 0;

	for(k = 0; k < LAYER_NUM; k++)
	{
		for(i = 0; i < POLY_NUM; i++)
		{
			for(j = 0; j < (TERM_SIZE * TERM_SIZE); j++)
			{
				g_term_c[k][i][j] = 0xFF;
				g_term_0_y_c[k][i][j] = 0xFF;
				g_term_x_0_c[k][i][j] = 0xFF;
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
	unsigned int i = 0, j = 0;

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

int g_term_new_gen(unsigned char layer_idx, unsigned char tern_idx, unsigned char root_insert)
{
	unsigned int i = 0, j = 0, k = 0, l = 0, r = 0, s = 0;

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
			printf("g_term_init: %d %d | %x %x\n",
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
			printf("g_term_init: %d %d | %x %x\n",
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
		memcpy(g_term_c[layer_idx][tern_idx],
			   g_term_c[layer_idx][tern_idx - 1],
			   sizeof(unsigned char) * k);
	}

	for(k = 0; k < (TERM_SIZE * TERM_SIZE); k++)//for every term contain "y"
	{
		g_term_c[layer_idx + 1][tern_idx][k] = g_term_c[layer_idx][tern_idx][k];
		
		if((0 != g_term_y[k])
			&& (0xFF != g_term_c[layer_idx][tern_idx][k]))
		{
			printf("*********************\n");
			printf("g_term_have_y: %d | %d %d | %x\n",
				    tern_idx,
				    g_term_x[k],
				    g_term_y[k],
				    g_term_c[layer_idx][tern_idx][k]);
			
			for(l = 0; l < (g_term_y[k] - 0); l++)//for pow_cal, g*m => tmp
			{
				//printf("*********************\n");
				for(i = 0; i < (TERM_SIZE * TERM_SIZE); i++)//for every a
				{
					if(0xFF != g_term_c_expand[i])
					{
						printf("g_term_c_expand: %d %d | %x\n",
					    g_term_x[i],
					    g_term_y[i],
					    g_term_c_expand[i]);

						for(j = 0; j < (TERM_SIZE * TERM_SIZE); j++)//for every b
						{
							if(0xFF != mul_g_term_c_expand[j])
							{
								printf("mul_g_term_c_expand: %d %d | %x\n",
							    g_term_x[j],
							    g_term_y[j],
							    mul_g_term_c_expand[j]);

								for(r = 0; r < (TERM_SIZE * TERM_SIZE); r++)//find the place
								{
									if(((g_term_x[i] + g_term_x[j]) == g_term_x[r])
										&& ((g_term_y[i] + g_term_y[j]) == g_term_y[r]))//r <= place(i, j)
									{
										tmp_g_term_c_expand[r] = gf_add(tmp_g_term_c_expand[r], gf_multp(g_term_c_expand[i], mul_g_term_c_expand[j]));
										printf("tmp_g_term_c_expand: %d %d | %x\n",
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

					if(0xFF != g_term_c_expand[r])
					{
						printf("g_term_expand_update: %d %d | %x\n",
							    g_term_x[r],
							    g_term_y[r],
							    g_term_c_expand[r]);
					}
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
							printf("g_term_expand_store updating: %d + %d = %d, %d | %x %x %x\n",
								    g_term_x[r],
								    g_term_x[k],
							    	g_term_x[s],
							    	g_term_y[s],
							    	g_term_c_expand[r],
							    	g_term_c[layer_idx][tern_idx][k],
							    	g_term_c_expand_store[s]);
							//g_term_c_expand_store[s] = gf_add(g_term_c_expand_store[s], g_term_c_expand[r]);
							//g_term_c_expand_store[s] = gf_multp(g_term_c_expand_store[s], g_term_c[layer_idx][tern_idx][k]);
							tmp = gf_multp(g_term_c_expand[r], g_term_c[layer_idx][tern_idx][k]);
							g_term_c_expand_store[s] = gf_add(g_term_c_expand_store[s], tmp);
							printf("g_term_expand_store: %d %d | %x\n",
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
			g_term_c[layer_idx + 1][tern_idx][k] = 0xFF;

		}
	}

	printf("*********************\n");
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
		
		g_term_c[layer_idx + 1][tern_idx][r] = gf_add(g_term_c_expand_store[r], g_term_c[layer_idx + 1][tern_idx][r]);
		if(0xFF != g_term_c[layer_idx + 1][tern_idx][r])
		{
			printf("g_term_update: %d | %d %d | %x\n",
				    tern_idx,
				    g_term_x[r],
				    g_term_y[r],
				    g_term_c[layer_idx + 1][tern_idx][r]);
		}
	}

	printf("*********************\n");
	/*eliminate common factor*/
	unsigned char cmn_factor = 0xFF;//deal with x
	for(r = 0; r < (TERM_SIZE * TERM_SIZE); r++)
	{
		if(0xFF != g_term_c[layer_idx + 1][tern_idx][r])
		{
			if(cmn_factor > g_term_x[r])
			{
				cmn_factor = g_term_x[r];
			}
		}
	}
	printf("cmn_factor for x: %d\n", cmn_factor);
	if(0 != cmn_factor)
	{
		for(r = 0; r < (TERM_SIZE * TERM_SIZE); r++)
		{
			if(0xFF != g_term_c[layer_idx + 1][tern_idx][r])
			{
				for(s = 0; s < (TERM_SIZE * TERM_SIZE); s++)
				{
					if((g_term_x[s] == (g_term_x[r] - cmn_factor))
						&& (g_term_y[s] == g_term_y[r]))
					{
						g_term_c[layer_idx + 1][tern_idx][s] = g_term_c[layer_idx + 1][tern_idx][r];
						g_term_c[layer_idx + 1][tern_idx][r] = 0xFF;
						printf("g_term_update_x: %d | %d %d | %x\n",
							    tern_idx,
							    g_term_x[s],
							    g_term_y[s],
							    g_term_c[layer_idx + 1][tern_idx][s]);
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
	printf("cmn_factor for y: %d\n", cmn_factor);
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
						printf("g_term_update_y: %d | %d %d | %x\n",
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

	printf("*********************\n");
	for(r = 0; r < (TERM_SIZE * TERM_SIZE); r++)
	{
		if(0xFF != g_term_c[layer_idx + 1][tern_idx][r])
		{
			printf("g_term_update_OK: %d | %d %d | %x\n",
				    tern_idx,
				    g_term_x[r],
				    g_term_y[r],
				    g_term_c[layer_idx + 1][tern_idx][r]);
		}
	}
	
	return 0;
}

int g_term_0_y_cal(unsigned char layer_idx, unsigned char tern_idx)
{
	unsigned int i = 0, j = 0;

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
				printf("g_term_set_to_zero: %d | %d %d | %x %x\n",
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
					printf("g_term_0_y: %d | %d %d | %x\n",
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
		g_term_0_y_c[layer_idx][tern_idx][j] = 0xFF;
	}

	for(j = 0; j < (TERM_SIZE * TERM_SIZE); j++)
	{
		if((0 != g_term_x[j])
			&& (0xFF != g_term_c[layer_idx][tern_idx][j]))
		{
			g_term_0_y_c[layer_idx][tern_idx][j] = 0xFF;
			printf("g_term_set_to_zero: %d | %d %d | %x %x\n",
				    i,
				    g_term_x[j],
				    g_term_y[j],
				    g_term_c[layer_idx][tern_idx][j],
				    g_term_0_y_c[layer_idx][tern_idx][j]);
		}
		else
		{
			g_term_0_y_c[layer_idx][tern_idx][j] = g_term_c[layer_idx][tern_idx][j];
			if(0xFF != g_term_0_y_c[layer_idx][tern_idx][j])
			{
				printf("g_term_0_y: %d | %d %d | %x\n",
						i,
						g_term_x[j],
						g_term_y[j],
						g_term_0_y_c[layer_idx][tern_idx][j]);
			}
		}
	}
#endif

	return 0;
}

unsigned char g_term_x_0_cal(unsigned char layer_idx, unsigned char tern_idx)
{
	unsigned char val = 0xFF;

	unsigned int i = 0, j = 0;

	/*clear*/
	for(j = 0; j < (TERM_SIZE * TERM_SIZE); j++)
	{
		g_term_x_0_c[layer_idx][tern_idx][j] = 0xFF;
	}

	for(j = 0; j < (TERM_SIZE * TERM_SIZE); j++)
	{
		if((0 != g_term_y[j])
			&& (0xFF != g_term_c[layer_idx][tern_idx][j]))
		{
			g_term_x_0_c[layer_idx][tern_idx][j] = 0xFF;
			printf("g_term_set_to_zero: %d | %d %d | %x %x\n",
				    i,
				    g_term_x[j],
				    g_term_y[j],
				    g_term_c[layer_idx][tern_idx][j],
				    g_term_x_0_c[layer_idx][tern_idx][j]);
		}
		else
		{
			g_term_x_0_c[layer_idx][tern_idx][j] = g_term_c[layer_idx][tern_idx][j];
			if(0xFF != g_term_x_0_c[layer_idx][tern_idx][j])
			{
				printf("g_term_0_y: %d | %d %d | %x\n",
						i,
						g_term_x[j],
						g_term_y[j],
						g_term_x_0_c[layer_idx][tern_idx][j]);
			}
		}

		val = gf_add(val, g_term_x_0_c[layer_idx][tern_idx][j]);
	}
	
	return val;
}

int chien_searching_for_g_0_y(unsigned char layer_idx, unsigned char tern_idx, unsigned char root_test)
{
	int is_root = 0;
	unsigned int k = 0;
	unsigned char tmp = 0xFF, tmp_sum = 0xFF;
	
	for(k = 0; k < (TERM_SIZE * TERM_SIZE); k++)
	{
		if(0xFF != g_term_0_y_c[layer_idx][tern_idx][k])
		{
			tmp = 0xFF;
			tmp = gf_pow_cal(root_test, g_term_y[k]);
			//printf("tmp: %x | %x^%d\n", tmp, root_test, g_term_y[k]);
			//printf("tmp: %x | %x*%x\n", gf_multp(tmp, g_term_0_y_c[tern_idx][k]), tmp, g_term_0_y_c[tern_idx][k]);
			tmp = gf_multp(tmp, g_term_0_y_c[layer_idx][tern_idx][k]);
			//printf("tmp_sum: %x | %x+%x\n", gf_add(tmp_sum, tmp), tmp_sum, tmp);
			tmp_sum = gf_add(tmp_sum, tmp);
		}
	}

	if(0xFF == tmp_sum)
	{
		is_root = 1;
		printf("is_root: %d %x\n", tern_idx, root_test);
	}
	else
	{
		is_root = 0;
	}

	return is_root;
}

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

int rr_factorization()
{
	int i = 0, j = 0, k = 0, l = 0, r = 0, s = 0;
	int is_root = 0;

	//test_factorization_init();
	f_root_init();

	//chien_searching_for_g_0_y(1, power_polynomial_table[0][0]);//test
	//g_term_new_gen(1, power_polynomial_table[0][0]);//test

	for(s = 0; s < MESSAGE_LEN; s++)//layer
	{
		printf("-------------------------------------------\n");
		l = 0;
		for(r = 0; r < f_root_cnt[s]; r++)//roots per layer
		{
			g_term_0_y_cal(s, r);

			for(k = 0; k < GF_FIELD; k++)//test roots
			{
				is_root = chien_searching_for_g_0_y(s, r, power_polynomial_table[k][0]);
				if(1 == is_root)
				{
					f_root_val[s][l] = power_polynomial_table[k][0];
					f_root_prev[s][l] = r;
					f_root_cnt[s + 1] = f_root_cnt[s + 1] + 1;

					printf("f_root: %d %d | %x\n", s, l, f_root_val[s][l]);

					g_term_new_gen(s, l, f_root_val[s][l]);

					l = l + 1;
				}
			}
		}
	}

	for(r = 0; r < f_root_cnt[MESSAGE_LEN]; r++)
	{
		if(0xFF == g_term_x_0_cal(MESSAGE_LEN, r))
		{
			printf("Decoding OK!\n");

			printf("Message: %d %d | %x\n", MESSAGE_LEN - 1, r, f_root_val[MESSAGE_LEN - 1][r]);
			for(s = MESSAGE_LEN - 2; s >= 0; s--)
			{
				printf("Message: %d %d | %x\n", s, r, f_root_val[s][f_root_prev[s + 1][r]]);
			}
		}
		else
		{
			printf("Decoding Fail!\n");
		}
	}

	return 0;
}

int as_decoding()
{
	unsigned int i = 0;

	printf("Received:\n");
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		printf("%x ", received_polynomial[i]);
	}
	printf("\n");

	g_term_init();

	koetter_interpolation();

	rr_factorization();

	return 0;
}
