#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "gf_cal.h"
#include "encoding.h"
#include "as_decoding.h"

#define S_MUL	10
#define K_M		8

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
	unsigned int i = 0, j = 0;

	for(i = 0; i < d_x; i++)
	{
		for(j = 0; j < d_y; j++)
		{
			*((unsigned int *)lex_table + i * d_y + j) = 1 * i + (MESSAGE_LEN - 1) * j;
		}
	}
}

int koetter_interpolation()
{
	int i = 0, j = 0, k = 0, m = 0, n = 0;
	unsigned int a = 0, b = 0, v = 0;
	unsigned int tmp_sum = 0;

	unsigned int d_x = 0, d_y = 0, c = 0, lm = 0;
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

	unsigned int lex_order_table[d_x][d_y];
	lex_order((unsigned int **)lex_order_table, d_x, d_y);
	printf("lex_order_table:\n");
	for(i = 0; i < d_x; i++)
	{
		for(j = 0; j < d_y; j++)
		{
			printf("%d ", lex_order_table[i][j]);
			if(K_M >= lex_order_table[i][j])
			{
				lm = lm + 1;
				if(d_x_max < i)
				{
					d_x_max = i;
				}
				if(d_y_max < j)
				{
					d_y_max = j;
				}
			}
		}
		printf("\n");
	}
	printf("constraint: %x %x %x %x %x %x\n", c, d_x, d_y, lm, d_x_max, d_y_max);
	unsigned char g_table_c[lm];
	unsigned char g_table_x[lm];
	unsigned char g_table_y[lm];
	memset(g_table_c, 0, sizeof(unsigned char) * lm);
	memset(g_table_x, 0, sizeof(unsigned char) * lm);
	memset(g_table_y, 0, sizeof(unsigned char) * lm);
	k = 0;
	for(tmp_sum = 0; tmp_sum <= K_M; tmp_sum++)
	{
		for(i = d_x_max; i >= 0; i--)
		{
			for(j = d_y_max; j >= 0; j--)
			{
				if(tmp_sum == lex_order_table[i][j])
				{
					g_table_c[k] = 0;
					g_table_x[k] = i;
					g_table_y[k] = j;
					printf("g_table: %x %x %x %x\n", tmp_sum, k, g_table_x[k], g_table_y[k]);
					k = k + 1;
				}
			}
		}
	}
	
	unsigned char g_term[d_y_max + 1];//j, degree of x, degree of y
	for(i = 0; i <= d_y_max; i++)
	{
		for(j = 0; j < lm; j++)
		{
			if((0 == g_table_x[j]) && (i == g_table_y[j]))
			{
				break;
			}
		}
		g_term[i] = j;
		printf("g_term: %d %d\n", i, g_term[i]);
	}
	
	unsigned char discrepancy[d_y_max];
	memset(discrepancy, 0xFF, sizeof(unsigned char) * (d_y_max + 1));
	for(i = 0; i < (CODEWORD_LEN + 1); i++)
	{
		for(j = 0; j < CODEWORD_LEN; j++)
		{
			if(0 == mul_matrix[i][j])
			{
				continue;
			}

			for(a = 0; a < (mul_matrix[i][j] - 1); a++)
			{
				for(b = 0; b < (mul_matrix[i][j] - a - 1); b++)
				{
					memset(discrepancy, 0xFF, sizeof(unsigned char) * (d_y_max + 1));
					for(v = 0; v < (d_y_max + 1); v++)
					{
#if 0//need to be checked
						tmp_sum = 0;
						for(n = b; n < d_y_max; n++)
						{
							for(m = a; m < (d_x + d_y - n); m++)
							{
								tmp_sum = tmp_sum
										+ real_combine(m, a) * real_combine(n, b)
										* (power_polynomial_table[j + 1][0] * (m - a) % (GF_FIELD - 1))
										* (beta_matrix[i][j] * (n - b) % (GF_FIELD - 1))
										* g_table_c[g_term[v]];
							}
						}
#endif
						discrepancy[v] = tmp_sum;
					}
				}
			}
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

	koetter_interpolation();

	return 0;
}
