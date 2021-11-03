#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "gf_cal.h"
#include "as_decoding.h"

#define S_MUL	31

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
unsigned char rel_group_seq[MESSAGE_LEN];
unsigned char unrel_group_seq[CODEWORD_LEN - MESSAGE_LEN];
unsigned char syndrome[CODEWORD_LEN - MESSAGE_LEN];
unsigned char tao[CODEWORD_LEN - MESSAGE_LEN + 1];
unsigned char sigma[((CODEWORD_LEN - MESSAGE_LEN) + (CODEWORD_LEN - MESSAGE_LEN + 1) - 1) - (CODEWORD_LEN - MESSAGE_LEN)];

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
	reg[0] = unrel_group_seq[0];
	reg[1] = 0;
	unsigned char b[CODEWORD_LEN - MESSAGE_LEN];
	memset(b, 0xFF, sizeof(unsigned char) * (CODEWORD_LEN - MESSAGE_LEN));

	for(i = 1; i < (CODEWORD_LEN - MESSAGE_LEN); i++)
	{
		a[0] = unrel_group_seq[i];
		a[1] = 0;
		memcpy(b, reg, sizeof(unsigned char) * (1 + i));
		gf_multp_poly_hw(a, 2,
						  b, (1 + i),
						  reg, (2 + i));
#if 0
		printf("reg: ");
		for(j = 0; j < (2 + i); j++)
		{
			printf("%x ", reg[j]);
		}
		printf("\n");
#endif
	}

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

	return 0;
}

unsigned char poly_eva(unsigned char *poly, unsigned char poly_len, unsigned char input_val)
{
	unsigned int i = 0;
	unsigned char poly_val = 0xFF, tmp_product = 0;

	for(i = 0; i < poly_len; i++)
	{
		tmp_product = gf_multp(*(poly + i), (input_val * i) % (GF_FIELD - 1));
		poly_val = gf_add(tmp_product, poly_val);
	}
	
	return poly_val;
}

unsigned char coordinate_trans(unsigned char locator, unsigned char r, unsigned char r_hd)
{
	unsigned int i = 0, j = 0;
	unsigned char coordinate = 0xFF;

	unsigned char locator_product = 0, tmp_sum = 0xFF;
	for(i = 0; i < CODEWORD_LEN - MESSAGE_LEN; i++)
	{
		if(locator == unrel_group_seq[i])
		{
			continue;
		}
		tmp_sum = gf_add(locator, unrel_group_seq[i]);
		locator_product = gf_multp(locator_product, tmp_sum);
	}
	locator_product = gf_multp(locator_product, locator);
	tmp_sum = gf_add(r_hd, r);
	locator_product = gf_multp(tmp_sum, locator_product);

	tmp_sum = poly_eva(sigma, ((CODEWORD_LEN - MESSAGE_LEN) + (CODEWORD_LEN - MESSAGE_LEN + 1) - 1) - (CODEWORD_LEN - MESSAGE_LEN), locator);

	coordinate = gf_add(locator_product, tmp_sum);
	
	return coordinate;
}

int re_encoding()
{
	unsigned int i = 0, j = 0, k = 0, l = 0;
	unsigned int tmp = 0;
	memset(syndrome, 0xFF, sizeof(unsigned char) * (CODEWORD_LEN - MESSAGE_LEN));

	syndrome_cal(received_polynomial, syndrome,
				  CODEWORD_LEN, MESSAGE_LEN);

	rel_group();

	tao_cal();
	sigma_cal();

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

	return 0;
}
