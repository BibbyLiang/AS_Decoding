#include <stdio.h>
#include <string.h>
#include "gf_cal.h"

unsigned char power_polynomial_table[GF_FIELD][2] = 
{
	/*power <---> polynomial*/
	/*These are also the coefficients of generator polynomials.*/
	{0xFF, 0x0},
	{0x0, 0x1},
	{0x1, 0x2},
	{0x2, 0x4},
	{0x3, 0x3},
	{0x4, 0x6},
	{0x5, 0x7},
	{0x6, 0x5}
};

unsigned char gf_pow2poly(unsigned char val_in_pow)
{
	unsigned char val_in_poly = 0;

	return power_polynomial_table[val_in_pow + 1][1];
}

unsigned char gf_poly2pow(unsigned char val_in_poly)
{
	unsigned char i = 0;
	unsigned char val_in_pow = 0;

	for(i = 0; i < GF_FIELD; i++)
	{
		if(power_polynomial_table[i][1] == val_in_poly)
		{
			val_in_pow = power_polynomial_table[i][0];
			break;
		}
	}

	return val_in_pow;
}

unsigned char gf_location(unsigned char val)
{
	unsigned char val_location = power_polynomial_table[val + 1][0];

	return val_location;
}

unsigned char gf_add(unsigned char a, unsigned char b)
{
	unsigned char i = 0;
	unsigned char sum_in_pow = 0;
	
	unsigned char sum_in_poly = gf_pow2poly(a) ^ gf_pow2poly(b);
	sum_in_pow = gf_poly2pow(sum_in_poly);

	return sum_in_pow;
}

unsigned char gf_multp(unsigned char a, unsigned char b)
{
	if((0xFF == a) || (0xFF == b))
	{
		return 0xFF;
	}

	unsigned char product_in_pow = (a + b) % (GF_FIELD - 1);

	return product_in_pow;
}

unsigned char gf_div(unsigned char a, unsigned char b)
{
	if(0xFF == a)
	{
		return 0xFF;
	}
	if(0xFF == b)
	{
		printf("div err.\n");
		return 0xFF;
	}

	//printf("div: %x %x\n", a, b);
	unsigned char quotient_in_pow = 0;
	if(a >= b)
	{
		quotient_in_pow = (a - b) % (GF_FIELD - 1);
	}
	else
	{
		quotient_in_pow = ((b / (GF_FIELD - 1) + 1) * (GF_FIELD - 1) + a - b) % (GF_FIELD - 1);
	}

	return quotient_in_pow;
}

unsigned char gf_mod_single_term(unsigned char a, unsigned char b)
{
	unsigned char remainder_in_pow = 0xFF;
	if(a >= b)
	{
		/*remove small term*/
		remainder_in_pow = 0xFF;
	}
	else
	{
		/*keep big term*/
		remainder_in_pow = 0;
	}

	return remainder_in_pow;
}

unsigned char gf_degree(unsigned char* a, unsigned char len_a)
{
	unsigned char i = 0;

	for(i = len_a - 1; i >= 0; i--)
	{
		if(0xFF != a[i])
		{
			break;
		}
	}

	return i;
}

unsigned char gf_div_q_r(unsigned char* dividend, unsigned char len_dividend,
							unsigned char* divisor, unsigned char len_divisor,
							unsigned char* quotien, unsigned char len_quotien,
							unsigned char* remainder, unsigned char len_remainder)
{
	unsigned char i = 0, j = 0, k = 0;
	unsigned char locator = 0, factor = 0, locator_rmd = 0, factor_rmd = 0;
	unsigned char dividend_tmp[len_dividend], remainder_tmp[len_remainder];
	memset(dividend_tmp, 0xFF, sizeof(unsigned char) * len_dividend);
	memset(remainder_tmp, 0xFF, sizeof(unsigned char) * len_remainder);

	if(gf_degree(divisor, len_divisor) > gf_degree(dividend, len_dividend))
	{
		for(i = 0; i < len_remainder; i++)
		{
			remainder[i] = dividend[i];
		}
		//printf("quotien is zero: %d %d\n", gf_degree(divisor, len_divisor), gf_degree(dividend, len_dividend));

		return 0;
	}

	memcpy(dividend_tmp, dividend, sizeof(unsigned char) * len_dividend);
	for(i = 0; i < len_dividend; i++)
	{
		locator = gf_degree(dividend_tmp, len_dividend) - gf_degree(divisor, len_divisor);
		factor = gf_div(dividend_tmp[gf_degree(dividend_tmp, len_dividend)], divisor[gf_degree(divisor, len_divisor)]);

		quotien[locator] = factor;
		//printf("quotien: %x %d %x\n", quotien[locator], locator, factor);

		for(j = 0; j < len_divisor; j++)
		{
			factor_rmd = gf_multp(factor, divisor[gf_degree(divisor, len_divisor) - j]);
			locator_rmd = locator + gf_degree(divisor, len_divisor) - j;
			remainder_tmp[locator_rmd] = gf_add(dividend_tmp[locator_rmd], factor_rmd);
		}
#if 0
		printf("remainder_tmp:\n");
		for(k = 0; k < len_remainder; k++)
		{
			printf("%x ", remainder_tmp[k]);
		}
		printf("\n");
#endif
		if(gf_degree(divisor, len_divisor) > gf_degree(remainder_tmp, len_remainder))
		{
			for(k = 0; k < len_remainder; k++)
			{
				remainder[k] = remainder_tmp[k];
			}
			break;
		}
		else
		{
			memcpy(dividend_tmp, remainder_tmp, sizeof(unsigned char) * len_dividend);
			memset(remainder_tmp, 0xFF, sizeof(unsigned char) * len_remainder);
		}
	}

	return 0;
}

unsigned char gf_multp_poly(unsigned char* a, unsigned char len_a,
								unsigned char* b, unsigned char len_b,
								unsigned char* product, unsigned char len_product)
{
	unsigned char i = 0, j = 0, idx = 0;

	for(i = 0; i < len_a; i++)
	{
		for(j = 0; j < len_b; j++)
		{
			idx = i + j;
			if(len_product <= idx)
			{
				//printf("product len err: %d\n", idx);
				continue;
			}
			product[idx] = gf_add(product[idx], gf_multp(a[i], b[j]));
		}
	}
}

int gf_multp_poly_hw(unsigned char* a, unsigned char len_a,
				 		   unsigned char* b, unsigned char len_b,
				 		   unsigned char* product, unsigned char len_product)
{
	unsigned int i = 0, j = 0, idx = len_product - 1;
	unsigned char reg[len_a - 1];
	memset(reg, 0xFF, sizeof(unsigned char) * (len_a - 1));
	unsigned char pd_tmp = 0xFF;

	/*high -> low*/
	for(i = 0; i < len_b; i++)
	{
		pd_tmp = gf_multp(*(a + len_a - 1), *(b + (len_b - 1 - i)));
		product[idx] = gf_add(reg[len_a - 2], pd_tmp);
		//printf("unrel_group_seq: %x\n", a[0]);
#if 0
		printf("%x %x %x %x %x %x\n", 
									idx,
									product[idx],
					  	   	 		reg[len_a - 2],
					       	 		pd_tmp,
					       	 		*(a + len_a - 1),
					       	 		*(b + (len_b - 1 - i)));
#endif
		if(0 >= idx)
		{
			printf("product len err: %d\n", idx);
			break;
		}
		idx = idx - 1;

		for(j = 1; j < (len_a - 1); j++)
		{
			pd_tmp = gf_multp(*(a + len_a - 1 - j), *(b + (len_b - 1 - i)));
			reg[len_a - 1 - j] = gf_add(reg[len_a - 1 - j - 1], pd_tmp);
#if 0
			printf("idx: %d\n", len_a - 1 - j);
			printf("%x %x %x %x %x\n", reg[len_a - 1 - j], reg[len_a - 1 - j - 1], pd_tmp,
								 	   *(a + len_a - 1 - j), *(b + (len_b - 1 - i)));
#endif
		}
		reg[0] = gf_multp(*(a + 0), *(b + (len_b - 1 - i)));
		//printf("%x %x %x\n", reg[0], *(a + 0), *(b + (len_b - 1 - i)));
	}

	for(i = 0; i < len_a; i++)
	{
		pd_tmp = gf_multp(*(a + len_a - 1), 0xFF);
		product[idx] = gf_add(reg[len_a - 2], pd_tmp);
#if 0
		printf("%x %x %x %x %x %x\n", 
									idx,
									product[idx],
					  	   	 		reg[len_a - 2],
					       	 		pd_tmp,
					       	 		*(a + len_a - 1),
					       	 		0xFF);
#endif
		if(0 >= idx)
		{
			printf("product len err: %d\n", idx);
			break;
		}
		idx = idx - 1;

		for(j = 1; j < (len_a - 1); j++)
		{
			pd_tmp = gf_multp(*(a + len_a - 1 - j - 1), 0xFF);
			reg[len_a - 1 - j] = gf_add(reg[len_a - 1 - j - 1], pd_tmp);
		}
		reg[0] = gf_multp(*(a + 0), 0xFF);
	}
	
	return 0;
}
