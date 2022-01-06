#include <stdio.h>
#include <string.h>
#include "debug_info.h"
#include "gf_cal.h"

#if (3 == GF_Q)
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
#endif

#if (4 == GF_Q)
unsigned char power_polynomial_table[GF_FIELD][2] = 
{
	/*power <---> polynomial*/
	/*These are also the coefficients of generator polynomials.*/
	{0xFF, 0x0},
	{0x0, 0x1},
	{0x1, 0x2},
	{0x2, 0x4},
	{0x3, 0x8},
	{0x4, 0x3},
	{0x5, 0x6},
	{0x6, 0xc},
	{0x7, 0xb},
	{0x8, 0x5},
	{0x9, 0xa},
	{0xa, 0x7},
	{0xb, 0xe},
	{0xc, 0xf},
	{0xd, 0xd},
	{0xe, 0x9}
};
#endif

#if (6 == GF_Q)
unsigned char power_polynomial_table[GF_FIELD][2] = 
{
	/*power <---> polynomial*/
	/*These are also the coefficients of generator polynomials.*/
	{0xFF, 0},
	{0, 1},
	{1, 2},
	{2, 4},
	{3, 8},
	{4, 16},
	{5, 32},
	{6, 3},
	{7, 6},
	{8, 12},
	{9, 24},
	{10, 48},
	{11, 35},
	{12, 5},
	{13, 10},
	{14, 20},
	{15, 40},
	{16, 19},
	{17, 38},
	{18, 15},
	{19, 30},
	{20, 60},
	{21, 59},
	{22, 53},
	{23, 41},
	{24, 17},
	{25, 34},
	{26, 7},
	{27, 14},
	{28, 28},
	{29, 56},
	{30, 51},
	{31, 37},
	{32, 9},
	{33, 18},
	{34, 36},
	{35, 11},
	{36, 22},
	{37, 44},
	{38, 27},
	{39, 54},
	{40, 47},
	{41, 29},
	{42, 58},
	{43, 55},
	{44, 45},
	{45, 25},
	{46, 50},
	{47, 39},
	{48, 13},
	{49, 26},
	{50, 52},
	{51, 43},
	{52, 21},
	{53, 42},
	{54, 23},
	{55, 46},
	{56, 31},
	{57, 62},
	{58, 63},
	{59, 61},
	{60, 57},
	{61, 49},
	{62, 33}
};
#endif

unsigned char gf_pow2poly(unsigned char val_in_pow)
{
	unsigned char val_in_poly = 0;
	if(0xFF == val_in_pow)
	{
		return power_polynomial_table[0][1];
	}
	else
	{
		return power_polynomial_table[val_in_pow + 1][1];
	}
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
	unsigned char val_location = 0xFF;

	if(0xFF == val)
	{
		val_location = power_polynomial_table[0][0];
	}
	else
	{
		val_location = power_polynomial_table[val + 1][0];
	}

	return val_location;
}

unsigned char gf_add(unsigned char a, unsigned char b)
{
	unsigned char i = 0;
	unsigned char sum_in_pow = 0;
	
	unsigned char sum_in_poly = gf_pow2poly(a) ^ gf_pow2poly(b);
#if 0	
	if(26 == b)
	{
		DEBUG_NOTICE("sum_in_poly: %x = %x + %x | %x ^ %x\n", sum_in_poly, a, b, gf_pow2poly(a), gf_pow2poly(b));
	}
#endif	
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
		DEBUG_NOTICE("div err.\n");
		return 0xFF;
	}

	//DEBUG_NOTICE("div: %x %x\n", a, b);
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
		//DEBUG_NOTICE("quotien is zero: %d %d\n", gf_degree(divisor, len_divisor), gf_degree(dividend, len_dividend));

		return 0;
	}

	memcpy(dividend_tmp, dividend, sizeof(unsigned char) * len_dividend);
	for(i = 0; i < len_dividend; i++)
	{
		locator = gf_degree(dividend_tmp, len_dividend) - gf_degree(divisor, len_divisor);
		factor = gf_div(dividend_tmp[gf_degree(dividend_tmp, len_dividend)], divisor[gf_degree(divisor, len_divisor)]);

		quotien[locator] = factor;
		//DEBUG_NOTICE("quotien: %x %d %x\n", quotien[locator], locator, factor);

		for(j = 0; j < len_divisor; j++)
		{
			factor_rmd = gf_multp(factor, divisor[gf_degree(divisor, len_divisor) - j]);
			locator_rmd = locator + gf_degree(divisor, len_divisor) - j;
			remainder_tmp[locator_rmd] = gf_add(dividend_tmp[locator_rmd], factor_rmd);
		}
#if 0
		DEBUG_NOTICE("remainder_tmp:\n");
		for(k = 0; k < len_remainder; k++)
		{
			DEBUG_NOTICE("%x ", remainder_tmp[k]);
		}
		DEBUG_NOTICE("\n");
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
				//DEBUG_NOTICE("product len err: %d\n", idx);
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
		//DEBUG_NOTICE("unrel_group_seq: %x\n", a[0]);
#if 0
		DEBUG_NOTICE("%x %x %x %x %x %x\n", 
									idx,
									product[idx],
					  	   	 		reg[len_a - 2],
					       	 		pd_tmp,
					       	 		*(a + len_a - 1),
					       	 		*(b + (len_b - 1 - i)));
#endif
		if(0 >= idx)
		{
			//DEBUG_NOTICE("product len err: %d\n", idx);
			break;
		}
		idx = idx - 1;

		for(j = 1; j < (len_a - 1); j++)
		{
			pd_tmp = gf_multp(*(a + len_a - 1 - j), *(b + (len_b - 1 - i)));
			reg[len_a - 1 - j] = gf_add(reg[len_a - 1 - j - 1], pd_tmp);
#if 0
			DEBUG_NOTICE("idx: %d\n", len_a - 1 - j);
			DEBUG_NOTICE("%x %x %x %x %x\n", reg[len_a - 1 - j], reg[len_a - 1 - j - 1], pd_tmp,
								 	   *(a + len_a - 1 - j), *(b + (len_b - 1 - i)));
#endif
		}
		reg[0] = gf_multp(*(a + 0), *(b + (len_b - 1 - i)));
		//DEBUG_NOTICE("%x %x %x\n", reg[0], *(a + 0), *(b + (len_b - 1 - i)));
	}

	for(i = 0; i < len_a; i++)
	{
		pd_tmp = gf_multp(*(a + len_a - 1), 0xFF);
		product[idx] = gf_add(reg[len_a - 2], pd_tmp);
#if 0
		DEBUG_NOTICE("%x %x %x %x %x %x\n", 
									idx,
									product[idx],
					  	   	 		reg[len_a - 2],
					       	 		pd_tmp,
					       	 		*(a + len_a - 1),
					       	 		0xFF);
#endif
		if(0 >= idx)
		{
			//DEBUG_NOTICE("product len err: %d\n", idx);
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

unsigned long long real_combine(unsigned long long n, unsigned long long k)
{
	unsigned long long combine_num = 0;

#if 0//it is useless when values are too large
	int i = 0;
	long tmp_n = 1, tmp_k = 1, tmp_n_k = 1;

	for(i = 1; i < (n + 1); i++)
	{
		tmp_n = tmp_n * i;
	}

	for(i = 1; i < (k + 1); i++)
	{
		tmp_k = tmp_k * i;
	}

	for(i = 1; i < (n - k + 1); i++)
	{
		tmp_n_k = tmp_n_k * i;
	}

	combine_num = tmp_n / tmp_k / tmp_n_k;
#else//fast calculation for finite field
	if(k == (n & k))
	{
		combine_num = 1;
	}
	else
	{
		combine_num = 2;
	}
#endif

	return combine_num;
}

unsigned char gf_real_mutp_ff(unsigned long long n, unsigned char ff)
{
	unsigned char val = 0xFF;

	if(0 != (n % 2))
	{
		val = ff;
	}
	else
	{
		val = 0xFF;
	}

	return val;
}

unsigned char gf_pow_cal(unsigned char ff, unsigned long long n)
{
	unsigned char val = 0xFF;
	if(0xFF == ff)
	{
		if(0 != n)
		{
			return 0xFF;
		}
		else
		{
			return 0x0;
		}
	}

	if(0 <= n)
	{
		val = (ff * n) % (GF_FIELD - 1);
	}
	else
	{
		val = (power_polynomial_table[-n + 1][0]) % (GF_FIELD - 1);
		val = (ff * val) % (GF_FIELD - 1);
	}

	return val;
}
