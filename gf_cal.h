#ifndef GF_CAL_H
#define GF_CAL_H

#define GF_Q			 4
#if (3 == GF_Q)
#define GF_FIELD        8
#define MESSAGE_LEN     3
#endif
#if (4 == GF_Q)
#define GF_FIELD        16
#define MESSAGE_LEN     7
#endif
#if (6 == GF_Q)
#define GF_FIELD        64
#define MESSAGE_LEN     21
#endif
#define CODEWORD_LEN    (GF_FIELD - 1)   

extern unsigned char power_polynomial_table[GF_FIELD][2];

extern unsigned char gf_pow2poly(unsigned char val_in_pow);
extern unsigned char gf_poly2pow(unsigned char val_in_poly);
extern unsigned char gf_location(unsigned char val);
extern unsigned char gf_add(unsigned char a, unsigned char b);
extern unsigned char gf_multp(unsigned char a, unsigned char b);
extern unsigned char gf_div(unsigned char a, unsigned char b);
extern unsigned char gf_mod_single_term(unsigned char a, unsigned char b);
extern unsigned char gf_degree(unsigned char* a, unsigned char len_a);
extern unsigned char gf_div_q_r(unsigned char* dividend, unsigned char len_dividend,
								   unsigned char* divisor, unsigned char len_divisor,
								   unsigned char* quotien, unsigned char len_quotien,
								   unsigned char* remainder, unsigned char len_remainder);
extern unsigned char gf_multp_poly(unsigned char* a, unsigned char len_a,
									   unsigned char* b, unsigned char len_b,
									   unsigned char* product, unsigned char len_product);

extern int gf_multp_poly_hw(unsigned char* a, unsigned char len_a,
				 				  unsigned char* b, unsigned char len_b,
				 				  unsigned char* product, unsigned char len_product);
extern unsigned long long real_combine(unsigned long long n, unsigned long long k);
extern unsigned char gf_real_mutp_ff(unsigned long long n, unsigned char ff);
unsigned char gf_pow_cal(unsigned char ff, unsigned long long n);
#endif
