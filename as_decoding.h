#ifndef AS_DECODING_H
#define AS_DECODING_H

#define RELEX_ORDER				1
#define SIMPLE_ASD				0
#define RE_ENCODING				0
#define RECUR_RR				1

#define S_MUL					1
//#define K_M						S_MUL * (S_MUL - 1) * MESSAGE_LEN
#define LAYER_NUM				(MESSAGE_LEN + 1)
//#define LAYER_NUM				2
/*approximate and sufficient size allocation*/
#define TERM_SIZE				((5 * (S_MUL + 1) * CODEWORD_LEN) >> 2)//REAL_SIZE = TERM_SIZE^2
#if (0 == RECUR_RR)
#define POLY_NUM				TERM_SIZE * TERM_SIZE * LAYER_NUM// notice there may be dangerous! POLY_NUM = ROOT_CNT_LAST_LAYER * GF_FIELD
#define ROOT_SIZE				POLY_NUM
#else
#define POLY_NUM				1// notice there may be dangerous! POLY_NUM = ROOT_CNT_LAST_LAYER * GF_FIELD
#define ROOT_SIZE				1
#endif
#define LEX_TABLE_EXPAND_SIZE	2

extern unsigned char received_polynomial[CODEWORD_LEN];
extern unsigned char output_polynomial[CODEWORD_LEN];
extern unsigned char decoded_codeword[CODEWORD_LEN];
extern unsigned char decoded_message[MESSAGE_LEN];

extern unsigned long long err_num;
extern unsigned char decoding_ok_flag;
extern unsigned long long weight_stored;
extern unsigned long long hamm_distance_debug;
extern unsigned long long rr_err;

extern int chanl_rel_init();
extern int chnl_rel_cal();
extern int mul_assign();
extern int tao_cal();
extern unsigned char poly_eva(unsigned char *poly, unsigned long long poly_len, unsigned char input_val);
extern int re_encoding();
extern int as_decoding();
extern int g_term_malloc();
extern int g_term_destroy();
extern int dfs_rr_recur(unsigned char *g_c_q,
					unsigned char *g_c_0_y,
					unsigned long long l);
extern unsigned long long hamm_distance_code_cal(unsigned char *a,
									  					  unsigned char *b,
									  					  unsigned long long len);
extern unsigned long long hamm_distance_bit_cal(unsigned char *a,
									  			  unsigned char *b,
									  			  unsigned long long len);

#endif
