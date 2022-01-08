#ifndef AS_DECODING_H
#define AS_DECODING_H

#define RELEX_ORDER				1

#define S_MUL					1
//#define K_M						S_MUL * (S_MUL - 1) * MESSAGE_LEN
#define LAYER_NUM				(MESSAGE_LEN + 1)
//#define LAYER_NUM				2
/*approximate and sufficient size allocation*/
#define TERM_SIZE				((5 * (S_MUL + 1) * CODEWORD_LEN) >> 2)//REAL_SIZE = TERM_SIZE^2
#define POLY_NUM				TERM_SIZE// notice there may be dangerous! POLY_NUM = ROOT_CNT_LAST_LAYER * GF_FIELD
#define ROOT_SIZE				POLY_NUM
#define LEX_TABLE_EXPAND_SIZE	3 / 2

extern unsigned char received_polynomial[CODEWORD_LEN];
extern unsigned char output_polynomial[CODEWORD_LEN];
extern unsigned char decoded_codeword[CODEWORD_LEN];
extern unsigned char decoded_message[MESSAGE_LEN];

extern int chanl_rel_init();
extern int mul_assign();
extern int tao_cal();
extern unsigned char poly_eva(unsigned char *poly, unsigned char poly_len, unsigned char input_val);
extern int re_encoding();
extern int as_decoding();
extern int g_term_malloc();
extern int g_term_destroy();
#endif
