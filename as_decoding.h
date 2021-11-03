#ifndef AS_DECODING_H
#define AS_DECODING_H

extern unsigned char received_polynomial[CODEWORD_LEN];
extern unsigned char output_polynomial[CODEWORD_LEN];

extern int chanl_rel_init();
extern int mul_assign();
extern int tao_cal();
extern int re_encoding();
extern int as_decoding();
#endif
