#ifndef MOD_H
#define MOD_H

#define BITS_PER_SYMBOL_BPSK	1

extern int bpsk_mod(unsigned char *input_seq, 
				 		unsigned int input_len,
				 		float **output_seq,
				 		unsigned int output_len);

extern int bpsk_demod(float **input_seq, 
				    		unsigned int input_len,
				    		unsigned char *output_seq,
				    		unsigned int output_len);

#endif