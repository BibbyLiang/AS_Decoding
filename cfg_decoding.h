#ifndef CFG_DECODING_H
#define CFG_DECODING_H

#define GF_Q			 8
#define GF_CAL_COUNT	 0

#define SYS_ENC		1

#define TEST_MODE	0

#define EARLY_TERMINATION		1
#define EARLY_TERMINATION_NUM	200
#define OUTPUT_LOG				1

#define RELEX_ORDER				1
#define SIMPLE_ASD				0
#define RE_ENCODING				0
#define RECUR_RR				1
#define DYNAMIC_MEM				1
#define REDUNDANT_SIZE			0
#define DYNAMIC_TERM			1
#if (1 == DYNAMIC_TERM)
#define DYNAMIC_TERM_ITER		1
#define DYNAMIC_TERM_X			9/8
#define DYNAMIC_TERM_Y			10/9
#endif

#define S_MUL					5
#define LEX_TABLE_EXPAND_SIZE	4

#if (1 == RE_ENCODING)
#define Y_WEIGHT				(-1)
#else
#define Y_WEIGHT				(MESSAGE_LEN - 1)
#endif

#endif
