#ifndef DEBUG_INFO_H
#define DEBUG_INFO_H

#include <stdio.h>

#define CFG_DEBUG_SYS			1
#define CFG_DEBUG_IMPOTANT		0
#define CFG_DEBUG_INFO			0
#define CFG_DEBUG_NOTICE		0

#if (1 == CFG_DEBUG_SYS)
#define DEBUG_SYS			printf
#else
#define DEBUG_SYS			//printf
#endif
#if (1 == CFG_DEBUG_IMPOTANT)
#define DEBUG_IMPOTANT			printf
#else
#define DEBUG_IMPOTANT			//printf
#endif
#if (1 == CFG_DEBUG_INFO)
#define DEBUG_INFO			printf
#else
#define DEBUG_INFO			//printf
#endif
#if (1 == CFG_DEBUG_NOTICE)
#define DEBUG_NOTICE			printf
#else
#define DEBUG_NOTICE			//printf
#endif

#endif
