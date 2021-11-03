#ifndef ENCODING_H
#define ENCODING_H

extern unsigned char generator_polynomial[CODEWORD_LEN - MESSAGE_LEN + 1];
extern unsigned char message_polynomial[MESSAGE_LEN];
extern unsigned char error_polynomial[CODEWORD_LEN];
extern unsigned char encoded_polynomial[CODEWORD_LEN];

extern unsigned char evaluation_encoding();
extern unsigned char systematic_encoding();

#endif