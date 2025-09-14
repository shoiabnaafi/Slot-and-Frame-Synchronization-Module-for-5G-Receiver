#ifndef CORR_H
#define CORR_H

#include "sync_config.h"

// Cross-correlation functions
int correlation_init(void);
void correlation_cleanup(void);

// Time-domain cross-correlation (for shorter sequences)
void cross_correlate_time(const complex_t* signal, const complex_t* reference, 
                         int signal_len, int ref_len, float* correlation);

// Frequency-domain cross-correlation using FFT
void cross_correlate_fft(const complex_t* signal, const complex_t* reference,
                        int signal_len, int ref_len, float* correlation);

// Sliding window correlation for real-time processing
int sliding_window_correlate(const complex_t* rx_buffer, int buffer_len,
                            const complex_t* ref_seq, int ref_len,
                            float* correlation_output, int window_step);

// Utility functions for correlation
void complex_multiply(const complex_t* a, const complex_t* b, complex_t* result);
void complex_conjugate(const complex_t* input, complex_t* output, int len);
float complex_magnitude_squared(const complex_t* c);

#endif // CORR_H