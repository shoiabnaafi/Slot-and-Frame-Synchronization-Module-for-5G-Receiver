#include "corr.h"
#include <stdlib.h>
#include <string.h>

static int correlation_initialized = 0;

// Initialize correlation engine
int correlation_init(void) {
    if (correlation_initialized) {
        return SYNC_SUCCESS;
    }
    
    correlation_initialized = 1;
    return SYNC_SUCCESS;
}

// Cleanup correlation engine
void correlation_cleanup(void) {
    correlation_initialized = 0;
}

// Time-domain cross-correlation implementation
void cross_correlate_time(const complex_t* signal, const complex_t* reference, 
                         int signal_len, int ref_len, float* correlation) {
    
    if (!signal || !reference || !correlation) {
        return;
    }
    
    // Calculate correlation for each possible shift
    int num_lags = signal_len - ref_len + 1;
    
    for (int lag = 0; lag < num_lags; lag++) {
        float real_sum = 0.0f;
        float imag_sum = 0.0f;
        
        // Calculate correlation at this lag
        for (int n = 0; n < ref_len; n++) {
            // Complex correlation: signal * conj(reference)
            real_sum += signal[lag + n].real * reference[n].real + 
                       signal[lag + n].imag * reference[n].imag;
            imag_sum += signal[lag + n].imag * reference[n].real - 
                       signal[lag + n].real * reference[n].imag;
        }
        
        // Store magnitude squared for peak detection
        correlation[lag] = real_sum * real_sum + imag_sum * imag_sum;
    }
}

// Simplified FFT-based correlation (placeholder for full FFT implementation)
void cross_correlate_fft(const complex_t* signal, const complex_t* reference,
                        int signal_len, int ref_len, float* correlation) {
    
    // For a production implementation, you would use a proper FFT library
    // like FFTW, KissFFT, or CMSIS-DSP. Here we fall back to time-domain
    // for simplicity in this demonstration.
    
    cross_correlate_time(signal, reference, signal_len, ref_len, correlation);
}

// Sliding window correlation for real-time processing
int sliding_window_correlate(const complex_t* rx_buffer, int buffer_len,
                            const complex_t* ref_seq, int ref_len,
                            float* correlation_output, int window_step) {
    
    if (!rx_buffer || !ref_seq || !correlation_output) {
        return SYNC_INVALID_PARAMS;
    }
    
    if (buffer_len < ref_len) {
        return SYNC_INVALID_PARAMS;
    }
    
    int num_windows = (buffer_len - ref_len + window_step) / window_step;
    
    for (int w = 0; w < num_windows; w++) {
        int window_start = w * window_step;
        
        if (window_start + ref_len > buffer_len) {
            break;
        }
        
        // Correlate reference with current window
        float real_sum = 0.0f;
        float imag_sum = 0.0f;
        
        for (int n = 0; n < ref_len; n++) {
            real_sum += rx_buffer[window_start + n].real * ref_seq[n].real + 
                       rx_buffer[window_start + n].imag * ref_seq[n].imag;
            imag_sum += rx_buffer[window_start + n].imag * ref_seq[n].real - 
                       rx_buffer[window_start + n].real * ref_seq[n].imag;
        }
        
        correlation_output[w] = sqrtf(real_sum * real_sum + imag_sum * imag_sum);
    }
    
    return num_windows;
}

// Complex number multiplication
void complex_multiply(const complex_t* a, const complex_t* b, complex_t* result) {
    if (!a || !b || !result) {
        return;
    }
    
    result->real = a->real * b->real - a->imag * b->imag;
    result->imag = a->real * b->imag + a->imag * b->real;
}

// Complex conjugate
void complex_conjugate(const complex_t* input, complex_t* output, int len) {
    if (!input || !output) {
        return;
    }
    
    for (int i = 0; i < len; i++) {
        output[i].real = input[i].real;
        output[i].imag = -input[i].imag;
    }
}

// Complex magnitude squared
float complex_magnitude_squared(const complex_t* c) {
    if (!c) {
        return 0.0f;
    }
    
    return c->real * c->real + c->imag * c->imag;
}

// Utility function for normalized cross-correlation
void normalized_cross_correlate(const complex_t* signal, const complex_t* reference,
                               int signal_len, int ref_len, float* correlation) {
    
    if (!signal || !reference || !correlation) {
        return;
    }
    
    // Calculate reference power for normalization
    float ref_power = 0.0f;
    for (int i = 0; i < ref_len; i++) {
        ref_power += complex_magnitude_squared(&reference[i]);
    }
    
    if (ref_power == 0.0f) {
        memset(correlation, 0, (signal_len - ref_len + 1) * sizeof(float));
        return;
    }
    
    int num_lags = signal_len - ref_len + 1;
    
    for (int lag = 0; lag < num_lags; lag++) {
        float real_sum = 0.0f;
        float imag_sum = 0.0f;
        float signal_power = 0.0f;
        
        // Calculate correlation and signal power at this lag
        for (int n = 0; n < ref_len; n++) {
            real_sum += signal[lag + n].real * reference[n].real + 
                       signal[lag + n].imag * reference[n].imag;
            imag_sum += signal[lag + n].imag * reference[n].real - 
                       signal[lag + n].real * reference[n].imag;
            signal_power += complex_magnitude_squared(&signal[lag + n]);
        }
        
        // Normalized correlation
        float corr_magnitude = sqrtf(real_sum * real_sum + imag_sum * imag_sum);
        if (signal_power > 0.0f) {
            correlation[lag] = corr_magnitude / sqrtf(signal_power * ref_power);
        } else {
            correlation[lag] = 0.0f;
        }
    }
}