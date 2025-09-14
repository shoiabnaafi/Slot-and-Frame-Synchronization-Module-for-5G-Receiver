#include "sss.h"
#include "corr.h"
#include "peak.h"
#include <stdlib.h>
#include <string.h>

static int sss_initialized = 0;

// SSS m-sequence lookup table for faster generation
static int sss_m0_table[31];
static int sss_m1_table[31];

// Initialize SSS detector
int sss_init(void) {
    if (sss_initialized) {
        return SYNC_SUCCESS;
    }
    
    // Pre-compute SSS m-sequence tables
    generate_m_sequence(sss_m0_table, 0x1f); // x^5 + x^2 + 1
    generate_m_sequence(sss_m1_table, 0x1f); // x^5 + x^3 + 1
    
    sss_initialized = 1;
    return SYNC_SUCCESS;
}

// Cleanup SSS detector
void sss_cleanup(void) {
    sss_initialized = 0;
    memset(sss_m0_table, 0, sizeof(sss_m0_table));
    memset(sss_m1_table, 0, sizeof(sss_m1_table));
}

// Generate SSS m-sequence
void generate_m_sequence(int* x, int init_state) {
    if (!x) return;
    
    // Initialize shift register
    x[4] = (init_state >> 0) & 1;
    x[3] = (init_state >> 1) & 1;
    x[2] = (init_state >> 2) & 1;
    x[1] = (init_state >> 3) & 1;
    x[0] = (init_state >> 4) & 1;
    
    // Generate 31-length m-sequence
    for (int i = 0; i < 31; i++) {
        x[i + 5] = (x[i + 2] + x[i]) & 1; // x^5 + x^2 + 1
    }
}

// Generate SSS sequence based on 3GPP TS 38.211
int generate_sss_sequence(int cell_id_1, int cell_id_2, complex_t* sss_seq) {
    if (cell_id_1 < 0 || cell_id_1 >= 336 || cell_id_2 < 0 || 
        cell_id_2 >= NUM_PSS_SEQUENCES || !sss_seq) {
        return SYNC_INVALID_PARAMS;
    }
    
    if (!sss_initialized) {
        return SYNC_INVALID_PARAMS;
    }
    
    // SSS generation parameters
    int q_prime = cell_id_1 / 30;
    int q = (cell_id_1 + q_prime * (q_prime + 1) / 2) / 30;
    int m_prime = cell_id_1 + q * (q + 1) / 2;
    int m0 = m_prime % 31;
    int m1 = (m0 + (m_prime / 31) + 1) % 31;
    
    // Generate x0 and x1 sequences
    int x0[SSS_LENGTH];
    int x1[SSS_LENGTH];
    
    // x0(n) = x_tilde((n + m0) mod 31)
    for (int n = 0; n < SSS_LENGTH; n++) {
        x0[n] = sss_m0_table[(n + m0) % 31];
    }
    
    // x1(n) = x_tilde((n + m1) mod 31)  
    for (int n = 0; n < SSS_LENGTH; n++) {
        x1[n] = sss_m1_table[(n + m1) % 31];
    }
    
    // Generate SSS sequence: d_sss(n) = (1 - 2*x0(n)) * (1 - 2*x1(n))
    for (int n = 0; n < SSS_LENGTH; n++) {
        float s0 = 1.0f - 2.0f * x0[n];
        float s1 = 1.0f - 2.0f * x1[n];
        sss_seq[n].real = s0 * s1;
        sss_seq[n].imag = 0.0f; // SSS is BPSK modulated
    }
    
    return SYNC_SUCCESS;
}

// Detect SSS in received samples around PSS position
int detect_sss(const complex_t* rx_samples, int pss_position, int num_samples,
               float* correlation, int* sss_peak_index, int* cell_id_1) {
    
    if (!rx_samples || !correlation || !sss_peak_index || !cell_id_1) {
        return SYNC_INVALID_PARAMS;
    }
    
    if (!sss_initialized) {
        return SYNC_INVALID_PARAMS;
    }
    
    // Search window around expected SSS position
    // SSS is typically located before PSS in the same symbol
    int search_start = pss_position - SSS_LENGTH - 32; // Guard samples
    int search_end = pss_position + SSS_LENGTH + 32;
    
    if (search_start < 0) search_start = 0;
    if (search_end >= num_samples) search_end = num_samples - SSS_LENGTH;
    
    float best_correlation = 0.0f;
    int best_sss_index = -1;
    int best_cell_id_1 = -1;
    
    // Try different cell_id_1 values (simplified search)
    for (int test_cell_id_1 = 0; test_cell_id_1 < 168; test_cell_id_1 += 8) {
        complex_t test_sss[SSS_LENGTH];
        
        // Generate test SSS sequence (assume cell_id_2 = 0 for now)
        if (generate_sss_sequence(test_cell_id_1, 0, test_sss) != SYNC_SUCCESS) {
            continue;
        }
        
        // Correlate in search window
        for (int pos = search_start; pos <= search_end; pos++) {
            float corr_sum = 0.0f;
            
            for (int n = 0; n < SSS_LENGTH; n++) {
                float real_corr = rx_samples[pos + n].real * test_sss[n].real +
                                rx_samples[pos + n].imag * test_sss[n].imag;
                corr_sum += real_corr;
            }
            
            if (corr_sum > best_correlation) {
                best_correlation = corr_sum;
                best_sss_index = pos;
                best_cell_id_1 = test_cell_id_1;
            }
        }
    }
    
    // Validate detection
    float noise_floor = estimate_noise_mean(correlation, search_end - search_start + 1);
    float threshold = noise_floor * THRESHOLD_FACTOR;
    
    if (best_correlation > threshold) {
        *sss_peak_index = best_sss_index;
        *cell_id_1 = best_cell_id_1;
        
        // Fill correlation array for analysis
        for (int i = 0; i <= search_end - search_start; i++) {
            correlation[i] = (i == best_sss_index - search_start) ? 
                           best_correlation : noise_floor;
        }
        
        return SYNC_SUCCESS;
    }
    
    return SYNC_NO_PEAK_FOUND;
}

// Refine timing using SSS correlation
int refine_timing_with_sss(const complex_t* rx_samples, int pss_pos, 
                          int* refined_offset) {
    
    if (!rx_samples || !refined_offset) {
        return SYNC_INVALID_PARAMS;
    }
    
    // Simple refinement: look for SSS peak around expected position
    int search_window = 16; // +/- 8 samples
    int best_offset = 0;
    float best_metric = 0.0f;
    
    for (int offset = -search_window/2; offset <= search_window/2; offset++) {
        int test_pos = pss_pos + offset;
        if (test_pos < 0 || test_pos + SSS_LENGTH >= BUFFER_SIZE) {
            continue;
        }
        
        // Calculate simple correlation metric
        float metric = 0.0f;
        for (int n = 0; n < SSS_LENGTH/4; n++) { // Use subset for speed
            metric += rx_samples[test_pos + n].real * rx_samples[test_pos + n].real +
                     rx_samples[test_pos + n].imag * rx_samples[test_pos + n].imag;
        }
        
        if (metric > best_metric) {
            best_metric = metric;
            best_offset = offset;
        }
    }
    
    *refined_offset = best_offset;
    return SYNC_SUCCESS;
}

// Calculate cell_id_1 from SSS correlation analysis
int calculate_cell_id_1(const float* sss_correlation, int peak_pos) {
    if (!sss_correlation || peak_pos < 0) {
        return -1;
    }
    
    // Simplified cell_id_1 calculation
    // In practice, this would involve more sophisticated analysis
    // of the SSS correlation pattern
    
    // Placeholder implementation
    return 0; // Default cell_id_1
}