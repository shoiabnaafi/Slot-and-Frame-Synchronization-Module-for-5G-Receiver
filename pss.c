#include "pss.h"
#include "corr.h"
#include "peak.h"
#include <stdlib.h>
#include <string.h>

// Global PSS sequences storage
complex_t pss_sequences[NUM_PSS_SEQUENCES][PSS_LENGTH];
static int pss_initialized = 0;

// PSS m-sequence initial states for NID2 = 0, 1, 2
static const int pss_init_states[NUM_PSS_SEQUENCES] = {0x01, 0x02, 0x04};

// Initialize PSS detector
int pss_init(void) {
    if (pss_initialized) {
        return SYNC_SUCCESS;
    }
    
    // Generate all 3 PSS sequences
    for (int nid2 = 0; nid2 < NUM_PSS_SEQUENCES; nid2++) {
        if (generate_pss_sequence(nid2, pss_sequences[nid2]) != SYNC_SUCCESS) {
            return SYNC_INVALID_PARAMS;
        }
    }
    
    pss_initialized = 1;
    return SYNC_SUCCESS;
}

// Cleanup PSS detector
void pss_cleanup(void) {
    pss_initialized = 0;
    memset(pss_sequences, 0, sizeof(pss_sequences));
}

// Generate PSS sequence based on 3GPP TS 38.211
int generate_pss_sequence(int cell_id_2, complex_t* pss_seq) {
    if (cell_id_2 < 0 || cell_id_2 >= NUM_PSS_SEQUENCES || !pss_seq) {
        return SYNC_INVALID_PARAMS;
    }
    
    // PSS m-sequence generation
    int x[PSS_LENGTH + 7] = {0}; // Extra space for shift register
    
    // Initialize shift register based on cell_id_2
    x[6] = (pss_init_states[cell_id_2] >> 0) & 1;
    x[5] = (pss_init_states[cell_id_2] >> 1) & 1;
    x[4] = (pss_init_states[cell_id_2] >> 2) & 1;
    x[3] = (pss_init_states[cell_id_2] >> 3) & 1;
    x[2] = (pss_init_states[cell_id_2] >> 4) & 1;
    x[1] = (pss_init_states[cell_id_2] >> 5) & 1;
    x[0] = (pss_init_states[cell_id_2] >> 6) & 1;
    
    // Generate m-sequence: x(i+7) = (x(i+4) + x(i)) mod 2
    for (int i = 0; i < PSS_LENGTH; i++) {
        x[i + 7] = (x[i + 4] + x[i]) & 1;
    }
    
    // Generate PSS sequence: d_pss(n) = 1 - 2*x(m) where m = (n + 43*NID2) mod 127
    for (int n = 0; n < PSS_LENGTH; n++) {
        int m = (n + 43 * cell_id_2) % PSS_LENGTH;
        pss_seq[n].real = 1.0f - 2.0f * x[m + 7];
        pss_seq[n].imag = 0.0f; // PSS is BPSK modulated
    }
    
    return SYNC_SUCCESS;
}

// Detect PSS in received samples
int detect_pss(const complex_t* rx_samples, int num_samples, 
               float* correlation, int* peak_index, float* peak_magnitude) {
    
    if (!rx_samples || !correlation || !peak_index || !peak_magnitude) {
        return SYNC_INVALID_PARAMS;
    }
    
    if (!pss_initialized) {
        return SYNC_INVALID_PARAMS;
    }
    
    float best_correlation = 0.0f;
    int best_peak_idx = -1;
    int best_cell_id_2 = -1;
    
    // Test correlation with all 3 PSS sequences
    for (int nid2 = 0; nid2 < NUM_PSS_SEQUENCES; nid2++) {
        float* temp_correlation = (float*)malloc((num_samples - PSS_LENGTH + 1) * sizeof(float));
        if (!temp_correlation) {
            return SYNC_INVALID_PARAMS;
        }
        
        // Perform cross-correlation
        cross_correlate_time(rx_samples, pss_sequences[nid2], 
                           num_samples, PSS_LENGTH, temp_correlation);
        
        // Find peak in this correlation
        int temp_peak_idx;
        float temp_peak_mag;
        float noise_floor = estimate_noise_floor(temp_correlation, num_samples - PSS_LENGTH + 1);
        float threshold = calculate_adaptive_threshold(temp_correlation, 
                                                     num_samples - PSS_LENGTH + 1,
                                                     noise_floor, THRESHOLD_FACTOR);
        
        int result = find_peak_with_threshold(temp_correlation, 
                                            num_samples - PSS_LENGTH + 1,
                                            threshold, &temp_peak_idx, &temp_peak_mag);
        
        if (result == SYNC_SUCCESS && temp_peak_mag > best_correlation) {
            best_correlation = temp_peak_mag;
            best_peak_idx = temp_peak_idx;
            best_cell_id_2 = nid2;
            
            // Copy the best correlation result
            memcpy(correlation, temp_correlation, 
                   (num_samples - PSS_LENGTH + 1) * sizeof(float));
        }
        
        free(temp_correlation);
    }
    
    if (best_peak_idx >= 0) {
        *peak_index = best_peak_idx;
        *peak_magnitude = best_correlation;
        return best_cell_id_2; // Return cell_id_2 as success indicator
    }
    
    return SYNC_NO_PEAK_FOUND;
}

// Estimate noise floor from correlation output
float estimate_noise_floor(const float* correlation, int len) {
    if (!correlation || len <= 0) {
        return 0.0f;
    }
    
    // Use median-based estimation for robustness
    float* temp_array = (float*)malloc(len * sizeof(float));
    if (!temp_array) {
        return 0.0f;
    }
    
    memcpy(temp_array, correlation, len * sizeof(float));
    
    // Calculate absolute values
    for (int i = 0; i < len; i++) {
        temp_array[i] = fabsf(temp_array[i]);
    }
    
    float noise_floor = estimate_noise_median(temp_array, len);
    
    free(temp_array);
    return noise_floor;
}

// Get cell_id_2 from peak analysis
int get_pss_cell_id_2(int peak_index, const float* correlation) {
    // This is a simplified version - in practice, you would 
    // need to track which PSS sequence gave the best correlation
    // For now, return based on correlation pattern analysis
    
    if (!correlation || peak_index < 0) {
        return -1;
    }
    
    // Placeholder implementation - would need more sophisticated analysis
    return 0; // Default to NID2 = 0
}