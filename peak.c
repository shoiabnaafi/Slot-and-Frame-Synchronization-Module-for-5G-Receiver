#include "peak.h"
#include <stdlib.h>
#include <string.h>

static int peak_detection_initialized = 0;

// Initialize peak detection
int peak_detection_init(void) {
    if (peak_detection_initialized) {
        return SYNC_SUCCESS;
    }
    
    peak_detection_initialized = 1;
    return SYNC_SUCCESS;
}

// Cleanup peak detection
void peak_detection_cleanup(void) {
    peak_detection_initialized = 0;
}

// Main peak detection function with thresholding
int find_peak_with_threshold(const float* correlation, int len, 
                           float threshold, int* peak_index, 
                           float* peak_magnitude) {
    
    if (!correlation || !peak_index || !peak_magnitude || len <= 0) {
        return SYNC_INVALID_PARAMS;
    }
    
    float max_value = 0.0f;
    int max_index = -1;
    
    // Find maximum value above threshold
    for (int i = 0; i < len; i++) {
        if (correlation[i] > threshold && correlation[i] > max_value) {
            max_value = correlation[i];
            max_index = i;
        }
    }
    
    if (max_index >= 0) {
        *peak_index = max_index;
        *peak_magnitude = max_value;
        return SYNC_SUCCESS;
    }
    
    return SYNC_NO_PEAK_FOUND;
}

// Calculate adaptive threshold based on noise statistics
float calculate_adaptive_threshold(const float* data, int len, 
                                 float noise_floor, float factor) {
    
    if (!data || len <= 0) {
        return 0.0f;
    }
    
    // Use noise floor with multiplicative factor
    float threshold = noise_floor * factor;
    
    // Optional: add mean-based component
    float mean = estimate_noise_mean(data, len);
    float std = estimate_noise_std(data, len, mean);
    
    // Adaptive threshold: noise_floor * factor + mean + N*std
    threshold = fmaxf(threshold, mean + 2.0f * std);
    
    return threshold;
}

// Estimate noise mean
float estimate_noise_mean(const float* data, int len) {
    if (!data || len <= 0) {
        return 0.0f;
    }
    
    double sum = 0.0;
    for (int i = 0; i < len; i++) {
        sum += fabsf(data[i]); // Use absolute values
    }
    
    return (float)(sum / len);
}

// Estimate noise standard deviation
float estimate_noise_std(const float* data, int len, float mean) {
    if (!data || len <= 1) {
        return 0.0f;
    }
    
    double sum_sq_diff = 0.0;
    for (int i = 0; i < len; i++) {
        float diff = fabsf(data[i]) - mean;
        sum_sq_diff += diff * diff;
    }
    
    return sqrtf((float)(sum_sq_diff / (len - 1)));
}

// Estimate noise median (more robust than mean)
float estimate_noise_median(const float* data, int len) {
    if (!data || len <= 0) {
        return 0.0f;
    }
    
    // Create a copy for sorting
    float* temp_data = (float*)malloc(len * sizeof(float));
    if (!temp_data) {
        return 0.0f;
    }
    
    // Copy absolute values
    for (int i = 0; i < len; i++) {
        temp_data[i] = fabsf(data[i]);
    }
    
    // Sort the array
    sort_array(temp_data, len);
    
    float median;
    if (len % 2 == 0) {
        median = (temp_data[len/2 - 1] + temp_data[len/2]) / 2.0f;
    } else {
        median = temp_data[len/2];
    }
    
    free(temp_data);
    return median;
}

// Validate peak (check for spurious detections)
int validate_peak(const float* correlation, int peak_index, 
                 int min_distance, float min_ratio) {
    
    if (!correlation || peak_index < 0) {
        return 0; // Invalid
    }
    
    float peak_value = correlation[peak_index];
    
    // Check nearby values to ensure this is a true peak
    int start = (peak_index >= min_distance) ? peak_index - min_distance : 0;
    int end = peak_index + min_distance;
    
    for (int i = start; i <= end; i++) {
        if (i != peak_index && correlation[i] > peak_value / min_ratio) {
            return 0; // Nearby peak too strong, might be spurious
        }
    }
    
    return 1; // Valid peak
}

// Find multiple peaks above threshold
int find_multiple_peaks(const float* correlation, int len, 
                       float threshold, int* peak_indices, 
                       float* peak_magnitudes, int max_peaks) {
    
    if (!correlation || !peak_indices || !peak_magnitudes || 
        len <= 0 || max_peaks <= 0) {
        return 0;
    }
    
    int num_peaks = 0;
    
    for (int i = 1; i < len - 1 && num_peaks < max_peaks; i++) {
        // Check if this is a local maximum above threshold
        if (correlation[i] > threshold &&
            correlation[i] > correlation[i-1] &&
            correlation[i] > correlation[i+1]) {
            
            // Validate the peak
            if (validate_peak(correlation, i, MIN_PEAK_DISTANCE, 2.0f)) {
                peak_indices[num_peaks] = i;
                peak_magnitudes[num_peaks] = correlation[i];
                num_peaks++;
            }
        }
    }
    
    // Sort peaks by magnitude (descending)
    for (int i = 0; i < num_peaks - 1; i++) {
        for (int j = i + 1; j < num_peaks; j++) {
            if (peak_magnitudes[j] > peak_magnitudes[i]) {
                // Swap magnitudes
                float temp_mag = peak_magnitudes[i];
                peak_magnitudes[i] = peak_magnitudes[j];
                peak_magnitudes[j] = temp_mag;
                
                // Swap indices
                int temp_idx = peak_indices[i];
                peak_indices[i] = peak_indices[j];
                peak_indices[j] = temp_idx;
            }
        }
    }
    
    return num_peaks;
}

// Simple quicksort implementation for array sorting
void sort_array(float* arr, int len) {
    if (!arr || len <= 1) {
        return;
    }
    
    // Simple bubble sort for small arrays (demo purposes)
    for (int i = 0; i < len - 1; i++) {
        for (int j = 0; j < len - i - 1; j++) {
            if (arr[j] > arr[j + 1]) {
                float temp = arr[j];
                arr[j] = arr[j + 1];
                arr[j + 1] = temp;
            }
        }
    }
}

// Calculate SNR estimate from peak and noise floor
float calculate_snr_estimate(float peak_magnitude, float noise_floor) {
    if (noise_floor <= 0.0f) {
        return 0.0f;
    }
    
    // SNR in dB: 10 * log10(signal_power / noise_power)
    float snr_linear = peak_magnitude / noise_floor;
    float snr_db = 10.0f * log10f(snr_linear);
    
    return snr_db;
}

// Advanced peak detection with hysteresis
int find_peak_with_hysteresis(const float* correlation, int len,
                             float high_threshold, float low_threshold,
                             int* peak_index, float* peak_magnitude) {
    
    if (!correlation || !peak_index || !peak_magnitude || len <= 0) {
        return SYNC_INVALID_PARAMS;
    }
    
    if (high_threshold < low_threshold) {
        return SYNC_INVALID_PARAMS;
    }
    
    int in_peak = 0;
    int current_peak_start = -1;
    float current_peak_max = 0.0f;
    int current_peak_max_idx = -1;
    
    float best_peak_magnitude = 0.0f;
    int best_peak_index = -1;
    
    for (int i = 0; i < len; i++) {
        if (!in_peak && correlation[i] > high_threshold) {
            // Start of new peak
            in_peak = 1;
            current_peak_start = i;
            current_peak_max = correlation[i];
            current_peak_max_idx = i;
        } else if (in_peak) {
            // Update peak maximum
            if (correlation[i] > current_peak_max) {
                current_peak_max = correlation[i];
                current_peak_max_idx = i;
            }
            
            // Check for end of peak
            if (correlation[i] < low_threshold) {
                // End of peak - evaluate it
                if (current_peak_max > best_peak_magnitude) {
                    best_peak_magnitude = current_peak_max;
                    best_peak_index = current_peak_max_idx;
                }
                in_peak = 0;
            }
        }
    }
    
    // Handle case where peak extends to end of data
    if (in_peak && current_peak_max > best_peak_magnitude) {
        best_peak_magnitude = current_peak_max;
        best_peak_index = current_peak_max_idx;
    }
    
    if (best_peak_index >= 0) {
        *peak_index = best_peak_index;
        *peak_magnitude = best_peak_magnitude;
        return SYNC_SUCCESS;
    }
    
    return SYNC_NO_PEAK_FOUND;
}