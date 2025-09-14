#ifndef PEAK_H
#define PEAK_H

#include "sync_config.h"

// Peak detection and thresholding functions
int peak_detection_init(void);
void peak_detection_cleanup(void);

// Main peak detection function
int find_peak_with_threshold(const float* correlation, int len, 
                           float threshold, int* peak_index, 
                           float* peak_magnitude);

// Adaptive threshold calculation
float calculate_adaptive_threshold(const float* data, int len, float noise_floor, float threshold_factor);

// Noise floor estimation methods
float estimate_noise_mean(const float* data, int len);
float estimate_noise_std(const float* data, int len, float mean);
float estimate_noise_median(const float* data, int len);
float estimate_noise_floor(const float* correlation, int len);

// Peak validation and filtering
int validate_peak(const float* correlation, int peak_index, 
                 int min_distance, float min_ratio);
int find_multiple_peaks(const float* correlation, int len, 
                       float threshold, int* peak_indices, 
                       float* peak_magnitudes, int max_peaks);

// Statistical utilities
void sort_array(float* arr, int len);
float calculate_snr_estimate(float peak_magnitude, float noise_floor);

#endif // PEAK_H