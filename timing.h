#ifndef TIMING_H
#define TIMING_H

#include "sync_config.h"

// Timing correction and synchronization functions
int timing_init(void);
void timing_cleanup(void);

// Main timing correction function
int correct_timing_offset(const complex_t* input_samples, 
                          complex_t* output_samples,
                          int num_samples, int timing_offset);

// Fractional timing offset correction (for sub-sample precision)
int fractional_timing_correction(const complex_t* input_samples,
                                complex_t* output_samples,
                                int num_samples, float fractional_offset);

// Timing offset estimation refinement
int estimate_fine_timing_offset(const complex_t* rx_samples,
                               const complex_t* reference_seq,
                               int coarse_offset, int search_window,
                               float* fine_offset);

// Slot and frame boundary detection
int detect_slot_boundary(const complex_t* samples, int pss_position,
                         int* slot_start_index);
int detect_frame_boundary(const complex_t* samples, int sss_position,
                          int* frame_start_index);

// Interpolation functions for fractional delay
void linear_interpolate(const complex_t* input, complex_t* output,
                        int len, float delay);
void cubic_interpolate(const complex_t* input, complex_t* output,
                       int len, float delay);

// Timing tracking and adjustment
typedef struct {
    float accumulated_offset;
    int sample_count;
    float drift_estimate;
    int last_correction_sample;
} timing_tracker_t;

int init_timing_tracker(timing_tracker_t* tracker);
int update_timing_tracker(timing_tracker_t* tracker, float measured_offset);
int get_timing_correction(timing_tracker_t* tracker, int current_sample);

float calculate_snr_estimate(float peak_magnitude, float noise_floor);
float calculate_fractional_correlation(const complex_t* signal, const complex_t* reference, int signal_len, int ref_len, float fractional_offset);

#endif // TIMING_H