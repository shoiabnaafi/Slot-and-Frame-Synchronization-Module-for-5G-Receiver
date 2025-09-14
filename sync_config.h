#ifndef SYNC_CONFIG_H
#define SYNC_CONFIG_H

#include <stdint.h>
#include <math.h>

// 5G NR Synchronization Configuration
#define PSS_LENGTH              127
#define SSS_LENGTH              127
#define BUFFER_SIZE             4096
#define CORRELATION_WINDOW      1024
#define MAX_CELL_ID             1007
#define NUM_PSS_SEQUENCES       3

// Thresholding parameters
#define NOISE_ESTIMATION_LEN    512
#define THRESHOLD_FACTOR        3.5f
#define MIN_PEAK_DISTANCE       64

// Complex number structure
typedef struct {
    float real;
    float imag;
} complex_t;

// Synchronization result structure
typedef struct {
    int is_synchronized;
    int pss_peak_index;
    int sss_peak_index;
    int cell_id;
    float peak_magnitude;
    float snr_estimate;
    int timing_offset;
} sync_result_t;

// Function return codes
#define SYNC_SUCCESS            0
#define SYNC_NO_PEAK_FOUND      -1
#define SYNC_THRESHOLD_NOT_MET  -2
#define SYNC_INVALID_PARAMS     -3

#endif // SYNC_CONFIG_H