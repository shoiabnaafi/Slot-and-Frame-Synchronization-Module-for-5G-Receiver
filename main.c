#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "sync_config.h"
#include "pss.h"
#include "sss.h"
#include "corr.h"
#include "peak.h"
#include "timing.h"

// Function prototypes
int initialize_sync_system(void);
void cleanup_sync_system(void);
int generate_test_signal(complex_t* buffer, int buffer_len);
int add_noise_to_signal(complex_t* buffer, int buffer_len, float noise_power);
void print_sync_results(const sync_result_t* result);
int run_synchronization_test(void);

// Main synchronization processing function
int process_sync_frame(const complex_t* rx_samples, int num_samples, 
                      sync_result_t* result) {
    
    if (!rx_samples || !result || num_samples < PSS_LENGTH + SSS_LENGTH) {
        return SYNC_INVALID_PARAMS;
    }
    
    // Initialize result structure
    memset(result, 0, sizeof(sync_result_t));
    
    // Allocate correlation buffer
    int correlation_len = num_samples - PSS_LENGTH + 1;
    float* pss_correlation = (float*)malloc(correlation_len * sizeof(float));
    if (!pss_correlation) {
        return SYNC_INVALID_PARAMS;
    }
    
    // Step 1: PSS Detection
    printf("Step 1: Detecting PSS...\n");
    int pss_peak_index;
    float pss_peak_magnitude;
    
    int pss_result = detect_pss(rx_samples, num_samples, pss_correlation, 
                               &pss_peak_index, &pss_peak_magnitude);
    
    if (pss_result < 0) {
        printf("PSS detection failed: %d\n", pss_result);
        free(pss_correlation);
        return pss_result;
    }
    
    printf("PSS detected at index %d with magnitude %.3f (Cell ID 2: %d)\n", 
           pss_peak_index, pss_peak_magnitude, pss_result);
    
    // Store PSS results
    result->pss_peak_index = pss_peak_index;
    result->peak_magnitude = pss_peak_magnitude;
    int cell_id_2 = pss_result;
    
    // Step 2: SSS Detection
    printf("Step 2: Detecting SSS...\n");
    float* sss_correlation = (float*)malloc(correlation_len * sizeof(float));
    if (!sss_correlation) {
        free(pss_correlation);
        return SYNC_INVALID_PARAMS;
    }
    
    int sss_peak_index;
    int cell_id_1;
    
    int sss_result = detect_sss(rx_samples, pss_peak_index, num_samples,
                               sss_correlation, &sss_peak_index, &cell_id_1);
    
    if (sss_result == SYNC_SUCCESS) {
        printf("SSS detected at index %d (Cell ID 1: %d)\n", 
               sss_peak_index, cell_id_1);
        result->sss_peak_index = sss_peak_index;
        result->cell_id = cell_id_1 * 3 + cell_id_2; // Combined cell ID
    } else {
        printf("SSS detection failed: %d\n", sss_result);
        // Continue with PSS-only synchronization
        result->sss_peak_index = -1;
        result->cell_id = cell_id_2;
    }
    
    // Step 3: Timing Offset Correction
    printf("Step 3: Calculating timing offset...\n");
    int timing_offset;
    int slot_start;
    
    if (detect_slot_boundary(rx_samples, pss_peak_index, &slot_start) == SYNC_SUCCESS) {
        timing_offset = pss_peak_index - slot_start;
        printf("Slot boundary detected. Timing offset: %d samples\n", timing_offset);
    } else {
        timing_offset = 0;
        printf("Using PSS position as timing reference\n");
    }
    
    result->timing_offset = timing_offset;
    
    // Step 4: SNR Estimation
    float noise_floor = estimate_noise_floor(pss_correlation, correlation_len);
    result->snr_estimate = calculate_snr_estimate(pss_peak_magnitude, noise_floor);
    
    printf("Estimated SNR: %.2f dB\n", result->snr_estimate);
    
    // Mark as synchronized if we have good PSS detection
    result->is_synchronized = 1;
    
    free(pss_correlation);
    free(sss_correlation);
    
    return SYNC_SUCCESS;
}

// Generate test signal with PSS/SSS
int generate_test_signal(complex_t* buffer, int buffer_len) {
    if (!buffer || buffer_len < BUFFER_SIZE) {
        return SYNC_INVALID_PARAMS;
    }
    
    // Initialize buffer with zeros
    memset(buffer, 0, buffer_len * sizeof(complex_t));
    
    // Insert PSS at position 1000
    int pss_position = 1000;
    if (pss_position + PSS_LENGTH < buffer_len) {
        // Use cell_id_2 = 1 for test
        complex_t test_pss[PSS_LENGTH];
        generate_pss_sequence(1, test_pss);
        
        memcpy(&buffer[pss_position], test_pss, PSS_LENGTH * sizeof(complex_t));
        printf("Test PSS inserted at position %d\n", pss_position);
    }
    
    // Insert SSS nearby (simplified placement)
    int sss_position = pss_position - 200;
    if (sss_position >= 0 && sss_position + SSS_LENGTH < buffer_len) {
        complex_t test_sss[SSS_LENGTH];
        generate_sss_sequence(100, 1, test_sss); // Test cell_id_1 = 100, cell_id_2 = 1
        
        memcpy(&buffer[sss_position], test_sss, SSS_LENGTH * sizeof(complex_t));
        printf("Test SSS inserted at position %d\n", sss_position);
    }
    
    return SYNC_SUCCESS;
}

// Add AWGN to test signal
int add_noise_to_signal(complex_t* buffer, int buffer_len, float noise_power) {
    if (!buffer || buffer_len <= 0) {
        return SYNC_INVALID_PARAMS;
    }
    
    srand((unsigned int)time(NULL));
    
    for (int i = 0; i < buffer_len; i++) {
        // Generate Gaussian noise (Box-Muller transform)
        static int has_spare = 0;
        static double spare;
        
        if (has_spare) {
            has_spare = 0;
            buffer[i].imag += (float)(spare * sqrt(noise_power));
        } else {
            has_spare = 1;
            double u = ((double)rand() / RAND_MAX) * 2.0 - 1.0;
            double v = ((double)rand() / RAND_MAX) * 2.0 - 1.0;
            double mag = u * u + v * v;
            
            if (mag < 1.0 && mag > 0.0) {
                mag = sqrt(-2.0 * log(mag) / mag);
                spare = v * mag;
                buffer[i].real += (float)(u * mag * sqrt(noise_power));
            }
        }
    }
    
    return SYNC_SUCCESS;
}

// Print synchronization results
void print_sync_results(const sync_result_t* result) {
    if (!result) {
        return;
    }
    
    printf("\n=== Synchronization Results ===\n");
    printf("Synchronized: %s\n", result->is_synchronized ? "YES" : "NO");
    
    if (result->is_synchronized) {
        printf("PSS Peak Index: %d\n", result->pss_peak_index);
        printf("PSS Peak Magnitude: %.3f\n", result->peak_magnitude);
        
        if (result->sss_peak_index >= 0) {
            printf("SSS Peak Index: %d\n", result->sss_peak_index);
        }
        
        printf("Cell ID: %d\n", result->cell_id);
        printf("Timing Offset: %d samples\n", result->timing_offset);
        printf("SNR Estimate: %.2f dB\n", result->snr_estimate);
    }
    printf("===============================\n");
}

// Initialize all synchronization subsystems
int initialize_sync_system(void) {
    printf("Initializing synchronization system...\n");
    
    if (pss_init() != SYNC_SUCCESS) {
        printf("Failed to initialize PSS detector\n");
        return SYNC_INVALID_PARAMS;
    }
    
    if (sss_init() != SYNC_SUCCESS) {
        printf("Failed to initialize SSS detector\n");
        return SYNC_INVALID_PARAMS;
    }
    
    if (correlation_init() != SYNC_SUCCESS) {
        printf("Failed to initialize correlation engine\n");
        return SYNC_INVALID_PARAMS;
    }
    
    if (peak_detection_init() != SYNC_SUCCESS) {
        printf("Failed to initialize peak detection\n");
        return SYNC_INVALID_PARAMS;
    }
    
    if (timing_init() != SYNC_SUCCESS) {
        printf("Failed to initialize timing correction\n");
        return SYNC_INVALID_PARAMS;
    }
    
    printf("Synchronization system initialized successfully\n");
    return SYNC_SUCCESS;
}

// Cleanup all synchronization subsystems
void cleanup_sync_system(void) {
    printf("Cleaning up synchronization system...\n");
    
    pss_cleanup();
    sss_cleanup();
    correlation_cleanup();
    peak_detection_cleanup();
    timing_cleanup();
    
    printf("Synchronization system cleaned up\n");
}

// Run complete synchronization test
int run_synchronization_test(void) {
    printf("\n=== Running 5G Slot/Frame Synchronization Test ===\n");
    
    // Allocate test buffer
    complex_t* test_buffer = (complex_t*)malloc(BUFFER_SIZE * sizeof(complex_t));
    if (!test_buffer) {
        printf("Failed to allocate test buffer\n");
        return SYNC_INVALID_PARAMS;
    }
    
    // Generate test signal with PSS and SSS
    if (generate_test_signal(test_buffer, BUFFER_SIZE) != SYNC_SUCCESS) {
        printf("Failed to generate test signal\n");
        free(test_buffer);
        return SYNC_INVALID_PARAMS;
    }
    
    // Add noise to make it realistic
    float snr_db = 10.0f; // 10 dB SNR
    float noise_power = powf(10.0f, -snr_db / 10.0f);
    add_noise_to_signal(test_buffer, BUFFER_SIZE, noise_power);
    printf("Added noise for %.1f dB SNR\n", snr_db);
    
    // Run synchronization
    sync_result_t result;
    int sync_status = process_sync_frame(test_buffer, BUFFER_SIZE, &result);
    
    if (sync_status == SYNC_SUCCESS) {
        print_sync_results(&result);
        
        // Test timing correction
        if (result.timing_offset != 0) {
            complex_t* corrected_buffer = (complex_t*)malloc(BUFFER_SIZE * sizeof(complex_t));
            if (corrected_buffer) {
                printf("\nApplying timing correction...\n");
                correct_timing_offset(test_buffer, corrected_buffer, 
                                    BUFFER_SIZE, result.timing_offset);
                printf("Timing correction applied\n");
                free(corrected_buffer);
            }
        }
        
    } else {
        printf("Synchronization failed with error: %d\n", sync_status);
    }
    
    free(test_buffer);
    return sync_status;
}

// Main function
int main(int argc, char* argv[]) {
    printf("5G NR Slot and Frame Synchronization Module\n");
    printf("==========================================\n");
    
    // Initialize system
    if (initialize_sync_system() != SYNC_SUCCESS) {
        printf("System initialization failed\n");
        return EXIT_FAILURE;
    }
    
    // Run test
    int result = run_synchronization_test();
    
    // Cleanup
    cleanup_sync_system();
    
    if (result == SYNC_SUCCESS) {
        printf("\nTest completed successfully!\n");
        return EXIT_SUCCESS;
    } else {
        printf("\nTest failed with error: %d\n", result);
        return EXIT_FAILURE;
    }
}