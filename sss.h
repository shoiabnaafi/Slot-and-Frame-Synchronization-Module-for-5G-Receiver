#ifndef SSS_H
#define SSS_H

#include "sync_config.h"

// SSS detection functions
int sss_init(void);
void sss_cleanup(void);
void generate_m_sequence(int* x, int init_state);
int generate_sss_sequence(int cell_id_1, int cell_id_2, complex_t* sequence);
int detect_sss(const complex_t* rx_samples, int pss_peak_index, int num_samples,
               float* correlation_output, int* sss_peak_index, int* cell_id_1);

#endif // SSS_H