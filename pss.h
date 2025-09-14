#ifndef PSS_H
#define PSS_H

#include "sync_config.h"

// PSS detection functions
int pss_init(void);
void pss_cleanup(void);
int generate_pss_sequence(int n_id_2, complex_t* sequence);
int detect_pss(const complex_t* rx_samples, int num_samples, float* correlation_output,
               int* pss_peak_index, float* pss_peak_magnitude);

#endif // PSS_H