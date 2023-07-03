#include <iostream>

uint64_t* parallel_sub(uint64_t* res, uint64_t *a, uint64_t* b, int size) {
    for(int i = 0; i < size; i++) {
        res[i] = a[i] - b[i];
    }
    return res;
}

uint64_t* parallel_sub_from_zero(uint64_t* res, uint64_t *a, int start_ind, int end_ind) {
    for(int i = start_ind; i < end_ind; i++) {
        res[i] = (uint64_t) 0 - a[i];
    }
    return res;
}

void parallel_pos_neg_sign_assign(uint64_t *pos_e_contributions, uint64_t *neg_e_contributions, uint64_t *repeated_msb_a,
                                  uint64_t *one_contributions, uint64_t *a, uint64_t *abs_a, uint64_t *msb_a,
                                  uint64_t *pec, uint64_t *nec, int p_role, int n_bits, int start_ind, int end_ind) {
    for(uint32_t i = start_ind; i < end_ind; i++) {
        pos_e_contributions[i * (n_bits + 1)] = a[i];
        neg_e_contributions[i * (n_bits + 1)] = abs_a[i];
        repeated_msb_a[i * (n_bits + 1)] = msb_a[i];
        for (int bi = 0; bi < n_bits; bi++) {
            pos_e_contributions[(i * (n_bits + 1)) + bi + 1] = pec[bi];
            neg_e_contributions[(i * (n_bits + 1)) + bi + 1] = nec[bi];
            one_contributions[(i * n_bits) + bi] = p_role * (((uint64_t) 1) << FRAC);
            repeated_msb_a[(i * (n_bits + 1)) + bi + 1] = msb_a[i];
        }
    }
}

void parallel_exp_sub(uint64_t *pos_e_contributions, uint64_t *neg_e_contributions, uint64_t *repeated_msb_a,
                      uint64_t *one_contributions, uint64_t *a, uint64_t *abs_a, uint64_t *msb_a,
                      int neg_n_bits, int p_role, int n_bits, int start_ind, int end_ind) {
    for(int i = start_ind; i < end_ind; i++) {
        abs_a[i] = (uint64_t) 0 - a[i];
    }

    // compute the possible contribution of positive and negative values
    uint64_t* pec = new uint64_t[n_bits];
    uint64_t* nec = new uint64_t[n_bits];
    if(p_role == P2) {
        for (int i = n_bits - 1; i >= 0; i--) {
            pec[n_bits - i - 1] = convert2uint64(exp(pow(2, i - FRAC)));
            if (i > neg_n_bits - 1) {
                nec[n_bits - i - 1] = (((uint64_t) 1) << FRAC);
            } else {
                nec[n_bits - i - 1] = convert2uint64(1.0 / exp(pow(2, i - FRAC)));
            }
        }
    }
    else {
        for (int i = n_bits - 1; i >= 0; i--) {
            pec[n_bits - i - 1] = 0;
            nec[n_bits - i - 1] = 0;
        }
    }

    for(uint32_t i = start_ind; i < end_ind; i++) {
        pos_e_contributions[i * (n_bits + 1)] = a[i];
        neg_e_contributions[i * (n_bits + 1)] = abs_a[i];
        repeated_msb_a[i * (n_bits + 1)] = msb_a[i];
        for (int bi = 0; bi < n_bits; bi++) {
            pos_e_contributions[(i * (n_bits + 1)) + bi + 1] = pec[bi];
            neg_e_contributions[(i * (n_bits + 1)) + bi + 1] = nec[bi];
            one_contributions[(i * n_bits) + bi] = p_role * (((uint64_t) 1) << FRAC);
            repeated_msb_a[(i * (n_bits + 1)) + bi + 1] = msb_a[i];
        }
    }
}

