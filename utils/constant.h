//
// Created by Mete Akgun on 03.07.20.
//

#ifndef PML_CONSTANT_H
#define PML_CONSTANT_H

#define L_BIT 64
#define LP 67
#define FRAC 0
#define RING_N 0xffffffffffffffff  // ring size
#define N1_MASK 0x7fffffffffffffff
#define N1 0x8000000000000000
#define EVEN_MASK 0xfffffffffffffff7

#define MAXMUL 16384

#define PRECISION 100
#define MAX_SAMPLE 0xfffff
#define MAX_SAMPLE_MASK 0x7ffff

#define MAX_MATRIX_SIZE 20000
#define BUFFER_SIZE 40000000

#define DEBUG_FLAG 0


// constants for INVSQRT
#define ORTHMASK 0x3e000
#define MAX_DELTA 10
#define MIN_DELTA 1
#define MAXSCALAR 0xfffff
#define MAXA 0x3fffff

enum role {
    P1, P2, HELPER
};

enum op {
    // Core
    CORE_MMSB,CORE_END,CORE_MUL,CORE_MMUL,CORE_MUX,CORE_MMUX,CORE_MMC,CORE_MC,CORE_MCMP,CORE_CMP,CORE_MSB,
    CORE_EXP, CORE_MEXP,CORE_DP,CORE_MDP,CORE_MATMATMUL,CORE_MMATMATMUL,CORE_MATVECMUL,CORE_MMATVECMUL, CORE_DIV,
    // AUC
    AUC_MSB,AUC_TDIV,AUC_MROU,AUC_MDIV,
    // CNN                                 CL1: 1 image, several kernel; CL2: several images, 1 kernel
    CNN_MAX, CNN_MMAX, CNN_RELU, CNN_DRLU, CNN_CL1, CNN_CL2, CNN_MDI,
    // RKN
    RKN_EIG, RKN_MEIG, RKN_GM2KM, RKN_INVSQRT, RKN_MINVSQRT
};


#endif //PML_CONSTANT_H
