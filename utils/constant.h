//
// Created by Mete Akgun on 03.07.20.
//

#ifndef PML_CONSTANT_H
#define PML_CONSTANT_H


#define L 64
#define LP 67
#define FRAC 0
#define N 0xffffffffffffffff
#define N1_MASK 0x7fffffffffffffff
#define N1 0x8000000000000000

#define MAXMUL 16384

#define PRECISION 100
#define MAX_SAMPLE 0xfffff
#define MAX_SAMPLE_MASK 0x7ffff

#define MAX_MATRIX_SIZE 1048576


#define DEBUG_FLAG 0


enum role {
    P1, P2, HELPER
};

enum op {
    CORE_MMSB,CORE_END,CORE_MUL,CORE_MMUL,CORE_MUX,CORE_MMUX,CORE_MMC,CORE_MC,CORE_MCMP,CORE_CMP,CORE_MSB, // Core
    AUC_MSB,AUC_TDIV,AUC_MROU,AUC_MDIV //AUC
};


#endif //PML_CONSTANT_H
