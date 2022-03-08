//
// Created by Mete Akgun on 04.07.20.
//

#ifndef PML_FLIB_H
#define PML_FLIB_H

#include <stdlib.h>
#include <stdint.h>

uint64_t convert2Long(unsigned char **ptr){
    uint64_t val = 0;
    for (int i=56;i>=0;i-=8){
        val=val+((uint64_t)(**ptr)<<i);
        (*ptr)++;
    }
    return val;
}
uint32_t convert2Int(unsigned char **ptr){
    uint32_t val = 0;
    for (int i=24;i>=0;i-=8){
        val=val+((uint64_t)(**ptr)<<i);
        (*ptr)++;
    }
    return val;
}
uint8_t convert2uint8(unsigned char **ptr){
    uint8_t val = (uint8_t)(**ptr);
    (*ptr)++;
    return val;
}
void convert2Array(unsigned char **ptr, uint8_t arr[],int sz){
    for (int i=0;i<sz;i++) {
        arr[i] = (**ptr);
        (*ptr)++;
    }
}
void convert2Array(unsigned char **ptr, uint64_t arr[], int sz){
    for (int i=0;i<sz;i++) {
        arr[i] = (**ptr);
        (*ptr)++;
    }
}
void addVal2CharArray(uint64_t val,unsigned char **ptr){
    for (int i=56;i>=0;i-=8){
        (**ptr)=(val>>i)&0xff;
        (*ptr)++;
    }
}
void addVal2CharArray(uint32_t val,unsigned char **ptr){
    for (int i=24;i>=0;i-=8){
        (**ptr)=(val>>i)&0xff;
        (*ptr)++;
    }
}
void addVal2CharArray(uint8_t val,unsigned char **ptr){
    (**ptr)=(val)&0xff;
    (*ptr)++;
}
void addVal2CharArray(uint8_t val[],unsigned char **ptr, int sz){
    for (int i=0;i<sz;i++){
        (**ptr)=(val[i])&0xff;
        (*ptr)++;
    }
}
void addVal2CharArray(uint64_t *val, unsigned char **ptr, int sz){
    for (int i=0;i<sz;i++){
        addVal2CharArray(val[i],ptr);
    }

    // debora ? did you add this function ?
    /*for (int i=0;i<sz;i++){
        (**ptr)=(val[i])&0xff;
        (*ptr)++;
    }*/
}
uint8_t bit(uint64_t val,uint8_t ind){
    return (val>>ind)&0x1;
}
uint8_t mod(int k, int n) {
    return ((k %= n) < 0) ? k+n : k;
}

double convert2double(uint64_t x, int precision=FRAC) {
    double tmp = 1 << precision;
    if ((int) (x >> 63) == 1) {
        return -1 * ((double) (~x + 1) / tmp);
    } else {
        return ((double) x / tmp);
    }
}
uint64_t convert2uint64(double x, int precision=FRAC) {
    if (x < 0) {
        return (uint64_t) 0 - (uint64_t) floor(abs(x * (1 << precision)));
    } else {
        return (uint64_t) floor(x * (1 << precision));
    }
}


double *convert2double(uint64_t *x, uint32_t size, int precision=FRAC) {
    double *res = new double[size];
    double tmp = 1 << precision;
    for (int i = 0; i < size; i++) {
        if ((int) (x[i] >> 63) == 1) {
            res[i] = -1 * ((double) (~x[i] + 1) / tmp);
        } else {
            res[i] = ((double) x[i] / tmp);
        }
    }
    return res;
}
uint64_t* convert2uint64(double* x, uint32_t size, int precision=FRAC) {
    uint64_t *res = new uint64_t[size];
    double tmp = 1 << precision;
    for (int i = 0; i < size; i++) {
        if (x[i] < 0) {
            res[i] = (uint64_t) 0 - (uint64_t) floor(abs(x[i] * (1 << precision)));
        } else {
            res[i] = (uint64_t) floor(x[i] * (1 << precision));
        }
    }
    return res;
}

/**
 * Get the MoDular Inverse (MDI) of a given number a with specified modulo. For the resulting/returned value b must hold
 *      ab mod(modulo) are congruent to 1.
 * @param a the value for which the modular inverse shall be calculated.
 * The modulo under which a and the inverse are multiplied equal to 1 will always be the ring size.
 * @return the modular inverse of a under the ring size of 16.
 */
uint64_t getModularInverse(uint64_t a){
    // start with 1 because 0 does not have an inverse value.
    for (uint64_t x = 1; x < N; x++){
        // use bitwise and instead of %: faster, and because N = 0xf...
        if (((a & N) * (x & N)) & N == 1) {
            return x;
        }
    }
    return 0;
}

#endif //PML_FLIB_H
