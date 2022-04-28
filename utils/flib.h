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
        arr[i] = convert2Long(ptr);
    }
}
void convert2Array(unsigned char **ptr, uint64_t *&arr, uint32_t size){
    // Recover a one dimensional dynamic array
    arr = new uint64_t[size];
    for( uint32_t i = 0; i < size; i++) {
        arr[i] = convert2Long(ptr);
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
}
void addArray2CharArray( uint64_t **val, unsigned char **ptr, uint32_t n_row, uint32_t n_col){
    // Add uint64_t vals in **val to the buffer to send
    for( uint32_t i = 0; i < n_row; i++) {
        for( uint32_t j = 0; j < n_col; j++) {
            addVal2CharArray(val[i][j], &*ptr);
        }
    }
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

void convert22DArray(unsigned char **ptr, uint64_t **&arr, uint32_t n_row, uint32_t n_col){
    // Recover a two dimensional dynamic array from a straightened two dimensional array
    arr = new uint64_t*[n_row];
    for( uint32_t i = 0; i < n_row; i++) {
        arr[i] = new uint64_t[n_col];
        for( uint32_t j = 0; j < n_col; j++) {
            arr[i][j] = convert2Long(ptr);
        }
    }
}

void convert23DArray(unsigned char **ptr, uint64_t ***&arr, uint32_t n_arrs, uint32_t n_row, uint32_t n_col){
    // Recover a two dimensional dynamic array from a straightened two dimensional array
    arr = new uint64_t**[n_arrs];
    for( uint32_t g = 0; g < n_arrs; g++) {
        arr[g] = new uint64_t*[n_row];
        for( uint32_t i = 0; i < n_row; i++) {
            arr[g][i] = new uint64_t[n_col];
            for( uint32_t j = 0; j < n_col; j++) {
                arr[g][i][j] = convert2Long(ptr);
            }
        }
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

double **convert2double(uint64_t **x, uint32_t n_row, uint32_t n_col, int precision=FRAC) {
    /*
     * Convert two-dimensional uint64 array to two-dimensional double array
     *
     * Input(s):
     *  - x: two-dimensional uint64 array
     *  - n_rows: number of rows in matrix
     *  - n_col: number of columns in matrix
     *  - precision: precision of the number format of uint64 array
     *
     * Output(s):
     *  - res: two-dimensional double array
     */
    double **res = new double*[n_row];
    double tmp = 1 << precision;
    for (uint32_t i = 0; i < n_row; i++) {
        res[i] = new double[n_col];
        for(uint32_t j = 0; j < n_col; j++) {
            if ((int) (x[i][j] >> 63) == 1) { // negative value
                res[i][j] = -1 * ((double) (~x[i][j] + 1) / tmp);
            } else { // positive value
                res[i][j] = ((double) x[i][j] / tmp);
            }
        }
    }
    return res;
}

uint64_t **convert2uint64(double** x, uint32_t n_row, uint32_t n_col, int precision=FRAC) {
    /*
     * Convert two-dimensional double array to two-dimensional uint64 array
     *
     * Input(s):
     *  - x: two-dimensional double array
     *  - n_rows: number of rows in matrix
     *  - n_col: number of columns in matrix
     *  - precision: precision of the number format of desired uint64 array
     *
     * Output(s):
     *  - res: two-dimensional uint64 array
     */
    uint64_t **res = new uint64_t*[n_row];
    double tmp = 1 << precision;
    for (uint32_t i = 0; i < n_row; i++) {
        res[i] = new uint64_t[n_col];
        for( uint32_t j = 0; j < n_col; j++) {
            if (x[i][j] < 0) { // negative values
                res[i][j] = (uint64_t) 0 - (uint64_t) floor(abs(x[i][j] * (1 << precision)));
            } else { // positive values
                res[i][j] = (uint64_t) floor(x[i][j] * (1 << precision));
            }
        }
    }
    return res;
}

uint64_t *straighten2DArray(uint64_t** x, uint32_t n_row, uint32_t n_col) {
    /*
     * Straighten two-dimensional uint64 array
     *
     * Input(s):
     *  - x: two-dimensional uint64 array
     *  - n_rows: number of rows in matrix
     *  - n_col: number of columns in matrix
     *
     * Output(s):
     *  - str_x: one-dimensional uint64 vector
     */
    uint64_t* str_x = new uint64_t[n_row * n_col];
    for(uint32_t i = 0; i < n_row; i++) {
        for( uint32_t j = 0; j < n_col; j++) {
            str_x[i * n_col + j] = x[i][j];
        }
    }
    return str_x;
}

double *straighten2DArray(double** x, uint32_t n_row, uint32_t n_col) {
    /*
     * Straighten two-dimensional double array
     *
     * Input(s):
     *  - x: two-dimensional double array
     *  - n_rows: number of rows in matrix
     *  - n_col: number of columns in matrix
     *
     * Output(s):
     *  - str_x: one-dimensional double vector
     */
    double* str_x = new double[n_row * n_col];
    for(uint32_t i = 0; i < n_row; i++) {
        for( uint32_t j = 0; j < n_col; j++) {
            str_x[i * n_col + j] = x[i][j];
        }
    }
    return str_x;
}

uint64_t getModularInverse(uint64_t a){
    /**
     * Get the MoDular Inverse (MDI) of a given number a with specified modulo. For the resulting/returned value b must hold
     *      ab mod(modulo) are congruent to 1.
     * @param a the value for which the modular inverse shall be calculated.
     * The modulo under which a and the inverse are multiplied equal to 1 will always be the ring size.
     * @return the modular inverse of a under the ring size of 16.
     */
    uint64_t r = a;
    for (int i = 0; i < 6; i++) {// (n = 6) for (2^64)
        r = r * (2 - r * a); // ignore overflow.
    }
    return r;
}





// Local functions which does not require security and works with secret shared values
uint64_t local_MUL(uint64_t a, uint64_t b) {
    /*
     * Input(s)
     * a: the first multiplicand in our number format- uint64_t
     * b: the second multiplicand in our number format - uint64_t
     *
     * Output(s)
     * Returns the multiplication of a and b - uint64_t
     */
    uint64_t z = a * b;
    // restore the fractional part - refer to SecureNN for more details
    if ((z >> 63) == 0) {
        z = z >> FRAC;
    } else {
        z = -1 * ((-1 * z) >> FRAC);
    }
    return z;
}

uint64_t* local_MUL(uint64_t *a, uint64_t *b, uint32_t size) {
    /*
     * Input(s)
     * a: one of the vectors of the multiplicands - uint64_t vector
     * b: the other vector of the multiplicands - uint64_t vector
     * size: the size of the vectors a and b
     *
     * Output(s)
     * Returns an uint64_t vector containing the result of the multiplication
     */
    uint64_t *result = new uint64_t[size];
    for (uint32_t i = 0; i < size; i++) {
        result[i] = local_MUL(a[i], b[i]);
    }
    return result;
}

uint64_t** local_MATMATMUL(uint64_t **a, uint64_t **b, uint32_t a_row, uint32_t a_col, uint32_t b_col) {
    /*
     * Perform multiplication of matrices a and b. The function assumes that the number of columns of a equals to
     * the number of rows of b.
     *
     * Input(s)
     * a: two dimensional matrix of size a_row-by-a_col
     * b: two dimensional matrix of size a_col-by-b_col
     *
     * Output(s)
     * Returns a matrix of size a_row-by-b_col
     */
    uint64_t **result = new uint64_t *[a_row];
    uint64_t tmp_sum = 0;
    for (uint32_t i = 0; i < a_row; i++) {
        result[i] = new uint64_t[b_col];
        for (uint32_t j = 0; j < b_col; j++) {
            tmp_sum = 0;
            for (uint32_t k = 0; k < a_col; k++) {
                tmp_sum += local_MUL(a[i][k], b[k][j]);
            }
            result[i][j] = tmp_sum;
        }
    }
    return result;
}

uint64_t*** local_MATMATMUL(uint64_t ***a, uint64_t ***b, uint32_t n_mats, uint32_t a_row, uint32_t a_col, uint32_t b_col) {
    /*
     * Perform several multiplication of matrices a and b. The function assumes that the number of columns of a equals to
     * the number of rows of b.
     *
     * Input(s)
     * a: three dimensional matrix of size n_mats-by-a_row-by-a_col
     * b: three dimensional matrix of size n_mats-by-a_col-by-b_col
     *
     * Output(s)
     * Returns a matrix of size n_mats-by-a_row-by-b_col
     */
    uint64_t ***result = new uint64_t **[n_mats];
    for(uint32_t g = 0; g < n_mats; g++) {
        result[g] = new uint64_t*[a_row];
        uint64_t tmp_sum = 0;
        for (uint32_t i = 0; i < a_row; i++) {
            result[g][i] = new uint64_t[b_col];
            for (uint32_t j = 0; j < b_col; j++) {
                tmp_sum = 0;
                for (uint32_t k = 0; k < a_col; k++) {
                    tmp_sum += local_MUL(a[g][i][k], b[g][k][j]);
                }
                result[g][i][j] = tmp_sum;
            }
        }
    }
    return result;
}

#endif //PML_FLIB_H
