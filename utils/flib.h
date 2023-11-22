//
// Created by Mete Akgun on 04.07.20.
//

#ifndef PML_FLIB_H
#define PML_FLIB_H

#include <stdlib.h>
#include <stdint.h>
#include <thread>

uint64_t b_mask[9]= {0,0xff,0xffff,0xffffff,0xffffffff,0xffffffffff,0xffffffffffff,0xffffffffffffff,0xffffffffffffffff};

uint64_t ConvertToLong(unsigned char **ptr){
    uint64_t val = (*((uint64_t *)(*ptr)));
    (*ptr)+=8;
    return val;
}

uint64_t ConvertToLong(unsigned char **ptr, int bsz){
    uint64_t val = (*((uint64_t *)(*ptr)))&b_mask[bsz];
    (*ptr)+=bsz;
    return val;
}

uint32_t ConvertToInt(unsigned char **ptr){
    uint32_t val = 0;
    for (int i=24;i>=0;i-=8){
        val=val+((uint64_t)(**ptr)<<i);
        (*ptr)++;
    }
    return val;
}

uint8_t ConvertToUint8(unsigned char **ptr){
    uint8_t val = (uint8_t)(**ptr);
    (*ptr)++;
    return val;
}

void ConvertToArray(unsigned char **ptr, uint8_t arr[], int sz){
    for (int i=0;i<sz;i++) {
        arr[i] = (**ptr);
        (*ptr)++;
    }
}

void ConvertToArray(unsigned char **ptr, uint64_t arr[], int sz){
    for (int i=0;i<sz;i++) {
        arr[i] = ConvertToLong(ptr);
    }
}

void ConvertToArray(unsigned char **ptr, uint64_t *&arr, uint32_t size){
    // Recover a one dimensional dynamic array
    arr = new uint64_t[size];
    for( uint32_t i = 0; i < size; i++) {
        arr[i] = ConvertToLong(ptr);
    }
}

void AddBitToCharArray(uint8_t val, uint8_t **ptr, uint8_t *bit_index){
    (**ptr)=(**ptr)^(val<<(*bit_index));
    if ((*bit_index) == 0){
        (*bit_index) = 7;
        (*ptr)++;
        (**ptr)=0;
    }
    else
        (*bit_index)-=1;
}


//This is compatible with AddBitToCharArray
uint8_t ConvertToByte(uint8_t **ptr, uint8_t *bit_index){
    uint8_t val = ((**ptr)>>(*bit_index))&0x1;
    if ((*bit_index) == 0){
        (*bit_index) = 7;
        (*ptr)++;
    }
    else
        (*bit_index)-=1;
    return val;
}

void AddValueToCharArray(uint64_t val,unsigned char **ptr, int bsz=8){
    *((uint64_t *)(*ptr)) = val;
    (*ptr)+=bsz;
}

void AddValueToCharArray(uint32_t val, unsigned char **ptr){
    for (int i=24;i>=0;i-=8){
        (**ptr)=(val>>i)&0xff;
        (*ptr)++;
    }
}
void AddValueToCharArray(uint8_t val, unsigned char **ptr){
    (**ptr)=(val)&0xff;
    (*ptr)++;
}

void AddValueToCharArray(uint8_t val[], unsigned char **ptr, int sz){
    for (int i=0;i<sz;i++){
        (**ptr)=(val[i])&0xff;
        (*ptr)++;
    }
}
void AddValueToCharArray(uint64_t *val, unsigned char **ptr, int sz){
    for (int i=0;i<sz;i++){
        AddValueToCharArray(val[i], ptr);
    }

}
void AddArrayToCharArray(uint64_t **val, unsigned char **ptr, uint32_t n_row, uint32_t n_col){
    // Add uint64_t vals in **val to the buffer to send
    for( uint32_t i = 0; i < n_row; i++) {
        for( uint32_t j = 0; j < n_col; j++) {
            AddValueToCharArray(val[i][j], &*ptr);
        }
    }
}
uint8_t Bit(uint64_t val, uint8_t ind){
    return (val>>ind)&0x1;
}

// what does this do exactly? TODO give more descriptive name
uint8_t Mod(int k, int n) {
    return ((k %= n) < 0) ? k+n : k;
}

uint64_t MersenneMod(uint64_t k, uint64_t p, uint8_t s) {
    uint64_t i = (k & p) + (k >> s);
    return (i >= p) ? i - p : i;
}

double ConvertToDouble(uint64_t x, int precision= FRACTIONAL_BITS) {
    double tmp = (double)((uint64_t) 1 << precision);
    if ((int) (x >> 63) == 1) {
        return -1 * ((double) (~x + 1) / tmp);
    } else {
        return ((double) x / tmp);
    }
}

uint64_t ConvertToUint64(double x, int precision = FRACTIONAL_BITS) {
    if (x < 0) {
        return (uint64_t) 0 - (uint64_t) floor(abs(x * (((uint64_t) 1) << precision)));
    } else {
        return (uint64_t) floor(x * (((uint64_t) 1) << precision));
    }
}

void ConvertTo2dArray(unsigned char **ptr, uint64_t **&arr, uint32_t n_row, uint32_t n_col){
    // Recover a two-dimensional dynamic array from a straightened two-dimensional array
    arr = new uint64_t*[n_row];
    for( uint32_t i = 0; i < n_row; i++) {
        arr[i] = new uint64_t[n_col];
        for( uint32_t j = 0; j < n_col; j++) {
            arr[i][j] = ConvertToLong(ptr);
        }
    }
}

void ConvertTo3dArray(unsigned char **ptr, uint64_t ***&arr, uint32_t n_arrs, uint32_t n_row, uint32_t n_col){
    // Recover a three-dimensional dynamic array from a straightened three-dimensional array
    arr = new uint64_t**[n_arrs];
    for( uint32_t g = 0; g < n_arrs; g++) {
        arr[g] = new uint64_t*[n_row];
        for( uint32_t i = 0; i < n_row; i++) {
            arr[g][i] = new uint64_t[n_col];
            for( uint32_t j = 0; j < n_col; j++) {
                arr[g][i][j] = ConvertToLong(ptr);
            }
        }
    }
}

double *ConvertToDouble(uint64_t *x, uint32_t size, int precision= FRACTIONAL_BITS) {
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

uint64_t* ConvertToUint64(double* x, uint32_t size, int precision= FRACTIONAL_BITS) {
    uint64_t *res = new uint64_t[size];
    double tmp = 1 << precision;
    for (int i = 0; i < size; i++) {
        if (x[i] < 0) {
            res[i] = (uint64_t) 0 - (uint64_t) floor(abs(x[i] * tmp));
        } else {
            res[i] = (uint64_t) floor(x[i] * tmp);
        }
    }
    return res;
}

double **ConvertToDouble(uint64_t **x, uint32_t n_row, uint32_t n_col, int precision= FRACTIONAL_BITS) {
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

double ***ConvertToDouble(
    uint64_t ***values, size_t matrix_count, size_t row_count, size_t column_count, int precision = FRACTIONAL_BITS
) {
    double*** result = new double**[matrix_count];
    for (int i = 0; i < matrix_count; i++) {
        result[i] = ConvertToDouble(values[i], row_count, column_count, precision);
    }
    return result;
}

uint64_t **ConvertToUint64(double** x, uint32_t n_row, uint32_t n_col, int precision= FRACTIONAL_BITS) {
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
                res[i][j] = (uint64_t) 0 - (uint64_t) floor(abs(x[i][j] * tmp));
            } else { // positive values
                res[i][j] = (uint64_t) floor(x[i][j] * tmp);
            }
        }
    }
    return res;
}

uint64_t *Straighten2dArray(uint64_t** x, uint32_t n_row, uint32_t n_col) {
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

double *Straighten2dArray(double** x, uint32_t n_row, uint32_t n_col) {
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
long long GetModularInverseN(long long a, long long m) {
    long long m0 = m;
    long long y = 0, x = 1;

    if (m == 1)
        return 0;

    while (a > 1) {
        long long q = a / m;
        long long t = m;
        m = a % m, a = t;
        t = y;
        y = x - q * y;
        x = t;

    }

    if (x < 0)
        x += m0;
    return x;
}

/** Get modular multiplication of given numbers
 * @param a, b are multiplicants
 * @param m modulus
 * @return a*b%n
 *
 * Modulo operation is done by eMod Method
 * 3rd parameter of MersenneMod should be log(m). In this case it is 61 since 2^61-1
 * */
long long MultMod(long long x, long long y, long long m) {
    long long res = 0;
    x = MersenneMod(x, m, 61);//x % m;
    while (y > 0) {
        if (y % 2 == 1)
            res = (res + x) % m;
        x = (x * 2) % m;
        y /= 2;
    }
    return res % m;
}
uint64_t GetModularInverse(uint64_t a){
    /**
     * Get the Modular Inverse (ModularInverse) of a given number a with specified modulo. For the resulting/returned value b must hold
     *      ab Mod(modulo) are congruent to 1.
     * @param a the value for which the modular inverse shall be calculated.
     * The modulo under which a and the inverse are multiplied equal to 1 will always be the ring size.
     * @return the modular inverse of a under the ring size of 16.
     */
    uint64_t r = a;
    for (int i = 0; i < 6; i++) {// (n = 6) because 2^6 is 64
        r = r * (2 - r * a); // ignore overflow.
    }
    return r;
}


/*
 * Arithmetic shift defined in SecureNN: we fill the significant bits with the most significant bit.
 */
uint64_t ArithmeticShift(uint64_t z, int n_shift = FRACTIONAL_BITS) {
    z = static_cast<uint64_t>( static_cast<int64_t>(z) >> n_shift);
    return z;
}


// Local functions which does not require security and works with secret shared values
uint64_t LocalMultiply(uint64_t a, uint64_t b, int shift = FRACTIONAL_BITS) {
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
    // v1
    if ((z >> 63) == 0) {
        z = z >> shift;
    } else {
        z = -1 * ((-1 * z) >> shift);
    }
    return z;
}

uint64_t* LocalMultiply(uint64_t *a, uint64_t *b, uint32_t size) {
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
        result[i] = LocalMultiply(a[i], b[i]);
    }
    return result;
}

uint64_t** LocalMatrixMatrixMultiply(uint64_t **a, uint64_t **b, uint32_t a_row, uint32_t a_col, uint32_t b_col, int shift = FRACTIONAL_BITS) {
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
    for (uint32_t i = 0; i < a_row; i++) {
        result[i] = new uint64_t[b_col];
        for(uint32_t j = 0; j < b_col; j++) {
            result[i][j] = 0;
        }
        for (uint32_t j = 0; j < a_col; j++) {
            for (uint32_t k = 0; k < b_col; k++) {
                result[i][k] += LocalMultiply(a[i][j], b[j][k], shift);
            }
        }
    }


    return result;
}

uint64_t*** LocalMatrixMatrixMultiply(
        const uint64_t *const *const *const a,
        const uint64_t *const *const *const b,
        uint32_t n_mats,
        uint32_t a_row,
        uint32_t a_col,
        uint32_t b_col,
        int shift = FRACTIONAL_BITS
        ) {
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
        for (uint32_t i = 0; i < a_row; i++) {
            result[g][i] = new uint64_t[b_col];
            for(uint32_t j = 0; j < b_col; j++) {
                result[g][i][j] = 0;
            }
            for (uint32_t j = 0; j < a_col; j++) {
                for (uint32_t k = 0; k < b_col; k++) {
                    result[g][i][k] += LocalMultiply(a[g][i][j], b[g][j][k], shift);
                }
            }
        }
    }


    return result;
}

void Delete2dMatrix(uint64_t** matrix, size_t row_count) {
    for (int i = 0; i < row_count; i++) {
        delete[] matrix[i];
    }
    delete[] matrix;
}

void Delete2dMatrix(double** matrix, size_t row_count) {
    for (int i = 0; i < row_count; i++) {
        delete[] matrix[i];
    }
    delete[] matrix;
}

void Delete3dMatrix(uint64_t*** matrix, size_t matrix_count, size_t row_count) {
    for (int i = 0; i < matrix_count; i++) {
        for (int ii = 0; ii < row_count; ii++) {
            delete[] matrix[i][ii];
        }
        delete[] matrix[i];
    }
    delete[] matrix;
}

void Delete3dMatrix(double*** matrix, size_t matrix_count, size_t row_count) {
    for (int i = 0; i < matrix_count; i++) {
        for (int ii = 0; ii < row_count; ii++) {
            delete[] matrix[i][ii];
        }
        delete[] matrix[i];
    }
    delete[] matrix;
}

void AddToBuffer(const uint64_t *const val, unsigned char *ptr, int size, uint32_t bsz= 8){
    unsigned char *ptr_tmp = ptr;
    for (int i = 0; i<size-2; i++) {
        AddValueToCharArray(*(val+i), &ptr_tmp, bsz);
    }
    for (int i = 0; i<2*bsz; i++) {
        *(ptr_tmp+i) = 0;
    }
    *((uint64_t *)(ptr_tmp)) += *(val+size-2);
    ptr_tmp+=bsz;
    *((uint64_t *)(ptr_tmp)) += *(val+size-1);
}

void WriteToBuffer(const uint64_t *const val, unsigned char *ptr, int size, uint32_t bsz= 8){
    thread thr[SOCKET_NUMBER];
    int block_size = (int)ceil(size*1.0/SOCKET_NUMBER);
    if (block_size<50)
        block_size = size;
    uint32_t written = 0;
    uint32_t thr_num = 0;
    for (int i=0;i<SOCKET_NUMBER;i++){
        if ((size-written)<=block_size){
            block_size = size-written;
            thr[i] = thread(AddToBuffer, val + written, ptr + (written * bsz), block_size, bsz);
            thr_num++;
            break;
        }else{
            thr[i] = thread(AddToBuffer, val + written, ptr + (written * bsz), block_size, bsz);
            thr_num++;
            written+=block_size;
        }
    }
    for (int i=0;i<thr_num;i++){
        thr[i].join();
    }
}

void ReadBufferSingular(uint64_t *val, unsigned char *ptr, int size, uint32_t bsz= 8){
    unsigned char *ptr_tmp = ptr;
    for (int i = 0; i < size; i++) {
        *(val+i) = ConvertToLong(&ptr_tmp, bsz);
    }
}

void ReadBuffer(uint64_t *val, unsigned char *ptr, int size, uint32_t bsz= 8){
    thread thr[SOCKET_NUMBER];
    int block_size = (int)ceil(size*1.0/SOCKET_NUMBER);
    uint32_t read = 0;
    uint32_t thr_num = 0;
    for (int i=0;i<SOCKET_NUMBER;i++){
        if ((size-read)<=block_size){
            block_size = size-read;
            thr[i] = thread(ReadBufferSingular, val + read, ptr + (read * bsz), block_size, bsz);
            thr_num++;
            break;
        }else{
            thr[i] = thread(ReadBufferSingular, val + read, ptr + (read * bsz), block_size, bsz);
            thr_num++;
            read+=block_size;
        }
    }
    for (int i=0;i<thr_num;i++){
        thr[i].join();
    }
}

#endif //PML_FLIB_H
