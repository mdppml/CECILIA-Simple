//
// Created by Debora Jutz on 03.02.22.
//

#ifndef CECILIA_TEST_FUNCTIONS_H
#define CECILIA_TEST_FUNCTIONS_H

#include "../core/Party.h"
#include "flib.h"
#include <fstream>
#include <tuple>

#endif //CECILIA_TEST_FUNCTIONS_H

/**
 *
 * @param proxy
 * @param size size of the vector to be created
 * @param max_num highest random number to be contained in the returned data, default is 10
 * @param neg_flag allow negative values as default, set to false if they shall not be allowed
 * @return
 */
static double* random_1D_data(Party *proxy, int size, int max_num=10, bool neg_flag=true) {
    double* mat_data = new double[size];
    for (int i = 0; i < size; i++) {
        mat_data[i] = max_num * (static_cast <float> (rand()) / static_cast <float> (RAND_MAX));
        if (neg_flag == true && rand() % 2 == 0) {
            mat_data[i] *= -1;
        }
    }
    return mat_data;
}

/**
 * Generate a random matrix whose values are stored in a simple vector. Matrix values as well as size are random but
 * some constraints exist: (1) matrix size will be high enough so that at least 4 windows fit in total;
 *                         (2) window size and matrix size match, meaning that
 *                         matrix row size = x * window row size and matrix column size = y * window column size
 *                         (3) row and column size will be greater 1.
 * @param proxy
 * @param matrix_number specifies the number of matrices of the random shape that shall be generated.
 *                      default is 2, so the matrix and a tmp vector are provided to generated secret shares.
 * @param max_num maximum value to be allowed in the matrix
 * @param neg_flag specifies if negative values are allowed.
 * @return a tuple containing two pointers:
 *          pointer of type uint64_t -> pointing to the generated matrix data
 *          pointer of type uint32_t -> pointing to four values in the given order:
 *                  matrix column size
 *                  matrix row size
 *                  window column size
 *                  window row size
 *         (matrix column size * matrix row size) will be the dimension of the matrix (m_size) and the generated matrix data will contain (m_size * matrix_number) values
 */
static double* random_window_matrix(Party *proxy, uint8_t matrix_number=2, double max_num=255.99, bool neg_flag=true){
    cout << "Generate matrix with random values... " << endl;
    uint32_t mColSize = 0, mRowSize = 0;
    uint32_t w_dim = 0; // w_rows = 0;
    while (w_dim < 1){
        w_dim = proxy->generateCommonRandom() % (L / 4);  // divide by 2 so that a valid mColSize can be found
    }
    //TODO adapt to asymmetric window sizes
    /**
    while (w_rows < 1){
        w_rows = proxy->generateCommonRandom() % (L/4);
    }*/
    while (mColSize <= w_dim || (mColSize % w_dim) != 0){ // matrix contains at least 2 windows, windows fit completely
        mColSize = proxy->generateCommonRandom() % L;
    }
    while (mRowSize <= w_dim || (mRowSize % w_dim) != 0){
        mRowSize = proxy->generateCommonRandom() % L;
    }

    //cout << "mR=" << mRowSize << "; mC=" << mColSize << "; window=" << w_dim << " x " << w_dim << endl;
    uint64_t mSize = mRowSize * mColSize;
    //                            for the 4 random parts: matrix A and B, mTmp A and B
    double *randData = new double [mSize*matrix_number + 4]; // +4 for the 4 values specifying matrix + window dimensions
    randData = random_1D_data(proxy, mSize*matrix_number, max_num, neg_flag);
    randData[-4] = mColSize;
    randData[-3] = mRowSize;
    randData[-2] = w_dim;
    randData[-1] = 0; //w_rows;
    return randData;
}
