//
// Created by Debora Jutz on 03.02.22.
//

#ifndef CECILIA_TEST_FUNCTIONS_H
#define CECILIA_TEST_FUNCTIONS_H

#include "../core/core.h"
#include "flib.h"
#include <fstream>
#include <tuple>
#include <random>

#endif //CECILIA_TEST_FUNCTIONS_H

// Random matrix generation functions
/**
 *
 * @param proxy
 * @param size size of the vector to be created
 * @param max_num highest random number to be contained in the returned data, default is 10
 * @param neg_flag allow negative values as default, set to false if they shall not be allowed
 * @return
 */
static double* random_1D_data(Party *proxy, int size, double max_num=10, bool neg_flag=true) {
    double* mat_data = new double[size];
    for (int i = 0; i < size; i++) {
        mat_data[i] = max_num * (static_cast <float> (rand()) / static_cast <float> (RAND_MAX));
        if (neg_flag == true && rand() % 2 == 0) {
            mat_data[i] *= -1;
        }
    }
    return mat_data;
}

static double* random_1D_data(Party *proxy, int size, double min_num, double max_num) {
    double* mat_data = new double[size];
    random_device rd; // obtain a random number from hardware
    mt19937 gen(rd()); // seed the generator
    uniform_real_distribution<> distr(min_num, max_num); // define the range
    for (int i = 0; i < size; i++) {
        mat_data[i] = distr(gen);
    }
    return mat_data;
}

static double** random_2D_data(Party *proxy, int n_row, int n_col, double min_num, double max_num) {
    double d_tmp1;
    uint64_t tmp1;
    double** mat_data = new double *[n_row];
    random_device rd; // obtain a random number from hardware
    mt19937 gen(rd()); // seed the generator
    uniform_real_distribution<> distr(min_num, max_num); // define the range
    for (int i = 0; i < n_row; i++) {
        mat_data[i] = new double[n_col];
        for(int j = 0; j < n_col; j++) {
            mat_data[i][j] = distr(gen);
        }
    }
    return mat_data;
}

static uint64_t** random_gram_matrix(Party *proxy, int n_row, int n_col) {
    double d_tmp1;
    uint64_t tmp1, a;
    int role = proxy->getPRole();

    // data matrix generation
    double **rand_matrix = new double *[n_row];
    for(int i = 0; i < n_row; i++) {
        rand_matrix[i] = new double[n_col];

        double tmp_sum = 0;
        for(int j = 0; j < n_col; j++) {
            d_tmp1 = 10 * (static_cast <float> (rand()) / static_cast <float> (RAND_MAX));
            tmp_sum += pow(d_tmp1, 2);
//            if (rand() % 2 == 0) {
//                d_tmp1 *= -1;
//            }
            rand_matrix[i][j] = d_tmp1;
        }

//        tmp_sum = sqrt(tmp_sum);
//        for( int j = 0; j < n_col; j++) {
//            rand_matrix[i][j] /= tmp_sum;
//        }
    }

    // gram matrix computation
    cout << "\nGram matrix: " << endl;
    double** gram_matrix = new double*[n_row];
    for(int i = 0; i < n_row; i++) {
        gram_matrix[i] = new double[n_row];
        for(int j = 0; j < n_row; j++) {
            double tmp_sum = 0;
            for(int k = 0; k < n_col; k++) {
                tmp_sum += rand_matrix[i][k] * rand_matrix[j][k];
            }
            gram_matrix[i][j] = tmp_sum;
            cout << gram_matrix[i][j] << "\t";
        }
        cout << endl;
    }
    cout << endl;


    // share generation
    uint64_t** invsqrt_data = new uint64_t*[n_row];
    for (int i = 0; i < n_row; i++) {
        invsqrt_data[i] = new uint64_t[n_row];
        for(int j = 0; j < n_row; j++) {
            // generate shares
            uint64_t val = 0;
            for (int k = 3; k >= 0; k -= 1) {
                a = rand() & 0xffff;
                val = val ^ (a << (k * 16));
            }

            if (role == 0) {
                invsqrt_data[i][j] = val;
            } else {
                invsqrt_data[i][j] = convert2uint64(gram_matrix[i][j]) - val;
            }
        }
    }

    return invsqrt_data;
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
        w_dim = proxy->generateCommonRandom() % (L_BIT / 4);  // divide by 2 so that a valid mColSize can be found
    }
    //TODO adapt to asymmetric window sizes
    /**
    while (w_rows < 1){
        w_rows = proxy->generateCommonRandom() % (L/4);
    }*/
    while (mColSize <= w_dim || (mColSize % w_dim) != 0){ // matrix contains at least 2 windows, windows fit completely
        mColSize = proxy->generateCommonRandom() % L_BIT;
    }
    while (mRowSize <= w_dim || (mRowSize % w_dim) != 0){
        mRowSize = proxy->generateCommonRandom() % L_BIT;
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


/*
 * Printing functions for debugging purposes. There are functions to print out secret shared scalar, 1D and
 * 2D arrays either in secret shared form or after reconstruction. There are similar functions for plaintext data as well.
 * Moreoever, for 1D and 2D arrays, there is an option that one can set to print out horizontally or vertically.
 */
// Printing the data of type uint64_t for debugging purposes
void print2DArray(string const &str1, uint64_t** x, uint32_t n_row, uint32_t n_col, bool horizontal=true) {
    // horizontal: true if the resulting print-out is desired row-by-col
    cout << "======================= " << str1 << " =======================" << endl;
    if(horizontal) {
        for(uint32_t i = 0; i < n_row; i++) {
            for(uint32_t j = 0; j < n_col; j++) {
                cout << x[i][j] << "\t";
            }
            cout << endl;
        }
    }
    else {
        for(uint32_t i = 0; i < n_col; i++) {
            cout << i << ".\t";
            for(uint32_t j = 0; j < n_row; j++) {
                cout << x[j][i] << "\t";
            }
            cout << endl;
        }
    }
    cout << "==============================================================" << endl;
}

void print1DArray(string const &str1, uint64_t* x, uint32_t size, bool horizontal=true) {
    cout << "======================= " << str1 << " =======================" << endl;
    if(horizontal) {
        for (uint32_t i = 0; i < size; i++) {
            cout << x[i] << "\t";
        }
        cout << endl;
    }
    else {
        for(uint32_t i = 0; i < size; i++) {
            cout << i << ". " << x[i] << endl;
        }
    }
    cout << "==============================================================" << endl;
}

void printValue(string const &str1, uint64_t x) {
    cout << "======================= " << str1 << " =======================" << endl;
    cout << x << endl;
    cout << "==============================================================" << endl;
}

// Printing the data of type uint64_t for debugging purposes after reconstructing and converting
/*
void print2DArrayRecAndConv(Party* proxy, string const &str1, uint64_t** x, uint32_t n_row, uint32_t n_col,
                            bool horizontal=true) {
    // horizontal: true if the resulting print-out is desired row-by-col
    double **d_x = convert2double(REC(proxy, x, n_row, n_col), n_row, n_col);
    cout << "======================= " << str1 << " =======================" << endl;
    if (horizontal) {
        for (uint32_t i = 0; i < n_row; i++) {
            for (uint32_t j = 0; j < n_col; j++) {
                cout << d_x[i][j] << "\t";
            }
            cout << endl;
        }
    }
    else {
        for(uint32_t i = 0; i < n_col; i++) {
            cout << i << ".\t";
            for(uint32_t j = 0; j < n_row; j++) {
                cout << d_x[j][i] << "\t";
            }
            cout << endl;
        }
    }
    cout << "==============================================================" << endl;
}


void print1DArrayRecAndConv(Party *proxy, string const &str1, uint64_t* x, uint32_t size, bool horizontal=true) {
    double * d_x = convert2double(REC( proxy, x, size), size);
    cout << "======================= " << str1 << " =======================" << endl;
    if(horizontal) {
        for(uint32_t i = 0; i < size; i++) {
            cout << d_x[i] << "\t";
        }
        cout << endl;
    }
    else {
        for(uint32_t i = 0; i < size; i++) {
            cout << i << ". " << d_x[i] << endl;
        }
    }
    cout << "==============================================================" << endl;
}

void printValueRecAndConv(Party *proxy, string const &str1, uint64_t x) {
    double d_x = convert2double(REC(proxy, x));
    cout << "======================= " << str1 << " =======================" << endl;
    cout << d_x << endl;
    cout << "==============================================================" << endl;
}
*/

// Printing the data of type double for debugging purposes
void print2DArray(string const &str1, double** x, uint32_t n_row, uint32_t n_col, bool horizontal=true) {
    // horizontal: true if the resulting print-out is desired row-by-col
    cout << "======================= " << str1 << " =======================" << endl;
    if(horizontal) {
        for(uint32_t i = 0; i < n_row; i++) {
            for(uint32_t j = 0; j < n_col; j++) {
                cout << x[i][j] << "\t";
            }
            cout << endl;
        }
    }
    else {
        for(uint32_t i = 0; i < n_col; i++) {
            cout << i << ".\t";
            for(uint32_t j = 0; j < n_row; j++) {
                cout << x[j][i] << "\t";
            }
            cout << endl;
        }
    }
    cout << "==============================================================" << endl;
}

void print1DArray(string const &str1, double* x, uint32_t size, bool horizontal=true) {
    cout << "======================= " << str1 << " =======================" << endl;
    if(horizontal) {
        for(uint32_t i = 0; i < size; i++) {
            cout << x[i] << "\t";
        }
        cout << endl;
    }
    else {
        for(uint32_t i = 0; i < size; i++) {
            cout << i << ". " << x[i] << endl;
        }
    }
    cout << "==============================================================" << endl;
}

void printValue(string const &str1, double x) {
    cout << "======================= " << str1 << " =======================" << endl;
    cout << x << endl;
    cout << "==============================================================" << endl;
}

// Matrix operations
static double** multiply_matrices(double** m1, double** m2, int m1_row, int m1_col, int m2_col) {
    double tmp;
    double** res = new double*[m1_row];
    for( int i = 0; i < m1_row; i++) {
        res[i] = new double[m2_col];
        for(int j = 0; j < m2_col; j++) {
            tmp = 0;
            for(int k = 0; k < m1_col; k++) {
                tmp += m1[i][k] * m2[k][j];
            }
            res[i][j] = tmp;
        }
    }
    return res;
}

static double* multiply_matrice_vector(double** m1, double* m2, int m1_row, int m1_col) {
    double tmp;
    double* res = new double[m1_row];
    for( int i = 0; i < m1_row; i++) {
        tmp = 0;
        for(int k = 0; k < m1_col; k++) {
            tmp += m1[i][k] * m2[k];
        }
        res[i] = tmp;
    }
    return res;
}

double*** double_MATMATMUL(double ***a, double ***b, uint32_t n_mats, uint32_t a_row, uint32_t a_col, uint32_t b_col) {
    /*
     * Perform several multiplication of double matrices a and b. The function assumes that the number of columns of a equals to
     * the number of rows of b.
     *
     * Input(s)
     * a: three dimensional matrix of size n_mats-by-a_row-by-a_col
     * b: three dimensional matrix of size n_mats-by-a_col-by-b_col
     *
     * Output(s)
     * Returns a matrix of size n_mats-by-a_row-by-b_col
     */
    if(DEBUG_FLAG == 1)
        cout << "************************************************************\ndouble_MATMATMUL is called" << endl;
    double ***result = new double **[n_mats];
    for(uint32_t g = 0; g < n_mats; g++) {
        result[g] = new double*[a_row];
        double tmp_sum = 0;
        for (uint32_t i = 0; i < a_row; i++) {
            result[g][i] = new double[b_col];
            for (uint32_t j = 0; j < b_col; j++) {
                tmp_sum = 0;
                for (uint32_t k = 0; k < a_col; k++) {
                    tmp_sum += a[g][i][k] * b[g][k][j];
                }
                result[g][i][j] = tmp_sum;
            }
        }
    }
    if(DEBUG_FLAG == 1)
        cout << "Returning from double_MATMATMUL...\n************************************************************" << endl;
    return result;
}

static void bubbleSort(double x[], int n) {
    bool exchanges;
    do {
        exchanges = false;  // assume no exchanges
        for (int i=0; i<n-1; i++) {
            if (x[i] > x[i+1]) {
                double temp = x[i]; x[i] = x[i+1]; x[i+1] = temp;
                exchanges = true;  // after exchange, must look again
            }
        }
    } while (exchanges);
}