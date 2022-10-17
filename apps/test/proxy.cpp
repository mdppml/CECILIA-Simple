#include <cstdlib>
#include <iostream>
#include <deque>
#include <chrono>
#include <tuple>
#include <iomanip>
#include "../../core/core.h"
#include "../../core/cnn.h"
#include "../../core/rkn.h"
#include "../../utils/flib.h"
#include "../cnn/model_parser.h"
#include "../cnn/mnist_data/mnist_reader.hpp"
#include <bitset>


using namespace std;
//using namespace Eigen;

constexpr int MIN_VAL = -150;
constexpr int MAX_VAL = 1400;
constexpr int sz = 1000;
constexpr int WSZ = 10;

bool MUL_Test(Party *proxy) {
    cout << setfill('*') << setw(50) << "Calling MUL";
    cout << setfill('*') << setw(49) << "*" << endl;

    double xd = MIN_VAL + (double)(proxy->generateCommonRandom() & RAND_MAX) / ((double)(RAND_MAX / (MAX_VAL - MIN_VAL)));
    double yd = MIN_VAL + (double)(proxy->generateCommonRandom() & RAND_MAX) / ((double)(RAND_MAX / (MAX_VAL - MIN_VAL)));
    uint64_t x = proxy->createShare(xd);
    uint64_t y = proxy->createShare(yd);
    proxy->SendBytes(CORE_MUL);
    uint64_t r = MUL(proxy, x, y);
    // checking the result
    double rd = convert2double(REC(proxy, r));
    double rcd = (xd * yd);
    cout << "diff: " << rd - rcd << endl;
    if ((int) (rd - rcd) == 0)
        cout << "MUL works correctly" << endl;
    else {
        cout << "MUL works incorrectly" << endl;
    }
    return ((int) (rd - rcd) == 0);
}

bool MMUL_Test(Party *proxy) {
    cout << setfill('*') << setw(50) << "Calling MMUL";
    cout << setfill('*') << setw(49) << "*" << endl;
    uint64_t x[sz], y[sz], z[sz];
    uint32_t *params = new uint32_t[1];
    params[0] = sz;
    for (int i = 0; i < sz; i++) {
        double xd = MIN_VAL +
                    (double) (proxy->generateCommonRandom() & RAND_MAX) / ((double) (RAND_MAX / (MAX_VAL - MIN_VAL)));
        double yd = MIN_VAL +
                    (double) (proxy->generateCommonRandom() & RAND_MAX) / ((double) (RAND_MAX / (MAX_VAL - MIN_VAL)));
        x[i] = proxy->createShare(xd);
        y[i] = proxy->createShare(yd);
    }
    proxy->SendBytes(CORE_MMUL, params, 1);
    uint64_t *r = MUL(proxy, x, y, sz);
    // checking the result
    bool flag = true;
    for (int i = 0; i < sz; i++) {
        double xd = convert2double(REC(proxy, x[i]));
        double yd = convert2double(REC(proxy, y[i]));
        double rd = convert2double(REC(proxy, r[i]));
        double rcd = (xd * yd);
        if ((int) (rd - rcd) != 0) {
            flag = false;
            cout << "Absolute difference of multiplication of " << xd << " and " << yd << ": " << abs(rd - rcd) << endl;
            break;
        }
    }
    if (flag)
        cout << "MMUL works correctly" << endl;
    else
        cout << "MMUL works incorrectly" << endl;
    return flag;
}

void MOC_Test(Party *proxy) {
    cout << setfill('*') << setw(50) << "Calling MOC";
    cout << setfill('*') << setw(49) << "*" << endl;
    uint64_t x = proxy->generateRandom() & N1_MASK;
    proxy->SendBytes(CORE_MC);
    uint64_t r = MOC(proxy, x);
    // checking the result
    uint64_t x_reconstructed = REC(proxy, x, N1_MASK);
    uint64_t r_reconstructed = REC(proxy, r);
    if (x_reconstructed == r_reconstructed)
        cout << "MOC works correctly" << endl;
    else
        cout << "MOC works incorrectly" << endl;
}

void MMOC_Test(Party *proxy) {
    cout << setfill('*') << setw(50) << "Calling Vectorized MOC";
    cout << setfill('*') << setw(49) << "*" << endl;
    uint64_t x[sz];
    uint32_t *params = new uint32_t[1];
    params[0] = sz;

    for (int i = 0; i < sz; i++)
        x[i] = proxy->generateRandom() & N1_MASK;
    proxy->SendBytes(CORE_MMC, params, 1);
    uint64_t *r = MOC(proxy, x, sz);
    // checking the result
    uint64_t *x_reconstructed = REC(proxy, x, sz, N1_MASK);
    uint64_t *r_reconstructed = REC(proxy, r, sz);
    bool flag = true;
    for (int i = 0; i < sz; i++) {
        if (x_reconstructed[i] != r_reconstructed[i]) {
            flag = false;
            break;
        }
    }
    if (flag)
        cout << "Vectorized MOC works correctly" << endl;
    else
        cout << "Vectorized MOC works incorrectly" << endl;

}

void MSB_Test(Party *proxy) {
    cout << setfill('*') << setw(50) << "Calling MSB";
    cout << setfill('*') << setw(49) << "*" << endl;
    uint64_t x = 3 - 5; //proxy->generateRandom();
    proxy->SendBytes(CORE_MSB);
    uint64_t r = MSB(proxy, x);
    // checking the result
    uint64_t x_reconstructed = REC(proxy, x);
    uint64_t r_reconstructed = REC(proxy, r);
    uint64_t r_computed = (x_reconstructed >> (L_BIT - 1)) << FRAC;
    if (r_reconstructed == r_computed) {
        cout << "MSB works correctly" << endl;
    } else {
        cout << "MSB works incorrectly" << endl;
        cout << "computed MSB = " << bitset<64>(r_computed) << "reconstructed MSB = " << bitset<64>(r_reconstructed)
             << "for original X = " << bitset<64>(x_reconstructed) << endl;
    }
}

void MMSB_Test(Party *proxy) {
    cout << setfill('*') << setw(50) << "Calling Vectorized MSB";
    cout << setfill('*') << setw(49) << "*" << endl;
    uint32_t *params = new uint32_t[1];
    params[0] = sz;

    uint64_t x[sz];
    for (int i = 0; i < sz; i++) {
        x[i] = proxy->generateRandom();
    }
    proxy->SendBytes(CORE_MMSB, params, 1);
    uint64_t *r = MSB(proxy, x, sz);
    //uint64_t *r = MSBv2(proxy,x,sz);
    // checking the result
    uint64_t *x_reconstructed = REC(proxy, x, sz);
    uint64_t *r_reconstructed = REC(proxy, r, sz);
    bool flag = true;
    for (int i = 0; i < sz; i++) {
        uint64_t r_computed = (x_reconstructed[i] >> (L_BIT - 1)) << FRAC;
        if (r_computed != r_reconstructed[i]) {
            flag = false;
            break;
        }
    }
    if (flag) {
        cout << "Vectorized MSB works correctly" << endl;
    } else {
        cout << "Vectorized MSB works incorrectly" << endl;
    }
}

void CMP_Test(Party *proxy) {
    cout << setfill('*') << setw(50) << "Calling CMP";
    cout << setfill('*') << setw(49) << "*" << endl;
    double xd =
            MIN_VAL + (double) (proxy->generateCommonRandom() & RAND_MAX) / ((double) (RAND_MAX / (MAX_VAL - MIN_VAL)));
    double yd =
            MIN_VAL + (double) (proxy->generateCommonRandom() & RAND_MAX) / ((double) (RAND_MAX / (MAX_VAL - MIN_VAL)));
    uint64_t x = proxy->createShare(xd);
    uint64_t y = proxy->createShare(yd);
    proxy->SendBytes(CORE_CMP);
    uint64_t r = CMP(proxy, x, y);
    // checking the result
    xd = convert2double(REC(proxy, x));
    yd = convert2double(REC(proxy, y));
    uint64_t rd = REC(proxy, r);
    uint64_t r_computed = (xd >= yd) << FRAC;
    if (rd == r_computed)
        cout << "CMP works correctly" << endl;
    else
        cout << "CMP works incorrectly" << endl;
}

void MCMP_Test(Party *proxy) {
    cout << setfill('*') << setw(50) << "Calling Vectorized CMP";
    cout << setfill('*') << setw(49) << "*" << endl;
    uint64_t x[sz], y[sz];
    uint32_t *params = new uint32_t[1];
    params[0] = sz;

    for (int i = 0; i < sz; i++) {
        double xd = MIN_VAL +
                    (double) (proxy->generateCommonRandom() & RAND_MAX) / ((double) (RAND_MAX / (MAX_VAL - MIN_VAL)));
        double yd = MIN_VAL +
                    (double) (proxy->generateCommonRandom() & RAND_MAX) / ((double) (RAND_MAX / (MAX_VAL - MIN_VAL)));
        x[i] = proxy->createShare(xd);
        y[i] = proxy->createShare(yd);
    }
    proxy->SendBytes(CORE_MCMP, params, 1);
    uint64_t *r = CMP(proxy, x, y, sz);
    // checking the result
    uint64_t *x_reconstructed = REC(proxy, x, sz);
    uint64_t *y_reconstructed = REC(proxy, y, sz);
    uint64_t *r_reconstructed = REC(proxy, r, sz);
    bool flag = true;
    for (int i = 0; i < sz; i++) {
        uint64_t r_computed = (convert2double(x_reconstructed[i]) >= convert2double(y_reconstructed[i])) << FRAC;
        if (r_computed != r_reconstructed[i]) {
            flag = false;
            break;
        }
    }
    if (flag)
        cout << "Vectorized CMP works correctly" << endl;
    else
        cout << "Vectorized CMP works incorrectly" << endl;

}

void MUX_Test(Party *proxy) {
    cout << setfill('*') << setw(50) << "Calling MUX";
    cout << setfill('*') << setw(49) << "*" << endl;
    double xd =
            MIN_VAL + (double) (proxy->generateCommonRandom() & RAND_MAX) / ((double) (RAND_MAX / (MAX_VAL - MIN_VAL)));
    double yd =
            MIN_VAL + (double) (proxy->generateCommonRandom() & RAND_MAX) / ((double) (RAND_MAX / (MAX_VAL - MIN_VAL)));
    double zd = (double) (proxy->generateCommonRandom() & 0x1);
    uint64_t x = proxy->createShare(xd);
    uint64_t y = proxy->createShare(yd);
    uint64_t z = proxy->createShare(zd);

    proxy->SendBytes(CORE_MUX);
    uint64_t r = MUX(proxy, x, y, z);
    // checking the result
    uint64_t x_reconstructed = REC(proxy, x);
    uint64_t y_reconstructed = REC(proxy, y);
    uint64_t z_reconstructed = REC(proxy, z);
    uint64_t r_reconstructed = REC(proxy, r);
    uint64_t r_computed;
    if (z_reconstructed == 0)
        r_computed = x_reconstructed;
    else if (z_reconstructed == (1 << FRAC))
        r_computed = y_reconstructed;
    if (r_reconstructed == r_computed)
        cout << "MUX works correctly" << endl;
    else
        cout << "MUX works incorrectly" << endl;
}

void MMUX_Test(Party *proxy) {
    cout << setfill('*') << setw(50) << "Calling Vectorized MUX";
    cout << setfill('*') << setw(49) << "*" << endl;
    uint64_t x[sz], y[sz], z[sz];

    uint32_t *params = new uint32_t[1];
    params[0] = sz;

    for (int i = 0; i < sz; i++) {
        double xd = MIN_VAL +
                    (double) (proxy->generateCommonRandom() & RAND_MAX) / ((double) (RAND_MAX / (MAX_VAL - MIN_VAL)));
        double yd = MIN_VAL +
                    (double) (proxy->generateCommonRandom() & RAND_MAX) / ((double) (RAND_MAX / (MAX_VAL - MIN_VAL)));
        double zd = (double) (proxy->generateCommonRandom() & 0x1);
        x[i] = proxy->createShare(xd);
        y[i] = proxy->createShare(yd);
        z[i] = proxy->createShare(zd);
    }
    proxy->SendBytes(CORE_MMUX, params, 1);
    uint64_t *r = MUX(proxy, x, y, z, sz);
    // checking the result
    double *x_reconstructed = convert2double(REC(proxy, x, sz), sz);
    double *y_reconstructed = convert2double(REC(proxy, y, sz), sz);
    double *z_reconstructed = convert2double(REC(proxy, z, sz), sz);
    double *r_reconstructed = convert2double(REC(proxy, r, sz), sz);
    bool flag = true;
    double r_computed;
    for (int i = 0; i < sz; i++) {
        if (z_reconstructed[i] == 0) {
            r_computed = x_reconstructed[i];
        } else if (z_reconstructed[i] == 1) {
            r_computed = y_reconstructed[i];
        }
        if (r_reconstructed[i] != r_computed) {
            flag = false;
            break;
        }
    }
    if (flag)
        cout << "Vectorized MUX works correctly" << endl;
    else {
        cout << "Vectorized MUX works incorrectly" << endl;
    }
}

void MAX_Test(Party *proxy) {
    cout << setfill('*') << setw(50) << "Calling MAX";
    cout << setfill('*') << setw(49) << "*" << endl;

    uint32_t mRows = WSZ * WSZ;
    uint32_t mCols = WSZ * WSZ;
    uint64_t mSize = mCols * mRows;
    uint32_t *params = new uint32_t[1];
    params[0] = mSize;

    uint64_t *shareOfMatrix = proxy->createShare(random_1D_data(proxy, mSize), mSize);

    proxy->SendBytes(CNN_MAX, params, 1);
    uint64_t max = MAX(proxy, shareOfMatrix, mSize);

    // checking the result
    double computed_max = -1;
    for (uint32_t position = 0; position < mSize; position++) {
        double matrixVal = convert2double(REC(proxy, shareOfMatrix[position]));
        if (matrixVal > computed_max) {
            computed_max = matrixVal;
        }
    }

    double pp_result = convert2double(REC(proxy, max));
    if (computed_max == pp_result) {
        cout << "MAX works correctly" << endl;
    } else {
        cout << "MAX works incorrectly" << endl;
        cout << "computed: " << pp_result << " should be: " << computed_max << endl;
    }

}

void MMAX_Test(Party *proxy) {
    cout << setfill('*') << setw(50) << "Calling vectorized MAX";
    cout << setfill('*') << setw(49) << "*" << endl;
    int precision = 3;

    // INIT PARAMETER
    uint32_t mmaxParams[4];
    uint32_t n_matrices = 3;
    mmaxParams[0] = n_matrices * WSZ; // matrix Rows
    mmaxParams[1] = WSZ; // matrix Columns
    uint64_t mSize = mmaxParams[1] * mmaxParams[0];
    mmaxParams[2] = 2; // window rows
    mmaxParams[3] = 2; // window cols

    uint64_t *shareOfMatrix = proxy->createShare(random_1D_data(proxy, mSize), mSize);
    // PERFORMING MMAX
    proxy->SendBytes(CNN_MMAX, mmaxParams, 4);
    uint64_t *max = MAX(proxy, shareOfMatrix, mmaxParams[0], mmaxParams[1], mmaxParams[2], mmaxParams[3]);

    // TESTING
    uint32_t window_length = mmaxParams[2] * mmaxParams[3];
    uint32_t number_of_windows = mSize / window_length; // in total, for all of the n_matrices matrices
    uint64_t *reconstructed_max = REC(proxy, max, number_of_windows);
    bool flag = true;
    uint64_t resorted[mSize];
    RST(shareOfMatrix, mmaxParams[1], mmaxParams[0], mmaxParams[3], mmaxParams[2], resorted);
    double *d_matrix = convert2double(REC(proxy, resorted, mSize), mSize);
    double computed_max[number_of_windows];

    for (uint32_t win = 0; win < number_of_windows; win++) {
        computed_max[win] = d_matrix[window_length * win];
        for (uint32_t win_element = 1; win_element < window_length; win_element++) {
            double matrixVal = d_matrix[window_length * win + win_element];
            if (matrixVal > computed_max[win]) {
                computed_max[win] = matrixVal;
            }
        }
        if (abs(computed_max[win] - convert2double(reconstructed_max[win]) > 0.001)) {
            flag = false;
        }
    }

    if (flag) {
        cout << "Vectorized MAX works correctly" << endl;
    } else {
        cout << "Vectorized MAX works incorrectly" << endl;
        if (mmaxParams[0] * mmaxParams[1] < 200) { //otherwise too many values
            print1DMatrixByWindows("resorted Matrix: ", d_matrix, mmaxParams[0], mmaxParams[1], 1,
                                   mmaxParams[2] * mmaxParams[3]);
        }
        print1DMatrixByWindows("computed max values (test): ", computed_max, 1, number_of_windows, 1, 1);
        print1DMatrixByWindows("VS result from method: ", convert2double(reconstructed_max, number_of_windows), 1,
                               number_of_windows, 1, 1);
    }
}

void ARGMAX_Test(Party *proxy) {
    cout << setfill('*') << setw(50) << "Calling ARGMAX";
    cout << setfill('*') << setw(49) << "*" << endl;

    uint32_t mRows = WSZ / 2;
    uint32_t mCols = WSZ / 2;
    uint64_t mSize = mCols * mRows;
    uint32_t *params = new uint32_t[1];
    params[0] = mSize;

    uint64_t *shareOfMatrix = proxy->createShare(random_1D_data(proxy, mSize), mSize);

    proxy->SendBytes(CNN_ARGMAX, params, 1);
    uint64_t argmax = ARGMAX(proxy, shareOfMatrix, mSize);

    // checking the result
    double computed_max = -1;
    double computed_argmax = -1;
    for (uint32_t position = 0; position < mSize; position++) {
        double matrixVal = convert2double(REC(proxy, shareOfMatrix[position]));
        if (matrixVal > computed_max) {
            computed_max = matrixVal;
            computed_argmax = position;
        }
    }

    double pp_result = convert2double(REC(proxy, argmax));
    if (computed_argmax == pp_result) {
        cout << "ARGMAX works correctly" << endl;
    } else {
        cout << "ARGMAX works incorrectly" << endl;
        cout << "computed: " << pp_result << " should be: " << computed_argmax << endl;
    }
}


void RST_Test(Party *proxy) {
    cout << setfill('*') << setw(50) << "Calling RST";
    cout << setfill('*') << setw(49) << "*" << endl;

    uint32_t mRows = WSZ * 10;
    uint32_t mCols = WSZ * 10;
    uint64_t mSize = mCols * mRows;

    uint32_t wRows = WSZ;
    uint32_t wCols = WSZ;
    uint64_t *shareOfMatrix = proxy->createShare(random_1D_data(proxy, mSize, 1000.0), mSize);
    auto *resorted = new uint64_t[sz];
    //print1DMatrixByWindows("RST original", convert2double(REC(proxy, shareOfMatrix, mSize), mSize), mRows, mCols, wRows, wCols);
    RST(shareOfMatrix, mCols, mRows, wCols, wRows, resorted);

    uint64_t wElements = wRows * wCols;
    uint64_t numberOfWins = mSize / wElements;
    if (mSize < 100) {
        print1DMatrixByWindows("RST finished: resorted", convert2double(REC(proxy, resorted, mSize), mSize),
                               numberOfWins, wElements, 1, 1);
    }
}

void RELU_Test(Party *proxy) {
    cout << setfill('*') << setw(50) << "Calling RELU";
    cout << setfill('*') << setw(49) << "*" << endl;

    uint64_t x = proxy->createShare(convert2double(proxy->generateCommonRandom()));
    proxy->SendBytes(CNN_RELU);
    uint64_t relu = RELU(proxy, x);
    uint64_t reconstructed_relu = REC(proxy, relu);

    // checking the result
    double computed_relu = -1;
    double originalX = convert2double(REC(proxy, x));
    if (originalX >= 0) {
        computed_relu = originalX;
    } else {
        computed_relu = 0;
    }

    double pp_result = convert2double(reconstructed_relu);
    if (computed_relu == pp_result) {
        cout << "RELU works correctly" << endl;
    } else {
        cout << "RELU works incorrectly" << endl;
        cout << "computed: " << pp_result << " should be: " << computed_relu << endl;
    }

}

void MRELU_Test(Party *proxy) {
    cout << setfill('*') << setw(50) << "Calling MRELU";
    cout << setfill('*') << setw(49) << "*" << endl;

    uint64_t size = WSZ * WSZ;
    uint32_t *params = new uint32_t[1];
    params[0] = size;
    uint64_t *x = proxy->createShare(random_1D_data(proxy, size, 255.0), size);
    proxy->SendBytes(CNN_MRELU, params, 1);
    uint64_t *mrelu = RELU(proxy, x, size);
    uint64_t *rec_mrelu = REC(proxy, mrelu, size);

    // checking the result
    double *correct_relu = new double[size];
    double *originalX = convert2double(REC(proxy, x, size), size);
    double *pp_result = convert2double(rec_mrelu, size);
    bool foundIncorrect = false;
    for (uint64_t i = 0; i < size; i++) {
        if (originalX[i] >= 0) {
            correct_relu[i] = originalX[i];
        } else {
            correct_relu[i] = 0;
        }
        if (correct_relu[i] != pp_result[i]) {
            foundIncorrect = true;
        }
    }

    if (!foundIncorrect) {
        cout << "MRELU works correctly" << endl;
    } else {
        cout << "MRELU works incorrectly" << endl;
        print1DArray("Original Values:", originalX, size);
        print1DArray("Computed RELU:", pp_result, size);
        print1DArray("VS Correct RELU:", correct_relu, size);
    }

}


void DRLU_Test(Party *proxy) {
    cout << setfill('*') << setw(50) << "Calling DRLU";
    cout << setfill('*') << setw(49) << "*" << endl;

    uint64_t x = proxy->createShare(convert2double(proxy->generateRandom()));

    proxy->SendBytes(CNN_DRLU);
    uint64_t drelu = DRELU(proxy, x);
    uint64_t reconstructed_drelu = REC(proxy, drelu);

    // checking the result
    double originalX = convert2double(REC(proxy, x));
    uint64_t computed_drelu = 0;
    if (originalX > 0)
        computed_drelu = 1;

    if (computed_drelu == reconstructed_drelu) {
        cout << "DRLU works correctly" << endl;
    } else {
        cout << "DRLU works incorrectly" << endl;
        cout << "computed: " << reconstructed_drelu << " should be: " << computed_drelu << endl;
    }

}

void MDRLU_Test(Party *proxy) {
    cout << setfill('*') << setw(50) << "Calling MDRLU";
    cout << setfill('*') << setw(49) << "*" << endl;

    uint64_t size = 10;
    uint32_t *params = new uint32_t[1];
    params[0] = size;
    uint64_t *x = proxy->createShare(random_1D_data(proxy, size), size);

    proxy->SendBytes(CNN_MDRLU, params, 1);
    uint64_t *drelu = DRELU(proxy, x, size);
    double *rec_drelu = convert2double(REC(proxy, drelu, size), size);

    // checking the result
    double *originalX = convert2double(REC(proxy, x, size), size);
    double *correct_drelu = new double[size];
    bool allCorrect = true;
    for (uint64_t i = 0; i < size; i++) {
        if (originalX[i] > 0) {
            correct_drelu[i] = 1;
        } else {
            correct_drelu[i] = 0;
        }
        if (abs(correct_drelu[i] - rec_drelu[i]) > 0.0001) {
            allCorrect = false;
        }
    }
    if (allCorrect) {
        cout << "MDRLU works correctly" << endl;
    } else {
        cout << "MDRLU works incorrectly" << endl;
        print1DArray("Original Values:", originalX, size);
        print1DArray("Computed DRELU:", rec_drelu, size);
        print1DArray("VS Correct DRELU:", correct_drelu, size);
    }

}


void PAD_Test(Party *proxy) {
    cout << setfill('*') << setw(50) << "Calling PAD (padding)";
    cout << setfill('*') << setw(49) << "*" << endl;

    uint32_t row = 28;
    uint32_t col = 28;
    uint32_t padding_value = 0;
    uint32_t padding_size = 2;
    uint64_t **x = proxy->createShare(random_2D_data(proxy, row, col, 0.0, 255.0), row, col);
    //print2DArray("Original matrix: ", convert2double(REC(proxy, x, row, col), row, col), row, col);
    uint64_t **padded = PAD(x, row, col, padding_value, padding_size);
    double **rec_padded = convert2double(REC(proxy, padded, row + 2 * padding_size, col + 2 * padding_size),
                                         row + 2 * padding_size, col + 2 * padding_size);

    double **computed_padding = new double *[row + 2 * padding_size];
    for (uint32_t p = 0; p < padding_size; p++) {
        computed_padding[p] = new double[col + 2 * padding_size];
        computed_padding[p + padding_size + row] = new double[col + 2 * padding_size];

        for (uint32_t c = 0; c < (col + 2 * padding_size); c++) {
            computed_padding[p][c] = padding_value;                      //top
            computed_padding[p + padding_size + row][c] = padding_value; //bottom
        }
    }
    for (uint64_t r = 0; r < row; r++) {
        computed_padding[padding_size + r] = new double[col + 2 * padding_size];
        for (uint32_t p = 0; p < padding_size; p++) {
            computed_padding[padding_size + r][p] = padding_value;
            computed_padding[padding_size + r][padding_size + col + p] = padding_value;
        }
        for (uint32_t c = 0; c < col; c++) {
            computed_padding[padding_size + r][c + padding_size] = convert2double(REC(proxy, x[r][c]));
        }
    }

    // checking the result
    bool allCorrect = true;
    for (int i = 0; i < (2 * padding_size + row); ++i) {
        for (int j = 0; j < (2 * padding_size + col); ++j) {
            if (abs(computed_padding[i][j] - rec_padded[i][j]) > 0.0001) {
                allCorrect = false;
            }
        }
    }
    if (allCorrect) {
        cout << "PAD works correctly" << endl;
        print2DArray("Computed Padding", rec_padded, row + 2 * padding_size, col + 2 * padding_size);
    } else {
        cout << "PAD works incorrectly" << endl;
        print2DArray("Computed Padding", rec_padded, row + 2 * padding_size, col + 2 * padding_size);
        print2DArray("Correct Padding", computed_padding, row + 2 * padding_size, col + 2 * padding_size);
    }

}

void CL_Test(Party *proxy) {
    cout << setfill('*') << setw(50) << "Calling CL (convolutional layer)";
    cout << setfill('*') << setw(49) << "*" << endl;
    // set params
    uint32_t row = 8;
    uint32_t col = 8;
    uint32_t i_channel = 1;
    uint32_t k_number = 2;
    uint32_t k_dim = 5;
    uint32_t stride = 1;
    uint32_t max_width = 2;
    uint32_t max_height = 2;
    bool flatten_result = true;

    // generate secret shares of random data
    auto ***x = new uint64_t **[i_channel];
    auto ***kernel = new uint64_t **[k_number];
    for (uint64_t i = 0; i < i_channel; i++) {
        x[i] = proxy->createShare(random_2D_data(proxy, row, col, -10.0, 10.0), row, col);
        //print2DArray("input", convert2double(REC(proxy, x[i], row, col), row, col), row, col);
        kernel[i] = proxy->createShare(random_2D_data(proxy, k_number, k_dim * k_dim, -5.0, 5.0), k_number,
                                       k_dim * k_dim);
        //print2DArray("kernel", convert2double(REC(proxy, kernel[i], k_number, k_dim*k_dim), k_number, k_dim*k_dim), k_number, k_dim*k_dim);
    }
    uint64_t conv_width = floor((col - k_dim + 1) / stride);
    uint64_t conv_height = floor((row - k_dim + 1) / stride);

    double *rec_bias = random_1D_data(proxy, k_number);
    uint64_t *bias = proxy->createShare(rec_bias, k_number);
    // send params
    uint32_t params[9];
    params[0] = i_channel;
    params[1] = row;
    params[2] = col;
    params[3] = k_dim;      // kernel size
    params[4] = k_number;   // kernel number = output channel
    params[5] = stride;
    params[6] = max_height;
    params[7] = max_width;
    params[8] = flatten_result;
    proxy->SendBytes(CNN_CL, params, 9);

    uint64_t ***conv = CL(proxy, x, i_channel, row, col, kernel, k_dim, k_number, stride, max_height, max_width,
                          bias, flatten_result);
    delete[] bias;
    uint64_t lastpos = row - k_dim + 1;
    uint32_t out_height;
    if (max_height > 1) {
        out_height = conv_height / max_height;
    }
    uint32_t out_width;
    if (max_width > 1) {
        out_width = conv_width / max_width;
    }
    //print2DArray("CL result, first output channel", convert2double(REC(proxy, conv[0], out_height, out_width), out_height, out_width), out_height, out_width);
    //print2DArray("CL result, second output channel", convert2double(REC(proxy, conv[1], out_height, out_width), out_height, out_width), out_height, out_width);

    // checking the result
    double ***correct_conv = new double **[k_number];
    double ***pooled_conv = new double **[k_number];
    bool allCorrect = true;
    for (uint32_t kern = 0; kern < k_number; kern++) {    // for each output channel
        correct_conv[kern] = new double *[conv_height];                // init kernels conv result
        for (uint32_t cr = 0; cr < lastpos; cr += stride) {
            correct_conv[kern][cr / stride] = new double[conv_width];   // init row of conv result
            for (uint32_t cc = 0; cc < lastpos; cc += stride) {
                double dot_product = 0;
                for (uint32_t channel = 0; channel < i_channel; channel++) {
                    for (uint32_t k_value = 0; k_value < (k_dim * k_dim); k_value++) {
                        double v = convert2double(
                                REC(proxy, x[channel][cr + k_value / k_dim][cc + k_value % k_dim]));
                        double weight = convert2double(REC(proxy, kernel[channel][kern][k_value]));
                        dot_product += v * weight;
                    }
                }
                dot_product += rec_bias[kern];
                // Activation: RELU
                double relu = 0.0;
                if (dot_product > 0) {
                    relu = dot_product;
                }
                correct_conv[kern][cr / stride][cc / stride] = relu;
                if (max_height <= 1 and max_width <= 1) {
                    // no maxpooling --> check if correct
                    double comp_value = 0;
                    if (flatten_result) {
                        comp_value = convert2double(REC(proxy, conv[0][0][kern * conv_height * conv_width +
                                                                          (cr / stride * conv_width) +
                                                                          cc / stride]));
                    } else {
                        comp_value = convert2double(REC(proxy, conv[kern][cr / stride][cc / stride]));
                    }
                    if (abs(relu - comp_value) > 0.0001) {
                        allCorrect = false;
                    }
                }
            }
        }
        //print2DArray("activated conv", correct_conv[kern], conv_height, conv_width);
        if (max_height > 0 and max_width > 0 and (max_height + max_width) > 2) {
            //MAXPOOLING
            if (flatten_result) {
                pooled_conv[kern] = new double *[0];
                pooled_conv[kern][0] = new double[out_height * out_width];
            } else {
                pooled_conv[kern] = new double *[out_height];
            }
            for (int r = 0; r < (conv_height - max_height + 1); r += max_height) {
                if (!flatten_result) {
                    pooled_conv[kern][r / max_height] = new double[out_width];
                }
                for (int c = 0; c < (conv_width - max_width + 1); c += max_width) {
                    //find max of window:
                    double max = correct_conv[kern][r][c]; // first value in window
                    for (int max_r = 0; max_r < max_height; ++max_r) {
                        for (int max_c = 0; max_c < max_width; ++max_c) {
                            double next_value = correct_conv[kern][r + max_r][c + max_c];
                            if (next_value > max) {
                                max = next_value;
                            }
                        }
                    }
                    double comp_max = 0;
                    if (flatten_result) {
                        pooled_conv[kern][0][r / max_height * max_width + c / max_width] = max;
                        comp_max = convert2double(REC(proxy, conv[0][0][kern * max_height * max_width +
                                                                        (r / max_height * max_width) + c / max_width]));
                    } else {
                        pooled_conv[kern][r / max_height][c / max_width] = max;
                        comp_max = convert2double(REC(proxy, conv[kern][r / max_height][c / max_width]));
                    }
                    if (abs(max - comp_max) > 0.001) {
                        allCorrect = false;
                    }
                }
            }
        }
    }
    if (allCorrect) {
        cout << "CL works correctly" << endl;
    } else {
        cout << "CL works incorrectly" << endl;
        /*cout << setfill('=') << setw(25) << "GIVEN PARAMETERS";
        cout << setfill('=') << setw(24) << "=" << endl;
        cout << "Input channel of X" << endl;
        for (int i = 0; i < i_channel; ++i) {
            print2DArray("x", convert2double(REC(proxy, x[i], row, col), row, col), row, col);
            cout << "Output channel of kernel " << endl;
            print2DArray("Channel of kernel",
                         convert2double(REC(proxy, kernel[i], k_number, k_dim * k_dim), k_number, k_dim * k_dim), k_number, k_dim * k_dim);
        }*/
        cout << setfill('=') << setw(25) << "OUTPUT PARAMETERS";
        cout << setfill('*') << setw(24) << "*" << endl;
        if (flatten_result) {
            if (max_width > 1 and max_height > 1) {
                print1DArray("Flattened, maxpooled computed convolution",
                             convert2double(REC(proxy, conv[0][0], out_height * out_width * k_number),
                                            out_height * out_width * k_number), out_height * out_width * k_number);
                for (int k = 0; k < k_number; ++k) {
                    print1DArray("Correct, pooled convolution", pooled_conv[k][0], out_height * out_width);
                }
            } else {
                print1DArray("Flattened computed convolution",
                             convert2double(REC(proxy, conv[0][0], conv_height * conv_width * k_number),
                                            conv_height * conv_width * k_number), conv_height * conv_width * k_number);
                for (int k = 0; k < k_number; ++k) {
                    print1DArray("Correct, pooled convolution", pooled_conv[k][0], conv_height * conv_width);
                }
            }
        } else {
            for (uint64_t k = 0; k < k_number; k++) {
                cout << "Output channel " << k << endl;
                if (max_width > 0 and max_height > 0) {
                    print2DArray("Computed convolution",
                                 convert2double(REC(proxy, conv[k], out_height, out_width), out_height, out_width),
                                 out_height, out_width);
                    print2DArray("Correct, pooled convolution", pooled_conv[k], out_height, out_width);
                } else {
                    print2DArray("Computed convolution",
                                 convert2double(REC(proxy, conv[k], conv_height, conv_width), conv_height, conv_width),
                                 conv_height, conv_width);
                    print2DArray("Correct convolution", correct_conv[k], conv_height, conv_width);
                }
            }
        }

    }

    for (int i = 0; i < i_channel; ++i) {
        for (int r = 0; r < row; ++r) {
            delete[] x[i][r];
        }
        for (int k = 0; k < k_number; ++k) {
            delete[] kernel[i][k];
        }
        delete[] x[i];
        delete[] kernel[i];
    }
    for (int k = 0; k < k_number; ++k) {
        for (int c = 0; c < conv_height; ++c) {
            delete[] correct_conv[k][c];
        }
        if (max_width > 0 and max_width > 0) {
            for (int m = 0; m < (conv_height - max_height) / max_height; ++m) {
                delete[] pooled_conv[k][m];
            }
        }
        delete[] correct_conv[k];
        delete[] pooled_conv[k];
    }
    delete[] kernel;
    delete[] correct_conv;
    delete[] pooled_conv;

}

void FCL_Test(Party *proxy) {
    cout << setfill('*') << setw(50) << "Calling FCL (fully connected layer)";
    cout << setfill('*') << setw(49) << "*" << endl;

    // init
    uint32_t row = 4;
    uint32_t col = 4;
    uint32_t i_channel = 50;
    uint32_t i_nodes = row * col * i_channel;
    uint32_t o_nodes = 2;

    bool allCorrect = true;
    while (allCorrect) {
        uint64_t *x = proxy->createShare(random_1D_data(proxy, i_nodes, 0.0, 255.0), i_nodes);
        uint64_t **weights = proxy->createShare(random_2D_data(proxy, i_nodes, o_nodes, -5.0, 5.0), i_nodes, o_nodes);
        uint64_t *bias = proxy->createShare(random_1D_data(proxy, o_nodes), o_nodes);

        // send params
        uint32_t params[2];
        params[0] = i_nodes;
        params[1] = o_nodes;
        proxy->SendBytes(CNN_FCL, params, 2);

        uint64_t *res = FCL(proxy, x, i_nodes, weights, o_nodes, bias);
        double *reconstructed_res = convert2double(REC(proxy, res, o_nodes), o_nodes);

        // checking the result
        double *rec_input = convert2double(REC(proxy, x, i_nodes), i_nodes);
        //print1DArray("reconstructed input: ", rec_input, i_nodes);
        double **rec_weights = convert2double(REC(proxy, weights, i_nodes, o_nodes), i_nodes, o_nodes);
        //print2DArray("reconstructed weights: ", rec_weights, i_nodes, o_nodes);
        double *rec_bias = convert2double(REC(proxy, bias, o_nodes), o_nodes);
        //print1DArray("reconstructed bias: ", rec_bias, o_nodes);

        double *original_res = new double[o_nodes];
        for (uint64_t o = 0; o < o_nodes; o++) {
            double value = 0;
            for (int i = 0; i < i_nodes; i++) {
                value += rec_input[i] * rec_weights[o][i]; //weights would be tansposed in FCL
            }
            //ReLU activation
            if (value < 0) {
                value = 0;
            }
            original_res[o] = value + rec_bias[o];
            if (original_res[o] - reconstructed_res[o] > 0.0001) {
                allCorrect = false;
            }
        }

        if (allCorrect) {
            cout << "FCL works correctly" << endl;
        } else {
            cout << "FCL works incorrectly" << endl;
            print1DArray("Computed Result: ", reconstructed_res, o_nodes);
            print1DArray("Correct Result", original_res, o_nodes);
            print1DArray("input", convert2double(REC(proxy, x, i_nodes), i_nodes), i_nodes);
            print2DArray("weights", convert2double(REC(proxy, weights, i_nodes, o_nodes), i_nodes, o_nodes), i_nodes,
                         o_nodes);
            print1DArray("bias", convert2double(REC(proxy, bias, o_nodes), o_nodes), o_nodes);
        }
    }
}

void FLT_Test(Party *proxy) {
    cout << setfill('*') << setw(50) << "Calling FLT (flattening)";
    cout << setfill('*') << setw(49) << "*" << endl;

    uint32_t row = 5;
    uint32_t col = 5;
    uint32_t i_channel = 3;

    uint64_t ***x = new uint64_t **[i_channel];
    for (uint64_t i = 0; i < i_channel; i++) {
        x[i] = proxy->createShare(random_2D_data(proxy, row, col, 0, 255), row, col);
    }
    uint32_t i_nodes = row * col * i_channel;

    uint64_t *res = FLT(x, row, col, i_channel);
    double *reconstructed_res = convert2double(REC(proxy, res, i_nodes), i_nodes);

    // checking the result
    bool allCorrect = true;
    double *original_res = new double[i_nodes];
    for (int channel = 0; channel < i_channel; ++channel) {
        for (int r = 0; r < row; ++r) {
            for (int c = 0; c < col; ++c) {
                uint64_t position = channel * row * col + r * col + c;
                original_res[position] = convert2double(REC(proxy, x[channel][r][c]));
                if ((original_res[position] - reconstructed_res[position]) > 0.0001) {
                    allCorrect = false;
                }
            }
        }
    }

    if (allCorrect) {
        cout << "FLT works correctly" << endl;
    } else {
        cout << "FLT works incorrectly" << endl;
        print1DArray("Computed Result: ", reconstructed_res, i_nodes);
        print1DArray("Correct Result", original_res, i_nodes);
    }
}

void INC_Test(Party *proxy) {
    cout << setfill('*') << setw(50) << "Calling INC (increasing input matrix for conv)";
    cout << setfill('*') << setw(49) << "*" << endl;

    uint32_t row = 28;
    uint32_t col = 28;
    uint32_t i_channel = 1;
    uint64_t ***x = new uint64_t **[i_channel];
    for (uint64_t i = 0; i < i_channel; i++) {
        x[i] = proxy->createShare(random_2D_data(proxy, row, col, 0, 255), row, col);
    }
    //print2DArray("Original X: ", convert2double(REC(proxy, x[0], row, col), row, col), row, col);

    uint32_t k_number = 5;
    uint32_t k_dim = 5;
    uint32_t k_len = k_dim * k_dim;
    uint32_t stride = 1;

    uint32_t lastpos = row - k_dim + 1;
    uint32_t conv_height = static_cast<uint32_t>(floor(lastpos / stride));
    uint32_t conv_width = static_cast<uint32_t>(floor((col - k_dim + 1) / stride));

    uint64_t ***stretchedX = INC(x, i_channel, row, col, k_dim, stride);
    // checking the result
    if (row * col < 100) {
        print2DArray("RESULTED STRETCHED MATRIX (first channel)",
                     convert2double(REC(proxy, stretchedX[0], conv_height * conv_width, k_len),
                                    conv_height * conv_width, k_len), conv_height * conv_width, k_len);
    }
}

void EXP_Test(Party *proxy) {
    cout << setfill('*') << setw(50) << "Calling EXP";
    cout << setfill('*') << setw(49) << "*" << endl;

    cout << "Min power: " << proxy->getMinPower() << endl;
    cout << "Max power: " << proxy->getMaxPower() << endl;

    uint64_t x = proxy->createShare(random_1D_data(proxy, 1, proxy->getMinPower(), proxy->getMaxPower())[0]);

    proxy->SendBytes(CORE_EXP);
    uint64_t shr_exp = EXP(proxy, x);
    uint64_t reconstructed_exp = REC(proxy, shr_exp);
    double rec_exp = convert2double(reconstructed_exp);

    // checking the result
    double originalX = convert2double(REC(proxy, x));
    double true_exp = exp(originalX);

    if (abs(true_exp - rec_exp) <= 0.1) {
        cout << "EXP works correctly" << endl;
        cout << "power: " << originalX << " -- computed: " << rec_exp << " should be: " << true_exp << endl;
    } else {
        cout << "EXP works incorrectly" << endl;
        cout << "power: " << originalX << " -- computed: " << rec_exp << " should be: " << true_exp << endl;
    }

}

void MEXP_Test(Party *proxy) {
    cout << setfill('*') << setw(50) << "Calling MEXP";
    cout << setfill('*') << setw(49) << "*" << endl;

    uint32_t params[1];
    int n_samples = 10;
    params[0] = n_samples;
    uint64_t *x = proxy->createShare(random_1D_data(proxy, n_samples, proxy->getMinPower(), proxy->getMaxPower()),
                                     n_samples);
    proxy->SendBytes(CORE_MEXP, params, 1);
    uint64_t *shr_exp = EXP(proxy, x, n_samples);
    uint64_t *reconstructed_exp = REC(proxy, shr_exp, n_samples);
    double *rec_exp = convert2double(reconstructed_exp, n_samples);

    // checking the result
    double *originalX = convert2double(REC(proxy, x, n_samples), n_samples);
    double *true_exp = new double[n_samples];
    for (int i = 0; i < n_samples; i++) {
        true_exp[i] = exp(originalX[i]);
    }

    bool flag = true;
    for (int i = 0; i < n_samples; i++) {
        cout << fixed << "power: " << originalX[i] << " -- computed: " << rec_exp[i] << " - expected: " << true_exp[i]
             <<
             " - absolute difference: " << abs(rec_exp[i] - true_exp[i]) << endl;
        if (abs(true_exp[i] - rec_exp[i]) >= 0.1) {
            flag = false;
        }
    }
    if (flag) {
        cout << "EXP works correctly" << endl;
    } else {
        cout << "EXP works incorrectly" << endl;
    }

}

void DP_Test(Party *proxy) {
    cout << setfill('*') << setw(50) << "Calling DP";
    cout << setfill('*') << setw(49) << "*" << endl;

    uint32_t *params = new uint32_t[1];
    uint32_t size = 20; // size of the vector
    params[0] = size;
    double min_val = -10;
    double max_val = 10;

    // generate a vector of random values
    uint64_t *vec1 = proxy->createShare(random_1D_data(proxy, size, min_val, max_val), size);
    uint64_t *vec2 = proxy->createShare(random_1D_data(proxy, size, min_val, max_val), size);

    // call DP
    proxy->SendBytes(CORE_DP, params, 1);
    uint64_t res = DP(proxy, vec1, vec2, size);

    // check the results
    double gt = 0;
    double *rec_vec1 = convert2double(REC(proxy, vec1, size), size);
    double *rec_vec2 = convert2double(REC(proxy, vec2, size), size);
    for (int i = 0; i < size; i++) {
        gt += rec_vec1[i] * rec_vec2[i];
    }

    double rec_res = convert2double(REC(proxy, res));

    printValue("Computed dot product", rec_res);
    printValue("GT dot product", gt);
}

void MDP_Test(Party *proxy) {
    cout << setfill('*') << setw(50) << "Calling MDP";
    cout << setfill('*') << setw(49) << "*" << endl;

    uint32_t *params = new uint32_t[1];
    uint32_t size = 20; // size of the total vector
    params[0] = size;

    uint32_t d = 5; // size of each individual vector in the main vector
    double min_val = -10;
    double max_val = 10;

    if (size % d != 0) {
        throw invalid_argument("DP_Test: The size must be divisible by d.");
    }

    // generate a vector of random values
    uint64_t *vec1 = proxy->createShare(random_1D_data(proxy, size, min_val, max_val), size);
    uint64_t *vec2 = proxy->createShare(random_1D_data(proxy, size, min_val, max_val), size);

    // call DP
    proxy->SendBytes(CORE_MDP, params, 1);
    uint64_t *res = DP(proxy, vec1, vec2, size, d);

    // check the results
    double *gt = new double[size / d];
    double *rec_vec1 = convert2double(REC(proxy, vec1, size), size);
    double *rec_vec2 = convert2double(REC(proxy, vec2, size), size);
    for (int i = 0; i < size; i += d) {
        double tmp_sum = 0;
        for (int j = i; j < i + d; j++) {
            tmp_sum += rec_vec1[j] * rec_vec2[j];
        }
        gt[i / d] = tmp_sum;
    }

    double *rec_res = convert2double(REC(proxy, res, size / d), size / d);

    print1DArray("Computed dot product", rec_res, size / d);
    print1DArray("GT dot product", gt, size / d);
}

bool MATMATMUL_Test(Party *proxy) {
    cout << setfill('*') << setw(50) << "Calling MATMATMUL";
    cout << setfill('*') << setw(49) << "*" << endl;

    // setting
    int a_row = 3;
    int a_col = 10;
    int b_col = 2;

    uint32_t *params = new uint32_t[1];
    uint32_t size = a_row * a_col * b_col; // size of the vector
    params[0] = size;
    double min_val = -0.784378;
    double max_val = 1481.76;

    uint64_t **mat1 = proxy->createShare(random_2D_data(proxy, a_row, a_col, min_val, max_val), a_row, a_col);
    uint64_t **mat2 = proxy->createShare(random_2D_data(proxy, a_col, b_col, min_val, max_val), a_col, b_col);

    proxy->SendBytes(CORE_MATMATMUL, params, 1);
    uint64_t **res = MATMATMUL(proxy, mat1, mat2, a_row, a_col, b_col);
    double **rec_res = convert2double(REC(proxy, res, a_row, b_col), a_row, b_col);

    double **rec_mat1 = convert2double(REC(proxy, mat1, a_row, a_col), a_row, a_col);
    double **rec_mat2 = convert2double(REC(proxy, mat2, a_col, b_col), a_col, b_col);

//    proxy->print2DArray("mat1", rec_mat1, a_row, a_col);
//    proxy->print2DArray("mat2", rec_mat2, a_col, b_col);

    double **gt = multiply_matrices(rec_mat1, rec_mat2, a_row, a_col, b_col);

    double tmp = 0;
    for (int i = 0; i < a_row; i++) {
        for (int j = 0; j < b_col; j++) {
            tmp += abs(rec_res[i][j] - gt[i][j]);
        }
    }

    if (tmp <= 0.1) {
        cout << "MATMATMUL works correctly" << endl;
        //cout << "Total absolute difference: " << tmp << endl;
    } else {
        cout << "MATMATMUL works incorrectly" << endl;
        cout << "Total absolute difference: " << tmp << endl;

        print2DArray("mat", rec_mat1, a_row, a_col);
        print2DArray("mat", rec_mat2, a_col, b_col);

        print2DArray("Computed matrix multiplication", rec_res, a_row, b_col);
        print2DArray("GT matrix multiplication", gt, a_row, b_col);
    }
    return (tmp <= 0.1);
}

void MMATMATMUL_Test(Party *proxy) {
    // In this function, we test several matrix multiplications of random matrices with the same size.
    cout << setfill('*') << setw(50) << "Calling MMATMATMUL";
    cout << setfill('*') << setw(49) << "*" << endl;

    // setting
    int n_matrices = 3;
    int a_row = 5;
    int a_col = 6;
    int b_col = 3;
    uint32_t *params = new uint32_t[1];
    uint32_t size = n_matrices * a_row * a_col * b_col; // size of the vector
    params[0] = size;
    double min_val = -10;
    double max_val = 10;

    uint64_t ***mat1 = new uint64_t **[n_matrices];
    uint64_t ***mat2 = new uint64_t **[n_matrices];
    for (int i = 0; i < n_matrices; i++) {
        mat1[i] = proxy->createShare(random_2D_data(proxy, a_row, a_col, min_val, max_val), a_row, a_col);
        mat2[i] = proxy->createShare(random_2D_data(proxy, a_col, b_col, min_val, max_val), a_col, b_col);
    }

    proxy->SendBytes(CORE_MMATMATMUL, params, 1);
    uint64_t ***res = MATMATMUL(proxy, mat1, mat2, n_matrices, a_row, a_col, b_col);

    double ***rec_mat1 = new double **[n_matrices];
    double ***rec_mat2 = new double **[n_matrices];
    for (int i = 0; i < n_matrices; i++) {
        rec_mat1[i] = convert2double(REC(proxy, mat1[i], a_row, a_col), a_row, a_col);
        rec_mat2[i] = convert2double(REC(proxy, mat2[i], a_col, b_col), a_col, b_col);
//        proxy->print2DArray("mat1's reconstructed matrix " + to_string(i), rec_mat1[i], a_row, a_col);
//        proxy->print2DArray("mat2's reconstructed matrix " + to_string(i), rec_mat2[i], a_col, b_col);
    }

//    for(int i = 0; i < n_matrices; i++) {
//        print2DArray("Computed matrix multiplication - " + to_string(i), convert2double(
//                REC(proxy, res[i], a_row, b_col), a_row, b_col), a_row, b_col);
//        double** gt = multiply_matrices(rec_mat1[i], rec_mat2[i], a_row, a_col, b_col);
//        print2DArray("GT matrix multiplication - " + to_string(i), gt, a_row, b_col);
//    }

    double *diffs = new double[n_matrices];
    bool flag = true;
    for (int m = 0; m < n_matrices; m++) {
        double **rec_res = convert2double(REC(proxy, res[m], a_row, b_col), a_row, b_col);
        double **gt = multiply_matrices(rec_mat1[m], rec_mat2[m], a_row, a_col, b_col);

        diffs[m] = 0;
        for (int i = 0; i < a_row; i++) {
            for (int j = 0; j < b_col; j++) {
                diffs[m] += abs(rec_res[i][j] - gt[i][j]);
            }
        }

        if (diffs[m] >= 0.1) {
            flag = false;
        }
    }

    if (flag) {
        cout << "MMATMATMUL works correctly" << endl;
        cout << "Total absolute differences: ";
        for (int m = 0; m < n_matrices; m++) {
            cout << diffs[m] << "\t";
        }
        cout << endl;
    } else {
        cout << "MMATMATMUL works incorrectly" << endl;
        cout << "Total absolute differences: ";
        for (int m = 0; m < n_matrices; m++) {
            cout << diffs[m] << "\t";
        }
        cout << endl;
    }
}

bool MATVECMUL_Test(Party *proxy) {
    /* In this function, we test the matrix multiplication of two random matrices. We first generate two
     * matrices of random values such that the number of column of the first matrix equals to the number
     * of rows of the second matrix.
     */
    cout << setfill('*') << setw(50) << "Calling MATVECMUL";
    cout << setfill('*') << setw(49) << "*" << endl;
    cout << "calling this one" << endl;
    // setting
    int a_row = 3;
    int a_col = 10;
    uint32_t *params = new uint32_t[1];
    uint32_t size = a_row * a_col; // size of the vector
    params[0] = size;
    double min_val = -0.784378;
    double max_val = 120.76;

    uint64_t **mat = proxy->createShare(random_2D_data(proxy, a_row, a_col, min_val, max_val), a_row, a_col);
    uint64_t *vec = proxy->createShare(random_1D_data(proxy, a_col, min_val, max_val), a_col);

    proxy->SendBytes(CORE_MATVECMUL, params, 1);
    uint64_t *res = MATVECMUL(proxy, mat, vec, a_row, a_col);
    double *rec_res = convert2double(REC(proxy, res, a_row), a_row);

    double **rec_mat = convert2double(REC(proxy, mat, a_row, a_col), a_row, a_col);
    double *rec_vec = convert2double(REC(proxy, vec, a_col), a_col);

    double *gt = multiply_matrice_vector(rec_mat, rec_vec, a_row, a_col);

    double tmp = 0;
    for (int i = 0; i < a_row; i++) {
        tmp += abs(rec_res[i] - gt[i]);
    }

    if (tmp <= 0.1) {
        cout << "MATVECMUL works correctly" << endl;
        //cout << "Total absolute difference: " << tmp << endl;
    } else {
        cout << "MATVECMUL works incorrectly" << endl;
        cout << "Total absolute difference: " << tmp << endl;

        print2DArray("mat", rec_mat, a_row, a_col);
        print1DArray("vec", rec_vec, a_col);

        print1DArray("Computed matrix-vector multiplication", convert2double(REC(proxy, res, a_row), a_row), a_row);
        print1DArray("GT matrix-vector multiplication", gt, a_row);
    }
    return tmp <= 0.1;
}

void MMATVECMUL_Test(Party *proxy) {
    // In this function, we test the matrix multiplication of n_matrices matrices and vectors.
    cout << setfill('*') << setw(50) << "Calling MMATVECMUL";
    cout << setfill('*') << setw(49) << "*" << endl;

    // setting
    int n_matrices = 3;
    int a_row = 5;
    int a_col = 6;
    int b_col = 3;
    uint32_t *params = new uint32_t[1];
    uint32_t size = n_matrices * a_row * a_col; // size of the vector
    params[0] = size;
    double min_val = -0.5;
    double max_val = 99.0;

    uint64_t ***mat = new uint64_t **[n_matrices];
    uint64_t **vec = new uint64_t *[n_matrices];
    for (int i = 0; i < n_matrices; i++) {
        mat[i] = proxy->createShare(random_2D_data(proxy, a_row, a_col, min_val, max_val), a_row, a_col);
        vec[i] = proxy->createShare(random_1D_data(proxy, a_col, min_val, max_val), a_col);
    }

    proxy->SendBytes(CORE_MMATVECMUL, params, 1);
    uint64_t **res = MATVECMUL(proxy, mat, vec, n_matrices, a_row, a_col);

    double ***rec_mat = new double **[n_matrices];
    double **rec_vec = new double *[n_matrices];
    for (int i = 0; i < n_matrices; i++) {
        rec_mat[i] = convert2double(REC(proxy, mat[i], a_row, a_col), a_row, a_col);
        rec_vec[i] = convert2double(REC(proxy, vec[i], a_col), a_col);
    }

//    for(int i = 0; i < n_matrices; i++) {
//        print1DArray("Computed matrix-vector multiplication - " + to_string(i), convert2double(REC(proxy, res[i], a_row), a_row), a_row);
//        double* gt = multiply_matrice_vector(rec_mat[i], rec_vec[i], a_row, a_col);
//        print1DArray("Ground truth matrix-vector multiplication - " + to_string(i), gt, a_row);
//    }

    double *diffs = new double[n_matrices];
    bool flag = true;
    for (int m = 0; m < n_matrices; m++) {
        double *rec_res = convert2double(REC(proxy, res[m], a_row), a_row);
        double *gt = multiply_matrice_vector(rec_mat[m], rec_vec[m], a_row, a_col);

        diffs[m] = 0;
        for (int i = 0; i < a_row; i++) {
            diffs[m] += abs(rec_res[i] - gt[i]);
        }

        if (diffs[m] >= 0.1) {
            flag = false;
        }
    }

    if (flag) {
        cout << "MMATVECMUL works correctly" << endl;
        cout << "Total absolute difference: " << endl;
        for (int m = 0; m < n_matrices; m++) {
            cout << diffs[m] << "\t";
        }
        cout << endl;
    } else {
        cout << "MMATVECMUL works incorrectly" << endl;
        cout << "Total absolute difference: " << endl;
        for (int m = 0; m < n_matrices; m++) {
            cout << diffs[m] << "\t";
        }
        cout << endl;
    }
}

//void INVSQRT_Test(Party* proxy) {
//    /* In this function, we test the computation of the inverse square root of a Gram matrix.
//     * We first generate a random Gram matrix by first generating a random data matrix D and
//     * then computing D^T * D.
//     */
//
//    cout<<setfill ('*')<<setw(50)<<"Calling INVSQRT";
//    cout<<setfill ('*')<<setw(49)<<"*"<<endl;
//    // setting
//    int n_row = 4;
//    int n_col = 5;
//
//    // generate a Gram matrix
//    uint64_t **invsqrt_data = random_gram_matrix(proxy, n_row, n_col);
//
//    double** tmp = convert2double(REC(proxy, invsqrt_data, n_row, n_row), n_row, n_row);
//
//    proxy->SendBytes(RKN_INVSQRT, n_row);
//    uint64_t** invsqrt_G = INVSQRT(proxy, invsqrt_data, n_row);
//
//    double** rec_invsqrt_G = convert2double(REC(proxy, invsqrt_G, n_row, n_row), n_row, n_row);
//
//    print2DArray("The inverse square root of the Gram matrix", rec_invsqrt_G, n_row, n_row);
//
//    double* straighten_invsqrt_G = new double[n_row * n_row];
//    for(uint32_t i = 0; i < n_row * n_row; i++) {
//        straighten_invsqrt_G[i] = tmp[i % n_row][i / n_row];
//    }
//
//    print2DArray("Gram matrix", tmp, n_row, n_row);
//    print1DArray("Straighten Gram matrix", straighten_invsqrt_G, n_row * n_row);
//
//    EigenSolver<Matrix<double, Dynamic, Dynamic>> ges;
//    Map<Matrix<double, Dynamic, Dynamic>> matrix_G(straighten_invsqrt_G, n_row, n_row);
//    ges.compute(matrix_G);
//    Matrix<double, Dynamic, Dynamic> eig_vecs = ges.eigenvectors().real();
//    Matrix<double, Dynamic, 1> eig_vals = ges.eigenvalues().real();
//
//    cout << "============= GT the eigenvalues ======================" << endl;
//    cout << eig_vals << endl;
//    cout << "============================================================================" << endl;
//
//    cout << "============= GT Inverse square root of the eigenvalues ======================" << endl;
//    Matrix<double, Dynamic, Dynamic> vals = eig_vals;
//    cout << vals.cwiseSqrt().cwiseInverse() << endl;
//    cout << "============================================================================" << endl;
//
//    cout << "============= GT Inverse square root of the Gram matrix ======================" << endl;
//    cout << eig_vecs * vals.cwiseSqrt().cwiseInverse().asDiagonal() * Transpose(eig_vecs) << endl;
//    cout << "============================================================================" << endl;
//
//    print2DArray("The inverse of the Gram matrix",
//                        multiply_matrices(rec_invsqrt_G, rec_invsqrt_G, n_row, n_row, n_row), n_row, n_row);
//
//    cout << "============= GT Inverse of the Gram matrix ======================" << endl;
//    cout << matrix_G.inverse() << endl;
//    cout << "============================================================================" << endl;
//
//}

//void MINVSQRT_Test(Party* proxy) {
//    /* In this function, we test the computation of the inverse square root of a Gram matrix.
//     * We first generate a random Gram matrix by first generating a random data matrix D and
//     * then computing D^T * D.
//     */
//    cout<<setfill ('*')<<setw(50)<<"Calling MINVSQRT";
//    cout<<setfill ('*')<<setw(49)<<"*"<<endl;
//    // setting
//    int n_row = 4;
//    int n_col = 5;
//    int n_gms = 3;
//
//    // generate a Gram matrix
//    uint64_t ***invsqrt_data = new uint64_t**[n_gms];
//    double ***Gs = new double**[n_gms];
//    for(int i = 0; i < n_gms; i++) {
//        invsqrt_data[i] = random_gram_matrix(proxy, n_row, n_col);
//        Gs[i] = convert2double(REC(proxy, invsqrt_data[i], n_row, n_row), n_row, n_row);
//    }
//
////    double** tmp = convert2double(REC(proxy, invsqrt_data, n_row, n_row), n_row, n_row);
//
//    proxy->SendBytes(RKN_MINVSQRT, n_gms,n_row);
//    uint64_t*** invsqrt_G = INVSQRT(proxy, invsqrt_data, n_gms, n_row);
//
//    double*** rec_invsqrt_G = new double**[n_gms];
//    for(int g = 0; g < n_gms; g++) {
//        rec_invsqrt_G[g] = convert2double(REC(proxy, invsqrt_G[g], n_row, n_row), n_row, n_row);
//        print2DArray("The inverse square root of the Gram matrix", rec_invsqrt_G[g], n_row, n_row);
//
//        double* straighten_invsqrt_G = new double[n_row * n_row];
//        for(uint32_t i = 0; i < n_row * n_row; i++) {
//            straighten_invsqrt_G[i] = Gs[g][i % n_row][i / n_row];
//        }
//        EigenSolver<Matrix<double, Dynamic, Dynamic>> ges;
//        Map<Matrix<double, Dynamic, Dynamic>> matrix_G(straighten_invsqrt_G, n_row, n_row);
//        ges.compute(matrix_G);
//        Matrix<double, Dynamic, Dynamic> eig_vecs = ges.eigenvectors().real();
//        Matrix<double, Dynamic, 1> eig_vals = ges.eigenvalues().real();
////        cout << eig_vals << endl;
//
//        cout << "============= GT Inverse square root of the Gram matrix ======================" << endl;
//        Matrix<double, Dynamic, Dynamic> vals = eig_vals;
//        cout << eig_vecs * vals.cwiseSqrt().cwiseInverse().asDiagonal() * Transpose(eig_vecs) << endl;
//        cout << "============================================================================" << endl;
//
//        print2DArray("The inverse of the Gram matrix",
//                     multiply_matrices(rec_invsqrt_G[g], rec_invsqrt_G[g], n_row, n_row, n_row), n_row, n_row);
//
//        cout << "============= GT Inverse of the Gram matrix ======================" << endl;
//        cout << matrix_G.inverse() << endl;
//        cout << "============================================================================" << endl;
//    }
//
//}

void DIV_Test(Party *proxy) {
    cout << setfill('*') << setw(50) << "Calling DIV";
    cout << setfill('*') << setw(49) << "*" << endl;

    uint64_t x = proxy->createShare(3);
    uint64_t y = proxy->createShare(10);

    proxy->SendBytes(CORE_DIV);
    uint64_t div = DIV(proxy, x, y);
    double reconstructed_div = convert2double(REC(proxy, div));

    // checking the result
    double originalX = convert2double(REC(proxy, x));
    double originalY = convert2double(REC(proxy, y));
    cout << "X: " << originalX << " Y: " << originalY << endl;
    double computed_div = originalX / originalY;


    if (computed_div == reconstructed_div) {
        cout << "DIV works correctly" << endl;
    } else {
        cout << "DIV works incorrectly" << endl;
        cout << "computed: " << reconstructed_div << " (" << bitset<L_BIT>(reconstructed_div) << ") but should be: "
             << computed_div << " (" << bitset<L_BIT>(computed_div) << ")" << endl;
    }

}


void ADD_Test(Party *proxy) {
    cout << setfill('*') << setw(50) << "Calling ADD functions";
    cout << setfill('*') << setw(49) << "*" << endl;

    uint64_t ***x = new uint64_t **[3];
    for (int i = 0; i < 3; ++i) {
        x[i] = proxy->createShare(random_2D_data(proxy, 5, 5, -4999.0, 5000.0), 5, 5);
    }
    uint64_t **res2 = new uint64_t *[5];
    for (int j = 0; j < 5; ++j) {
        uint64_t *combined_rows[3];
        for (int i = 0; i < 3; ++i) {
            combined_rows[i] = x[i][j];
        }
        res2[j] = ADD(proxy, combined_rows, 3, 5);
    }
    print2DArray("reconstruced result from 2D ADD", convert2double(REC(proxy, res2, 5, 5), 5, 5), 5, 5);

    uint64_t **res = ADD(proxy, x, 3, 5, 5);
    double **recon_res = convert2double(REC(proxy, res, 5, 5), 5, 5);
    print2DArray("reconstruced result from 3D ADD", recon_res, 5, 5);
    // checking the result
    double ***originalX = new double **[3];
    for (int i = 0; i < 3; ++i) {
        originalX[i] = convert2double(REC(proxy, x[i], 5, 5), 5, 5);
        print2DArray("ORIGINAL X: ", originalX[i], 5, 5);
    }

    double **correct_sum = new double *[5];
    bool allCorrect = true;
    for (uint64_t r = 0; r < 5; r++) {
        correct_sum[r] = new double[5];
        for (uint64_t c = 0; c < 5; c++) {
            double s = 0;
            for (int i = 0; i < 3; ++i) {
                s += originalX[i][r][c];
            }
            correct_sum[r][c] = s;
            if (abs(s - recon_res[r][c]) > 0.0001) {
                allCorrect = false;
            }
        }

    }
    if (allCorrect) {
        cout << "ADD works correctly" << endl;
    } else {
        cout << "ADD works incorrectly" << endl;
        print2DArray("Computed SUM: ", recon_res, 5, 5);
        print2DArray("Correct SUM: ", correct_sum, 5, 5);
    }

}

//void ppRKN_ITER_Test(Party* proxy) {
//    /*
//     * Test a single iteration of RKN which excludes the inverse square root of Gram matrix
//     */
//
//    // setup
//    int n_anc = 16; // number of anchor points
//    int n_dim = 20; // number of dimensionality of one-hot encoding
//    int k_mer = 10; // k-mer length
//    int length = 1; // length of the sequence
//
//    bool random_flag = false;
//    int size = k_mer * n_anc * n_dim;
//    int size2 = k_mer * n_anc;
//    double lambda = 0.9;
//    double alpha = 0.6;
//
//    // sequence
//    uint64_t** all_x = new uint64_t*[length];
//
//    // generate a random anchor points
//    uint64_t*** anchor_points = new uint64_t**[k_mer];
//    if(random_flag) {
//        cout << "Generate anchor points..." << endl;
//        for(int i = 0; i < k_mer; i++) {
//            anchor_points[i] = proxy->createShare(random_2D_data(proxy, n_anc, n_dim, 1, false), n_anc, n_dim);
//        }
//    }
//    else {
//        cout << "Reading anchor points..." << endl;
//        // right now, the order of the characters of the anchor points in each layer (i.e. for each k-mer) is ...
//        for(int i = 0; i < k_mer; i++) {
//            anchor_points[i] = read_2D_array(proxy, "/home/aburak/Projects/rkn_tcml/params/layer" + to_string(i) + "_k" + to_string(k_mer) +
//                                            "_anc" + to_string(n_anc) + "_dim" + to_string(n_dim), n_anc, n_dim, k_mer);
//        }
//    }
//
//    // generate a random data to represent the output of the previous time point at the same layer
//    cout << "Generate ct1..." << endl;
//    uint64_t* ct = zero_1D_data(proxy, size2 + n_anc);
//    uint64_t* initial_ct = zero_1D_data(proxy, size2 + n_anc);
//    for(int i = 0; i < n_anc; i++) {
//        ct[i] = proxy->getPRole() * ((uint64_t) 1 << FRAC);
//        initial_ct[i] = ct[i];
//    }
////    proxy->print1DArray("Initial ct", proxy->Mconvert2double(proxy->MReconstruct(ct, size2 + n_anc), size2 + n_anc), size2 + n_anc);
//
//    for(int s = 0; s < length; s++) {
//        // generate a random data
//        cout << "s: " << s << " - Generate sample data..." << endl;
//        all_x[s] = proxy->createShare(random_1D_data(proxy, n_dim, 1, false), n_dim);
//
//        // b part of the ppRKN
//        uint64_t* str_z = new uint64_t[size];
//
//        for(int i = 0; i < k_mer; i++) {
//            for(int j = 0; j < n_anc; j++) {
//                for(int k = 0; k < n_dim; k++) {
////                    cout << "check i: " << i << "\tj: " << j << "\tk: " << k << endl;
//                    str_z[(i * n_anc * n_dim) + (j * n_dim) + k] = anchor_points[i][j][k];
//                }
//            }
//        }
//        cout << "iteration " << s << endl;
//        proxy->SendBytes(RKN_ITER, size, size2);
//        print1DArray("all_x[s]", convert2double(REC(proxy, all_x[s], n_dim), n_dim), n_dim);
//        print1DArray("before ct", convert2double(REC(proxy, ct, size2), size2), size2);
//        uint64_t* tmp_ct = RKN_ITERATION(proxy, all_x[s], str_z, ct, n_dim, n_anc, k_mer, lambda, alpha);
//        copy(tmp_ct, tmp_ct + size2, ct + n_anc);
//        print1DArray("after ct", convert2double(REC(proxy, ct, size2), size2), size2);
//
//        delete [] str_z;
//    }
//    cout << "Initial mapping is done!" << endl;
////    proxy->print1DArray("c[t]", proxy->Mconvert2double(proxy->MReconstruct(ct, size2 + n_anc), size2 + n_anc), size2 + n_anc);
//
//
//    // ********************************************************************
//    // ********************************************************************
//    // ********************************************************************
//
//
//    // Ground truth computation
//    cout << "GT: reconstructing anchor points..." << endl;
//    double*** rec_anc_points = new double**[k_mer];
//    for(int i = 0; i < k_mer; i++) {
//        rec_anc_points[i] = convert2double(REC(proxy, anchor_points[i], n_anc, n_dim), n_anc, n_dim);
//    }
//
////    cout << "GT: x, t1k1 and t1k..." << endl;
//    double** rec_all_x = convert2double(REC(proxy, all_x, length, n_dim), length, n_dim);
//    double* rec_ct = convert2double(REC(proxy, initial_ct, size2 + n_anc), size2 + n_anc);
//
//    for(int iter = 0; iter < length; iter++) {
//        print1DArray("rec_all_x", rec_all_x[iter], n_dim);
//        print1DArray("before rec_ct[t]", rec_ct, size2 + n_anc);
//        // Ground truth: b part
//        // dot product
////        cout << "GT: computing dot product and exponential..." << endl;
//        double** gt_dp = new double*[k_mer];
//        double** exp_gt_dp = new double*[k_mer];
//        for(int k = 0; k < k_mer; k++) {
//            gt_dp[k] = new double[n_anc];
//            exp_gt_dp[k] = new double[n_anc];
//            for(int i = 0; i < n_anc; i++) {
//                double tmp_sum = 0;
//                for(int j = 0; j < n_dim; j++) {
//                    tmp_sum += rec_all_x[iter][j] * rec_anc_points[k][i][j];
//                }
//                gt_dp[k][i] = tmp_sum;
//                exp_gt_dp[k][i] = exp(alpha * (tmp_sum - 1));
//            }
//        }
//
//        print2DArray("exp_gt_dp " + to_string(iter), exp_gt_dp, k_mer, n_anc);
//
//        // Ground truth: c_{k-1}[t-1] * b_{l}[t]
////        cout << "GT: computing skt..." << endl;
//        double** gt_skt = new double*[k_mer];
//        for(int i = 0; i < k_mer; i++) {
//            gt_skt[i] = new double[n_anc];
//            for(int j = 0; j < n_anc; j++) {
//                gt_skt[i][j] = exp_gt_dp[i][j] * rec_ct[i * n_anc + j];
//            }
//        }
//
//        print2DArray("gt_skt " + to_string(iter), gt_skt, k_mer, n_anc);
//
//        // Ground truth: lambda * c_{k}[t-1] + s_{k}[t]
////        cout << "GT: computing ckt..." << endl;
//        double** gt_ckt = new double*[k_mer];
//        for(int i = 0; i < k_mer; i++) {
//            gt_ckt[i] = new double[n_anc];
//            for(int j = 0; j < n_anc; j++) {
//                gt_ckt[i][j] = lambda * gt_skt[i][j] + (1 - lambda) * rec_ct[(i + 1) * n_anc + j];
//            }
//        }
//
////        proxy->print2DArray("gt_ckt", gt_ckt, k_mer, n_anc);
//
//        // update c[t] based on the result of the mappings in each k-mer
//        for(int i = 1; i < k_mer + 1; i++) {
//            for(int j = 0; j < n_anc; j++) {
//                rec_ct[i * n_anc + j] = gt_ckt[i - 1][j];
//            }
//        }
//        print1DArray("after rec_ct[t]", rec_ct, size2 + n_anc);
//
//        // delete dynamically allocated arrays
//        for(int i = 0; i < k_mer; i++) {
//            delete [] gt_dp[i];
//            delete [] exp_gt_dp[i];
//            delete [] gt_skt[i];
//            delete [] gt_ckt[i];
//        }
//        delete [] gt_dp;
//        delete [] exp_gt_dp;
//        delete [] gt_skt;
//        delete [] gt_ckt;
//    }
//
//    cout << "Deleting the dynamically allocated arrays..." << endl;
//    for(int i = 0; i < length; i++) {
//        delete [] all_x[i];
//    }
//    delete [] all_x;
//
//    cout << "The computed mapping: " << endl;
//    print1DArray("c[t]", convert2double(REC(proxy, ct, size2 + n_anc), size2 + n_anc), size2 + n_anc);
//
//    cout << "Ground truth: " << endl;
//    print1DArray("GT c[t]", rec_ct, size2 + n_anc);
//}

//void ppRKN_PREDICTION_Test(Party* proxy) {
//    /*
//     * Test the whole prediction process of RKN including the inverse square root of Gram matrix
//     */
//
//    // setup
//    int n_layer = 1; // number of layers -- so far, we have only one layer
//    double reg = 0.1; // I do not remember this?
//    int n_anc = 16; // number of anchor points
//    int n_dim = 20; // number of dimensionality of one-hot encoding
//    int k_mer = 8; // k-mer length
//    double lambda = 0.5; // adjust the combination of ck[t-1] and ck[t]
//    double sigma = 0.4; // implicitly used in similarity computation
//    double alpha = 1.0 / (pow(sigma, 2) * k_mer);
//    string pooling = "gmp"; // pooling -- which is canceled and has no effect
//    string tfid = "a.101.1"; // sample id
//    string enc = "one_hot"; // encoding type
//    string eps = "_eps"; // do not remember?
//    int s_ind = 1; // test sample index
//    bool random_flag = false; // whether to use random values or a real example
//    double epsilon = 0.01; // epsilon added on top of eigenvalues for numeric problems - to replicate RKN
//
//    ostringstream oss;
//    oss << setprecision(1) << noshowpoint << lambda;
//    std::string str_lmb = oss.str();
//    ostringstream oss2;
//    oss2 << setprecision(1) << noshowpoint << sigma;
//    std::string str_sigma = oss2.str();
//    ostringstream oss3;
//    oss3 << setprecision(1) << noshowpoint << reg;
//    std::string str_reg = oss3.str();
//
//    int length;
//    uint64_t** all_x;
//    uint64_t*** anchor_points = new uint64_t**[k_mer];
//    uint64_t*** tr_anchor_points = new uint64_t**[k_mer]; // transpose of the anchor points in each layer
//    uint64_t* weights;
//    if(random_flag) { // random values
//        length = 20; // length of the synthetic sequence
//        all_x = new uint64_t*[length];
//        cout << "Generating data..." << endl;
//        for(int s = 0; s < length; s++) {
//            all_x[s] = proxy->createShare(random_1D_data(proxy, n_dim, 1, false), n_dim);
//        }
//
//        // generate a random anchor points
//        cout << "Generating anchor points..." << endl;
//        for(int i = 0; i < k_mer; i++) {
//            anchor_points[i] = proxy->createShare(random_2D_data(proxy, n_anc, n_dim, 1, false), n_anc, n_dim);
//
//            tr_anchor_points[i] = new uint64_t*[n_dim];
//            for(int r = 0; r < n_dim; r++) {
//                tr_anchor_points[i][r] = new uint64_t[n_anc];
//                for(int c = 0; c < n_anc; c++) {
//                    tr_anchor_points[i][r][c] = anchor_points[i][c][r];
//                }
//            }
//        }
//
//        // linear layer for the classification
//        weights = proxy->createShare(random_1D_data(proxy, n_anc + 1, 0.0, 1.0), n_anc + 1);
////        print1DArray("Weights", convert2double(REC(proxy, weights, n_anc), n_anc), n_anc);
//    }
//    else { // real values
//        // sequence
//        string folder_name = to_string(n_layer) + "_[" + to_string(n_anc) + "]_[" + to_string(k_mer) + "]_[" +
//                             str_lmb + "]_[" + str_sigma + "]_" + str_reg;
//        string base_fn = "/home/aburak/Projects/Framework/rkn_results/" +  pooling + "/" + enc + "/" + folder_name + "/" + tfid;
//        cout << "Base folder name: " << base_fn << endl;
//        string seq = recover_seq(base_fn + "/test_samples.csv", s_ind);
//        length = seq.length(); // length of the sequence
//        cout << "Sequence with length " << length << " :" << endl;
//        for(int i = 0; i < seq.length(); i++) {
//            cout << seq[i];
//        }
//        cout << endl;
//
//        all_x = encode_sequence(proxy, seq);
//
//        cout << "Reading anchor points..." << endl;
//        for(int i = 0; i < k_mer; i++) {
//            anchor_points[i] = read_2D_array(proxy, base_fn + "/layer" + to_string(i) + "_k" + to_string(k_mer) +
//                                    "_anc" + to_string(n_anc) + "_dim" + to_string(n_dim), n_anc, n_dim, k_mer);
//
//            tr_anchor_points[i] = new uint64_t*[n_dim];
//            for(int r = 0; r < n_dim; r++) {
//                tr_anchor_points[i][r] = new uint64_t[n_anc];
//                for(int c = 0; c < n_anc; c++) {
//                    tr_anchor_points[i][r][c] = anchor_points[i][c][r];
//                }
//            }
//        }
//
//        // linear layer for the classification
//        weights = read_1D_array(proxy, base_fn + "/linear_layer_k" + to_string(k_mer) + "_anc" + to_string(n_anc) +
//                                "_dim" + to_string(n_dim), n_anc + 1);
////        print1DArray("Weights", convert2double(REC(proxy, weights, n_anc), n_anc), n_anc);
//    }
//
//    int size = k_mer * n_anc * n_dim;
//    int size2 = k_mer * n_anc;
//
////    proxy->SendBytes(RKN_PRE);
//
//    // generate a random data to represent the output of the previous time point at the same layer
//    cout << "Generate ct1..." << endl;
//    uint64_t* ct = zero_1D_data(proxy, size2 + n_anc);
//    uint64_t* initial_ct = zero_1D_data(proxy, size2 + n_anc);
//    for(int i = 0; i < n_anc; i++) {
//        ct[i] = proxy->getPRole() * ((uint64_t) 1 << FRAC);
//        initial_ct[i] = ct[i];
//    }
//
//    // generate random sequence data
//    for(int s = 0; s < length; s++) {
//        // b part of the ppRKN
//        uint64_t* str_z = new uint64_t[size];
//
//        for(int i = 0; i < k_mer; i++) {
//            for(int j = 0; j < n_anc; j++) {
//                for(int k = 0; k < n_dim; k++) {
////                    cout << "check i: " << i << "\tj: " << j << "\tk: " << k << endl;
//                    str_z[(i * n_anc * n_dim) + (j * n_dim) + k] = anchor_points[i][j][k];
//                }
//            }
//        }
//        cout << "iteration " << s << endl;
//        proxy->SendBytes(RKN_ITER, size, size2);
////        print1DArray("all_x[s]", convert2double(REC(proxy, all_x[s], n_dim), n_dim), n_dim);
////        print1DArray("before ct", convert2double(REC(proxy, ct, size2), size2), size2);
//        uint64_t* tmp_ct = RKN_ITERATION(proxy, all_x[s], str_z, ct, n_dim, n_anc, k_mer, lambda, alpha);
//        copy(tmp_ct, tmp_ct + size2, ct + n_anc);
////        print1DArray("after ct", convert2double(REC(proxy, ct, size2), size2), size2);
//
//        uint64_t** mat_ct = new uint64_t *[k_mer];
//        for(int i = 0; i < k_mer; i++) {
//            mat_ct[i] = new uint64_t[n_anc];
//            for(int j = 0; j < n_anc; j++) {
//                mat_ct[i][j] = ct[n_anc + (i * n_anc) + j];
//            }
//        }
////        print2DArray("Mappings at " + to_string(s), convert2double(REC(proxy, mat_ct, k_mer, n_anc), k_mer, n_anc),
////                     k_mer, n_anc, false);
//
//        delete [] str_z;
//    }
//
//    cout << "Initial mapping is done!" << endl;
////    proxy->print1DArray("c[t]", proxy->Mconvert2double(proxy->MReconstruct(ct, size2 + n_anc), size2 + n_anc), size2 + n_anc);
//
//    // convert c[t] to matrix
//    uint64_t** mat_ct = new uint64_t *[k_mer];
//    for(int i = 0; i < k_mer; i++) {
//        mat_ct[i] = new uint64_t[n_anc];
//        for(int j = 0; j < n_anc; j++) {
//            mat_ct[i][j] = ct[n_anc + (i * n_anc) + j];
//        }
//    }
//
//    // Gram matrices of the anchor points
//    proxy->SendBytes(CORE_MMATMATMUL, k_mer * n_anc * n_dim * n_anc);
//    uint64_t*** gms = MATMATMUL(proxy, anchor_points, tr_anchor_points, k_mer, n_anc, n_dim, n_anc);
//
//    proxy->SendBytes(RKN_GM2KM, k_mer, n_anc);
//    uint64_t*** kmer_kms = GM2KM(proxy, gms, convert2uint64(alpha), k_mer, n_anc);
//
//    double*** rec_kmer_kms = new double**[k_mer];
//    for(int g = 0; g < k_mer; g++) {
//        rec_kmer_kms[g] = convert2double(REC(proxy, kmer_kms[g], n_anc, n_anc), n_anc, n_anc);
//    }
//
//    proxy->SendBytes(RKN_MINVSQRT, k_mer, n_anc);
//    uint64_t*** invsqrt_gms = INVSQRT(proxy, kmer_kms, k_mer, n_anc, epsilon);
//
//    proxy->SendBytes(CORE_MMATVECMUL, k_mer * n_anc * n_anc);
//    uint64_t** x_mapping = MATVECMUL(proxy, invsqrt_gms, mat_ct, k_mer, n_anc, n_anc);
//
//    double** rec_x_mapping = convert2double(REC(proxy, x_mapping, k_mer, n_anc), k_mer, n_anc);
//
//    proxy->SendBytes(CORE_DP, n_anc);
//    uint64_t prediction = DP(proxy, weights, x_mapping[k_mer - 1], n_anc);
//
//    for(int i = 0; i < k_mer; i++) {
//        delete [] mat_ct[i];
//        delete [] x_mapping[i];
//        for(int j = 0; j < n_anc; j++) {
//            delete [] gms[i][j];
//            if(i != 0)
//            delete [] kmer_kms[i][j];
//            delete [] invsqrt_gms[i][j];
//        }
//        delete [] gms[i];
//        if(i != 0)
//        delete [] kmer_kms[i];
//        delete [] invsqrt_gms[i];
//    }
//    delete [] mat_ct;
//    delete [] x_mapping;
//    delete [] gms;
//    delete [] kmer_kms;
//    delete [] invsqrt_gms;
//
//
//
//    // ********************************************************************
//    // ********************************************************************
//    // ********************************************************************
//
//
//    // Ground truth computation
//    cout << "Ground truth computation starts..." << endl;
//    double*** rec_anc_points = new double**[k_mer];
//    for(int i = 0; i < k_mer; i++) {
//        rec_anc_points[i] = convert2double(REC(proxy, anchor_points[i], n_anc, n_dim), n_anc, n_dim);
//    }
//
//    double** rec_all_x = convert2double(REC(proxy, all_x, length, n_dim), length, n_dim);
//    double* rec_ct = convert2double(REC(proxy, initial_ct, size2 + n_anc), size2 + n_anc);
//
//    double** gt_dp = new double*[k_mer];
//    double** exp_gt_dp = new double*[k_mer];
//    double** gt_skt = new double*[k_mer];
//    double** gt_ckt = new double*[k_mer];
//    for(int k = 0; k < k_mer; k++) {
//        gt_dp[k] = new double[n_anc];
//        exp_gt_dp[k] = new double[n_anc];
//        gt_skt[k] = new double[n_anc];
//        gt_ckt[k] = new double[n_anc];
//    }
//
//    for(int iter = 0; iter < length; iter++) {
//        // Ground truth: b part
//        // dot product
//        for(int k = 0; k < k_mer; k++) {
//            for(int i = 0; i < n_anc; i++) {
//                double tmp_sum = 0;
//                for(int j = 0; j < n_dim; j++) {
//                    tmp_sum += rec_all_x[iter][j] * rec_anc_points[k][i][j];
//                }
//                gt_dp[k][i] = tmp_sum;
//                exp_gt_dp[k][i] = exp(alpha * (tmp_sum - 1));
//            }
//        }
//
//        // Ground truth: c_{k-1}[t-1] * b_{l}[t]
//        for(int i = 0; i < k_mer; i++) {
//            for(int j = 0; j < n_anc; j++) {
//                gt_skt[i][j] = exp_gt_dp[i][j] * rec_ct[i * n_anc + j];
//            }
//        }
//
//        // Ground truth: lambda * c_{k}[t-1] + s_{k}[t]
//        for(int i = 0; i < k_mer; i++) {
//            for(int j = 0; j < n_anc; j++) {
//                gt_ckt[i][j] = lambda * gt_skt[i][j] + (1 - lambda) * rec_ct[(i + 1) * n_anc + j];
//            }
//        }
//
//        // update c[t] based on the result of the mappings in each k-mer
//        for(int i = 1; i < k_mer + 1; i++) {
//            for(int j = 0; j < n_anc; j++) {
//                rec_ct[i * n_anc + j] = gt_ckt[i - 1][j];
//            }
//        }
//    }
//
//    // delete dynamically allocated arrays
//    for(int i = 0; i < k_mer; i++) {
//        delete [] gt_dp[i];
//        delete [] exp_gt_dp[i];
//        delete [] gt_skt[i];
//        delete [] gt_ckt[i];
//    }
//    delete [] gt_dp;
//    delete [] exp_gt_dp;
//    delete [] gt_skt;
//    delete [] gt_ckt;
//
//    // ----------------------------------------------------------------------------------------------
//    // generate Gram matrices
//    double*** gt_gms = new double**[k_mer];
//    gt_gms[0] = inplace_dp(rec_anc_points[0], rec_anc_points[0], n_anc, n_dim);
//    for(int j = 0; j < n_anc; j++) {
//        for(int k = j; k < n_anc; k++) {
//            gt_gms[0][j][k] = exp(alpha * (gt_gms[0][j][k] - 1));
//            gt_gms[0][k][j] = gt_gms[0][j][k];
//        }
//    }
//
//    // initialize the rest of the gt_gms array
//    for(int i = 1; i < k_mer; i++) {
//        gt_gms[i] = new double*[n_anc];
//        for(int j = 0; j < n_anc; j++) {
//            gt_gms[i][j] = new double[n_anc];
//        }
//    }
//
//    //    proxy->print2DArray("GT Gram matrix 0", gt_gms[0], n_anc, n_anc);
//    for(int i = 1; i < k_mer; i++) {
//        double** tmp_gt_gms = inplace_dp(rec_anc_points[i], rec_anc_points[i], n_anc, n_dim);
//        for(int j = 0; j < n_anc; j++) {
//            for(int k = j; k < n_anc; k++) {
//                gt_gms[i][j][k] = exp(alpha * (tmp_gt_gms[j][k] - 1)) * gt_gms[i - 1][j][k];
//                gt_gms[i][k][j] = gt_gms[i][j][k];
//            }
//        }
//    }
//
//    // compute kernel matrices
//    double*** gt_kms = new double**[k_mer];
//    for(int i = 0; i < k_mer; i++) {
//        gt_kms[i] = new double*[n_anc];
//        for(int j = 0; j < n_anc; j++) {
//            gt_kms[i][j] = new double[n_anc];
//            for(int k = 0; k < n_anc; k++) {
//                gt_kms[i][j][k] = gt_gms[i][j][k];
//            }
//        }
//    }
//
//    double** gt_res = new double*[k_mer];
//    double** gt_eigvals = new double*[k_mer];
//    double** AT_gt_eigvals = new double*[k_mer];
//    for(int g = 0; g < k_mer; g++) {
//        double* straighten_G = new double[n_anc * n_anc];
//        double* AT_straighten_G = new double[n_anc * n_anc];
//        for(uint32_t i = 0; i < n_anc * n_anc; i++) {
//            straighten_G[i] = gt_kms[g][i % n_anc][i / n_anc];
//            AT_straighten_G[i] = rec_kmer_kms[g][i % n_anc][i / n_anc];
//        }
//
//        // ****************************************************************************************************
//        EigenSolver<Matrix<double, Dynamic, Dynamic, RowMajor>> AT_ges;
//        Map<Matrix<double, Dynamic, Dynamic, RowMajor>> AT_matrix_G(AT_straighten_G, n_anc, n_anc);
//        AT_ges.compute(AT_matrix_G);
//        Matrix<double, Dynamic, 1> AT_eig_vals = AT_ges.eigenvalues().real();
//        AT_gt_eigvals[g] = new double[n_anc];
//        Map<Matrix<double, Dynamic, 1>>(AT_gt_eigvals[g], n_anc) = AT_eig_vals;
//        // ****************************************************************************************************
//
//        EigenSolver<Matrix<double, Dynamic, Dynamic, RowMajor>> ges;
//        Map<Matrix<double, Dynamic, Dynamic, RowMajor>> matrix_G(straighten_G, n_anc, n_anc);
//        ges.compute(matrix_G);
//        Matrix<double, Dynamic, Dynamic, RowMajor> eig_vecs = ges.eigenvectors().real();
//        Matrix<double, Dynamic, 1> eig_vals = ges.eigenvalues().real();
//
//        gt_eigvals[g] = new double[n_anc];
//        Map<Matrix<double, Dynamic, 1>>(gt_eigvals[g], n_anc) = eig_vals;
//
//        //        cout << "********************************************\nGT eigenvalues of gram matrix " << g << ":\n" << eig_vals << endl;
//
//        Matrix<double, Dynamic, Dynamic, RowMajor> vals = eig_vals;
//
//        //        cout << "GT reconstructed inverse square root of the Gram matrix " << g << ":\n" <<
//        //        eig_vecs * vals.cwiseSqrt().cwiseInverse().asDiagonal() * Transpose(eig_vecs) << endl;
//
//        double* tmp_str_invsqrt = new double[n_anc * n_anc];
//        Map<Matrix<double, Dynamic, Dynamic, RowMajor>>(tmp_str_invsqrt, n_anc, n_anc) =
//                eig_vecs * (vals.cwiseSqrt().array() + epsilon).matrix().cwiseInverse().asDiagonal() * Transpose(eig_vecs);
//        double** tmp_invsqrt_gm = new double*[n_anc];
//        for(int at = 0; at < n_anc; at++) {
//            tmp_invsqrt_gm[at] = new double[n_anc];
//            for(int kafa = 0; kafa < n_anc; kafa++) {
//                tmp_invsqrt_gm[at][kafa] = tmp_str_invsqrt[at * n_anc + kafa];
//            }
//        }
//
//        gt_res[g] = multiply_matrice_vector(tmp_invsqrt_gm, &rec_ct[(g + 1) * n_anc], n_anc, n_anc);
//
//        // deleting dynamically allocated arrays
//        delete [] straighten_G;
//        delete [] tmp_str_invsqrt;
//        for(int d = 0; d < n_anc; d++) {
//            delete [] tmp_invsqrt_gm[d];
//        }
//        delete [] tmp_invsqrt_gm;
//    }
//
//    double* rec_weights = convert2double(REC(proxy, weights, n_anc), n_anc);
//    double gt_prediction = multiply_vector_vector(gt_res[k_mer - 1], rec_weights, n_anc);
//
//    double* total_diff = new double[k_mer];
//    for(int i = 0; i < k_mer; i++) {
//        total_diff[i] = 0;
//    }
//
//    double **diff = new double*[n_anc];
//    for(int i = 0; i < n_anc; i++) {
//        diff[i] = new double[k_mer];
//        for(int j = 0; j < k_mer; j++) {
//            diff[i][j] = gt_res[j][i] - rec_x_mapping[j][i];
//            total_diff[j] += abs(diff[i][j]);
//        }
//    }
//
//    print2DArray("Differences between mappings", diff, n_anc, k_mer, true);
//    print1DArray("Total differences between mappings", total_diff, k_mer);
//
//    printValue("Prediction", convert2double(REC(proxy, prediction)));
//    printValue("GT Prediction", gt_prediction);
//    printValue("|Prediction - GT Prediction|", abs(convert2double(REC(proxy, prediction)) - gt_prediction));
//
////    MbubbleSort(gt_eigvals, k_mer, n_anc);
////    MbubbleSort(AT_gt_eigvals, k_mer, n_anc);
//
////    proxy->print2DArray("GT Eigenvalues", gt_eigvals, k_mer, n_anc, false);
////    proxy->print2DArray("AT GT Eigenvalues", AT_gt_eigvals, k_mer, n_anc, false);
//
//    for(int i = 0; i < k_mer; i++) {
//        for(int j = 0; j < n_anc; j++) {
//            delete [] rec_anc_points[i][j];
//            delete [] anchor_points[i][j];
//            delete [] gt_gms[i][j];
//            delete [] gt_kms[i][j];
//        }
//        delete [] rec_anc_points[i];
//        delete [] anchor_points[i];
//        delete [] gt_gms[i];
//        delete [] gt_kms[i];
//    }
//    delete [] rec_anc_points;
//    delete [] anchor_points;
//    delete [] gt_gms;
//    delete [] gt_kms;
//    delete [] rec_ct;
//
//    // ---------------------------------
//    for(int i = 0; i < length; i++) {
//        delete [] rec_all_x[i];
//        delete [] all_x[i];
//    }
//    delete [] rec_all_x;
//    delete [] all_x;
//
//    for(int i = 0; i < n_anc; i++) {
//        delete [] diff[i];
//    }
//    delete [] diff;
//}



bool NETWORK_TEST(Party *proxy) {
    cout << setfill('*') << setw(50) << "Calling LeNet Test";
    cout << setfill('*') << setw(49) << "*" << endl;

    const int length = 28;
    double image[length][length] = {{0, 0, 0, 0, 0, 0, 0, 0,  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0},
                                    {0, 0, 0, 0, 0, 0, 0, 0,  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0},
                                    {0, 0, 0, 0, 0, 0, 0, 0,  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0},
                                    {0, 0, 0, 0, 0, 0, 0, 0,  0,   0,   116, 125, 171, 255, 255, 150, 93,  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0},
                                    {0, 0, 0, 0, 0, 0, 0, 0,  0,   169, 253, 253, 253, 253, 253, 253, 218, 30,  0,   0,   0,   0,   0,   0,   0,   0,   0,   0},
                                    {0, 0, 0, 0, 0, 0, 0, 0,  169, 253, 253, 253, 213, 142, 176, 253, 253, 122, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0},
                                    {0, 0, 0, 0, 0, 0, 0, 52, 250, 253, 210, 32,  12,  0,   6,   206, 253, 140, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0},
                                    {0, 0, 0, 0, 0, 0, 0, 77, 251, 210, 25,  0,   0,   0,   122, 248, 253, 65,  0,   0,   0,   0,   0,   0,   0,   0,   0,   0},
                                    {0, 0, 0, 0, 0, 0, 0, 0,  31,  18,  0,   0,   0,   0,   209, 253, 253, 65,  0,   0,   0,   0,   0,   0,   0,   0,   0,   0},
                                    {0, 0, 0, 0, 0, 0, 0, 0,  0,   0,   0,   0,   0,   117, 247, 253, 198, 10,  0,   0,   0,   0,   0,   0,   0,   0,   0,   0},
                                    {0, 0, 0, 0, 0, 0, 0, 0,  0,   0,   0,   0,   76,  247, 253, 231, 63,  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0},
                                    {0, 0, 0, 0, 0, 0, 0, 0,  0,   0,   0,   0,   128, 253, 253, 144, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0},
                                    {0, 0, 0, 0, 0, 0, 0, 0,  0,   0,   0,   176, 246, 253, 159, 12,  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0},
                                    {0, 0, 0, 0, 0, 0, 0, 0,  0,   0,   25,  234, 253, 233, 35,  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0},
                                    {0, 0, 0, 0, 0, 0, 0, 0,  0,   0,   198, 253, 253, 141, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0},
                                    {0, 0, 0, 0, 0, 0, 0, 0,  0,   78,  248, 253, 189, 12,  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0},
                                    {0, 0, 0, 0, 0, 0, 0, 0,  19,  200, 253, 253, 141, 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0},
                                    {0, 0, 0, 0, 0, 0, 0, 0,  134, 253, 253, 173, 12,  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0},
                                    {0, 0, 0, 0, 0, 0, 0, 0,  248, 253, 253, 25,  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0},
                                    {0, 0, 0, 0, 0, 0, 0, 0,  248, 253, 253, 43,  20,  20,  20,  20,  5,   0,   5,   20,  20,  37,  150, 150, 150, 147, 10,  0},
                                    {0, 0, 0, 0, 0, 0, 0, 0,  248, 253, 253, 253, 253, 253, 253, 253, 168, 143, 166, 253, 253, 253, 253, 253, 253, 253, 123, 0},
                                    {0, 0, 0, 0, 0, 0, 0, 0,  174, 253, 253, 253, 253, 253, 253, 253, 253, 253, 253, 253, 249, 247, 247, 169, 117, 117, 57,  0},
                                    {0, 0, 0, 0, 0, 0, 0, 0,  0,   118, 123, 123, 123, 166, 253, 253, 253, 155, 123, 123, 41,  0,   0,   0,   0,   0,   0,   0},
                                    {0, 0, 0, 0, 0, 0, 0, 0,  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0},
                                    {0, 0, 0, 0, 0, 0, 0, 0,  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0},
                                    {0, 0, 0, 0, 0, 0, 0, 0,  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0},
                                    {0, 0, 0, 0, 0, 0, 0, 0,  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0},
                                    {0, 0, 0, 0, 0, 0, 0, 0,  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0}};

    uint64_t ***secret = new uint64_t **[1];
    // image in shape: channel x rows x columns
    secret[0] = new uint64_t *[length];
    for (int i = 0; i < length; ++i) {
        secret[0][i] = proxy->createShare(image[i], length);
    }

    double kernel[5][5] = {{-0.09602198,  -0.066380806, 0.004841719,    -0.13074008, -0.21051669},
                           {-0.15994059,  0.061118018,  0.19219273,     0.06320523,  0.15816787},
                           {0.0129997665, 0.19746655,   -0.00050751516, -0.17033894, -0.070527315},
                           {0.054051124,  -0.20655234,  -0.1440479,     -0.21508642, 0.21944945},
                           {0.16889668,   -0.062005207, 0.14000067,     0.19843176,  0.11537137}};
    // kernel in shape: channel x num_of_kernel x length_of_kernel
    uint64_t ***sec_kernel = new uint64_t **[1];
    sec_kernel[0] = new uint64_t *[2];
    sec_kernel[0][0] = new uint64_t[25];
    sec_kernel[0][1] = new uint64_t[25];
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < 5; ++j) {
            sec_kernel[0][0][i * 5 + j] = proxy->createShare(kernel[i][j]);
            sec_kernel[0][1][i * 5 + j] = proxy->createShare(kernel[i][(j + 1) % 5]);
        }
    }

    double bias = 0.1175425;
    uint64_t *sec_bias = new uint64_t[0];
    sec_bias[0] = proxy->createShare(bias);

    uint64_t max_width = floor((length - 5 + 1) / 2);
    uint64_t max_height = floor((length - 5 + 1) / 2);

    // send params
    uint32_t params[9];
    params[0] = 1;
    params[1] = length;
    params[2] = length;
    params[3] = 5;      // kernel size
    params[4] = 2;   // kernel number = output channel
    params[5] = 1;
    params[6] = 2;
    params[7] = 2;
    params[8] = false;
    proxy->SendBytes(CNN_CL, params, 9);

    print2DArray("image: ", convert2double(REC(proxy, secret[0], length, length), length, length), length, length);
    uint64_t ***weights_secure = CL(proxy, secret, 1, length, length, sec_kernel, 5, 2, 1, 2, 2, sec_bias, false);
    // checking the result
    double **recon_res = convert2double(REC(proxy, weights_secure[0], max_height, max_width), max_height, max_width);

    uint64_t conv_size = length - 5 + 1;
    double **correct_res = new double *[conv_size];
    bool allCorrect = true;
    for (uint32_t cr = 0; cr < conv_size; cr++) {
        correct_res[cr] = new double[conv_size];   // init row of conv result
        for (uint32_t cc = 0; cc < conv_size; cc++) {
            double dot_product = 0;
            for (uint32_t kr = 0; kr < 5; kr++) {
                for (int kc = 0; kc < 5; ++kc) {
                    double v = image[cr + kr][cc + kc];
                    double weight = kernel[kr][kc];
                    dot_product += v * weight;
                }
            }
            dot_product += bias;
            // Activation: RELU
            double relu = 0.0;
            if (dot_product > 0) {
                relu = dot_product;
            }
            correct_res[cr][cc] = relu;
        }
    }
    //print2DArray("activated conv", correct_res, conv_size, conv_size);
    //MAXPOOLING
    double **pooled_conv = new double *[max_height];

    for (int r = 0; r < (conv_size - 2 + 1); r += 2) {
        pooled_conv[r / 2] = new double[max_width];
        for (int c = 0; c < (conv_size - 2 + 1); c += 2) {
            //find max of window:
            double max = correct_res[r][c]; // first value in window
            for (int max_r = 0; max_r < 2; ++max_r) {
                for (int max_c = 0; max_c < 2; ++max_c) {
                    double next_value = correct_res[r + max_r][c + max_c];
                    if (next_value > max) {
                        max = next_value;
                    }
                }
            }
            pooled_conv[r / 2][c / 2] = max;
            if (abs(max - recon_res[r / 2][c / 2]) > 0.1) {
                cout << r / 2 << " " << c / 2 << ": " << max << " (computed: " << recon_res[r / 2][c / 2] << ")"
                     << endl;
                allCorrect = false;
            }
        }
    }


    for (uint32_t cr = 0; cr < conv_size; cr++) {
        correct_res[cr] = new double[conv_size];   // init row of conv result
        for (uint32_t cc = 0; cc < conv_size; cc++) {
            double dot_product = 0;
            for (uint32_t kr = 0; kr < 5; kr++) {
                for (int kc = 0; kc < 5; ++kc) {
                    double v = image[cr + kr][cc + kc];
                    double weight = kernel[kr][kc];
                    dot_product += v * weight;
                }
            }
            dot_product += bias;
            // Activation: RELU
            double relu = 0.0;
            if (dot_product > 0) {
                relu = dot_product;
            }
            correct_res[cr][cc] = relu;
        }
    }
    //print2DArray("activated conv", correct_res, conv_size, conv_size);
    //MAXPOOLING
    pooled_conv = new double *[max_height];

    for (int r = 0; r < (conv_size - 2 + 1); r += 2) {
        pooled_conv[r / 2] = new double[max_width];
        for (int c = 0; c < (conv_size - 2 + 1); c += 2) {
            //find max of window:
            double max = correct_res[r][c]; // first value in window
            for (int max_r = 0; max_r < 2; ++max_r) {
                for (int max_c = 0; max_c < 2; ++max_c) {
                    double next_value = correct_res[r + max_r][c + max_c];
                    if (next_value > max) {
                        max = next_value;
                    }
                }
            }
            pooled_conv[r / 2][c / 2] = max;
            if (abs(max - recon_res[r / 2][c / 2]) > 0.1) {
                cout << r / 2 << " " << c / 2 << ": " << max << " (computed: " << recon_res[r / 2][c / 2] << ")"
                     << endl;
                allCorrect = false;
            }
        }
    }

    if (allCorrect) {
        cout << "networks first conv works correctly" << endl;
    } else {
        cout << "conv works incorrectly" << endl;
        print2DArray("Computed first conv: ", recon_res, max_height, max_width);
        print2DArray("Correct first conv: ", pooled_conv, max_height, max_width);
    }
    return allCorrect;
}

bool NETWORK_M_INPUTS_TEST(Party *proxy) {

    const int rows = 12;
    const int cols = 12;
    const int channels = 2;
    double input[channels][rows][cols] = {{{0,        0,        0,        0,        0,        0,        0,        0,        0,        0,        0,        0},
                                                  {0,        5.44205, 0,       0,       0,        3.11207,  0,        0,        0,        0,        0,        0},
                                                  {0,        4.51266, 0,       0,       0, 0, 0, 0, 0, 3.28263, 5.4091, 0},
                                                  {0,        118.905, 228.533, 292.015, 284.143, 208.377, 150.79,  142.398, 10.2446, 0,      0, 0},
                                                  {0,        58.1105, 124.143, 217.863, 259.282, 290.512, 290.146, 265.871, 205.192, 179.348, 42.3026, 0},
                                                  {0,        0,        0,        0,        9.73383, 18.3992, 27.8476, 34.2367, 115.436, 37.1317, 54.0832, 0},
                                                  {0,        0,        0,        0,        0,        7.12205,  9.38115, 119.877, 114.577, 71.8713, 18.4409,  0},
                                                  {0,        0,        0,        0,        1.1916,   4.85352, 34.714,  125.124, 49.0738, 56.4896, 0,        0},
                                                  {0,        0,        0,        0,        10.2245, 0,       101.164, 80.6467, 53.934, 8.44365,  0,        0},
                                                  {0,        0,        0,        4.97295,  3.08055, 113.148, 128.491, 71.8604, 55.6757,  0,        0,        0},
                                                  {0,        0,        0,        3.93804, 5.94114, 52.9342, 57.743,  76.2328, 5.07514,  0,        0,        0},
                                                  {0,        0,        0,        67.2699, 200.28,  158.916, 39.3759, 5.98702,  0,        0,        0,        0}},
                                          {{0.179874, 0.179874, 0.179874, 0.179874, 0.179874, 0.179874, 0.179874, 0.179874, 0.179874, 0.179874, 0.179874, 0.179874},
                                                  {0.179874, 10.1227, 39.0021, 15.1475, 0.179874, 0.179874, 0.179874, 0.179874, 0.179874, 0.179874, 0.179874, 0.179874},
                                                  {0.179874, 84.6016, 95.2082, 11.5725, 0, 0, 0, 0, 0, 0,       0,      0.179874},
                                                  {0.179874, 5.56343, 170.922, 92.3377, 51.0432, 31.77,   50.4084, 78.5217, 130.011, 30.587, 0, 0.179874},
                                                  {0.179874, 20.4813, 129.587, 63.875,  109.919, 82.8212, 87.0372, 140.246, 189.608, 21.0135, 0,       0.179874},
                                                  {0.179874, 0.179874, 0.179874, 0.179874, 17.1033, 10.8175, 53.67,   85.0075, 168.667, 0,       0,       0.179874},
                                                  {0.179874, 0.179874, 0.179874, 0.179874, 0.179874, 0.179874, 57.3658, 201.08,  167.161, 0,       0.179874, 0.179874},
                                                  {0.179874, 0.179874, 0.179874, 0.179874, 0.179874, 77.0553, 84.4292, 160.509, 0,       0,       0.179874, 0.179874},
                                                  {0.179874, 0.179874, 0.179874, 0.179874, 41.4585, 74.3464, 101.913, 107.129, 0,      0.179874, 0.179874, 0.179874},
                                                  {0.179874, 0.179874, 0.179874, 0.678918, 65.4045, 143.124, 129.084, 0,       0.179874, 0.179874, 0.179874, 0.179874},
                                                  {0.179874, 0.179874, 0.179874, 47.9556, 125.135, 213.215, 36.3252, 0,       0.179874, 0.179874, 0.179874, 0.179874},
                                                  {0.179874, 0.179874, 0.179874, 24.3328, 225.504, 217.689, 0,       0.179874, 0.179874, 0.179874, 0.179874, 0.179874}}};
    uint64_t ***secret = new uint64_t **[channels];
    for (int i = 0; i < channels; ++i) {
        secret[i] = new uint64_t *[rows];
        for (int j = 0; j < rows; ++j) {
            secret[i][j] = proxy->createShare(input[i][j], cols);
        }
    }
    const int n_kernel = 3;
    const int k_size = 5;
    double kernel[channels][n_kernel][k_size*k_size] = {{{0.0152697, 0.054374,  -0.0112883, 0.0392754, -0.0214642, -0.0312186, 0.0307207,  -0.00210262, 0.0411354, 0.037925,   -0.0441519, -0.0108989, -0.0450431, 0.018118,  0.020019,  0.0182398,  -0.0236182, -0.00400799, 0.0128652, -0.00035863, 0.0423267,  0.0190162,  -0.0499434, 0.0372975, -0.0153789},
                                               {-0.0120976, -0.0427761, 0.036385,   -0.0153403, -0.0346747,  -0.0203395, -0.0154173, 0.031334, 0.0262655,  0.00224767, 0.0172362, 0.0265062, -0.038035, 0.0301112,  -0.0243593,  0.0414339, 0.0451315,    0.0960202,  0.0864642,  0.0673048,   0.0375359,   0.0944946, 0.036669,  0.0578383, 0.0477528},
                                               {0.011804,  -0.0108315, -0.0148268, 0.0290907, -0.0321819, 0.0289785,  0.016873,   -0.0234544, 0.00653713, 0.0468413,   -0.00410596, 0.0145675,   0.00133672, 0.0338526, -0.00800818, 0.0219828,  0.0272399,  -0.037904,  0.000337112, -0.0222561, 0.0111167,  -0.02895,  0.0180652, 0.0372134,   -0.0268938}},
                                               {{0.0111452, 0.0214061, 0.0687376,  0.062753,  -0.0097982, -0.0045877, 0.00059125, 0.0382741,   0.0692399, -0.0240434, 0.00947216, 0.0385408,  0.0435236,  0.0615116, 0.0416076, -0.0390097, -0.0417169, 0.0151983,   0.0608362, 0.00883882,  -0.0401907, -0.0230137, 0.049076,   0.0608768, 0.0504488},
                                               {-0.0262289, -0.0297652, -0.0142716, 0.00283265, -0.00847504, -0.0262626, -0.0416967, -0.02454, -0.0123295, -0.0617562, 0.0243651, 0.0153887, 0.0177689, -0.0422747, -0.00707854, 0.0623879, -0.000610727, 0.00134396, -0.0315004, -0.00346995, 0.000966561, 0.0323579, 0.0769279, 0.0570635, 0.0289574},
                                               {0.0454316, -0.0314481, 0.00686538, 0.0484373, -0.0315532, -0.0090278, 0.00539585, 0.00143481, -0.0152588, -0.00149502, -0.00734612, -0.00654482, 0.0232091,  0.0353256, -0.0278187,  -0.0390783, 0.00789026, 0.00284264, -0.0160283,  -0.0279816, -0.0237621, 0.0404292, 0.0598701, 0.000395513, 0.000991492}}};
    /*double kernel[channels][n_kernel][k_size*k_size] = {{{0.0152697, 0.054374,  -0.0112883, 0.0392754, -0.0214642, -0.0312186, 0.0307207,  -0.00210262, 0.0411354, 0.037925,   -0.0441519, -0.0108989, -0.0450431, 0.018118,  0.020019,  0.0182398,  -0.0236182, -0.00400799, 0.0128652, -0.00035863, 0.0423267,  0.0190162,  -0.0499434, 0.0372975, -0.0153789},
                                                                {-0.0120976, -0.0427761, 0.036385,   -0.0153403, -0.0346747,  -0.0203395, -0.0154173, 0.031334, 0.0262655,  0.00224767, 0.0172362, 0.0265062, -0.038035, 0.0301112,  -0.0243593,  0.0414339, 0.0451315,    0.0960202,  0.0864642,  0.0673048,   0.0375359,   0.0944946, 0.036669,  0.0578383, 0.0477528},
                                                                {0.011804,  -0.0108315, -0.0148268, 0.0290907, -0.0321819, 0.0289785,  0.016873,   -0.0234544, 0.00653713, 0.0468413,   -0.00410596, 0.0145675,   0.00133672, 0.0338526, -0.00800818, 0.0219828,  0.0272399,  -0.037904,  0.000337112, -0.0222561, 0.0111167,  -0.02895,  0.0180652, 0.0372134,   -0.0268938},
                                                                {0.0391821,  -0.0394245, -0.00421925, 0.0145944, -0.0166062, 0.0360033,  -0.0389615,   -0.00759398, -0.0100245,  -0.00673693, -0.0279515, 0.024065,  -0.037399,  0.0372219,  -0.00597352, 0.00321439,  -0.0249007, 0.0385082,    0.0317517,  -0.0147183, -0.0444428,  -0.0268456, 0.0408177,  -0.00474454, 0.0186405},
                                                                {0.0200414,  0.000466362, -0.0390842, -0.0235088, 0.0226218,  -0.00499683, -0.00579998, 0.00352646, 0.0229466,   -0.0344984, -0.0195158, 0.0264518, 0.00946806, -0.00541125, 0.00877971, 0.0320743, -0.0134945, 0.00798949, -0.0239902, -0.0411824,  -0.0457198, 0.00470364, -0.021372,  0.0478193,  0.00360444},
                                                                {-0.043034,  -0.0155584, 0.01389,   -0.0425494, -0.0136172,  0.0104553,  -0.00727323, -0.0382959, -0.0401407, 0.000288676, -0.00164152, 0.00100106, -0.00428327, -0.00714999, -0.0308735, 0.0182834,  0.0444369,  0.0115464,  -0.0411644,  -0.033078, 0.0141336,  -0.0249119, -0.0280341, 0.0206403, 0.0228296},
                                                                {0.030832,   0.0233214,    0.00397517, -0.0328605, 0.00906477, 0.0111248,  -0.00125119, -0.00299969, 0.00562291, 0.000493974, 0.0130933,  0.0433407, 0.0198971, 0.0024231, 0.00805547, 0.0100079, 0.0295246,  0.00608106,        -0.0271838, 0.0122188, 0.0400462, 0.00690566, -0.0502115, -0.0173953,  -0.0395861},
                                                                {-0.0228776, 0.0127145, 0.00537419, -0.00203479, 0.000142038, 0.00437803, 0.0174998,  -0.0274117, -0.0283462,  -0.0121211, 0.0456775, -0.0151794, 0.018327,  -0.0197859, 0.0141643, 0.0144876,  0.000815491, -0.0187011, -0.0348825,  0.0122389, -0.0194639,  -0.00645966, -0.00384462, -0.0188769, -0.00730397},
                                                                {-0.00798795, 0.0361338, 0.0211809,  -0.028986,  0.0435885,  0.00586178,  0.0714831, 0.0103799, 0.0893413, 0.0683302,  -0.00516072, 0.0163892, 0.00637049, -0.00697006, -0.0108477,   -0.00527127, 0.0021458,  -0.059419,   0.0179344,  -0.064917, 0.014252,   -0.00922859, -0.040268, 0.0289978,  -0.0131975},
                                                                {0.0149063, -0.00347616, -0.0408025, -0.0354836, -0.0493283, -0.00397038, 0.0395944, 0.0132072,  0.0369898, -0.0175434, 0.00519899, 0.0678152, 0.0626374,   -0.00976002, -0.0263756, 0.0264579, 0.03205,   -0.00176096, -0.0200349, 0.00567852, 0.0135142, -0.0252004, 0.00880843, 0.000574279, -0.0203634},
                                                                {-0.0230446, 0.00985364, 0.0212001,  -0.0196618,  0.00780229, -0.0323597, -0.0269155, 0.0145217,  0.0157026, -0.00739551, -0.00997144, -0.00634764, -0.0263921, -0.0432097, -0.0056069,  -0.0382243,  -0.00191331, -0.0279346, -0.0390441, 0.00318105, 0.0231178, -0.0510458, -0.0175983, 0.0470836, 0.00624285},
                                                                {-0.0448882,  -0.0358049, 0.0497897,  -0.0108792, -0.0270764, -0.03076,    0.0383663,  0.0217393,  0.00974446, -0.0108536, -0.0141754, 0.040067,  0.0279887,  -0.000699077, -0.0522192, -0.015762,  -0.003283,  -0.00854266, 0.0420529,  -0.0243678,  0.0429316,  0.0219718,  0.0250591,  -0.0225186, -0.039197},
                                                                {-0.0241347, 0.0173971, 0.0260786, 0.0151541, 0.0206457,   0.00713279, 0.0311899, 0.00955618, -0.0200317, -0.036952,   -0.0470167,  -0.0456353, -0.0333356, 0.000449326, -0.02412,   -0.0435136, 0.00244766, -0.0481601, 0.0259174,  0.0278652,  -0.0565282, -0.0128141, -0.015584, -0.0424491, -0.0366166},
                                                                {0.0110689,   0.0242046,  -0.0466518, 0.0177097,  -0.00479673,  0.0438076,  0.0195632,  -0.00822297, -0.00139199, 0.0236992,  0.0164488,  -0.0458919, -0.0301612, -0.00419701, 0.0219087, -0.0379356, -0.0468354, -0.0476783, -0.0328077, -0.0235427, -0.0352294, -0.0603228, 0.00969842, -0.014057,  -0.0230743},
                                                                {-0.0289203, 0.0398778, 0.00575217, 0.0257819,  -0.043041,  -0.0140566, -0.0185023,  -0.0265948, -0.00995332, 0.00400524,  -0.0184137, 0.00152011, 0.0201326,  -0.000611862, -0.0523959, 0.0345174, 0.0042873,  0.00284178,        0.0175376,  -0.0065791, -0.000166531, 0.00709658, 0.00880513, 0.0636001,  0.0159178},
                                                                {-0.0386404, -0.0425945, -0.0421754, -0.0172307,  0.00646786, 0.0397402,  0.0074409, -0.003731,  -0.00748909, -0.0031684, 0.0181335,  0.0172275,  0.0119325, 0.000621515, 0.0215294, 0.0460379,  -0.0374972, -0.00729005, 0.0160222, 0.00188377, 0.0323744,  -0.0380576, 0.0166708,  -0.00233671, 0.00592996},
                                                                {-0.0220373, -0.0285451, -0.0070975, 0.0127433,  -0.0239843, -0.00012684, 0.0203211,  -0.00860936,  0.0184402, -0.043504, 0.0273562,  0.0273329,  -0.0174899, -0.0124525, 0.00428926, 0.0277281,  0.00968005, -0.00460388, 0.00129189, -0.0294741, 0.0413946, -0.0194698, 0.0126256,  0.0236837,  -0.0318845},
                                                                {-0.0293658, 0.0404151,   0.0148683,  0.0532676,  0.0468546, 0.00656924, -0.0309491, -0.0197955, -0.0373579, 0.0147139, -0.0148166, 0.0169157,  -0.0024255, -0.0445677, -0.00407154,  -0.0117758, 0.0193231,  -0.0491882,  -0.0182584, -0.000270851, 0.0335786,  -0.0236811,  0.0233153,  -0.0185489, -0.0361295},
                                                                {-0.0119443,  -0.0170581, -0.0428174, -0.00820237, 0.0320459, -0.0104259,  0.00306205, 0.0010499, -0.0528613, -0.0245712, 0.0177581,  -0.00433131, -0.0279069, 0.00606032, -0.00556103, 0.0268547,  -0.0335537, -0.0248975, -0.0304636,  -0.0014182, 0.0299969, 0.0467459, 0.0583735, -0.00546816, -0.0305165},
                                                                {0.018892,   0.0313354,  -0.0414115, 0.00634449, -0.0103772, -0.0303332, 0.0307468,  -0.0429189, 0.00722902, 0.00201799, -0.022127,  -0.0380311, -0.0503551, 0.00920168, 0.00147819, -0.0103358, 0.0107836, -0.0180369, 0.0295644, 0.019449,   -0.0274645, -0.0332418, -0.00863405, -0.0139857, 0.0347327},
                                                                {0.0358038,  -0.0461728, -0.0299941,  0.0136563, 0.0251586,   0.0363407,  0.0396898,  -0.032974, -0.0536901, 0.0134123, -0.0463866, 0.00113573, 0.00533407,  -0.0463077, -0.0409405, -0.00047124, 0.012465,  -0.0193692, -0.0559089, -0.0329343, 0.0228868, -0.0232691, -0.0360485, -0.0173055, 0.0302332},
                                                                {-0.0236541, -0.033996,  -0.0180238, 0.0185029,  -0.0254945, -0.0287216, -0.00749095, 0.0078843, 0.0292374, -0.0463648, -0.0672033, -0.0272002, -0.0127069,  0.0179428, -0.0217023, 0.0267252, 0.011766,  -0.0248935, 0.0194846,  -0.0496919, 0.00651729, 0.032636,    -0.0117236, 0.030754,  -0.00185802},
                                                                {0.0220525,  0.0335789,  -0.00600046, 0.0335769, 0.0558175, 0.0342919,   -0.0193606, 0.0181103,  -0.0253193, -0.0208243,  0.0308181, -0.00863886, -0.0167984, -0.0316913, -0.0347017, 0.00589598, 8.95151e-06, -0.0338229, -0.0139525, 0.00651146, 0.00256258, 0.0164423, -0.0452407, -0.00462715, -0.0343737},
                                                                {0.0387321, -0.0100077, -0.0273543, 0.00421041, -0.0392451, 0.0218486,  -0.00365168, 0.0161956,  0.00324022, -0.0153046, 0.05226,    0.0372115, 0.00434048, 0.00691883, 0.0110492, -0.0274617, 0.03837,   -0.0167165, 0.0213024, -0.00268492, -0.0146594, 0.0378255, -0.0243652, 0.00454103, 0.0364164},
                                                                {0.0188908, 0.0396095,   -0.0228412,  0.0175988, 0.0335067,   -0.00222021, 0.0452796, 0.0219438,   -0.0152897, -0.0153121, -0.0381136, 0.00160986, 0.0460428,  0.0149103,  0.012072,   -0.020803,  -0.00180168, 0.0215293, 0.0298198,  0.0199506,  -0.0111147, -0.0216928, 0.000643242, 0.0193776, 0.00709579},
                                                                {0.0211861,  -0.0455798,  -0.0234709, -0.0449939, -0.0226419, 0.0102192,    -0.0426667, 0.0144101,   -0.0150897, 0.0305073, 0.0378711, -0.0393478, -0.0358986, 0.0308827, -0.000453334, 0.00395561, -0.032646,  -0.0536718, 0.0271696, 0.00985262, -0.0343754, 0.0150585,  -0.0405455, 0.0260934, -0.00839432},
                                                                {-0.00405847, -0.0216137, -0.0284036, -0.00349631, 0.0322282,  0.0562579,  0.0280847, 0.0141361,  0.0461761,  0.0316113, -0.0153152, -0.0132464, 0.0241652, 0.0341555,  0.0220608, -0.0217795, -0.0480004, -0.0463083, -0.0286769, 0.0309273, -0.0520326, -0.0344642, -0.0240315,  0.00476425, -0.0377886},
                                                                {0.0137054,   0.0467618,  -0.0105536, -0.000859622, 0.0199925,  0.0357268,   -0.00668708, 0.0299802, -0.0187628,  0.0521298, -0.0200462, 0.0317907, 0.0222946,  0.0593796, -0.0206713, -0.0114117, -0.0298722,  0.0351698,  0.0104105,  0.0387852, 0.00569732, -0.0248471, 0.0275976, -0.00390932, -0.0418406},
                                                                {-0.0488407, -0.055779,  -0.0482501, -0.00144554, 0.000842305, 0.0322646,   0.0153163, -0.00383202, 0.0391403,  0.0231164,  0.0272481, 0.103729,  0.0836727, 0.0237131, 0.0230123,  -0.018825,  0.0194232,  0.0667909,  -0.0126696, 0.00814122, 0.0191094,  -0.0519054, -0.0523414, -0.0169977, -0.0474668},
                                                                {0.0113022,  -0.0111541, -0.0264198, -0.0371764, -0.0118721, 0.00034396, -0.0389479, -0.0467772, 0.0377191, 0.0270101, -0.0317914, -0.00406681, -0.00780517, 0.027396, -0.0211436, 0.0280218, 0.0135839,  -0.0116187, -0.00904291, -0.0287109, -0.043211,    0.0335861,  -0.0349338, -0.00525466, -0.00152422},
                                                                {0.0423882, -0.0255072, -0.00502006, 0.0398017, -0.0427315, 0.0128578, 0.0432658, 0.0335412, 0.029192,   0.0306435,  -0.0317784, -0.0253629, -0.00542368, -0.0448831,  0.0316099, 0.0154443, -0.0117677, 0.0123438, -0.0422633, -0.0125835, 0.035643,   -0.0407636, -0.0373229, 0.0371261,  0.0239131},
                                                                {0.0184353,  -0.0174708, -0.0310311,  -0.0389449, -0.0287567,  -0.0291059, -0.00986319, -0.0247067, -0.0203216, 0.0185252, -0.0243477, 0.0271498, -0.0150879, -0.0419058, 0.0271249,  0.0239446, -0.0406384, -0.0297722, 0.016713,   0.0287004,   0.0131561, -0.0375068, 0.0121983,   -0.0327358,  0.0332652},
                                                                {-0.0321208, -0.0406959, 0.0227666,  0.0234774,  -0.0310603, -0.0084431, 0.0325153,  0.00726795, -0.00192832, -0.0427938, 0.0351807, -0.0166339, 0.0525522,  0.0364239,  -0.000586823, 0.0301032, 0.0386101, 0.0168853, 0.0445353,  -0.00901288, 0.00119475, 0.0291047, 0.0188715, -0.0264925, -0.0230797},
                                                                {-0.0297932, -0.030524, 0.0317254,  0.011159,   -0.0390071, -0.0564094, 0.00425511, 0.00790463, 0.00074114, -0.0308029, 0.00584595, 0.0439219, 0.0615915,  0.0186295,  0.0462435,  0.0427643,   0.0296231,  -0.00807032, -0.0217291, 0.0293974,  0.0194275, 0.0442447, -0.0363749, -0.0455779, -0.0109932},
                                                                {0.0376329, 0.02114,   -0.0207287, -0.0332071, 0.0391149,   -0.0266982,  0.0146344, 0.00315235, 0.0342845,  -0.023625,  -0.0452676,  -0.037589,  -0.0264766, 0.0167646, -0.00270109, -0.0225167, 0.0158131, -0.0411788, 0.0296316,   0.0164202,  0.0085697, -0.00053501, 0.0365701,  0.021246, -0.0349039},
                                                                {-0.00412691, -0.0318361,  0.00609924, 0.0408805, -0.0313765, -0.00773796, -0.0466985, -0.0527277, 0.0266871,  -0.00655214, -0.0436253,   -0.0658782, -0.0388605, 0.0315675,  0.0348028, -0.0552401, -0.0570824, 0.00798213, 0.0228847, -0.0179765, -0.0438445, -0.0375925, -0.0360353, -0.0122165, -0.0114959},
                                                                {0.0341611,  0.0313014,  0.0103648,   0.0112398,  -0.0285909, -0.0205428, 0.0120612, 0.0354749, -0.0137108, 0.024269,  0.0110727,   -0.0321526, 0.0421073,  -0.0447471, -0.037726, 0.017205,   0.0237076,  -0.0036042, -0.0388371, -0.0213912, 0.0223702, -0.0352585,  0.00311928, -0.0209426, -0.0187599},
                                                                {0.0029443, -0.0435967, 0.00811902, 0.0324785,  0.0304622, 0.0171143,  -0.0307008, 0.00925542, 0.00808559, 0.0321488,  -0.0354378, 0.0413961, -0.015934, 0.033418,    0.0568193,   -0.0174234, 0.00240194, 0.0107034, -0.010155, 0.0405467, -0.00664929, -0.0411937, -0.0454384, 0.00853229, 0.00214266},
                                                                {-0.0073464, -0.026202, -0.0207623, -0.0425809, 0.0202833,  0.0228551, 0.0299314,   0.00399659, 0.0200871, -0.0217129, 0.00361728,  -0.0111652, -0.0236218, 0.0454659, 0.0312724, 0.0231763,  0.0250114,  0.0345299,  -0.0117073, -0.0246173, -0.0475959, -0.0512378, -0.0422496, 0.00776245, 0.00172571},
                                                                {0.0427185,  0.0305539,  0.0162719,  -0.0182851, 0.0191428,  -0.0099904, -0.0179902, -0.00360675, -0.0368078, -0.00588811, -0.0434607, -0.0181829, -0.0379488, -0.00219687, -0.0278822, -0.00571475, 0.00297709, 0.0177353, 0.0382817,  0.0136851,  -0.0341578, -0.00669893,  -0.0126966, -0.0049522, -0.0298855},
                                                                {-0.0163267, -0.0349787, 0.0189835,  -0.0228238, -0.0310472, -0.0320215, -0.0023646, -0.026673,  -0.00369588, -0.031883, 0.0454067,  0.0205591,  -0.00769644, -0.0188449, 0.0318185, -0.0301191, -0.0429963, -0.00584844, -0.00749585, -0.0473836, -0.0271292, -0.0268211, -0.0333486,  -0.00459809, 0.0348291},
                                                                {-0.0443317, 0.0214588,  0.0393616,  0.0252526,  0.0187325,  -0.0362443, -0.035316, -0.00221464, 0.0410806,  0.0152988,  0.0203929,  0.0556836, 0.00574289, 0.0537931,  -0.0156834,  -0.00807516, 0.0164043, 0.0727797, 0.0814804, 0.00927312, -0.0168872, -0.0313365, 0.0633258,  -0.0207063,  0.0295728},
                                                                {0.034896,    -0.0238712, 0.0121059,  -0.0398509, -0.0317776, 0.0206706,  0.00201212, 0.00198384, -0.0295415, -0.0003371, 0.0351158, 0.0284745, 0.0406333,  0.00798849, -0.0390773, 0.00176297,   -0.0070067, -0.0181014,  0.000336216, 0.0233977,  0.00422754, -0.0272089, 0.0499933,  -0.0129372, 0.0281913},
                                                                {0.0216921,  0.0445015, 0.000168543, 0.00625171, 0.0114243, 0.0197006, -0.0342946, 0.0412882,  -0.000508577, 0.0280645,  -0.000728538, 0.0200743, -0.028671, 0.00806389, 0.00833313, -0.0393026, -0.00511164, -0.0295207, -0.0403186, 0.0357299, -0.0307141, -0.0456501, -0.0407076, 0.0400227,  -0.0400412},
                                                                {0.0421602, 0.031921,   0.0149275, -0.00760219, 0.0113893,  0.0189283, -0.0394624,  0.0169157,  0.0404885,  0.0300481, 0.00302595, -0.020342, 0.0482973,  -0.00172933, 0.0311519,   0.0232336, 0.0169641,   0.0358628,  -0.0226057, -0.0249838, -0.0124767, -0.0389264, -0.00167821, 0.0127125,  -0.0402431},
                                                                {-0.00307082, 0.0375079, -0.02261,  0.0241309, 0.0421152,   -0.0114934, 0.0384512, 0.0378949, 0.0415564, -0.0129287, -0.0033225, 0.0127063,  -0.0176174, 0.000890074, 0.0426573,  0.0140132, -0.0160256, 0.00378939, 0.0217366, -0.0273393, -0.00477793, 0.0207298,  0.0354092,  -0.0293594, -0.0186128},
                                                                {-0.0367555, 0.02757,   -0.0347377, 0.0478517,   0.0141144,  0.0155965, -0.000784785, -0.0237877, -0.00267585, 0.0173004, 0.0458214, 0.0090384, -0.00986296, -6.62658e-05, 0.026162,  0.0289767,  0.0262571, 0.00766926,  0.00552154, -0.00764805, -0.0139594, 0.0252822, 0.0679049, 0.0621685,  0.0283238},
                                                                {-0.023624, -0.0471896, -0.0146308, 0.0157374,  -0.0185129, -0.012242,  -0.0175303, 0.0205201,  0.0469567,  0.0308388,   -0.0442018, 0.0394344, -0.0108842, 0.0394109, 0.0277099, 0.0169094,   0.0487091, 0.0585808, 0.00148553, 0.00542836,  -0.0545618, -0.0031425, 0.0331201, -0.00867209, -0.024615},
                                                                {-0.0187278, -0.0270348, 0.00187145, -0.00931263, -0.0304743, 0.0453622, 0.0341908,  -0.00831162, 0.00383209, -0.0317976, -0.0331837, 0.0233167,  -0.0232171, -0.0259483, -0.0379926, -0.0391898, 0.00799124, -0.00516822, -0.0158196, -0.00469465, -0.0500213, -0.0100694, 0.034041,  -0.0117255, 0.00113898},
                                                                {-0.0199597,  -0.0302331, -0.0470626, 0.0241145,  0.00836413, 0.0169402, -0.0402565, -0.0226507, 0.0304263,  -0.0329044, 0.0129732,  -0.0393084, -0.0358482, -0.0336125, 0.0150114, -0.0157294, 0.00725665,  -0.0439306, 0.0342664,  -0.0269005, 0.00300051, -0.0165155, 0.0324295,  -0.0106851,  0.0424262}},
                                                        {{0.0111452, 0.0214061, 0.0687376,  0.062753,  -0.0097982, -0.0045877, 0.00059125, 0.0382741,   0.0692399, -0.0240434, 0.00947216, 0.0385408,  0.0435236,  0.0615116, 0.0416076, -0.0390097, -0.0417169, 0.0151983,   0.0608362, 0.00883882,  -0.0401907, -0.0230137, 0.049076,   0.0608768, 0.0504488},
                                                                {-0.0262289, -0.0297652, -0.0142716, 0.00283265, -0.00847504, -0.0262626, -0.0416967, -0.02454, -0.0123295, -0.0617562, 0.0243651, 0.0153887, 0.0177689, -0.0422747, -0.00707854, 0.0623879, -0.000610727, 0.00134396, -0.0315004, -0.00346995, 0.000966561, 0.0323579, 0.0769279, 0.0570635, 0.0289574},
                                                                {0.0454316, -0.0314481, 0.00686538, 0.0484373, -0.0315532, -0.0090278, 0.00539585, 0.00143481, -0.0152588, -0.00149502, -0.00734612, -0.00654482, 0.0232091,  0.0353256, -0.0278187,  -0.0390783, 0.00789026, 0.00284264, -0.0160283,  -0.0279816, -0.0237621, 0.0404292, 0.0598701, 0.000395513, 0.000991492},
                                                                {-0.0299046, 0.0289346,  -0.00445837, 0.0154693, -0.0380455, 0.00302158, -0.000872039, 0.0303523,   -0.00820767, -0.032161,   -0.0011792, 0.0104731, -0.0409921, -0.0296581, -0.0205781,  -0.00628431, 0.0296354,  -0.000804274, -0.0332409, -0.0368767, -0.00786672, -0.0100082, -0.0156731, -0.0278286,  -0.0258468},
                                                                {-0.0268991, 0.0086293,   0.0233972,  -0.0399138, -0.0470906, 0.00736785,  -0.0136902,  0.00479192, -0.00558729, 0.0118972,  0.0481933,  0.0282175, -0.0338202, -0.0488259,  -0.0333209, 0.0701404, 0.0600975,  -0.0193266, 0.029491,   -0.00673269, 0.0255243,  -0.0176082, -0.0216489, -0.0317515, -0.0189168},
                                                                {0.00283216, 0.0117937,  0.0323666, -0.0221766, -0.00630659, -0.0343983, -0.0267821,  -0.0297195, 0.0179823,  -0.0158665,  -0.00766797, -0.0279834, -0.0489096,  0.0192914,   -0.0419595, -0.0284768, 0.00715147, 0.00444432, -0.00807851, -0.032398, -0.0304195, 0.0126221,  0.00661658, 0.0159156, -0.00916186},
                                                                {-0.0193211, -0.000547358, -0.0157824, -0.0129677, 0.0310979,  0.00331432, -0.00254168, 0.00302796,  -0.0210725, -0.0201162,  0.00845827, 0.051185,  0.0492263, -0.012428, 0.00726103, 0.0111661, -0.0164616, -0.0315427, -0.0489084, 0.0113544, 0.0340786, 0.0247148,  0.0127547,  -0.00419607, -0.00235878},
                                                                {-0.025662,  0.0212055, 0.0189266,  -0.0213349,  0.0061761,   0.0344473,  0.00989348, -0.0204797, -0.00573309, -0.0274657, 0.00306,   0.028072,   0.0187681, 0.0496673,  0.0505504, -0.0281962, 0.0328507,   0.0254329,  0.000113956, 0.0159563, -0.00899269, 0.0225859,   -0.0198967,  0.00943356, 0.0262406},
                                                                {0.00959323,  0.02588,   -0.0489082, -0.0334411, -0.0292887, -0.00535428, 0.0103725, 0.0161699, 0.0497832, -0.0187308, 0.0525653,   0.0349993, 0.0125323,  0.0508338,   -0.000935133, 0.0474764,   -0.0101837, -0.00230091, -0.0187122, 0.018446,  -0.0153773, -0.0474248,  0.0111279, -0.0332261, -0.0127984},
                                                                {0.0263951, 0.0371203,   0.0328889,  -0.0689752, -0.0572537, 0.0186033,   0.0576943, -0.0310508, 0.0112133, -0.0369863, -0.0177304, 0.0247459, -0.00414027, -0.0119849,  0.0513438,  0.0308422, 0.0311621, 0.0316195,   0.033816,   -0.034425,  0.024102,  -0.0383942, -0.0189766, -0.00400827, 0.0181765},
                                                                {0.0284531,  0.0255357,  -0.0357535, -0.00396116, -0.0189407, 0.0308016,  -0.030733,  -0.0146711, 0.0273068, -0.0107088,  0.0198198,   -0.0133166,  -0.0149172, -0.0136359, 0.000570898, 0.000680849, 0.0345873,   0.00760996, 0.0345747,  -0.0356431, 0.0281828, -0.0491642, -0.0234967, -0.034027, 0.00684164},
                                                                {-0.00159193, -0.0415886, 0.00656466, -0.0274305, 0.0272829,  -0.00032264, -0.0101163, -0.0195991, 0.045267,   -0.046772,  0.0185086,  0.0153234, -0.0335858, 0.0404656,    -0.0052835, -0.0322108, -0.0144035, -0.00710319, -0.0157839, -0.00585947, -0.0131159, -0.0339709, 0.00500658, -0.0331117, 0.0343231},
                                                                {0.0431001,  0.0571341, 0.0296366, 0.0265962, -0.00197533, 0.0548786,  0.0126628, 0.0214142,  -0.0147852, -0.00705603, -0.00230175, 0.0679948,  0.0276556,  -0.0395533,  -0.0338753, 0.0330785,  0.0161559,  -0.014544,  -0.0105823, -0.0155031, 0.016141,   0.0280229,  0.0221144, -0.015458,  -0.0431418},
                                                                {-0.00951111, -0.0162203, -0.0132255, -0.0119282, -0.000668327, -0.0400885, -0.0341753, 0.0470505,   0.0235046,   -0.0206587, -0.0128239, -0.0126904, 0.049587,   -0.0101674,  0.0279503, -0.0602121, 0.00362046, 0.0122701,  0.00553195, -0.0106028, 0.0436626,  0.0121015,  0.0385112,  -0.0567563, 0.0197721},
                                                                {0.039183,   0.0370663, 0.0025068,  -0.0373269, 0.00824782, 0.0301259,  -0.00866739, 0.00563494, -0.0057687,  -0.00102343, 0.0360258,  0.0119704,  -0.0342076, -0.044364,    -0.0210854, -0.0258521,-0.0240587, -0.0317274, -0.0448867, -0.0405472, -0.00468012,  0.00393374, -0.0427435, -0.0395067, -0.0342314},
                                                                {-0.0174629, 0.0387413,  0.00462218, -0.00496565, -0.0402058, -0.0288904, 0.0290261, -0.0301044, 0.0520717,   0.012239,   -0.0448436, -0.0119689, -0.021131, 0.0101275,   0.0173376, -0.0269224, 0.0204113,  0.0115975,   0.0571399, 0.0582118,  -0.0157864, -0.0280393, -0.0565889, 0.048871,    0.0613245},
                                                                {0.0343614,  0.0128673,  -0.0113658, -0.0244249, 0.0452343,  0.0371514,   -0.0425103, -0.000130784, 0.0323683, 0.0365266, -0.0160755, 0.00186808, 0.0193998,  0.0306265,  0.039748,   -0.0218154, 0.00118873, 0.00876941,  0.0295887,  0.00826758, -0.029455, 0.0188712,  -0.0199358, -0.0184353, 0.0230551},
                                                                {-0.0199034, -0.00757668, 0.00936449, 0.00388894, 0.0415092, -0.0399883, 0.0144226,  -0.0600156, -0.0431923, 0.0195241, -0.0162577, 0.00557196, -0.0337653, -0.0288888, -0.000791223, -0.0421934, -0.0607949, -0.00700759, -0.0180358, 0.0282742,    -0.0107624, -0.00671286, 0.00454409, 0.0519776,  0.0146275},
                                                                {-0.00527168, -0.0261813, 0.0434279,  -0.0433571,  0.0036021, -0.00673802, 0.0301819,  0.0463499, 0.00741173, 0.028798,   -0.0408672, 0.0301552,   -0.0256224, -0.0483527, -0.0488236,  -0.0170272, 0.0153385,  0.0151407,  -0.00899206, 0.0403727,  0.0310164, 0.0591481, 0.0474933, 0.0373783,   -0.0119572},
                                                                {0.00662188, -0.0513499, 0.022221,   -0.0358734, 0.0396327,  -0.0143499, -0.0364621, -0.0257765, 0.0366919,  0.0139978,  -0.0115729, 0.0106868,  -0.0166122, -0.0352178, 0.0199969,  -0.0327593, 0.0279175, 0.0465126,  0.0177942, -0.0168576, -0.0195465, 0.0255866,  -0.0107609,  0.00213476, -0.0113748},
                                                                {-0.0271109, 0.0169626,  -0.00907441, 0.011862,  -0.00509303, -0.0169598, -0.0150831, 0.0272795, 0.00854549, 0.0366177, 0.00130756, -0.0646169, -0.00159759, 0.0915029,  0.0110538,  -0.0219049,  0.0091155, 0.0135896,  0.0665597,  0.0403909,  0.0209017, 0.00216972, -0.0176056, 0.0360126,  0.017171},
                                                                {-0.0182101, -0.0120828, -0.0206816, -0.0242444, -0.0126865, -0.0137138, 0.041468,    0.0436763, 0.0186151, 0.0215389,  0.0629368,  0.031119,   -0.00400534, 0.0156596, -0.0366142, 0.0943117, 0.0111441, 0.0562237,  -0.0321628, 0.0183965,  0.107646,   0.000146582, 0.00764901, 0.0259787, 0.0300986},
                                                                {0.00835903, -0.0126845, 0.0312573,   0.0382055, 0.051263,  -0.00584904, 0.0667478,  -0.0185292, 0.0150497,  -0.00888862, 0.0192203, 0.0370774,   0.0340956,  0.00651099, -0.0457359, -0.0361057, 0.0113699,   0.0235073,  -0.0170129, -0.0464816, -0.0403509, 0.0418298, -0.0338514, -0.0218917,  -0.0304264},
                                                                {-0.035616, 0.0480467,  0.0387901,  -0.0373967, -0.0383844, 0.00562777, 0.0132326,   0.00507553, -0.0305659, -0.0383709, -0.0255991, -0.027306, -0.0101793, 0.0225409,  0.0291646, -0.0287813, 0.0335458, 0.0124209,  0.02949,   0.00409783,  0.0295057,  0.0369492, 0.040998,   -0.0417652, 0.011094},
                                                                {0.0116365, -0.00725318, -0.00612309, 0.0118113, -0.00885055, -0.00344479, 0.0287655, -0.00312306, 0.024219,   -0.0269113, 0.00643188, -0.0122493, -0.0131843, 0.00504313, 0.00501371, -0.0228229, 0.0378537,   0.0010786, -0.0159995, -0.0308154, 0.0182755,  0.0294065,  -0.0172994,  0.0250657, 0.00391557},
                                                                {-0.0189869, 0.000980851, -0.0216779, -0.0212062, 0.0439679,  -0.000405666, -0.0412976, 0.000219678, 0.00993249, 0.0454954, 0.0352018, -0.0428802, 0.0483599,  0.065787,  0.0027998,    -0.0446154, -0.0790229, 0.0282665,  0.125071,  0.0479316,  -0.0388288, -0.0498337, 0.0584282,  0.115887,  0.055302},
                                                                {0.0408286,   -0.0132532, 0.0306442,  0.023036,    -0.0182449, -0.0319386, -0.037887, -0.0089088, 0.00811697, 0.0338139, -0.0282731, 0.0154161,  0.0138779, -0.0307008, 0.0173334, 0.0125785,  -0.024623,  -0.0631282, 0.0447809,  0.0186619, 0.00114453, -0.0692566, -0.00406697, 0.0380954,  0.0633586},
                                                                {-0.00172364, -0.0287413, -0.0156258, 0.0254611,    -0.0331937, -0.00788512, 0.00964132,  0.031717,  -0.00426464, 0.012118,  -0.0386688, 0.0105546, -0.0171042, -0.037885, -0.0129349, -0.0417157, -0.00292557, -0.0409747, 0.00360375, 0.013908,  0.0103072,  -0.0121163, 0.0271786, -0.0500795,  -0.0061068},
                                                                {0.0694844,  0.00719632, 0.0272572,  0.0055685,   0.0218418,   -0.00108906, 0.0528378, 0.0210115,   -0.0311304, -0.0192641, 0.066582,  0.0208346, 0.0315786, 0.0398836, -0.0158965, -0.0404795, 0.00064669, -0.0458137, 0.0420712,  0.0214754,  -0.0294227, -0.0716142, 0.00798443, -0.0380188, -0.038087},
                                                                {-0.0286567, 0.00510436, 0.00472005, 0.0396853,  -0.0054875, 0.033812,   -0.0336524, 0.00324664, 0.0512558, 0.0366155, 0.0154339,  0.0190149,   0.0530778,   0.061568, -0.015941,  0.0425509, -0.0163808, -0.0249069, 0.0504514,   -0.0340801, -0.000843235, -0.0596294, 0.0287889,  0.0155008,   0.0185284},
                                                                {0.0389951, 0.060115,   0.0256978,   0.019685,  0.00374602, -0.014087, 0.0280188, 0.0131402, -0.0427097, -0.0408594, -0.0489571, 0.0734447,  0.0103367,   -0.00918153, 0.0047435, 0.0345159, 0.0840146,  0.0635886, -0.0443392, -0.0232791, -0.0131222, 0.0791887,  0.00667727, -0.0335962, 0.0273652},
                                                                {0.00576462, -0.0235738, 0.000998885, -0.0266601, -0.00561583, -0.0187478, -0.0491152,  0.0275791,  0.00923532, 0.0257761, -0.0260993, 0.0263672, 0.00350195, -0.0343259, -0.0155265, 0.012814,  -0.0283393, 0.0045276,  -0.0206577, 0.000597529, 0.0520244, 0.0128925,  -0.00703738, -0.00206791, -0.00866661},
                                                                {0.0292015,  -0.014307,  0.00163265, -0.0188393, -0.0335035, 0.0306784,  -0.0229826, -0.0249523, -0.02658,    0.0340654,  0.0512657, 0.068819,   0.00387513, -0.0360479, 0.0152875,    0.0708423, 0.0720872, 0.0103996, -0.0173319, 0.033484,    -0.0333445, 0.0137646, 0.0548227, 0.01364,    -0.0363533},
                                                                {0.0185011,  0.0310581, -0.0159287, -0.0313342, -0.0260454, -0.0180349, 0.027101,   -0.0403271, 0.0104721,  -0.041282,  -0.0424533, 0.0115199, -0.0134404, -0.0552404, -0.0293282, -0.00440231, 0.00218685, -0.031483,   -0.0192044, -0.0453239, 0.0402021, 0.058951,  0.0238424,  -0.0533671, -0.0611196},
                                                                {0.0412366, 0.0212415, 0.010863,   0.0359989,  0.000937785, -0.00538879, 0.0142973, 0.00636223, -0.0360028, -0.0452367, -0.00175794, -0.0334509, -0.0326031, 0.0113963, -0.0416716,  0.0396711,  0.0200257, -0.0375914, -0.00671215, -0.0418601, 0.0508628, 0.0178893,   -0.0439536, 0.012544, 0.00946594},
                                                                {-0.0112303,  0.000516836, 0.0715293,  0.0493259, 0.00436093, -0.0149268,  0.0151691,  0.0931838,  -0.0266886, 0.050837,    -6.01837e-05, 0.0430865,  0.0840167,  -0.0135191, 0.0121925, -0.0246385, 0.00559652, 0.0968811,  0.0319041, -0.0121523, -0.0328117, 0.05472,    0.0748125,  0.0565409,  0.00215958},
                                                                {-0.0340892, -0.0336653, -0.00513862, -0.0039206, -0.0272237, 0.0015047,  0.0440929, 0.0160913, -0.0332609, 0.0328768, -0.00572704, 0.0065004,  -0.0133354, 0.0364126,  -0.036195, -0.0397307, -0.0296595, 0.0346104,  0.0144716,  -0.0259771, 0.0560061, -0.00412903, 0.0177051,  0.0414276,  0.000622505},
                                                                {0.0420108, 0.0483931,  0.034963,   0.00447661, 0.0424173, 0.00742184, 0.041809,   -0.0203197, 0.00334741, -0.0516262, 0.058048,   0.0235585, 0.0242302, -0.00425054, -0.00245889, -0.0268727, 0.038372,   0.0398711, 0.0513349, 0.0131632, -0.0267638,  0.0681429,  -0.0176363, -0.0400695, -0.0535873},
                                                                {0.0355464,  0.0429512, 0.0200018,  -0.012221,  -0.0162071, 0.041165,  -0.00299502, 0.0643522,  0.0406552, 0.0020538,  -0.00291054, -0.0446626, 0.0185599,  0.0174461, 0.0223844, -0.0133634, -0.0190212, -0.0285087, 0.0182802,  -0.0112059, -0.0219915, -0.071335,  -0.0415287, 0.0140583,  0.0157937},
                                                                {-0.0183127, -0.0442626, -0.0427302, 0.0172328,  0.00902011, 0.0345708,  -0.0398323, -0.0212236,  0.0120569,  0.0116306,   0.0178826,  -0.027211,  -0.0164032, 0.0332274,   0.0067394,  -0.00254157, -0.0364203, 0.0206249, -0.0305141, 0.00507685, -0.0221643, -0.000834887, 0.0106812,  -0.0259973, -0.0100871},
                                                                {0.0327011,  0.0571782,  -0.0382239, 0.0169591,  -0.0289104, 0.0338171,  -0.0104924, -0.0170535, -0.0327096,  0.0100476, -0.0166532, -0.0130156, 0.0569673,   0.0285412,  0.0258438, -0.023394,  -0.0138645, 0.0213111,   0.0154546,   -0.0182666, 0.0116268,  0.0432447,  -0.00984466, 0.0323731,   0.0263512},
                                                                {-0.0227396, -0.0182448, -0.0246106, -0.0265962, -0.0370221, -0.0398526, 0.045771,  0.0435462,   -0.0432432, -0.0357719, 0.00526256, 0.0559469, -0.0300626, -0.0407137, -0.00110145, -0.0167063,  0.0116857, 0.0158856, 0.0124826, 0.00798978, -0.0053303, -0.0203056, 0.00752481, -0.00173675, -0.0594815},
                                                                {-0.00971684, -0.0562297, 0.00912522, -0.0183643, 0.0313984,  -0.0454661, 0.0188142,  0.00926486, -0.0245789, -0.0153713, 0.0262025, 0.0146859, 0.00881377, -0.0231112, -0.0423737, -0.000492567, 0.0274332,  0.000746618, -0.0207927,  -0.0189245, 0.0345237,  0.0457388,  -0.0398928, -0.0116106, 0.0230937},
                                                                {-0.0176953, -0.026058, -0.0362481,  -0.0338871, 0.0212024, 0.0366254, 0.0536979,  -0.0243561, 0.0405597,    0.00431846, 0.0253358,    -0.017441, 0.0475548, 0.0437661,  0.00455736, 0.0427877,  -0.0109489,  0.0251903,  0.00197734, 0.0487086, -0.0292899, 0.0158722,  0.00197222, -0.0183227, 0.0160844},
                                                                {0.0193123, 0.00108713, 0.0462581, -0.0421578,  -0.0119699, 0.0105379, 0.000343623, 0.00218918, -0.0490092, 0.010437,  -0.0187061, 0.0134311, 0.00287931, -0.034573,   -0.00559516, 0.0260874, -0.00346018, -0.0141633, 0.00663126, -0.0327329, -0.0176909, -0.0341432, 0.0310592,   -0.0106786, -0.0338937},
                                                                {0.0240399,   0.0484935, 0.0118163, 0.0192529, 0.000887559, 0.0382235,  0.0323101, 0.0363976, 0.0485372, -0.010019,  0.0221361,  0.00942272, 0.0251012,  0.0388019,   0.00175635, 0.0236711, 0.0233248,  0.0521483,  0.0871396, 0.0789742,  -0.0142468,  -0.0271494, -0.0444581, 0.025525,   0.0276568},
                                                                {0.0291279,  0.0180686, -0.0155264, -0.00265822, -0.0315797, 0.022872,  -0.0244181,   -0.0214506, 0.0191469,   0.0372042, 0.0640437, 0.0622148, 0.00933544,  0.000519634,  0.0474548, 0.00406494, 0.0690197, -0.00056783, -0.0283944, -0.0320343,  -0.0147218, 0.0604753, 0.0486756, 0.00503297, 0.0274816},
                                                                {0.0428296, 0.0201027,  -0.0562257, -0.0524113, -0.0183919, -0.0222528, 0.0802219,  -0.0119269, -0.0158851, -0.00501718, 0.0289741,  0.0258415, 0.0891881,  0.0145919, 0.0324552, -0.00534373, 0.0575143, 0.0998186, 0.0332356,  0.000910856, -0.0210284, -0.0448606, 0.0594259, 0.0831678,   0.0433262},
                                                                {0.0032089,  0.0126572,  -0.0125621, -0.0425796,  -0.0411645, 0.0211815, 0.00442325, -0.00662896, 0.0304516,  -0.033878,  0.00793653, -0.0202221, 0.0277867,  -0.0257056, -0.0200151, -0.0113685, -0.0124259, 0.0409443,   0.0350205,  0.0239098,   0.0133941,  -0.0449913, 0.0187963, 0.0684132,  -0.00277541},
                                                                {-0.00264468, -0.031152,  -0.0144575, -0.0219479, -0.033037,  -0.021668, -0.0264561, -0.0295309, -0.0175333, -0.0231822, -0.0306009, 0.0283153,  0.00463261, -0.0116295, 0.0393127, -0.0324565, -0.00727183, -0.0191237, 0.00401866, -0.0204299, 0.0472385,  0.0534394,  -0.0325408, -0.00896093, -0.0245041}}};
*/
    uint64_t ***sec_kernel = new uint64_t **[channels];
    for (int i = 0; i < channels; ++i) {
        sec_kernel[i] = new uint64_t * [n_kernel];
        for (int j = 0; j < n_kernel; ++j) {
            sec_kernel[i][j] = proxy->createShare(kernel[i][j], k_size*k_size);
        }
    }
    double bias[n_kernel] = {0.006893986,0.028841767,-0.011976943}; //,0.034752205,0.026028743,0.045701645,-0.014377818,-0.0099255275,-0.02820568,0.038818415,0.017578261,-0.037074637,-0.033514146,0.006196057,-0.013816738,-0.03915469,0.039954044,0.05106463,-0.032943178,-0.023218332,-0.01375858,0.011314128,-0.0073599,-0.012863117,0.0317652,-0.022406086,0.0480716,0.04324196,-0.018559806,0.044528794,-0.038379807,-0.035707023,-0.011929117,-0.03329016,0.03991615,0.013486342,-0.035355117,0.010999725,0.019694611,-0.0106660845,0.021526804,0.0036931,-0.015429371,-0.00035950931,-0.026158424,-0.023419743,-0.03008074,0.052684106,0.02619821,-0.012335571};
    uint64_t * sec_bias = proxy->createShare(bias, n_kernel);

    uint32_t params[9];
    params[0] = channels;
    params[1] = rows;
    params[2] = cols;
    params[3] = k_size;
    params[4] = n_kernel;
    params[5] = 1; //stride
    params[6] = 2; //max
    params[7] = 2;
    params[8] = false;
    proxy->SendBytes(CNN_CL, params, 9);
    cout << "call CL" << endl;
    uint64_t ***conv = CL(proxy, secret, channels, rows, cols, sec_kernel, k_size, n_kernel, 1, 2, 2, sec_bias, false);
    cout << "returned from CL" << endl;
    double *** mul = new double ** [n_kernel];

    int mul_row = (rows - k_size + 1);
    int mul_col = (cols - k_size + 1);
    for (int i = 0; i < channels; ++i) {
        for (int k = 0; k < n_kernel; ++k) {
            if(i == 0)
                mul[k] = new double *[mul_row];
            for (int k_pos_r = 0; k_pos_r < mul_row; ++k_pos_r) {
                if(i == 0)
                    mul[k][k_pos_r] = new double [mul_col];
                for (int k_pos_c = 0; k_pos_c < mul_col; ++k_pos_c) {
                    if(i == 0) {
                        mul[k][k_pos_r][k_pos_c] = 0;
                    }
                    for (int k_value = 0; k_value < k_size*k_size; ++k_value) {
                        double kernel_value = kernel[i][k][k_value];
                        int i_row = k_pos_r + k_value/k_size;
                        int i_col = k_pos_c + k_value%k_size;
                        double sec_value = input[i][i_row][i_col];

                        mul[k][k_pos_r][k_pos_c] += kernel_value*sec_value;
                    }
                    mul[k][k_pos_r][k_pos_c] += bias[k];
                    //RELU
                    if((i == (channels - 1)) and (mul[k][k_pos_r][k_pos_c] < 0)){
                        mul[k][k_pos_r][k_pos_c] = 0;
                    }
                }
            }
        }
    }

    //maxpool
    bool correct = true;
    double *** corr_conv = new double ** [n_kernel];
    int conv_row = mul_row/2; //because of maxpool
    int conv_col = mul_col/2;
    double ***rec_conv;
        if(!params[8]){
            rec_conv = new double **[n_kernel];
        }
        else{
            rec_conv = new double **[1];
            rec_conv[0] = new double *[1];
            rec_conv[0][0] = convert2double(REC(proxy, conv[0][0], n_kernel*conv_row*conv_col), n_kernel*conv_row*conv_col);
        }
        for (int k = 0; k < n_kernel; ++k) {
            if(!params[8]){
                rec_conv[k] = convert2double(REC(proxy, conv[k], conv_row, conv_col), conv_row, conv_col);
            }
            //print2DArray("computed conv kernel 1", rec_conv[k], conv_row, conv_col);
            //print2DArray("correct conv output", mul[k], mul_row, mul_col);
            corr_conv[k] = new double *[conv_row];
            for (int r = 0; r < mul_row; r += 2) {
                corr_conv[k][r / 2] = new double[conv_col];
                for (int c = 0; c < mul_col; c += 2) {
                    //find max of window:
                    double max = mul[k][r][c]; // first value in window
                    for (int max_r = 0; max_r < 2; ++max_r) {
                        for (int max_c = 0; max_c < 2; ++max_c) {
                            double next_value = mul[k][r + max_r][c + max_c];
                            //cout << "cmp: " << max << " " << next_value << endl;
                            if (next_value > max) {
                                max = next_value;
                            }
                        }
                    }
                    corr_conv[k][r / 2][c / 2] = max;
                    double cmp_value = 0;
                    if(params[8]){ //flattened
                        cmp_value = rec_conv[0][0][k*conv_row*conv_col + (r / 2 * conv_col) + c / 2];
                    }
                    else{
                        cmp_value = rec_conv[k][r / 2][c / 2];
                    }
                    if (abs(max - cmp_value) > 0.1) {
                        cout << r / 2 << " " << c / 2 << ": " << max << " (computed: " << cmp_value << ")" << endl;
                        correct = false;
                    }
                }
            }
        }
    if (!correct) {
        for (int i = 0; i < channels; ++i) {
            print2DArray("secret", convert2double(REC(proxy, secret[i], rows, cols), rows, cols), rows, cols);
            print2DArray("kernels", convert2double(REC(proxy, sec_kernel[i], n_kernel, k_size*k_size), n_kernel, k_size*k_size), n_kernel, k_size*k_size);

            for (int j = 0; j < rows; ++j) {
                delete secret[i][j];
            }
            for (int j = 0; j < n_kernel; ++j) {
                delete[] sec_kernel[i][j];
            }
            delete[] secret[i];
            delete[] sec_kernel[i];
        }

        for (int o = 0; o < n_kernel; ++o) {
            if(!params[8]){
                print2DArray("computed conv: ", rec_conv[o], conv_row, conv_col);
                delete[] rec_conv[o];
            }
            print2DArray("actual result: ", corr_conv[o], conv_row, conv_col);
            delete[] corr_conv[o];
        }
        if(params[8]){
            print1DArray("computed conv: ", rec_conv[0][0], n_kernel*conv_row*conv_col);
            delete[] rec_conv[0][0];
            delete[] rec_conv[0];
            delete[] rec_conv;
        }
    }
    return correct;

}

// Main function to run the experiments
int main(int argc, char *argv[]) {
    uint8_t role = atoi(argv[1]);
    uint16_t cport = atoi(argv[2]);
    string caddress(argv[3]);
    uint16_t hport = atoi(argv[4]);
    string haddress(argv[5]);

    Party *proxy;
    if (role == 0)
        proxy = new Party(P1, hport, haddress, cport, caddress);
    else
        proxy = new Party(P2, hport, haddress, cport, caddress);

/*
    ADD_Test(proxy);
    MUL_Test(proxy);


//    MOC_Test(proxy);
//    MMOC_Test(proxy);

    MSB_Test(proxy);
    MMSB_Test(proxy);

//    CMP_Test(proxy);
//    MCMP_Test(proxy);

//    MUX_Test(proxy);
//    MMUX_Test(proxy);

    MAX_Test(proxy);
    //MMAX_Test(proxy);
    ARGMAX_Test(proxy);

    RST_Test(proxy);
    RELU_Test(proxy);
    MRELU_Test(proxy);

    DRLU_Test(proxy);
*/
    //MDRLU_Test(proxy);//TODO
    /*DIV_Test(proxy);

    INC_Test(proxy);
    FLT_Test(proxy);
    FCL_Test(proxy);
    PAD_Test(proxy);*/
    //CL_Test(proxy);
    //MEMCPY_TEST(proxy);
/*
    EXP_Test(proxy);
    MEXP_Test(proxy);

    DP_Test(proxy);
    MDP_Test(proxy);*/
    bool all_correct0 = true;
    int counter = 0;
    while (all_correct0 and counter < 1000) {
        all_correct0 = NETWORK_M_INPUTS_TEST(proxy);
        counter++;
        cout << counter << endl;
    }
    //MAX_Test(proxy);
//    MMATMATMUL_Test(proxy);

//    INVSQRT_Test(proxy);
//    MINVSQRT_Test(proxy);

//    ppRKN_ITER_Test(proxy);
//    ppRKN_PREDICTION_Test(proxy);

    proxy->SendBytes(CORE_END);
    proxy->PrintBytes();
    return 0;
}

