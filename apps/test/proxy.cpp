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

constexpr int MIN_VAL = -100;
constexpr int MAX_VAL = 100;
constexpr int sz = 1000;
constexpr int WSZ = 10;

void MUL_Test(Party *proxy){
    cout<<setfill ('*')<<setw(50)<<"Calling MUL";
    cout<<setfill ('*')<<setw(49)<<"*"<<endl;

    double xd = MIN_VAL + (double)(proxy->generateCommonRandom() & RAND_MAX) / ((double)(RAND_MAX / (MAX_VAL - MIN_VAL)));
    double yd = MIN_VAL + (double)(proxy->generateCommonRandom() & RAND_MAX) / ((double)(RAND_MAX / (MAX_VAL - MIN_VAL)));
    uint64_t x = proxy->createShare(xd);
    uint64_t y = proxy->createShare(yd);
    proxy->SendBytes(CORE_MUL);
    uint64_t r = MUL(proxy,x, y);
    // checking the result
    xd = convert2double(REC(proxy,x));
    yd = convert2double(REC(proxy,y));
    double rd = convert2double(REC(proxy,r));
    double rcd = (xd*yd);
    if ((int)(rd - rcd) == 0)
        cout<<"MUL works correctly"<<endl;
    else {
        cout<<"MUL works incorrectly"<<endl;
    }
}

void MMUL_Test(Party *proxy){
    cout<<setfill ('*')<<setw(50)<<"Calling MMUL";
    cout<<setfill ('*')<<setw(49)<<"*"<<endl;
    uint64_t x[sz], y[sz],z[sz];
    uint32_t* params = new uint32_t[1];
    params[0] = sz;
    cout << "stored params" << endl;
    for (int i=0;i<sz;i++){
        double xd = MIN_VAL + (double)(proxy->generateCommonRandom() & RAND_MAX) / ((double)(RAND_MAX / (MAX_VAL - MIN_VAL)));
        double yd = MIN_VAL + (double)(proxy->generateCommonRandom() & RAND_MAX) / ((double)(RAND_MAX / (MAX_VAL - MIN_VAL)));
        x[i] = proxy->createShare(xd);
        y[i] = proxy->createShare(yd);
    }
    proxy->SendBytes(CORE_MMUL,params, 1);
    uint64_t *r = MUL(proxy,x, y,sz);
    // checking the result
    bool flag = true;
    for (int i=0;i<sz;i++) {
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
        cout<<"MMUL works correctly"<<endl;
    else
        cout<<"MMUL works incorrectly"<<endl;
}

void MOC_Test(Party *proxy){
    cout<<setfill ('*')<<setw(50)<<"Calling MOC";
    cout<<setfill ('*')<<setw(49)<<"*"<<endl;
    uint64_t x = proxy->generateRandom()&N1_MASK;
    proxy->SendBytes(CORE_MC);
    uint64_t r = MOC(proxy,x);
    // checking the result
    uint64_t x_reconstructed = REC(proxy,x,N1_MASK);
    uint64_t r_reconstructed = REC(proxy,r);
    if (x_reconstructed == r_reconstructed)
        cout<<"MOC works correctly"<<endl;
    else
        cout<<"MOC works incorrectly"<<endl;
}

void MMOC_Test(Party *proxy){
    cout<<setfill ('*')<<setw(50)<<"Calling Vectorized MOC";
    cout<<setfill ('*')<<setw(49)<<"*"<<endl;
    uint64_t x[sz];
    uint32_t* params = new uint32_t[1];
    params[0] = sz;

    for (int i=0;i<sz;i++)
        x[i] = proxy->generateRandom()&N1_MASK;
    proxy->SendBytes(CORE_MMC,params, 1);
    uint64_t *r = MOC(proxy,x,sz);
    // checking the result
    uint64_t *x_reconstructed = REC(proxy,x,sz,N1_MASK);
    uint64_t *r_reconstructed = REC(proxy,r,sz);
    bool flag = true;
    for (int i=0;i<sz;i++){
        if (x_reconstructed[i] != r_reconstructed[i]){
            flag = false;
            break;
        }
    }
    if (flag)
        cout<<"Vectorized MOC works correctly"<<endl;
    else
        cout<<"Vectorized MOC works incorrectly"<<endl;

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

void MMSB_Test(Party *proxy){
    cout<<setfill ('*')<<setw(50)<<"Calling Vectorized MSB";
    cout<<setfill ('*')<<setw(49)<<"*"<<endl;
    uint32_t* params = new uint32_t[1];
    params[0] = sz;

    uint64_t x[sz];
    for (int i=0;i<sz;i++) {
        x[i] = proxy->generateRandom();
    }
    proxy->SendBytes(CORE_MMSB,params, 1);
//    uint64_t *r = MSB(proxy,x,sz);
    uint64_t *r = MSBv2(proxy,x,sz);
    // checking the result
    uint64_t *x_reconstructed = REC(proxy,x,sz);
    uint64_t *r_reconstructed = REC(proxy,r,sz);
    bool flag = true;
    for (int i=0;i<sz;i++){
        uint64_t r_computed = (x_reconstructed[i]>>(L_BIT - 1)) << FRAC;
        if (r_computed != r_reconstructed[i]){
            flag = false;
            break;
        }
    }
    if (flag){
        cout<<"Vectorized MSB works correctly"<<endl;
    }
    else {
        cout << "Vectorized MSB works incorrectly" << endl;
    }
}

void CMP_Test(Party *proxy){
    cout<<setfill ('*')<<setw(50)<<"Calling CMP";
    cout<<setfill ('*')<<setw(49)<<"*"<<endl;
    double xd = MIN_VAL + (double)(proxy->generateCommonRandom() & RAND_MAX) / ((double)(RAND_MAX / (MAX_VAL - MIN_VAL)));
    double yd = MIN_VAL + (double)(proxy->generateCommonRandom() & RAND_MAX) / ((double)(RAND_MAX / (MAX_VAL - MIN_VAL)));
    uint64_t x = proxy->createShare(xd);
    uint64_t y = proxy->createShare(yd);
    proxy->SendBytes(CORE_CMP);
    uint64_t r = CMP(proxy,x,y);
    // checking the result
    xd = convert2double(REC(proxy,x));
    yd = convert2double(REC(proxy,y));
    uint64_t rd = REC(proxy,r);
    uint64_t r_computed = (xd>=yd)<<FRAC;
    if (rd == r_computed)
        cout<<"CMP works correctly"<<endl;
    else
        cout<<"CMP works incorrectly"<<endl;
}

void MCMP_Test(Party *proxy){
    cout<<setfill ('*')<<setw(50)<<"Calling Vectorized CMP";
    cout<<setfill ('*')<<setw(49)<<"*"<<endl;
    uint64_t x[sz], y[sz];
    uint32_t* params = new uint32_t[1];
    params[0] = sz;

    for (int i=0;i<sz;i++){
        double xd = MIN_VAL + (double)(proxy->generateCommonRandom() & RAND_MAX) / ((double)(RAND_MAX / (MAX_VAL - MIN_VAL)));
        double yd = MIN_VAL + (double)(proxy->generateCommonRandom() & RAND_MAX) / ((double)(RAND_MAX / (MAX_VAL - MIN_VAL)));
        x[i] = proxy->createShare(xd);
        y[i] = proxy->createShare(yd);
    }
    proxy->SendBytes(CORE_MCMP,params, 1);
    uint64_t *r = CMP(proxy,x,y,sz);
    // checking the result
    uint64_t *x_reconstructed = REC(proxy,x,sz);
    uint64_t *y_reconstructed = REC(proxy,y,sz);
    uint64_t *r_reconstructed = REC(proxy,r,sz);
    bool flag = true;
    for (int i=0;i<sz;i++){
        uint64_t r_computed = (convert2double(x_reconstructed[i])>=convert2double(y_reconstructed[i]))<<FRAC;
        if (r_computed != r_reconstructed[i]){
            flag = false;
            break;
        }
    }
    if (flag)
        cout<<"Vectorized CMP works correctly"<<endl;
    else
        cout<<"Vectorized CMP works incorrectly"<<endl;

}

void MUX_Test(Party *proxy){
    cout<<setfill ('*')<<setw(50)<<"Calling MUX";
    cout<<setfill ('*')<<setw(49)<<"*"<<endl;
    double xd = MIN_VAL + (double)(proxy->generateCommonRandom() & RAND_MAX) / ((double)(RAND_MAX / (MAX_VAL - MIN_VAL)));
    double yd = MIN_VAL + (double)(proxy->generateCommonRandom() & RAND_MAX) / ((double)(RAND_MAX / (MAX_VAL - MIN_VAL)));
    double zd = (double)(proxy->generateCommonRandom()&0x1);
    uint64_t x = proxy->createShare(xd);
    uint64_t y = proxy->createShare(yd);
    uint64_t z = proxy->createShare(zd);

    proxy->SendBytes(CORE_MUX);
    uint64_t r = MUX(proxy,x, y, z);
    // checking the result
    uint64_t x_reconstructed = REC(proxy,x);
    uint64_t y_reconstructed = REC(proxy,y);
    uint64_t z_reconstructed = REC(proxy,z);
    uint64_t r_reconstructed = REC(proxy,r);
    uint64_t r_computed;
    if (z_reconstructed == 0)
        r_computed = x_reconstructed;
    else if (z_reconstructed == (1<<FRAC))
        r_computed = y_reconstructed;
    if (r_reconstructed == r_computed)
        cout<<"MUX works correctly"<<endl;
    else
        cout<<"MUX works incorrectly"<<endl;
}

void MMUX_Test(Party *proxy){
    cout<<setfill ('*')<<setw(50)<<"Calling Vectorized MUX";
    cout<<setfill ('*')<<setw(49)<<"*"<<endl;
    uint64_t x[sz], y[sz], z[sz];

    uint32_t* params = new uint32_t[1];
    params[0] = sz;

    for (int i=0;i<sz;i++){
        double xd = MIN_VAL + (double)(proxy->generateCommonRandom() & RAND_MAX) / ((double)(RAND_MAX / (MAX_VAL - MIN_VAL)));
        double yd = MIN_VAL + (double)(proxy->generateCommonRandom() & RAND_MAX) / ((double)(RAND_MAX / (MAX_VAL - MIN_VAL)));
        double zd = (double)(proxy->generateCommonRandom()&0x1);
        x[i] = proxy->createShare(xd);
        y[i] = proxy->createShare(yd);
        z[i] = proxy->createShare(zd);
    }
    proxy->SendBytes(CORE_MMUX,params, 1);
    uint64_t *r = MUX(proxy,x, y, z, sz);
    // checking the result
    double *x_reconstructed = convert2double(REC(proxy, x, sz), sz);
    double *y_reconstructed = convert2double(REC(proxy, y, sz), sz);
    double *z_reconstructed = convert2double(REC(proxy, z, sz), sz);
    double *r_reconstructed = convert2double(REC(proxy, r, sz), sz);
    bool flag = true;
    double r_computed;
    for (int i=0;i<sz;i++){
        if (z_reconstructed[i] == 0) {
            r_computed = x_reconstructed[i];
        }
        else if (z_reconstructed[i] == 1){
            r_computed = y_reconstructed[i];
        }
        if (r_reconstructed[i] != r_computed){
            flag = false;
            break;
        }
    }
    if (flag)
        cout<<"Vectorized MUX works correctly"<<endl;
    else {
        cout << "Vectorized MUX works incorrectly" << endl;
    }
}

void MAX_Test(Party *proxy){
    cout<<setfill ('*')<<setw(50)<<"Calling MAX";
    cout<<setfill ('*')<<setw(49)<<"*"<<endl;

    uint32_t mRows = WSZ*15;
    uint32_t mCols = WSZ*15;
    uint64_t mSize = mCols*mRows;
    uint32_t* params = new uint32_t[1];
    params[0] = mSize;

    uint64_t *shareOfMatrix = proxy->createShare(random_1D_data(proxy, mSize), mSize);

    proxy->SendBytes(CNN_MAX, params, 1);
    uint64_t max = MAX(proxy, shareOfMatrix, mSize);

    // checking the result
    double computed_max = -1;
    for(uint32_t position = 0; position < mSize; position++){
        double matrixVal = convert2double(REC(proxy, shareOfMatrix[position]));
        if (matrixVal > computed_max){
            computed_max = matrixVal;
        }
    }

    double pp_result = convert2double(REC(proxy, max));
    if(computed_max == pp_result){
        cout<<"MAX works correctly"<<endl;
    }
    else{
        cout<<"MAX works incorrectly" <<endl;
        cout<< "computed: " << pp_result << " should be: " << computed_max << endl;
    }

}

void MMAX_Test(Party *proxy){
    cout<<setfill ('*')<<setw(50)<<"Calling vectorized MAX";
    cout<<setfill ('*')<<setw(49)<<"*"<<endl;
    int precision = 3;

    // INIT PARAMETER
    uint32_t mmaxParams[4];
    mmaxParams[0] = WSZ*3; // matrix Rows
    mmaxParams[1] = WSZ*3; // matrix Columns
    uint64_t mSize = mmaxParams[1] * mmaxParams[0];
    mmaxParams[2] = WSZ; // window rows
    mmaxParams[3] = WSZ; // window cols

    uint64_t *shareOfMatrix = proxy->createShare(random_1D_data(proxy, mSize), mSize);
    // PERFORMING MMAX
    proxy->SendBytes(CNN_MMAX, mmaxParams, 4);
    uint64_t *max = MAX(proxy, shareOfMatrix, mmaxParams[0], mmaxParams[1], mmaxParams[2], mmaxParams[3]);

    // TESTING
    uint32_t window_length = mmaxParams[2] * mmaxParams[3];
    uint32_t number_of_windows = mSize / window_length;
    uint64_t* reconstructed_max = REC(proxy, max, number_of_windows);
    bool flag = true;
    uint64_t resorted [mSize];
    RST(shareOfMatrix, mmaxParams[1], mmaxParams[0], mmaxParams[3], mmaxParams[2], resorted);
    double *d_matrix = convert2double(REC(proxy, resorted, mSize), mSize);
    double computed_max[number_of_windows];

    for(uint32_t win = 0; win < number_of_windows; win++){
        computed_max[win] = d_matrix[window_length*win];
        for(uint32_t win_element = 1; win_element < window_length; win_element++){
            double matrixVal = d_matrix[window_length*win + win_element];
            if (matrixVal > computed_max[win]){
                computed_max[win] = matrixVal;
            }
        }
        if (abs(computed_max[win] - convert2double(reconstructed_max[win]) > 0.001) ){
            flag = false;
        }
    }

    if(flag){
        cout<<"Vectorized MAX works correctly"<<endl;
    }
    else{
        cout<<"Vectorized MAX works incorrectly"<<endl;
        if(mmaxParams[0]*mmaxParams[1] < 200){ //otherwise too many values
            print1DMatrixByWindows("resorted Matrix: ", d_matrix, mmaxParams[0], mmaxParams[1], 1, mmaxParams[2]*mmaxParams[3]);
        }
        print1DMatrixByWindows("computed max values (test): ", computed_max, 1, number_of_windows, 1, 1);
        print1DMatrixByWindows("VS result from method: ", convert2double(reconstructed_max, number_of_windows), 1, number_of_windows, 1, 1);
    }
}

void ARGMAX_Test(Party *proxy){
    cout<<setfill ('*')<<setw(50)<<"Calling ARGMAX";
    cout<<setfill ('*')<<setw(49)<<"*"<<endl;

    uint32_t mRows = WSZ*15;
    uint32_t mCols = WSZ*15;
    uint64_t mSize = mCols*mRows;
    uint32_t* params = new uint32_t[1];
    params[0] = mSize;

    uint64_t *shareOfMatrix = proxy->createShare(random_1D_data(proxy, mSize), mSize);

    proxy->SendBytes(CNN_ARGMAX, params, 1);
    uint64_t argmax = ARGMAX(proxy, shareOfMatrix, mSize);

    // checking the result
    double computed_max = -1;
    double computed_argmax = -1;
    for(uint32_t position = 0; position < mSize; position++){
        double matrixVal = convert2double(REC(proxy, shareOfMatrix[position]));
        if (matrixVal > computed_max){
            computed_max = matrixVal;
            computed_argmax = position;
        }
    }

    double pp_result = convert2double(REC(proxy, argmax));
    if(computed_argmax == pp_result){
        cout<<"ARGMAX works correctly"<<endl;
    }
    else{
        cout<<"ARGMAX works incorrectly" <<endl;
        cout<< "computed: " << pp_result << " should be: " << computed_argmax << endl;
    }

}


void RST_Test(Party *proxy){
    cout<<setfill ('*')<<setw(50)<<"Calling RST";
    cout<<setfill ('*')<<setw(49)<<"*"<<endl;

    uint32_t mRows = WSZ*10;
    uint32_t mCols = WSZ*10;
    uint64_t mSize = mCols*mRows;

    uint32_t wRows = WSZ;
    uint32_t wCols = WSZ;
    uint64_t *shareOfMatrix = proxy->createShare(random_1D_data(proxy, mSize, 1000.0), mSize);
    auto *resorted = new uint64_t [sz];
    //print1DMatrixByWindows("RST original", convert2double(REC(proxy, shareOfMatrix, mSize), mSize), mRows, mCols, wRows, wCols);
    RST(shareOfMatrix, mCols, mRows, wCols, wRows, resorted);

    uint64_t wElements = wRows * wCols;
    uint64_t numberOfWins = mSize / wElements;
    //print1DMatrixByWindows("RST finished: resorted", convert2double(REC(proxy, resorted, mSize), mSize), numberOfWins, wElements, 1, 1);
}

void RELU_Test(Party *proxy){
    cout<<setfill ('*')<<setw(50)<<"Calling RELU";
    cout<<setfill ('*')<<setw(49)<<"*"<<endl;

    uint64_t x = proxy->createShare(convert2double(proxy->generateCommonRandom()));
    proxy->SendBytes(CNN_RELU);
    uint64_t relu = RELU(proxy, x);
    uint64_t reconstructed_relu = REC(proxy, relu);

    // checking the result
    double computed_relu = -1;
    double originalX = convert2double(REC(proxy, x));
    if (originalX >= 0){
        computed_relu = originalX;
    }
    else{
        computed_relu = 0;
    }

    double pp_result = convert2double(reconstructed_relu);
    if(computed_relu == pp_result){
        cout<<"RELU works correctly"<<endl;
    }
    else{
        cout<<"RELU works incorrectly" <<endl;
        cout<< "computed: " << pp_result << " should be: " << computed_relu << endl;
    }

}

void MRELU_Test(Party *proxy){
    cout<<setfill ('*')<<setw(50)<<"Calling MRELU";
    cout<<setfill ('*')<<setw(49)<<"*"<<endl;

    uint64_t size = sz;
    uint32_t* params = new uint32_t[1];
    params[0] = size;
    uint64_t* x = proxy->createShare(random_1D_data(proxy, size, 255.0), size);
    proxy->SendBytes(CNN_MRELU, params, 1);
    uint64_t* mrelu = RELU(proxy, x, size);
    uint64_t* rec_mrelu = REC(proxy, mrelu, size);

    // checking the result
    double* correct_relu = new double [size];
    double* originalX = convert2double(REC(proxy, x, size), size);
    double* pp_result = convert2double(rec_mrelu, size);
    bool foundIncorrect = false;
    for(uint64_t i = 0; i<size; i++){
        if (originalX[i] >= 0){
            correct_relu[i] = originalX[i];
        }
        else{
            correct_relu[i] = 0;
        }
        if(correct_relu[i] != pp_result[i]) {
            foundIncorrect = true;
        }
    }

    if(!foundIncorrect){
        cout<<"MRELU works correctly"<<endl;
    }
    else{
        cout<<"MRELU works incorrectly" <<endl;
        print1DArray("Original Values:", originalX, size);
        print1DArray("Computed RELU:", pp_result, size);
        print1DArray("VS Correct RELU:", correct_relu, size);
    }

}


void DRLU_Test(Party *proxy){
    cout<<setfill ('*')<<setw(50)<<"Calling DRLU";
    cout<<setfill ('*')<<setw(49)<<"*"<<endl;

    uint64_t x = proxy->createShare(convert2double(proxy->generateRandom()));

    proxy->SendBytes(CNN_DRLU);
    uint64_t drelu = DRELU(proxy, x);
    uint64_t reconstructed_drelu = REC(proxy, drelu);

    // checking the result
    double originalX = convert2double(REC(proxy, x));
    uint64_t computed_drelu = 0;
    if (originalX > 0)
        computed_drelu = 1;

    if(computed_drelu == reconstructed_drelu){
        cout<<"DRLU works correctly"<<endl;
    }
    else{
        cout<<"DRLU works incorrectly" <<endl;
        cout<< "computed: " << reconstructed_drelu << " should be: " << computed_drelu << endl;
    }

}

void MDRLU_Test(Party *proxy){
    cout<<setfill ('*')<<setw(50)<<"Calling MDRLU";
    cout<<setfill ('*')<<setw(49)<<"*"<<endl;

    uint64_t size = 10;
    uint32_t* params = new uint32_t[1];
    params[0] = size;
    uint64_t* x = proxy->createShare(random_1D_data(proxy, size), size);

    proxy->SendBytes(CNN_MDRLU, params, 1);
    uint64_t* drelu = DRELU(proxy, x, size);
    double* rec_drelu = convert2double(REC(proxy, drelu, size), size);

    // checking the result
    double* originalX = convert2double(REC(proxy, x, size), size);
    double* correct_drelu = new double [size];
    bool allCorrect = true;
    for(uint64_t i = 0; i<size; i++){
        if (originalX[i] > 0) {
            correct_drelu[i] = 1;
        }
        else{
            correct_drelu[i] = 0;
        }
        if(abs(correct_drelu[i] - rec_drelu[i]) > 0.0001){
            allCorrect = false;
        }
    }
    if(allCorrect){
        cout<<"MDRLU works correctly"<<endl;
    }
    else{
        cout<<"MDRLU works incorrectly" <<endl;
        print1DArray("Original Values:", originalX, size);
        print1DArray("Computed DRELU:", rec_drelu, size);
        print1DArray("VS Correct DRELU:", correct_drelu, size);
    }

}


void PAD_Test(Party *proxy){
    cout<<setfill ('*')<<setw(50)<<"Calling PAD (padding)";
    cout<<setfill ('*')<<setw(49)<<"*"<<endl;

    uint32_t row = 28;
    uint32_t col = 28;
    uint32_t padding_value = 0;
    uint32_t padding_size = 2;
    uint64_t** x = proxy->createShare(random_2D_data(proxy, row, col, 0.0, 255.0), row, col);
    //print2DArray("Original matrix: ", convert2double(REC(proxy, x, row, col), row, col), row, col);
    uint64_t** padded = PAD(x, row, col, padding_value, padding_size);
    double** rec_padded = convert2double(REC(proxy, padded, row+2*padding_size, col+2*padding_size), row+2*padding_size, col+2*padding_size);

    double** computed_padding = new double *[row + 2*padding_size];
    for (uint32_t p = 0; p<padding_size; p++){
        computed_padding[p] = new double [col + 2*padding_size];
        computed_padding[p + padding_size + row] = new double [col + 2*padding_size];

        for(uint32_t c = 0; c<(col + 2*padding_size); c++){
            computed_padding[p][c] = padding_value;                      //top
            computed_padding[p + padding_size + row][c] = padding_value; //bottom
        }
    }
    for(uint64_t r = 0; r < row; r++){
        computed_padding[padding_size + r] = new double [col + 2*padding_size];
        for(uint32_t p = 0; p<padding_size;p++){
            computed_padding[padding_size + r][p] = padding_value;
            computed_padding[padding_size + r][padding_size + col + p] = padding_value;
        }
        for(uint32_t c = 0; c<col; c++) {
            computed_padding[padding_size + r][c + padding_size] = convert2double(REC(proxy, x[r][c]));
        }
    }

    // checking the result
    bool allCorrect = true;
    for (int i = 0; i < (2 * padding_size + row); ++i) {
        for (int j = 0; j < (2 * padding_size + col); ++j) {
            if (abs(computed_padding[i][j] - rec_padded[i][j]) > 0.0001){
                allCorrect = false;
            }
        }
    }
    print2DArray("Computed Padding", rec_padded, row + 2*padding_size, col + 2*padding_size);
    if(allCorrect){
        cout<<"PAD works correctly"<<endl;
    }
    else{
        cout<<"PAD works incorrectly" <<endl;
        print2DArray("Computed Padding", rec_padded, row + 2*padding_size, col + 2*padding_size);
        print2DArray("Correct Padding", computed_padding, row + 2*padding_size, col + 2*padding_size);
    }

}

void CL_Test(Party *proxy){
    cout<<setfill ('*')<<setw(50)<<"Calling CL (convolutional layer)";
    cout<<setfill ('*')<<setw(49)<<"*"<<endl;
                                                                                            // set params
    uint32_t row = 6;
    uint32_t col = 6;
    uint32_t i_channel = 1;
    uint32_t k_number = 1;
    uint32_t k_dim = 3;
    uint32_t stride = 1;
    uint32_t maxpool_window_dim = 0;
                                                                                            // generate secret shares of random data
    uint64_t*** x = new uint64_t** [i_channel];
    for(uint64_t i = 0; i<i_channel; i++){
        x[i] = proxy->createShare(random_2D_data(proxy, row, col, 0, 255), row, col);
    }
    uint64_t*** kernel = new uint64_t **[k_number];
    for(uint32_t o = 0; o<k_number; o++){
        kernel[o] = proxy->createShare(random_2D_data(proxy, i_channel, k_dim * k_dim, -5, 5), i_channel, k_dim * k_dim);
        double** rec_kernel = convert2double(REC(proxy, kernel[o], i_channel, k_dim * k_dim), i_channel, k_dim * k_dim);
        print2DArray("recon kernel: ", rec_kernel, i_channel, k_dim * k_dim);
    }
    uint64_t *bias = proxy->createShare(random_1D_data(proxy, k_number), k_number);
    double* rec_bias = convert2double(REC(proxy, bias, k_number), k_number);
    print1DArray("recon bias: ", rec_bias, k_number);
    delete[] bias;
                                                                                            // send params
    uint32_t params[7];
    params[0] = i_channel;
    params[1] = row;
    params[2] = col;
    params[3] = k_dim;      // kernel size
    params[4] = k_number;   // kernel number = output channel
    params[5] = stride;
    params[6] = maxpool_window_dim;
    params[7] = false;
    proxy->SendBytes(CNN_CL, params, 8);

    uint64_t*** conv = CL(proxy, x, i_channel, row, col, kernel, k_dim, k_number, stride, maxpool_window_dim, bias, false);
    uint64_t conv_size = floor((row - k_dim + 1) / stride);
    uint64_t lastpos = row - k_dim + 1;
    uint64_t out_size = conv_size;
    if(maxpool_window_dim > 0){
        out_size = conv_size/maxpool_window_dim;
    }

    // checking the result
    double*** conv_unreduced = new double **[k_number]; // no reduction with MAX
    double*** correct_conv = new double **[k_number];
    bool allCorrect = true;
    for(uint32_t kern = 0; kern < k_number; kern++){    // for each output channel
        conv_unreduced[kern] = new double *[conv_size];                // init kernels conv result
        correct_conv[kern] = new double *[out_size];
        for(uint32_t cr = 0; cr < lastpos; cr+=stride){
            conv_unreduced[kern][cr / stride] = new double [conv_size];   // init row of conv result
            if(maxpool_window_dim <= 0){
                correct_conv[kern][cr / stride] = new double [out_size];
            }
            for(uint32_t cc = 0; cc < lastpos; cc+=stride){
                double dot_product = 0;
                for(uint32_t channel = 0; channel<i_channel; channel++){
                    for(uint32_t k_value = 0; k_value<(k_dim*k_dim); k_value++){
                        double v = convert2double(REC(proxy, x[channel][cr + k_value / k_dim][cc + k_value % k_dim]));
                        double weight = convert2double(REC(proxy, kernel[kern][channel][k_value]));
                        dot_product += v * weight;
                    }
                }
                // Activation: RELU
                double relu = 0.0;
                if(dot_product > 0) {
                    relu = dot_product;
                }
                relu += rec_bias[kern];
                conv_unreduced[kern][cr / stride][cc / stride] = relu; // does not matter if added here or after maxpool
                if(maxpool_window_dim <= 0){
                    correct_conv[kern][cr / stride][cc / stride] = relu;
                    if(abs(relu - convert2double(REC(proxy, conv[kern][cr / stride][cc / stride]))) > 0.0001) {
                        allCorrect = false;
                    }
                }
            }
        }
        if(maxpool_window_dim > 0){
            // MAXPOOL reduction for window size 2x2
            for(uint32_t mr = 0; mr < out_size; mr++){
                correct_conv[kern][mr] = new double [out_size];
                for(uint32_t mc = 0; mc < out_size; mc++){
                    uint32_t conv_r = mr * 2, conv_c = mc * 2;
                    double a = std::max(conv_unreduced[kern][conv_r][conv_c], conv_unreduced[kern][conv_r][conv_c+1]);
                    double b = std::max(conv_unreduced[kern][conv_r+1][conv_c], conv_unreduced[kern][conv_r+1][conv_c+1]);
                    double m = std::max(a,b);
                    correct_conv[kern][mr][mc] = m;
                    if (abs(m - convert2double(REC(proxy, conv[kern][mr][mc]))) > 0.0001){
                        allCorrect = false;
                    }
                }
            }
        }
    }
    if(allCorrect){
        cout<<"CL works correctly"<<endl;
    }
    else{
        cout<<"CL works incorrectly" <<endl;
        print2DArray("First channel of x", convert2double(REC(proxy, x[0], row, col), row, col), row, col);
        print1DArray("kernel 0", convert2double(REC(proxy, kernel[0][0], k_dim*k_dim), k_dim*k_dim), k_dim*k_dim);
        for (uint64_t k = 0; k<k_number; k++){
            cout << "Channel " << k << endl;
            print2DArray("Channel of computed convolution", convert2double(REC(proxy, conv[k], out_size, out_size), out_size, out_size), out_size, out_size);
            print2DArray("Channel of correct convolution", correct_conv[k], out_size, out_size);
        }
    }

}


void FCL_Test(Party *proxy){
    cout<<setfill ('*')<<setw(50)<<"Calling FCL (fully connected layer)";
    cout<<setfill ('*')<<setw(49)<<"*"<<endl;

    // init
    uint32_t row = 2;
    uint32_t col = 2;
    uint32_t i_channel = 1;
    uint32_t i_nodes = row*col*i_channel;
    uint32_t o_nodes = 2;

    uint64_t* x = proxy->createShare(random_1D_data(proxy, i_nodes, 0.0, 255.0), i_nodes);
    uint64_t** weights = proxy->createShare(random_2D_data(proxy, i_nodes, o_nodes, -5.0, 5.0), i_nodes, o_nodes);
    uint64_t *bias = proxy->createShare(random_1D_data(proxy, o_nodes), o_nodes);

    // send params
    uint32_t params[2];
    params[0] = i_nodes;
    params[1] = o_nodes;
    proxy->SendBytes(CNN_FCL, params, 2);

    uint64_t *res = FCL(proxy, x, i_nodes, weights, o_nodes, bias);
    double* reconstructed_res = convert2double(REC(proxy, res, o_nodes), o_nodes);

    // checking the result
    bool allCorrect = true;
    double *rec_input = convert2double(REC(proxy, x, i_nodes), i_nodes);
    print1DArray("reconstructed input: ", rec_input, i_nodes);
    double **rec_weights = convert2double(REC(proxy, weights, i_nodes, o_nodes), i_nodes, o_nodes);
    print2DArray("reconstructed weights: ", rec_weights, i_nodes, o_nodes);
    double* rec_bias = convert2double(REC(proxy, bias, o_nodes), o_nodes);
    print1DArray("reconstructed bias: ", rec_bias, o_nodes);

    double * original_res = new double [o_nodes];
    for(uint64_t o = 0; o<o_nodes; o++){
        double value = 0;
        for (int i = 0; i < i_nodes; i++) {
            value += rec_input[i] * rec_weights[o][i]; //weights would be tansposed in FCL
        }
        //ReLU activation
        if(value < 0){
            value = 0;
        }
        original_res[o] = value + rec_bias[o];
        if (original_res[o] - reconstructed_res[o] > 0.0001){
            allCorrect = false;
        }
    }

    if(allCorrect){
        cout<<"FCL works correctly"<<endl;
    }
    else{
        cout<<"FCL works incorrectly" <<endl;
        print1DArray("Computed Result: ", reconstructed_res, o_nodes);
        print1DArray("Correct Result", original_res, o_nodes);
    }
}

void FLT_Test(Party *proxy){
    cout<<setfill ('*')<<setw(50)<<"Calling FLT (flattening)";
    cout<<setfill ('*')<<setw(49)<<"*"<<endl;

    uint32_t row = 5;
    uint32_t col = 5;
    uint32_t i_channel = 3;

    uint64_t*** x = new uint64_t **[i_channel];
    for(uint64_t i = 0; i<i_channel; i++){
        x[i] = proxy->createShare(random_2D_data(proxy, row, col, 0, 255), row, col);
    }
    uint32_t i_nodes = row*col*i_channel;

    uint64_t *res = FLT(x, row, col, i_channel);
    double* reconstructed_res = convert2double(REC(proxy, res, i_nodes), i_nodes);

    // checking the result
    bool allCorrect = true;
    double *original_res = new double [i_nodes];
    for (int channel = 0; channel < i_channel; ++channel) {
        for (int r = 0; r < row; ++r) {
            for (int c = 0; c < col; ++c) {
                uint64_t position = channel*row*col + r*col + c;
                original_res[position] = convert2double(REC(proxy, x[channel][r][c]));
                if ((original_res[position] - reconstructed_res[position]) > 0.0001){
                    allCorrect = false;
                }
            }
        }
    }

    if(allCorrect){
        cout<<"FLT works correctly"<<endl;
    }
    else{
        cout<<"FLT works incorrectly" <<endl;
        print1DArray("Computed Result: ", reconstructed_res, i_nodes);
        print1DArray("Correct Result", original_res, i_nodes);
    }
}

void INC_Test(Party *proxy){
    cout<<setfill ('*')<<setw(50)<<"Calling INC (increasing input matrix for conv)";
    cout<<setfill ('*')<<setw(49)<<"*"<<endl;

    uint32_t row = 28;
    uint32_t col = 28;
    uint32_t i_channel = 1;
    uint64_t *** x = new uint64_t **[i_channel];
    for(uint64_t i = 0; i<i_channel; i++){
        x[i] = proxy->createShare(random_2D_data(proxy, row, col, 0, 255), row, col);
    }
    //print2DArray("Original X: ", convert2double(REC(proxy, x[0], row, col), row, col), row, col);

    uint32_t k_number = 5;
    uint32_t k_dim = 5;
    uint32_t k_len = k_dim*k_dim;
    uint32_t stride = 1;

    uint32_t lastpos = row - k_dim + 1;
    uint32_t conv_height = static_cast<uint32_t>(floor(lastpos / stride));
    uint32_t conv_width = static_cast<uint32_t>(floor((col - k_dim + 1) / stride));

    uint64_t*** stretchedX = INC(x, i_channel, row, col, k_dim, stride);
    // checking the result
    print2DArray("RESULTED STRETCHED MATRIX (first channel)", convert2double(REC(proxy, stretchedX[0], conv_height*conv_width, k_len), conv_height*conv_width, k_len), conv_height*conv_width, k_len);
}

void EXP_Test(Party *proxy){
    cout<<setfill ('*')<<setw(50)<<"Calling EXP";
    cout<<setfill ('*')<<setw(49)<<"*"<<endl;

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

    if(abs(true_exp - rec_exp) <= 0.1 ){
        cout<<"EXP works correctly"<<endl;
        cout<< "power: " << originalX << " -- computed: " << rec_exp << " should be: " << true_exp << endl;
    }
    else{
        cout<<"EXP works incorrectly" <<endl;
        cout<< "power: " << originalX << " -- computed: " << rec_exp << " should be: " << true_exp << endl;
    }

}

void MEXP_Test(Party *proxy){
    cout<<setfill ('*')<<setw(50)<<"Calling MEXP";
    cout<<setfill ('*')<<setw(49)<<"*"<<endl;

    uint32_t params[1];
    int n_samples = 10;
    params[0] = n_samples;
    uint64_t *x = proxy->createShare(random_1D_data(proxy, n_samples, proxy->getMinPower(), proxy->getMaxPower()), n_samples);
    proxy->SendBytes(CORE_MEXP, params, 1);
    uint64_t* shr_exp = EXP(proxy, x, n_samples);
    uint64_t* reconstructed_exp = REC(proxy, shr_exp, n_samples);
    double* rec_exp = convert2double(reconstructed_exp, n_samples);

    // checking the result
    double* originalX = convert2double(REC(proxy, x, n_samples), n_samples);
    double* true_exp = new double[n_samples];
    for(int i = 0; i < n_samples; i++) {
        true_exp[i] = exp(originalX[i]);
    }

    bool flag = true;
    for(int i = 0; i < n_samples; i++) {
        cout<< fixed << "power: " << originalX[i] << " -- computed: " << rec_exp[i] << " - expected: " << true_exp[i] <<
        " - absolute difference: " << abs(rec_exp[i] - true_exp[i])<< endl;
        if(abs(true_exp[i] - rec_exp[i]) >= 0.1 ){
            flag = false;
        }
    }
    if(flag){
        cout<<"EXP works correctly"<<endl;
    }
    else{
        cout<<"EXP works incorrectly" <<endl;
    }

}

void DP_Test(Party* proxy) {
    cout<<setfill ('*')<<setw(50)<<"Calling DP";
    cout<<setfill ('*')<<setw(49)<<"*"<<endl;

    uint32_t* params = new uint32_t[1];
    uint32_t size = 20; // size of the vector
    params[0] = size;
    double min_val = -10;
    double max_val = 10;

    // generate a vector of random values
    uint64_t* vec1 = proxy->createShare(random_1D_data(proxy, size, min_val, max_val), size);
    uint64_t* vec2 = proxy->createShare(random_1D_data(proxy, size, min_val, max_val), size);

    // call DP
    proxy->SendBytes(CORE_DP, params, 1);
    uint64_t res = DP(proxy, vec1, vec2, size);

    // check the results
    double gt = 0;
    double* rec_vec1 = convert2double(REC(proxy, vec1, size), size);
    double* rec_vec2 = convert2double(REC(proxy, vec2, size), size);
    for(int i = 0; i < size; i++) {
        gt += rec_vec1[i] * rec_vec2[i];
    }

    double rec_res = convert2double(REC(proxy, res));

    printValue("Computed dot product", rec_res);
    printValue("GT dot product", gt);
}

void MDP_Test(Party* proxy) {
    cout<<setfill ('*')<<setw(50)<<"Calling MDP";
    cout<<setfill ('*')<<setw(49)<<"*"<<endl;

    uint32_t* params = new uint32_t[1];
    uint32_t size = 20; // size of the total vector
    params[0] = size;

    uint32_t d = 5; // size of each individual vector in the main vector
    double min_val = -10;
    double max_val = 10;

    if(size % d != 0) {
        throw invalid_argument("DP_Test: The size must be divisible by d.");
    }

    // generate a vector of random values
    uint64_t* vec1 = proxy->createShare(random_1D_data(proxy, size, min_val, max_val), size);
    uint64_t* vec2 = proxy->createShare(random_1D_data(proxy, size, min_val, max_val), size);

    // call DP
    proxy->SendBytes(CORE_MDP, params, 1);
    uint64_t* res = DP(proxy, vec1, vec2, size, d);

    // check the results
    double* gt = new double[size / d];
    double* rec_vec1 = convert2double(REC(proxy, vec1, size), size);
    double* rec_vec2 = convert2double(REC(proxy, vec2, size), size);
    for(int i = 0; i < size; i += d) {
        double tmp_sum = 0;
        for(int j = i; j < i + d; j++) {
            tmp_sum += rec_vec1[j] * rec_vec2[j];
        }
        gt[i / d] = tmp_sum;
    }

    double *rec_res = convert2double(REC(proxy, res, size / d), size / d);

    print1DArray("Computed dot product", rec_res, size / d);
    print1DArray("GT dot product", gt, size / d);
}

void MATMATMUL_Test(Party* proxy) {
    cout<<setfill ('*')<<setw(50)<<"Calling MATMATMUL";
    cout<<setfill ('*')<<setw(49)<<"*"<<endl;

    // setting
    int a_row = 5;
    int a_col = 6;
    int b_col = 3;

    uint32_t* params = new uint32_t[1];
    uint32_t size = a_row * a_col * b_col; // size of the vector
    params[0] = size;
    double min_val = -10;
    double max_val = 10;

    uint64_t** mat1 = proxy->createShare(random_2D_data(proxy, a_row, a_col, min_val, max_val), a_row, a_col);
    uint64_t** mat2 = proxy->createShare( random_2D_data(proxy, a_col, b_col, min_val, max_val), a_col, b_col);

    proxy->SendBytes(CORE_MATMATMUL, params, 1);
    uint64_t** res = MATMATMUL(proxy, mat1, mat2, a_row, a_col, b_col);
    double** rec_res = convert2double(REC(proxy, res, a_row, b_col), a_row, b_col);

    double** rec_mat1 = convert2double(REC(proxy, mat1, a_row, a_col), a_row, a_col);
    double** rec_mat2 = convert2double(REC(proxy, mat2, a_col, b_col), a_col, b_col);

//    proxy->print2DArray("mat1", rec_mat1, a_row, a_col);
//    proxy->print2DArray("mat2", rec_mat2, a_col, b_col);

    double** gt = multiply_matrices(rec_mat1, rec_mat2, a_row, a_col, b_col);

    double tmp = 0;
    for( int i = 0; i < a_row; i++) {
        for(int j = 0; j < b_col; j++) {
            tmp += abs(rec_res[i][j] - gt[i][j]);
        }
    }

    if(tmp <= 0.1) {
        cout << "MATMATMUL works correctly" << endl;
        cout << "Total absolute difference: " << tmp << endl;
    }
    else {
        cout << "MATMATMUL works incorrectly" << endl;
        cout << "Total absolute difference: " << tmp << endl;
    }

//    print2DArray("Computed matrix multiplication", rec_res, a_row, b_col);
//    print2DArray("GT matrix multiplication", gt, a_row, b_col);
}

void MMATMATMUL_Test(Party *proxy) {
    // In this function, we test several matrix multiplications of random matrices with the same size.
    cout<<setfill ('*')<<setw(50)<<"Calling MMATMATMUL";
    cout<<setfill ('*')<<setw(49)<<"*"<<endl;

    // setting
    int n_matrices = 3;
    int a_row = 5;
    int a_col = 6;
    int b_col = 3;
    uint32_t* params = new uint32_t[1];
    uint32_t size = n_matrices * a_row * a_col * b_col; // size of the vector
    params[0] = size;
    double min_val = -10;
    double max_val = 10;

    uint64_t*** mat1 = new uint64_t**[n_matrices];
    uint64_t*** mat2 = new uint64_t**[n_matrices];
    for(int i = 0; i < n_matrices; i++) {
        mat1[i] = proxy->createShare(random_2D_data(proxy, a_row, a_col, min_val, max_val), a_row, a_col);
        mat2[i] = proxy->createShare(random_2D_data(proxy, a_col, b_col, min_val, max_val), a_col, b_col);
    }

    proxy->SendBytes(CORE_MMATMATMUL, params, 1);
    uint64_t*** res = MATMATMUL(proxy, mat1, mat2, n_matrices, a_row, a_col, b_col);

    double*** rec_mat1 = new double**[n_matrices];
    double*** rec_mat2 = new double**[n_matrices];
    for(int i = 0; i < n_matrices; i++) {
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

    double* diffs = new double[n_matrices];
    bool flag = true;
    for(int m = 0; m < n_matrices; m++) {
        double** rec_res = convert2double(REC(proxy, res[m], a_row, b_col), a_row, b_col);
        double** gt = multiply_matrices(rec_mat1[m], rec_mat2[m], a_row, a_col, b_col);

        diffs[m] = 0;
        for( int i = 0; i < a_row; i++) {
            for(int j = 0; j < b_col; j++) {
                diffs[m] += abs(rec_res[i][j] - gt[i][j]);
            }
        }

        if(diffs[m] >= 0.1) {
            flag = false;
        }
    }

    if(flag) {
        cout << "MMATMATMUL works correctly" << endl;
        cout << "Total absolute differences: ";
        for(int m = 0; m < n_matrices; m++) {
            cout << diffs[m] << "\t";
        }
        cout << endl;
    }
    else {
        cout << "MMATMATMUL works incorrectly" << endl;
        cout << "Total absolute differences: ";
        for(int m = 0; m < n_matrices; m++) {
            cout << diffs[m] << "\t";
        }
        cout << endl;
    }
}

void MATVECMUL_Test(Party *proxy) {
    /* In this function, we test the matrix multiplication of two random matrices. We first generate two
     * matrices of random values such that the number of column of the first matrix equals to the number
     * of rows of the second matrix.
     */
    cout<<setfill ('*')<<setw(50)<<"Calling MATVECMUL";
    cout<<setfill ('*')<<setw(49)<<"*"<<endl;

    // setting
    int a_row = 5;
    int a_col = 6;
    int b_col = 3;
    uint32_t* params = new uint32_t[1];
    uint32_t size = a_row * a_col; // size of the vector
    params[0] = size;
    double min_val = -10;
    double max_val = 10;

    uint64_t** mat = proxy->createShare(random_2D_data(proxy, a_row, a_col, min_val, max_val), a_row, a_col);
    uint64_t* vec = proxy->createShare(random_1D_data(proxy, a_col, min_val, max_val), a_col);

    proxy->SendBytes(CORE_MATVECMUL, params, 1);
    uint64_t* res = MATVECMUL(proxy, mat, vec, a_row, a_col);
    double* rec_res = convert2double(REC(proxy, res, a_row), a_row);

    double** rec_mat = convert2double(REC(proxy, mat, a_row, a_col), a_row, a_col);
    double* rec_vec = convert2double(REC(proxy, vec, a_col), a_col);

//    proxy->print2DArray("mat", rec_mat, a_row, a_col);
//    proxy->print1DArray("vec", rec_vec, a_col);

    double* gt = multiply_matrice_vector(rec_mat, rec_vec, a_row, a_col);

    double tmp = 0;
    for( int i = 0; i < a_row; i++) {
        tmp += abs(rec_res[i] - gt[i]);
    }

    if(tmp <= 0.1) {
        cout << "MATVECMUL works correctly" << endl;
        cout << "Total absolute difference: " << tmp << endl;
    }
    else {
        cout << "MATVECMUL works incorrectly" << endl;
        cout << "Total absolute difference: " << tmp << endl;
    }

//    print1DArray("Computed matrix-vector multiplication", convert2double(REC(proxy, res, a_row), a_row), a_row);
//    print1DArray("GT matrix-vector multiplication", gt, a_row);
}

void MMATVECMUL_Test(Party *proxy) {
    // In this function, we test the matrix multiplication of n_matrices matrices and vectors.
    cout<<setfill ('*')<<setw(50)<<"Calling MMATVECMUL";
    cout<<setfill ('*')<<setw(49)<<"*"<<endl;

    // setting
    int n_matrices = 3;
    int a_row = 5;
    int a_col = 6;
    int b_col = 3;
    uint32_t* params = new uint32_t[1];
    uint32_t size = n_matrices * a_row * a_col; // size of the vector
    params[0] = size;
    double min_val = -10;
    double max_val = 10;

    uint64_t*** mat = new uint64_t**[n_matrices];
    uint64_t** vec = new uint64_t*[n_matrices];
    for(int i = 0; i < n_matrices; i++) {
        mat[i] = proxy->createShare(random_2D_data(proxy, a_row, a_col, min_val, max_val), a_row, a_col);
        vec[i] = proxy->createShare(random_1D_data(proxy, a_col, min_val, max_val), a_col);
    }

    proxy->SendBytes(CORE_MMATVECMUL, params, 1);
    uint64_t** res = MATVECMUL(proxy, mat, vec, n_matrices, a_row, a_col);

    double*** rec_mat = new double**[n_matrices];
    double** rec_vec = new double*[n_matrices];
    for(int i = 0; i < n_matrices; i++) {
        rec_mat[i] = convert2double(REC(proxy, mat[i], a_row, a_col), a_row, a_col);
        rec_vec[i] = convert2double(REC(proxy, vec[i], a_col), a_col);
    }

//    for(int i = 0; i < n_matrices; i++) {
//        print1DArray("Computed matrix-vector multiplication - " + to_string(i), convert2double(REC(proxy, res[i], a_row), a_row), a_row);
//        double* gt = multiply_matrice_vector(rec_mat[i], rec_vec[i], a_row, a_col);
//        print1DArray("Ground truth matrix-vector multiplication - " + to_string(i), gt, a_row);
//    }

    double* diffs = new double[n_matrices];
    bool flag = true;
    for(int m = 0; m < n_matrices; m++) {
        double* rec_res = convert2double(REC(proxy, res[m], a_row), a_row);
        double* gt = multiply_matrice_vector(rec_mat[m], rec_vec[m], a_row, a_col);

        diffs[m] = 0;
        for( int i = 0; i < a_row; i++) {
            diffs[m] += abs(rec_res[i] - gt[i]);
        }

        if(diffs[m] >= 0.1) {
            flag = false;
        }
    }

    if(flag) {
        cout << "MMATVECMUL works correctly" << endl;
        cout << "Total absolute difference: " << endl;
        for(int m = 0; m < n_matrices; m++) {
            cout << diffs[m] << "\t";
        }
        cout << endl;
    }
    else {
        cout << "MMATVECMUL works incorrectly" << endl;
        cout << "Total absolute difference: " << endl;
        for(int m = 0; m < n_matrices; m++) {
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

void DIV_Test(Party *proxy){
    cout<<setfill ('*')<<setw(50)<<"Calling DIV";
    cout<<setfill ('*')<<setw(49)<<"*"<<endl;

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


    if(computed_div == reconstructed_div){
        cout<<"DIV works correctly"<<endl;
    }
    else{
        cout<<"DIV works incorrectly" <<endl;
        cout<< "computed: " << reconstructed_div << " (" << bitset<L_BIT>(reconstructed_div) << ") but should be: " << computed_div << " (" << bitset<L_BIT>(computed_div) << ")" << endl;
    }

}


void ADD_Test(Party *proxy){
    cout<<setfill ('*')<<setw(50)<<"Calling ADD for one 2D vector";
    cout<<setfill ('*')<<setw(49)<<"*"<<endl;

    uint64_t **x = proxy->createShare(random_2D_data(proxy, 3, 10, 0, 5), 3, 10);
    uint64_t *res = ADD(proxy, x, 3, 10);

    double *recon_res = convert2double(REC(proxy, res, 10), 10);

    // checking the result
    double **originalX = convert2double(REC(proxy, x, 3, 10), 3, 10);
    double *correct_sum = new double [10];
    print2DArray("ORIGINAL X: ", originalX, 3, 10);
    bool allCorrect = true;
    for(uint64_t v = 0; v<10; v++){
        uint64_t s = 0;
        for(uint64_t c = 0; c<3; c++){
            s += originalX[c][v];
        }
        correct_sum[v] = s;
        if((s - recon_res[v]) > 0.0001){
            allCorrect = false;
        }
    }
    if(allCorrect){
        cout<<"ADD works correctly"<<endl;
    }
    else{
        cout<<"ADD works incorrectly" <<endl;
        print1DArray("Computed SUM: ", recon_res, 10);
        print1DArray("Correct SUM: ", correct_sum, 10);
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



// Main function to run the experiments
int main(int argc, char* argv[]) {
    uint8_t role = atoi(argv[1]);
    uint16_t cport = atoi(argv[2]);
    string caddress(argv[3]);
    uint16_t hport = atoi(argv[4]);
    string haddress(argv[5]);

    Party *proxy;
    if (role==0)
        proxy = new Party(P1,hport, haddress, cport, caddress);
    else
        proxy = new Party(P2,hport, haddress, cport, caddress);

    /*
    ADD_Test(proxy);
    MUL_Test(proxy);
    MMUL_Test(proxy);

//    MOC_Test(proxy);
//    MMOC_Test(proxy);

    MSB_Test(proxy);
    MMSB_Test(proxy);

//    CMP_Test(proxy);
//    MCMP_Test(proxy);

//    MUX_Test(proxy);
//    MMUX_Test(proxy);
*/
    MAX_Test(proxy);
    MMAX_Test(proxy);/*

    RST_Test(proxy);
    RELU_Test(proxy);
    MRELU_Test(proxy);

    DRLU_Test(proxy);
    ARGMAX_Test(proxy);*/
    //MDRLU_Test(proxy);//TODO
    //DIV_Test(proxy);
/*
    INC_Test(proxy);
    FLT_Test(proxy);
    FCL_Test(proxy);
    PAD_Test(proxy);
    CL_Test(proxy);

    EXP_Test(proxy);
    MEXP_Test(proxy);

    DP_Test(proxy);
    MDP_Test(proxy);
    MATMATMUL_Test(proxy);
    MMATMATMUL_Test(proxy);

    MATVECMUL_Test(proxy);*/

//    INVSQRT_Test(proxy);
//    MINVSQRT_Test(proxy);

//    ppRKN_ITER_Test(proxy);
//    ppRKN_PREDICTION_Test(proxy);

    proxy->SendBytes(CORE_END);
    proxy->PrintBytes();
    return 0;
}
