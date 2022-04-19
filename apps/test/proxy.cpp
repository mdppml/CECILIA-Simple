#include <cstdlib>
#include <iostream>
#include <deque>
#include <chrono>
#include <tuple>
#include <iomanip>
#include "../../core/core.h"
#include <bitset>

#include "../../core/cnn.h"
#include "../../core/rkn.h"
#include "../../utils/flib.h"
//#include <Eigen/Eigen>

using namespace std;
//using namespace Eigen;

constexpr int MIN_VAL = -100;
constexpr int MAX_VAL = 100;
constexpr int sz = 1000;
constexpr int WSZ = 4;
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

    for (int i=0;i<sz;i++){
        double xd = MIN_VAL + (double)(proxy->generateCommonRandom() & RAND_MAX) / ((double)(RAND_MAX / (MAX_VAL - MIN_VAL)));
        double yd = MIN_VAL + (double)(proxy->generateCommonRandom() & RAND_MAX) / ((double)(RAND_MAX / (MAX_VAL - MIN_VAL)));
        x[i] = proxy->createShare(xd);
        y[i] = proxy->createShare(yd);
    }
    proxy->SendBytes(CORE_MMUL,sz);
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
    for (int i=0;i<sz;i++)
        x[i] = proxy->generateRandom()&N1_MASK;
    proxy->SendBytes(CORE_MMC,sz);
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

void MSB_Test(Party *proxy){
    cout<<setfill ('*')<<setw(50)<<"Calling MSB";
    cout<<setfill ('*')<<setw(49)<<"*"<<endl;
    uint64_t x = proxy->generateRandom();
    proxy->SendBytes(CORE_MSB);
    uint64_t r = MSB(proxy,x);
    // checking the result
    uint64_t x_reconstructed = REC(proxy,x);
    uint64_t r_reconstructed = REC(proxy,r);
    uint64_t r_computed = (x_reconstructed>>(L_BIT - 1)) << FRAC;
    if (r_reconstructed == r_computed)
        cout<<"MSB works correctly"<<endl;
    else
        cout<<"MSB works incorrectly"<<endl;
}

void MMSB_Test(Party *proxy){
    cout<<setfill ('*')<<setw(50)<<"Calling Vectorized MSB";
    cout<<setfill ('*')<<setw(49)<<"*"<<endl;
    uint64_t x[sz];
    for (int i=0;i<sz;i++)
        x[i] = proxy->generateRandom();
    proxy->SendBytes(CORE_MMSB,sz);
    uint64_t *r = MSB(proxy,x,sz);
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
    if (flag)
        cout<<"Vectorized MSB works correctly"<<endl;
    else
        cout<<"Vectorized MSB works incorrectly"<<endl;
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
    for (int i=0;i<sz;i++){
        double xd = MIN_VAL + (double)(proxy->generateCommonRandom() & RAND_MAX) / ((double)(RAND_MAX / (MAX_VAL - MIN_VAL)));
        double yd = MIN_VAL + (double)(proxy->generateCommonRandom() & RAND_MAX) / ((double)(RAND_MAX / (MAX_VAL - MIN_VAL)));
        x[i] = proxy->createShare(xd);
        y[i] = proxy->createShare(yd);
    }
    proxy->SendBytes(CORE_MCMP,sz);
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
    for (int i=0;i<sz;i++){
        double xd = MIN_VAL + (double)(proxy->generateCommonRandom() & RAND_MAX) / ((double)(RAND_MAX / (MAX_VAL - MIN_VAL)));
        double yd = MIN_VAL + (double)(proxy->generateCommonRandom() & RAND_MAX) / ((double)(RAND_MAX / (MAX_VAL - MIN_VAL)));
        double zd = (double)(proxy->generateCommonRandom()&0x1);
        x[i] = proxy->createShare(xd);
        y[i] = proxy->createShare(yd);
        z[i] = proxy->createShare(zd);
    }
    proxy->SendBytes(CORE_MMUX,sz);
    uint64_t *r = MUX(proxy,x, y, z, sz);
    // checking the result
    uint64_t *x_reconstructed = REC(proxy,x,sz);
    uint64_t *y_reconstructed = REC(proxy,y,sz);
    uint64_t *z_reconstructed = REC(proxy,z,sz);
    uint64_t *r_reconstructed = REC(proxy,r,sz);
    bool flag = true;
    uint64_t r_computed;
    for (int i=0;i<sz;i++){
        if (z_reconstructed[i] == 0)
            r_computed = x_reconstructed[i];
        else if (z_reconstructed[i] == (1<<FRAC))
            r_computed = y_reconstructed[i];
        if (r_reconstructed[i] != r_computed){
            flag = false;
            break;
        }
    }
    if (flag)
        cout<<"Vectorized MUX works correctly"<<endl;
    else
        cout<<"Vectorized MUX works incorrectly"<<endl;

}

void MAX_Test(Party *proxy){
    cout<<setfill ('*')<<setw(50)<<"Calling MAX";
    cout<<setfill ('*')<<setw(49)<<"*"<<endl;

    uint32_t mRows = WSZ*15;
    uint32_t mCols = WSZ*15;
    uint64_t mSize = mCols*mRows;

    uint64_t *shareOfMatrix = proxy->createShare(random_1D_data(proxy, mSize), mSize);

    proxy->SendBytes(CNN_MAX, mSize);
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

void RST_Test(Party *proxy){
    cout<<setfill ('*')<<setw(50)<<"Calling RST";
    cout<<setfill ('*')<<setw(49)<<"*"<<endl;

    uint32_t mRows = 6;
    uint32_t mCols = 4;
    uint64_t mSize = mCols*mRows;

    uint32_t wRows = 2;
    uint32_t wCols = 2;
    uint64_t *shareOfMatrix = proxy->createShare(random_1D_data(proxy, mSize, 32), mSize);
    auto *resorted = new uint64_t [sz];
    print1DMatrixByWindows("RST original", convert2double(REC(proxy, shareOfMatrix, mSize), mSize), mRows, mCols, wRows, wCols);
    RST(shareOfMatrix, mCols, mRows, wCols, wRows, resorted);

    uint64_t wElements = wRows * wCols;
    uint64_t numberOfWins = mSize / wElements;
    print1DMatrixByWindows("RST finished: resorted", convert2double(REC(proxy, resorted, mSize), mSize), numberOfWins, wElements, 1, 1);
}

void MMAX_Test(Party *proxy){
    cout<<setfill ('*')<<setw(50)<<"Calling vectorized MAX";
    cout<<setfill ('*')<<setw(49)<<"*"<<endl;
    int precision = 3;

    // INIT PARAMETER
    uint64_t mmaxParams[4];
    mmaxParams[0] = WSZ*3; // matrix Rows
    mmaxParams[1] = WSZ*3; // matrix Columns
    uint64_t mSize = mmaxParams[1] * mmaxParams[0];

    mmaxParams[2] = WSZ; // window rows
    mmaxParams[3] = WSZ;  // window columns

    uint64_t *shareOfMatrix = proxy->createShare(random_1D_data(proxy, mSize), mSize);

    // PERFORMING MMAX
    proxy->SendBytes(CNN_MMAX);

    unsigned char *ptr_out = proxy->getBuffer1();
    addVal2CharArray(mmaxParams, &ptr_out, 4);
    Send(proxy->getSocketHelper(), proxy->getBuffer1(), 4 * 8);

    uint64_t *max = MAX(proxy, shareOfMatrix, mmaxParams[0], mmaxParams[1], mmaxParams[3]);

    // TESTING
    uint32_t window_length = mmaxParams[2] * mmaxParams[3];
    uint32_t number_of_windows = mSize / window_length;
    uint64_t* reconstructed_max = REC(proxy, max, number_of_windows);

    bool flag = true;
    uint64_t resorted [mSize];
    RST(shareOfMatrix, mmaxParams[1], mmaxParams[0], mmaxParams[3], mmaxParams[3], resorted);
    double *d_matrix = convert2double(REC(proxy, resorted, mSize), mSize);
    double computed_max[number_of_windows];

    for(uint32_t win = 0; win < number_of_windows; win++){
        for(uint32_t win_element = 0; win_element < window_length; win_element++){
            double matrixVal = d_matrix[window_length*win + win_element];
            cout << matrixVal << " > " << computed_max[win] << " ?" << endl;
            if (matrixVal > computed_max[win]){
                computed_max[win] = matrixVal;
                cout << "Y" << endl;
            }
        }
        if (computed_max[win] != convert2double(reconstructed_max[win])){
            flag = false;
            break;
        }
    }

    if(flag){
        cout<<"Vectorized MAX works correctly"<<endl;
    }
    else{
        cout<<"Vectorized MAX works incorrectly"<<endl;
        print1DMatrixByWindows("Matrix: ", d_matrix, mmaxParams[0], mmaxParams[1], mmaxParams[2], mmaxParams[3]);
        print1DMatrixByWindows("computed max values (test): ", computed_max, 1, number_of_windows, 1, 1);
        print1DMatrixByWindows("VS result from method: ", convert2double(reconstructed_max, number_of_windows), 1, number_of_windows, 1, 1);
    }
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


    //double pp_result = convert2double(reconstructed_drelu, 0); dont need to convert to double relu is integer value delete this line
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

    int n_samples = 10;
    uint64_t *x = proxy->createShare(random_1D_data(proxy, n_samples, proxy->getMinPower(), proxy->getMaxPower()), n_samples);
    proxy->SendBytes(CORE_MEXP, n_samples);
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
    uint32_t size = 20; // size of the vector
    double min_val = -10;
    double max_val = 10;

    // generate a vector of random values
    uint64_t* vec1 = proxy->createShare(random_1D_data(proxy, size, min_val, max_val), size);
    uint64_t* vec2 = proxy->createShare(random_1D_data(proxy, size, min_val, max_val), size);

    // call DP
    proxy->SendBytes(CORE_DP, size);
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
    uint32_t size = 20; // size of the total vector
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
    proxy->SendBytes(CORE_MDP, size);
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
    double min_val = -10;
    double max_val = 10;

    uint64_t** mat1 = proxy->createShare(random_2D_data(proxy, a_row, a_col, min_val, max_val), a_row, a_col);
    uint64_t** mat2 = proxy->createShare( random_2D_data(proxy, a_col, b_col, min_val, max_val), a_col, b_col);

    proxy->SendBytes(CORE_MATMATMUL, a_row * a_col * b_col);
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
    double min_val = -10;
    double max_val = 10;

    uint64_t*** mat1 = new uint64_t**[n_matrices];
    uint64_t*** mat2 = new uint64_t**[n_matrices];
    for(int i = 0; i < n_matrices; i++) {
        mat1[i] = proxy->createShare(random_2D_data(proxy, a_row, a_col, min_val, max_val), a_row, a_col);
        mat2[i] = proxy->createShare(random_2D_data(proxy, a_col, b_col, min_val, max_val), a_col, b_col);
    }

    proxy->SendBytes(CORE_MMATMATMUL, n_matrices * a_row * a_col * b_col);
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
    double min_val = -10;
    double max_val = 10;

    uint64_t** mat = proxy->createShare(random_2D_data(proxy, a_row, a_col, min_val, max_val), a_row, a_col);
    uint64_t* vec = proxy->createShare(random_1D_data(proxy, a_col, min_val, max_val), a_col);

    proxy->SendBytes(CORE_MATVECMUL, a_row * a_col);
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
    double min_val = -10;
    double max_val = 10;

    uint64_t*** mat = new uint64_t**[n_matrices];
    uint64_t** vec = new uint64_t*[n_matrices];
    for(int i = 0; i < n_matrices; i++) {
        mat[i] = proxy->createShare(random_2D_data(proxy, a_row, a_col, min_val, max_val), a_row, a_col);
        vec[i] = proxy->createShare(random_1D_data(proxy, a_col, min_val, max_val), a_col);
    }

    proxy->SendBytes(CORE_MMATVECMUL, n_matrices * a_row * a_col);
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

void INVSQRT_Test(Party* proxy) {
    /* In this function, we test the computation of the inverse square root of a Gram matrix.
     * We first generate a random Gram matrix by first generating a random data matrix D and
     * then computing D^T * D.
     */

/**    cout<<setfill ('*')<<setw(50)<<"Calling INVSQRT";
    cout<<setfill ('*')<<setw(49)<<"*"<<endl;
    // setting
    int n_row = 4;
    int n_col = 5;

    // generate a Gram matrix
    uint64_t **invsqrt_data = random_gram_matrix(proxy, n_row, n_col);

    double** tmp = convert2double(REC(proxy, invsqrt_data, n_row, n_row), n_row, n_row);

    proxy->SendBytes(RKN_INVSQRT, n_row);
    uint64_t** invsqrt_G = INVSQRT(proxy, invsqrt_data, n_row);

    double** rec_invsqrt_G = convert2double(REC(proxy, invsqrt_G, n_row, n_row), n_row, n_row);

    print2DArray("The inverse square root of the Gram matrix", rec_invsqrt_G, n_row, n_row);

    double* straighten_invsqrt_G = new double[n_row * n_row];
    for(uint32_t i = 0; i < n_row * n_row; i++) {
        straighten_invsqrt_G[i] = tmp[i % n_row][i / n_row];
    }
    EigenSolver<Matrix<double, Dynamic, Dynamic>> ges;
    Map<Matrix<double, Dynamic, Dynamic>> matrix_G(straighten_invsqrt_G, n_row, n_row);
    ges.compute(matrix_G);
    Matrix<double, Dynamic, Dynamic> eig_vecs = ges.eigenvectors().real();
    Matrix<double, Dynamic, 1> eig_vals = ges.eigenvalues().real();
    cout << eig_vals << endl;

    cout << "============= GT Inverse square root of the Gram matrix ======================" << endl;
    Matrix<double, Dynamic, Dynamic> vals = eig_vals;
    cout << eig_vecs * vals.cwiseSqrt().cwiseInverse().asDiagonal() * Transpose(eig_vecs) << endl;
    cout << "============================================================================" << endl;

    print2DArray("The inverse of the Gram matrix",
                        multiply_matrices(rec_invsqrt_G, rec_invsqrt_G, n_row, n_row, n_row), n_row, n_row);

    cout << "============= GT Inverse of the Gram matrix ======================" << endl;
    cout << matrix_G.inverse() << endl;
    cout << "============================================================================" << endl;
*/
}

void MINVSQRT_Test(Party* proxy) {
    /* In this function, we test the computation of the inverse square root of a Gram matrix.
     * We first generate a random Gram matrix by first generating a random data matrix D and
     * then computing D^T * D.
     */
/*    cout<<setfill ('*')<<setw(50)<<"Calling MINVSQRT";
    cout<<setfill ('*')<<setw(49)<<"*"<<endl;
    // setting
    int n_row = 4;
    int n_col = 5;
    int n_gms = 3;

    // generate a Gram matrix
    uint64_t ***invsqrt_data = new uint64_t**[n_gms];
    double ***Gs = new double**[n_gms];
    for(int i = 0; i < n_gms; i++) {
        invsqrt_data[i] = random_gram_matrix(proxy, n_row, n_col);
        Gs[i] = convert2double(REC(proxy, invsqrt_data[i], n_row, n_row), n_row, n_row);
    }

//    double** tmp = convert2double(REC(proxy, invsqrt_data, n_row, n_row), n_row, n_row);

    proxy->SendBytes(RKN_MINVSQRT, n_gms,n_row);
    uint64_t*** invsqrt_G = INVSQRT(proxy, invsqrt_data, n_gms, n_row);

    double*** rec_invsqrt_G = new double**[n_gms];
    for(int g = 0; g < n_gms; g++) {
        rec_invsqrt_G[g] = convert2double(REC(proxy, invsqrt_G[g], n_row, n_row), n_row, n_row);
        print2DArray("The inverse square root of the Gram matrix", rec_invsqrt_G[g], n_row, n_row);

        double* straighten_invsqrt_G = new double[n_row * n_row];
        for(uint32_t i = 0; i < n_row * n_row; i++) {
            straighten_invsqrt_G[i] = Gs[g][i % n_row][i / n_row];
        }
        EigenSolver<Matrix<double, Dynamic, Dynamic>> ges;
        Map<Matrix<double, Dynamic, Dynamic>> matrix_G(straighten_invsqrt_G, n_row, n_row);
        ges.compute(matrix_G);
        Matrix<double, Dynamic, Dynamic> eig_vecs = ges.eigenvectors().real();
        Matrix<double, Dynamic, 1> eig_vals = ges.eigenvalues().real();
//        cout << eig_vals << endl;

        cout << "============= GT Inverse square root of the Gram matrix ======================" << endl;
        Matrix<double, Dynamic, Dynamic> vals = eig_vals;
        cout << eig_vecs * vals.cwiseSqrt().cwiseInverse().asDiagonal() * Transpose(eig_vecs) << endl;
        cout << "============================================================================" << endl;

        print2DArray("The inverse of the Gram matrix",
                     multiply_matrices(rec_invsqrt_G[g], rec_invsqrt_G[g], n_row, n_row, n_row), n_row, n_row);

        cout << "============= GT Inverse of the Gram matrix ======================" << endl;
        cout << matrix_G.inverse() << endl;
        cout << "============================================================================" << endl;
    }
*/
}

void DIV_Test(Party *proxy){
    cout<<setfill ('*')<<setw(50)<<"Calling DIV";
    cout<<setfill ('*')<<setw(49)<<"*"<<endl;

    uint64_t x = proxy->createShare(16.0);
    uint64_t y = proxy->createShare(4.0);

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
        cout<< "computed: " << reconstructed_div << " should be: " << computed_div << endl;
    }

}

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

    MUL_Test(proxy);
    MMUL_Test(proxy);

    MOC_Test(proxy);
    MMOC_Test(proxy);

    MSB_Test(proxy);
    MMSB_Test(proxy);

    CMP_Test(proxy);
    MCMP_Test(proxy);

    MUX_Test(proxy);
    MMUX_Test(proxy);

    MAX_Test(proxy);
//    MMAX_Test(proxy); //TODO adapt to asymmetric window size

    //RST_Test(proxy); // works (needs much space in console as it prints matrices)

//    RELU_Test(proxy);
//    DRLU_Test(proxy);

    DIV_Test(proxy);

//    EXP_Test(proxy);
//    MEXP_Test(proxy);

//    DP_Test(proxy);
//    MDP_Test(proxy);

//    MATMATMUL_Test(proxy);
//    MMATMATMUL_Test(proxy);

//    MATVECMUL_Test(proxy);
//    MMATVECMUL_Test(proxy);

//    INVSQRT_Test(proxy);
    MINVSQRT_Test(proxy);

    proxy->SendBytes(CORE_END);
    proxy->PrintBytes();
    return 0;
}
