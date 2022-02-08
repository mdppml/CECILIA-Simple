#include <cstdlib>
#include <iostream>
#include <deque>
#include <chrono>
#include <tuple>
#include <iomanip>
#include "../../core/core.h"

#include "../../core/cnn.h"
#include "../../utils/test_functions.h"
using namespace std;

constexpr int MIN_VAL = -100;
constexpr int MAX_VAL = 100;
constexpr int sz = 1000;
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
    else
        cout<<"MUL works incorrectly"<<endl;
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
    uint64_t r_computed = (x_reconstructed>>(L-1))<<FRAC;
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
        uint64_t r_computed = (x_reconstructed[i]>>(L-1))<<FRAC;
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
void MAX_Test(Party *proxy, double* random_matrix, uint32_t mCols, uint32_t mRows){
    cout<<setfill ('*')<<setw(50)<<"Calling MAX";
    cout<<setfill ('*')<<setw(49)<<"*"<<endl;

    uint64_t mSize = mCols * mRows;
    uint64_t *shareOfMatrix = new uint64_t [mSize];
    shareOfMatrix = proxy->createShare(random_matrix, mSize);
    // print1DMatrixByWindows("Ground truth matrix for MAX", shareOfMatrix, mRows, mCols, mRows, mCols);

    proxy->SendBytes(CNN_MAX, mSize);

    uint64_t max = MAX(proxy, shareOfMatrix, mSize);
    uint64_t reconstructed_max_value = REC(proxy, max);
    //cout << "RECONSTRUCTED MAX = " << reconstructed_max_value << endl;

    // checking the result
    double computed_max_value = MIN_VAL;
    for(uint32_t i = 0; i<mSize; i++){
        double matrixVal = random_matrix[i];
        //cout << to_string(matrixVal) << " > " << to_string(refMax) << " ? --> ";
        if (matrixVal > computed_max_value){
            computed_max_value = matrixVal;
        }
    }
    // TODO check for precision > 0
    u_int64_t comp_result = convert2uint64(computed_max_value, 0);
    if(comp_result == reconstructed_max_value){
        cout<<"MAX works correctly"<<endl;
    }
    else{
        cout<<"MAX works incorrectly" <<endl;
        cout<< "computed: " << reconstructed_max_value << " should be: " << comp_result << endl;
    }

}
void MMAX_Test(Party *proxy, double* random_matrix, uint32_t mCols, uint32_t mRows, uint32_t wCols, uint32_t wRows){
    cout<<setfill ('*')<<setw(50)<<"Calling vectorized MAX";
    cout<<setfill ('*')<<setw(49)<<"*"<<endl;

    uint32_t mSize = mCols * mRows;
    uint64_t *shareOfMatrix = new uint64_t [mSize];
    shareOfMatrix = proxy->createShare(random_matrix, mSize);
    //print1DMatrixByWindows("Ground truth matrix for vectorized MAX", shareOfMatrix, mRows, mCols, wRows, wCols);

    uint32_t mmaxParams;
    //contains matrix_size bitwise at 16 MSBs and window size at the 16 LSBs
    mmaxParams = ((uint64_t)mSize << 16) + (uint64_t) wCols;
    cout << "It should be: mSize = " << mSize << " and wCols = " << wCols << endl;
    proxy->SendBytes(CNN_MMAX, mmaxParams);
    uint64_t *max = MAX(proxy, shareOfMatrix, mSize, wCols);

    uint64_t number_of_windows = floor(mSize/(wCols*wRows));
    uint64_t* reconstructed_max = REC(proxy, max, number_of_windows);

    // checking the result
    bool flag = true;
    uint64_t *resorted = new uint64_t [mSize];
    resorted = RST(shareOfMatrix, mCols, mRows, wCols, wRows, resorted);

    // TODO check for precision > 0
    double *d_matrix = convert2double(REC(proxy, resorted, mSize), mSize, 0);
    double *computed_max = new double[number_of_windows];

    for(uint32_t win = 0; win < mRows; win++){
        for(uint32_t win_element = 0; win_element < mCols; win_element++){
            double matrixVal = d_matrix[win_element];
            if (matrixVal > computed_max[win]){
                computed_max[win] = matrixVal;
            }
        }
        // TODO check for precision > 0
        if (computed_max[win] != convert2double(reconstructed_max[win], 0)){
            flag = false;
            break;
        }
    }

    if(flag){
        cout<<"Vectorized MAX works correctly"<<endl;
    }
    else{
        cout<<"Vectorized MAX works incorrectly"<<endl;
    }
}

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

    // generate random matrix for MAX tests:
    double * random_matrix;
    uint32_t * matrix_dimensions = new uint32_t [4]; // this will always contain 4 values: mColSize, mRowSize, wColSize, wRowSize
    cout << "randMatrix: " << random_matrix << " dims: " << matrix_dimensions << "; matrix_dims: " << matrix_dimensions << endl;
    std::tuple(random_matrix, matrix_dimensions) = random_window_matrix(proxy, 2);
    cout << "created tuple. Accessing values..." << "; matrix_dims: " << matrix_dimensions << endl;
    uint32_t mColSize = matrix_dimensions[0];
    uint32_t mRowSize = matrix_dimensions[1];
    uint32_t wColSize = matrix_dimensions[2];
    uint32_t wRowSize = matrix_dimensions[3];
    cout << "dimensions are: " << mColSize << " x " << mRowSize << " with windows " << wColSize << " x " << wColSize << " " << endl;

    MAX_Test(proxy, random_matrix, mColSize, mRowSize);
    //MMAX_Test(proxy, random_matrix, mColSize, mRowSize, wColSize, wColSize); //TODO adapt to asymmetric window size


    proxy->SendBytes(CORE_END);
    proxy->PrintBytes();
    return 0;
}
