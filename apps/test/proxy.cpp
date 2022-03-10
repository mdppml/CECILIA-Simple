#include <cstdlib>
#include <iostream>
#include <deque>
#include <chrono>
#include <tuple>
#include <iomanip>
#include "../../core/core.h"

#include "../../core/cnn.h"
#include "../../utils/test_functions.h"
#include "../../utils/flib.h"
using namespace std;

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
void MAX_Test(Party *proxy){
    cout<<setfill ('*')<<setw(50)<<"Calling MAX";
    cout<<setfill ('*')<<setw(49)<<"*"<<endl;

    uint32_t mRows = WSZ*10;
    uint32_t mCols = WSZ*10;
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

    double pp_result = convert2double(REC(proxy, max), 0);
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

    uint32_t mRows = WSZ*10;
    uint32_t mCols = WSZ*10;
    uint64_t mSize = mCols*mRows;

    uint32_t wRows = WSZ;
    uint32_t wCols = WSZ;
    uint64_t *shareOfMatrix = proxy->createShare(random_1D_data(proxy, mSize, 32), mSize);
    auto *resorted = new uint64_t [sz];
    print1DMatrixByWindows("RST original", convert2double(REC(proxy, shareOfMatrix, mSize), mSize), mCols, mRows, wRows, wCols);
    RST(shareOfMatrix, mCols, mRows, wCols, wRows, resorted);

    print1DMatrixByWindows("RST finished: resorted", convert2double(REC(proxy, resorted, mSize), mSize), mCols, mRows, wRows, wCols);
}

void MMAX_Test(Party *proxy){
    cout<<setfill ('*')<<setw(50)<<"Calling vectorized MAX";
    cout<<setfill ('*')<<setw(49)<<"*"<<endl;

    auto *mmaxParams = new uint64_t [4];
    uint32_t mRows = WSZ*10;
    mmaxParams[0] = mRows;
    uint32_t mCols = WSZ*10;
    mmaxParams[1] = mCols;
    uint64_t mSize = mCols*mRows;

    uint32_t wRows = WSZ;
    mmaxParams[2] = wRows;
    uint32_t wCols = WSZ;
    mmaxParams[3] = wCols;

    uint64_t *shareOfMatrix = proxy->createShare(random_1D_data(proxy, mSize), mSize);

    uint8_t * buffer = proxy->getBuffer1();
    int socket_helper = proxy->getSocketHelper();

    proxy->SendBytes(CNN_MMAX);
    unsigned char *ptr_out = &buffer[0];
    addVal2CharArray(mmaxParams, &ptr_out, 4);
    Send(socket_helper, buffer, 4 * 8);

    //contains mRows bitwise at 16 MSBs, then mCols at next 16 MSBs and window size at the 16 LSBs (16 are still reserved for different win size)
    //mmaxParams = ((uint64_t)mRows << 48) + ((uint64_t)mCols << 32) + (uint64_t) wCols;
    //TODO send mmaxParams in buffer!!
    uint64_t *max = MAX(proxy, shareOfMatrix, mRows, mCols, wCols);

    uint64_t number_of_windows = floor(mSize/(wCols*wRows));
    uint64_t* reconstructed_max = REC(proxy, max, number_of_windows);

    // checking the result
    bool flag = true;
    //uint64_t *resorted = new uint64_t [mSize];
    //resorted = RST(shareOfMatrix, mCols, mRows, wCols, wRows, resorted);

    // TODO check for precision > 0
    double *d_matrix = convert2double(REC(proxy, shareOfMatrix, mSize), mSize, 0);
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
    if (originalX > 0){
        computed_relu = originalX;
    }
    else{
        //cout << "X = " << to_string(computed_relu) << " --> RELU = 0." << endl;
        computed_relu = 0;
    }

    double pp_result = convert2double(reconstructed_relu, 0);
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


    //double pp_result = convert2double(reconstructed_drelu, 0); dont need to convert to double relu is integer value delete this line
    if(computed_drelu == reconstructed_drelu){
        cout<<"DRLU works correctly"<<endl;
    }
    else{
        cout<<"DRLU works incorrectly" <<endl;
        cout<< "computed: " << reconstructed_drelu << " should be: " << computed_drelu << endl;
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

    MAX_Test(proxy);

    //RST_Test(proxy);

    //MMAX_Test(proxy); //TODO adapt to asymmetric window size

    RELU_Test(proxy);
    DRLU_Test(proxy);

    proxy->SendBytes(CORE_END);
    proxy->PrintBytes();
    return 0;
}
