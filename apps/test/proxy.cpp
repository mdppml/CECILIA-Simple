#include <cstdlib>
#include <iostream>
#include <deque>
#include <chrono>
#include <iomanip>
#include "../../core/core.h"
using namespace std;

constexpr int MIN = -100;
constexpr int MAX = 100;
constexpr int sz = 1000;
void MUL_Test(Party *proxy){
    cout<<setfill ('*')<<setw(50)<<"Calling MUL";
    cout<<setfill ('*')<<setw(49)<<"*"<<endl;

    double xd = MIN + (double)(proxy->generateCommonRandom()&RAND_MAX) / ((double)(RAND_MAX/(MAX - MIN)));
    double yd = MIN + (double)(proxy->generateCommonRandom()&RAND_MAX) / ((double)(RAND_MAX/(MAX - MIN)));
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
        double xd = MIN + (double)(proxy->generateCommonRandom()&RAND_MAX) / ((double)(RAND_MAX/(MAX - MIN)));
        double yd = MIN + (double)(proxy->generateCommonRandom()&RAND_MAX) / ((double)(RAND_MAX/(MAX - MIN)));
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
    double xd = MIN + (double)(proxy->generateCommonRandom()&RAND_MAX) / ((double)(RAND_MAX/(MAX - MIN)));
    double yd = MIN + (double)(proxy->generateCommonRandom()&RAND_MAX) / ((double)(RAND_MAX/(MAX - MIN)));
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
        double xd = MIN + (double)(proxy->generateCommonRandom()&RAND_MAX) / ((double)(RAND_MAX/(MAX - MIN)));
        double yd = MIN + (double)(proxy->generateCommonRandom()&RAND_MAX) / ((double)(RAND_MAX/(MAX - MIN)));
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
    double xd = MIN + (double)(proxy->generateCommonRandom()&RAND_MAX) / ((double)(RAND_MAX/(MAX - MIN)));
    double yd = MIN + (double)(proxy->generateCommonRandom()&RAND_MAX) / ((double)(RAND_MAX/(MAX - MIN)));
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
        double xd = MIN + (double)(proxy->generateCommonRandom()&RAND_MAX) / ((double)(RAND_MAX/(MAX - MIN)));
        double yd = MIN + (double)(proxy->generateCommonRandom()&RAND_MAX) / ((double)(RAND_MAX/(MAX - MIN)));
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

int main(int argc, char* argv[]) {
    uint8_t role = atoi(argv[1]);

    Party *proxy;
    if (role==0)
        proxy = new Party(P1,7777, "127.0.0.1", 8888, "127.0.0.1");
    else
        proxy = new Party(P2,7777, "127.0.0.1", 8888, "127.0.0.1");

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


    proxy->SendBytes(CORE_END);
    proxy->PrintBytes();
    return 0;
}
