#include <cstdlib>
#include <iostream>
#include <deque>
#include <chrono>
#include <tuple>
#include <iomanip>
#include <ctype.h>
#include "../../utils/parse_options.h"
#include "../../core/core.h"
#include "../../core/cnn.h"
#include "../../core/rkn.h"
#include "../../utils/flib.h"
#include "../cnn/mnist_reader.hpp"
#include <bitset>
#include <unordered_map>


using namespace std;
//using namespace Eigen;

constexpr int MIN_VAL = -100;
constexpr int MAX_VAL = static_cast<int>((uint64_t) 1 << 43);
constexpr int sz = 10;
constexpr int WSZ = 3;
constexpr int num_repetition = 1;

// ************************************ Ali  ***********************************************
bool MUL_Test_v2(Party *proxy, int i, unordered_map<string, int> &cases, int &cnt){
//    cout << setfill('*') << setw(50) << "Calling MultiplyNarrow";
//    cout << setfill('*') << setw(49) << "*" << endl;

    double xd = MIN_VAL + (double)(proxy->GenerateCommonRandom() & RAND_MAX) / ((double)(RAND_MAX / (MAX_VAL - MIN_VAL)));
    double yd = MIN_VAL + (double)(proxy->GenerateCommonRandom() & RAND_MAX) / ((double)(RAND_MAX / (MAX_VAL - MIN_VAL)));
    uint64_t x = proxy->CreateShare(xd);
    uint64_t y = proxy->CreateShare(yd);
    proxy->SendBytes(coreMultiply);
    uint64_t r = Multiply(proxy, x, y);
    // checking the result
    uint64_t rec_x = Reconstruct(proxy, x);
    uint64_t rec_y = Reconstruct(proxy, y);
    uint64_t rec_r = Reconstruct(proxy, r);
    xd = ConvertToDouble(rec_x);
    yd = ConvertToDouble(rec_y);
    double rd = ConvertToDouble(rec_r);
    double rcd = (xd*yd);

    string key = to_string((rec_x > x)) + to_string(rec_y > y) + to_string(rec_r > r);
    if(cases.find(key) == cases.end()) {
        cases[key] = 1;
    }
    else {
        cases[key]++;
    }
//    if(key == "001" || key == "010" || key == "100" || key == "111") {
//        cout << "-----------------------------------------" << endl;
//        cout << "x: " << xd << "\ny: " << yd << "\nComputed r: " << rd << "\nGT r: " << rcd << endl;
//        cout << "Bitwise computed r: " << bitset<64>(rec_r) << endl;
//        cout << "-----------------------------------------" << endl;
//    }
    if ((int)(rd - rcd) == 0) {
//        cout<<"MultiplyNarrow works correctly"<<endl;
//        cout << "x: " << xd << "\ny: " << yd << "\nComputed r: " << rd << "\nGT r: " << rcd << endl;
//        cout << "Bitwise computed r: " << bitset<64>(rec_r) << endl;
//        cout << "-----------------------------------------" << endl;
        return true;
    }
    else {
//        cout << "----------------" << i << "-------------------------" << endl;
//        cout << "x: " << rec_x << endl;
//        cout << "Share of x: " << x << endl;
//        cout << "rec_x >? share_x: " << (rec_x > x) << endl;
//        cout << "y: " << rec_y << endl;
//        cout << "Share of y: " << y << endl;
//        cout << "rec_y >? share_y: " << (rec_y > y) << endl;
//        cout << "r: " << rec_r << endl;
//        cout << "Share of r: " << r << endl;
//        cout << "rec_r >? share_r: " << (rec_r > r) << endl;
//        cout << "-----------------------------------------" << endl;
//        cout<<"MultiplyNarrow works incorrectly"<<endl;
//        cout << "x: " << xd <""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""< "\ny: " << yd << "\nComputed r: " << rd << "\nGT r: " << rcd << endl;
//        cout << "Bitwise computed r: " << bitset<64>(rec_r) << endl;
//        cnt++;
//        cout << "-----------------------------------------" << endl;
    }
    return ((int) (rd - rcd) == 0);
}

bool MMUL_Test_v2(Party *proxy){
    cout << setfill('*') << setw(50) << "Calling MMUL";
    cout << setfill('*') << setw(49) << "*" << endl;
    uint64_t x[sz], y[sz], z[sz];
    uint32_t *params = new uint32_t[1];
    params[0] = sz;
    for (int i = 0; i < sz; i++) {
        double xd = MIN_VAL +
                    (double) (proxy->GenerateCommonRandom() & RAND_MAX) / ((double) (RAND_MAX / (MAX_VAL - MIN_VAL)));
        double yd = MIN_VAL +
                    (double) (proxy->GenerateCommonRandom() & RAND_MAX) / ((double) (RAND_MAX / (MAX_VAL - MIN_VAL)));
        x[i] = proxy->CreateShare(xd);
        y[i] = proxy->CreateShare(yd);
    }
    proxy->SendBytes(coreVectorisedMultiply, params, 1);
    uint64_t *r = Multiply(proxy, x, y, sz);
    // checking the result
    unordered_map<string, int> umap;
    bool flag = true;
    for (int i=0;i<sz;i++) {
        uint64_t rec_x = Reconstruct(proxy, x[i]);
        uint64_t rec_y = Reconstruct(proxy, y[i]);
        uint64_t rec_r = Reconstruct(proxy, r[i]);
        double xd = ConvertToDouble(rec_x);
        double yd = ConvertToDouble(rec_y);
        double rd = ConvertToDouble(rec_r);
        double rcd = (xd * yd);
        cout << "----------------" << i << "-------------------------" << endl;
        cout << "x: " << rec_x << endl;
        cout << "Share of x: " << x[i] << endl;
        cout << "rec_x >? share_x: " << (rec_x > x[i]) << endl;
        cout << "y: " << rec_y << endl;
        cout << "Share of y: " << y[i] << endl;
        cout << "rec_y >? share_y: " << (rec_y > y[i]) << endl;
        cout << "r: " << rec_r << endl;
        cout << "Share of r: " << r[i] << endl;
        cout << "rec_r >? share_r: " << (rec_r > r[i]) << endl;
        cout << "-----------------------------------------" << endl;
        string key = to_string((rec_x > x[i])) + to_string(rec_y > y[i]) + to_string(rec_r > r[i]);
        if(umap.find(key) == umap.end()) {
            umap[key] = 1;
        }
        else {
            umap[key]++;
        }
        if ((int) (rd - rcd) != 0) {
            flag = false;
            cout << "Absolute difference of multiplication of " << xd << " and " << yd << ": " << abs(rd - rcd) << endl;
            for (auto& it: umap) {
                // Do stuff
                cout << it.first << ": " << it.second << endl;
            }
            break;
        }
    }
    if (flag) {
        cout<<"MMUL works correctly"<<endl;
        return true;
    }
    else {
        cout<<"MMUL works incorrectly"<<endl;
        return false;
    }
}
// ****************************************************************************************

bool TRUNCATE_Test(Party *proxy, int &cnt, unordered_map<string, int> &map_cnt) {
//    cout<<setfill ('*')<<setw(50)<<"Calling Truncate";
//    cout<<setfill ('*')<<setw(49)<<"*"<<endl;

    // negative values
//    uint64_t x = 0xf000000000000000; // small negative value
//    uint64_t x = 0x8000000000000000; // large negative value
//    uint64_t x = (uint64_t) 1 | ((uint64_t) 1 << (L_BIT - 1));

    // large positive value
//    uint64_t x = 0x7fffffffffffffff;
//    uint64_t x = ((uint64_t) 0 - 1) - 0xffff000000000001;
//    uint64_t x = ((uint64_t) 0 - 1) - 0xffffff0000000000;
    uint64_t x = proxy->GenerateCommonRandom();

    uint64_t gt_trun_x;
    if((x >> (L_BIT - 1)) == 0x0) {
        gt_trun_x = x >> FRACTIONAL_BITS;
    }
    else {
        gt_trun_x = ((((uint64_t) 1 << FRACTIONAL_BITS) - 1) << (L_BIT - FRACTIONAL_BITS)) | (x >> FRACTIONAL_BITS);
    }

//    cout << "x: " << bitset<64>(x) << " => " << ConvertToDouble(x) << endl;

    uint64_t xi = proxy->CreateShare(x);
//    uint64_t rec_x = Reconstruct(proxy, xi);
//    cout << "x: " << bitset<64>(rec_x) << endl;
//    cout << "xi: " << bitset<64>(xi) << " => " << ConvertToDouble(xi) << endl;

/*
    if(proxy->getPRole() == P1) {
        uint64_t tmp = xi >> FRACTIONAL_BITS;
        cout << "xi >> FRACTIONAL_BITS: " << bitset<64>(tmp) << endl;
    }
    else if(proxy->getPRole() == P2) {
        uint64_t tmp = -1 * xi;
        cout << "-1 * xi: " << bitset<64>(tmp) << endl;
        tmp = tmp >> FRACTIONAL_BITS;
        cout << "(-1 * xi) >> FRACTIONAL_BITS: " << bitset<64>(tmp) << endl;
        tmp = -1 * tmp;
        cout << "-1 * ((-1 * xi) >> FRACTIONAL_BITS): " << bitset<64>(tmp) << endl;
    }
*/
    uint64_t trun_x = Truncate(proxy, xi);
    uint64_t rec_trun_x = Reconstruct(proxy, trun_x);

//    cout << "GT: " << bitset<64>(gt_trun_x) << " => " << ConvertToDouble(gt_trun_x) << endl;
//    cout << "CM: " << bitset<64>(rec_trun_x) << " => " << ConvertToDouble(rec_trun_x) << endl;

    uint64_t r = xi - x;
    uint64_t neg_x = (uint64_t) 0 - x;

    cout << "------------------------------------------" << endl;
    cout << " x: " << bitset<64>(x) << " => " << ConvertToDouble(x) << endl;
//    cout << "-x: " << bitset<64>(neg_x) << " => " << ConvertToDouble(neg_x) << endl;
//    cout << " r: " << bitset<64>(r) << " => " << ConvertToDouble(r) << endl;
//    cout << "xi: " << bitset<64>(xi) << " => " << ConvertToDouble(xi) << endl;
//    if(proxy->getPRole() == P1) {
//        uint64_t tmp = xi >> FRACTIONAL_BITS;
//        cout << "xi >> FRACTIONAL_BITS: " << bitset<64>(tmp) << endl;
//    }
//    else if(proxy->getPRole() == P2) {
//        uint64_t tmp = -1 * xi;
//        cout << "-1 * xi: " << bitset<64>(tmp) << endl;
//        tmp = tmp >> FRACTIONAL_BITS;
//        cout << "(-1 * xi) >> FRACTIONAL_BITS: " << bitset<64>(tmp) << endl;
//        tmp = -1 * tmp;
//        cout << "-1 * ((-1 * xi) >> FRACTIONAL_BITS): " << bitset<64>(tmp) << endl;
//    }
//    cout << "GT: " << bitset<64>(gt_trun_x) << " => " << ConvertToDouble(gt_trun_x) << endl;
//    cout << "CM: " << bitset<64>(rec_trun_x) << " => " << ConvertToDouble(rec_trun_x) << endl;

    cout << (neg_x > r) << " -- " << (r > x) << endl;
    string key = to_string((x >> (L_BIT - 1)) & 0x1) + to_string(neg_x > r) + "" + to_string(r > x);
    if(map_cnt.find(key) == map_cnt.end()) {
        map_cnt[key] = 1;
    }
    else {
        map_cnt[key]++;
    }
    if(rec_trun_x - gt_trun_x <= 0x1) {
        cout << "Truncate works correctly" << endl;
        return true;
    }
    else {
        cnt++;
        cout << "Truncate works incorrectly" << endl;
        return false;
    }
}

bool MUL_Test(Party *proxy){
    cout<<setfill ('*')<<setw(50)<<"Calling Multiply";
    cout<<setfill ('*')<<setw(49)<<"*"<<endl;

    double xd = MIN_VAL + (double)(proxy->GenerateCommonRandom() & RAND_MAX) / ((double)(RAND_MAX / (MAX_VAL - MIN_VAL)));
    double yd = MIN_VAL + (double)(proxy->GenerateCommonRandom() & RAND_MAX) / ((double)(RAND_MAX / (MAX_VAL - MIN_VAL)));
    uint64_t x = proxy->CreateShare(xd);
    uint64_t y = proxy->CreateShare(yd);
    proxy->SendBytes(coreMultiply);
    uint64_t r = Multiply(proxy, x, y);
    // checking the result
    uint64_t rec_x = Reconstruct(proxy, x);
    uint64_t rec_y = Reconstruct(proxy, y);
    uint64_t rec_r = Reconstruct(proxy, r);
    xd = ConvertToDouble(rec_x);
    yd = ConvertToDouble(rec_y);
    double rd = ConvertToDouble(rec_r);
    double rcd = (xd*yd);

    if ((int)(rd - rcd) == 0) {
        cout<<"Multiply works correctly"<<endl;
        cout << "x: " << xd << "\ny: " << yd << "\nComputed r: " << rd << "\nGT r: " << rcd << endl;
        return true;
    }
    else {
        cout<<"Multiply works incorrectly"<<endl;
        cout << "x: " << xd << "\ny: " << yd << "\nComputed r: " << rd << "\nGT r: " << rcd << endl;
        return false;
    }
    return ((int) (rd - rcd) == 0);
}

bool MMUL_Test(Party *proxy, double &exe_time, bool only_timing = false){
    cout<<setfill ('*')<<setw(50)<<"Calling MMUL";
    cout<<setfill ('*')<<setw(49)<<"*"<<endl;
    uint64_t *x = new uint64_t[sz];
    uint64_t *y = new uint64_t[sz];
    uint32_t* params = new uint32_t[1];
    params[0] = sz;
    for (int i=0;i<sz;i++){
        double xd = MIN_VAL + (double)(proxy->GenerateCommonRandom() & MAX_MULTIPLY) / ((double)(MAX_MULTIPLY / (MAX_VAL - MIN_VAL)));
        double yd = MIN_VAL + (double)(proxy->GenerateCommonRandom() & MAX_MULTIPLY) / ((double)(MAX_MULTIPLY / (MAX_VAL - MIN_VAL)));
        x[i] = proxy->CreateShare(xd);
        y[i] = proxy->CreateShare(yd);
    }
    proxy->SendBytes(coreMultiply);
    auto r = Multiply(proxy, x[0], y[0]);

    bool flag = true;
    unordered_map<string, int> umap;
    uint64_t rec_x = Reconstruct(proxy, x[0]);
    uint64_t rec_y = Reconstruct(proxy, y[0]);
    uint64_t rec_r = Reconstruct(proxy, r);
    double xd = ConvertToDouble(rec_x);
    double yd = ConvertToDouble(rec_y);
    double rd = ConvertToDouble(rec_r);
    double rcd = (xd * yd);
    if ((int) (rd - rcd) != 0) {
        flag = false;
        cout << "Absolute difference of multiplication of " << xd << " and " << yd << ": " << abs(rd - rcd) << endl;
        for (auto& it: umap) {
            // Do stuff
            cout << it.first << ": " << it.second << endl;
        }
    }
    delete [] x;
    delete [] y;
    delete [] params;

    if (flag) {
        cout<<"MMUL works correctly"<<endl;
        return true;
    }
    else {
        cout<<"MMUL works incorrectly"<<endl;
        return false;
    }
}

bool MMUL2_Test(Party *proxy){
    cout<<setfill ('*')<<setw(50)<<"Calling MMUL";
    cout<<setfill ('*')<<setw(49)<<"*"<<endl;
    uint64_t *x = new uint64_t[sz];
    uint64_t *y = new uint64_t[sz];
    uint32_t* params = new uint32_t[2];
    cout<<"SIZE:\t" << sz << endl;
    auto bsz = 23;
    params[0] = sz;
    params[1] = bsz;
    for (int i=0;i<sz;i++){
        x[i] = 0x7fffff-i;
        y[i] = i;
    }
    uint64_t *r;
    auto start = chrono::high_resolution_clock::now();
    for (int i = 0; i < 2; ++i) {
        proxy->SendBytes(coreVectorisedMultiply2, params, 2);
        r = Multiply(proxy, x, y, sz, bsz);
    }
    auto end = chrono::high_resolution_clock::now();
    double time_taken = chrono::duration_cast<chrono::nanoseconds>(end - start).count();
    time_taken *= 1e-9;
    cout<<"MultiplyNarrow Time:\t" << fixed
        << time_taken << setprecision(9) << " sec" << endl;
    // checking the result
    auto rec_r = ReconstructNarrow(proxy, r, sz, 23);
    for (int i = 0; i < sz; ++i) {
        cout << rec_r[i] << endl;
    }
    unordered_map<string, int> umap;
    bool flag = true;

    delete [] x;
    delete [] y;
    delete [] params;
    delete [] r;

    if (flag) {
        cout<<"MMUL works correctly"<<endl;
        return true;
    }
    else {
        cout<<"MMUL works incorrectly"<<endl;
        return false;
    }
}

void MOC_Test(Party *proxy) {
    cout << setfill('*') << setw(50) << "Calling ModularConversion";
    cout << setfill('*') << setw(49) << "*" << endl;
    uint64_t x = proxy->GenerateRandom() & N1_MASK;
    proxy->SendBytes(coreModularConversion);
    uint64_t r = ModularConversion(proxy, x);
    // checking the result
    uint64_t x_reconstructed = Reconstruct(proxy, x, N1_MASK);
    uint64_t r_reconstructed = Reconstruct(proxy, r);
    if (x_reconstructed == r_reconstructed)
        cout << "ModularConversion works correctly" << endl;
    else
        cout << "ModularConversion works incorrectly" << endl;
}

bool MMOC_Test(Party *proxy, double &exe_time, bool only_timing = false) {
    cout << setfill('*') << setw(50) << "Calling Vectorized ModularConversion";
    cout << setfill('*') << setw(49) << "*" << endl;
    uint64_t *x = new uint64_t[sz];
    uint32_t *params = new uint32_t[1];
    params[0] = sz;

    for (int i = 0; i < sz; i++)
        x[i] = proxy->GenerateRandom() & N1_MASK;

    proxy->SendBytes(coreVectorisedModularConversion, params, 1);
    auto start = chrono::high_resolution_clock::now();
    uint64_t *r = ModularConversion(proxy, x, sz);
    auto end = chrono::high_resolution_clock::now();
    double time_taken = chrono::duration_cast<chrono::nanoseconds>(end - start).count();
    time_taken *= 1e-9;
    cout << "Vectorized MOC Time:\t" << fixed << time_taken << setprecision(9) << " sec" << endl;
    exe_time = time_taken;

    // checking the result
    bool flag = true;
    if(!only_timing) {
        uint64_t *x_reconstructed = Reconstruct(proxy, x, sz, N1_MASK);
        uint64_t *r_reconstructed = Reconstruct(proxy, r, sz);
        for (int i = 0; i < sz; i++) {
            if (x_reconstructed[i] != r_reconstructed[i]) {
                flag = false;
                break;
            }
        }
    }
    if (flag)
        cout << "Vectorized ModularConversion works correctly" << endl;
    else
        cout << "Vectorized ModularConversion works incorrectly" << endl;
    return flag;
}

void MSB_Test(Party *proxy) {
    cout << setfill('*') << setw(50) << "Calling MostSignificantBit";
    cout << setfill('*') << setw(49) << "*" << endl;
    uint64_t x = 3 - 5; //proxy->GenerateRandom();
    proxy->SendBytes(coreMostSignificantBit);
    uint64_t r = MostSignificantBit(proxy, x);
    // checking the result
    uint64_t x_reconstructed = Reconstruct(proxy, x);
    uint64_t r_reconstructed = Reconstruct(proxy, r);
    uint64_t r_computed = (x_reconstructed >> (L_BIT - 1)) << FRACTIONAL_BITS;
    if (r_reconstructed == r_computed) {
        cout << "MostSignificantBit works correctly" << endl;
    } else {
        cout << "MostSignificantBit works incorrectly" << endl;
        cout << "computed MostSignificantBit = " << bitset<64>(r_computed) << "reconstructed MostSignificantBit = " << bitset<64>(r_reconstructed)
             << "for original X = " << bitset<64>(x_reconstructed) << endl;
    }
}

void MMSB_Test(Party *proxy, double &exe_time, bool only_timing = false) {
    cout << setfill('*') << setw(50) << "Calling Vectorized MostSignificantBit";
    cout << setfill('*') << setw(49) << "*" << endl;
    uint32_t *params = new uint32_t[1];
    params[0] = sz;

    uint64_t x[sz];
    for (int i = 0; i < sz; i++) {
        x[i] = proxy->GenerateRandom();
    }
    uint64_t *r;
    proxy->SendBytes(coreVectorisedMostSignificantBit, params, 1);
    auto start = chrono::high_resolution_clock::now();
    r = MostSignificantBit(proxy, x, sz);
    auto end = chrono::high_resolution_clock::now();
    double time_taken = chrono::duration_cast<chrono::nanoseconds>(end - start).count();
    time_taken *= 1e-9;
    cout << "Vectorized MostSignificantBit Time:\t" << fixed << time_taken << setprecision(9) << " sec" << endl;
    exe_time = time_taken;

    // checking the result
    bool flag = true;
    if(!only_timing) {
        uint64_t *x_reconstructed = Reconstruct(proxy, x, sz);
        uint64_t *r_reconstructed = Reconstruct(proxy, r, sz);
        for (int i = 0; i < sz; i++) {
            uint64_t r_computed = (x_reconstructed[i] >> (L_BIT - 1)) << FRACTIONAL_BITS;
            if (r_computed != r_reconstructed[i]) {
                flag = false;
                break;
            }
        }
    }
    if (flag) {
        cout << "Vectorized MostSignificantBit works correctly" << endl;
    } else {
        cout << "Vectorized MostSignificantBit works incorrectly" << endl;
    }
}

void CMP_Test(Party *proxy) {
    cout << setfill('*') << setw(50) << "Calling Compare";
    cout << setfill('*') << setw(49) << "*" << endl;
    double xd =
            MIN_VAL + (double) (proxy->GenerateCommonRandom() & RAND_MAX) / ((double) (RAND_MAX / (MAX_VAL - MIN_VAL)));
    double yd =
            MIN_VAL + (double) (proxy->GenerateCommonRandom() & RAND_MAX) / ((double) (RAND_MAX / (MAX_VAL - MIN_VAL)));
    uint64_t x = proxy->CreateShare(xd);
    uint64_t y = proxy->CreateShare(yd);
    proxy->SendBytes(coreCompare);
    uint64_t r = Compare(proxy, x, y);
    // checking the result
    xd = ConvertToDouble(Reconstruct(proxy, x));
    yd = ConvertToDouble(Reconstruct(proxy, y));
    uint64_t rd = Reconstruct(proxy, r);
    uint64_t r_computed = (xd >= yd) << FRACTIONAL_BITS;
    if (rd == r_computed)
        cout << "Compare works correctly" << endl;
    else
        cout << "Compare works incorrectly" << endl;
}

bool MCMP_Test(Party *proxy, double &exe_time, bool only_timing = false) {
    cout << setfill('*') << setw(50) << "Calling Vectorized Compare";
    cout << setfill('*') << setw(49) << "*" << endl;
    uint64_t *x = new uint64_t[sz];
    uint64_t *y = new uint64_t[sz];
    uint32_t *params = new uint32_t[1];
    params[0] = sz;

    for (int i = 0; i < sz; i++) {
        double xd = MIN_VAL +
                    (double) (proxy->GenerateCommonRandom() & RAND_MAX) / ((double) (RAND_MAX / (MAX_VAL - MIN_VAL)));
        double yd = MIN_VAL +
                    (double) (proxy->GenerateCommonRandom() & RAND_MAX) / ((double) (RAND_MAX / (MAX_VAL - MIN_VAL)));
        x[i] = proxy->CreateShare(xd);
        y[i] = proxy->CreateShare(yd);
    }
    proxy->SendBytes(coreVectorisedCompare, params, 1);
    auto start = chrono::high_resolution_clock::now();
    uint64_t *r = Compare(proxy, x, y, sz);
    auto end = chrono::high_resolution_clock::now();
    double time_taken = chrono::duration_cast<chrono::nanoseconds>(end - start).count();
    time_taken *= 1e-9;
    cout << "Vectorized CMP Time:\t" << fixed << time_taken << setprecision(9) << " sec" << endl;
    exe_time = time_taken;

    // checking the result
    bool flag = true;
    if(!only_timing) {
        uint64_t *x_reconstructed = Reconstruct(proxy, x, sz);
        uint64_t *y_reconstructed = Reconstruct(proxy, y, sz);
        uint64_t *r_reconstructed = Reconstruct(proxy, r, sz);

        for (int i = 0; i < sz; i++) {
            uint64_t r_computed = (ConvertToDouble(x_reconstructed[i]) >= ConvertToDouble(y_reconstructed[i])) << FRACTIONAL_BITS;
            if (r_computed != r_reconstructed[i]) {
                flag = false;
                break;
            }
        }
    }
    if (flag)
        cout << "Vectorized Compare works correctly" << endl;
    else
        cout << "Vectorized Compare works incorrectly" << endl;
    return flag;
}

void MUX_Test(Party *proxy, int &cnt) {
//    cout << setfill('*') << setw(50) << "Calling Multiplex";
//    cout << setfill('*') << setw(49) << "*" << endl;
    double xd =
            MIN_VAL + (double) (proxy->GenerateCommonRandom() & RAND_MAX) / ((double) (RAND_MAX / (MAX_VAL - MIN_VAL)));
    double yd =
            MIN_VAL + (double) (proxy->GenerateCommonRandom() & RAND_MAX) / ((double) (RAND_MAX / (MAX_VAL - MIN_VAL)));
    double zd = (double) (proxy->GenerateCommonRandom() & 0x1);
    uint64_t x = proxy->CreateShare(xd);
    uint64_t y = proxy->CreateShare(yd);
    uint64_t z = proxy->CreateShare(zd);

    proxy->SendBytes(coreMultiplex);
    uint64_t r = Multiplex(proxy, x, y, z);
    // checking the result
    uint64_t x_reconstructed = Reconstruct(proxy, x);
    uint64_t y_reconstructed = Reconstruct(proxy, y);
    uint64_t z_reconstructed = Reconstruct(proxy, z);
    uint64_t r_reconstructed = Reconstruct(proxy, r);
    uint64_t r_computed;
    if (z_reconstructed == 0)
        r_computed = x_reconstructed;
    else if (z_reconstructed == (1 << FRACTIONAL_BITS))
        r_computed = y_reconstructed;
    if (r_reconstructed == r_computed) {
//        cout << "Multiplex works correctly" << endl;
    }
    else {
        cnt++;
//        cout << "Multiplex works incorrectly" << endl;
    }
}

bool MMUX_Test(Party *proxy, double &exe_time, bool only_timing = false) {
    cout << setfill('*') << setw(50) << "Calling Vectorized Multiplex";
    cout << setfill('*') << setw(49) << "*" << endl;
    uint64_t *x = new uint64_t[sz];
    uint64_t *y = new uint64_t[sz];
    uint64_t *z = new uint64_t[sz];


    uint32_t *params = new uint32_t[1];
    params[0] = sz;

    for (int i = 0; i < sz; i++) {
        double xd = MIN_VAL +
                    (double) (proxy->GenerateCommonRandom() & RAND_MAX) / ((double) (RAND_MAX / (MAX_VAL - MIN_VAL)));
        double yd = MIN_VAL +
                    (double) (proxy->GenerateCommonRandom() & RAND_MAX) / ((double) (RAND_MAX / (MAX_VAL - MIN_VAL)));
        double zd = (double) (proxy->GenerateCommonRandom() & 0x1);
        x[i] = proxy->CreateShare(xd);
        y[i] = proxy->CreateShare(yd);
        z[i] = proxy->CreateShare(zd);
    }

    proxy->SendBytes(coreVectorisedMultiplex, params, 1);
    auto start = chrono::high_resolution_clock::now();
    uint64_t *r = Multiplex(proxy, x, y, z, sz);
    auto end = chrono::high_resolution_clock::now();
    double time_taken = chrono::duration_cast<chrono::nanoseconds>(end - start).count();
    time_taken *= 1e-9;
    cout << "Vectorized MUX Time:\t" << fixed << time_taken << setprecision(9) << " sec" << endl;
    exe_time = time_taken;

    // checking the result
    bool flag = true;
    if(!only_timing) {
        uint64_t *tmp = Reconstruct(proxy, x, sz);
        double *x_reconstructed = ConvertToDouble(tmp, sz);
        delete [] tmp;
        tmp = Reconstruct(proxy, y, sz);
        double *y_reconstructed = ConvertToDouble(tmp, sz);
        delete [] tmp;
        tmp = Reconstruct(proxy, z, sz);
        double *z_reconstructed = ConvertToDouble(tmp, sz);
        delete [] tmp;
        tmp = Reconstruct(proxy, r, sz);
        double *r_reconstructed = ConvertToDouble(tmp, sz);
        delete [] tmp;
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
        delete [] x_reconstructed;
        delete [] y_reconstructed;
        delete [] z_reconstructed;
        delete [] r_reconstructed;
    }
    if (flag)
        cout << "Vectorized Multiplex works correctly" << endl;
    else {
        cout << "Vectorized Multiplex works incorrectly" << endl;
    }
    delete [] params;
    delete [] r;
    return flag;
}

void MAX_Test(Party *proxy) {
    cout << setfill('*') << setw(50) << "Calling Max";
    cout << setfill('*') << setw(49) << "*" << endl;

    uint32_t mRows = WSZ * WSZ;
    uint32_t mCols = WSZ * WSZ;
    uint64_t mSize = mCols * mRows;
    uint32_t *params = new uint32_t[1];
    params[0] = mSize;

    uint64_t *shareOfMatrix = proxy->CreateShare(Random1dData(proxy, mSize), mSize);

    proxy->SendBytes(cnnMax, params, 1);
    uint64_t max = Max(proxy, shareOfMatrix, mSize);

    // checking the result
    double computed_max = -1;
    for (uint32_t position = 0; position < mSize; position++) {
        double matrixVal = ConvertToDouble(Reconstruct(proxy, shareOfMatrix[position]));
        if (matrixVal > computed_max) {
            computed_max = matrixVal;
        }
    }

    double pp_result = ConvertToDouble(Reconstruct(proxy, max));
    if (computed_max == pp_result) {
        cout << "Max works correctly" << endl;
    } else {
        cout << "Max works incorrectly" << endl;
        cout << "computed: " << pp_result << " should be: " << computed_max << endl;
    }

}

void MAX_Specific_Test(Party *proxy) {
    cout << setfill('*') << setw(50) << "Calling Max with specific input";
    cout << setfill('*') << setw(49) << "*" << endl;

    const int kernel = 50; //20;
    const int mRows = 8; //24;
    const int mCols = 8; //24;
    const int channel_size = mCols * mRows;
    const int full_size = kernel * channel_size;

    uint32_t *params = new uint32_t[4];
    params[0] = kernel*mRows;
    params[1] = mCols;
    params[2] = 2;
    params[3] = 2;

    double matrix[kernel][channel_size] = {{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.227232,0,0,0,0,0,0,0,0},
                                           {0,0,0,0,0.0859022,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                           {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                           {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                           {0,0,0,0,0,0,0,0,0,0,0,0,0,0.0752745,0.264879,0.0158701,0.158383,0,0,0,0,0,0,0,0,0,0,0,0,0.100305,0.20839,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0764523,0.0910606,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                           {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                           {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                           {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                           {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                           {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.362823,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.10355,0.355088,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0935564,0.178823,0,0,0,0,0,0,0,0},
                                           {0,0,0,0,0,0,0,0.0194263,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                           {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                           {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                           {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0563087,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.027648,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                           {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                           {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                           {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                           {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                           {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                           {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                            {0,0,0,0,0,0,0,0,0,0,0,0,0.494861,0.118445,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.142336,0.0122185,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0960207,0,0,0},
                                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.198657,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.121045,0.442252,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.169669,0,0,0,0,0,0,0},
                                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0416222,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.2063,0.0495777,0,0,0,0,0,0,0,0,0,0,0,0},
                                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                            {0,0,0,0,0.134719,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0545073,0,0,0,0,0,0,0,0,0,0,0,0.0783863,0,0,0,0.0188293,0,0,0},
                                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0509005,0,0,0,0,0,0,0.235451,0.015522,0,0,0,0,0,0,0.21901,0.0371971,0,0,0,0,0},
                                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.01964,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0757437,0.183009,0.124008,0,0,0,0,0,0,0,0},
                                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                            {0,0.0187225,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.118331,0,0,0,0,0,0,0.232381,0,0,0,0,0,0,0,0.204376,0.483225,0,0,0,0,0,0,0.162704,0,0,0,0,0,0,0,0.494878,0,0,0,0,0,0,0},
                                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.025219,0,0,0,0,0},
                                            {0,0,0,0,0,0,0,0.0165958,0,0,0,0,0,0,0,0.0812492,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}};

    /*double matrix[kernel][channel_size] = {{0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,          0,        0,        0,          0,          0,         0,          0,         0,         0,         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0432825, 0.737419, 1.20078, 0.926585, 0.567595, 0.430356, 0.0985632, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.266747, 0.761726, 1.15512, 1.19953, 1.259, 1.33864, 1.35749, 1.33956, 1.25385, 1.33286, 1.24358, 0.524811, 0.116308, 0.566143, 0.900823, 0.484933, 0.110026, 0.0415068, 0, 0, 0, 0, 0, 0, 0.0496149, 0.298505, 0.423376, 0.433193, 0.560034, 0.81949, 0.932854, 0.937939, 1.0021, 1.00747, 0.514659, 0, 0, 0.593557, 1.18136, 1.09555, 0.56113, 0.121064, 0, 0, 0, 0, 0, 0, 0, 0, 0.0584536, 0.0162401, 0, 0.271658, 0.286316, 0.388371, 0.432861, 0.170568, 0, 0, 0, 0.617371, 1.29066, 0.941217, 0.564678, 0.0604715, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.00190067, 0, 0, 0, 0, 0, 0.725859, 1.19235, 0.914327, 0.357224, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.387172, 1.20801, 1.04574, 0.787889, 0.147707, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.99095, 1.24617, 0.936285, 0.350068, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.480506, 1.18821, 0.922119, 0.67782, 0.136627, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.826264, 1.09702, 0.821839, 0.274472, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.210215, 1.09141, 0.888957, 0.606357, 0.0578222, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.943862, 1.10342, 0.816594, 0.259459, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.551895, 1.31267, 1.04617, 0.583283, 0.0493488, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.192739, 1.38934, 1.32576, 0.856151, 0.220819, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.903785, 1.25367, 0.95167, 0.438142, 0.032342, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.285908, 1.12757, 0.632457, 0.459846, 0.11523, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.729205, 0.850974, 0.368022, 0.153832, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.647827, 1.31563, 1.15062, 0.662708, 0.0616331, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.500298, 2.10337, 2.29506, 1.45795, 0.59522, 0.0153494, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                           {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0734491, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.362925, 0.716253, 0.619985, 0.350005, 0.102894, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.057806, 0.892962, 1.07456, 1.24886, 1.24168, 1.1096, 0.774299, 0.504557, 0.49774, 0.507096, 0.503409, 0.476536, 0.455185, 0.50011, 0.423846, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.234496, 0.457969, 0.802769, 1.17083, 1.26638, 0.953167, 0.951771, 1.02423, 1.00445, 0.9438, 0.882361, 1.28464, 1.56168, 0.91546, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.00835609, 0, 0.629375, 1.53825, 1.23936, 0.182053, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.485014, 0.864035, 0.196854, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.395518, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.268579, 0.206752, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.391879, 0.109731, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.103093, 0.129759, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.186019, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2506, 0.127074, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.00219822, 0, 0, 0, 0.31146, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.247535, 0.141356, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.289079, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.165562, 0.156967, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0153313, 0, 0, 0, 0, 0.0767155, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.340187, 0.574853, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.621376, 0.964777, 0.0631924, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.589972, 0.770303, 0.178358, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                           {0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,          0,        0,        0,          0,          0,         0,          0,         0,         0,         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.237984, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.228874, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.310358, 0, 0, 0, 0, 0, 0, 0, 0, 0.46898, 0.511594, 0, 0, 0, 0, 0, 0, 0, 0.0456448, 0, 0, 0, 0, 0.806241, 0.471326, 0, 0, 0, 0, 0, 0, 0, 0.120534, 0.219846, 0.0505295, 0, 0, 0, 0, 0, 0, 0, 0.0285769, 0, 0, 0, 0.158561, 1.13608, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.12045, 0.84389, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.73989, 1.84284, 0.260638, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.89313, 1.30266, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.05871, 1.89211, 0.0854731, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.93603, 1.18415, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.09261, 1.87932, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.96585, 1.02396, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.53205, 1.72829, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.750181, 1.99137, 0.714845, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.01195, 1.29429, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.53819, 1.91344, 0.0337887, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.206527, 1.85473, 0.873343, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.957108, 1.44047, 0.0771542, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.36761, 0.897438, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.970666, 1.31086, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                           {0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,          0,        0,        0,          0,          0,         0,          0,         0,         0,         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.158424, 0.145867, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0164242, 0.334898, 0.203192, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.140348, 0.189475, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.136819, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.333672, 0.965915, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.0282, 0.563186, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.557622, 1.103, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.988375, 0.485306, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.536626, 0.972125, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.02458, 0.319715, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.793219, 0.83905, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.467673, 1.07768, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0164528, 1.09895, 0.6079, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.800079, 0.880129, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.22402, 0.944175, 0.0867357, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.62048, 0.70813, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.313625, 0.929912, 0.17405, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0700369, 0.832417, 0.970158, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                           {0.280945,  0.280945,  0.280945,  0.280945,  0.280945,  0.280945,  0.280945,  0.280945,  0.280945,  0.280945,  0.280945,  0.280945,  0.280945,  0.280945,  0.280945,  0.280945,  0.280945,  0.280945,  0.280945,  0.280945,  0.280945,  0.280945,  0.280945,  0.280945,  0.280945,  0.280945,  0.280945,  0.280945,  0.280945,  0.280945,  0.280945,  0.280945,  0.280945,  0.280945,  0.280945,  0.280945,  0.280945,  0.280945,  0.280945,  0.280945,  0.280945,  0.280945,  0.280945,  0.280945,  0.280945,  0.280945,  0.280945,  0.280945,  0.280945,  0.280945,  0.280945,  0.280945,  0.280945,  0.280945,  0.280945,  0.280945,  0.280945,  0.280945,  0.280945,  0.280945,  0.280945,  0.280945,  0.280945,  0.280945,  0.280945,  0.280945,  0.280945,  0.280945,  0.280945,  0.280945,  0.280945,  0.280945,  0.280945,  0.280945,  0.207191,  0.0537281, 0,         0.00669956, 0.115641, 0.239457, 0.28947,    0.321052,   0.297474,  0.291881,   0.280945,  0.280945,  0.280945,  0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.0333185, 0, 0, 0, 0, 0, 0, 0, 0.00253105, 0, 0, 0, 0, 0, 0.0136518, 0.138852, 0.277794, 0.328836, 0.330134, 0.296741, 0.280945, 0.280945, 0.280945, 0.280945, 0.0289803, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.069665, 0.300183, 0.338274, 0.316589, 0.280945, 0.280945, 0.280945, 0.280945, 0.0405521, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.227547, 0.375062, 0.310232, 0.280945, 0.280945, 0.280945, 0.280945, 0.11235, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.171496, 0.276066, 0.280945, 0.280945, 0.280945, 0.280945, 0.29567, 0.132909, 0.0365934, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.170278, 0.280945, 0.280945, 0.280945, 0.280945, 0.298882, 0.284643, 0.244589, 0.270036, 0.228269, 0.17215, 0.133272, 0.0651722, 0.0058403, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.109022, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.285496, 0.291809, 0.255734, 0.284547, 0.262274, 0.164466, 0, 0, 0, 0, 0, 0, 0, 0, 0.19263, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.12715, 0, 0, 0, 0, 0, 0, 0, 0.086153, 0.268204, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.273043, 0, 0, 0, 0, 0, 0, 0, 0, 0.189719, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.164668, 0, 0, 0, 0, 0, 0, 0, 0.0711136, 0.249801, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.215094, 0, 0, 0, 0, 0, 0, 0, 0, 0.233369, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.264263, 0.0251904, 0, 0, 0, 0, 0, 0, 0, 0.106101, 0.277406, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.278311, 0.0884714, 0, 0, 0, 0, 0, 0, 0, 0, 0.23989, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.245697, 0, 0, 0, 0, 0, 0, 0, 0, 0.110049, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.253726, 0.0345955, 0, 0, 0, 0, 0, 0, 0, 0.00280762, 0.240598, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.144717, 0, 0, 0, 0, 0, 0, 0, 0, 0.137613, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.227386, 0, 0, 0, 0, 0, 0, 0, 0, 0.0716982, 0.256171, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.13643, 0, 0, 0, 0, 0, 0, 0, 0.0412064, 0.225703, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.0596828, 0, 0, 0, 0, 0, 0, 0, 0.171558, 0.280237, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.0867596, 0, 0, 0, 0, 0, 0, 0, 0.217778, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945, 0.280945},
                                           {0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0.190792,  0.417209,   0.566278, 0.468496, 0.286523,   0.0902662,  0,         0,          0,         0,         0,         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0953951, 0.521488, 0.847225, 0.921467, 0.99016, 0.922757, 0.896797, 0.782505, 0.743255, 0.71077, 0.71077, 0.710768, 0.71077, 0.721519, 0.730688, 0.554837, 0.271, 0.034914, 0, 0, 0, 0, 0, 0, 0, 0.143544, 0, 0, 0, 0.235899, 0.597982, 0.797321, 0.938406, 0.971693, 0.992268, 1.01062, 0.967701, 0.957336, 0.994376, 0.905666, 0.506709, 0.193145, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0749311, 0.190907, 0.0842648, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0885544, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.352369, 0.246889, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0029335, 0, 0, 0, 0, 0, 0, 0.284, 0.804563, 0.0792618, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.170414, 0.654415, 0.52401, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.334496, 0.64537, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0271978, 0.514347, 0.312236, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.302699, 0.637236, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.202598, 0.742551, 0.502673, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.104168, 0.515867, 0.801261, 0.0185318, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0714512, 0.408633, 0.806725, 0.315039, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.133413, 0.691304, 0.421447, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.11162, 0.528804, 0.603087, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.279576, 0.802347, 0.301196, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.165013, 0.583327, 0.780985, 0.0161362, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.290513, 0.605277, 0.375111, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.331939, 0.380505, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                           {0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0.0306339,  0,        0,        0,          0,          0,         0,          0,         0,         0,         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0977488, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0245686, 0.0135546, 0, 0, 0, 0, 0, 0, 0.131175, 0.455553, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.00699043, 0, 0.0303907, 0.0130301, 0, 0, 0, 0.551286, 0.473157, 0, 0, 0, 0, 0, 0, 0, 0.0764084, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.245277, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.150772, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.11565, 0.277762, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.440934, 0.0162191, 0, 0, 0.0222292, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.133646, 0.335653, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.388079, 0.119604, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.193846, 0.23859, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.380293, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.316046, 0.0897198, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.167913, 0.343784, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.44462, 0, 0, 0, 0.0245399, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.346041, 0.152405, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.229223, 0.00677967, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.139061, 0.0957308, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.174007, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0553904, 0.0657215, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                           {0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,          0,        0,        0,          0,          0,         0,          0,         0,         0,         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0886116, 0.133471, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0710545, 0.0765257, 0.174959, 0.218726, 0.112304, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.102848, 0.128772, 0.0220127, 0, 0, 0, 0, 0, 0.181618, 0.227471, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0159473, 0.293757, 0.452892, 0.219244, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.00955486, 0.150662, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.123599, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0955944, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.163246, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.087513, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.127883, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0254564, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.110393, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0509348, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0853271, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0528021, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.223194, 0.377392, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0623112, 0.353663, 0.128427, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0746231, 0.292293, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                           {0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0302534, 0.105946,  0.247674,  0.375164,   0.410138, 0.238834, 0.167532,   0.0291662,  0.0335112, 0.00920105, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0, 0.00174141, 0.203443, 0.606244, 0.867716, 1.02344, 0.976373, 0.806461, 0.579286, 0.459189, 0.425822, 0.425822, 0.425823, 0.425823, 0.425637, 0.399853, 0.276884, 0.128291, 0, 8.96454e-05, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0, 0, 0, 0, 0.118505, 0.488564, 0.907069, 1.1144, 1.08857, 0.997608, 0.863682, 0.910763, 0.890536, 0.869487, 0.906591, 1.04154, 1.0512, 0.78002, 0.368882, 0.0335455, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0262098, 0.417497, 1.03152, 1.07208, 0.738652, 0.247338, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.237531, 0.190769, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0430632, 0.343282, 0.258451, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.397781, 0.749886, 0.858565, 0.603208, 0.210557, 0.0788279, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.126471, 0.292295, 0.382054, 0.405916, 0.382632, 0.573469, 0.817734, 0.707501, 0.585099, 0.438447, 0.387509, 0.329747, 0, 0, 0.0016222, 0, 0, 0, 0, 0, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0542545, 0.149878, 0.167134, 0.219129, 0.188626, 0.219624, 0.128932, 0, 0, 0.0778389, 0.0345478, 0, 0, 0, 0, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0, 0, 0, 0.0489826, 0, 0, 0, 0, 0, 0.00786781, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297594, 0, 0, 0, 0.132257, 0.049449, 0, 0, 0, 0, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0202694, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0301933, 0, 0, 0.149218, 0.25722, 0, 0, 0, 0, 0, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0298262, 0, 0, 0, 0.32396, 0.185249, 0.0293932, 0, 0, 0, 0.0236349, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297203, 0.0120564, 0, 0, 0.0213566, 0.105277, 0.105862, 0, 0, 0, 0, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0265303, 0, 0, 0, 0.116514, 0.0949287, 0, 0, 0, 0, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0299044, 0.0123777, 0, 0, 0, 0.00301552, 0, 0, 0, 0, 0, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0, 0, 0, 0.0695648, 0.0697069, 0.0451555, 0, 0, 0, 0, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0301018, 0, 0, 0, 0.382325, 0.476505, 0.262918, 0, 0, 0, 0, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0, 0, 0, 0, 0.110424, 0.260411, 0.123774, 0, 0, 0, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0284872, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003, 0.0297003},
                                           {0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0.0449495, 0.18644,   0.160052,  0.0618372,  0,        0,        0,          0,          0,         0,          0,         0,         0,         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0712719, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.155187, 0.637244, 1.14965, 1.48843, 1.1999, 0.838457, 0.247337, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.278148, 0.652989, 1.11822, 1.53788, 1.83262, 1.97777, 2.11531, 2.24905, 2.13914, 1.98822, 1.89069, 1.97707, 2.0512, 1.29699, 0.160268, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0261002, 0.127216, 0.193717, 0.351787, 0.490805, 0.655335, 0.895895, 1.12346, 1.33667, 1.46972, 1.49926, 1.69915, 1.25309, 0.325438, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0365677, 0.0302486, 0.138096, 0.179729, 0.289762, 0.295687, 0, 0, 0, 0, 0, 0, 0.0127974, 0.0184717, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.007761, 0, 0, 0, 0, 0, 0, 0, 0.019392, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0346136, 0.0818844, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0793848, 0, 0, 0, 0, 0, 0, 0.070138, 0.141627, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0344944, 0.0430784, 0, 0, 0, 0, 0, 0, 0.0438213, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0899315, 0, 0, 0, 0, 0, 0, 0, 0.0673685, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.153747, 0, 0, 0, 0, 0, 0, 0, 0.0324097, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0156021, 0.0581999, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.150402, 0, 0, 0, 0, 0, 0, 0.0943346, 0.0666065, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0524235, 0, 0, 0, 0, 0, 0, 0.161052, 0.157246, 0.0557833, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0182304, 0, 0, 0, 0, 0, 0, 0, 0.154524, 0.116697, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0190077, 0.00422573, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                           {0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0.0758972, 0.29908,   0.446904,  0.42829,    0.241879, 0.147053, 0.00268841, 0.00120926, 0,         0,          0,         0,         0,         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.159139, 0.34116, 0.611259, 0.609203, 0.894455, 1.0011, 0.89355, 0.75139, 0.597558, 0.572004, 0.522465, 0.522469, 0.522466, 0.522466, 0.487383, 0.307316, 0.0850143, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0200644, 0.218307, 0.49372, 0.844885, 0.945918, 1.02231, 0.945604, 0.927444, 0.919522, 0.886961, 0.864299, 0.86841, 0.913664, 0.930615, 0.650671, 0.333348, 0.162891, 0.0643635, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.107982, 0.258506, 0.419777, 0.354527, 0.36553, 0.55779, 0.790424, 0.949646, 0.704148, 0.46075, 0.361305, 0.163948, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0170355, 0.0402374, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0672741, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0743093, 0.163032, 0.200128, 0.0583572, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0923967, 0.179338, 0.0615358, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.118915, 0.20257, 0.142535, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.06462, 0.236382, 0.28235, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.189349, 0.276633, 0.23079, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.208173, 0.283091, 0.347896, 0.0250111, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0150557, 0.114587, 0.29909, 0.179275, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.00948715, 0.246238, 0.281779, 0.1172, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.1042, 0.19825, 0.354316, 0.141526, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0470781, 0.200229, 0.300705, 0.270118, 0.13003, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0571375, 0.155724, 0.296813, 0.232208, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0787201, 0.120833, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                           {0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,          0,        0,        0,          0,          0,         0,          0,         0,         0,         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                           {0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,          0,        0,        0,          0,          0,         0,          0,         0,         0,         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.152178, 0.764385, 1.38179, 1.2346, 0.99877, 0.85547, 0.527638, 0.126223, 0, 0, 0, 0, 0, 0, 0, 0.0900545, 0, 0, 0, 0, 0, 0, 0, 0.356476, 0.297223, 0.118483, 0.66028, 0.898064, 1.13605, 1.41496, 1.89443, 2.10715, 1.92975, 1.89167, 1.99256, 1.94951, 1.77269, 1.60839, 1.50616, 1.32703, 0.618596, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.229382, 0.35519, 0.317621, 0.473087, 0.379195, 0.20658, 0.0655146, 0.488178, 0.794012, 0.507361, 0.1328, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0229082, 0.00973606, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0368643, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.070797, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0817547, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.165371, 0.0433149, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.00865936, 0.345485, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0901413, 0.143704, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0101328, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.105902, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0488892, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0357666, 0, 0, 0, 0, 0.375725, 0.124421, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.075491, 0, 0, 0, 0.0132027, 0.0968084, 0.114005, 0.0502071, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                           {0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,          0,        0,        0,          0,          0,         0,          0,         0,         0,         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                           {0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,          0,        0,        0,          0,          0,         0,          0,         0,         0,         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.00023365, 0.214173, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.062993, 0.490165, 0.00245857, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.000724792, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.265717, 0.587636, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0510654, 0.0926476, 0.0728989, 0, 0, 0.0150318, 0, 0, 0, 0, 0, 0, 0, 0.626134, 0.532009, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.00258732, 0, 0, 0, 0.0526285, 0.917048, 0.358248, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.569925, 0.747371, 0.0406618, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.802263, 0.358312, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.659881, 0.66299, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.120478, 0.911771, 0.316689, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.697849, 0.654089, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.150175, 0.834767, 0.255295, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.680448, 0.582118, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.344605, 0.719445, 0.200043, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.683311, 0.415438, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.726001, 0.608833, 0.0271015, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.388103, 1.04659, 0.277339, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.797542, 0.949734, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.758975, 0.504703, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.132224, 0.520643, 0.16678, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                           {0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,          0,        0,        0,          0,          0,         0,          0,         0,         0,         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0396738, 0.184656, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.172009, 0.576, 0.259283, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.743513, 0.66107, 0.582356, 0.49578, 0.435237, 0.178647, 0, 0, 0.0282049, 0.0156765, 0, 0.0530138, 0.258425, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.391929, 0.431149, 0.262116, 0.412843, 0.645967, 0.719496, 0.404282, 0.424668, 0.533831, 0.491512, 0.552042, 0.667781, 0.93711, 0.482172, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.276676, 0.152375, 0.12008, 0.236005, 0.296032, 0.486886, 0.500675, 0.765926, 0.261073, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0722322, 0.0237045, 0.0350819, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0987711, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.15524, 0.128613, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.28529, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.145818, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.044652, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0138063, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0931358, 0.219022, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.357976, 0.524323, 0.197508, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.762316, 1.02435, 0.740483, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                           {0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0.164753,  0.0905457, 0,          0,        0,        0,          0,          0,         0,          0,         0,         0,         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.302049, 0.632489, 0.400491, 0.120374, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.406887, 0.830869, 0.602783, 0.171552, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.456203, 0.688563, 0.199932, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.236612, 0.119971, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.292718, 0.975502, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.422639, 1.21495, 0.901145, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.118719, 0.955182, 1.34018, 0.0928049, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.464827, 1.38467, 0.840444, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0360556, 0.933207, 1.31554, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.466526, 1.23895, 0.730035, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.29744, 0.965719, 1.16358, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.179699, 0.804631, 1.18646, 0.522422, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.582205, 1.27272, 0.915565, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.262005, 1.11709, 1.11469, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0789385, 0.708958, 1.4098, 0.496603, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.497803, 1.18945, 1.30817, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.10163, 0.958932, 1.43877, 0.807294, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.307568, 1.3496, 1.26516, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.30358, 0.984871, 0.530167, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                           {0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,          0,        0,        0,          0,          0,         0,          0,         0,         0,         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.120149, 0.400809, 0.376833, 0.14692, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.146519, 0.443768, 0.628616, 0.762361, 0.557635, 0.565026, 0.536958, 0.531481, 0.392584, 0.313382, 0.352374, 0.355478, 0.362147, 0.316507, 0.285995, 0.247098, 0.092598, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.190767, 0.346801, 0.59898, 0.763471, 0.730837, 0.686392, 0.679904, 0.721581, 0.682769, 0.418087, 0.543624, 0.688566, 0.368727, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.269508, 0.482572, 0.220734, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.095892, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0355043, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.220806, 0.504333, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.136187, 0.414672, 0.277226, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.202203, 0.460056, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.321074, 0.129833, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.297476, 0.447349, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.246963, 0.526723, 0.179194, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.107104, 0.491141, 0.407704, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.395029, 0.544072, 0.0770311, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.232516, 0.501597, 0.170622, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0554829, 0.361138, 0.375785, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.350862, 0.37549, 0.098794, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.112092, 0.313046, 0.383846, 0.0441227, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.109595, 0.331562, 0.3494, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0319843, 0.179198, 0.0964727, 0, 0, 0, 0, 0.0191956, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                           {0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0.100121,  0.287907,   0.377955, 0.290541, 0.155945,   0.0130482,  0,         0,          0,         0,         0,         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.288069, 0.767823, 0.862352, 0.553749, 0.223339, 0.20607, 0.212403, 0.338199, 0.38498, 0.446012, 0.446012, 0.446012, 0.446012, 0.45488, 0.470009, 0.356818, 0.138635, 0, 0, 0, 0, 0, 0, 0.063962, 0.488465, 0.571764, 0.267329, 0, 0, 0, 0, 0.207896, 0.364509, 0.502231, 0.526629, 0.533298, 0.503569, 0.471123, 0.377291, 0.0943632, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0568361, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.373775, 0.680672, 0.227518, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.237681, 0.692734, 0.542728, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.372824, 0.627514, 0.00121784, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0903645, 0.576621, 0.540049, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.504531, 0.791862, 0.0847921, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.368144, 0.831078, 0.564584, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.201797, 0.826106, 0.809855, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0122833, 0.628487, 0.806569, 0.197739, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.329732, 0.817646, 0.567723, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.157513, 0.71932, 0.725593, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.529183, 0.761038, 0.471131, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.222958, 0.651875, 0.65927, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.264426, 0.525858, 0.367188, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.133455, 0.106062, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                           {0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,          0,        0,        0,          0,          0,         0,          0,         0,         0,         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2033, 1.06223, 1.19903, 0.611937, 0.00995922, 0, 0, 0, 0, 0, 0, 0.0200768, 0.0300989, 0.0191116, 0.238425, 0.380352, 0.184207, 0, 0, 0, 0, 0, 0, 0.128551, 0.882641, 0.991519, 0.198293, 0, 0, 0, 0.372183, 0.63732, 0.739592, 0.8561, 0.909415, 0.981788, 0.730843, 0.262082, 0.305978, 0.278606, 0, 0, 0, 0, 0, 0, 0, 0.0288038, 0.16503, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.553515, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.50291, 0.796584, 0, 0, 0, 0.000760078, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.06933, 1.43157, 0, 0, 0, 0.113768, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.308579, 1.59907, 0.412732, 0, 0, 0, 0.0452404, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.05658, 1.04179, 0, 0, 0, 0.0774679, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.117063, 1.19715, 0, 0, 0, 0, 0.182841, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.896158, 1.0377, 0, 0, 0, 0.135352, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.384721, 1.55718, 0.217612, 0, 0, 0, 0.0880127, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.41497, 1.0814, 0, 0, 0, 0.135665, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.12927, 1.61058, 0, 0, 0, 0, 0.0763063, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.544815, 1.67397, 0.304942, 0, 0, 0, 0.165591, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.41103, 0.890423, 0, 0, 0, 0.147182, 0.0505705, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.720549, 1.11311, 0, 0, 0, 0.0218458, 0.162202, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.226832, 1.0343, 0.344746, 0, 0, 0, 0.155346, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.736975, 0.892434, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.56616, 1.22317, 0.771925, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                      };*/
    uint64_t *secret = new uint64_t [full_size];
    for (int i = 0; i < full_size; ++i) {
        secret[i] = proxy->CreateShare(matrix[i / channel_size][i % channel_size]);
    }

    proxy->SendBytes(cnnVectorisedMax, params, 4);
    uint64_t *max = Max(proxy, secret, kernel * mRows, mCols, 2, 2);

    // TESTING
    uint32_t window_length = params[2] * params[3];
    uint32_t number_of_windows = full_size / window_length;
    uint64_t *reconstructed_max = Reconstruct(proxy, max, number_of_windows);

    bool flag = true;
    uint64_t resorted[full_size];
    Print1dArray("matrix 6th row: ", matrix[5], channel_size);
    Print1dArray("matrix 23th row: ", matrix[22], channel_size);
    Resort(secret, params[1], params[0], params[3], params[2], resorted);
    double *d_matrix = ConvertToDouble(Reconstruct(proxy, resorted, channel_size), channel_size);
    cout << "resorted matrix 6th row: " << endl;
    for (int i = 0; i < channel_size; ++i) {
        cout << d_matrix[i + 5*channel_size] << " ";
    }
    cout << endl;
    cout << "resorted matrix 23th row: " << endl;
    for (int i = 0; i < channel_size; ++i) {
        cout << d_matrix[i + 22*channel_size] << " ";
    }
    cout << endl;
    double computed_max[number_of_windows];

    for (uint32_t win = 0; win < number_of_windows; win++) {
        computed_max[win] = d_matrix[window_length * win];
        for (uint32_t win_element = 1; win_element < window_length; win_element++) {
            double matrixVal = d_matrix[window_length * win + win_element];
            if (matrixVal > computed_max[win]) {
                computed_max[win] = matrixVal;
            }
        }
        if (computed_max[win] - ConvertToDouble(reconstructed_max[win]) > 0.001) {
            flag = false;
        }
    }


    if (flag) {
        cout << "Vectorized Max works correctly" << endl;
    } else {
        cout << "Vectorized Max works incorrectly" << endl;
        if (params[0] * params[1] < 200) { //otherwise too many values
            Print1dMatrixByWindows("resorted Matrix: ", d_matrix, params[0], params[1], 1,
                                   params[2] * params[3]);
        }
        Print1dMatrixByWindows("computed max values (test): ", computed_max, 1, number_of_windows, 1, 1);
        Print1dMatrixByWindows("VS result from method: ", ConvertToDouble(reconstructed_max, number_of_windows), 1,
                               number_of_windows, 1, 1);
    }

}

void MMAX_Test(Party *proxy) {
    cout << setfill('*') << setw(50) << "Calling vectorized Max";
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

    uint64_t *shareOfMatrix = proxy->CreateShare(Random1dData(proxy, mSize), mSize);
    // PERFORMING MMAX
    proxy->SendBytes(cnnVectorisedMax, mmaxParams, 4);
    uint64_t *max = Max(proxy, shareOfMatrix, mmaxParams[0], mmaxParams[1], mmaxParams[2], mmaxParams[3]);

    // TESTING
    uint32_t window_length = mmaxParams[2] * mmaxParams[3];
    uint32_t number_of_windows = mSize / window_length; // in total, for all of the n_matrices matrices
    uint64_t *reconstructed_max = Reconstruct(proxy, max, number_of_windows);
    bool flag = true;
    uint64_t resorted[mSize];
    Resort(shareOfMatrix, mmaxParams[1], mmaxParams[0], mmaxParams[3], mmaxParams[2], resorted);
    double *d_matrix = ConvertToDouble(Reconstruct(proxy, resorted, mSize), mSize);
    double computed_max[number_of_windows];

    for (uint32_t win = 0; win < number_of_windows; win++) {
        computed_max[win] = d_matrix[window_length * win];
        for (uint32_t win_element = 1; win_element < window_length; win_element++) {
            double matrixVal = d_matrix[window_length * win + win_element];
            if (matrixVal > computed_max[win]) {
                computed_max[win] = matrixVal;
            }
        }
        if (computed_max[win] - ConvertToDouble(reconstructed_max[win]) > 0.001) {
            flag = false;
        }
    }

    if (flag) {
        cout << "Vectorized Max works correctly" << endl;
    } else {
        cout << "Vectorized Max works incorrectly" << endl;
        if (mmaxParams[0] * mmaxParams[1] < 200) { //otherwise too many values
            Print1dMatrixByWindows("resorted Matrix: ", d_matrix, mmaxParams[0], mmaxParams[1], 1,
                                   mmaxParams[2] * mmaxParams[3]);
        }
        Print1dMatrixByWindows("computed max values (test): ", computed_max, 1, number_of_windows, 1, 1);
        Print1dMatrixByWindows("VS result from method: ", ConvertToDouble(reconstructed_max, number_of_windows), 1,
                               number_of_windows, 1, 1);
    }
}

void ARGMAX_Test(Party *proxy) {
    cout << setfill('*') << setw(50) << "Calling ArgMax";
    cout << setfill('*') << setw(49) << "*" << endl;

    uint32_t mRows = WSZ / 2;
    uint32_t mCols = WSZ / 2;
    uint64_t mSize = mCols * mRows;
    uint32_t *params = new uint32_t[1];
    params[0] = mSize;

    uint64_t *shareOfMatrix = proxy->CreateShare(Random1dData(proxy, mSize), mSize);

    proxy->SendBytes(cnnArgMax, params, 1);
    uint64_t argmax = ArgMax(proxy, shareOfMatrix, mSize);

    // checking the result
    double computed_max = -1;
    double computed_argmax = -1;
    for (uint32_t position = 0; position < mSize; position++) {
        double matrixVal = ConvertToDouble(Reconstruct(proxy, shareOfMatrix[position]));
        if (matrixVal > computed_max) {
            computed_max = matrixVal;
            computed_argmax = position;
        }
    }

    double pp_result = ConvertToDouble(Reconstruct(proxy, argmax));
    if (computed_argmax == pp_result) {
        cout << "ArgMax works correctly" << endl;
    } else {
        cout << "ArgMax works incorrectly" << endl;
        cout << "computed: " << pp_result << " should be: " << computed_argmax << endl;
    }
}

void RST_Test(Party *proxy) {
    cout << setfill('*') << setw(50) << "Calling ReSort";
    cout << setfill('*') << setw(49) << "*" << endl;

    const int kernel = 50;
    const int mRows = 8; //WSZ * 10;
    const int mCols = 8; // WSZ * 10;
    const int mSize = mCols * mRows;

    uint32_t wRows = 2; //WSZ;
    uint32_t wCols = 2; //WSZ;
    double matrix[kernel][mSize] = {{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.227232,0,0,0,0,0,0,0,0},
                                           {0,0,0,0,0.0859022,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                           {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                           {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                           {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                           {0,0,0,0,0,0,0,0,0,0,0,0,0,0.0752745,0.264879,0.0158701,0.158383,0,0,0,0,0,0,0,0,0,0,0,0,0.100305,0.20839,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0764523,0.0910606,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                           {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                           {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                           {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                           {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                           {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.362823,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.10355,0.355088,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0935564,0.178823,0,0,0,0,0,0,0,0},
                                           {0,0,0,0,0,0,0,0.0194263,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                           {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                           {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                           {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0563087,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.027648,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                           {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                           {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                           {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                           {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                           {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                           {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                           {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                           {0,0,0,0,0,0,0,0,0,0,0,0,0.494861,0.118445,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.142336,0.0122185,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                           {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                           {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0960207,0,0,0},
                                           {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                           {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                           {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                           {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.198657,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.121045,0.442252,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.169669,0,0,0,0,0,0,0},
                                           {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                           {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0416222,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.2063,0.0495777,0,0,0,0,0,0,0,0,0,0,0,0},
                                           {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                           {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                           {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                           {0,0,0,0,0.134719,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0545073,0,0,0,0,0,0,0,0,0,0,0,0.0783863,0,0,0,0.0188293,0,0,0},
                                           {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                           {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                           {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0509005,0,0,0,0,0,0,0.235451,0.015522,0,0,0,0,0,0,0.21901,0.0371971,0,0,0,0,0},
                                           {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                           {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.01964,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0757437,0.183009,0.124008,0,0,0,0,0,0,0,0},
                                           {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                           {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                           {0,0.0187225,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.118331,0,0,0,0,0,0,0.232381,0,0,0,0,0,0,0,0.204376,0.483225,0,0,0,0,0,0,0.162704,0,0,0,0,0,0,0,0.494878,0,0,0,0,0,0,0},
                                           {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                           {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                           {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.025219,0,0,0,0,0},
                                           {0,0,0,0,0,0,0,0.0165958,0,0,0,0,0,0,0,0.0812492,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                           {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                           {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                           {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}};
    cout << "create shares... " << endl;
    uint64_t *shareOfMatrix = new uint64_t [kernel*mSize];
    for (int i = 0; i < kernel; ++i) {
        uint64_t *share = proxy->CreateShare(matrix[i], mSize);
        for (int j = 0; j < mSize; ++j) {
            shareOfMatrix[i*mSize + j] = share[j];
        }
        delete[] share;
    }

    auto *resorted = new uint64_t[kernel*mSize];
    //Print1dMatrixByWindows("ReSort original", ConvertToDouble(Reconstruct(proxy, shareOfMatrix, mSize), mSize), mRows, mCols, wRows, wCols);
    cout << "resort" << endl;
    Resort(shareOfMatrix, mCols, kernel * mRows, wCols, wRows, resorted);

    uint64_t wElements = wRows * wCols;
    uint64_t numberOfWins = kernel * mSize / wElements;
    Print1dMatrixByWindows("ReSort finished: resorted", ConvertToDouble(Reconstruct(proxy, resorted, mSize), mSize),
                           numberOfWins, wElements, 1, 1);
    if (mSize < 100) {
        Print1dMatrixByWindows("ReSort finished: resorted", ConvertToDouble(Reconstruct(proxy, resorted, mSize), mSize),
                               numberOfWins, wElements, 1, 1);
    }
}

void RELU_Test(Party *proxy) {
    cout << setfill('*') << setw(50) << "Calling Relu";
    cout << setfill('*') << setw(49) << "*" << endl;

    uint64_t x = proxy->CreateShare(ConvertToDouble(proxy->GenerateCommonRandom()));
    proxy->SendBytes(cnnRelu);
    uint64_t relu = ReLU(proxy, x);
    uint64_t reconstructed_relu = Reconstruct(proxy, relu);

    // checking the result
    double computed_relu = -1;
    double originalX = ConvertToDouble(Reconstruct(proxy, x));
    if (originalX >= 0) {
        computed_relu = originalX;
    } else {
        computed_relu = 0;
    }

    double pp_result = ConvertToDouble(reconstructed_relu);
    if (computed_relu == pp_result) {
        cout << "Relu works correctly" << endl;
    } else {
        cout << "Relu works incorrectly" << endl;
        cout << "computed: " << pp_result << " should be: " << computed_relu << endl;
    }

}

void MRELU_Test(Party *proxy) {
    cout << setfill('*') << setw(50) << "Calling MRELU";
    cout << setfill('*') << setw(49) << "*" << endl;

    uint64_t size = WSZ * WSZ;
    uint32_t *params = new uint32_t[1];
    params[0] = size;
    uint64_t *x = proxy->CreateShare(Random1dData(proxy, size, 255.0), size);
    proxy->SendBytes(cnnVectorisedRelu, params, 1);
    uint64_t *mrelu = Relu(proxy, x, size);
    uint64_t *rec_mrelu = Reconstruct(proxy, mrelu, size);

    // checking the result
    double *correct_relu = new double[size];
    double *originalX = ConvertToDouble(Reconstruct(proxy, x, size), size);
    double *pp_result = ConvertToDouble(rec_mrelu, size);
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
        Print1dArray("Original Values:", originalX, size);
        Print1dArray("Computed Relu:", pp_result, size);
        Print1dArray("VS Correct Relu:", correct_relu, size);
    }

}

void DRLU_Test(Party *proxy) {
    cout << setfill('*') << setw(50) << "Calling DRLU";
    cout << setfill('*') << setw(49) << "*" << endl;

    uint64_t x = proxy->CreateShare(ConvertToDouble(proxy->GenerateRandom()));

    proxy->SendBytes(cnnDerivativeRelu);
    uint64_t drelu = DerivativeRelu(proxy, x);
    uint64_t reconstructed_drelu = Reconstruct(proxy, drelu);

    // checking the result
    double originalX = ConvertToDouble(Reconstruct(proxy, x));
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
    uint64_t *x = proxy->CreateShare(Random1dData(proxy, size), size);

    proxy->SendBytes(cnnVectorisedDerivativeRelu, params, 1);
    uint64_t *drelu = DerivativeRelu(proxy, x, size);
    double *rec_drelu = ConvertToDouble(Reconstruct(proxy, drelu, size), size);

    // checking the result
    double *originalX = ConvertToDouble(Reconstruct(proxy, x, size), size);
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
        Print1dArray("Original Values:", originalX, size);
        Print1dArray("Computed DerivativeRelu:", rec_drelu, size);
        Print1dArray("VS Correct DerivativeRelu:", correct_drelu, size);
    }

}

void PAD_Test(Party *proxy) {
    cout << setfill('*') << setw(50) << "Calling Pad (padding)";
    cout << setfill('*') << setw(49) << "*" << endl;

    uint32_t row = 28;
    uint32_t col = 28;
    uint32_t padding_value = 0;
    uint32_t padding_size = 2;
    uint64_t **x = proxy->CreateShare(Random2dData(proxy, row, col, 0.0, 255.0), row, col);
    //Print2dArray("Original matrix: ", ConvertToDouble(Reconstruct(proxy, x, row, col), row, col), row, col);

    uint64_t **padded = Pad(x, row, col, padding_value, padding_size);
    double **rec_padded = ConvertToDouble(Reconstruct(proxy, padded, row + 2 * padding_size, col + 2 * padding_size),
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
            computed_padding[padding_size + r][c + padding_size] = ConvertToDouble(Reconstruct(proxy, x[r][c]));
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
        cout << "Pad works correctly" << endl;
        //Print2dArray("Computed Padding", rec_padded, row + 2 * padding_size, col + 2 * padding_size);
    } else {
        cout << "Pad works incorrectly" << endl;
        Print2dArray("Computed Padding", rec_padded, row + 2 * padding_size, col + 2 * padding_size);
        Print2dArray("Correct Padding", computed_padding, row + 2 * padding_size, col + 2 * padding_size);
    }

}

void CL_Test(Party *proxy) {
    cout << setfill('*') << setw(50) << "Calling ConvolutionalLayer (convolutional layer)";
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
        x[i] = proxy->CreateShare(Random2dData(proxy, row, col, -10.0, 10.0), row, col);
        //Print2dArray("input", ConvertToDouble(Reconstruct(proxy, x[i], row, col), row, col), row, col);
        kernel[i] = proxy->CreateShare(Random2dData(proxy, k_number, k_dim * k_dim, -5.0, 5.0), k_number,
                                       k_dim * k_dim);
        //Print2dArray("kernel", ConvertToDouble(Reconstruct(proxy, kernel[i], k_number, k_dim*k_dim), k_number, k_dim*k_dim), k_number, k_dim*k_dim);
    }
    uint64_t conv_width = floor((col - k_dim + 1) / stride);
    uint64_t conv_height = floor((row - k_dim + 1) / stride);

    double *rec_bias = Random1dData(proxy, k_number);
    uint64_t *bias = proxy->CreateShare(rec_bias, k_number);
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
    proxy->SendBytes(cnnConvolutionalLayer, params, 9);

    uint64_t ***conv = ConvolutionalLayer(proxy, x, i_channel, row, col, kernel, k_dim, k_number, stride, max_height,
                                          max_width,
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
    //Print2dArray("ConvolutionalLayer result, first output channel", ConvertToDouble(Reconstruct(proxy, conv[0], out_height, out_width), out_height, out_width), out_height, out_width);
    //Print2dArray("ConvolutionalLayer result, second output channel", ConvertToDouble(Reconstruct(proxy, conv[1], out_height, out_width), out_height, out_width), out_height, out_width);

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
                        double v = ConvertToDouble(
                                Reconstruct(proxy, x[channel][cr + k_value / k_dim][cc + k_value % k_dim]));
                        double weight = ConvertToDouble(Reconstruct(proxy, kernel[channel][kern][k_value]));
                        dot_product += v * weight;
                    }
                }
                dot_product += rec_bias[kern];
                // Activation: Relu
                double relu = 0.0;
                if (dot_product > 0) {
                    relu = dot_product;
                }
                correct_conv[kern][cr / stride][cc / stride] = relu;
                if (max_height <= 1 and max_width <= 1) {
                    // no maxpooling --> check if correct
                    double comp_value = 0;
                    if (flatten_result) {
                        comp_value = ConvertToDouble(Reconstruct(proxy, conv[0][0][kern * conv_height * conv_width +
                                                                                  (cr / stride * conv_width) +
                                                                                  cc / stride]));
                    } else {
                        comp_value = ConvertToDouble(Reconstruct(proxy, conv[kern][cr / stride][cc / stride]));
                    }
                    if (abs(relu - comp_value) > 0.0001) {
                        allCorrect = false;
                    }
                }
            }
        }
        //Print2dArray("activated conv", correct_conv[kern], conv_height, conv_width);
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
                        comp_max = ConvertToDouble(Reconstruct(proxy, conv[0][0][kern * max_height * max_width +
                                                                                (r / max_height * max_width) +
                                                                                c / max_width]));
                    } else {
                        pooled_conv[kern][r / max_height][c / max_width] = max;
                        comp_max = ConvertToDouble(Reconstruct(proxy, conv[kern][r / max_height][c / max_width]));
                    }
                    if (abs(max - comp_max) > 0.001) {
                        allCorrect = false;
                    }
                }
            }
        }
    }
    if (allCorrect) {
        cout << "ConvolutionalLayer works correctly" << endl;
    } else {
        cout << "ConvolutionalLayer works incorrectly" << endl;
        /*cout << setfill('=') << setw(25) << "GIVEN PARAMETERS";
        cout << setfill('=') << setw(24) << "=" << endl;
        cout << "Input channel of X" << endl;
        for (int i = 0; i < i_channel; ++i) {
            Print2dArray("x", ConvertToDouble(Reconstruct(proxy, x[i], row, col), row, col), row, col);
            cout << "Output channel of kernel " << endl;
            Print2dArray("Channel of kernel",
                         ConvertToDouble(Reconstruct(proxy, kernel[i], k_number, k_dim * k_dim), k_number, k_dim * k_dim), k_number, k_dim * k_dim);
        }*/
        cout << setfill('=') << setw(25) << "OUTPUT PARAMETERS";
        cout << setfill('*') << setw(24) << "*" << endl;
        if (flatten_result) {
            if (max_width > 1 and max_height > 1) {
                Print1dArray("Flattened, maxpooled computed convolution",
                             ConvertToDouble(Reconstruct(proxy, conv[0][0], out_height * out_width * k_number),
                                            out_height * out_width * k_number), out_height * out_width * k_number);
                for (int k = 0; k < k_number; ++k) {
                    Print1dArray("Correct, pooled convolution", pooled_conv[k][0], out_height * out_width);
                }
            } else {
                Print1dArray("Flattened computed convolution",
                             ConvertToDouble(Reconstruct(proxy, conv[0][0], conv_height * conv_width * k_number),
                                            conv_height * conv_width * k_number), conv_height * conv_width * k_number);
                for (int k = 0; k < k_number; ++k) {
                    Print1dArray("Correct, pooled convolution", pooled_conv[k][0], conv_height * conv_width);
                }
            }
        } else {
            for (uint64_t k = 0; k < k_number; k++) {
                cout << "Output channel " << k << endl;
                if (max_width > 0 and max_height > 0) {
                    Print2dArray("Computed convolution",
                                 ConvertToDouble(Reconstruct(proxy, conv[k], out_height, out_width), out_height, out_width),
                                 out_height, out_width);
                    Print2dArray("Correct, pooled convolution", pooled_conv[k], out_height, out_width);
                } else {
                    Print2dArray("Computed convolution",
                                 ConvertToDouble(Reconstruct(proxy, conv[k], conv_height, conv_width), conv_height, conv_width),
                                 conv_height, conv_width);
                    Print2dArray("Correct convolution", correct_conv[k], conv_height, conv_width);
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

bool FCL_Test(Party *proxy) {
    cout << setfill('*') << setw(50) << "Calling FullyConnectedLayer (fully connected layer)";
    cout << setfill('*') << setw(49) << "*" << endl;

    // init
    uint32_t row = WSZ*2;
    uint32_t col = WSZ*2;
    uint32_t i_channel = 50;
    uint32_t i_nodes = row * col * i_channel;
    uint32_t o_nodes = WSZ;

    bool allCorrect = true;
    uint64_t *x = proxy->CreateShare(Random1dData(proxy, i_nodes, 0.0, 255.0), i_nodes);
    uint64_t **weights = proxy->CreateShare(Random2dData(proxy, i_nodes, o_nodes, -5.0, 5.0), i_nodes, o_nodes);
    uint64_t *bias = proxy->CreateShare(Random1dData(proxy, o_nodes), o_nodes);

    // send params
    uint32_t params[2];
    params[0] = i_nodes;
    params[1] = o_nodes;
    proxy->SendBytes(cnnFullyConnectedLayer, params, 2);

    uint64_t *res = FullyConnectedLayer(proxy, x, i_nodes, weights, o_nodes, bias);
    double *reconstructed_res = ConvertToDouble(Reconstruct(proxy, res, o_nodes), o_nodes);

    // checking the result
    double *rec_input = ConvertToDouble(Reconstruct(proxy, x, i_nodes), i_nodes);
    delete[] x;
    double **rec_weights = ConvertToDouble(Reconstruct(proxy, weights, i_nodes, o_nodes), i_nodes, o_nodes);
    for (int i = 0; i < i_nodes; ++i) {
        delete[] weights[i];
    }
    delete[] weights;
    double *rec_bias = ConvertToDouble(Reconstruct(proxy, bias, o_nodes), o_nodes);
    delete[] bias;

    double *original_res = new double[o_nodes];
    for (uint64_t o = 0; o < o_nodes; o++) {
        double value = 0;
        for (int i = 0; i < i_nodes; i++) {
            value += rec_input[i] * rec_weights[o][i]; //weights would be tansposed in FullyConnectedLayer
        }
        //Relu activation
        if (value < 0) {
            value = 0;
        }
        original_res[o] = value + rec_bias[o];
        if (original_res[o] - reconstructed_res[o] > 0.0001) {
            allCorrect = false;
        }
    }

    if (allCorrect) {
        cout << "FullyConnectedLayer works correctly" << endl;
    } else {
        cout << "FullyConnectedLayer works incorrectly" << endl;
        Print1dArray("Computed Result: ", reconstructed_res, o_nodes);
        Print1dArray("Correct Result", original_res, o_nodes);
        Print1dArray("input", rec_input, i_nodes);
        Print2dArray("weights", rec_weights, i_nodes, o_nodes);
        Print1dArray("bias", rec_bias, o_nodes);
    }
    delete[] original_res;
    delete[] reconstructed_res;
    for (int i = 0; i < i_nodes; ++i) {
        delete[] weights[i];
    }
    delete[] rec_weights;
    delete[] bias;
    return allCorrect;
}

bool FLT_Test(Party *proxy) {
    cout << setfill('*') << setw(50) << "Calling Flatten (flattening)";
    cout << setfill('*') << setw(49) << "*" << endl;

    uint32_t row = 14;
    uint32_t col = 14;
    uint32_t i_channel = WSZ/2;

    uint64_t ***x = new uint64_t **[i_channel];
    for (uint64_t i = 0; i < i_channel; i++) {
        x[i] = proxy->CreateShare(Random2dData(proxy, row, col, 0, 255), row, col);
    }
    uint32_t i_nodes = row * col * i_channel;

    uint64_t *res = Flatten(x, row, col, i_channel);
    double *reconstructed_res = ConvertToDouble(Reconstruct(proxy, res, i_nodes), i_nodes);
    delete[] res;

    // checking the result
    bool allCorrect = true;
    double *original_res = new double[i_nodes];
    for (int channel = 0; channel < i_channel; ++channel) {
        for (int r = 0; r < row; ++r) {
            for (int c = 0; c < col; ++c) {
                uint64_t position = channel * row * col + r * col + c;
                original_res[position] = ConvertToDouble(Reconstruct(proxy, x[channel][r][c]));
                if ((original_res[position] - reconstructed_res[position]) > 0.0001) {
                    allCorrect = false;
                }
            }
        }
    }

    if (allCorrect) {
        cout << "Flatten works correctly" << endl;
    } else {
        cout << "Flatten works incorrectly" << endl;
        Print1dArray("Computed Result: ", reconstructed_res, i_nodes);
        Print1dArray("Correct Result", original_res, i_nodes);
    }
    delete[] original_res;
    delete[] reconstructed_res;
    for (int i = 0; i < i_channel; ++i) {
        delete[] x[i];
    }
    delete[] x;
    return allCorrect;
}

void INC_Test(Party *proxy) {
    cout << setfill('*') << setw(50) << "Calling Increase (increasing input matrix for conv)";
    cout << setfill('*') << setw(49) << "*" << endl;

    uint32_t row = 24;
    uint32_t col = 24;
    uint32_t i_channel = 20;
    uint64_t ***x = new uint64_t **[i_channel];
    for (uint64_t i = 0; i < i_channel; i++) {
        x[i] = proxy->CreateShare(Random2dData(proxy, row, col, 0, 255), row, col);
    }
    //Print2dArray("Original X: ", ConvertToDouble(Reconstruct(proxy, x[0], row, col), row, col), row, col);

    uint32_t k_number = 5;
    uint32_t k_dim = 5;
    uint32_t k_len = k_dim * k_dim;
    uint32_t stride = 1;

    uint32_t lastpos = row - k_dim + 1;
    uint32_t conv_height = static_cast<uint32_t>(floor(lastpos / stride));
    uint32_t conv_width = static_cast<uint32_t>(floor((col - k_dim + 1) / stride));

    uint64_t ***stretchedX = Increase(x, i_channel, row, col, k_dim, stride);
    // checking the result
    if (row * col < 100) {
        Print2dArray("RESULTED STRETCHED MATRIX (first channel)",
                     ConvertToDouble(Reconstruct(proxy, stretchedX[0], conv_height * conv_width, k_len),
                                    conv_height * conv_width, k_len), conv_height * conv_width, k_len);
    }
}

void EXP_Test(Party *proxy) {
    cout << setfill('*') << setw(50) << "Calling Exp";
    cout << setfill('*') << setw(49) << "*" << endl;

    cout << "Min power: " << proxy->GetMinPower() << endl;
    cout << "Max power: " << proxy->GetMaxPower() << endl;

    uint64_t x = proxy->CreateShare(Random1dData(proxy, 1, proxy->GetMinPower(), proxy->GetMaxPower())[0]);

    proxy->SendBytes(coreExp);
    uint64_t shr_exp = Exp(proxy, x);
    uint64_t reconstructed_exp = Reconstruct(proxy, shr_exp);
    double rec_exp = ConvertToDouble(reconstructed_exp);

    // checking the result
    double originalX = ConvertToDouble(Reconstruct(proxy, x));
    double true_exp = exp(originalX);

    if (abs(true_exp - rec_exp) <= 0.1) {
        cout << "Exp works correctly" << endl;
        cout << "power: " << originalX << " -- computed: " << rec_exp << " should be: " << true_exp << endl;
    } else {
        cout << "Exp works incorrectly" << endl;
        cout << "power: " << originalX << " -- computed: " << rec_exp << " should be: " << true_exp << endl;
    }

}

bool MEXP_Test(Party *proxy, double &exe_time, int &cnt, bool only_timing = false) {
    cout << setfill('*') << setw(50) << "Calling MEXP";
    cout << setfill('*') << setw(49) << "*" << endl;

    uint32_t params[1];
    int n_samples = sz;
    params[0] = n_samples;
    // random values
    uint64_t *x = proxy->CreateShare(Random1dData(proxy, n_samples, proxy->GetMinPower() + 10, proxy->GetMaxPower() - 10), n_samples);

    // specific values
//    int n_samples = 6;
//    params[0] = n_samples;
//    double tmp_test_values[6] = {-2.08377, -3.29118, -2.71972, -3.9096, -3.63482, -2.61542};
//    uint64_t *x = proxy->CreateShare(tmp_test_values, n_samples);

    proxy->SendBytes(coreVectorisedExp, params, 1);
    auto start = chrono::high_resolution_clock::now();
    uint64_t *shr_exp = Exp(proxy, x, n_samples);
    auto end = chrono::high_resolution_clock::now();
    double time_taken = chrono::duration_cast<chrono::nanoseconds>(end - start).count();
    time_taken *= 1e-9;
    cout << "Vectorized EXP Time:\t" << fixed << time_taken << setprecision(9) << " sec" << endl;
    exe_time = time_taken;

    bool flag = true;
    if(!only_timing) {
        uint64_t *reconstructed_exp = Reconstruct(proxy, shr_exp, n_samples);
        double *rec_exp = ConvertToDouble(reconstructed_exp, n_samples);

        // checking the result
        double *originalX = ConvertToDouble(Reconstruct(proxy, x, n_samples), n_samples);
        double *true_exp = new double[n_samples];
        for (int i = 0; i < n_samples; i++) {
            true_exp[i] = exp(originalX[i]);
        }

        double max_diff = MIN_VAL;
        double tmp_power;
        double tmp_true_exp;
        double tmp_rec_exp;
        double tmp;
        for (int i = 0; i < n_samples; i++) {
            tmp = abs(rec_exp[i] - true_exp[i]);
            cout << fixed << "power: " << originalX[i] << " -- computed: " << rec_exp[i] << " - expected: " << true_exp[i]
                 <<
                 " - absolute difference: " << tmp << endl;
            if(tmp > max_diff) {
                max_diff = tmp;
                tmp_power = originalX[i];
                tmp_true_exp = true_exp[i];
                tmp_rec_exp = rec_exp[i];
            }
            if (abs((true_exp[i] - rec_exp[i]) * 100.0 / abs(true_exp[i])) >= 1) {
                flag = false;
                cnt++;
            }
        }
        if (flag) {
            cout << "EXP works correctly" << endl;
        } else {
            cout << "EXP works incorrectly" << endl;
        }
        cout << "Maximum absolute difference is " << max_diff << " and obtained for " << tmp_power << endl;
        cout << fixed << "Computed: " << tmp_rec_exp << " - expected: " << tmp_true_exp << endl;
    }
    return flag;
}

void DP_Test(Party *proxy) {
    cout << setfill('*') << setw(50) << "Calling DotProduct";
    cout << setfill('*') << setw(49) << "*" << endl;

    uint32_t *params = new uint32_t[1];
    uint32_t size = 20; // size of the vector
    params[0] = size;
    double min_val = -10;
    double max_val = 10;

    // generate a vector of random values
    uint64_t *vec1 = proxy->CreateShare(Random1dData(proxy, size, min_val, max_val), size);
    uint64_t *vec2 = proxy->CreateShare(Random1dData(proxy, size, min_val, max_val), size);

    // call DotProduct
    proxy->SendBytes(coreDotProduct, params, 1);
    uint64_t res = DotProduct(proxy, vec1, vec2, size);

    // check the results
    double gt = 0;
    double *rec_vec1 = ConvertToDouble(Reconstruct(proxy, vec1, size), size);
    double *rec_vec2 = ConvertToDouble(Reconstruct(proxy, vec2, size), size);
    for (int i = 0; i < size; i++) {
        gt += rec_vec1[i] * rec_vec2[i];
    }

    double rec_res = ConvertToDouble(Reconstruct(proxy, res));

    PrintValue("Computed dot product", rec_res);
    PrintValue("GT dot product", gt);
}

bool MDP_Test(Party *proxy, double &exe_time, bool only_timing = false) {
    cout << setfill('*') << setw(50) << "Calling MDP";
    cout << setfill('*') << setw(49) << "*" << endl;

    uint32_t *params = new uint32_t[1];
    params[0] = sz;

    uint32_t d = 5; // size of each individual vector in the main vector
    double min_val = -10;
    double max_val = 10;

    if (sz % d != 0) {
        throw invalid_argument("DP_Test: The size must be divisible by d.");
    }

    // generate a vector of random values
    uint64_t *vec1 = new uint64_t[sz];
    uint64_t *vec2 = new uint64_t[sz];
    for (int i = 0; i < sz; i++) {
        double xd = MIN_VAL + (double) (proxy->GenerateCommonRandom() & RAND_MAX) / ((double) (RAND_MAX / (MAX_VAL - MIN_VAL)));
        double yd = MIN_VAL + (double) (proxy->GenerateCommonRandom() & RAND_MAX) / ((double) (RAND_MAX / (MAX_VAL - MIN_VAL)));
        double zd = (double) (proxy->GenerateCommonRandom() & 0x1);
        vec1[i] = proxy->CreateShare(xd);
        vec2[i] = proxy->CreateShare(yd);
    }

    // call DotProduct
    proxy->SendBytes(coreVectorisedDotProduct, params, 1);
    auto start = chrono::high_resolution_clock::now();
    uint64_t *res = DotProduct(proxy, vec1, vec2, sz, d);
    auto end = chrono::high_resolution_clock::now();
    double time_taken = chrono::duration_cast<chrono::nanoseconds>(end - start).count();
    time_taken *= 1e-9;
    cout << "Vectorized DP Time:\t" << fixed << time_taken << setprecision(9) << " sec" << endl;
    exe_time = time_taken;

    // check the results
    bool flag = true;
    if(!only_timing) {
        double *gt = new double[sz / d];
        double *rec_vec1 = ConvertToDouble(Reconstruct(proxy, vec1, sz), sz);
        double *rec_vec2 = ConvertToDouble(Reconstruct(proxy, vec2, sz), sz);
        for (int i = 0; i < sz; i += d) {
            double tmp_sum = 0;
            for (int j = i; j < i + d; j++) {
                tmp_sum += rec_vec1[j] * rec_vec2[j];
            }
            gt[i / d] = tmp_sum;
        }

        double *rec_res = ConvertToDouble(Reconstruct(proxy, res, sz / d), sz / d);

        Print1dArray("Computed dot product", rec_res, sz / d);
        Print1dArray("GT dot product", gt, sz / d);
    }
    return flag;
}

bool MATMATMUL_Test(Party *proxy) {
    cout << setfill('*') << setw(50) << "Calling MatrixMatrixMultiply";
    cout << setfill('*') << setw(49) << "*" << endl;

    // setting
    int a_row = WSZ*3;
    int a_col = WSZ*2;
    int b_col = WSZ*4;

    uint32_t *params = new uint32_t[1];
    uint32_t size = a_row * a_col * b_col; // size of the vector
    params[0] = size;
    double min_val = -0.784378;
    double max_val = 1481.76;
//    double max_val = 1000;

    double **data1 = Random2dData(proxy, a_row, a_col, min_val, max_val);
    double **data2 = Random2dData(proxy, a_col, b_col, min_val, max_val);
    uint64_t **mat1 = proxy->CreateShare(data1, a_row, a_col);
    uint64_t **mat2 = proxy->CreateShare(data2, a_col, b_col);

    proxy->SendBytes(coreMatrixMatrixMultiply, params, 1);
    uint64_t **res = MatrixMatrixMultiply(proxy, mat1, mat2, a_row, a_col, b_col);
    uint64_t **tmp_rec = Reconstruct(proxy, res, a_row, b_col);
    double **rec_res = ConvertToDouble(tmp_rec, a_row, b_col);

    uint64_t** tmp_rec_mat1 = Reconstruct(proxy, mat1, a_row, a_col);
    uint64_t** tmp_rec_mat2 = Reconstruct(proxy, mat2, a_col, b_col);
    double **rec_mat1 = ConvertToDouble(tmp_rec_mat1, a_row, a_col);
    double **rec_mat2 = ConvertToDouble(tmp_rec_mat2, a_col, b_col);

//    proxy->Print2dArray("mat1", rec_mat1, a_row, a_col);
//    proxy->Print2dArray("mat2", rec_mat2, a_col, b_col);

    double **gt = MultiplyMatrices(rec_mat1, rec_mat2, a_row, a_col, b_col);

    double tmp = 0;
    for (int i = 0; i < a_row; i++) {
        for (int j = 0; j < b_col; j++) {
            tmp += abs(rec_res[i][j] - gt[i][j]);
        }
    }

    if (tmp <= 0.1) {
        cout << "MatrixMatrixMultiply works correctly" << endl;
        //cout << "Total absolute difference: " << tmp << endl;
    } else {
        cout << "MatrixMatrixMultiply works incorrectly" << endl;
        cout << "Total absolute difference: " << tmp << endl;

        Print2dArray("mat", rec_mat1, a_row, a_col);
        Print2dArray("mat", rec_mat2, a_col, b_col);

        Print2dArray("Computed matrix multiplication", rec_res, a_row, b_col);
        Print2dArray("GT matrix multiplication", gt, a_row, b_col);
    }

    delete [] params;
    for( int i = 0; i < a_row; i++) {
        delete [] mat1[i];
        delete [] res[i];
        delete [] rec_res[i];
        delete [] rec_mat1[i];
        delete [] gt[i];
        delete [] tmp_rec[i];
        delete [] tmp_rec_mat1[i];
        delete [] data1[i];
    }
    delete [] mat1;
    delete [] res;
    delete [] rec_res;
    delete [] rec_mat1;
    delete [] gt;
    delete [] tmp_rec;
    delete [] tmp_rec_mat1;
    delete [] data1;
    for( int i = 0; i < a_col; i++) {
        delete [] mat2[i];
        delete [] rec_mat2[i];
        delete [] tmp_rec_mat2[i];
        delete [] data2[i];
    }
    delete [] mat2;
    delete [] rec_mat2;
    delete [] tmp_rec_mat2;
    delete [] data2;

    return (tmp <= 0.1);
}

void MMATMATMUL_Test(Party *proxy) {
    // In this function, we test several matrix multiplications of random matrices with the same size.
    cout << setfill('*') << setw(50) << "Calling MMATMATMUL";
    cout << setfill('*') << setw(49) << "*" << endl;

    // setting
    int n_matrices = 3;
    int a_row = WSZ*5;
    int a_col = WSZ*6;
    int b_col = WSZ*3;
    uint32_t *params = new uint32_t[1];
    uint32_t size = n_matrices * a_row * a_col * b_col; // size of the vector
    params[0] = size;
    double min_val = -100;
    double max_val = 100;

    uint64_t ***mat1 = new uint64_t **[n_matrices];
    uint64_t ***mat2 = new uint64_t **[n_matrices];
    for (int i = 0; i < n_matrices; i++) {
        mat1[i] = proxy->CreateShare(Random2dData(proxy, a_row, a_col, min_val, max_val), a_row, a_col);
        mat2[i] = proxy->CreateShare(Random2dData(proxy, a_col, b_col, min_val, max_val), a_col, b_col);
    }

    proxy->SendBytes(coreVectorisedMatrixMatrixMultiply, params, 1);
    uint64_t ***res = MatrixMatrixMultiply(proxy, mat1, mat2, n_matrices, a_row, a_col, b_col);

    double ***rec_mat1 = new double **[n_matrices];
    double ***rec_mat2 = new double **[n_matrices];
    for (int i = 0; i < n_matrices; i++) {
        rec_mat1[i] = ConvertToDouble(Reconstruct(proxy, mat1[i], a_row, a_col), a_row, a_col);
        rec_mat2[i] = ConvertToDouble(Reconstruct(proxy, mat2[i], a_col, b_col), a_col, b_col);
//        proxy->Print2dArray("mat1's reconstructed matrix " + to_string(i), rec_mat1[i], a_row, a_col);
//        proxy->Print2dArray("mat2's reconstructed matrix " + to_string(i), rec_mat2[i], a_col, b_col);
    }

//    for(int i = 0; i < n_matrices; i++) {
//        Print2dArray("Computed matrix multiplication - " + to_string(i), ConvertToDouble(
//                Reconstruct(proxy, res[i], a_row, b_col), a_row, b_col), a_row, b_col);
//        double** gt = MultiplyMatrices(rec_mat1[i], rec_mat2[i], a_row, a_col, b_col);
//        Print2dArray("GT matrix multiplication - " + to_string(i), gt, a_row, b_col);
//    }

    double *diffs = new double[n_matrices];
    bool flag = true;
    for (int m = 0; m < n_matrices; m++) {
        double **rec_res = ConvertToDouble(Reconstruct(proxy, res[m], a_row, b_col), a_row, b_col);
        double **gt = MultiplyMatrices(rec_mat1[m], rec_mat2[m], a_row, a_col, b_col);

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
    cout << setfill('*') << setw(50) << "Calling MatrixVectorMultiply";
    cout << setfill('*') << setw(49) << "*" << endl;
    cout << "calling this one" << endl;
    // setting
    int a_row = WSZ*3;
    int a_col = WSZ*WSZ;
    uint32_t *params = new uint32_t[1];
    uint32_t size = a_row * a_col; // size of the vector
    params[0] = size;
    double min_val = -0.784378;
    double max_val = 120.76;

    uint64_t **mat = proxy->CreateShare(Random2dData(proxy, a_row, a_col, min_val, max_val), a_row, a_col);
    uint64_t *vec = proxy->CreateShare(Random1dData(proxy, a_col, min_val, max_val), a_col);

    proxy->SendBytes(coreMatrixVectorMultiply, params, 1);
    uint64_t *res = MatrixVectorMultiply(proxy, mat, vec, a_row, a_col);
    double *rec_res = ConvertToDouble(Reconstruct(proxy, res, a_row), a_row);

    double **rec_mat = ConvertToDouble(Reconstruct(proxy, mat, a_row, a_col), a_row, a_col);
    double *rec_vec = ConvertToDouble(Reconstruct(proxy, vec, a_col), a_col);
    double *gt = MultiplyMatrixVector(rec_mat, rec_vec, a_row, a_col);

    double tmp = 0;
    for (int i = 0; i < a_row; i++) {
        tmp += abs(rec_res[i] - gt[i]);
    }

    if (tmp <= 0.1) {
        cout << "MatrixVectorMultiply works correctly" << endl;
        //cout << "Total absolute difference: " << tmp << endl;
    } else {
        cout << "MatrixVectorMultiply works incorrectly" << endl;
        cout << "Total absolute difference: " << tmp << endl;

        Print2dArray("mat", rec_mat, a_row, a_col);
        Print1dArray("vec", rec_vec, a_col);

        Print1dArray("Computed matrix-vector multiplication", ConvertToDouble(Reconstruct(proxy, res, a_row), a_row), a_row);
        Print1dArray("GT matrix-vector multiplication", gt, a_row);
    }
    return tmp <= 0.1;
}

void MMATVECMUL_Test(Party *proxy) {
    // In this function, we test the matrix multiplication of n_matrices matrices and vectors.
    cout << setfill('*') << setw(50) << "Calling MMATVECMUL";
    cout << setfill('*') << setw(49) << "*" << endl;

    // setting
    int n_matrices = 3;
    int a_row = WSZ*5;
    int a_col = WSZ*6;
    int b_col = WSZ*3;
    uint32_t *params = new uint32_t[1];
    uint32_t size = n_matrices * a_row * a_col; // size of the vector
    params[0] = size;
    double min_val = -0.5;
    double max_val = 99.0;

    uint64_t ***mat = new uint64_t **[n_matrices];
    uint64_t **vec = new uint64_t *[n_matrices];
    for (int i = 0; i < n_matrices; i++) {
        mat[i] = proxy->CreateShare(Random2dData(proxy, a_row, a_col, min_val, max_val), a_row, a_col);
        vec[i] = proxy->CreateShare(Random1dData(proxy, a_col, min_val, max_val), a_col);
    }

    proxy->SendBytes(coreVectorisedMatrixVectorMultiply, params, 1);
    uint64_t **res = MatrixVectorMultiply(proxy, mat, vec, n_matrices, a_row, a_col);

    double ***rec_mat = new double **[n_matrices];
    double **rec_vec = new double *[n_matrices];
    for (int i = 0; i < n_matrices; i++) {
        rec_mat[i] = ConvertToDouble(Reconstruct(proxy, mat[i], a_row, a_col), a_row, a_col);
        rec_vec[i] = ConvertToDouble(Reconstruct(proxy, vec[i], a_col), a_col);
    }

//    for(int i = 0; i < n_matrices; i++) {
//        Print1dArray("Computed matrix-vector multiplication - " + to_string(i), ConvertToDouble(Reconstruct(proxy, res[i], a_row), a_row), a_row);
//        double* gt = multiply_matrice_vector(rec_mat[i], rec_vec[i], a_row, a_col);
//        Print1dArray("Ground truth matrix-vector multiplication - " + to_string(i), gt, a_row);
//    }

    double *diffs = new double[n_matrices];
    bool flag = true;
    for (int m = 0; m < n_matrices; m++) {
        double *rec_res = ConvertToDouble(Reconstruct(proxy, res[m], a_row), a_row);
        double *gt = MultiplyMatrixVector(rec_mat[m], rec_vec[m], a_row, a_col);

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
//    cout<<setfill ('*')<<setw(50)<<"Calling InverseSqrt";
//    cout<<setfill ('*')<<setw(49)<<"*"<<endl;
//    // setting
//    int n_row = 4;
//    int n_col = 5;
//
//    // generate a Gram matrix
//    uint64_t **invsqrt_data = random_gram_matrix(proxy, n_row, n_col);
//
//    double** tmp = ConvertToDouble(Reconstruct(proxy, invsqrt_data, n_row, n_row), n_row, n_row);
//
//    proxy->SendBytes(RKN_INVSQRT, n_row);
//    uint64_t** invsqrt_G = InverseSqrt(proxy, invsqrt_data, n_row);
//
//    double** rec_invsqrt_G = ConvertToDouble(Reconstruct(proxy, invsqrt_G, n_row, n_row), n_row, n_row);
//
//    Print2dArray("The inverse square root of the Gram matrix", rec_invsqrt_G, n_row, n_row);
//
//    double* straighten_invsqrt_G = new double[n_row * n_row];
//    for(uint32_t i = 0; i < n_row * n_row; i++) {
//        straighten_invsqrt_G[i] = tmp[i % n_row][i / n_row];
//    }
//
//    Print2dArray("Gram matrix", tmp, n_row, n_row);
//    Print1dArray("Straighten Gram matrix", straighten_invsqrt_G, n_row * n_row);
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
//    Print2dArray("The inverse of the Gram matrix",
//                        MultiplyMatrices(rec_invsqrt_G, rec_invsqrt_G, n_row, n_row, n_row), n_row, n_row);
//
//    cout << "============= GT Inverse of the Gram matrix ======================" << endl;
//    cout << matrix_G.inverse() << endl;
//    cout << "============================================================================" << endl;
//
//}

bool MINVSQRT_Test(Party* proxy, double &exe_time, bool only_timing = false) {
    /* In this function, we test the computation of the inverse square root of a Gram matrix.
     * We first generate a random Gram matrix by first generating a random data matrix D and
     * then computing D^T * D.
     */
    cout<<setfill ('*')<<setw(50)<<"Calling MINVSQRT";
    cout<<setfill ('*')<<setw(49)<<"*"<<endl;
    // setting
    int n_row = 100;
    int n_col = 100;
    int n_gms = 10;

    // generate a Gram matrix
    uint64_t ***invsqrt_data = new uint64_t**[n_gms];
    for(int i = 0; i < n_gms; i++) {
        invsqrt_data[i] = RandomGramMatrix(proxy, n_row, n_col);
    }

//    double** tmp = ConvertToDouble(REC(proxy, invsqrt_data, n_row, n_row), n_row, n_row);

    uint32_t *params = new uint32_t[2];
    params[0] = n_gms;
    params[1] = n_row;

    proxy->SendBytes(rknVectorisedInverseSqrt, params, 2);
    auto start = chrono::high_resolution_clock::now();
    uint64_t*** invsqrt_G = InverseSqrt(proxy, invsqrt_data, n_gms, n_row);
    auto end = chrono::high_resolution_clock::now();
    double time_taken = chrono::duration_cast<chrono::nanoseconds>(end - start).count();
    time_taken *= 1e-9;
    cout << "Vectorized INVSQRT Time:\t" << fixed << time_taken << setprecision(9) << " sec" << endl;
    exe_time = time_taken;

    bool flag = true;
    if(!only_timing) {
        double ***Gs = new double**[n_gms];
        for(int i = 0; i < n_gms; i++) {
            Gs[i] = ConvertToDouble(Reconstruct(proxy, invsqrt_data[i], n_row, n_row), n_row, n_row);
        }

        double*** rec_invsqrt_G = new double**[n_gms];
        for(int g = 0; g < n_gms; g++) {
            rec_invsqrt_G[g] = ConvertToDouble(Reconstruct(proxy, invsqrt_G[g], n_row, n_row), n_row, n_row);
            Print2dArray("The inverse square root of the Gram matrix", rec_invsqrt_G[g], n_row, n_row);

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
            cout << eig_vecs * vals.cwiseSqrt().cwiseInverse().asDiagonal() * Eigen::Transpose(eig_vecs) << endl;
            cout << "============================================================================" << endl;

            Print2dArray("The inverse of the Gram matrix",
                         MultiplyMatrices(rec_invsqrt_G[g], rec_invsqrt_G[g], n_row, n_row, n_row), n_row, n_row);

            cout << "============= GT Inverse of the Gram matrix ======================" << endl;
            cout << matrix_G.inverse() << endl;
            cout << "============================================================================" << endl;
        }
    }
    return flag;
}

double fRand(double fMin, double fMax)
{
    double f = (double) rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

void DIV_Test(Party *proxy, int &cnt, bool verbose = false){
    if(verbose)
        cout << setfill('*') << setw(50) << "Calling Divide" << setfill('*') << setw(49) << "*" << endl;


//    int y_int = (rand() % 10 + 1) / 3.0;
//    int result = (rand() % 10 + 1) / 4.0;
//    uint64_t y = proxy->CreateShare(y_int);
//    uint64_t x = proxy->CreateShare(y_int * result + (rand() % 10));
    double x_d = fRand(-100, 100);
    double y_d = fRand(-100, 100);
    uint64_t y = proxy->CreateShare(y_d);
    uint64_t x = proxy->CreateShare(x_d);

    proxy->SendBytes(coreDivide);
    uint64_t div = Divide(proxy, x, y);
    uint64_t rec_dev = Reconstruct(proxy, div);
    double reconstructed_div = ConvertToDouble(rec_dev);

    // checking the result
    double originalX = ConvertToDouble(Reconstruct(proxy, x));
    double originalY = ConvertToDouble(Reconstruct(proxy, y));
    double computed_div = originalX / originalY;

    if(verbose) {
        cout << " =========================================== " << endl;
        cout << "X: " << originalX << " Y: " << originalY << endl;
    }
    if(abs(computed_div - reconstructed_div) < 0.0001){
        if(verbose) {
            cout<<"Divide works correctly"<<endl;
            cout<< "computed: " << reconstructed_div << " -- ground truth: " << computed_div << endl;
        }
    }
    else{
        cnt++;
        if(verbose) {
            cout<<"Divide works incorrectly" <<endl;
            cout<< "computed: " << reconstructed_div << " -- but should be: " << computed_div << endl;
        }
    }
    if(verbose) {
        cout << "Bitwise computed result: " << bitset<L_BIT>(rec_dev) << endl;
        cout << " =========================================== " << endl;
    }

}

void MDIV_Test(Party *proxy, int &cnt, bool verbose = false){
    if(verbose)
        cout << setfill ('*') << setw(50) << "Calling MDIV" << setfill ('*') << setw(49) << "*" << endl;

    double *x_d = new double[sz];
    double *y_d = new double[sz];
    uint64_t *x = new uint64_t[sz];
    uint64_t *y = new uint64_t[sz];
    for(int i = 0; i < sz; i++) {
        x_d[i] = fRand(-100, 100);
        y_d[i] = fRand(-100, 100);
        x[i] = proxy->CreateShare(x_d[i]);
        y[i] = proxy->CreateShare(y_d[i]);
    }

    uint32_t* params = new uint32_t[1];
    params[0] = sz;
    proxy->SendBytes(coreVectorisedDivide, params, 1);
    uint64_t *div = Divide(proxy, x, y, sz);
    uint64_t *rec_dev = Reconstruct(proxy, div, sz);
    double *reconstructed_div = ConvertToDouble(rec_dev, sz);

    // checking the result
    double *originalX = ConvertToDouble(Reconstruct(proxy, x, sz), sz);
    double *originalY = ConvertToDouble(Reconstruct(proxy, y, sz), sz);
    double *computed_div = new double[sz];
    for(int i = 0; i < sz; i++) {
        computed_div[i] = originalX[i] / originalY[i];
    }

    if(verbose)
        cout << " =========================================== " << endl;
    for(int i = 0; i < sz; i++) {
        if(abs(computed_div[i] - reconstructed_div[i]) < 0.0001){
//            cout << " --------------------------------------- " << endl;
//            cout<<"Divide works correctly"<<endl;
//            cout << "X: " << originalX[i] << " Y: " << originalY[i] << endl;
//            cout<< "computed: " << reconstructed_div[i] << " -- ground truth: " << computed_div[i] << endl;
//            cout << " --------------------------------------- " << endl;
        }
        else{
            cnt++;
            if(verbose) {
                cout << " --------------------------------------- " << endl;
                cout<<"Divide works incorrectly" <<endl;
                cout << "X: " << originalX[i] << " Y: " << originalY[i] << endl;
                cout<< "computed: " << reconstructed_div[i] << " -- but should be: " << computed_div[i] << endl;
                cout << " --------------------------------------- " << endl;
            }
        }
//        cout << "Bitwise computed result: " << bitset<L_BIT>(rec_dev) << endl;
    }
    if(verbose)
        cout << " =========================================== " << endl;

    delete [] x_d;
    delete [] y_d;
    delete [] x;
    delete [] y;
    delete [] div;
    delete [] rec_dev;
    delete [] reconstructed_div;
    delete [] originalX;
    delete [] originalY;
}

void NORM_Test(Party *proxy, int &cnt, bool verbose = false){
    if(verbose)
        cout << setfill ('*') << setw(50) << "Calling MNORM" << setfill ('*') << setw(49) << "*" << endl;

    double *x_d = new double[sz];
    double *y_d = new double[sz];
    uint64_t *x = new uint64_t[sz];
    uint64_t *y = new uint64_t[sz];
    for(int i = 0; i < sz; i++) {
        x_d[i] = fRand(0, 100);
        y_d[i] = fRand(100, 200);
        x[i] = proxy->CreateShare(x_d[i]);
        y[i] = proxy->CreateShare(y_d[i]);
    }

    uint32_t* params = new uint32_t[1];
    params[0] = sz;
    proxy->SendBytes(coreNormalise, params, 1);
    uint64_t *div = Normalise(proxy, x, y, sz);
    uint64_t *rec_dev = Reconstruct(proxy, div, sz);
    double *reconstructed_div = ConvertToDouble(rec_dev, sz);

    // checking the result
    double *originalX = ConvertToDouble(Reconstruct(proxy, x, sz), sz);
    double *originalY = ConvertToDouble(Reconstruct(proxy, y, sz), sz);
    double *computed_div = new double[sz];
    for(int i = 0; i < sz; i++) {
        computed_div[i] = originalX[i] / originalY[i];
    }

    if(verbose)
        cout << " =========================================== " << endl;
    for(int i = 0; i < sz; i++) {
        if(abs(computed_div[i] - reconstructed_div[i]) < 0.0001){
            cout << " --------------------------------------- " << endl;
            cout<<"MNORM works correctly"<<endl;
            cout << "X: " << originalX[i] << " Y: " << originalY[i] << endl;
            cout<< "computed: " << reconstructed_div[i] << " -- ground truth: " << computed_div[i] << endl;
            cout << " --------------------------------------- " << endl;
        }
        else{
            cnt++;
            if(verbose) {
                cout << " --------------------------------------- " << endl;
                cout<<"MNORM works incorrectly" <<endl;
                cout << "X: " << originalX[i] << " Y: " << originalY[i] << endl;
                cout<< "computed: " << reconstructed_div[i] << " -- but should be: " << computed_div[i] << endl;
                cout << " --------------------------------------- " << endl;
            }
        }
//        cout << "Bitwise computed result: " << bitset<L_BIT>(rec_dev) << endl;
    }
    if(verbose)
        cout << " =========================================== " << endl;

    delete [] x_d;
    delete [] y_d;
    delete [] x;
    delete [] y;
    delete [] div;
    delete [] rec_dev;
    delete [] reconstructed_div;
    delete [] originalX;
    delete [] originalY;

}

void LocalMultiply_Test(Party *proxy, int &cnt, bool verbose = false){
    if(verbose)
        cout << setfill('*') << setw(50) << "Calling LocalMultiply" << setfill('*') << setw(49) << "*" << endl;

    // secret data
    double x_d = fRand(-100, 100);
    uint64_t x = proxy->CreateShare(x_d);

    // public value
    double y_d = fRand(-100, 100);
    uint64_t y = ConvertToUint64(y_d);

    uint64_t res = LocalMultiply(x, y);
    uint64_t rec_res = Reconstruct(proxy, res);
    double reconstructed_res = ConvertToDouble(rec_res);

    // checking the result
    double originalX = ConvertToDouble(Reconstruct(proxy, x));
    double originalY = ConvertToDouble(y);
    double computed_mul = originalX * originalY;

    if(verbose) {
        cout << " =========================================== " << endl;
        cout << "X: " << originalX << " (" << x_d << ") " << " Y: " << originalY << " (" << y_d << ")" << endl;
    }
    if(abs(computed_mul - reconstructed_res) < 0.0001){
        if(verbose) {
            cout<<"LocalMultiply works correctly"<<endl;
            cout<< "computed: " << reconstructed_res << " -- ground truth: " << computed_mul << endl;
        }
    }
    else{
        cnt++;
        if(verbose) {
            cout<<"LocalMultiply works incorrectly" <<endl;
            cout<< "computed: " << reconstructed_res << " -- but should be: " << computed_mul << endl;
        }
    }
    if(verbose) {
        cout << "Bitwise computed result: " << bitset<L_BIT>(rec_res) << endl;
        cout << " =========================================== " << endl;
    }

}

void local_MMUL_Test(Party *proxy, int &cnt, bool verbose = false){
    if(verbose)
        cout << setfill ('*') << setw(50) << "Calling local_MMUL" << setfill ('*') << setw(49) << "*" << endl;

    // secret values
    double *x_d = new double[sz];
    uint64_t *x = new uint64_t[sz];

    // public values
    double *y_d = new double[sz];
    uint64_t *y = new uint64_t[sz];

    for(int i = 0; i < sz; i++) {
        x_d[i] = fRand(-100, 100);
        y_d[i] = fRand(-100, 100);
        x[i] = proxy->CreateShare(x_d[i]);
        y[i] = ConvertToUint64(y_d[i]);
    }

    uint64_t *res = LocalMultiply(x, y, sz);
    uint64_t *rec_res = Reconstruct(proxy, res, sz);
    double *reconstructed_res = ConvertToDouble(rec_res, sz);

    // checking the result
    uint64_t *rec_x = Reconstruct(proxy, x, sz);
    double *originalX = ConvertToDouble(rec_x, sz);
    double *originalY = ConvertToDouble(y, sz);
    double *computed_mul = new double[sz];
    for(int i = 0; i < sz; i++) {
        computed_mul[i] = originalX[i] * originalY[i];
    }

    if(verbose)
        cout << " =========================================== " << endl;
    for(int i = 0; i < sz; i++) {
        if(abs(computed_mul[i] - reconstructed_res[i]) < 0.0001){
            if(verbose) {
                cout << " --------------------------------------- " << endl;
                cout<<"local_MMUL works correctly"<<endl;
                cout << "X: " << originalX[i] << " Y: " << originalY[i] << endl;
                cout<< "computed: " << reconstructed_res[i] << " -- ground truth: " << computed_mul[i] << endl;
                cout << " --------------------------------------- " << endl;
            }
        }
        else{
            cnt++;
            if(verbose) {
                cout << " --------------------------------------- " << endl;
                cout<<"local_MMUL works incorrectly" <<endl;
                cout << "X: " << originalX[i] << " Y: " << originalY[i] << endl;
                cout<< "computed: " << reconstructed_res[i] << " -- but should be: " << computed_mul[i] << endl;
                cout << " --------------------------------------- " << endl;
            }
        }
    }
    if(verbose)
        cout << " =========================================== " << endl;

    delete [] x_d;
    delete [] x;
    delete [] y_d;
    delete [] y;
    delete [] computed_mul;
    delete [] res;
    delete [] rec_res;
    delete [] reconstructed_res;
    delete [] rec_x;
    delete [] originalX;
    delete [] originalY;
}

void ADD_Test(Party *proxy) {
    cout << setfill('*') << setw(50) << "Calling Add functions";
    cout << setfill('*') << setw(49) << "*" << endl;

    const int channel = 3;
    uint64_t ***x = new uint64_t **[channel];
    for (int i = 0; i < channel; ++i) {
        x[i] = proxy->CreateShare(Random2dData(proxy, 5, 5, -4999.0, 5000.0), 5, 5);
    }

    uint64_t **res = Add(proxy, x, channel, 5, 5);
    double **recon_res = ConvertToDouble(Reconstruct(proxy, res, 5, 5), 5, 5);

    // checking the result
    double ***originalX = new double **[channel];
    for (int i = 0; i < channel; ++i) {
        originalX[i] = ConvertToDouble(Reconstruct(proxy, x[i], 5, 5), 5, 5);
        delete[] x[i];
    }
    delete[] x;
    for (int r = 0; r < 5; ++r) {
        delete[] res[r];
    }
    delete[] res;

    double **correct_sum = new double *[5];
    bool allCorrect = true;
    for (uint64_t r = 0; r < 5; r++) {
        correct_sum[r] = new double[5];
        for (uint64_t c = 0; c < 5; c++) {
            double s = 0;
            for (int i = 0; i < channel; ++i) {
                s += originalX[i][r][c];
            }
            correct_sum[r][c] = s;
            if (abs(s - recon_res[r][c]) > 0.0001) {
                allCorrect = false;
            }
        }

    }
    if (allCorrect) {
        cout << "Add works correctly" << endl;
    } else {
        cout << "Add works incorrectly" << endl;
        Print2dArray("Computed SUM: ", recon_res, 5, 5);
        Print2dArray("Correct SUM: ", correct_sum, 5, 5);
    }
    for (int r = 0; r < 5; ++r) {
        delete[] correct_sum[r];
        delete[] recon_res[r];
    }
    delete[] correct_sum;
    delete[] recon_res;
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
//            anchor_points[i] = proxy->CreateShare(Random2dData(proxy, n_anc, n_dim, 1, false), n_anc, n_dim);
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
//        ct[i] = proxy->getPRole() * ((uint64_t) 1 << FRACTIONAL_BITS);
//        initial_ct[i] = ct[i];
//    }
////    proxy->Print1dArray("Initial ct", proxy->MConvertToDouble(proxy->MReconstruct(ct, size2 + n_anc), size2 + n_anc), size2 + n_anc);
//
//    for(int s = 0; s < length; s++) {
//        // generate a random data
//        cout << "s: " << s << " - Generate sample data..." << endl;
//        all_x[s] = proxy->CreateShare(Random1dData(proxy, n_dim, 1, false), n_dim);
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
//        Print1dArray("all_x[s]", ConvertToDouble(Reconstruct(proxy, all_x[s], n_dim), n_dim), n_dim);
//        Print1dArray("before ct", ConvertToDouble(Reconstruct(proxy, ct, size2), size2), size2);
//        uint64_t* tmp_ct = RknIteration(proxy, all_x[s], str_z, ct, n_dim, n_anc, k_mer, lambda, alpha);
//        copy(tmp_ct, tmp_ct + size2, ct + n_anc);
//        Print1dArray("after ct", ConvertToDouble(Reconstruct(proxy, ct, size2), size2), size2);
//
//        delete [] str_z;
//    }
//    cout << "Initial mapping is done!" << endl;
////    proxy->Print1dArray("c[t]", proxy->MConvertToDouble(proxy->MReconstruct(ct, size2 + n_anc), size2 + n_anc), size2 + n_anc);
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
//        rec_anc_points[i] = ConvertToDouble(Reconstruct(proxy, anchor_points[i], n_anc, n_dim), n_anc, n_dim);
//    }
//
////    cout << "GT: x, t1k1 and t1k..." << endl;
//    double** rec_all_x = ConvertToDouble(Reconstruct(proxy, all_x, length, n_dim), length, n_dim);
//    double* rec_ct = ConvertToDouble(Reconstruct(proxy, initial_ct, size2 + n_anc), size2 + n_anc);
//
//    for(int iter = 0; iter < length; iter++) {
//        Print1dArray("rec_all_x", rec_all_x[iter], n_dim);
//        Print1dArray("before rec_ct[t]", rec_ct, size2 + n_anc);
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
//        Print2dArray("exp_gt_dp " + to_string(iter), exp_gt_dp, k_mer, n_anc);
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
//        Print2dArray("gt_skt " + to_string(iter), gt_skt, k_mer, n_anc);
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
////        proxy->Print2dArray("gt_ckt", gt_ckt, k_mer, n_anc);
//
//        // update c[t] based on the result of the mappings in each k-mer
//        for(int i = 1; i < k_mer + 1; i++) {
//            for(int j = 0; j < n_anc; j++) {
//                rec_ct[i * n_anc + j] = gt_ckt[i - 1][j];
//            }
//        }
//        Print1dArray("after rec_ct[t]", rec_ct, size2 + n_anc);
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
//    Print1dArray("c[t]", ConvertToDouble(Reconstruct(proxy, ct, size2 + n_anc), size2 + n_anc), size2 + n_anc);
//
//    cout << "Ground truth: " << endl;
//    Print1dArray("GT c[t]", rec_ct, size2 + n_anc);
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
//            all_x[s] = proxy->CreateShare(Random1dData(proxy, n_dim, 1, false), n_dim);
//        }
//
//        // generate a random anchor points
//        cout << "Generating anchor points..." << endl;
//        for(int i = 0; i < k_mer; i++) {
//            anchor_points[i] = proxy->CreateShare(Random2dData(proxy, n_anc, n_dim, 1, false), n_anc, n_dim);
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
//        weights = proxy->CreateShare(Random1dData(proxy, n_anc + 1, 0.0, 1.0), n_anc + 1);
////        Print1dArray("Weights", ConvertToDouble(Reconstruct(proxy, weights, n_anc), n_anc), n_anc);
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
////        Print1dArray("Weights", ConvertToDouble(Reconstruct(proxy, weights, n_anc), n_anc), n_anc);
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
//        ct[i] = proxy->getPRole() * ((uint64_t) 1 << FRACTIONAL_BITS);
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
////        Print1dArray("all_x[s]", ConvertToDouble(Reconstruct(proxy, all_x[s], n_dim), n_dim), n_dim);
////        Print1dArray("before ct", ConvertToDouble(Reconstruct(proxy, ct, size2), size2), size2);
//        uint64_t* tmp_ct = RknIteration(proxy, all_x[s], str_z, ct, n_dim, n_anc, k_mer, lambda, alpha);
//        copy(tmp_ct, tmp_ct + size2, ct + n_anc);
////        Print1dArray("after ct", ConvertToDouble(Reconstruct(proxy, ct, size2), size2), size2);
//
//        uint64_t** mat_ct = new uint64_t *[k_mer];
//        for(int i = 0; i < k_mer; i++) {
//            mat_ct[i] = new uint64_t[n_anc];
//            for(int j = 0; j < n_anc; j++) {
//                mat_ct[i][j] = ct[n_anc + (i * n_anc) + j];
//            }
//        }
////        Print2dArray("Mappings at " + to_string(s), ConvertToDouble(Reconstruct(proxy, mat_ct, k_mer, n_anc), k_mer, n_anc),
////                     k_mer, n_anc, false);
//
//        delete [] str_z;
//    }
//
//    cout << "Initial mapping is done!" << endl;
////    proxy->Print1dArray("c[t]", proxy->MConvertToDouble(proxy->MReconstruct(ct, size2 + n_anc), size2 + n_anc), size2 + n_anc);
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
//    uint64_t*** gms = MatrixMatrixMultiply(proxy, anchor_points, tr_anchor_points, k_mer, n_anc, n_dim, n_anc);
//
//    proxy->SendBytes(RKN_GM2KM, k_mer, n_anc);
//    uint64_t*** kmer_kms = GaussianKernel(proxy, gms, convert2uint64(alpha), k_mer, n_anc);
//
//    double*** rec_kmer_kms = new double**[k_mer];
//    for(int g = 0; g < k_mer; g++) {
//        rec_kmer_kms[g] = ConvertToDouble(Reconstruct(proxy, kmer_kms[g], n_anc, n_anc), n_anc, n_anc);
//    }
//
//    proxy->SendBytes(RKN_MINVSQRT, k_mer, n_anc);
//    uint64_t*** invsqrt_gms = InverseSqrt(proxy, kmer_kms, k_mer, n_anc, epsilon);
//
//    proxy->SendBytes(CORE_MMATVECMUL, k_mer * n_anc * n_anc);
//    uint64_t** x_mapping = MatrixVectorMultiply(proxy, invsqrt_gms, mat_ct, k_mer, n_anc, n_anc);
//
//    double** rec_x_mapping = ConvertToDouble(Reconstruct(proxy, x_mapping, k_mer, n_anc), k_mer, n_anc);
//
//    proxy->SendBytes(CORE_DP, n_anc);
//    uint64_t prediction = DotProduct(proxy, weights, x_mapping[k_mer - 1], n_anc);
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
//        rec_anc_points[i] = ConvertToDouble(Reconstruct(proxy, anchor_points[i], n_anc, n_dim), n_anc, n_dim);
//    }
//
//    double** rec_all_x = ConvertToDouble(Reconstruct(proxy, all_x, length, n_dim), length, n_dim);
//    double* rec_ct = ConvertToDouble(Reconstruct(proxy, initial_ct, size2 + n_anc), size2 + n_anc);
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
//    //    proxy->Print2dArray("GT Gram matrix 0", gt_gms[0], n_anc, n_anc);
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
//    double* rec_weights = ConvertToDouble(Reconstruct(proxy, weights, n_anc), n_anc);
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
//    Print2dArray("Differences between mappings", diff, n_anc, k_mer, true);
//    Print1dArray("Total differences between mappings", total_diff, k_mer);
//
//    PrintValue("Prediction", ConvertToDouble(Reconstruct(proxy, prediction)));
//    PrintValue("GT Prediction", gt_prediction);
//    PrintValue("|Prediction - GT Prediction|", abs(ConvertToDouble(Reconstruct(proxy, prediction)) - gt_prediction));
//
////    MbubbleSort(gt_eigvals, k_mer, n_anc);
////    MbubbleSort(AT_gt_eigvals, k_mer, n_anc);
//
////    proxy->Print2dArray("GT Eigenvalues", gt_eigvals, k_mer, n_anc, false);
////    proxy->Print2dArray("AT GT Eigenvalues", AT_gt_eigvals, k_mer, n_anc, false);
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
        secret[0][i] = proxy->CreateShare(image[i], length);
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
            sec_kernel[0][0][i * 5 + j] = proxy->CreateShare(kernel[i][j]);
            sec_kernel[0][1][i * 5 + j] = proxy->CreateShare(kernel[i][(j + 1) % 5]);
        }
    }

    double bias = 0.1175425;
    uint64_t *sec_bias = new uint64_t[0];
    sec_bias[0] = proxy->CreateShare(bias);

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
    proxy->SendBytes(cnnConvolutionalLayer, params, 9);

    Print2dArray("image: ", ConvertToDouble(Reconstruct(proxy, secret[0], length, length), length, length), length, length);
    uint64_t ***weights_secure = ConvolutionalLayer(proxy, secret, 1, length, length, sec_kernel, 5, 2, 1, 2, 2,
                                                    sec_bias, false);
    // checking the result
    double **recon_res = ConvertToDouble(Reconstruct(proxy, weights_secure[0], max_height, max_width), max_height, max_width);

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
            // Activation: Relu
            double relu = 0.0;
            if (dot_product > 0) {
                relu = dot_product;
            }
            correct_res[cr][cc] = relu;
        }
    }
    //Print2dArray("activated conv", correct_res, conv_size, conv_size);
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
            // Activation: Relu
            double relu = 0.0;
            if (dot_product > 0) {
                relu = dot_product;
            }
            correct_res[cr][cc] = relu;
        }
    }
    //Print2dArray("activated conv", correct_res, conv_size, conv_size);
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
        Print2dArray("Computed first conv: ", recon_res, max_height, max_width);
        Print2dArray("Correct first conv: ", pooled_conv, max_height, max_width);
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
            secret[i][j] = proxy->CreateShare(input[i][j], cols);
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
            sec_kernel[i][j] = proxy->CreateShare(kernel[i][j], k_size*k_size);
        }
    }
    double bias[n_kernel] = {0.006893986,0.028841767,-0.011976943}; //,0.034752205,0.026028743,0.045701645,-0.014377818,-0.0099255275,-0.02820568,0.038818415,0.017578261,-0.037074637,-0.033514146,0.006196057,-0.013816738,-0.03915469,0.039954044,0.05106463,-0.032943178,-0.023218332,-0.01375858,0.011314128,-0.0073599,-0.012863117,0.0317652,-0.022406086,0.0480716,0.04324196,-0.018559806,0.044528794,-0.038379807,-0.035707023,-0.011929117,-0.03329016,0.03991615,0.013486342,-0.035355117,0.010999725,0.019694611,-0.0106660845,0.021526804,0.0036931,-0.015429371,-0.00035950931,-0.026158424,-0.023419743,-0.03008074,0.052684106,0.02619821,-0.012335571};
    uint64_t * sec_bias = proxy->CreateShare(bias, n_kernel);

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
    proxy->SendBytes(cnnConvolutionalLayer, params, 9);
    cout << "call ConvolutionalLayer" << endl;
    uint64_t ***conv = ConvolutionalLayer(proxy, secret, channels, rows, cols, sec_kernel, k_size, n_kernel, 1, 2, 2,
                                          sec_bias, false);
    cout << "returned from ConvolutionalLayer" << endl;
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
                    //Relu
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
            rec_conv[0][0] = ConvertToDouble(Reconstruct(proxy, conv[0][0], n_kernel * conv_row * conv_col), n_kernel * conv_row * conv_col);
        }
        for (int k = 0; k < n_kernel; ++k) {
            if(!params[8]){
                rec_conv[k] = ConvertToDouble(Reconstruct(proxy, conv[k], conv_row, conv_col), conv_row, conv_col);
            }
            //Print2dArray("computed conv kernel 1", rec_conv[k], conv_row, conv_col);
            //Print2dArray("correct conv output", mul[k], mul_row, mul_col);
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
            Print2dArray("secret", ConvertToDouble(Reconstruct(proxy, secret[i], rows, cols), rows, cols), rows, cols);
            Print2dArray("kernels", ConvertToDouble(Reconstruct(proxy, sec_kernel[i], n_kernel, k_size * k_size), n_kernel, k_size * k_size), n_kernel, k_size * k_size);

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
                Print2dArray("computed conv: ", rec_conv[o], conv_row, conv_col);
                delete[] rec_conv[o];
            }
            Print2dArray("actual result: ", corr_conv[o], conv_row, conv_col);
            delete[] corr_conv[o];
        }
        if(params[8]){
            Print1dArray("computed conv: ", rec_conv[0][0], n_kernel*conv_row*conv_col);
            delete[] rec_conv[0][0];
            delete[] rec_conv[0];
            delete[] rec_conv;
        }
    }
    return correct;

}


// Main function to run the experiments
int main(int argc, char *argv[]) {
    // MAPS FUNCTIONS TO NUMBERS SO THAT WE CAN USE SWITCH-CASE TO DETERMINE WHICH FUNCTION WILL BE RUN
    unordered_map<int, string> function_mapping;
    function_mapping[0] = "All (default)";
    function_mapping[1] = "Vectorized MUL";
    function_mapping[2] = "Vectorized MUX";
    function_mapping[3] = "Vectorized DP";
    function_mapping[4] = "Vectorized MOC";
    function_mapping[5] = "Vectorized MSB";
    function_mapping[6] = "Vectorized CMP";
    function_mapping[7] = "Vectorized EXP";
    function_mapping[8] = "Vectorized INVSQRT";

    // DEFINES THE FLAGS AND THE OPTIONS
    unordered_map<string, string> flag_mapping;
    flag_mapping["role"] = "-r";
    flag_mapping["proxy-port"] = "-cp";
    flag_mapping["proxy-ip"] = "-ci";
    flag_mapping["helper-port"] = "-hp";
    flag_mapping["helper-ip"] = "-hi";
    flag_mapping["function"] = "-f";
    flag_mapping["save-exe-time-result"] = "-s";
    flag_mapping["only-timing"] = "-t";
    flag_mapping["exp-id"] = "-id";
    flag_mapping["network"] = "-n";

    // CHECK IF AT LEAST THE ROLE IS GIVEN
    if(argc == 1) {
        cout << "Too less argument!" << endl;
        cout << "Run \"proxy_test -h\" to see the required arguments and other options" << endl;
        return -1;
    }

    // CHECK IF THERE ARE TOO MANY ARGUMENTS
    if(argc > (flag_mapping.size() + 1)) {
        cout << "Too many arguments!" << endl;
        cout << "Run \"proxy_test -h\" to see the required arguments and other options" << endl;
        return -14;
    }

    // HELP PAGE PRINT OUTS
    if(argc == 2 && strcmp(argv[1], "-h") == 0) {
        cout << "proxy_test [flag] [flag-value]" << endl;
        cout << flag_mapping["role"] << "\t\t" << "(REQUIRED) 0 or 1 determining the role, i.e. whether the proxy acts as P1 or P2, respectively" << endl;
        cout << flag_mapping["proxy-port"] << "\t\t" << "Port for the connection between proxies - default 9999" << endl;
        cout << flag_mapping["proxy-ip"] << "\t\t" << "IP address for proxy connection - default 127.0.0.1" << endl;
        cout << flag_mapping["helper-port"] << "\t\t" << "Port for the connection to the helper - default 7777" << endl;
        cout << flag_mapping["helper-ip"] << "\t\t" << "IP address for helper connection - default 127.0.0.1" << endl;
        cout << flag_mapping["function"] << "\t\t" << "Function: determines which functions are going to be tested" << endl;
        for(int i = 0; i < function_mapping.size(); i++) {
            cout << "\t\t\t" << i << ": " << function_mapping[i] << endl;
        }
        cout << flag_mapping["only-timing"] << "\t\t" << "Indicate whether the ground truths are going to be computed as well - default 0" << endl;
        cout << "\t\t\t" << "0: Compute ground truths and check correctness" << endl;
        cout << "\t\t\t" << "1: Testing only timing and skip ground truth computations" << endl;
        cout << flag_mapping["save-exe-time-result"] << "\t\t" << "Determines if the result will be saved - default \"False\"" << endl;
        cout << "\t\t\t" << "0: False - Do not save results" << endl;
        cout << "\t\t\t" << "1: True - Save the results" << endl;
        cout << flag_mapping["exp-id"] << "\t\t" << "ID of the experiment that is going to be used in the name of the file in case the output wants to be saved - default \"DEFAULT\"" << endl;
        cout << flag_mapping["network"] << "\t\t" << "Network type (LAN or WAN) that is going to be used in the name of the file in case the output wants to be saved - default \"NA\"" << endl;
        return 0;
    }

    uint8_t role;
    uint16_t cport = 9999;
    string caddress = "127.0.0.1";
    uint16_t hport = 7777;
    string haddress = "127.0.0.1";
    int function = 0;
    bool timing_exp = false;
    bool save_flag = false;
    string exp_id = "DEFAULT";
    string network = "NA";

    bool role_set = false;

    // CHECK THE GIVEN ARGUMENTS AND ASSIGN THEM ACCORDINGLY IF THEY ARE CORRECT
    for(int i = 1; i < argc; i += 2) {
        char* option = argv[i];
        // CHECK IF THERE IS A NEXT ARGUMENT
        if(argc < i + 2) {
            cout << "Argument is not given for the flag " << option << "!" << endl;
            return -15;
        }

        if(option == flag_mapping["role"]) {
            if(isNumberCheck(argv[i + 1]) && (atoi(argv[i + 1]) == 0 || atoi(argv[i + 1]) == 1)) {
                role = atoi(argv[i + 1]);
                role_set = true;
            }
            else {
                cout << "Role needs to be 0 or 1!" << endl;
                return -2;
            }
        }
        else if(option == flag_mapping["proxy-port"]) {
            if(isNumberCheck(argv[i + 1]))
                cport = atoi(argv[i + 1]);
            else {
                cout << "Proxy port needs to be number!" << endl;
                return -3;
            }
        }
        else if(option == flag_mapping["proxy-ip"]) {
            if(validateIP(argv[i + 1]))
                caddress = argv[i + 1];
            else {
                cout << "Proxy IP is not valid!" << endl;
                return -4;
            }
        }
        else if(option == flag_mapping["helper-port"]) {
            if(isNumberCheck(argv[i + 1]))
                hport = atoi(argv[i + 1]);
            else {
                cout << "Helper port needs to be number!" << endl;
                return -5;
            }
        }
        else if(option == flag_mapping["helper-ip"]) {
            if(validateIP(argv[i + 1]))
                haddress = argv[i + 1];
            else {
                cout << "Helper IP is not valid!" << endl;
                return -6;
            }
        }
        else if(option == flag_mapping["function"]) {
            if(isNumberCheck(argv[i + 1]) && atoi(argv[i + 1]) < function_mapping.size() && atoi(argv[i + 1]) >= 0)
                function = atoi(argv[i + 1]);
            else {
                cout << "Function is not valid!" << endl;
                return -7;
            }
        }
        else if(option == flag_mapping["only-timing"]) {
            if(isNumberCheck(argv[i + 1])) {
                if(atoi(argv[i + 1]) == 1) {
                    timing_exp = true;
                }
            }
            else {
                cout << flag_mapping["only-timing"] << " needs to be an integer!" << endl;
                return -8;
            }
        }
        else if(option == flag_mapping["save-exe-time-result"]) {
            if(isNumberCheck(argv[i + 1])) {
                if(atoi(argv[i + 1]) == 1) {
                    save_flag = true;
                }
            }
            else {
                cout << flag_mapping["save-exe-time-result"] << " needs to be an integer!" << endl;
                return -9;
            }
        }
        else if(option == flag_mapping["exp-id"]) {
            exp_id = argv[i + 1];
        }
        else if(option == flag_mapping["network"]) {
            if(!isNumberCheck(argv[i + 1])) {
                for (int c = 0; c < strlen(argv[i + 1]); c++)
                    argv[i + 1][c] = toupper(argv[i + 1][c]);
                if(strcmp(argv[i + 1], "LAN") == 0 || strcmp(argv[i + 1], "WAN") == 0) {
                    network = argv[i + 1];
                }
                else {
                    cout << flag_mapping["network"] << " needs to be either LAN or WAN!" << endl;
                    return -11;
                }
            }
            else {
                cout << flag_mapping["network"] << " needs to be string and either LAN or WAN!" << endl;
                return -12;
            }
        }
        else {
            cout << "There is no flag with name " << option << "!" << endl;
            return -13;
        }
    }

    if(!role_set) {
        cout << "Role must be set!" << endl;
        return -16;
    }




//    uint8_t role = atoi(argv[1]);
//    uint16_t cport = atoi(argv[2]);
//    string caddress(argv[3]);
//    uint16_t hport = atoi(argv[4]);
//    string haddress(argv[5]);

    // function to run
//    int function = 0;
//    if(argc >= 7) {
//        function = atoi(argv[6]);
//    }
//    cout << "function: " << function << endl;

    // true if the experiments are only for timing analysis -- skipping ground truths of the functions
//    bool timing_exp = false;
//    if(argc >= 8 && atoi(argv[7]) == 1) {
//        timing_exp = true;
//    }
//    cout << "timing experiment: " << timing_exp << endl;

    // experiment ID
//    string exp_id = "default";
//    if(argc >= 9) {
//        exp_id = argv[8];
//    }
//    cout << "experiment id: " << exp_id << endl;

    // experiment ID
//    string network = "not-specified";
//    if(argc >= 10) {
//        network = argv[9];
//    }
//    cout << "network: " << network << endl;

//    if(argc > 10) {
//        cout << "Too many arguments!" << endl;
//        return -1;
//    }

    Party *proxy;
    if (role == 0)
        proxy = new Party(proxy1, hport, haddress, cport, caddress);
    else
        proxy = new Party(proxy2, hport, haddress, cport, caddress);

    srand((unsigned) time(NULL));

    bool result = true;
    int ind = 0;
    int cnt = 0;
    unordered_map<string, int> umap;
    double *exe_time = new double[num_repetition];
    double tmp_time;
    auto start = chrono::high_resolution_clock::now();
    exe_time[ind] = 0;
    while (ind < num_repetition ) { // && result
        cout << "Iteration: " << ind << endl;
        switch(function) {
            case 1:
                MMUL_Test(proxy, exe_time[ind], timing_exp);
                break;
            case 2:
                MMUX_Test(proxy, exe_time[ind], timing_exp);
                break;
            case 3:
                MDP_Test(proxy, exe_time[ind], timing_exp);
                break;
            case 4:
                MMOC_Test(proxy, exe_time[ind], timing_exp);
                break;
            case 5:
                MMSB_Test(proxy, exe_time[ind], timing_exp);
                break;
            case 6:
                MCMP_Test(proxy, exe_time[ind], timing_exp);
                break;
            case 7:
                MEXP_Test(proxy, exe_time[ind], cnt, timing_exp);
                break;
            case 8:
                MINVSQRT_Test(proxy, exe_time[ind], timing_exp);
                break;
            default: // all functions
                MMUL_Test(proxy, exe_time[ind], timing_exp);
                MMOC_Test(proxy, exe_time[ind], timing_exp);
                MMSB_Test(proxy, exe_time[ind], timing_exp);
                MCMP_Test(proxy, exe_time[ind], timing_exp);
                MMUX_Test(proxy, exe_time[ind], timing_exp);
                MDP_Test(proxy, exe_time[ind], timing_exp);
                MEXP_Test(proxy, exe_time[ind], cnt, timing_exp);
                MINVSQRT_Test(proxy, exe_time[ind], timing_exp);
        }

        // **************************** test cases for Ali ************************************
//        result = MUL_Test_v2(proxy, ind, umap, cnt);
        // ************************************************************************************

//        LocalMultiply_Test(proxy, cnt);
//        local_MMUL_Test(proxy, cnt);

//        result = TRUNCATE_Test(proxy, cnt, umap);
//
//        ADD_Test(proxy);
//
//        MUL_Test(proxy);
//        MMUL_Test(proxy, exe_time[ind], timing_exp);
//
//        MAX_Specific_Test(proxy);
//
//        MOC_Test(proxy);
//        MMOC_Test(proxy, exe_time[ind], timing_exp);
//
//        MSB_Test(proxy);
//        MMSB_Test(proxy, exe_time[ind], timing_exp);
//
//        CMP_Test(proxy);
//        MCMP_Test(proxy, exe_time[ind], timing_exp);
//
//        MUX_Test(proxy, cnt);
//        MMUX_Test(proxy, exe_time[ind], timing_exp);
//
//        MAX_Test(proxy);
//        MMAX_Test(proxy);
//
//        RST_Test(proxy);
//        RELU_Test(proxy);
//        MRELU_Test(proxy);
//
//        DRLU_Test(proxy);
//        ARGMAX_Test(proxy);
//        MDRLU_Test(proxy); //TODO
//
//        DIV_Test(proxy, cnt);
//        MDIV_Test(proxy, cnt);
//
//        NORM_Test(proxy, cnt, true);
//
//        INC_Test(proxy);
//        FLT_Test(proxy);
//        FCL_Test(proxy);
//        PAD_Test(proxy);
//        CL_Test(proxy);
//
//        EXP_Test(proxy);
//        MEXP_Test(proxy, exe_time[ind], timing_exp);
//
//        DP_Test(proxy);
//        MDP_Test(proxy, exe_time[ind], timing_exp);
//        MATMATMUL_Test(proxy);
//        MMATMATMUL_Test(proxy);
//
//        MATVECMUL_Test(proxy);
//
//        INVSQRT_Test(proxy);
//        MINVSQRT_Test(proxy, exe_time[ind], timing_exp);
//
//        ppRKN_ITER_Test(proxy);
//        ppRKN_PREDICTION_Test(proxy);
        ind++;
    }
    cout << "Number of wrong operations: " << cnt << endl;

    double total_exe = 0;
    for(int i = 0; i < num_repetition; i++) {
        cout << exe_time[i] << endl;
        total_exe += exe_time[i];
    }
    double average_runtime = total_exe / num_repetition;
    cout << "Average execution time: " << average_runtime << endl;

    // write the average runtime to a file if it is P1
//    if(proxy->getPRole() == P1) {
//        ofstream output_file;
//        output_file.open ("/Users/aliburak/Projects/CECILIA/exp_runners/rkn_experiments/benchmarking/function_" +
//            to_string(function) + "_" + exp_id + "_" + network + ".txt");
//        output_file << fixed << setprecision(3) << average_runtime;
//        output_file.close();
//    }



//    bool all_correct0 = true;
//    int counter = 0;
//    while (all_correct0 and counter < 1000) {
//        all_correct0 = FCL_Test(proxy); //NETWORK_M_INPUTS_TEST(proxy);
//        counter++;
//        cout << counter << endl;
//    }



    proxy->SendBytes(coreEnd);
    PrintBytes();

    auto end = chrono::high_resolution_clock::now();
//    double time_taken =
//            chrono::duration_cast<chrono::nanoseconds>(end - start).count();
//    cout<<"Total Time\t"<<fixed
//        << time_taken << setprecision(9) << " sec" << endl;

    delete proxy;

    return 0;
}

