//
// Created by Seyma Selcan on 02.12.22.
//
#include <cstdlib>
#include <iostream>
#include <chrono>
#include <tuple>
#include <iomanip>
#include "../../core/core.h"
#include "../../core/sort.h"

using namespace std;

constexpr int sz = 100;
constexpr int cols = 2;
constexpr int ringbits = 20;

void VSORT_test(Party *proxy){
    ofstream txt;
    cout<<setfill ('*')<<setw(50)<<"Calling SORT_test";
    cout<<setfill ('*')<<setw(49)<<"*"<<endl;
    uint32_t size = sz;
    cout<<"Size: "<<size<<endl;
    auto** a =new uint64_t*[cols];
    for(int i = 0; i<cols; i++) {
        a[i] =new uint64_t[sz];
        for (int j = 0; j < sz; ++j) {
            if (proxy->getPRole() == P1){
                a[i][j] = size-j-1;
            }else{
                a[i][j] = 0;
            }
        }

    }

    uint32_t params[2];
    params[0] = size;
    params[1] = cols;
    proxy->SendBytes(CORE_VSORT, params, 2);
    cout << "Calling VSORT..\n";
    auto start = chrono::high_resolution_clock::now();
    uint64_t** s = Sort(proxy, a, size, cols, 0);
    auto end = chrono::high_resolution_clock::now();
    double totaltime =
            chrono::duration_cast<chrono::nanoseconds>(end - start).count()*1e-9;
    cout<<totaltime<<endl;

    cout << "Callng Reconstruct..\n";
    uint64_t** sorted = Reconstruct(proxy, s, cols, size);

    for(int i = 0;i<10;i++){
        cout <<  sorted[0][i] << sorted[1][i]<< endl;
    }

    cout<<"Array successfully sorted"<<endl;

    delete []a;
}

void SORT_test(Party *proxy){
    ofstream txt;
    cout<<setfill ('*')<<setw(50)<<"Calling SORT_test";
    cout<<setfill ('*')<<setw(49)<<"*"<<endl;
    uint32_t size = sz;
    cout<<"Size: "<<size<<endl;
    auto* a =new uint64_t[size];
    for(int i = 0; i<size; i++) {
        //a[i] = proxy->generateCommonRandom();
        if (proxy->getPRole() == P1){
            a[i] = size-i-1;
        }else{
            a[i] = 0;
        }
    }


    uint32_t params[2];
    params[0] = size;
    params[1] = ringbits;
    proxy->SendBytes(CORE_SORT, params, 1);
    //proxy->SendBytes(CORE_SORT2, params, 2);
    cout << "Calling Sort..\n";
    auto start = chrono::high_resolution_clock::now();
    uint64_t* s = Sort(proxy, a, size);
    //uint64_t* s = SortNarrow(proxy, a, size,ringbits);
    auto end = chrono::high_resolution_clock::now();
    double totaltime =
            chrono::duration_cast<chrono::nanoseconds>(end - start).count()*1e-9;
    cout<<totaltime<<endl;

    cout<<"A2X Time:\t"<<a2x_time<<endl;
    cout<<"GRP Time:\t"<<grp_time<<endl;
    cout<<"APP Time:\t"<<app_time<<endl;
    cout<<"X2A Time:\t"<<x2a_time<<endl;
    cout<<"GEP Time:\t"<<gp_time<<endl;
    cout<<"Multiply Time:\t"<<mul_time<<endl;
    cout<<"Reconstruct Time:\t"<<recn_time<<endl;
    cout<<"CCM Time:\t"<<comp_time<<endl;
    cout<<"TOT Time:\t"<<t_time<<endl;

    cout << "Callng Reconstruct..\n";
    //uint64_t* sorted = ReconstructNarrow(proxy,s,size, ringbits);
    uint64_t* sorted = Reconstruct(proxy, s, size);

    for(int i = 0;i<20;i++){
        cout <<  sorted[i]<< endl;
    }

    cout<<"Array successfully sorted"<<endl;

    delete []a;
}

int main(int argc, char* argv[]) {
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



    auto start = chrono::high_resolution_clock::now();
    VSORT_test(proxy);
    //SORT_test(proxy);
    auto end = chrono::high_resolution_clock::now();
    double time_taken =
            chrono::duration_cast<chrono::nanoseconds>(end - start).count();
    time_taken *= 1e-9;
    cout<<time_taken<<endl;


    proxy->SendBytes(CORE_END);
    //proxy->PrintBytes();

    return 0;
}
