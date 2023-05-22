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

constexpr int sz = 5;
constexpr int ringbits = 23;

void SORT_test(Party *proxy){
    ofstream txt;
    cout<<setfill ('*')<<setw(50)<<"Calling SORT_test";
    cout<<setfill ('*')<<setw(49)<<"*"<<endl;
    uint32_t size = sz;
    cout<<"Size: "<<size<<endl;
    auto* a =new uint64_t[size];
    for(int i = 0; i<size; i++) {
        //a[i] = proxy->generateCommonRandom();
        a[i] = size-i-1;
        //cout << a[i] << " ";
    }

    cout << "Creating shares...\n";
    auto* x = new uint64_t[size];
    for (int i=0;i<size;i++) {
        x[i] = proxy->createShare(a[i]);
    }

    cout << "Calling SendBytes..\n";
    uint32_t params[2];
    params[0] = size;
    params[1] = ringbits;
    //proxy->SendBytes(CORE_SORT, params, 1);
    proxy->SendBytes(CORE_SORT2, params, 2);
    cout << "Calling SORT..\n";
    auto start = chrono::high_resolution_clock::now();
    uint64_t* s = SORT(proxy, x, size,ringbits);
    auto end = chrono::high_resolution_clock::now();
    double totaltime =
            chrono::duration_cast<chrono::nanoseconds>(end - start).count()*1e-9;
    cout<<totaltime<<endl;

   cout << "Callng REC..\n";
    uint64_t* sorted = RECN(proxy,s,size, ringbits);

    for(int i = 0;i<size;i++){
        cout <<  sorted[i]<< endl;
    }

    cout<<"Array successfully sorted"<<endl;

   delete []a;
   delete []x;


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
    SORT_test(proxy);
    auto end = chrono::high_resolution_clock::now();
    double time_taken =
            chrono::duration_cast<chrono::nanoseconds>(end - start).count();
    time_taken *= 1e-9;
    cout<<time_taken<<endl;


    proxy->SendBytes(CORE_END);
    //proxy->PrintBytes();

    return 0;
}
