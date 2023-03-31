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

constexpr int sz = 1000000;

void SORT_test(Party *proxy){
    ofstream txt;
    cout<<setfill ('*')<<setw(50)<<"Calling SORT_test";
    cout<<setfill ('*')<<setw(49)<<"*"<<endl;
    uint32_t size = sz;
    cout<<"Size: "<<size<<endl;
    auto* a =new uint64_t[size];
    for(int i = 0; i < size; i++) {
        a[i] = proxy->generateCommonRandom()*proxy->generateCommonRandom();
        //cout << a[i] << " ";
    }

    cout << "Creating shares...\n";
    auto* x = new uint64_t[size];
    for (int i=0;i<size;i++) {
        x[i] = proxy->createShare(a[i]);
    }

    cout << "Calling SendBytes..\n";
    uint32_t params[1];
    params[0] = size;
    proxy->SendBytes(CORE_SORT, params, 1);
    cout << "Calling SORT..\n";
    auto start = chrono::high_resolution_clock::now();
    uint64_t* s = SORT(proxy, x, size);
    auto end = chrono::high_resolution_clock::now();
    double totaltime =
            chrono::duration_cast<chrono::nanoseconds>(end - start).count()*1e-9;
    cout<<totaltime<<endl;
    txt.open("sorting-timings.txt", ios_base::app);
    //if (txt.is_open()){
    txt << sz << "\t" << totaltime << "\n";
    txt.close();
    //}
    //else cout << "Unable to open file";
    cout << "Callng REC..\n";
    uint64_t* sorted = REC(proxy,s,size);

    for(int i = 1;i<size;i++){
        if(sorted[i]<=sorted[i-1]){
            cout<<"Sort failed"<<endl;
            break;
        }
    }
    cout<<"Array successfully sorted"<<endl;


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
