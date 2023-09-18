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
constexpr int cols = 2;
constexpr int ringbits = 20;

void WriteResults(Party *proxy){

    string fileLine;
    vector<string> lines;
    string fileName = "/Users/sesame/CECILIA_Refactored/results1.txt";
    if(proxy->GetPRole()==proxy2) fileName = "/Users/sesame/CECILIA_Refactored/results2.txt";
    ifstream readFile(fileName);
    if (readFile.is_open()){
        while (getline(readFile, fileLine)){
            lines.push_back(fileLine);
        }
    }
    readFile.close();
    ofstream newFile(fileName, ofstream::out);

    newFile <<lines[0] << " " << a2x_time << endl;
    newFile <<lines[1] << " " << grp_time << endl;
    newFile <<lines[2] << " " << app_time << endl;
    newFile <<lines[3] << " " << x2a_time << endl;
    newFile <<lines[4] << " " << gp_time << endl;
    newFile <<lines[5] << " " << mul_time << endl;
    newFile <<lines[6] << " " << recn_time << endl;
    newFile <<lines[7] << " " << comp_time << endl;
    newFile <<lines[8] << " " << t_time << endl;

    newFile.close();
}

void VSortTest(Party *proxy){
    ofstream txt;
    cout<<setfill ('*')<<setw(50)<<"Calling SORT_test";
    cout<<setfill ('*')<<setw(49)<<"*"<<endl;
    uint32_t size = sz;
    cout<<"Size: "<<size<<endl;
    auto** a =new uint64_t*[cols];
    for(int i = 0; i<cols; i++) {
        a[i] =new uint64_t[sz];
        for (int j = 0; j < sz; ++j) {
            if (proxy->GetPRole() == proxy1){
                a[i][j] = size-j-1;
            }else{
                a[i][j] = 0;
            }
        }

    }

    uint32_t params[2];
    params[0] = size;
    params[1] = cols;
    proxy->SendBytes(coreVSort, params, 2);
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
void SortNarrowTest(Party *proxy){
    ofstream txt;
    cout<<setfill ('*')<<setw(50)<<"Calling SORT_test";
    cout<<setfill ('*')<<setw(49)<<"*"<<endl;
    uint32_t size = sz;
    cout<<"Size: "<<size<<endl;
    auto* a =new uint64_t[size];
    for(int i = 0; i<size; i++) {
        //a[i] = proxy->generateCommonRandom();
        if (proxy->GetPRole() == proxy1){
            a[i] = size-i-1;
        }else{
            a[i] = 0;
        }
    }


    uint32_t params[2];
    params[0] = size;
    params[1] = ringbits;
    proxy->SendBytes(coreSort2, params, 2);
    cout << "Calling Sort..\n";
    auto start = chrono::high_resolution_clock::now();
    uint64_t* s = SortNarrow(proxy, a, size,ringbits);
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

    WriteResults(proxy);

    cout << "Callng Reconstruct..\n";
    uint64_t* sorted = ReconstructNarrowPermutation(proxy,s,size, ringbits);

    for(int i = 0;i<10;i++){
        cout <<  sorted[i]<< endl;
    }

    cout<<"Array successfully sorted"<<endl;

    delete []a;
    delete []sorted;
    delete []s;}

void SortTest(Party *proxy){
    ofstream txt;
    cout<<setfill ('*')<<setw(50)<<"Calling SORT_test";
    cout<<setfill ('*')<<setw(49)<<"*"<<endl;
    uint32_t size = sz;
    cout<<"Size: "<<size<<endl;
    auto* a =new uint64_t[size];
    for(int i = 0; i<size; i++) {
        //a[i] = proxy->generateCommonRandom();
        if (proxy->GetPRole() == proxy1){
            a[i] = ConvertToUint64(double(size-i-1));
        }else{
            a[i] = ConvertToUint64(0.0);
        }
    }


    uint32_t params[1];
    params[0] = size;
    proxy->SendBytes(coreSort, params, 1);
    cout << "Calling Regular Sort..\n";
    auto start = chrono::high_resolution_clock::now();
    uint64_t* s = Sort(proxy, a, size);
    auto end = chrono::high_resolution_clock::now();
    double totaltime =
            chrono::duration_cast<chrono::nanoseconds>(end - start).count()*1e-9;
    cout<<totaltime<<endl;

    cout << "Callng Reconstruct..\n";
    uint64_t* sorted = Reconstruct(proxy, s, size);
    auto sorted2 = ConvertToDouble(sorted, size, FRACTIONAL_BITS);

    for(int i = 0;i<5;i++){
        cout <<  sorted2[i]<< endl;
    }

    cout<<"Array successfully sorted"<<endl;

    delete []a;
    delete []sorted;
    delete []s;
}

int main(int argc, char* argv[]) {
    uint8_t role = atoi(argv[1]);
    uint16_t cport = atoi(argv[2]);
    string caddress(argv[3]);
    uint16_t hport = atoi(argv[4]);
    string haddress(argv[5]);


    Party *proxy;
    if (role == 0)
        proxy = new Party(proxy1, hport, haddress, cport, caddress);
    else
        proxy = new Party(proxy2, hport, haddress, cport, caddress);



    auto start = chrono::high_resolution_clock::now();
    //VSortTest(proxy);
    //SortTest(proxy);
    SortNarrowTest(proxy);
    auto end = chrono::high_resolution_clock::now();
    double time_taken =
            chrono::duration_cast<chrono::nanoseconds>(end - start).count();
    time_taken *= 1e-9;
    cout<<time_taken<<endl;


    proxy->SendBytes(coreEnd);
    //proxy->PrintBytes();

    return 0;
}
