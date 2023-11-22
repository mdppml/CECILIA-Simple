//
// Created by Seyma Selcan on 02.12.22.
//
#include <cstdlib>
#include <iostream>
#include <math.h>
#include <chrono>
#include <tuple>
#include <iomanip>
#include "../../core/sort.h"

using namespace std;

constexpr int sz = 10000;
constexpr int cols = 2;
constexpr int ringbits = 20;

void WriteResults(Party *proxy){

    string fileLine;
    vector<string> lines;
    string fileName = "/home/ubuntu/CECILIA/results1.txt";
    if(proxy->GetPRole()==proxy2) fileName = "/home/ubuntu/CECILIA/results2.txt";
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

void SortTest(Party *proxy, int sz2, int ringbits2){
    ofstream txt;
    cout<<setfill ('*')<<setw(50)<<"Calling SORT_test";
    cout<<setfill ('*')<<setw(49)<<"*"<<endl;
    uint32_t size = sz2;
    cout<<"Size: "<<size << " \t ringbits: " << ringbits2 <<endl;
    auto* a =new uint64_t[size];
    for(int i = 0; i<size; i++) {
        //a[i] = proxy->generateCommonRandom();
        if (proxy->GetPRole() == proxy1){
            a[i] = size-i-1;
        }else{
            a[i] = 0;
        }
    }
    if (proxy->GetPRole() == proxy1){
        a[2] = 377;
    }


    uint32_t params[2];
    params[0] = size;
    params[1] = ringbits2;
    proxy->SendBytes(coreSort, params, 2);
    cout << "Calling Sort..\n";
    auto start = chrono::high_resolution_clock::now();
    uint32_t* s = Sort(proxy, a, size, ringbits2);
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
    cout<<"Total time:\t"<<t_time<<endl;
    //cout<<"f Time:\t"<< perm_f <<endl;
    //cout<<"s Time:\t"<< perm_s <<endl;
    //cout<<"p Time:\t"<< perm_p <<endl;
    //cout<<"REC transfer time:\t"<< rec_transfer <<endl;
    //cout<<"REC read time:\t"<< rec_read <<endl;
    //cout<<"MUL triple gen time:\t"<< mul_triple_gen <<endl;
    //cout<<"REC ef calc:\t"<< mul_ef_calc <<endl;
    
    //WriteResults(proxy);

    cout << "Callng Reconstruct..\n";
    uint32_t* sorted = ReconstructNarrowPermutation(proxy,s,size, ringbits2);

    for(int i = 0;i<10;i++){
        cout <<  sorted[i]<< endl;
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
    int sz2 = atoi(argv[6]);
    int ringbits2 = ceil(log2(sz2));
    Party *proxy;
    if (role == 0)
        proxy = new Party(proxy1, hport, haddress, cport, caddress);
    else
        proxy = new Party(proxy2, hport, haddress, cport, caddress);



    auto start = chrono::high_resolution_clock::now();
    SortTest(proxy, sz2, ringbits2);
    auto end = chrono::high_resolution_clock::now();
    double time_taken =
            chrono::duration_cast<chrono::nanoseconds>(end - start).count();
    time_taken *= 1e-9;
    cout<<time_taken<<endl;

    proxy->SendBytes(coreEnd);
    //proxy->PrintBytes();

    return 0;
}
