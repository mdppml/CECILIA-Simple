//
// Created by Seyma Selcan on 02.12.22.
//
#include <cstdlib>
#include <iostream>
#include <chrono>
#include <tuple>
#include <random>
#include <iomanip>
#include "../../core/core.h"
#include "../../core/sort.h"
#include "../../booleancore/core.h"

using namespace std;

constexpr int sz = 100;

void AND_test(Party *proxy){
    ofstream txt;
    cout<<setfill ('*')<<setw(50)<<"Calling Boolean AND";
    cout<<setfill ('*')<<setw(49)<<"*"<<endl;
    uint32_t size = sz;

    auto* x =new uint8_t[size];
    auto* y =new uint8_t[size];

    for(int i = 0; i < size; i++) {
        x[i] = proxy->generateRandomByte()&0x1;
    }

    for(int i = 0; i < size; i++) {
        y[i] = proxy->generateRandomByte()&0x1;
    }



    cout << "Calling SendBytes..\n";
    uint32_t params[1];
    params[0] = size;
    proxy->SendBytes(CORE_AND, params, 1);

    auto start = chrono::high_resolution_clock::now();
    uint8_t* s = AND(proxy, x, y, size);
    auto end = chrono::high_resolution_clock::now();
    double totaltime =
            chrono::duration_cast<chrono::nanoseconds>(end - start).count()*1e-9;
    cout<<totaltime<<endl;

    cout << "Callng REC..\n";
    uint8_t* anded = RECB(proxy,s,size);
    uint8_t* x_rec = RECB(proxy,x,size);
    uint8_t* y_rec = RECB(proxy,y,size);
    int numCorrect = sz;
    for(int i = 0;i<size;i++) {
        auto z_correct = x_rec[i] & y_rec[i];
        if ((int) anded[i] != (int) z_correct) {
            numCorrect--;
            cout << s[i] << " " << "calculated z: " << (int) anded[i] << "\tcorrect z: " << (int) z_correct
                 << "\tx_rec: " << (int) x_rec[i] << "\ty_rec: " << (int) y_rec[i] << "\tx: " << (int) x[i] << "\ty: "
                 << (int) y[i] << endl;
        }
    }

    cout  << numCorrect << " out of " << sz << " are calculated correctly" << endl;
    delete [] x;
    delete [] y;
    delete [] s;
    delete [] x_rec;
    delete [] anded;
    delete [] y_rec;
}
void SUB_test(Party *proxy){
    ofstream txt;
    cout<<setfill ('*')<<setw(50)<<"Calling Boolean SUB";
    cout<<setfill ('*')<<setw(49)<<"*"<<endl;
    uint32_t size = sz;

    auto* x =new uint64_t[size];
    auto* y =new uint64_t[size];

    for (int i = 0; i <size ; ++i) {
        x[i] = std::rand();     // This is not working with proxy->generateRandom() -> FIX?
        y[i] = std::rand();
//        if (proxy->getPRole()==P1){
//            x[i] = 1;
//            y[i] = 4;
//        }
//        else {
//            x[i] =2;
//            y[i] =2;
//        }
    }

    auto x_rec = RECB(proxy, x, size);
    auto y_rec = RECB(proxy, y, size);
    //cout << "recreated x and y " << x_rec[0] << " " << y_rec[0] << endl;
    std::uint64_t * z_correct =new std::uint64_t [size];
    for (int i = 0; i < size; ++i) {
        z_correct[i] =x_rec[i]-y_rec[i];
    }

    cout << "Calling SendBytes..\n";
    uint32_t params[1];
    params[0] = size;
    proxy->SendBytes(CORE_BSUB, params, 1);

    auto start = chrono::high_resolution_clock::now();

    uint64_t* s = BooleanSubstract(proxy, x, y, size);
    auto end = chrono::high_resolution_clock::now();
    double totaltime =
            chrono::duration_cast<chrono::nanoseconds>(end - start).count()*1e-9;
    cout<<totaltime<<endl;

    cout << "Callng REC..\n";
    auto sub = RECB(proxy, s, size);
    int numCorrect = sz;
    for(int i = 0;i<size;i++) {
        if ((uint64_t) sub[i] != (uint64_t) z_correct[i]) {
            numCorrect--;
            cout << "calculated z: " << (uint64_t) sub[i] << "\tcorrect z: " << (uint64_t) z_correct[i] << endl;
        }
    }
    cout  << numCorrect << " out of " << sz << " are calculated correctly" << endl;


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

    SUB_test(proxy);
    //AND_test(proxy);

    proxy->SendBytes(CORE_END);
    proxy->PrintBytes();

    return 0;
}
