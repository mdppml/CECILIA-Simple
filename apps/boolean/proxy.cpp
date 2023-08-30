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

constexpr int sz = 1000000;

void AND_test(Party *proxy){
    ofstream txt;
    cout<<setfill ('*')<<setw(50)<<"Calling Boolean And";
    cout<<setfill ('*')<<setw(49)<<"*"<<endl;
    uint32_t size = sz;
    cout << "size " << sz << endl;
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
    proxy->SendBytes(boolAnd, params, 1);

    auto start = chrono::high_resolution_clock::now();
    uint8_t* s = And(proxy, x, y, size);
    auto end = chrono::high_resolution_clock::now();
    double totaltime =
            chrono::duration_cast<chrono::nanoseconds>(end - start).count()*1e-9;
    cout<<totaltime<<endl;

    cout << "Callng Reconstruct..\n";
    uint8_t* anded = ReconstructBoolean(proxy, s, size);
    uint8_t* x_rec = ReconstructBoolean(proxy, x, size);
    uint8_t* y_rec = ReconstructBoolean(proxy, y, size);
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
    cout<<setfill ('*')<<setw(50)<<"Calling Boolean Subtract";
    cout<<setfill ('*')<<setw(49)<<"*"<<endl;
    uint32_t size = sz;

    auto* x =new uint64_t[size];
    auto* y =new uint64_t[size];

    for (int i = 0; i <size ; ++i) {
        x[i] = std::rand();     // This is not working with proxy->generateRandom() -> FIX?
        y[i] = std::rand();
//        if (proxy->getPRole()==P1){
//            x[i] = 9723820186448424737;
//            y[i] = 3;
//        }
//        else {
//            x[i] =9723820186448424749;
//            y[i] =1;
//        }
    }

    auto x_rec = ReconstructBoolean(proxy, x, size);
    auto y_rec = ReconstructBoolean(proxy, y, size);

    std::uint64_t * z_correct =new std::uint64_t [size];
    for (int i = 0; i < size; ++i) {
        z_correct[i] =x_rec[i]-y_rec[i];
    }

    cout << "Calling SendBytes..\n";
    uint32_t params[1];
    params[0] = size;
    proxy->SendBytes(boolSubtract, params, 1);

    auto start = chrono::high_resolution_clock::now();
    uint64_t* s = BooleanSubtract(proxy, x, y, size);
    auto end = chrono::high_resolution_clock::now();
    double totaltime =
            chrono::duration_cast<chrono::nanoseconds>(end - start).count()*1e-9;
    cout<<totaltime<<endl;

    cout << "Callng Reconstruct..\n";
    auto sub = ReconstructBoolean(proxy, s, size);
    int numCorrect = sz;
    for(int i = 0;i<size;i++) {
        if ((uint64_t) sub[i] != (uint64_t) z_correct[i]) {
            numCorrect--;
            cout << "calculated z: " << (uint64_t) sub[i] << "\tcorrect z: " << (uint64_t) z_correct[i] << endl;
        }
    }
    cout  << numCorrect << " out of " << sz << " are calculated correctly" << endl;

}
void Conversion_test(Party *proxy){
    ofstream txt;
    cout<<setfill ('*')<<setw(50)<<"Calling Arithmetic to XOR Conversion"<< setfill ('*')<<setw(49)<<"*"<<endl;
    uint32_t size = sz;

    auto* x =new uint64_t[size];
    cout << "size " << sz << endl;
    for (int i = 0; i <size ; ++i) {
        x[i] = 5;
    }

    uint32_t params[1];
    params[0] = size;
    proxy->SendBytes(boolArithmeticToXor, params, 1);

    auto start = chrono::high_resolution_clock::now();
    uint64_t* s = ArithmeticToXor(proxy, x, size);
    auto end = chrono::high_resolution_clock::now();
    double totaltime = chrono::duration_cast<chrono::nanoseconds>(end - start).count()*1e-9;
    cout<<totaltime<<endl;

    auto s_rec = ReconstructBoolean(proxy, s, size);
    auto x_rec = Reconstruct(proxy, x, size);
    int  counter = 0;
    for(int i = 0;i<size;i++) {
        if (s_rec[i]!= x_rec[i]) break;//cout << "A2B not working " << endl;
        else counter++;
    }
    if(counter == sz) cout << "A2B Conversion works correctly \n";
    else cout << "A2B Conversion is INCORRECT \n";



    cout<<setfill ('*')<<setw(50)<<"Calling XOR (64bit) to Arithmetic Conversion"<< setfill ('*')<<setw(49)<<"*"<<endl;
    cout << "Calling SendBytes..\n";
    params[0] = size;
    proxy->SendBytes(boolXorToArithmetic, params, 1);

    start = chrono::high_resolution_clock::now();
    uint64_t* t = XorToArithmetic(proxy, s, size);
    end = chrono::high_resolution_clock::now();
    totaltime = chrono::duration_cast<chrono::nanoseconds>(end - start).count()*1e-9;
    cout<<totaltime<<endl;

    cout << "Callng Reconstruct..\n";
    auto t_rec = Reconstruct(proxy, t, size);
    counter = 0;
    for(int i = 0;i<size;i++) {
        if (t_rec[i]!= x_rec[i]) cout << t_rec[i] << "\t " << x_rec[i] << endl;
        else counter++;
    }
    if(counter == sz) cout << "B2A Conversion works correctly \n";
    else cout << "Conversion is not correct!\n Only " << counter << " out of " << sz << " elements converted correctly\n";



    cout<<setfill ('*')<<setw(50)<<"Calling XOR (1Bit) to Arithmetic Conversion"<< setfill ('*')<<setw(49)<<"*"<<endl;
    cout << "Calling SendBytes..\n";
    params[0] = size;
    proxy->SendBytes(boolXorToArithmetic3, params, 1);
    int sz = size/8 +1;
    auto* b =new uint8_t[sz];
    for (int i = 0; i <sz ; ++i) {
        if(proxy->getPRole()==P1) b[i] = 5;
        else b[i] = 0;
    }

    start = chrono::high_resolution_clock::now();
    t = XorToArithmetic3(proxy, b, size);
    end = chrono::high_resolution_clock::now();
    totaltime = chrono::duration_cast<chrono::nanoseconds>(end - start).count()*1e-9;
    cout<<totaltime<<endl;

//    cout << "Callng Reconstruct..\n";
//    t_rec = Reconstruct(proxy, t, size);
//
//    counter = 0;
//    for(int i = 0;i<size;i++) {
//        cout << t[i] << endl;
//    }
//    if(counter == sz) cout << "B2A Conversion works correctly \n";
//    else cout << "Conversion is not correct!\n Only " << counter << " out of " << sz << " elements converted correctly\n";

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

    Conversion_test(proxy);
//    SUB_test(proxy);
//    AND_test(proxy);

    proxy->SendBytes(CORE_END);
    proxy->PrintBytes();

    return 0;
}
