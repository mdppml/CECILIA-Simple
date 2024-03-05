//
// Created by noah on 21/06/22.
//
#include <iostream>
#include "../core/core.h"
#include "../booleancore/core.h"


int main(int argc, char* argv[]) {
    int const default_argc = 2;
    char* default_args[] = {"", "127.0.0.1", "7777"};
    if(argc == 1 ){
            argc = default_argc;
            argv = default_args;
    }
    else if (argc != 3){
        clog << "Error: The program requires exactly two arguments; the IP address of the helper and the port it is listening on." << endl;
        clog << "In case of no arguments the program will use the default IP address (localhost) and port number (7777)." << endl;
        return 1;
    }
    else{
        clog << "IP address and port number received, helper is starting..." << endl;
    }
    string address(argv[1]);
    uint16_t port = strtol(argv[2], nullptr, 10);
    auto *proxy = new Party(helper, port, address);
    bool keep_looping = true;
    uint32_t sz, n_gms, n_matrices, size1, size2;
    uint64_t params [9];
    Operation operation;
    auto start = chrono::high_resolution_clock::now();
    while (keep_looping){
        operation = static_cast<Operation>(proxy->ReadByte());
        switch(operation) {
            case coreVectorisedMostSignificantBit:
                sz = proxy->ReadInt();
                MostSignificantBit(proxy, nullptr, sz);
                break;
            case coreMostSignificantBit:
                MostSignificantBit(proxy, 0);
                break;
            case coreModularConversion:
                ModularConversion(proxy, 0);
                break;
            case coreVectorisedModularConversion:
                sz = proxy->ReadInt();
                ModularConversion(proxy, nullptr, sz);
                break;
            case coreCompare:
                Compare(proxy, 0, 0);
                break;
            case coreVectorisedCompare:
                sz = proxy->ReadInt();
                Compare(proxy, nullptr, nullptr, sz);
                break;
            case coreMultiplex:
                Multiplex(proxy, 0, 0, 0);
                break;
            case coreVectorisedMultiplex:
                sz = proxy->ReadInt();
                Multiplex(proxy, nullptr, nullptr, nullptr, sz);
                break;
            case coreVectorisedMultiplex2:
                sz = proxy->ReadInt();
                size2 = proxy->ReadInt();
                MultiplexNarrow(proxy, nullptr, nullptr, nullptr, sz,size2);
                break;
            case coreMultiply:
                Multiply(proxy, 0, 0);
                break;
            case coreVectorisedMultiply:
                sz = proxy->ReadInt();
                Multiply(proxy, nullptr, nullptr, sz);
                break;
            case coreVectorisedMultiply2:
                sz = proxy->ReadInt();
                size2 = proxy->ReadInt();
                MultiplyNarrow(proxy, nullptr, nullptr, sz, size2);
                break;
            case coreDotProduct:
                sz = proxy->ReadInt();
                DotProduct(proxy, nullptr, nullptr, sz);
                break;
            case coreVectorisedDotProduct:
                sz = proxy->ReadInt();
                DotProduct(proxy, nullptr, nullptr, sz, 0);
                break;
            case coreExp:
                Exp(proxy, 0);
                break;
            case coreVectorisedExp:
                sz = proxy->ReadInt();
                Exp(proxy, nullptr, sz);
                break;
            case coreMatrixMatrixMultiply:
                sz = proxy->ReadInt();
                size1 = proxy->ReadInt();
                size2 = proxy->ReadInt();
                MatrixMatrixMultiply(proxy, nullptr, nullptr, sz, size1, size2);
                break;
            case coreVectorisedMatrixMatrixMultiply:
                n_matrices = proxy->ReadInt();
                sz = proxy->ReadInt();
                size1 = proxy->ReadInt();
                size2 = proxy->ReadInt();
                MatrixMatrixMultiply(proxy, nullptr, nullptr, n_matrices, sz, size1, size2);
                break;
            case coreMatrixVectorMultiply:
                size1 = proxy->ReadInt();
                size2 = proxy->ReadInt();
                MatrixVectorMultiply(proxy, nullptr, nullptr, size1, size2);
                break;
            case coreVectorisedMatrixVectorMultiply:
                n_matrices = proxy->ReadInt();
                size1 = proxy->ReadInt();
                size2 = proxy->ReadInt();
                MatrixVectorMultiply(proxy, nullptr, nullptr, n_matrices, size1, size2);
                break;
            case coreDivide:
                Divide(proxy, 0, 0);
                break;
            case coreVectorisedDivide:
                size1 = proxy->ReadInt();
                Divide(proxy, 0, 0, size1);
                break;
            case coreNormalise:
                size1 = proxy->ReadInt();
                Normalise(proxy, 0, 0, size1);
                break;
            case boolAnd:
                sz = proxy->ReadInt();
                And2(proxy, 0, 0, sz);
                break;
            case boolSubtract:
                sz = proxy->ReadInt();
                BooleanSubtract2(proxy, 0, 0, sz);
                break;
            case boolArithmeticToXor:
                sz = proxy->ReadInt();
                ArithmeticToXor(proxy, 0, sz);
                break;
            case boolXorToArithmetic:
                sz = proxy->ReadInt();
                XorToArithmetic(proxy, 0, sz);
                break;
            case boolXorToArithmetic2:
                sz = proxy->ReadInt();
                XorToArithmetic2(proxy, 0, sz);
                break;
            case boolXorToArithmetic3:
                sz = proxy->ReadInt();
                size1 = proxy->ReadInt();
                XorToArithmetic3(proxy, 0, sz, size1);
                break;
            case coreEnd:
                keep_looping = false;
                break;
        }
    }
    auto end = chrono::high_resolution_clock::now();
    PrintBytes();

    double time_taken = chrono::duration_cast<chrono::nanoseconds>(end - start).count();
    time_taken *= 1e-9;
    proxy->PrintPaperFriendly(time_taken);
    delete proxy;
    return 0;
}
