//
// Created by noah on 21/06/22.
//
#include <iostream>
#include "../core/core.h"
#include "../core/cnn.h"
#include "../core/rkn.h"
#include "../core/auc.h"
#include "../core/sort.h"
#include "../booleancore/core.h"


int main(int argc, char* argv[]) {
    if (argc != 3) {
        clog << "Error: The program requires exactly two arguments; the IP address of the helper and the port it is listening on." << endl;
        return 1;
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
            case cnnMax:
                sz = proxy->ReadInt();
                Max(proxy, nullptr, sz);
                break;
            case cnnVectorisedMax:
                params[0] = proxy->ReadInt();
                params[1] = proxy->ReadInt();
                params[2] = proxy->ReadInt();
                params[3] = proxy->ReadInt();
                if (
                        params[0] > 0
                        and params[1] > 0
                        and params[2] > 0
                        and params[2] <= params[0]
                        and params[3] <= params[1]
                    ){
                    Max(proxy, nullptr, params[0], params[1], params[2], params[3], true);
                }
                else{
                    cout << "ERROR: received mmax parameters were not in valid range..." << endl;
                }
                break;
            case cnnArgMax:
                sz = proxy->ReadInt();
                ArgMax(proxy, nullptr, sz);
                break;
            case cnnRelu:
                ReLU(proxy, 0);
                break;
            case cnnVectorisedRelu:
                sz = proxy->ReadInt();
                Relu(proxy, nullptr, sz);
                break;
            case cnnDerivativeRelu:
                DerivativeRelu(proxy, 0);
                break;
            case cnnVectorisedDerivativeRelu:
                sz = proxy->ReadInt();
                DerivativeRelu(proxy, nullptr, sz);
                break;
            case cnnConvolutionalLayer:
                params[0] = proxy->ReadInt();
                params[1] = proxy->ReadInt();
                params[2] = proxy->ReadInt();
                params[3] = proxy->ReadInt();
                params[4] = proxy->ReadInt();
                params[5] = proxy->ReadInt();
                params[6] = proxy->ReadInt();
                params[7] = proxy->ReadInt();
                params[8] = proxy->ReadInt();
                ConvolutionalLayer(proxy, nullptr, params[0], params[1], params[2], nullptr, params[3], params[4],
                                   params[5], params[6],
                                   params[7], nullptr, params[8]);
                break;
            case cnnFullyConnectedLayer:
                params[0] = proxy->ReadInt();
                params[1] = proxy->ReadInt();
                FullyConnectedLayer(proxy, nullptr, params[0], nullptr, params[1], nullptr);
                break;
            case rknGaussianKernel:
                n_gms = proxy->ReadInt();
                sz = proxy->ReadInt();
                GaussianKernel(proxy, nullptr, 0, n_gms, sz);
                break;
            case rknInverseSqrt:
                sz = proxy->ReadInt();
                InverseSqrt(proxy, nullptr, sz);
                break;
            case rknVectorisedInverseSqrt:
                n_gms = proxy->ReadInt();
                sz = proxy->ReadInt();
                InverseSqrt(proxy, nullptr, n_gms, sz);
                break;
            case rknIteration:
                size1 = proxy->ReadInt();
                size2 = proxy->ReadInt();
                RknIteration(proxy, nullptr, nullptr, nullptr, size1, size2, 0, 0, 0);
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
            case aucMostSignificantBit:
                sz = proxy->ReadInt();
                AucMostSignificantBit(proxy, nullptr, sz);
                break;
            case aucDivide:
                AucDivide(proxy, 0, 0);
                break;
            case aucVectorisedRound:
                sz = proxy->ReadInt();
                Round(proxy, nullptr, sz);
                break;
            case aucVectorisedDivide:
                sz = proxy->ReadInt();
                AucDivide(proxy, nullptr, nullptr, sz);
                break;
            case aucRocNoTie:
                sz = proxy->ReadInt();
                cout << "Sz: " << sz << endl;
                RocNoTie(proxy, 0, sz);
                break;
            case aucRocWithTie:
                sz = proxy->ReadInt();
                RocWithTie(proxy, 0, sz);
                break;
            case aucPrCurve:
                sz = proxy->ReadInt();
                PrCurve(proxy, 0, sz);
                break;
            case rknEigenDecomposition:
                sz = proxy->ReadInt();
                EigenDecomposition(proxy, sz);
                break;
            case rknVectorisedEigenDecomposition:
                n_gms = proxy->ReadInt();
                sz = proxy->ReadInt();
                EigenDecomposition(proxy, n_gms, sz);
                break;
            case coreSort:
                sz = proxy->ReadInt();
                Sort(proxy, 0, sz);
                break;
            case coreVSort:
                sz = proxy->ReadInt();
                size1 = proxy->ReadInt();
                Sort(proxy, 0, sz, size1, 0);
                break;
            case coreSort2:
                sz = proxy->ReadInt();
                size1 = proxy->ReadInt();
                SortNarrow(proxy, 0, sz, size1);
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
