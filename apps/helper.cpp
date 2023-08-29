//
// Created by noah on 21/06/22.
//
#include <iostream>
#include "../core/core.h"
#include "../core/cnn.h"
#include "../core/rkn.h"
#include "../core/auc.h"
#include "../core/sort.h"


int main(int argc, char* argv[]) {
    if (argc != 3) {
        clog << "Error: The program requires exactly two arguments; the IP address of the helper and the port it is listening on." << endl;
        return 1;
    }
    string address(argv[1]);
    uint16_t port = strtol(argv[2], nullptr, 10);

    auto *helper = new Party(helper, port, address);
    bool keep_looping = true;
    uint32_t sz, n_gms, size1, size2;
    uint64_t params [9];
    Operation operation;
    auto start = chrono::high_resolution_clock::now();
    while (keep_looping){
        operation = static_cast<Operation>(helper->ReadByte());
        switch(operation) {
            case coreVectorisedMostSignificantBit:
                sz = helper->ReadInt();
                MostSignificantBit(helper, nullptr, sz);
                break;
            case coreMostSignificantBit:
                MostSignificantBit(helper, 0);
                break;
            case coreModularConversion:
                ModularConversion(helper, 0);
                break;
            case coreVectorisedModularConversion:
                sz = helper->ReadInt();
                ModularConversion(helper, nullptr, sz);
                break;
            case coreCompare:
                Compare(helper, 0, 0);
                break;
            case coreVectorisedCompare:
                sz = helper->ReadInt();
                Compare(helper, nullptr, nullptr, sz);
                break;
            case coreMultiplex:
                Multiplex(helper, 0, 0, 0);
                break;
            case coreVectorisedMultiplex:
                sz = helper->ReadInt();
                Multiplex(helper, nullptr, nullptr, nullptr, sz);
                break;
            case coreMultiply:
                Multiply(helper, 0, 0);
                break;
            case coreVectorisedMultiply:
                sz = helper->ReadInt();
                Multiply(helper, nullptr, nullptr, sz);
                break;
            case coreDotProduct:
                sz = helper->ReadInt();
                DotProduct(helper, nullptr, nullptr, sz);
                break;
            case coreVectorisedDotProduct:
                sz = helper->ReadInt();
                DotProduct(helper, nullptr, nullptr, sz, 0);
                break;
            case coreExp:
                Exp(helper, 0);
                break;
            case coreVectorisedExp:
                sz = helper->ReadInt();
                Exp(helper, nullptr, sz);
                break;
            case coreMatrixMatrixMultiply:
                sz = helper->ReadInt();
                // note that a_row is the required size of the multiplication that will be performed in MatrixMatrixMultiply
                MatrixMatrixMultiply(helper, nullptr, nullptr, sz, 0, 0);
                break;
            case coreVectorisedMatrixMatrixMultiply:
                sz = helper->ReadInt();
                // note that a_row is the required size of the multiplication that will be performed in MatrixMatrixMultiply
                MatrixMatrixMultiply(helper, nullptr, nullptr, 0, sz, 0, 0);
                break;
            case coreMatrixVectorMultiply:
                sz = helper->ReadInt();
                // note that a_row is the required size of the multiplication that will be performed in MatrixVectorMultiply
                MatrixVectorMultiply(helper, nullptr, nullptr, sz, 0);
                break;
            case coreVectorisedMatrixVectorMultiply:
                sz = helper->ReadInt();
                // note that a_row is the required size of the multiplication that will be performed in MatrixVectorMultiply
                MatrixVectorMultiply(helper, nullptr, nullptr, 0, sz, 0);
                break;
            case cnnMax:
                sz = helper->ReadInt();
                Max(helper, nullptr, sz);
                break;
            case cnnVectorisedMax:
                params[0] = helper->ReadInt();
                params[1] = helper->ReadInt();
                params[2] = helper->ReadInt();
                params[3] = helper->ReadInt();
                if (
                        params[0] > 0
                        and params[1] > 0
                        and params[2] > 0
                        and params[2] <= params[0]
                        and params[3] <= params[1]
                    ){
                    Max(helper, nullptr, params[0], params[1], params[2], params[3], true);
                }
                else{
                    cout << "ERROR: received mmax parameters were not in valid range..." << endl;
                }
                break;
            case cnnArgMax:
                sz = helper->ReadInt();
                ArgMax(helper, nullptr, sz);
                break;
            case cnnReLU:
                ReLU(helper, 0);
                break;
            case cnnVectorisedReLU:
                sz = helper->ReadInt();
                ReLU(helper, nullptr, sz);
                break;
            case cnnDerivativeReLU:
                DerivativeReLU(helper, 0);
                break;
            case cnnVectorisedDerivativeReLU:
                sz = helper->ReadInt();
                DerivativeReLU(helper, nullptr, sz);
                break;
            case cnnConvolutionalLayer:
                params[0] = helper->ReadInt();
                params[1] = helper->ReadInt();
                params[2] = helper->ReadInt();
                params[3] = helper->ReadInt();
                params[4] = helper->ReadInt();
                params[5] = helper->ReadInt();
                params[6] = helper->ReadInt();
                params[7] = helper->ReadInt();
                params[8] = helper->ReadInt();
                ConvolutionalLayer(helper, nullptr, params[0], params[1], params[2], nullptr, params[3], params[4],
                                   params[5], params[6],
                                   params[7], nullptr, params[8]);
                break;
            case cnnFullyConnectedLayer:
                params[0] = helper->ReadInt();
                params[1] = helper->ReadInt();
                FullyConnectedLayer(helper, nullptr, params[0], nullptr, params[1], nullptr);
                break;
            case rknGaussianKernel:
                n_gms = helper->ReadInt();
                sz = helper->ReadInt();
                GaussianKernel(helper, nullptr, 0, n_gms, sz);
                break;
            case rknInverseSqrt:
                sz = helper->ReadInt();
                InverseSqrt(helper, nullptr, sz);
                break;
            case rknVectorisedInverseSqrt:
                n_gms = helper->ReadInt();
                sz = helper->ReadInt();
                InverseSqrt(helper, nullptr, n_gms, sz);
                break;
            case rknIteration:
                size1 = helper->ReadInt();
                size2 = helper->ReadInt();
                RknIteration(helper, nullptr, nullptr, nullptr, size1, size2, 0, 0, 0);
                break;
            case coreDivide:
                Divide(helper, 0, 0);
                break;
            case coreVectorisedDivide:
                size1 = helper->ReadInt();
                Divide(helper, 0, 0, size1);
                break;
            case coreNormalise:
                size1 = helper->ReadInt();
                Normalize(helper, 0, 0, size1);
                break;
            case aucMostSignificantBit:
                sz = helper->ReadInt();
                aucMostSignificantBit(helper, nullptr, sz);
                break;
            case aucDivide:
                aucDivide(helper, 0, 0);
                break;
            case aucVectorisedRound:
                sz = helper->ReadInt();
                Round(helper, nullptr, sz);
                break;
            case aucVectorisedDivide:
                sz = helper->ReadInt();
                aucDivide(helper, nullptr, nullptr, sz);
                break;
            case aucRocNoTie:
                sz = helper->ReadInt();
                cout << "Sz: " << sz << endl;
                RocNoTie(helper, 0, sz);
                break;
            case aucRocWithTie:
                sz = helper->ReadInt();
                RocWithTie(helper, 0, sz);
                break;
            case aucPrCurve:
                sz = helper->ReadInt();
                PrCurve(helper, 0, sz);
                break;
            case rknEigenDecomposition:
                sz = helper->ReadInt();
                EigenDecomposition(helper, sz);
                break;
            case rknVectorisedEigenDecomposition:
                n_gms = helper->ReadInt();
                sz = helper->ReadInt();
                EigenDecomposition(helper, n_gms, sz);
                break;
            case coreSort:
                sz = helper->ReadInt();
                Sort(helper, 0, sz);
                break;
            case coreEnd:
                keep_looping = false;
                break;
        }
    }
    auto end = chrono::high_resolution_clock::now();
    helper->PrintBytes();

    double time_taken = chrono::duration_cast<chrono::nanoseconds>(end - start).count();
    time_taken *= 1e-9;
    helper->PrintPaperFriendly(time_taken);
    delete helper;
    return 0;
}